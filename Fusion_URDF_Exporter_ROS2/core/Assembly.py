# -*- coding: utf-8 -*-
"""
Graph-based reader for nested Fusion assemblies.

Unlike the original exporter (which only looks at first-level components and
root-level joints), this module:

* reads every joint in the whole assembly tree (root.allJoints), in world
  coordinates, including joints inside linked / nested sub-assemblies;
* groups every occurrence that has bodies ("leaf") into rigid links using
  Rigid joints and rigid groups (union-find);
* attaches leaves that have no joint at all (screws, nuts...) to the nearest
  leaf of the same top-level sub-assembly (bounding-box distance);
* builds the kinematic tree from Revolute / Slider joints with a BFS from the
  base link, so the order in which components were picked does not matter;
* combines mass / center of mass / inertia of every leaf of a link;
* exports one STL per link, in world coordinates (mm), WITHOUT modifying the
  design (bodies are exported one by one and transformed in Python).

The output (joints_dict / inertial_dict) uses the same format the original
Write.py expects, so the URDF / xacro / launch writers are reused unchanged.
"""

import adsk, adsk.core, adsk.fusion
import math, os, re, struct, tempfile, unicodedata, json

# Fusion JointTypes enum order
_RIGID, _REVOLUTE, _SLIDER = 0, 1, 2
_JT_NAMES = ['Rigid', 'Revolute', 'Slider', 'Cylindrical', 'PinSlot', 'Planar', 'Ball']

GROUND = '__ground__'
ROOT_BODIES = '__root__'


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def clean_name(name):
    """'pé_protótipo v25:1' -> 'pe_prototipo'"""
    n = re.sub(r'\s+v\d+(:\d+)?$', '', name)      # version suffix
    n = re.sub(r':\d+$', '', n)                     # instance suffix
    n = unicodedata.normalize('NFKD', n).encode('ascii', 'ignore').decode()
    n = re.sub(r'[^A-Za-z0-9_]+', '_', n).strip('_')
    if not n or n[0].isdigit():
        n = 'l_' + n
    return n


def _bbox(entity):
    b = entity.boundingBox
    return (b.minPoint.x, b.minPoint.y, b.minPoint.z,
            b.maxPoint.x, b.maxPoint.y, b.maxPoint.z)


def _bbox_gap(a, b):
    d = [max(0.0, a[i] - b[i + 3], b[i] - a[i + 3]) for i in range(3)]
    return math.sqrt(sum(x * x for x in d))


def _bbox_overlap(a, b):
    v = 1.0
    for i in range(3):
        s = min(a[i + 3], b[i + 3]) - max(a[i], b[i])
        if s <= 0:
            return 0.0
        v *= s
    return v


def _top(path):
    return path.split('+')[0]


def _point(geom_or_origin):
    """JointGeometry or JointOrigin -> world point (cm)."""
    try:
        return geom_or_origin.origin.asArray()
    except Exception:
        return geom_or_origin.geometry.origin.asArray()


class _UF:
    def __init__(self):
        self.p = {}

    def add(self, k):
        self.p.setdefault(k, k)

    def find(self, k):
        self.add(k)
        while self.p[k] != k:
            self.p[k] = self.p[self.p[k]]
            k = self.p[k]
        return k

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


# --------------------------------------------------------------------------
# main analysis
# --------------------------------------------------------------------------
def load_mass_overrides(script_dir):
    """Optional mass_overrides.json next to the script: {"MX-106T": 0.153, ...}
    Keys are matched against the cleaned component name (no version)."""
    f = os.path.join(script_dir, 'mass_overrides.json')
    if not os.path.exists(f):
        return {}
    with open(f, encoding='utf-8') as fh:
        raw = json.load(fh)
    return {clean_name(k): float(v) for k, v in raw.items() if not k.startswith('_')}


def analyze(root, base_hint=None, mass_overrides=None):
    """
    Returns a model dict:
      links:   {link_name: {'leaves': [occ,...], 'group': gid}}
      joints:  ordered list of dicts (BFS order)
      base:    'base_link'
      warnings: [str]
    """
    warnings = []
    mass_overrides = mass_overrides or {}

    # ---- leaves (occurrences that own bodies) ---------------------------
    leaves = {}
    for occ in root.allOccurrences:
        if occ.bRepBodies.count > 0:
            leaves[occ.fullPathName] = occ
    has_root_bodies = root.bRepBodies.count > 0

    uf = _UF()
    for k in leaves:
        uf.add(k)
    if has_root_bodies:
        uf.add(ROOT_BODIES)

    def leaves_under(occ):
        out = []
        stack = [occ]
        while stack:
            o = stack.pop()
            if o.bRepBodies.count > 0:
                out.append(o.fullPathName)
            for c in o.childOccurrences:
                stack.append(c)
        return out

    def node_of(occ, jname):
        if occ is None:
            return ROOT_BODIES if has_root_bodies else GROUND
        if occ.bRepBodies.count > 0:
            return occ.fullPathName
        sub = leaves_under(occ)
        if not sub:
            warnings.append('Junta "%s": "%s" nao tem corpos; ignorada.' % (jname, occ.fullPathName))
            return None
        for s in sub[1:]:
            uf.union(sub[0], s)
        return sub[0]

    # ---- joints ---------------------------------------------------------
    touched = set()
    movable = []
    all_joints = list(root.allJoints) + list(root.allAsBuiltJoints)
    for j in all_joints:
        try:
            if j.isSuppressed:
                continue
        except Exception:
            pass
        jt = j.jointMotion.jointType
        n1 = node_of(j.occurrenceOne, j.name)
        n2 = node_of(j.occurrenceTwo, j.name)
        if n1 is None or n2 is None:
            continue
        uf.add(n1); uf.add(n2)
        touched.update((n1, n2))
        if jt == _RIGID:
            uf.union(n1, n2)
        elif jt in (_REVOLUTE, _SLIDER):
            movable.append((j, n1, n2))
        else:
            warnings.append('Junta "%s" (%s) nao existe em URDF; tratada como fixa.'
                            % (j.name, _JT_NAMES[jt]))
            uf.union(n1, n2)

    try:
        for g in root.allRigidGroups:
            occs = [o for o in g.occurrences]
            keys = []
            for o in occs:
                k = node_of(o, g.name)
                if k:
                    keys.append(k)
            for k in keys[1:]:
                uf.union(keys[0], k)
            touched.update(keys)
    except Exception:
        pass

    # ---- orphan leaves (no joint at all): attach to nearest leaf ---------
    bboxes = {k: _bbox(o) for k, o in leaves.items()}
    orphans = [k for k in leaves if k not in touched]
    anchored = [k for k in leaves if k in touched]
    n_orph = 0
    for k in orphans:
        same_top = [a for a in anchored if _top(a) == _top(k)] or anchored
        if not same_top:
            continue
        best = min(same_top, key=lambda a: (_bbox_gap(bboxes[k], bboxes[a]),
                                            -_bbox_overlap(bboxes[k], bboxes[a])))
        uf.union(best, k)
        n_orph += 1
    if n_orph:
        warnings.append('%d pecas sem nenhuma junta (parafusos etc.) foram presas a peca mais proxima.' % n_orph)

    # ---- base group -----------------------------------------------------
    def groups():
        g = {}
        for k in list(uf.p.keys()):
            g.setdefault(uf.find(k), []).append(k)
        return g

    base_key = None
    for k, o in leaves.items():
        if o.component.name == 'base_link':
            base_key = k
            break
    if base_key is None and base_hint:
        h = base_hint.lower()
        cands = [k for k in leaves if h in k.lower()]
        if cands:
            base_key = cands[0]
    if base_key is None and GROUND in uf.p:
        base_key = GROUND
    if base_key is None:
        # heaviest group
        masses = {}
        for k, o in leaves.items():
            r = uf.find(k)
            masses[r] = masses.get(r, 0.0) + o.getPhysicalProperties(
                adsk.fusion.CalculationAccuracy.LowCalculationAccuracy).mass
        base_key = max(masses, key=masses.get)
        warnings.append('Base escolhida automaticamente pelo grupo mais pesado.')
    if GROUND in uf.p:
        uf.union(base_key, GROUND)

    # ---- BFS over movable joints; attach unreachable groups --------------
    while True:
        base_g = uf.find(base_key)
        adj = {}
        for (j, n1, n2) in movable:
            g1, g2 = uf.find(n1), uf.find(n2)
            adj.setdefault(g1, []).append((j, n1, n2))
            adj.setdefault(g2, []).append((j, n1, n2))
        seen = {base_g}
        order = [base_g]
        tree = []           # (joint, parent_g, child_g, child_node, flipped)
        loops = []
        used = set()
        qi = 0
        while qi < len(order):
            g = order[qi]; qi += 1
            for (j, n1, n2) in adj.get(g, []):
                if id(j) in used:
                    continue
                g1, g2 = uf.find(n1), uf.find(n2)
                if g1 == g2:
                    used.add(id(j))
                    loops.append(j)
                    continue
                other, child_node, flipped = (g1, n1, False) if g2 == g else (g2, n2, True)
                used.add(id(j))
                if other in seen:
                    loops.append(j)
                    continue
                seen.add(other)
                order.append(other)
                tree.append((j, g, other, child_node, flipped))
        all_g = groups()
        unreached = [g for g in all_g if g not in seen and any(k in leaves for k in all_g[g])]
        if not unreached:
            break
        # attach each unreachable group rigidly to the nearest reached leaf
        reached_leaves = [k for g in seen for k in all_g.get(g, []) if k in leaves]
        if not reached_leaves:
            warnings.append('Base sem pecas; nao foi possivel montar a arvore.')
            break
        for g in unreached:
            gl = [k for k in all_g[g] if k in leaves]
            best = None
            for a in gl:
                for b in reached_leaves:
                    d = _bbox_gap(bboxes[a], bboxes[b])
                    if best is None or d < best[0]:
                        best = (d, a, b)
            if best:
                uf.union(best[2], best[1])
                warnings.append('Grupo sem junta com o resto do robo (%s, %d pecas) foi fixado em "%s".'
                                % (_short(gl[0]), len(gl), _short(best[2])))
        movable = [m for m in movable]  # recompute with new unions

    for j in loops:
        warnings.append('Junta "%s" fecha uma cadeia (loop) ou liga pecas ja rigidas entre si; ignorada.'
                        % _jname(j))

    # ---- names ----------------------------------------------------------
    all_g = groups()
    group_leaves = {g: [leaves[k] for k in ks if k in leaves] for g, ks in all_g.items()}
    link_of_group = {uf.find(base_key): 'base_link'}
    used_names = {'base_link'}
    raw_names = {}
    for (j, pg, cg, child_node, flipped) in tree:
        comp = clean_name(leaves[child_node].component.name) if child_node in leaves else 'link'
        # short prefix from the top-level sub-assembly: 'perna_direita_montagem_x' -> 'perna_direita'
        top = '_'.join(clean_name(_top(child_node)).split('_')[:2]) if '+' in child_node else ''
        raw_names[cg] = (comp, top)
    for (j, pg, cg, child_node, flipped) in tree:
        comp, top = raw_names[cg]
        name = '%s_%s' % (top, comp) if top else comp
        base_name, i = name, 2
        while name in used_names:
            name = '%s_%d' % (base_name, i); i += 1
        used_names.add(name)
        link_of_group[cg] = name

    # ---- joints_dict ----------------------------------------------------
    joints_dict = {}
    for (j, pg, cg, child_node, flipped) in tree:
        jm = j.jointMotion
        jt = jm.jointType
        d = {'parent': link_of_group[pg], 'child': link_of_group[cg],
             'axis': [0, 0, 0], 'upper_limit': 0.0, 'lower_limit': 0.0}
        xyz = _point(j.geometryOrOriginOne)
        d['xyz'] = [round(v / 100.0, 6) for v in xyz]
        sign = -1.0 if flipped else 1.0
        if jt == _REVOLUTE:
            ax = jm.rotationAxisVector
            ax.normalize()
            d['axis'] = [round(sign * v, 6) + 0.0 for v in ax.asArray()]
            lim = jm.rotationLimits
            cur = jm.rotationValue
            if lim.isMaximumValueEnabled and lim.isMinimumValueEnabled:
                lo, hi = lim.minimumValue - cur, lim.maximumValue - cur
                if flipped:
                    lo, hi = -hi, -lo
                d['type'] = 'revolute'
                d['lower_limit'], d['upper_limit'] = round(lo, 6), round(hi, 6)
            else:
                if lim.isMaximumValueEnabled or lim.isMinimumValueEnabled:
                    warnings.append('Junta "%s" tem so um limite; exportada como continuous.' % _jname(j))
                d['type'] = 'continuous'
        else:  # slider
            ax = jm.slideDirectionVector
            ax.normalize()
            d['axis'] = [round(sign * v, 6) + 0.0 for v in ax.asArray()]
            lim = jm.slideLimits
            cur = jm.slideValue
            d['type'] = 'prismatic'
            if lim.isMaximumValueEnabled and lim.isMinimumValueEnabled:
                lo, hi = (lim.minimumValue - cur) / 100.0, (lim.maximumValue - cur) / 100.0
                if flipped:
                    lo, hi = -hi, -lo
                d['lower_limit'], d['upper_limit'] = round(lo, 6), round(hi, 6)
            else:
                warnings.append('Slider "%s" sem limites; usando +-0.1 m.' % _jname(j))
                d['lower_limit'], d['upper_limit'] = -0.1, 0.1
        d['fusion_name'] = _jname(j)
        joints_dict[link_of_group[cg] + '_joint'] = d

    links = {link_of_group[g]: group_leaves[g] for g in link_of_group}
    inertial_dict = {name: _inertial(occs, mass_overrides) for name, occs in links.items()}

    return {'links': links, 'joints_dict': joints_dict, 'inertial_dict': inertial_dict,
            'warnings': warnings}


def _short(path):
    parts = path.split('+')
    return parts[0] if len(parts) == 1 else parts[0] + ' > ' + parts[-1]


def _jname(j):
    try:
        ctx = j.assemblyContext
        return (ctx.name + ' > ' + j.name) if ctx else j.name
    except Exception:
        return j.name


# --------------------------------------------------------------------------
# global re-orientation (ROS: X forward, Y left, Z up)
# --------------------------------------------------------------------------
FRONT_TO_YAW = {'+X': 0.0, '-X': 180.0, '+Y': -90.0, '-Y': 90.0}


def yaw_matrix(deg):
    a = math.radians(deg)
    c, s = round(math.cos(a), 12), round(math.sin(a), 12)
    return [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]


def _mv(R, v):
    return [R[i][0] * v[0] + R[i][1] * v[1] + R[i][2] * v[2] for i in range(3)]


def apply_yaw(model, deg):
    """Rotate the whole robot about world Z (positions, axes, COM, inertia)."""
    if not deg:
        return
    R = yaw_matrix(deg)
    for d in model['joints_dict'].values():
        d['xyz'] = [round(x, 6) + 0.0 for x in _mv(R, d['xyz'])]
        d['axis'] = [round(x, 6) + 0.0 for x in _mv(R, d['axis'])]
    for d in model['inertial_dict'].values():
        d['center_of_mass'] = [round(x, 6) + 0.0 for x in _mv(R, d['center_of_mass'])]
        xx, yy, zz, xy, yz, xz = d['inertia']
        I = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]
        RI = [[sum(R[i][k] * I[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        T = [[sum(RI[i][k] * R[j][k] for k in range(3)) for j in range(3)] for i in range(3)]
        d['inertia'] = [round(T[0][0], 9), round(T[1][1], 9), round(T[2][2], 9),
                        round(T[0][1], 9), round(T[1][2], 9), round(T[0][2], 9)]
    model['yaw'] = deg


# --------------------------------------------------------------------------
# inertia
# --------------------------------------------------------------------------
def _inertial(occs, mass_overrides):
    """Combine leaves: returns {'mass', 'center_of_mass' (m, world), 'inertia' (kg m^2 about COM)}."""
    acc = adsk.fusion.CalculationAccuracy.HighCalculationAccuracy
    parts = []
    for o in occs:
        if not o.isVisible:
            continue
        p = o.getPhysicalProperties(acc)
        m = p.mass
        if m <= 0:
            continue
        c = [v / 100.0 for v in p.centerOfMass.asArray()]   # m
        ok, i1, i2, i3 = p.getPrincipalMomentsOfInertia()
        ok2, a1, a2, a3 = p.getPrincipalAxes()
        I = [[0.0] * 3 for _ in range(3)]
        for mom, ax in ((i1, a1), (i2, a2), (i3, a3)):
            v = ax.asArray()
            n = math.sqrt(sum(x * x for x in v)) or 1.0
            v = [x / n for x in v]
            for r in range(3):
                for s in range(3):
                    I[r][s] += mom / 10000.0 * v[r] * v[s]     # kg cm^2 -> kg m^2
        key = clean_name(o.component.name)
        if key in mass_overrides:
            k = mass_overrides[key] / m
            m *= k
            I = [[x * k for x in row] for row in I]
        parts.append((m, c, I))

    M = sum(p[0] for p in parts)
    if M <= 0:
        return {'mass': 1e-6, 'center_of_mass': [0, 0, 0], 'inertia': [1e-9, 1e-9, 1e-9, 0, 0, 0]}
    C = [sum(p[0] * p[1][i] for p in parts) / M for i in range(3)]
    T = [[0.0] * 3 for _ in range(3)]
    for m, c, I in parts:
        d = [c[i] - C[i] for i in range(3)]
        dd = sum(x * x for x in d)
        for r in range(3):
            for s in range(3):
                T[r][s] += I[r][s] + m * ((dd if r == s else 0.0) - d[r] * d[s])
    return {'mass': round(M, 6),
            'center_of_mass': [round(x, 6) for x in C],
            'inertia': [round(T[0][0], 9), round(T[1][1], 9), round(T[2][2], 9),
                        round(T[0][1], 9), round(T[1][2], 9), round(T[0][2], 9)]}


# --------------------------------------------------------------------------
# STL export (no design modification)
# --------------------------------------------------------------------------
_REFINE = {
    'low': adsk.fusion.MeshRefinementSettings.MeshRefinementLow,
    'medium': adsk.fusion.MeshRefinementSettings.MeshRefinementMedium,
    'high': adsk.fusion.MeshRefinementSettings.MeshRefinementHigh,
}


def _read_stl(path):
    with open(path, 'rb') as f:
        data = f.read()
    if data[:5] == b'solid' and b'facet' in data[:400]:
        tris = []
        verts = [tuple(map(float, l.split()[1:4])) for l in data.decode(errors='ignore').splitlines()
                 if l.strip().startswith('vertex')]
        for i in range(0, len(verts) - 2, 3):
            tris.append(verts[i:i + 3])
        return tris
    n = struct.unpack('<I', data[80:84])[0]
    out = []
    for t in struct.iter_unpack('<12fH', data[84:84 + n * 50]):
        out.append(((t[3], t[4], t[5]), (t[6], t[7], t[8]), (t[9], t[10], t[11])))
    return out


def export_link_meshes(design, model, save_dir, refinement='low', progress=None):
    """One binary STL per link, world coordinates in mm (URDF uses scale 0.001)."""
    mesh_dir = os.path.join(save_dir, 'meshes')
    os.makedirs(mesh_dir, exist_ok=True)
    em = design.exportManager
    tmp = os.path.join(tempfile.gettempdir(), 'f2u_body.stl')
    errors = []
    total = sum(len(v) for v in model['links'].values())
    done = 0
    Y = yaw_matrix(model.get('yaw', 0.0))
    for name, occs in model['links'].items():
        tris_out = []
        for o in occs:
            done += 1
            if progress and progress(done, total) is False:
                raise RuntimeError('Exportacao cancelada.')
            if not o.isVisible:
                continue
            M = o.transform2.asArray()          # row-major 4x4, translation in cm
            R0 = [[M[0], M[1], M[2]], [M[4], M[5], M[6]], [M[8], M[9], M[10]]]
            R = [[sum(Y[i][k] * R0[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
            t = _mv(Y, [M[3] * 10.0, M[7] * 10.0, M[11] * 10.0])
            for bi in range(o.bRepBodies.count):
                body = o.bRepBodies.item(bi)
                if not body.isVisible:
                    continue
                try:
                    opt = em.createSTLExportOptions(body, tmp)
                    opt.isBinaryFormat = True
                    opt.sendToPrintUtility = False
                    opt.meshRefinement = _REFINE.get(refinement, _REFINE['low'])
                    if not em.execute(opt):
                        raise RuntimeError('export failed')
                    for tri in _read_stl(tmp):
                        tris_out.append(tuple(
                            (R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2] + t[0],
                             R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2] + t[1],
                             R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2] + t[2])
                            for v in tri))
                except Exception as e:
                    errors.append('%s / %s: %s' % (o.fullPathName, body.name, e))
        _write_stl(os.path.join(mesh_dir, name + '.stl'), tris_out, name)
    try:
        os.remove(tmp)
    except Exception:
        pass
    return errors


def _write_stl(path, tris, name):
    with open(path, 'wb') as f:
        f.write(name.encode('ascii', 'ignore')[:80].ljust(80, b' '))
        f.write(struct.pack('<I', len(tris)))
        for a, b, c in tris:
            u = [b[i] - a[i] for i in range(3)]
            v = [c[i] - a[i] for i in range(3)]
            n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
            ln = math.sqrt(sum(x * x for x in n)) or 1.0
            f.write(struct.pack('<12fH', n[0] / ln, n[1] / ln, n[2] / ln,
                                a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], 0))


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------
def report(model):
    lines = ['LINKS (%d):' % len(model['links'])]
    tot = 0.0
    for name, occs in model['links'].items():
        m = model['inertial_dict'][name]['mass']
        tot += m
        lines.append('  %-40s %3d pecas  %.3f kg' % (name, len(occs), m))
    lines.append('  massa total: %.3f kg' % tot)
    lines.append('')
    lines.append('JOINTS (%d):' % len(model['joints_dict']))
    for n, d in model['joints_dict'].items():
        lines.append('  %-40s %-10s %s -> %s  axis=%s  xyz=%s  [%s]' % (
            n, d['type'], d['parent'], d['child'], d['axis'], d['xyz'], d.get('fusion_name', '')))
    if model['warnings']:
        lines.append('')
        lines.append('AVISOS:')
        lines += ['  - ' + w for w in model['warnings']]
    return '\n'.join(lines)
