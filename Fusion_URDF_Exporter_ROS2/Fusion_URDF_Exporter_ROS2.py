#Author-syuntoku14, Dheena2k2, Lentin Joseph, MarceloJordao01
#Description-Generate URDF file from Fusion 360 (supports nested / linked sub-assemblies)

import adsk, adsk.core, adsk.fusion, traceback
import os
import json
from .utils import utils
from .core import Write, Assembly

"""
Units: Fusion internal cm -> URDF meters. Meshes are written in mm (scale 0.001).

What changed compared to upstream:
* links/joints are read from the WHOLE assembly tree (root.allJoints), so
  sub-assemblies (also linked ones) work and no restructuring is needed;
* Rigid joints merge parts into one link; Revolute -> revolute/continuous,
  Slider -> prismatic; parent/child is found with a BFS from the base, so the
  selection order inside Fusion does not matter;
* the base does not have to be called base_link (you are asked for it);
* the robot can be rotated on export so that it faces +X (ROS convention);
* the design is NOT modified (no copied components / renamed occurrences).
"""

TITLE = 'Fusion 360 -> ROS 2 URDF'
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
SETTINGS_FILE = os.path.join(SCRIPT_DIR, 'export_settings.json')


def _load_settings():
    try:
        with open(SETTINGS_FILE, encoding='utf-8') as f:
            return json.load(f)
    except Exception:
        return {}


def _save_settings(s):
    try:
        with open(SETTINGS_FILE, 'w', encoding='utf-8') as f:
            json.dump(s, f, indent=2, ensure_ascii=False)
    except Exception:
        pass


def _ask(ui, prompt, default):
    ret = ui.inputBox(prompt, TITLE, default)
    if isinstance(ret, (list, tuple)):
        value, cancelled = ret[0], ret[1]
    else:
        value, cancelled = ret, False
    if cancelled:
        return None
    return value.strip()


def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        design = adsk.fusion.Design.cast(app.activeProduct)
        if not design:
            ui.messageBox('No active Fusion design', TITLE)
            return

        root = design.rootComponent
        robot_name = Assembly.clean_name(root.name.split()[0])
        package_name = robot_name + '_description'
        settings = _load_settings().get(robot_name, {})

        # ---------------- questions ----------------
        has_base_link = any(o.component.name == 'base_link' for o in root.allOccurrences)
        base_hint = None
        if not has_base_link:
            base_hint = _ask(ui,
                'Nao existe componente "base_link".\n'
                'Digite o nome (ou parte do nome) de uma peca/submontagem que pertence a base\n'
                '(ex.: torax). Tudo que estiver rigidamente ligado a ela vira base_link:',
                settings.get('base_hint', ''))
            if base_hint is None:
                return

        front = _ask(ui,
            'Para qual eixo do Fusion a FRENTE do robo aponta?\n'
            '(+X, -X, +Y ou -Y). O modelo sera girado para a frente ficar em +X (padrao ROS).',
            settings.get('front', '+X'))
        if front is None:
            return
        front = front.upper().replace(' ', '')
        if front not in Assembly.FRONT_TO_YAW:
            ui.messageBox('Eixo invalido: %s' % front, TITLE)
            return

        quality = _ask(ui, 'Qualidade das malhas STL: low, medium ou high',
                       settings.get('quality', 'low'))
        if quality is None:
            return
        quality = quality.lower()

        all_settings = _load_settings()
        all_settings[robot_name] = {'base_hint': base_hint or '', 'front': front, 'quality': quality}
        _save_settings(all_settings)

        # ---------------- analysis ----------------
        mass_overrides = Assembly.load_mass_overrides(SCRIPT_DIR)
        model = Assembly.analyze(root, base_hint=base_hint, mass_overrides=mass_overrides)
        Assembly.apply_yaw(model, Assembly.FRONT_TO_YAW[front])
        rep = Assembly.report(model)

        summary = '%d links, %d juntas, massa total %.2f kg.\n' % (
            len(model['links']), len(model['joints_dict']),
            sum(d['mass'] for d in model['inertial_dict'].values()))
        if mass_overrides:
            summary += '(%d massas substituidas por mass_overrides.json)\n' % len(mass_overrides)
        if model['warnings']:
            summary += '\nAvisos:\n' + '\n'.join('- ' + w for w in model['warnings'][:12])
            if len(model['warnings']) > 12:
                summary += '\n... (%d avisos no export_report.txt)' % len(model['warnings'])
        summary += '\n\nContinuar e escolher a pasta de saida?'
        if ui.messageBox(summary, TITLE, adsk.core.MessageBoxButtonTypes.OKCancelButtonType) \
                != adsk.core.DialogResults.DialogOK:
            return

        save_dir = utils.file_dialog(ui)
        if save_dir is False:
            return
        save_dir = save_dir + '/' + package_name
        os.makedirs(save_dir, exist_ok=True)

        harmonic = ui.messageBox(
            'Voce usa Gazebo Harmonic (ROS 2 Jazzy)?\nSim = Harmonic, Nao = Gazebo Classic (Humble).',
            TITLE, adsk.core.MessageBoxButtonTypes.YesNoButtonType) == adsk.core.DialogResults.DialogYes

        joints_dict = model['joints_dict']
        inertial_dict = model['inertial_dict']
        links_xyz_dict = {}
        package_dir = SCRIPT_DIR + '/package/'

        # ---------------- write package ----------------
        if harmonic:
            Write.write_urdf_sim(joints_dict, links_xyz_dict, inertial_dict, package_name, robot_name, save_dir)
            Write.write_materials_xacro(joints_dict, links_xyz_dict, inertial_dict, package_name, robot_name, save_dir)
            Write.write_ros2control_xacro(joints_dict, links_xyz_dict, inertial_dict, package_name, robot_name, save_dir)
            Write.write_gazebo_sim_xacro(joints_dict, links_xyz_dict, inertial_dict, package_name, robot_name, save_dir)
            Write.write_display_launch(package_name, robot_name, save_dir)
            Write.write_gazebo_sim_launch(package_name, robot_name, save_dir)
        else:
            Write.write_urdf(joints_dict, links_xyz_dict, inertial_dict, package_name, robot_name, save_dir)
            Write.write_materials_xacro(joints_dict, links_xyz_dict, inertial_dict, package_name, robot_name, save_dir)
            Write.write_transmissions_xacro(joints_dict, links_xyz_dict, inertial_dict, package_name, robot_name, save_dir)
            Write.write_gazebo_xacro(joints_dict, links_xyz_dict, inertial_dict, package_name, robot_name, save_dir)
            Write.write_display_launch(package_name, robot_name, save_dir)
            Write.write_gazebo_launch(package_name, robot_name, save_dir)

        utils.create_package(package_name, save_dir, package_dir)
        utils.update_setup_py(save_dir, package_name)
        utils.update_setup_cfg(save_dir, package_name)
        utils.update_package_xml(save_dir, package_name)

        # ---------------- meshes ----------------
        total = sum(len(v) for v in model['links'].values())
        pd = ui.createProgressDialog()
        pd.isCancelButtonShown = True
        pd.show(TITLE, 'Exportando malhas: %v de %m pecas', 0, total)

        def progress(done, tot):
            pd.progressValue = done
            adsk.doEvents()
            return not pd.wasCancelled

        try:
            mesh_errors = Assembly.export_link_meshes(design, model, save_dir, quality, progress)
        finally:
            pd.hide()

        with open(os.path.join(save_dir, 'export_report.txt'), 'w', encoding='utf-8') as f:
            f.write(rep)
            if mesh_errors:
                f.write('\n\nERROS AO EXPORTAR MALHAS:\n' + '\n'.join('  - ' + e for e in mesh_errors))

        msg = 'Pacote gerado em:\n%s\n\n%d links, %d juntas (%s).' % (
            save_dir, len(model['links']), len(joints_dict),
            'Gazebo Harmonic' if harmonic else 'Gazebo Classic')
        if mesh_errors:
            msg += '\n\n%d corpos falharam na exportacao (ver export_report.txt).' % len(mesh_errors)
        ui.messageBox(msg, TITLE)

    except Exception:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()), TITLE)
