"""Open an isolated, ready-to-run street Import timing scene for the user.

Launch interactive Blender with this chat's profile, factory startup and
autoexec disabled. By default it waits for the user's click; an explicit
--auto-run-after flag is available for controlled window-focus experiments. Its temporary
Street Test panel invokes the normal Import operator and records light CPU and
foreground-window observations. Its default is system-managed scheduling;
--qos high explicitly scopes a test-process request around Import only.
"""
from collections import Counter
from contextlib import ExitStack
import argparse
import ctypes
from ctypes import wintypes
import importlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
import traceback
from unittest.mock import patch

import bpy
from mathutils import Quaternion, Vector

ROOT = Path(__file__).resolve().parents[1]
STATE = {}


class BLOSM_TEST_OT_Run(bpy.types.Operator):
    bl_idname = 'blosm_test.run_streets'
    bl_label = 'Run street import'
    bl_description = 'Run the normal 3D realistic Import and record time and foreground/core observations'
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        from blosm.tests.street_execution_diagnostics import ProcessorState, ProcessQoS
        processor = ProcessorState()
        user = ctypes.WinDLL('user32', use_last_error=True)
        user.GetForegroundWindow.restype = wintypes.HWND
        user.GetWindowThreadProcessId.argtypes = [wintypes.HWND, ctypes.POINTER(wintypes.DWORD)]

        def sample():
            active_pid = wintypes.DWORD()
            user.GetWindowThreadProcessId(user.GetForegroundWindow(), ctypes.byref(active_pid))
            cpu = processor.current()
            return (active_pid.value == os.getpid(), active_pid.value, cpu.get('cpu'), cpu.get('efficiency_class'))

        prefix = 'blosm.way.street_network.python_backend.'
        charts = importlib.import_module(prefix + 'quad_charts')
        original = charts.visibility_center
        counts = Counter()

        def observed(*args, **kwargs):
            counts[str(sample())] += 1
            try:
                return original(*args, **kwargs)
            finally:
                counts[str(sample())] += 1

        report = dict(checkout=str(ROOT), revision=STATE['revision'], blender=bpy.app.version_string,
            pid=os.getpid(), settings=settings(), topology=processor.topology,
            policy_before=processor.power(), current_before=sample(),
            priority_before=processor.current(), requested_qos=STATE.get('qos', 'system'),
            instrumentation='Visibility-kernel entry/exit observations only; no affinity override')
        before_collections = {c.as_pointer() for c in bpy.data.collections}
        request = None
        started = time.time()
        (STATE['output'] / 'started.json').write_text(json.dumps(dict(pid=os.getpid(),
            unix_time=started, requested_qos=STATE.get('qos', 'system'))), encoding='utf8')
        start, cpu_start = time.perf_counter(), time.process_time()
        try:
            with ExitStack() as stack:
                if STATE.get('qos', 'system') == 'high':
                    request = stack.enter_context(ProcessQoS('high'))
                report['policy_during'] = processor.power()
                for module in list(sys.modules.values()):
                    if (getattr(module, '__name__', '').startswith(prefix)
                            and getattr(module, 'visibility_center', None) is original):
                        stack.enter_context(patch.object(module, 'visibility_center', observed))
                try:
                    status = bpy.ops.blosm.import_data()
                finally:
                    # Also restore aliases imported lazily during generation.
                    for module in list(sys.modules.values()):
                        if (getattr(module, '__name__', '').startswith(prefix)
                                and getattr(module, 'visibility_center', None) is observed):
                            module.visibility_center = original
            report['operator_seconds'] = time.perf_counter() - start
            report['cpu_seconds'] = time.process_time() - cpu_start
            report['status'] = sorted(status)
            for collection in bpy.data.collections:
                if collection.as_pointer() not in before_collections and collection.get('blosm_street_report'):
                    report['importer'] = json.loads(bpy.data.texts[collection['blosm_street_report']].as_string())
            report['passed'] = status == {'FINISHED'} and report.get('importer', {}).get('coverage') == 'COMPLETE'
        except Exception as error:
            report.update(passed=False, error=str(error), traceback=traceback.format_exc(),
                          operator_seconds=time.perf_counter() - start)
        report.update(processor_foreground_checks=dict(counts), current_after=sample(),
                      policy_after=processor.power(), priority_after=processor.current(),
                      execution_policy=request.report if request else None)
        number = 1
        while (STATE['output'] / f'foreground_test_{number:03d}.json').exists():
            number += 1
        receipt = STATE['output'] / f'foreground_test_{number:03d}.json'
        receipt.write_text(json.dumps(report, indent=2), encoding='utf8')
        STATE['last'] = report
        STATE['receipt'] = str(receipt)
        self.report({'INFO' if report['passed'] else 'ERROR'},
                    f'Street test: {report["operator_seconds"]:.2f} seconds; report: {receipt.name}')
        print('STREET_USER_TEST_RESULT', str(receipt), report['operator_seconds'], flush=True)
        if STATE.get('save_result'):
            bpy.ops.wm.save_as_mainfile(filepath=str(STATE['output'] / 'street_result.blend'))
        if STATE.get('quit_after_run'):
            def quit_test():
                bpy.ops.wm.quit_blender()
                return None
            bpy.app.timers.register(quit_test, first_interval=.5)
        return {'FINISHED'} if report['passed'] else {'CANCELLED'}


class BLOSM_TEST_OT_Reset(bpy.types.Operator):
    bl_idname = 'blosm_test.reset_streets'
    bl_label = 'Reset test scene'
    bl_description = 'Reload this isolated test copy, removing its generated outputs before another measurement'

    def execute(self, context):
        bpy.ops.wm.open_mainfile(filepath=str(STATE['copy']), load_ui=False, use_scripts=False)
        return {'FINISHED'}


class BLOSM_TEST_PT_Performance(bpy.types.Panel):
    bl_idname = 'BLOSM_TEST_PT_performance'
    bl_label = 'Street import timing test'
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Street Test'

    def draw(self, context):
        layout = self.layout
        props = context.scene.blosm
        layout.label(text='Isolated streets development checkout')
        layout.label(text='Revision: ' + STATE['revision'][:8])
        layout.label(text='3D realistic | Python | streets only')
        layout.label(text=f'Terrain: {props.terrainObject} / {props.streetTerrainMode}')
        layout.label(text=f'Resolution {props.streetResolution} / cell {props.streetCellSize:.2f} m')
        layout.label(text='HighQoS during Import' if STATE.get('qos') == 'high' else 'System-managed CPU scheduling')
        layout.separator()
        layout.label(text='Keep this Blender active for the test.')
        row = layout.row()
        row.scale_y = 1.5
        row.operator('blosm_test.run_streets', icon='PLAY')
        layout.operator('blosm_test.reset_streets', icon='FILE_REFRESH')
        last = STATE.get('last')
        if last:
            layout.separator()
            layout.label(text=f'Last test: {last["operator_seconds"]:.2f} seconds')
            layout.label(text='Complete street coverage' if last['passed'] else 'See diagnostic report')
        layout.separator()
        layout.label(text='Normal settings are in the Blosm tab.')


def settings():
    props = bpy.context.scene.blosm
    keys = ('mode', 'dataType', 'osmSource', 'osmFilepath', 'buildings', 'highways', 'forests',
            'vegetation', 'aviation', 'railways', 'water', 'streetBackend', 'terrainObject',
            'streetCellSize', 'streetResolution', 'streetProfileRise', 'streetExaggeration',
            'streetTerrainMode', 'streetTerrainSampling', 'streetTerrainClearance', 'assetLibraryMode')
    result = {key: getattr(props, key) for key in keys}
    result['preferPerformanceDuringImport'] = bpy.context.preferences.addons['blosm'].preferences.preferPerformanceDuringImport
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('scene', 'osm', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--auto-run-after', type=float, default=0., help='Explicit automation only: delay after readiness before Import')
    parser.add_argument('--save-result', action='store_true')
    parser.add_argument('--quit-after-run', action='store_true')
    parser.add_argument('--qos', choices=('system', 'high'), default='system')
    parser.add_argument('--background-target-pid', type=int, help='Wait for this owned idle process to be foreground before automated Import')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    STATE.update(output=args.output.resolve(), revision=args.revision,
                 copy=args.output.resolve() / 'street_foreground_test.blend',
                 save_result=args.save_result, quit_after_run=args.quit_after_run, qos=args.qos)
    if STATE['copy'].exists():
        raise RuntimeError('Use a fresh output directory; do not overwrite the existing test copy')
    spec = importlib.util.spec_from_file_location('blosm', ROOT / '__init__.py',
                                                submodule_search_locations=[str(ROOT)])
    package = importlib.util.module_from_spec(spec)
    sys.modules['blosm'] = package
    spec.loader.exec_module(package)
    package.register()
    bpy.context.preferences.view.show_splash = False
    preferences = bpy.context.preferences.addons.new()
    preferences.module = package.__name__
    preferences.preferences.dataDir = str(STATE['output'])
    preferences.preferences.enableExperimentalFeatures = True
    # This harness compares external --qos modes; do not silently add the
    # production preference to the system-managed control.
    preferences.preferences.preferPerformanceDuringImport = False
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()), load_ui=False, use_scripts=False)
    props = bpy.context.scene.blosm
    for name, value in dict(mode='3Drealistic', dataType='osm', osmSource='file',
            osmFilepath=str(args.osm.resolve()), buildings=False, forests=False,
            vegetation=False, aviation=False, highways=True, railways=False,
            water=False, loadMissingMembers=False, setupScript='', streetBackend='python',
            streetResolution=1, streetCellSize=.7, streetProfileRise=.08, streetExaggeration=1.,
            streetTerrainMode='PROFILE', streetTerrainSampling=1., streetTerrainClearance=.02,
            streetTerrainGrade=.1, streetTerrainLift=25., streetMaterialFile='',
            assetLibraryMode='AUTO').items():
        setattr(props, name, value)
    props.streetCarriagewayMaterial = None
    assert props.terrainObject == 'Terrain' and bpy.context.scene.objects.get('Terrain')
    scene = bpy.context.scene
    scene.render.threads_mode = 'FIXED'
    scene.render.threads = 16
    for cls in (BLOSM_TEST_OT_Run, BLOSM_TEST_OT_Reset, BLOSM_TEST_PT_Performance):
        bpy.utils.register_class(cls)
    corners = [bpy.data.objects['Terrain'].matrix_world @ Vector(v) for v in bpy.data.objects['Terrain'].bound_box]
    center = sum(corners, Vector()) / len(corners)

    def finish_ui_setup():
        # Sidebar visibility triggers an editor update requiring a live window.
        # Never change this property on inactive screens during CLI startup.
        windows = list(bpy.context.window_manager.windows)
        if not windows:
            return .25
        for window in windows:
            for area in window.screen.areas:
                if area.type != 'VIEW_3D':
                    continue
                with bpy.context.temp_override(window=window, area=area):
                    space = area.spaces.active
                    space.show_region_ui = True
                    space.shading.type = 'SOLID'
                    space.region_3d.view_perspective = 'ORTHO'
                    space.region_3d.view_rotation = Quaternion((1., 0., 0., 0.))
                    space.region_3d.view_location = center
                    space.region_3d.view_distance = 1200.
                    for region in area.regions:
                        if region.type == 'UI' and hasattr(region, 'active_panel_category'):
                            try:
                                region.active_panel_category = 'Street Test'
                            except (AttributeError, TypeError, ValueError, RuntimeError):
                                pass
        bpy.context.preferences.filepaths.file_preview_type = 'NONE'
        with bpy.context.temp_override(window=windows[0]):
            bpy.ops.wm.save_as_mainfile(filepath=str(STATE['copy']))
        ready = dict(checkout=str(ROOT), revision=args.revision, pid=os.getpid(), scene=str(STATE['copy']),
                     settings=settings(), auto_import=args.auto_run_after > 0,
                     auto_run_after=args.auto_run_after, qos=args.qos,
                     background_target_pid=args.background_target_pid,
                     numeric_threads={key: os.environ.get(key) for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')})
        (STATE['output'] / 'ready.json').write_text(json.dumps(ready, indent=2), encoding='utf8')
        print('STREET_USER_TEST_READY', str(STATE['copy']), flush=True)
        if args.auto_run_after > 0:
            def start_test():
                if args.background_target_pid:
                    user = ctypes.WinDLL('user32', use_last_error=True)
                    user.GetForegroundWindow.restype = wintypes.HWND
                    user.GetWindowThreadProcessId.argtypes = [wintypes.HWND, ctypes.POINTER(wintypes.DWORD)]
                    active = wintypes.DWORD()
                    user.GetWindowThreadProcessId(user.GetForegroundWindow(), ctypes.byref(active))
                    if active.value != args.background_target_pid:
                        return .25
                window = bpy.context.window_manager.windows[0]
                area = next(a for a in window.screen.areas if a.type == 'VIEW_3D')
                with bpy.context.temp_override(window=window, area=area):
                    bpy.ops.blosm_test.run_streets()
                return None
            bpy.app.timers.register(start_test, first_interval=args.auto_run_after)
        return None

    bpy.app.timers.register(finish_ui_setup, first_interval=1.)


if __name__ == '__main__':
    main()
