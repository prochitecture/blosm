"""Native Blender qualification of the production Import performance scopes."""
import argparse
import importlib
import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
from unittest.mock import patch

import bpy

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = 'bl_ext.import_performance_check.blosm'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    for name in ('bl_ext', 'bl_ext.import_performance_check'):
        if name not in sys.modules:
            shell = ModuleType(name)
            shell.__path__ = []
            sys.modules[name] = shell
    spec = importlib.util.spec_from_file_location(PACKAGE, ROOT / '__init__.py',
                                                submodule_search_locations=[str(ROOT)])
    addon = importlib.util.module_from_spec(spec)
    sys.modules[PACKAGE] = addon
    spec.loader.exec_module(addon)
    addon.register()
    entry = bpy.context.preferences.addons.new()
    entry.module = PACKAGE
    prefs = entry.preferences
    prefs.dataDir = str(args.output)
    prefs.enableExperimentalFeatures = True
    module = importlib.import_module(PACKAGE + '.app.import_performance')
    diagnostics = importlib.import_module(PACKAGE + '.tests.street_execution_diagnostics')
    runtime = importlib.import_module(PACKAGE + '.progressive.blender')
    api = module._WindowsProcessPolicy()
    original = api.read()
    state = diagnostics.ProcessorState()
    priority = state.current()
    assert not (original[0] & original[1] & 1), 'Fresh child unexpectedly starts in explicit EcoQoS'
    requested = original[0] | 1, original[1] & ~1
    checks = []
    assert prefs.preferPerformanceDuringImport is True
    assert prefs.bl_rna.properties['preferPerformanceDuringImport'].default is True

    def operator_case(status, expected, error=False):
        observed = []
        def dispatch(self, context):
            observed.append(api.read())
            if error:
                raise RuntimeError('Deliberate production Import failure')
            return {status}
        with patch.object(addon.BLOSM_OT_ImportData, '_execute', dispatch):
            if error:
                try:
                    bpy.ops.blosm.import_data()
                except RuntimeError as exc:
                    assert 'Deliberate production Import failure' in str(exc)
                else:
                    raise AssertionError('Import error was swallowed')
            else:
                assert bpy.ops.blosm.import_data() == {status}
        assert observed == [expected], observed

    try:
        for result in ('FINISHED', 'CANCELLED'):
            operator_case(result, requested)
            assert api.read() == original
        operator_case('FINISHED', requested, error=True)
        assert api.read() == original
        checks.append('real bpy operator: enabled default, success, cancellation and raised failure restore policy')

        prefs.preferPerformanceDuringImport = False
        operator_case('FINISHED', original)
        assert api.read() == original
        prefs.preferPerformanceDuringImport = True
        checks.append('preference disables the request')
        with patch.object(module, '_WindowsProcessPolicy', side_effect=OSError('Unavailable test API')):
            operator_case('FINISHED', original)
        checks.append('unavailable API preserves normal Import result')

        for mode, expected in (('eco', (original[0] | 1, original[1] | 1)), ('high', requested)):
            with diagnostics.ProcessQoS(mode):
                operator_case('FINISHED', expected)
                assert api.read() == expected
            assert api.read() == original
        checks.append('existing explicit EcoQoS and HighQoS are respected')

        props = bpy.context.scene.blosm
        props.progressiveImport = False
        missing = args.output / 'intentionally_missing.osm'
        assert not missing.exists()
        for key, value in dict(mode='3Drealistic', dataType='osm', osmSource='file',
                osmFilepath=str(missing), buildings=False, forests=False, vegetation=False,
                aviation=False, highways=True, railways=False, water=False,
                loadMissingMembers=False, setupScript='', streetBackend='python').items():
            setattr(props, key, value)
        observed = []
        osm_runtime = importlib.import_module(PACKAGE + '.app.osm_import')
        preflight = osm_runtime.OsmImport.check_prerequisites
        def preflight_failure(self, *a, **kw):
            observed.append(api.read())
            return preflight(self, *a, **kw)
        with patch.object(osm_runtime.OsmImport, 'check_prerequisites', preflight_failure):
            try:
                bpy.ops.blosm.import_data()
            except RuntimeError as exc:
                assert 'No such file' in str(exc)
            else:
                raise AssertionError('Missing OSM did not fail')
        assert observed == [requested] and api.read() == original
        checks.append('real missing-file OSM Import restores policy after preflight failure')

        observed = []
        def overlay_step(context):
            observed.append(api.read())
            return {'RUNNING_MODAL'}
        overlay = SimpleNamespace(importNextTile=overlay_step)
        modal = addon.BLOSM_OT_ControlOverlay.modal
        assert modal(overlay, bpy.context, SimpleNamespace(type='TIMER')) == {'RUNNING_MODAL'}
        assert observed == [requested] and api.read() == original
        observed.clear()
        modal(overlay, bpy.context, SimpleNamespace(type='MOUSEMOVE'))
        assert observed == [] and api.read() == original
        def overlay_error(context):
            assert api.read() == requested
            raise RuntimeError('Deliberate overlay step failure')
        overlay.importNextTile = overlay_error
        try:
            modal(overlay, bpy.context, SimpleNamespace(type='TIMER'))
        except RuntimeError:
            pass
        else:
            raise AssertionError('Overlay error was swallowed')
        assert api.read() == original
        checks.append('overlay timer processing scopes request; non-timer events and failed steps restore normally')

        for fail in (False, True):
            observed = []
            job = SimpleNamespace(done=False, error=None)
            def step():
                observed.append(('step', api.read()))
                if fail:
                    raise RuntimeError('Deliberate progressive step failure')
            job.step = step
            job.poll = lambda: observed.append(('poll', api.read()))
            adapter = SimpleNamespace(scene=bpy.context.scene, view_layer=bpy.context.view_layer,
                                      revision=0, status='Test', job=job)
            def cancel():
                observed.append(('cancel', api.read()))
                job.done = True
            adapter.cancel = cancel
            adapter.finish = lambda silent=False: observed.append(('finish', api.read()))
            controller = runtime.start(adapter, bpy.context)
            runtime._tick()
            assert api.read() == original
            if not fail:
                assert runtime.active is controller
                adapter.view_layer = None  # Paused: poll only, no performance scope.
                controller.tick()
                assert ('poll', original) in observed and api.read() == original
                runtime.stop_for_scene_change()
            assert runtime.active is None and api.read() == original
            assert all(policy == requested for label, policy in observed if label != 'poll'), observed
            assert any(label == 'step' for label, _ in observed)
            assert any(label == 'finish' for label, _ in observed)
        checks.append('progressive steps, pauses, cancellation/finalization and failed steps restore between callbacks')

        assert all(state.current()[key] == priority[key] for key in ('process_priority_class', 'thread_priority'))
        report = dict(passed=True, blender=bpy.app.version_string, checkout=str(ROOT),
                      package=PACKAGE, original=original, requested=requested, restored=api.read(), checks=checks)
        (args.output / 'checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
        print('PRODUCTION_IMPORT_PERFORMANCE_OK', json.dumps(report), flush=True)
    finally:
        runtime.stop_for_scene_change()
        api.write(original)
        addon.unregister()


if __name__ == '__main__':
    main()
