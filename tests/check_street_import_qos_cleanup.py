"""Verify HighQoS restoration in Blender, including a real failed Import."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import bpy

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location('blosm', ROOT / '__init__.py',
                                                submodule_search_locations=[str(ROOT)])
    addon = importlib.util.module_from_spec(spec)
    sys.modules['blosm'] = addon
    spec.loader.exec_module(addon)
    addon.register()
    from blosm.tests.street_execution_diagnostics import ProcessorState, ProcessQoS
    state = ProcessorState()
    original = state.power()['process_power_policy']
    priority = state.current()
    checks = {}

    class ExpectedFailure(RuntimeError):
        pass

    try:
        with ProcessQoS('high') as failure:
            assert state.power()['process_power_policy']['state_mask'] & 1 == 0
            raise ExpectedFailure('Deliberate exception inside the performance scope')
    except ExpectedFailure:
        pass
    assert state.power()['process_power_policy'] == original
    failure.close()  # Repeated cleanup is harmless.
    checks['exception_and_repeated_cleanup'] = failure.report

    with ProcessQoS('eco') as outer:
        eco = state.power()['process_power_policy']
        assert eco['control_mask'] & 1 and eco['state_mask'] & 1
        with ProcessQoS('high') as inner:
            assert state.power()['process_power_policy']['state_mask'] & 1 == 0
        assert state.power()['process_power_policy'] == eco
        checks['nested_high_from_existing_eco'] = inner.report
    assert state.power()['process_power_policy'] == original
    checks['nested_outer_restored'] = outer.report

    entry = bpy.context.preferences.addons.new()
    entry.module = 'blosm'
    entry.preferences.dataDir = str(args.output)
    entry.preferences.enableExperimentalFeatures = True
    props = bpy.context.scene.blosm
    missing = args.output / 'intentionally_missing.osm'
    assert not missing.exists()
    for name, value in dict(mode='3Drealistic', dataType='osm', osmSource='file',
            osmFilepath=str(missing), buildings=False, forests=False, vegetation=False,
            aviation=False, highways=True, railways=False, water=False,
            loadMissingMembers=False, setupScript='', streetBackend='python').items():
        setattr(props, name, value)
    status, error = None, None
    try:
        with ProcessQoS('high') as failed_import:
            status = sorted(bpy.ops.blosm.import_data())
    except Exception as exc:
        error = type(exc).__name__ + ': ' + str(exc)
    assert error or status != ['FINISHED'], (status, error)
    assert state.power()['process_power_policy'] == original
    checks['actual_failed_import'] = dict(status=status, error=error, policy=failed_import.report)
    after = state.current()
    assert all(after[k] == priority[k] for k in ('process_priority_class', 'thread_priority'))
    report = dict(passed=True, blender=bpy.app.version_string, checkout=str(ROOT),
                  before=original, after=state.power()['process_power_policy'], checks=checks)
    (args.output / 'cleanup.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print('IMPORT_QOS_CLEANUP_OK', json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
