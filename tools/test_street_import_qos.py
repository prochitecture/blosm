"""Run serial, fresh-window system/HighQoS Import comparisons behind an idle Blender.

Use this chat's isolated profile and a fresh artifacts directory. Computer Use
may activate the returned idle test window; importing windows wait until its
PID is foreground. No normal user Blender window is used or changed.
"""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DIRECTORIES = {'BLENDER_USER_RESOURCES': 'resources', 'BLENDER_USER_CONFIG': 'config',
    'BLENDER_USER_SCRIPTS': 'scripts', 'BLENDER_USER_EXTENSIONS': 'extensions',
    'BLENDER_USER_DATAFILES': 'datafiles', 'TEMP': 'temp', 'TMP': 'temp'}


def save(path, data):
    temporary = path.with_suffix('.new')
    temporary.write_text(json.dumps(data, indent=2), encoding='utf8')
    for attempt in range(6):
        try:
            temporary.replace(path)
            return
        except PermissionError:
            if attempt == 5:
                raise
            time.sleep(.05)


def environment(profile):
    env = os.environ.copy()
    for key, folder in DIRECTORIES.items():
        path = profile.resolve() / folder
        path.mkdir(parents=True, exist_ok=True)
        env[key] = str(path)
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        env[key] = '1'
    env.pop('BLOSM_TEST_WORKING_SET_MIB', None)
    return env


def summarize(output):
    state = json.loads((output / 'progress.json').read_text())
    assert state['status'] == 'COMPLETE', state.get('error', state['status'])
    assert len(state['cases']) == 4 and all(
        case.get('verified_background') and case.get('restored_policy')
        and case.get('exit_code') == 0 for case in state['cases']), state['cases']
    runs = [json.loads(p.read_text()) for p in sorted(output.glob('case_*/foreground_test_001.json'))]
    assert len(runs) == 4 and all(r['passed'] for r in runs)
    assert [r['requested_qos'] for r in runs] == ['system', 'high', 'high', 'system']
    assert len({json.dumps(r['settings'], sort_keys=True) for r in runs}) == 1
    grouped = {mode: [r['operator_seconds'] for r in runs if r['requested_qos'] == mode]
               for mode in ('system', 'high')}
    assert all(len(values) == 2 for values in grouped.values())
    medians = {mode: statistics.median(values) for mode, values in grouped.items()}
    report = dict(passed=True, order=[r['requested_qos'] for r in runs], seconds=grouped,
        median_seconds=medians, reduction_percent=100.*(1.-medians['high']/medians['system']),
        qualification='Fresh GUI process per run, identical light observations, one idle foreground target; saved-mesh equivalence is a separate gate.')
    save(output / 'summary.json', report)
    print(json.dumps(report, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('scene', 'osm', 'output', 'profile-root'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--blender', type=Path, default=Path('C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'))
    parser.add_argument('--guard', type=Path, default=ROOT.parent / 'prochitecture_gn/experiments/streets/blosm_integration/bounded_blender.py')
    parser.add_argument('--case-timeout', type=float, default=900.)
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    if args.summarize_only:
        return summarize(args.output)
    args.output.mkdir(parents=True, exist_ok=False)
    revision = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    spec = importlib.util.spec_from_file_location('qos_resource_guard', args.guard)
    guard = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = guard
    spec.loader.exec_module(guard)
    monitor = guard.Monitor()
    cap, reserve = 3*2**30, 1024*2**20
    common = [str(args.blender), '--factory-startup', '--disable-autoexec', '--threads', '16', '--python-exit-code', '1']
    state = dict(status='STARTING', checkout=str(ROOT), revision=revision,
        scene=str(args.scene.resolve()), osm=str(args.osm.resolve()),
        scene_sha256=hashlib.sha256(args.scene.read_bytes()).hexdigest(),
        osm_sha256=hashlib.sha256(args.osm.read_bytes()).hexdigest(),
        memory_cap_bytes=cap, reserve_bytes=reserve, cases=[])
    progress = args.output / 'progress.json'
    save(progress, state)
    children = []

    def launch(command, folder, profile):
        if monitor.available() < reserve:
            raise RuntimeError('System memory reserve is unavailable before launch')
        folder.mkdir(parents=True, exist_ok=True)
        with (folder / 'stdout.log').open('wb') as out, (folder / 'stderr.log').open('wb') as err:
            proc = subprocess.Popen(command, cwd=ROOT, env=environment(profile), stdout=out, stderr=err,
                                    creationflags=subprocess.CREATE_NEW_PROCESS_GROUP)
        children.append(proc)
        return proc

    def wait_ready(proc, folder):
        started = time.monotonic()
        while not (folder / 'ready.json').exists():
            if proc.poll() is not None:
                raise RuntimeError(f'GUI startup exited with {proc.returncode}: {folder}')
            if time.monotonic() - started > 60.:
                raise RuntimeError('GUI startup timed out: ' + str(folder))
            time.sleep(.25)

    try:
        target_folder = args.output / 'target'
        target = launch(common + ['--python', str(ROOT / 'tests/street_performance_focus_target.py'), '--',
            '--output', str(target_folder), '--result', str(args.output / 'session_done.json')],
            target_folder, args.profile_root / 'qos_focus_target')
        wait_ready(target, target_folder)
        state['target_pid'] = target.pid
        save(progress, state)
        for number, mode in enumerate(('system', 'high', 'high', 'system'), 1):
            folder = args.output / f'case_{number:02d}_{mode}'
            command = common + ['--python', str(ROOT / 'tests/open_street_performance_scene.py'), '--',
                '--scene', str(args.scene.resolve()), '--osm', str(args.osm.resolve()),
                '--output', str(folder), '--revision', revision, '--auto-run-after', '1',
                '--background-target-pid', str(target.pid), '--qos', mode, '--save-result', '--quit-after-run']
            proc = launch(command, folder, args.profile_root)
            case = dict(number=number, qos=mode, pid=proc.pid, folder=str(folder),
                        peak_resident=0, peak_private=0, working_set_throttle=False)
            state['cases'].append(case)
            wait_ready(proc, folder)
            state.update(status='WAITING_FOR_FOREGROUND_TARGET', current_case=number, current_pid=proc.pid)
            save(progress, state)
            started, low = time.monotonic(), 0
            while proc.poll() is None:
                if target.poll() is not None:
                    raise RuntimeError('The owned foreground target exited during the comparison')
                if time.monotonic() - started > args.case_timeout:
                    raise RuntimeError('Case timed out: ' + str(folder))
                sample = monitor.sample(proc.pid)
                case['last_sample'] = sample
                guard.record_memory_sample(case, sample)
                low = low + 1 if sample['system_available'] < reserve else 0
                if guard.memory_limit_metric(sample, cap) or low >= 3:
                    raise RuntimeError('Resource limit reached: ' + str(folder))
                if (folder / 'started.json').exists():
                    state['status'] = 'IMPORT_RUNNING'
                save(progress, state)
                time.sleep(.5)
            case['exit_code'] = proc.returncode
            guard.record_memory_sample(case, monitor.sample_handle(int(proc._handle)))
            assert proc.returncode == 0, case
            result = json.loads((folder / 'foreground_test_001.json').read_text())
            assert result['passed'], result.get('error')
            before, after = result['policy_before']['process_power_policy'], result['policy_after']['process_power_policy']
            assert before == after and before['success'], (before, after)
            assert all(result['priority_before'][k] == result['priority_after'][k]
                       for k in ('process_priority_class', 'thread_priority'))
            checks = result['processor_foreground_checks']
            assert checks and all(not ast.literal_eval(k)[0] and ast.literal_eval(k)[1] == target.pid for k in checks), checks
            if mode == 'high':
                policy = result['execution_policy']
                assert policy['before'] == policy['after_restore']
                assert policy['during']['control_mask'] == before['control_mask'] | 1
                assert policy['during']['state_mask'] == before['state_mask'] & ~1
            else:
                assert result['execution_policy'] is None
            case.update(seconds=result['operator_seconds'], verified_background=True, restored_policy=True)
            save(progress, state)
        state['status'] = 'COMPLETE'
        save(progress, state)
        summarize(args.output)
    except BaseException as error:
        state.update(status='FAILED', error=repr(error))
        save(progress, state)
        raise
    finally:
        save(args.output / 'session_done.json', dict(status=state['status']))
        for proc in children:
            if proc.poll() is None:
                try:
                    proc.wait(timeout=3.)
                except subprocess.TimeoutExpired:
                    proc.terminate()
                    proc.wait()


if __name__ == '__main__':
    main()
