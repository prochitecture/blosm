"""Repeat monitored, factory-startup imports in a chat-specific Blender profile.

Run with host Python. The existing bounded runner supplies safety caps but no
working-set throttle is used: deliberately paged runs are not a timing baseline.
"""
import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
DIRECTORIES = {'BLENDER_USER_RESOURCES': 'resources', 'BLENDER_USER_CONFIG': 'config',
    'BLENDER_USER_SCRIPTS': 'scripts', 'BLENDER_USER_DATAFILES': 'datafiles',
    'BLENDER_USER_EXTENSIONS': 'extensions', 'TEMP': 'temp', 'TMP': 'temp'}


def summarize(output):
    runs = []
    for path in sorted(output.glob('run_*/timings.json')):
        result = json.loads(path.read_text())
        # Older unprofiled receipts predate the selectable profiler scope.
        if not result['settings'].get('diagnostic_profile'):
            result['settings'].setdefault('profile_scope', None)
        result['settings'].setdefault('construction_diagnostics', False)
        resources = json.loads((path.parent / 'process/resources.json').read_text())
        if not result['passed'] or not resources['completed']:
            raise RuntimeError('Incomplete run: ' + str(path))
        if resources['working_set_requested_bytes'] is not None:
            raise RuntimeError('A throttled run cannot enter this baseline')
        runs.append((result, resources))
    if not runs:
        raise RuntimeError('No completed runs in ' + str(output))
    assert len({r['input_sha256'] for r, _ in runs}) == 1
    assert len({r.get('source_scene_sha256') for r, _ in runs}) == 1
    assert len({r['revision'] for r, _ in runs}) == 1
    assert len({r.get('runtime_source_sha256') for r, _ in runs}) == 1
    assert len({m.get('affinity_cpu_requested') for _, m in runs}) == 1
    assert all(r['settings'] == runs[0][0]['settings'] for r, _ in runs)
    assert len({(r['meshes'], r['quads'], r['vertices']) for r, _ in runs}) == 1
    timing_keys = runs[0][0]['timings']
    report = dict(passed=True, run_count=len(runs), revision=runs[0][0]['revision'],
        input=runs[0][0]['input'], input_sha256=runs[0][0]['input_sha256'],
        source_scene_sha256=runs[0][0].get('source_scene_sha256'),
        runtime_source_sha256=runs[0][0].get('runtime_source_sha256'),
        blender=runs[0][0]['blender'], settings=runs[0][0]['settings'],
        meshes=runs[0][0]['meshes'], quads=runs[0][0]['quads'],
        generation_seconds=[r['timings']['street_generation_seconds'] for r, _ in runs],
        median_timings={key: statistics.median(r['timings'][key] for r, _ in runs) for key in timing_keys},
        median_task_seconds={key: statistics.median(r['task_seconds'].get(key, 0) for r, _ in runs)
                             for key in runs[0][0]['task_seconds']},
        peak_resident_bytes=max(r['peak_resident'] for _, r in runs),
        peak_sampled_private_bytes=max(r['peak_private'] for _, r in runs),
        exclusions=runs[0][0]['exclusions'],
        affinity_cpu_requested=runs[0][1].get('affinity_cpu_requested'),
        conditions=('Profiler-instrumented diagnostic; not an unprofiled timing baseline.'
            if runs[0][0]['settings'].get('diagnostic_profile') else
            'Construction attempt ledger; diagnostic timings include overhead.'
            if runs[0][0]['settings'].get('construction_diagnostics') else
            'Coarse wrappers and calibration kernels; diagnostic timings include overhead.'
            if runs[0][0]['settings'].get('execution_diagnostics') else
            'Fresh Blender processes, OS file cache not flushed, no working-set throttle; CPU affinity, if requested, is reported separately.'))
    (output / 'summary.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(report, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--osm', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--profile-root', type=Path, required=True)
    parser.add_argument('--addon-root', type=Path, help='Read-only addon snapshot for reference comparisons')
    parser.add_argument('--addon-revision', help='Commit archived into --addon-root; its runtime digest is also recorded')
    parser.add_argument('--blender', type=Path, default=Path('C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'))
    parser.add_argument('--guard', type=Path, default=ROOT.parent / 'prochitecture_gn/experiments/streets/blosm_integration/bounded_blender.py')
    parser.add_argument('--runs', type=int, default=3)
    parser.add_argument('--threads', type=int, default=min(64, os.cpu_count() or 1))
    parser.add_argument('--resolution', type=int, choices=(1, 2, 3, 4), default=1)
    parser.add_argument('--scene', type=Path, help='Use effective saved street settings and terrain instead of benchmark defaults')
    parser.add_argument('--without-terrain', action='store_true', help='Clear the saved terrain selection for a control run')
    parser.add_argument('--numeric-threads', type=int, choices=range(1, 65), default=1)
    parser.add_argument('--diagnostic-profile', action='store_true', help='Collect call profiles; do not compare these timings with unprofiled runs')
    parser.add_argument('--profile-scope', choices=('numerical', 'publication', 'all'), default='numerical',
                        help='With --diagnostic-profile, collect separate numerical/publication call profiles')
    parser.add_argument('--operator', action='store_true', help='Use the normal streets-only Import in 3D realistic mode')
    parser.add_argument('--execution-diagnostics', action='store_true', help='Record coarse calls and CPU/calibration evidence')
    parser.add_argument('--construction-diagnostics', action='store_true', help='Record chart attempts and bounded replay inputs')
    parser.add_argument('--affinity-cpu', type=int, help='Diagnostic: bind only the test child to this CPU using the resource guard')
    parser.add_argument('--qos', choices=('system', 'high', 'eco'), default='system', help='Explicit current-test-process scheduling policy; restored after timing')
    parser.add_argument('--import-performance', action='store_true', help='Enable the production performance preference for normal Import')
    parser.add_argument('--memory-gib', type=float, default=3.)
    parser.add_argument('--reserve-mib', type=int, default=1024)
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    if args.summarize_only:
        return summarize(args.output)
    if bool(args.addon_root) != bool(args.addon_revision):
        parser.error('--addon-root and --addon-revision must be supplied together')
    if args.without_terrain and not args.scene:
        parser.error('--without-terrain requires --scene')
    if args.operator and not args.scene:
        parser.error('--operator requires --scene')
    if args.import_performance and (not args.operator or args.qos != 'system'):
        parser.error('--import-performance requires --operator and no external --qos override')
    if args.execution_diagnostics and args.diagnostic_profile:
        parser.error('Keep coarse execution diagnostics separate from cProfile')
    if args.profile_scope != 'numerical' and not args.diagnostic_profile:
        parser.error('--profile-scope requires --diagnostic-profile')
    if args.affinity_cpu is not None and (args.threads != 1 or args.numeric_threads != 1):
        parser.error('--affinity-cpu requires --threads 1 and --numeric-threads 1')
    if not 1 <= args.runs <= 9:
        parser.error('Run count must be between 1 and 9')
    env = os.environ.copy()
    for key, folder in DIRECTORIES.items():
        path = args.profile_root.resolve() / folder
        path.mkdir(parents=True, exist_ok=True)
        env[key] = str(path)
    env.pop('BLOSM_TEST_WORKING_SET_MIB', None)
    revision = args.addon_revision or subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    for number in range(1, args.runs + 1):
        output = args.output / f'run_{number:02d}'
        if output.exists():
            parser.error('Use a fresh output directory: ' + str(output))
        command = [sys.executable, str(args.guard), '--blender', str(args.blender),
            '--script', str(ROOT / 'tests/benchmark_python_street_import.py'), '--output', str(output / 'process'),
            '--memory-gib', str(args.memory_gib), '--reserve-mib', str(args.reserve_mib),
            '--threads', str(args.threads), '--numeric-threads', str(args.numeric_threads)]
        if args.affinity_cpu is not None:
            command.extend(('--affinity-cpu', str(args.affinity_cpu)))
        command.extend(('--', '--osm', str(args.osm.resolve()), '--output', str(output),
                        '--revision', revision, '--resolution', str(args.resolution)))
        if args.addon_root:
            command.extend(('--addon-root', str(args.addon_root.resolve())))
        if args.scene:
            command.extend(('--scene', str(args.scene.resolve())))
        if args.without_terrain:
            command.append('--without-terrain')
        if args.diagnostic_profile:
            command.append('--diagnostic-profile')
            command.extend(('--profile-scope', args.profile_scope))
        if args.operator:
            command.append('--operator')
        if args.import_performance:
            command.append('--import-performance')
        if args.execution_diagnostics:
            command.append('--execution-diagnostics')
        if args.construction_diagnostics:
            command.append('--construction-diagnostics')
        if args.qos != 'system':
            command.extend(('--qos', args.qos))
        if number == args.runs:
            command.append('--save-scene')
        subprocess.run(command, env=env, cwd=ROOT, check=True,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0)
    summarize(args.output)


if __name__ == '__main__':
    main()
