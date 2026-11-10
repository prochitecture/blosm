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
        resources = json.loads((path.parent / 'process/resources.json').read_text())
        if not result['passed'] or not resources['completed']:
            raise RuntimeError('Incomplete run: ' + str(path))
        if resources['working_set_requested_bytes'] is not None:
            raise RuntimeError('A throttled run cannot enter this baseline')
        runs.append((result, resources))
    if not runs:
        raise RuntimeError('No completed runs in ' + str(output))
    assert len({r['input_sha256'] for r, _ in runs}) == 1
    assert len({r['revision'] for r, _ in runs}) == 1
    assert all(r['settings'] == runs[0][0]['settings'] for r, _ in runs)
    assert len({(r['meshes'], r['quads'], r['vertices']) for r, _ in runs}) == 1
    timing_keys = runs[0][0]['timings']
    report = dict(passed=True, run_count=len(runs), revision=runs[0][0]['revision'],
        input=runs[0][0]['input'], input_sha256=runs[0][0]['input_sha256'],
        blender=runs[0][0]['blender'], settings=runs[0][0]['settings'],
        meshes=runs[0][0]['meshes'], quads=runs[0][0]['quads'],
        generation_seconds=[r['timings']['street_generation_seconds'] for r, _ in runs],
        median_timings={key: statistics.median(r['timings'][key] for r, _ in runs) for key in timing_keys},
        median_task_seconds={key: statistics.median(r['task_seconds'].get(key, 0) for r, _ in runs)
                             for key in runs[0][0]['task_seconds']},
        peak_resident_bytes=max(r['peak_resident'] for _, r in runs),
        peak_sampled_private_bytes=max(r['peak_private'] for _, r in runs),
        exclusions=runs[0][0]['exclusions'],
        conditions='Fresh Blender processes, OS file cache not flushed, no CPU affinity or working-set throttle; monitored local wall-clock baseline.')
    (output / 'summary.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(report, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--osm', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--profile-root', type=Path, required=True)
    parser.add_argument('--blender', type=Path, default=Path('C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'))
    parser.add_argument('--guard', type=Path, default=ROOT.parent / 'prochitecture_gn/experiments/streets/blosm_integration/bounded_blender.py')
    parser.add_argument('--runs', type=int, default=3)
    parser.add_argument('--threads', type=int, default=min(64, os.cpu_count() or 1))
    parser.add_argument('--resolution', type=int, choices=(1, 2, 3, 4), default=1)
    parser.add_argument('--memory-gib', type=float, default=3.)
    parser.add_argument('--reserve-mib', type=int, default=1024)
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    if args.summarize_only:
        return summarize(args.output)
    if not 1 <= args.runs <= 9:
        parser.error('Run count must be between 1 and 9')
    env = os.environ.copy()
    for key, folder in DIRECTORIES.items():
        path = args.profile_root.resolve() / folder
        path.mkdir(parents=True, exist_ok=True)
        env[key] = str(path)
    env.pop('BLOSM_TEST_WORKING_SET_MIB', None)
    revision = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    for number in range(1, args.runs + 1):
        output = args.output / f'run_{number:02d}'
        if output.exists():
            parser.error('Use a fresh output directory: ' + str(output))
        command = [sys.executable, str(args.guard), '--blender', str(args.blender),
            '--script', str(ROOT / 'tests/benchmark_python_street_import.py'), '--output', str(output / 'process'),
            '--memory-gib', str(args.memory_gib), '--reserve-mib', str(args.reserve_mib),
            '--threads', str(args.threads), '--numeric-threads', '1', '--', '--osm', str(args.osm.resolve()),
            '--output', str(output), '--revision', revision, '--resolution', str(args.resolution)]
        if number == args.runs:
            command.append('--save-scene')
        subprocess.run(command, env=env, cwd=ROOT, check=True,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0)
    summarize(args.output)


if __name__ == '__main__':
    main()
