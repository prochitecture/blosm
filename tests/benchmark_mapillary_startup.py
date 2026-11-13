"""Cached real coverage + synthetic metadata through Blender's UI controller.

Launch only in an isolated factory-startup Blender UI process. Pass --artifacts
under workspace/artifacts; the script redirects Blender's recovery files there.
--revision 843c77c compares the committed reader/adapter with current code without
checking out files. Network time is deliberately excluded; optional fake response
delays are useful for checking whether metadata overlaps coverage processing.
--count-interval is a benchmark-only timer experiment during initial counting;
it does not change any production preference or installed Blender configuration.
"""
import argparse
import ast
from collections import deque
import importlib
import json
from pathlib import Path
import subprocess
import sys
import time
import types
from urllib.parse import parse_qs, urlsplit

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))

parser = argparse.ArgumentParser()
parser.add_argument('--artifacts', type=Path, required=True)
parser.add_argument('--fixture', type=Path,
                    default=Path('C:/dev/workspace/artifacts/mapillary_import/count_optimization'))
parser.add_argument('--revision', help='Read committed adapter/decoder instead of the working tree')
parser.add_argument('--complete', action='store_true', help='Continue through all cameras')
parser.add_argument('--workers', type=int, default=4)
parser.add_argument('--coverage-delay', type=float, default=0.)
parser.add_argument('--metadata-delay', type=float, default=0.)
parser.add_argument('--count-interval', type=float,
                    help='Experimental first-count timer interval, .001–.03 seconds')
parser.add_argument('--timeout', type=float, default=240.)
args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
if not 1 <= args.workers <= 8 or min(args.coverage_delay, args.metadata_delay) < 0 or args.timeout <= 0:
    parser.error('Use 1–8 workers, nonnegative fake delays and a positive timeout')
if args.count_interval is not None and not .001 <= args.count_interval <= .03:
    parser.error('--count-interval must be between .001 and .03 seconds')
args.artifacts.mkdir(parents=True, exist_ok=True)

import bpy
bpy.context.preferences.filepaths.temporary_directory = str(args.artifacts.resolve())
if bpy.app.background or bpy.context.window is None:
    raise RuntimeError('This benchmark requires an isolated Blender UI process')

addon = importlib.import_module(ROOT.name)
addon.register()
adapter = importlib.import_module(ROOT.name + '.mapillary.progressive')
runtime = importlib.import_module(ROOT.name + '.progressive.blender')
Projection = importlib.import_module(ROOT.name + '.util.transverse_mercator').TransverseMercator
runtime.show_result = lambda *_args: None


def revision_module(path, name, bindings=None):
    """Relative imports use this checkout; replace changed baseline dependencies below."""
    source = subprocess.check_output(['git', 'show', args.revision + ':' + path], cwd=ROOT)
    module = types.ModuleType(name)
    module.__package__ = ROOT.name + '.mapillary'
    module.__file__ = str(ROOT / path)
    tree = ast.parse(source, filename=module.__file__)
    if bindings is not None:
        # The current decoder may have renamed APIs. Do not import those when
        # reproducing a historical adapter; bind its historical dependencies.
        tree.body = [node for node in tree.body
                     if not (isinstance(node, ast.ImportFrom) and node.level == 1
                             and node.module in ('tiles', 'progress'))]
        module.__dict__.update(bindings)
    exec(compile(tree, module.__file__, 'exec'), module.__dict__)
    return module


if args.revision:
    old_tiles = revision_module('mapillary/tiles.py', 'startup_baseline_tiles')
    old_progress = revision_module('mapillary/progress.py', 'startup_baseline_progress')
    bindings = {'ImageProgress': old_progress.ImageProgress}
    for name in ('MAX_TILE_BYTES', 'contains', 'identifier', 'image_points',
                 'image_ids', 'tile_coordinates', 'validate_bounds'):
        bindings[name] = getattr(old_tiles, name)
    adapter = revision_module('mapillary/progressive.py', 'startup_baseline_adapter', bindings)

fixture = json.loads((args.fixture / 'baseline.json').read_text())
payload = (args.fixture / 'coverage.mvt').read_bytes()
bounds = fixture['bounds']
origin = ((bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2)
started = previous = None
events = {}
gaps = []
steps = []


def mark(name):
    events.setdefault(name, time.monotonic() - started)


class Transport:
    """Responses retain real request IDs; Job still enforces the shared window."""
    def __init__(self):
        self.max_requests = args.workers
        self.replies = deque()
        self.identity = 0
        self.coverage_requests = self.metadata_requests = 0

    def submit(self, request):
        self.identity += 1
        query = parse_qs(urlsplit(request.url).query)
        if 'image_ids' in query:
            mark('first_metadata_request')
            self.metadata_requests += 1
            entries = [{'id': identity, 'camera_type': 'perspective',
                        'width': 1920, 'height': 1080, 'exif_orientation': 1,
                        'camera_parameters': [.8, 0, 0], 'computed_compass_angle': 0}
                       for identity in query['image_ids'][0].split(',')]
            content = json.dumps({'data': entries}).encode()
            delay = args.metadata_delay
        else:
            self.coverage_requests += 1
            content, delay = payload, args.coverage_delay
        self.replies.append((time.monotonic() + delay, self.identity, content, None))
        return self.identity

    def poll(self):
        # Completion order, allowing a ready metadata response to bypass a slow
        # coverage reply when exercising a multi-tile fixture in the future.
        now = time.monotonic()
        for index, reply in enumerate(self.replies):
            if reply[0] <= now:
                del self.replies[index]
                return reply[1:]
        return None

    def close(self):
        self.replies.clear()


transport = Transport()
task = adapter.MapillaryImport(bpy.context, bounds, Projection(lat=origin[0], lon=origin[1]),
                               origin, False, 'offline-fixture-token', transport=transport)
original_step, original_commit = task.job.step, task.commit


def step(*args, **kwargs):
    turn = time.monotonic()
    try:
        return original_step(*args, **kwargs)
    finally:
        steps.append(time.monotonic() - turn)
        if task.progress.percent(False) is not None:
            mark('first_estimate')


def commit(*args, **kwargs):
    result = original_commit(*args, **kwargs)
    mark('first_camera')
    return result


task.job.step, task.commit = step, commit


def observe():
    global previous
    now = time.monotonic()
    if previous is not None:
        gaps.append(now - previous)
    previous = now
    timeout = now - started > args.timeout
    first_camera = task.count > 0 and not args.complete
    if not (task.job.done or timeout or first_camera):
        return .03

    error = str(task.job.error) if task.job.error else None
    complete = task.job.done and not task.job.cancelled and error is None
    if timeout:
        error = error or 'Benchmark deadline exceeded'
    if not task.count:
        error = error or 'No camera was created'
    created_ids = sorted(obj['mapillary_image_id'] for obj in task.collection.objects) if task.collection else []
    expected_ids = sorted(fixture['ids'])
    if complete and created_ids != expected_ids:
        error = error or 'Completed camera IDs differ from cached fixture baseline'
    report = {'status': 'failed' if error else 'passed', 'error': error,
              'revision': args.revision or 'working-tree', 'events_seconds': events,
              'elapsed_seconds': now - started, 'completed': complete,
              'cameras_created': task.count, 'created_ids': created_ids,
              'candidates': task.progress.candidates, 'expected_candidates': len(expected_ids),
              'metrics': task.metrics, 'max_job_turn': max(steps, default=0),
              'heartbeats': len(gaps), 'max_heartbeat_gap': max(gaps, default=0),
              'coverage_requests': transport.coverage_requests,
              'metadata_requests': transport.metadata_requests, 'errors': dict(task.errors),
              'import_status': task.status, 'blender': bpy.app.version_string,
              'global_undo': bpy.context.preferences.edit.use_global_undo,
              'experimental_count_interval': args.count_interval,
              'source': 'cached real coverage; synthetic metadata; no network or terrain'}
    if runtime.active is not None:
        runtime.active.close(silent=True)
    (args.artifacts / 'startup.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({key: value for key, value in report.items() if key != 'created_ids'}), flush=True)
    bpy.ops.wm.quit_blender()
    return None


if args.count_interval is not None:
    original_tick = runtime._tick

    def counting_tick():
        interval = original_tick()
        # Preserve shutdown and the shared default interval once a tile count
        # exists or cameras appear. Assign before start() so registration and
        # Controller.close() refer to this same callable when removing timers.
        if interval is not None and task.progress.counted_tiles == 0 and task.count == 0:
            return args.count_interval
        return interval

    runtime._tick = counting_tick

started = time.monotonic()
runtime.start(task, bpy.context)
bpy.app.timers.register(observe, first_interval=.05)
