"""Measure native buildings through the real GUI timer, including undo/redraw.

Run in an isolated GUI Blender with release_support arguments. By default stop
after 20 groups of Karl Marx Allee; --extended completes the full extract.
BLOSM_PADL_BENCHMARK_MODE=blocking selects the synchronous comparison.
BLOSM_PADL_BENCHMARK_BACKEND=python selects the Python backend comparison.
"""
from collections import Counter, defaultdict
from pathlib import Path
import functools
import importlib
import json
import os
import sys
import time
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bpy
from release_support import arguments, enable, configure_osm, catalog, PACKAGE

args = arguments()
enable(args)
driver = importlib.import_module(PACKAGE + '.progressive.blender')
runtime = importlib.import_module(PACKAGE + '.app.osm_import')
backends = importlib.import_module(PACKAGE + '.building.padl_integration.backends')
gn = importlib.import_module(PACKAGE + '.building.padl_integration.gn_runtime')
driver.show_result = lambda *a: None
mode = os.environ.get('BLOSM_PADL_BENCHMARK_MODE', 'progressive')
backend = os.environ.get('BLOSM_PADL_BENCHMARK_BACKEND', 'geometry_nodes')
limit = None if args.extended else 20
timings, updates = defaultdict(list), Counter()
commits, hearts = [], []
started = None


def measure(cls, name, label):
    original = getattr(cls, name)
    @functools.wraps(original)
    def measured(*a, **kw):
        start = time.perf_counter()
        try:
            return original(*a, **kw)
        finally:
            timings[label].append(time.perf_counter() - start)
    setattr(cls, name, measured)


for cls, name, label in (
    (backends.BuildingProgramLibrary, 'tree', 'verify_asset'),
    (backends._ProgressiveGraph, 'checked_tree', 'guard'),
    (backends.BuildingGenerationService, 'generate', 'generate'),
    (gn.StatusReader, '__init__', 'reader_create'),
    (gn.StatusReader, 'read', 'status_read'),
    (gn.StatusReader, '__exit__', 'reader_remove'),
    (driver.Controller, 'checkpoint', 'undo_checkpoint'),
    (runtime.OsmJob, 'step', 'job_step'),
):
    measure(cls, name, label)

original_updated = backends._ProgressiveGraph._updated
def updated(self, scene, graph):
    session = getattr(getattr(self.library, 'app', None), 'osm_session', None)
    if not self.busy and not (session and session.in_step):
        for update in graph.updates:
            item = update.id.original
            updates[(item.bl_rna.identifier, bool(update.is_updated_geometry),
                     bool(update.is_updated_transform), bool(update.is_updated_shading))] += 1
    return original_updated(self, scene, graph)
backends._ProgressiveGraph._updated = updated


def observe():
    if driver.active:
        count = driver.active.adapter.committed_units
        if count and (not commits or count != commits[-1]['groups']):
            commits.append(dict(groups=count, seconds=time.perf_counter() - started,
                                verifications=len(timings['verify_asset'])))
            print('PADL_PROGRESS', json.dumps(commits[-1]), flush=True)


def scenario():
    global started
    fixture = next(entry for entry in catalog() if entry['id'] == 'karl_marx_allee')
    settings = configure_osm(args, fixture, '3Dsimple', False)
    settings.progressiveImport = mode == 'progressive'
    settings.buildingBackend = backend
    settings.buildingSeed = 113
    started = time.perf_counter()
    assert bpy.ops.blosm.import_data('EXEC_DEFAULT', True) == {'FINISHED'}
    while driver.active:
        observe()
        if limit and driver.active.adapter.committed_units >= limit:
            bpy.ops.blosm.cancel_progressive()
            break
        yield
    result = runtime.last_result
    assert result['outcome'] == ('cancelled' if limit and mode == 'progressive' else 'success'), result
    assert not result['errors'], result


steps = scenario()
def tick():
    hearts.append(time.perf_counter())
    try:
        if started is not None and time.perf_counter() - started > 600:
            raise TimeoutError('PADL profiling exceeded ten minutes')
        next(steps)
        return .01
    except StopIteration:
        status, error = 'passed', None
    except Exception:
        status, error = 'failed', traceback.format_exc()
    elapsed = time.perf_counter() - started if started else 0.
    if driver.active:
        driver.active.close(silent=True)
    report = dict(status=status, error=error, mode=mode, backend=backend, limit=limit,
                  elapsed_seconds=elapsed, result=runtime.last_result, commits=commits,
                  heartbeats=len(hearts),
                  max_heartbeat_gap_seconds=max((b-a for a,b in zip(hearts, hearts[1:])), default=0),
                  timings={key: dict(calls=len(values), seconds=sum(values), maximum=max(values, default=0), samples=values)
                           for key, values in timings.items()},
                  external_updates=[dict(type=k[0], geometry=k[1], transform=k[2], shading=k[3], count=v)
                                    for k,v in updates.items()])
    intervals = [right['seconds'] - left['seconds'] for left, right in zip(commits, commits[1:])]
    report['publication_intervals'] = dict(
        count=len(intervals), maximum=max(intervals, default=0),
        mean=sum(intervals) / len(intervals) if intervals else 0,
        p95=sorted(intervals)[int(.95 * (len(intervals)-1))] if intervals else 0)
    (args.artifacts/'padl_performance.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('PADL_PERFORMANCE', json.dumps({key: value for key,value in report.items()
                                      if key not in ('timings','external_updates')}), flush=True)
    bpy.ops.wm.quit_blender()


bpy.app.timers.register(tick, first_interval=.5)
