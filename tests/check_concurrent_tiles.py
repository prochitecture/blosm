"""Blender integration and timing for serial/concurrent, blocking/progressive tiles.

Uses an isolated local HTTP fixture, including controlled network latency. Helsinki
payloads are optional benchmark inputs; results are not live-provider timings.
"""
import argparse
import base64
from collections import Counter
import importlib
import json
from math import radians
from pathlib import Path
import subprocess
import sys
import time
from urllib.request import urlopen

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gltf_test_support import bpy, REPO, PACKAGE, fixture, make_glb, make_b3dm, reset_scene, no_tile_writes_or_operator
import addon_utils

parser = argparse.ArgumentParser()
parser.add_argument('--artifacts', type=Path, required=True)
parser.add_argument('--helsinki', type=Path)
args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
args.artifacts.mkdir(parents=True, exist_ok=True)
package = PACKAGE
addon_utils.enable(package, default_set=True)
m = importlib.import_module(package + '.threed_tiles.manager')
r = importlib.import_module(package + '.threed_tiles.blender')
p = importlib.import_module(package + '.threed_tiles.progressive')
t = importlib.import_module(package + '.downloads.transport')

bounds = (24.93, 60.16, 24.95, 60.18)
volume = {'region': [radians(24.93), radians(60.16), radians(24.95), radians(60.18), -100, 100]}
center = m.BaseManager.fromGeographic(60.17, 24.94, 0)
doc, binary, glb = fixture()
b3dm = make_b3dm(glb, list(center))
doc['nodes'][0]['translation'] = [center.x, center.z, -center.y]
glb = make_glb(doc, binary)


def transform(x):
    return [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, -3, 2, 1]


def inputs(helsinki=False, failures=False):
    payloads, delays, root_children = {}, {}, []
    sources = ([(path.name, path.read_bytes()) for path in sorted(args.helsinki.glob('*.b3dm'))]
               if helsinki else [('tile_%d.%s' % (i, 'glb' if i % 2 == 0 else 'b3dm'),
                                  glb if i % 2 == 0 else b3dm) for i in range(4)])
    for branch in ('a', 'b'):
        children = []
        for index, (name, content) in enumerate(sources):
            path = branch + '/' + name
            payloads[path] = content
            delays[path] = .35 if index == 0 and branch == 'a' else .12
            children.append({'geometricError': 0, 'boundingVolume': volume,
                             'transform': transform(index * 10), 'content': {'uri': name}})
        if failures and branch == 'a':
            children.append({'geometricError': 0, 'boundingVolume': volume, 'content': {'uri': 'missing.glb'}})
            children.append({'geometricError': 0, 'boundingVolume': volume, 'content': {'uri': 'invalid.glb'}})
            payloads['a/invalid.glb'] = b'bad glTF'
        manifest = {'root': {'geometricError': 0, 'boundingVolume': volume, 'children': children}}
        payloads[branch + '/tileset.json'] = json.dumps(manifest).encode()
        root_children.append({'geometricError': 0, 'boundingVolume': volume,
                              'transform': transform(100 if branch == 'a' else 200),
                              'content': {'uri': branch + '/tileset.json?session=' + branch}})
    payloads['root.json'] = json.dumps({'root': {'geometricError': 0, 'boundingVolume': volume,
                                                'children': root_children}}).encode()
    return payloads, delays, 2 * len(sources)


def geometry():
    bpy.context.view_layer.update()
    return [tuple(obj.matrix_world @ v.co) for obj in bpy.context.scene.objects if obj.type == 'MESH'
            for v in obj.data.vertices]


def compare(expected, actual, join):
    from mathutils.kdtree import KDTree
    from math import frexp, ldexp
    assert len(expected) == len(actual), (len(expected), len(actual))
    tree = KDTree(len(expected))
    for index, point in enumerate(expected):
        tree.insert(point, index)
    tree.balance()
    magnitude = max(abs(v) for point in expected for v in point)
    tolerance = max(.001, ldexp(4., frexp(magnitude)[1] - 24)) if join else .001
    distance = max(tree.find(point)[2] for point in actual)
    assert distance <= tolerance, (distance, tolerance)
    return distance


def run(progressive, workers, join=False, helsinki=False, failures=False, truncate=False):
    reset_scene()
    payloads, delays, expected_count = inputs(helsinki, failures)
    server = subprocess.Popen([t.python_executable(), '-I', '-B', str(Path(__file__).with_name('progressive_http_fixture.py'))],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
    operation = None
    try:
        truncated = ['root.json', next(key for key in payloads if key.endswith(('.glb', '.b3dm')))] if truncate else []
        server.stdin.write((json.dumps({'payloads': {k: base64.b64encode(v).decode() for k, v in payloads.items()},
                                        'delays': delays, 'truncate_once': truncated}) + '\n').encode())
        server.stdin.close()
        port = int(server.stdout.readline())
        server.stdout.close()
        url = 'http://127.0.0.1:%s/' % port
        renderer = (p.ProgressiveRenderer if progressive else r.BlenderRenderer)('Concurrent test', join, '')
        manager = m.BaseManager(url + 'root.json', renderer)
        manager.centerLat, manager.centerLon = 60.17, 24.94
        manager.setGeometricError(1)
        manager.constantUriQuery = 'key=fixture'
        manager.concurrentDownloads = workers > 0
        manager.downloadWorkers = max(1, workers)
        consume = manager.consumeContent
        order, seconds = [], []
        def timed_consume(uri, *rest):
            start = time.perf_counter()
            consume(uri, *rest)
            seconds.append(time.perf_counter() - start)
            order.append(uri.split(url)[-1])
        manager.consumeContent = timed_consume
        start = time.perf_counter()
        with no_tile_writes_or_operator(no_dataset_files=True):
            if progressive:
                operation = p.TilesImport(manager, bounds, bpy.context, False)
                while not operation.job.done:
                    operation.job.step()
                    assert time.perf_counter() - start < 30
                    time.sleep(.001)
                assert operation.job.error is None, operation.job.error
                result = operation.job.result
            else:
                result = manager.render(*bounds)
        elapsed = time.perf_counter() - start
        assert result[0] == expected_count, result
        assert len(result[1]) == (2 if failures else 0), result
        with urlopen(url + 'stats') as response:
            stats = json.load(response)
        requested = stats['requests']
        for path, count in Counter(requested).items():
            assert count == (2 if path.split('?', 1)[0].lstrip('/') in truncated else 1), 'Unexpected retry count'
        for path in requested:
            if '/a/' in path or path.startswith('/a/'):
                assert 'session=a' in path and 'key=fixture' in path, path
            if '/b/' in path or path.startswith('/b/'):
                assert 'session=b' in path and 'key=fixture' in path, path
        assert stats['peak'] <= max(1, workers), stats
        if workers > 1:
            assert stats['peak'] > 1, stats
            if not helsinki:
                assert 'a/tile_0' not in order[0], order
        else:
            assert stats['peak'] == 1
        return {'mode': 'progressive' if progressive else 'blocking', 'workers': workers,
                'joined': join, 'dataset': 'Helsinki fixtures' if helsinki else 'synthetic',
                'seconds': elapsed, 'tile_processing_seconds': sum(seconds), 'peak_requests': stats['peak'],
                'errors': len(result[1]), 'truncated_resources_recovered': len(truncated), 'import_order': order}, geometry()
    finally:
        if operation is not None and not operation.job.done:
            operation.cancel()
        server.kill()
        server.wait(timeout=5)
        until = time.monotonic() + 5
        while t.reap() and time.monotonic() < until:
            time.sleep(.01)
        assert not t._live and not t.reap(), 'Helper leaked'


reports = []
try:
    for helsinki in ((False, True) if args.helsinki else (False,)):
        for join in (False, True):
            baseline, expected = run(False, 0, join, helsinki)
            reports.append(baseline)
            for progressive, workers in ((False, 1), (False, 4), (True, 1), (True, 4)):
                report, actual = run(progressive, workers, join, helsinki)
                report['max_geometry_error'] = compare(expected, actual, join)
                reports.append(report)
    for progressive in (False, True):
        report, _ = run(progressive, 4, failures=True)
        reports.append(report)
    for helsinki in ((False, True) if args.helsinki else (False,)):
        baseline, expected = run(False, 4, helsinki=helsinki)
        reports.append(baseline)
        for progressive in (False, True):
            report, actual = run(progressive, 4, helsinki=helsinki, truncate=True)
            report['max_geometry_error'] = compare(expected, actual, False)
            reports.append(report)
    report = {'blender': bpy.app.version_string, 'status': 'passed',
              'latency': 'Local HTTP fixture: first tile 350 ms; other tiles 120 ms', 'runs': reports}
    (args.artifacts / 'concurrent_checks.json').write_text(json.dumps(report, indent=2))
    print('CONCURRENT_CHECKS', json.dumps(report), flush=True)
finally:
    addon_utils.disable(package, default_set=True)
