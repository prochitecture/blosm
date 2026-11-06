"""Live-provider matrix in background Blender, without saving downloaded tiles.

Run with --artifacts DIR and --datasets helsinki aarhus ... . --inventory only
traverses actual manifests and reports the selected tile count, without importing.
Reports contain URL paths, hashes, counts, timings and failures, not payload bytes.
"""
import argparse
from collections import Counter
import hashlib
import importlib
import json
from math import frexp, ldexp
from pathlib import Path
import sys
import time
import traceback
from unittest.mock import patch

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gltf_test_support import bpy, REPO, reset_scene, no_tile_writes_or_operator
from live_tile_datasets import DATASETS
import addon_utils
import numpy as np
from mathutils.kdtree import KDTree

parser = argparse.ArgumentParser()
parser.add_argument('--artifacts', type=Path, required=True)
parser.add_argument('--datasets', nargs='+', choices=list(DATASETS), default=list(DATASETS))
parser.add_argument('--inventory', action='store_true')
parser.add_argument('--quick', action='store_true', help='Only blocking/progressive with four workers')
parser.add_argument('--max-tiles', type=int, default=80)
args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
args.artifacts.mkdir(parents=True, exist_ok=True)
bpy.context.preferences.filepaths.temporary_directory = str(args.artifacts)
package = REPO.name
addon_utils.enable(package, default_set=True)
m = importlib.import_module(package + '.threed_tiles.manager')
r = importlib.import_module(package + '.threed_tiles.blender')
p = importlib.import_module(package + '.threed_tiles.progressive')
t = importlib.import_module(package + '.downloads.transport')
core = importlib.import_module(package + '.progressive.core')
safe_url = importlib.import_module(package + '.downloads.diagnostics').redact_url
reports = []


def write_report():
    result = {'blender': bpy.app.version_string, 'live_network': True,
              'source': 'https://github.com/vvoovv/blosm/wiki/Data-sets-with-OGC-3D-Tiles',
              'runs': reports}
    (args.artifacts / 'live_checks.json').write_text(json.dumps(result, indent=2))


class MeasuredTransport(t.DownloadTransport):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.requests = {}
        self.records = []
        self.peak_pending = 0

    def submit(self, request):
        identity = super().submit(request)
        self.requests[identity] = request.url
        self.peak_pending = max(self.peak_pending, len(self.pending))
        return identity

    def poll(self):
        try:
            reply = super().poll()
        except Exception as exc:
            exc.live_probe = {'helper_exit': self.process.poll() if self.process else None,
                              'pending': [safe_url(url) for url in self.requests.values()],
                              'finished': self.records}
            raise
        if reply is not None:
            identity, data, failure = reply
            self.records.append({'url': safe_url(self.requests.pop(identity)),
                                 'bytes': len(data) if data is not None else 0,
                                 'sha256': hashlib.sha256(data).hexdigest() if data is not None else None,
                                 'error': str(failure) if failure else None})
        return reply


def manager(config, renderer):
    result = m.BaseManager(config['url'], renderer)
    x0, y0, x1, y1 = config['bounds']
    result.centerLat, result.centerLon = (y0 + y1) / 2, (x0 + x1) / 2
    result.setGeometricError(config['geometric_error'])
    result.concurrentDownloads = True
    result.errors = []
    return result


def inventory(config):
    tile_manager = manager(config, None)
    tile_manager.initAreaData(*config['bounds'])
    steps = tile_manager.iterContents()
    transport = MeasuredTransport(max_requests=1)
    selected = []
    value, failure = None, None
    started = time.monotonic()
    try:
        while True:
            try:
                command = steps.throw(failure) if failure else steps.send(value)
            except StopIteration as finished:
                assert not finished.value, finished.value
                break
            value, failure = None, None
            if isinstance(command, core.Download):
                selected.append(safe_url(command.request.url))
                assert len(selected) <= args.max_tiles, 'Live region exceeds tile budget'
            elif command is not None:
                transport.submit(command)
                while True:
                    reply = transport.poll()
                    if reply is not None:
                        _, value, failure = reply
                        break
                    assert time.monotonic() - started < 180, 'Manifest traversal timed out'
                    time.sleep(.005)
        assert not tile_manager.errors, tile_manager.errors
        return {'selected': len(selected), 'tile_urls': selected, 'manifests': transport.records,
                'seconds': time.monotonic() - started}
    finally:
        steps.close()
        transport.close()


def scene_data():
    bpy.context.view_layer.update()
    points, counts = [], Counter()
    for obj in bpy.context.scene.objects:
        if obj.type != 'MESH':
            continue
        mesh = obj.data
        counts.update(objects=1, vertices=len(mesh.vertices), faces=len(mesh.polygons), loops=len(mesh.loops))
        local = np.empty(len(mesh.vertices) * 3, dtype=np.float32)
        mesh.vertices.foreach_get('co', local)
        transform = np.array(obj.matrix_world, dtype=np.float64)
        points.append(local.reshape(-1, 3) @ transform[:3, :3].T + transform[:3, 3])
    textures = sorted(((list(image.size), hashlib.sha256(bytes(image.packed_file.data)).hexdigest()
                        if image.packed_file else None) for image in bpy.data.images), key=json.dumps)
    coords = np.concatenate(points) if points else np.empty((0, 3))
    assert np.isfinite(coords).all()
    return coords, dict(counts), textures


def run(config, progressive, workers, join):
    reset_scene()
    renderer = (p.ProgressiveRenderer if progressive else r.BlenderRenderer)('Live test', join, '')
    tile_manager = manager(config, renderer)
    tile_manager.downloadWorkers = workers
    transport = MeasuredTransport(max_requests=workers)
    process = transport.process
    operation = None
    started = time.monotonic()
    imported, commit_times, observed = [], [], []
    consume = tile_manager.consumeContent
    def commit(uri, *rest):
        before = time.monotonic()
        consume(uri, *rest)
        imported.append(safe_url(uri))
        commit_times.append(time.monotonic() - before)
        assert len(imported) <= args.max_tiles, 'Live import exceeds tile budget'
    tile_manager.consumeContent = commit
    try:
        with no_tile_writes_or_operator(no_dataset_files=True):
            if progressive:
                operation = p.TilesImport(tile_manager, config['bounds'], bpy.context, False, transport)
                while not operation.job.done:
                    operation.job.step()
                    observed.append(renderer.num_imported_tiles)
                    assert time.monotonic() - started < 240, 'Live import timed out'
                    time.sleep(.001)
                assert operation.job.error is None, operation.job.error
                result = operation.job.result
            else:
                with patch.object(m, 'DownloadTransport', return_value=transport):
                    result = tile_manager.render(*config['bounds'])
        assert len(result) == 2, result
        if result[0] and progressive:
            assert 1 in observed and result[0] in observed, 'No incremental tile commits'
        points, counts, textures = scene_data()
        assert not [c for c in bpy.data.collections if c.name.endswith('_import')]
        if result[0]:
            assert counts.get('vertices', 0) > 0
            assert counts['objects'] == 1 if join else counts['objects'] >= 1
        return {'mode': 'progressive' if progressive else 'blocking', 'workers': workers,
                'joined': join, 'tiles': result[0], 'errors': result[1],
                'seconds': time.monotonic() - started, 'commit_seconds': sum(commit_times),
                'longest_commit_seconds': max(commit_times, default=0),
                'peak_pending': transport.peak_pending, 'scene': counts,
                'textures': textures, 'requests': transport.records, 'import_order': imported}, points
    finally:
        if operation and not operation.job.done:
            operation.cancel()
        transport.close()
        process.wait(timeout=5)
        assert not t._live and not t.reap(), 'Download helper leaked'


try:
    for name in args.datasets:
        config = DATASETS[name]
        report = {'dataset': name, **config, 'runs': []}
        reports.append(report)
        try:
            if args.inventory:
                report.update(inventory(config), status='inventory')
            else:
                expected, tree, expected_counts, expected_textures, expected_hashes = None, None, None, None, None
                modes = [(False, 4, False), (True, 4, False)] if args.quick else [
                    (False, 1, False), (False, 4, False), (True, 1, False),
                    (True, 4, False), (False, 4, True), (True, 4, True)]
                for progressive, workers, join in modes:
                    result, points = run(config, progressive, workers, join)
                    report['runs'].append(result)
                    hashes = {item['url']: item['sha256'] for item in result['requests']}
                    if expected is None:
                        expected = points
                        expected_counts, expected_textures, expected_hashes = result['scene'], result['textures'], hashes
                        tree = KDTree(len(points))
                        for index, point in enumerate(points):
                            tree.insert(point, index)
                        tree.balance()
                    else:
                        assert hashes == expected_hashes, 'Live source payloads changed between runs'
                        assert result['textures'] == expected_textures, 'Texture contents changed'
                        for key in ('vertices', 'faces', 'loops'):
                            assert result['scene'].get(key, 0) == expected_counts.get(key, 0), key
                        assert result['tiles'] == report['runs'][0]['tiles']
                        assert len(result['errors']) == len(report['runs'][0]['errors'])
                        if len(points):
                            magnitude = max(1., float(np.abs(expected).max()))
                            tolerance = max(.002, ldexp(8., frexp(magnitude)[1] - 24))
                            distance = max(tree.find(point)[2] for point in points[::max(1, len(points) // 20000)])
                            result['sampled_max_geometry_error'] = distance
                            result['geometry_tolerance'] = tolerance
                            assert distance <= tolerance, (distance, tolerance)
                    write_report()
                    print('LIVE_RUN', name, progressive, workers, join, result['tiles'],
                          len(result['errors']), round(result['seconds'], 2), flush=True)
                expected_failure = config.get('expected_failure')
                if expected_failure:
                    assert report['runs'][0]['tiles'] == 0 and report['runs'][0]['errors']
                    assert all(expected_failure in error for error in report['runs'][0]['errors'])
                    report['status'] = 'passed_expected_failure'
                else:
                    assert report['runs'][0]['tiles'] > 0, 'No real tiles imported; see per-resource errors'
                    assert not report['runs'][0]['errors'], 'Dataset import has errors; see per-resource details'
                    report['status'] = 'passed'
        except Exception as exc:
            report['status'] = 'failed'
            report['error'] = str(exc)
            if hasattr(exc, 'live_probe'):
                report['transport_failure'] = exc.live_probe
            report['traceback'] = traceback.format_exc()
            print('LIVE_FAILURE', name, str(exc), flush=True)
        write_report()
        print('LIVE_DATASET', name, report['status'], report.get('selected', ''), flush=True)
finally:
    addon_utils.disable(package, default_set=True)

assert all(report['status'] in ('passed', 'passed_expected_failure', 'inventory') for report in reports), 'See live_checks.json'
