"""Shared, offline release checks. All output belongs in --artifacts."""
import argparse
from contextlib import contextmanager, ExitStack
import hashlib
import importlib
import json
import math
from pathlib import Path
import random
import socket
import sys
import time
import tomllib
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
from release_target import REPO, ADDON_ROOT, PACKAGE, PROFILE, assert_target
DEFAULT_ARTIFACTS = REPO.parent / "artifacts" / "blosm_release_tests"


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--artifacts', required=True, type=Path)
    parser.add_argument('--osm-root', type=Path, default=REPO.parent / 'osm_extracts')
    parser.add_argument('--fixture-root', type=Path, default=DEFAULT_ARTIFACTS / 'fixtures')
    parser.add_argument('--assets', type=Path)
    parser.add_argument('--candidate-zip', type=Path)
    parser.add_argument('--extended', action='store_true')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.artifacts.mkdir(parents=True, exist_ok=True)
    return args


def catalog(extended=False):
    with (REPO / 'tests/fixtures/osm_catalog.toml').open('rb') as handle:
        data = tomllib.load(handle)
    return [item for item in data['osm'] if extended or item['tier'] == 'core']


def verified(path, digest):
    if not path.is_file():
        raise FileNotFoundError(f'Required fixture missing: {path}')
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != digest:
        raise ValueError(f'Fixture changed: {path}; expected {digest}, got {actual}')
    return path


def osm_path(args, entry):
    path = verified(args.osm_root / entry['path'], entry['sha256'])
    if entry.get('companion_path'):
        verified(args.osm_root / entry['companion_path'], entry['companion_sha256'])
    return path


def real_fixture(args, name):
    with (REPO / 'tests/fixtures/local_catalog.toml').open('rb') as handle:
        entry = tomllib.load(handle)[name]
    return verified(args.fixture_root / entry['path'], entry['sha256'])


@contextmanager
def offline():
    def reject(*args, **kwargs):
        raise AssertionError('Unexpected network access in offline release test')
    with ExitStack() as stack:
        for target in ('socket.create_connection', 'socket.socket.connect',
                       'socket.socket.connect_ex', 'urllib.request.urlopen'):
            stack.enter_context(patch(target, side_effect=reject))
        yield


def enable(args):
    import bpy
    import addon_utils
    if PROFILE == 'general':
        sys.path.insert(0, str(REPO.parent))
    addon_utils.enable(PACKAGE, default_set=True)
    assert PACKAGE in bpy.context.preferences.addons, 'Addon registration failed'
    addon = importlib.import_module(PACKAGE)
    assert_target(addon)
    prefs = bpy.context.preferences.addons[PACKAGE].preferences
    (args.artifacts / 'data').mkdir(parents=True, exist_ok=True)
    prefs.dataDir = str(args.artifacts / 'data')
    if args.assets:
        prefs.assetsDir = str(args.assets)
    return addon


def reset(args):
    import bpy
    if bpy.context.object and bpy.context.object.mode != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    bpy.data.orphans_purge(do_recursive=True)
    scene = bpy.context.scene
    for key in ('lat', 'lon'):
        if key in scene:
            del scene[key]
    scene.property_unset('blosm')
    settings = scene.blosm
    settings.gnBlendFile2d = ''
    settings.gnSetup2d = '-'
    settings.setupScript = ''
    settings.terrainObject = ''
    settings.coordinatesAsFilter = False
    importlib.import_module(PACKAGE + '.gui').addDefaultLevels()
    random.seed(1729)
    return settings


def configure_osm(args, entry, mode='3Dsimple', single=True):
    settings = reset(args)
    settings.dataType = 'osm'
    settings.osmSource = 'file'
    settings.osmFilepath = str(osm_path(args, entry))
    settings.mode = mode
    settings.singleObject = single
    settings.loadMissingMembers = entry['load_missing_members']
    for layer in ('buildings', 'water', 'forests', 'vegetation', 'highways', 'railways'):
        setattr(settings, layer, layer == 'buildings')
    return settings


def meshes():
    import bpy
    return [obj for obj in bpy.context.scene.objects if obj.type == 'MESH']


def geometry(objects):
    """World-space face corners, independent of object names/origins and winding."""
    import bpy
    bpy.context.view_layer.update()
    faces = []
    for obj in objects:
        mesh = obj.data
        assert mesh.vertices and mesh.polygons, f'Empty mesh: {obj.name}'
        for vertex in mesh.vertices:
            assert all(math.isfinite(v) for v in vertex.co), f'Nonfinite mesh: {obj.name}'
        for face in mesh.polygons:
            assert face.area > 1e-9, f'Degenerate face: {obj.name}/{face.index}'
            corners = [tuple(obj.matrix_world @ mesh.vertices[i].co) for i in face.vertices]
            faces.append(sorted(corners))
    return sorted(faces)


def assert_geometry(test, expected, actual, tolerance=0.01):
    test.assertEqual(len(expected), len(actual), 'Face count changed with grouping')
    # Spatial matching avoids unstable lexicographic ordering of nearly tied
    # coordinates. Consume matches, preserving coincident duplicate face counts.
    from mathutils.kdtree import KDTree
    center = lambda face: tuple(sum(p[i] for p in face)/len(face) for i in range(3))
    tree = KDTree(len(actual))
    for i, face in enumerate(actual):
        tree.insert(center(face), i)
    tree.balance()
    consumed = set()
    for face in expected:
        candidates = [i for _, i, _ in tree.find_range(center(face), tolerance)
                      if i not in consumed and len(actual[i]) == len(face)]
        match = next((i for i in candidates if all(
            min(math.dist(p, q) for q in actual[i]) <= tolerance for p in face)), None)
        test.assertIsNotNone(match, f'No matching world-space face near {center(face)}')
        consumed.add(match)


class ReportResult(unittest.TextTestResult):
    def startTest(self, test):
        super().startTest(test)
        self.started = time.monotonic()

    def record(self, test, status, detail=None):
        self.records.append(dict(id=test.id(), status=status, detail=detail,
                                 seconds=round(time.monotonic() - self.started, 3)))

    def addSuccess(self, test):
        super().addSuccess(test)
        self.record(test, 'passed')

    def addFailure(self, test, err):
        super().addFailure(test, err)
        self.record(test, 'failed', self._exc_info_to_string(err, test))

    def addError(self, test, err):
        super().addError(test, err)
        self.record(test, 'error', self._exc_info_to_string(err, test))

    def addSkip(self, test, reason):
        super().addSkip(test, reason)
        self.record(test, 'skipped', reason)


def run(args, *classes):
    import bpy
    suite = unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(c) for c in classes)
    ReportResult.records = []
    with offline():
        result = unittest.TextTestRunner(verbosity=2, resultclass=ReportResult).run(suite)
    report = dict(blender=bpy.app.version_string, source=str(REPO),
                  profile=PROFILE, package=PACKAGE, runtime_root=str(ADDON_ROOT),
                  status='passed' if result.wasSuccessful() and not result.skipped else 'failed',
                  cases=result.records)
    (args.artifacts / 'results.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    if report['status'] != 'passed':
        raise RuntimeError('Release checks failed; see results.json')
