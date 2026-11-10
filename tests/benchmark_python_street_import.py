"""Time the normal local-OSM street import in a fresh factory Blender process.

No cProfile, rendering, terrain, buildings or network requests. Timings enclose
the existing importer and two component-level calls without changing geometry.
Always launch with this chat's isolated Blender profile. JSON goes to artifacts.
"""
import argparse
from collections import defaultdict
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

import bpy


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--osm', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--resolution', type=int, choices=(1, 2, 3, 4), default=1)
    parser.add_argument('--cell-size', type=float, default=.7)
    parser.add_argument('--save-scene', action='store_true')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    # Resolve this worktree, regardless of its directory/package name.
    spec = importlib.util.spec_from_file_location('blosm', ROOT / '__init__.py',
        submodule_search_locations=[str(ROOT)])
    package = importlib.util.module_from_spec(spec)
    sys.modules['blosm'] = package
    spec.loader.exec_module(package)
    from blosm.app import BaseApp
    from blosm.app.blender import BlenderApp
    from blosm.parse.osm import Osm
    from blosm.setup import Setup
    from blosm.action.generate_streets import StreetGenerator
    from blosm.way.renderer_streets import StreetRenderer
    from blosm.way.street_network.python_backend import VERSION
    from blosm.way.street_network.python_backend.backend import PythonMeshBackend
    from blosm.way.street_network.python_backend.scene import SceneWriter

    app = BlenderApp()
    app.mode = BaseApp.realistic
    app.loadMissingMembers = False
    app.straightAngleThreshold = 175.5
    app.streetBackend = 'python'
    app.osmFilepath = str(args.osm.resolve())
    app.assetsDir = ''
    app.pmlFilepathStreet = None
    app.trafficSide = 'auto'
    app.assetLibraryMode = 'APPEND'
    app.streetCellSize = args.cell_size
    app.streetResolution = args.resolution
    app.streetExaggeration = 1.
    app.streetProfileRise = .08
    BaseApp.initOsm(app)
    osm = Osm(app)
    setup = Setup(app, osm)
    setup.skipWays()
    setup.roadsAndPaths()
    manager = setup.getWayManager()
    manager.addAction(StreetGenerator(manager.styleStore, getStyle=manager.getStyle, doPixelCorrection=False))

    stages, numerical, publication = {}, [], []
    original_prepare, original_publish = PythonMeshBackend.prepare, SceneWriter.publish

    def prepared(self, *a, **kw):
        start = time.perf_counter()
        result = original_prepare(self, *a, **kw)
        numerical.append(dict(seconds=time.perf_counter() - start, timings=result.timings))
        return result

    def published(self, *a, **kw):
        start = time.perf_counter()
        result = original_publish(self, *a, **kw)
        publication.append(time.perf_counter() - start)
        return result

    PythonMeshBackend.prepare, SceneWriter.publish = prepared, published
    try:
        start = time.perf_counter()
        cpu = time.process_time()
        app.createLayers(osm)
        osm.parse(str(args.osm.resolve()), forceExtentCalculation=True)
        app.initLayers()
        app.process()
        stages['osm_preparation_seconds'] = time.perf_counter() - start
        start = time.perf_counter()
        renderer = StreetRenderer(app)
        renderer.prepare()
        context = renderer.render(manager, osm)
        renderer.finalize()
        renderer.cleanup()
        stages['street_generation_seconds'] = time.perf_counter() - start
        stages['prepare_and_generation_cpu_seconds'] = time.process_time() - cpu
        start = time.perf_counter()
        bpy.context.view_layer.update()
        stages['depsgraph_update_seconds'] = time.perf_counter() - start
    finally:
        PythonMeshBackend.prepare, SceneWriter.publish = original_prepare, original_publish

    stages['numerical_backend_seconds'] = sum(r['seconds'] for r in numerical)
    stages['blender_publication_seconds'] = sum(publication)
    stages['adaptation_materials_and_other_seconds'] = (stages['street_generation_seconds']
        - stages['numerical_backend_seconds'] - stages['blender_publication_seconds'])
    stages['total_ready_seconds'] = sum(stages[k] for k in
        ('osm_preparation_seconds', 'street_generation_seconds', 'depsgraph_update_seconds'))
    report = context.report
    meshes = [o for o in context.collection.all_objects if o.type == 'MESH' and o.get('st_py_key')]
    assert report['coverage'] == 'COMPLETE', report['preparation_issues']
    assert all(s['status'] == 'READY' for s in report['surfaces']), report['surfaces']
    assert meshes and all(o.data.uv_layers.get('UVMap') and o.data.materials for o in meshes)
    assert all(p.loop_total == 4 for obj in meshes for p in obj.data.polygons)
    tasks = defaultdict(float)
    for component in numerical:
        for key, seconds in component['timings']['tasks'].items():
            tasks[key.split(':', 1)[0]] += seconds
    resource_paths = {kind: bpy.utils.user_resource(kind) for kind in ('CONFIG', 'SCRIPTS', 'DATAFILES', 'EXTENSIONS')}
    profile = Path(os.environ['BLENDER_USER_CONFIG']).resolve().parent
    assert all(Path(p).resolve().is_relative_to(profile) for p in resource_paths.values())
    output = dict(passed=True, timings=stages, task_seconds=dict(tasks), importer=report,
        meshes=len(meshes), quads=sum(len(o.data.polygons) for o in meshes),
        vertices=sum(len(o.data.vertices) for o in meshes), components=len(report['surfaces']),
        input=str(args.osm.resolve()), input_sha256=hashlib.sha256(args.osm.read_bytes()).hexdigest(),
        revision=args.revision, implementation=VERSION, checkout=str(ROOT), blender=bpy.app.version_string,
        settings=dict(resolution=args.resolution, cell_size=args.cell_size, interpolation='smooth cubic',
            crown_rise=.08, exaggeration=1., terrain=False, materials='library defaults',
            blender_threads=bpy.context.scene.render.threads,
            numeric_thread_environment={k: os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')}),
        profile_paths=resource_paths,
        exclusions=['Blender startup', 'Python/module setup', 'validation assertions', 'file saving', 'rendering', 'UI redraw'])
    if args.save_scene:
        bpy.context.preferences.filepaths.file_preview_type = 'NONE'
        start = time.perf_counter()
        bpy.ops.wm.save_as_mainfile(filepath=str(args.output / 'street_benchmark.blend'))
        output['save_seconds'] = time.perf_counter() - start
    (args.output / 'timings.json').write_text(json.dumps(output, indent=2), encoding='utf8')
    print('STREET_BENCHMARK_OK', json.dumps(dict(timings=stages, task_seconds=dict(tasks),
        meshes=output['meshes'], quads=output['quads'], coverage=report['coverage'])), flush=True)


if __name__ == '__main__':
    main()
