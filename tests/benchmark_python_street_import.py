"""Time the normal local-OSM street import in a fresh factory Blender process.

By default no cProfile, rendering, buildings or network requests. Optional --scene uses
its effective street settings and selected terrain, without saving that file.
Timings enclose the existing importer and two component-level calls without
changing geometry. --diagnostic-profile adds explicit profiler overhead.
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


def runtime_digest():
    """Identify uncommitted runtime edits as well as the Git parent revision."""
    digest = hashlib.sha256()
    for path in sorted(ROOT.rglob('*.py')):
        relative = path.relative_to(ROOT)
        if relative.parts[0] in ('tests', 'tools'):
            continue
        digest.update(relative.as_posix().encode('utf8') + b'\0' + path.read_bytes() + b'\0')
    return digest.hexdigest()


def operator_import(package, app, args, renderer_class, osm_class):
    """Use the real 3D-realistic Import operator, with streets as its only feature."""
    from contextlib import ExitStack
    from unittest.mock import patch
    props = bpy.context.scene.blosm
    for name, value in dict(mode='3Drealistic', dataType='osm', osmSource='file',
            osmFilepath=str(args.osm.resolve()), buildings=False, forests=False,
            vegetation=False, aviation=False, highways=True, railways=False,
            water=False, loadMissingMembers=False, setupScript='').items():
        setattr(props, name, value)
    preferences = bpy.context.preferences.addons[package.__name__].preferences
    preferences.preferPerformanceDuringImport = args.import_performance
    data = args.output / 'download_cache'
    data.mkdir(parents=True, exist_ok=True)
    preferences.dataDir = str(data)
    preferences.enableExperimentalFeatures = True
    stages, contexts = defaultdict(float), []

    def timed(original, name, capture=False):
        def call(self, *a, **kw):
            start = time.perf_counter()
            try:
                result = original(self, *a, **kw)
                if capture:
                    contexts.append(result)
                return result
            finally:
                stages[name] += time.perf_counter() - start
        return call

    with ExitStack() as stack:
        for name in ('__init__', 'prepare', 'render', 'finalize', 'cleanup'):
            stack.enter_context(patch.object(renderer_class, name,
                timed(getattr(renderer_class, name), 'street_generation_seconds', name == 'render')))
        for name in ('createLayers', 'initLayers', 'process', 'setTerrain'):
            # Patch the effective class only: do not time its super calls twice.
            cls = type(app)
            stack.enter_context(patch.object(cls, name, timed(getattr(cls, name), 'osm_preparation_seconds')))
        stack.enter_context(patch.object(osm_class, 'parse', timed(osm_class.parse, 'osm_preparation_seconds')))
        start, cpu = time.perf_counter(), time.process_time()
        assert bpy.ops.blosm.import_data() == {'FINISHED'}
        stages['operator_seconds'] = time.perf_counter() - start
        stages['prepare_and_generation_cpu_seconds'] = time.process_time() - cpu
    stages['operator_other_seconds'] = (stages['operator_seconds'] - stages['street_generation_seconds']
                                         - stages['osm_preparation_seconds'])
    assert len(contexts) == 1, len(contexts)
    return dict(stages), contexts[0]


def main():
    global ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--osm', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--addon-root', type=Path, help='Explicit read-only reference snapshot instead of this checkout')
    parser.add_argument('--resolution', type=int, choices=(1, 2, 3, 4), default=1)
    parser.add_argument('--cell-size', type=float, default=.7)
    parser.add_argument('--save-scene', action='store_true')
    parser.add_argument('--scene', type=Path, help='Read saved street settings and terrain with embedded scripts disabled')
    parser.add_argument('--without-terrain', action='store_true', help='Control run with the scene terrain selection cleared')
    parser.add_argument('--diagnostic-profile', action='store_true', help='Collect call profiles; these timings include profiler overhead')
    parser.add_argument('--profile-scope', choices=('numerical', 'publication', 'all'), default='numerical')
    parser.add_argument('--operator', action='store_true', help='Use the normal Import operator in 3D realistic, streets-only mode')
    parser.add_argument('--execution-diagnostics', action='store_true', help='Record coarse numerical call counts, CPU placement and calibration timings')
    parser.add_argument('--construction-diagnostics', action='store_true', help='Record chart attempts and bounded kernel replay inputs; adds instrumentation overhead')
    parser.add_argument('--qos', choices=('system', 'high', 'eco'), default='system', help='Explicit current-test-process scheduling experiment; original policy restored')
    parser.add_argument('--import-performance', action='store_true', help='Enable the production Import preference; controlled baselines disable it')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    if args.addon_root is not None:
        ROOT = args.addon_root.resolve()
        if not (ROOT / '__init__.py').is_file():
            parser.error('The selected addon snapshot has no __init__.py')
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
    args.output.mkdir(parents=True, exist_ok=True)
    scene_sha256 = hashlib.sha256(args.scene.read_bytes()).hexdigest() if args.scene else None
    source_sha256 = runtime_digest()
    # Resolve this worktree, regardless of its directory/package name.
    spec = importlib.util.spec_from_file_location('blosm', ROOT / '__init__.py',
        submodule_search_locations=[str(ROOT)])
    package = importlib.util.module_from_spec(spec)
    sys.modules['blosm'] = package
    spec.loader.exec_module(package)
    if args.scene:
        package.register()
        if args.operator:
            # Registration is explicit for the selected worktree. addon_utils
            # would rediscover the alias "blosm" in the user's sibling checkout.
            item = bpy.context.preferences.addons.new()
            item.module = package.__name__
        bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()), load_ui=False, use_scripts=False)
    from blosm.app import BaseApp
    from blosm.app.blender import BlenderApp
    from blosm.parse.osm import Osm
    from blosm.setup import Setup
    from blosm.action.generate_streets import StreetGenerator
    from blosm.way.renderer_streets import StreetRenderer
    from blosm.way.street_network.python_backend import VERSION
    from blosm.way.street_network.python_backend.backend import PythonMeshBackend
    from blosm.way.street_network.python_backend.scene import SceneWriter

    app = package.blenderApp.app if args.operator else BlenderApp()
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
    scene_settings = {}
    if args.scene:
        props = bpy.context.scene.blosm
        for prop in props.bl_rna.properties:
            name = prop.identifier
            if name.startswith('street') or name in ('terrainObject', 'trafficSide', 'assetLibraryMode'):
                value = getattr(props, name)
                setattr(app, name, value)
                scene_settings[name] = value.name if isinstance(value, bpy.types.ID) else value
        if args.without_terrain:
            app.terrainObject = ''
            props.terrainObject = ''
        # This entry point always measures the Python street renderer only,
        # irrespective of other feature switches or the saved UI mode.
        app.streetBackend = 'python'
    BaseApp.initOsm(app)
    osm = Osm(app)
    setup = Setup(app, osm)
    setup.skipWays()
    setup.roadsAndPaths()
    manager = setup.getWayManager()
    manager.addAction(StreetGenerator(manager.styleStore, getStyle=manager.getStyle, doPixelCorrection=False))

    stages, numerical, publication = {}, [], []
    qos = None
    if args.qos != 'system':
        from blosm.tests.street_execution_diagnostics import ProcessQoS
        qos = ProcessQoS(args.qos)
    diagnostics = None
    if args.execution_diagnostics:
        from blosm.tests.street_execution_diagnostics import ExecutionDiagnostics
        diagnostics = ExecutionDiagnostics(args.output)
        diagnostics.install(package.__name__)
        diagnostics.calibrate('before_import')
    profilers = {}
    construction = None
    if args.construction_diagnostics:
        from blosm.tests.street_construction_diagnostics import ConstructionDiagnostics
        construction = ConstructionDiagnostics(args.output)
        construction.install(package.__name__)
    if args.diagnostic_profile:
        import cProfile
        scopes = ('numerical', 'publication') if args.profile_scope == 'all' else (args.profile_scope,)
        profilers = {scope: cProfile.Profile() for scope in scopes}
    original_prepare, original_publish = PythonMeshBackend.prepare, SceneWriter.publish

    def prepared(self, *a, **kw):
        if diagnostics:
            diagnostics.calibrate('before_numerical_' + str(len(numerical)))
        start = time.perf_counter()
        if 'numerical' in profilers:
            profilers['numerical'].enable()
        try:
            result = original_prepare(self, *a, **kw)
        finally:
            if 'numerical' in profilers:
                profilers['numerical'].disable()
        numerical.append(dict(seconds=time.perf_counter() - start, timings=result.timings))
        if diagnostics:
            diagnostics.calibrate('after_numerical_' + str(len(numerical) - 1))
        print('STREET_BENCHMARK_NUMERICAL', numerical[-1]['seconds'], flush=True)
        return result

    def published(self, *a, **kw):
        start = time.perf_counter()
        if 'publication' in profilers:
            profilers['publication'].enable()
        try:
            result = original_publish(self, *a, **kw)
        finally:
            if 'publication' in profilers:
                profilers['publication'].disable()
        publication.append(time.perf_counter() - start)
        print('STREET_BENCHMARK_PUBLICATION', publication[-1], flush=True)
        return result

    PythonMeshBackend.prepare, SceneWriter.publish = prepared, published
    try:
        if args.operator:
            stages, context = operator_import(package, app, args, StreetRenderer, Osm)
        else:
            start = time.perf_counter()
            cpu = time.process_time()
            app.createLayers(osm)
            if args.scene and 'lat' in bpy.context.scene and 'lon' in bpy.context.scene and props.relativeToInitialImport:
                osm.setProjection(bpy.context.scene['lat'], bpy.context.scene['lon'])
            osm.parse(str(args.osm.resolve()), forceExtentCalculation=True)
            if args.scene:
                app.setTerrain(bpy.context, createFlatTerrain=False, createBvhTree=True)
            app.initLayers()
            app.process()
            stages['osm_preparation_seconds'] = time.perf_counter() - start
            print('STREET_BENCHMARK_OSM', stages['osm_preparation_seconds'], flush=True)
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
        try:
            if construction:
                construction.close()
            if diagnostics:
                diagnostics.close()
        finally:
            if qos:
                qos.close()

    stages['numerical_backend_seconds'] = sum(r['seconds'] for r in numerical)
    stages['blender_publication_seconds'] = sum(publication)
    stages['adaptation_materials_and_other_seconds'] = (stages['street_generation_seconds']
        - stages['numerical_backend_seconds'] - stages['blender_publication_seconds'])
    stages['total_ready_seconds'] = sum(stages[k] for k in
        ('osm_preparation_seconds', 'street_generation_seconds', 'depsgraph_update_seconds'))
    if args.operator:
        stages['total_ready_seconds'] = stages['operator_seconds'] + stages['depsgraph_update_seconds']
    if args.diagnostic_profile:
        import pstats
        for scope, profiler in profilers.items():
            profiler.dump_stats(str(args.output / (scope + '.prof')))
            with (args.output / (scope + '_profile.txt')).open('w', encoding='utf8') as stream:
                stats = pstats.Stats(profiler, stream=stream).strip_dirs()
                stats.sort_stats('cumulative').print_stats(80)
                stats.sort_stats('tottime').print_stats(60)
    report = context.report
    stages['terrain_placement_seconds'] = sum(s.get('timings', {}).get('terrain', {}).get('seconds', 0.)
                                             for s in report['surfaces'])
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
    assert runtime_digest() == source_sha256, 'Runtime files changed during the measured import'
    profile = Path(os.environ['BLENDER_USER_CONFIG']).resolve().parent
    assert all(Path(p).resolve().is_relative_to(profile) for p in resource_paths.values())
    output = dict(passed=True, timings=stages, task_seconds=dict(tasks), importer=report,
        meshes=len(meshes), quads=sum(len(o.data.polygons) for o in meshes),
        vertices=sum(len(o.data.vertices) for o in meshes), components=len(report['surfaces']),
        input=str(args.osm.resolve()), input_sha256=hashlib.sha256(args.osm.read_bytes()).hexdigest(),
        revision=args.revision, implementation=VERSION, checkout=str(ROOT), blender=bpy.app.version_string,
        runtime_source_sha256=source_sha256,
        source_scene_sha256=scene_sha256,
        execution_policy=qos.report if qos else None,
        settings=dict(resolution=app.streetResolution, cell_size=app.streetCellSize, interpolation='smooth cubic',
            crown_rise=app.streetProfileRise, exaggeration=app.streetExaggeration,
            terrain=getattr(app, 'terrainObject', False) or False,
            scene=str(args.scene.resolve()) if args.scene else None, saved_street_settings=scene_settings,
            diagnostic_profile=args.diagnostic_profile,
            profile_scope=args.profile_scope if args.diagnostic_profile else None,
            entry_point='3D realistic Import, streets only' if args.operator else 'prepared street renderer',
            execution_diagnostics=args.execution_diagnostics,
            construction_diagnostics=args.construction_diagnostics,
            qos=args.qos,
            import_performance=args.import_performance,
            materials='explicit selection' if getattr(app, 'streetCarriagewayMaterial', None) else 'library defaults',
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
