"""OSM session, real extracts, separate-object/GN and cancellation contracts."""
from pathlib import Path
import importlib
import random
import time
import json
import os
import sys
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bpy
from release_target import PROFILE
from release_support import (arguments, enable, catalog, configure_osm, geometry,
                             assert_geometry, meshes, run, reset, real_fixture, PACKAGE)

args = arguments()
addon = enable(args)
runtime = importlib.import_module(PACKAGE + '.app.osm_import')
core = importlib.import_module(PACKAGE + '.progressive.core')
metrics = []


def peak_rss():
    if os.name != 'nt':
        return None
    import ctypes as c
    from ctypes import wintypes as w
    class Counters(c.Structure):
        _fields_ = [('size', w.DWORD), ('faults', w.DWORD)] + [
            (name, c.c_size_t) for name in ('peak', 'current', 'paged_peak', 'paged',
                                          'nonpaged_peak', 'nonpaged', 'pagefile', 'pagefile_peak')]
    value = Counters()
    value.size = c.sizeof(value)
    query = c.WinDLL('psapi').GetProcessMemoryInfo
    query.argtypes = [w.HANDLE, c.c_void_p, w.DWORD]
    query.restype = w.BOOL
    return value.peak if query(c.c_void_p(-1), c.byref(value), value.size) else None


class OsmProgressiveChecks(unittest.TestCase):
    def prepare(self, fixture='shared_edges', mode='3Dsimple'):
        entry = next(e for e in catalog() if e['id'] == fixture)
        settings = configure_osm(args, entry, mode, False)
        return settings

    def session(self, steps=None):
        session = runtime.OsmImport(bpy.context)
        if steps is not None:
            session.job.steps.close()
            session.job.steps = steps(session)
        return session

    def drain(self, session, limit=200000):
        try:
            for _ in range(limit):
                if session.job.done:
                    break
                session.job.step()
                if session.phase == 'Rendering export textures':
                    time.sleep(.01)
            self.assertTrue(session.job.done, 'Session failed to terminate')
            self.assertIsNone(session.job.error, str(session.job.error))
        finally:
            session.finish(silent=True)
        self.assertFalse(session.errors, session.errors)
        metrics.append(dict(test=self.id(), peak_parent_rss_bytes=peak_rss(), **runtime.last_result))

    def test_real_buildings_and_2d_parity(self):
        for mode in ('2D', '3Dsimple'):
            self.prepare(mode=mode)
            self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
            baseline = geometry(meshes())
            settings = self.prepare(mode=mode)
            settings.progressiveImport = True
            settings.singleObject = True
            selected = bpy.data.objects.new('User selection', None)
            bpy.context.scene.collection.objects.link(selected)
            selected.select_set(True)
            bpy.context.view_layer.objects.active = selected
            session = self.session()
            counts = []
            while not session.job.done:
                session.job.step()
                counts.append(len(meshes()))
                self.assertEqual(bpy.context.view_layer.objects.active, selected)
                self.assertTrue(selected.select_get())
            session.finish(silent=True)
            self.assertFalse(session.errors, session.errors)
            if mode == '3Dsimple' and hasattr(settings, 'buildingBackend'):
                # Native PADL publication retains its whole-manager fallback.
                self.assertEqual(session.committed_units, 1)
                self.assertEqual(len(meshes()), 3)
            else:
                self.assertEqual(session.committed_units, 3)
                self.assertIn(1, counts)
                self.assertIn(2, counts)
            self.assertTrue(settings.singleObject, 'Saved blocking preference changed')
            assert_geometry(self, baseline, geometry(meshes()))
            self.assertIn('lat', bpy.context.scene)

    def test_all_core_extracts_separate_mode_parity(self):
        for entry in catalog():
            if entry['id'] in ('frankfurter_allee', 'kuzmice_railways'):
                continue  # Dedicated curve test below.
            with self.subTest(fixture=entry['id']):
                configure_osm(args, entry, '3Dsimple', False)
                self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
                expected = geometry(meshes())
                settings = configure_osm(args, entry, '3Dsimple', True)
                settings.progressiveImport = True
                self.drain(self.session())
                assert_geometry(self, expected, geometry(meshes()))

    def test_cancel_retains_one_and_user_edits(self):
        settings = self.prepare(mode='2D')
        settings.progressiveImport = True
        session = self.session()
        while not session.committed_units and not session.job.done:
            session.job.step()
        self.assertIsNone(session.job.error)
        first = meshes()[0]
        first.name = 'Renamed by user'
        sentinel = bpy.data.objects.new('Unrelated edit', None)
        bpy.context.scene.collection.objects.link(sentinel)
        session.cancel()
        result = session.finish(silent=True)
        self.assertEqual(result['outcome'], 'cancelled')
        self.assertIn(first.name, bpy.data.objects)
        self.assertIn(sentinel.name, bpy.data.objects)
        self.assertEqual(len(meshes()), 1)
        self.assertIn('lat', bpy.context.scene)
        self.assertIsNone(session.finish(silent=True))

    def test_destination_delete_and_unrelated_sentinel(self):
        settings = self.prepare()
        settings.progressiveImport = True
        session = self.session()
        while not session.committed_units and not session.job.done:
            session.job.step()
        sentinel = bpy.data.objects.new('Unrelated edit', None)
        bpy.context.scene.collection.objects.link(sentinel)
        bpy.data.collections.remove(session.destination)
        session.job.step()
        self.assertTrue(session.job.done)
        self.assertIn('collection', str(session.job.error))
        session.finish(silent=True)
        self.assertIn(sentinel.name, bpy.data.objects)

    def test_local_work_does_not_launch_downloader(self):
        self.prepare()
        session = self.session()
        with patch.object(runtime, 'DownloadTransport', side_effect=AssertionError('Unexpected helper')):
            self.drain(session)

    def test_panel_and_storage_settings_are_captured(self):
        settings = self.prepare()
        settings.progressiveImport = True
        prefs = bpy.context.preferences.addons[PACKAGE].preferences
        saved_directory = prefs.dataDir
        session = self.session()
        settings.buildings = False
        settings.singleObject = True
        prefs.dataDir = str(args.artifacts/'changed_after_start_missing')
        try:
            self.drain(session)
            self.assertEqual(len(meshes()), 3)
            self.assertFalse(settings.buildings)
            self.assertEqual(runtime.app.dataDir, saved_directory)
        finally:
            prefs.dataDir = saved_directory

    def test_original_errors_reach_both_import_reports(self):
        prefs = bpy.context.preferences.addons[PACKAGE].preferences
        saved_directory = prefs.dataDir
        saved_assets = prefs.assetsDir
        failures = ['data_directory', 'invalid_data_directory', 'osm_file', 'setup_script', 'corrupt_osm']
        if PROFILE != 'base':
            failures += ['assets_directory', 'invalid_assets_directory']
        for failure in failures:
            with self.subTest(failure=failure):
                settings = self.prepare()
                if failure == 'data_directory':
                    prefs.dataDir = ''
                elif failure == 'invalid_data_directory':
                    prefs.dataDir = str(args.artifacts/'missing_original_data')
                elif failure in ('assets_directory', 'invalid_assets_directory'):
                    settings.mode = '3Drealistic'
                    settings.assetsDir = '' if failure == 'assets_directory' else str(args.artifacts/'missing_original_assets')
                    prefs.assetsDir = ''
                elif failure == 'osm_file':
                    settings.osmFilepath = str(args.artifacts/'missing_original.osm')
                elif failure == 'setup_script':
                    settings.setupScript = str(args.artifacts/'missing_original_setup.py')
                else:
                    corrupt = args.artifacts/'corrupt_original.osm'
                    corrupt.write_text('<osm><node')
                    settings.osmFilepath = str(corrupt)
                try:
                    with self.assertRaises(RuntimeError) as blocking:
                        bpy.ops.blosm.import_data()
                    original = runtime.last_result['errors'][0]
                    blocking_message = runtime.last_result['message']
                    self.assertIn(original, str(blocking.exception))
                    settings.progressiveImport = True
                    session = self.session()
                    while not session.job.done:
                        session.job.step()
                    with patch.object(runtime.driver, 'notify_result') as notify:
                        result = session.finish()
                    notify.assert_called_once_with('ERROR', blocking_message)
                    self.assertEqual(result['message'], blocking_message)
                    self.assertIn(original, result['message'])
                    if failure == 'corrupt_osm':
                        self.assertFalse(result['preliminary_failure'])
                        self.assertIn('OSM import failed: 0 completed, 1 errors.', result['message'])
                    else:
                        self.assertTrue(result['preliminary_failure'])
                        self.assertEqual(result['message'], original)
                        self.assertNotIn('OSM import failed:', result['message'])
                    self.assertIsNone(session.finish())
                finally:
                    prefs.dataDir = saved_directory
                    prefs.assetsDir = saved_assets

    def test_first_turn_checks_prerequisites_without_initializing(self):
        settings = self.prepare()
        settings.progressiveImport = True
        session = self.session()
        try:
            with patch.object(runtime.app, 'initOsm') as initialize, \
                 patch.object(runtime, 'DownloadTransport') as download:
                session.job.step(max_steps=1)
            self.assertIsNone(session.job.error)
            self.assertFalse(session.job.done)
            self.assertFalse(session.owns_app)
            self.assertFalse(session.initialized)
            self.assertFalse(session.ownership.pending)
            initialize.assert_not_called()
            download.assert_not_called()
        finally:
            session.cancel()
            session.finish(silent=True)

    def test_missing_relation_members_post_and_companion_reuse(self):
        source = args.artifacts/'missing_members.osm'
        source.write_text('''<osm version="0.6"><bounds minlat="52" minlon="13" maxlat="52.001" maxlon="13.001"/>
<node id="1" lat="52" lon="13"/><node id="2" lat="52" lon="13.0003"/>
<node id="3" lat="52.0002" lon="13.0003"/><node id="4" lat="52.0002" lon="13"/>
<relation id="300"><member type="way" ref="200" role="outer"/>
<tag k="type" v="multipolygon"/><tag k="building" v="yes"/></relation></osm>''')
        response = b'<osm version="0.6"><way id="200"><nd ref="1"/><nd ref="2"/><nd ref="3"/><nd ref="4"/><nd ref="1"/></way></osm>'
        requests = []
        def download(session, request, path):
            requests.append(request)
            Path(path).write_bytes(response)
            yield None
        def configure():
            settings = self.prepare()
            settings.osmFilepath = str(source)
            settings.loadMissingMembers = True
            settings.progressiveImport = True
        configure()
        with patch.object(runtime.OsmImport, '_download', download):
            self.drain(self.session())
        self.assertEqual(len(requests), 1)
        self.assertEqual(requests[0].method, 'POST')
        self.assertIn('200', requests[0].body)
        self.assertTrue(meshes())
        configure()
        with patch.object(runtime.OsmImport, '_download', side_effect=AssertionError('Companion not reused')):
            self.drain(self.session())
        self.assertTrue(meshes())

    def test_cancel_before_first_step_and_start_failure(self):
        self.prepare()
        session = self.session()
        session.cancel()
        self.assertEqual(session.finish(silent=True)['outcome'], 'cancelled')
        self.assertFalse(session.ownership.pending)
        driver = importlib.import_module(PACKAGE+'.progressive.blender')
        session = self.session()
        with patch.object(bpy.app.timers, 'register', side_effect=RuntimeError('Injected timer failure')):
            with self.assertRaisesRegex(RuntimeError, 'timer failure'):
                driver.start(session, bpy.context)
        self.assertIsNone(driver.active)
        self.assertTrue(session.finished)
        retained = bpy.data.node_groups.new('User retained group', 'GeometryNodeTree')
        retained.use_fake_user = True
        def pending(session):
            def create():
                group = bpy.data.node_groups.new('Pending unused group', 'GeometryNodeTree')
                group.use_fake_user = True
            session.call(create)
            yield core.YieldToUI()
        session = self.session(pending)
        session.job.step()
        session.cancel()
        session.finish(silent=True)
        self.assertNotIn('Pending unused group', bpy.data.node_groups)
        self.assertIn(retained.name, bpy.data.node_groups)

    def test_osm_ui_hides_single_option_and_keeps_2d_gn_controls(self):
        gui = importlib.import_module(PACKAGE+'.gui')
        class Layout:
            def __init__(self):
                self.props = []
            def prop(self, obj, name, **kwargs):
                self.props.append(name)
                return self
            def __getattr__(self, name):
                return lambda *a, **kw: self
        class Panel:
            def __init__(self, layout):
                self.layout = layout
            def __getattr__(self, name):
                return lambda *a, **kw: None
        settings = self.prepare()
        modes = ['2D', '3Dsimple'] + (['3Drealistic'] if PROFILE != 'base' else [])
        for mode in modes:
            settings.mode = mode
            for progressive in (False, True):
                settings.progressiveImport = progressive
                settings.singleObject = True
                layout = Layout()
                gui.BLOSM_PT_Settings.drawOsm(Panel(layout), bpy.context)
                self.assertEqual('singleObject' in layout.props, not progressive)
                if mode == '2D':
                    self.assertIn('gnBlendFile2d', layout.props)

    def test_terrain_is_not_mutated_and_external_edit_stops_import(self):
        from mathutils import Matrix
        settings = self.prepare()
        mesh = bpy.data.meshes.new('Existing terrain')
        mesh.from_pydata([(-10000, -10000, 5), (10000, -10000, 5),
                          (10000, 10000, 5), (-10000, 10000, 5)], [], [(0, 1, 2, 3)])
        terrain = bpy.data.objects.new('Existing terrain', mesh)
        bpy.context.scene.collection.objects.link(terrain)
        terrain.location = (2, 3, 4)
        settings.terrainObject = terrain.name
        settings.progressiveImport = True
        original = [tuple(v.co) for v in mesh.vertices]
        location = tuple(terrain.location)
        session = self.session()
        while not session.committed_units and not session.job.done:
            session.job.step()
        self.assertIsNone(session.job.error, session.job.error)
        self.assertEqual([tuple(v.co) for v in mesh.vertices], original)
        self.assertEqual(tuple(terrain.location), location)
        mesh.vertices[0].co.z += 1
        session.terrain_dirty = True
        session.job.step()
        self.assertIn('terrain changed', str(session.job.error))
        session.finish(silent=True)
        self.assertEqual(mesh.vertices[0].co.z, original[0][2]+1)

    def test_bulk_owned_commit_and_pending_failure(self):
        self.prepare()
        def producer(session):
            session.initialized = True
            session.set_origin = False
            runtime.app.managers = runtime.app.renderers = []
            yield None
            def commit():
                for name in ('Bulk A', 'Bulk B'):
                    obj = bpy.data.objects.new(name, bpy.data.meshes.new(name))
                    bpy.context.scene.collection.objects.link(obj)
            yield from session.unit(commit)
            def fail():
                obj = bpy.data.objects.new('Pending', bpy.data.meshes.new('Pending'))
                bpy.context.scene.collection.objects.link(obj)
                raise RuntimeError('Injected pending failure')
            session.call(fail)
        session = self.session(producer)
        while not session.committed_units:
            session.job.step()
        sentinel = bpy.data.objects.new('User sentinel', None)
        bpy.context.scene.collection.objects.link(sentinel)
        while not session.job.done:
            session.job.step()
        result = session.finish(silent=True)
        self.assertEqual(result['outcome'], 'partial')
        self.assertEqual(result['committed'], 1)
        self.assertIn('Injected pending failure', result['message'])
        self.assertIn('Bulk A', bpy.data.objects)
        self.assertIn('Bulk B', bpy.data.objects)
        self.assertIn('User sentinel', bpy.data.objects)
        self.assertNotIn('Pending', bpy.data.objects)
        self.assertNotIn('Pending', bpy.data.meshes)

    def test_new_backend_override_uses_atomic_fallback(self):
        self.prepare()
        calls = []
        class Original:
            def render(self):
                pass
            def render_steps(self, session):
                raise AssertionError('Inherited producer bypassed the backend override')
        class NewBackend(Original):
            def render(self):
                calls.append('new backend')
                obj = bpy.data.objects.new('Backend unit', bpy.data.meshes.new('Backend unit'))
                bpy.context.scene.collection.objects.link(obj)
        def steps(session):
            runtime.app.managers = runtime.app.renderers = []
            session.set_origin = False
            yield from session._hook(NewBackend(), 'render')
        session = self.session(steps)
        self.drain(session)
        self.assertEqual(calls, ['new backend'])
        self.assertEqual(session.committed_units, 1)

    def test_manager_preparation_order_even_without_way_output(self):
        self.prepare()
        sequence = []
        class Ways:
            def process(self):
                sequence.append('prepared ways')
        class Buildings:
            def process(self):
                self.assert_ready()
                sequence.append('prepared buildings')
            def assert_ready(self):
                if sequence != ['prepared ways']:
                    raise AssertionError('Building preparation ran before ways')
        def steps(session):
            runtime.app.managers = runtime.app.renderers = []
            yield from session._hook(Ways(), 'process')
            yield from session._hook(Buildings(), 'process')
        self.drain(self.session(steps))
        self.assertEqual(sequence, ['prepared ways', 'prepared buildings'])

    def test_random_stream_isolated_between_turns(self):
        self.prepare()
        values = []
        def producer(session):
            runtime.app.managers = runtime.app.renderers = []
            for _ in range(3):
                values.append(random.random())
                yield core.YieldToUI()
        random.seed(72)
        expected = [random.random() for _ in range(3)]
        random.seed(72)
        session = self.session(producer)
        while not session.job.done:
            session.job.step()
            random.random()
        session.finish(silent=True)
        self.assertEqual(values, expected)

    def test_background_entry_requires_explicit_blocking(self):
        settings = self.prepare()
        settings.progressiveImport = True
        with self.assertRaises(RuntimeError):
            bpy.ops.blosm.import_data()
        settings.commandLineMode = True
        settings.singleObject = True
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        self.assertEqual(len(meshes()), 3)
        settings.progressiveImport = settings.commandLineMode = False
        bpy.context.view_layer.objects.active = meshes()[0]
        meshes()[0].select_set(True)
        bpy.ops.object.mode_set(mode='EDIT')
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        self.assertEqual(bpy.context.mode, 'OBJECT')

    def test_native_gn_per_object_face_attributes(self):
        settings = self.prepare(mode='2D')
        settings.gnBlendFile2d = str(real_fixture(args, 'gn2d'))
        settings.gnSetup2d = 'building'
        settings.progressiveImport = True
        self.drain(self.session())
        layer_type = importlib.import_module(PACKAGE+'.building.gn_2d').GnBldg2dLayer
        unused = layer_type('unused_gn_layer', runtime.app)
        unused.finalize()  # A registered layer can legitimately receive no geometry.
        objects = meshes()
        self.assertEqual(len(objects), 3)
        groups = []
        for obj in objects:
            attribute = obj.data.attributes['building:levels']
            self.assertEqual(len(attribute.data), len(obj.data.polygons))
            self.assertTrue(all(value.value > 0 for value in attribute.data))
            groups.append(obj.modifiers[0].node_group)
        self.assertTrue(all(group == groups[0] for group in groups))

    def test_real_roads_and_multiline_parity(self):
        for fixture in ('frankfurter_allee', 'kuzmice_railways'):
            def config(progressive):
                settings = self.prepare(fixture)
                settings.buildings = False
                settings.highways = settings.railways = True
                settings.progressiveImport = progressive
            def points():
                bpy.context.view_layer.update()
                return sorted(tuple(round(c, 4) for c in obj.matrix_world @ p.co.xyz)
                              for obj in bpy.context.scene.objects if obj.type == 'CURVE'
                              for spline in obj.data.splines for p in spline.points)
            config(False)
            self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
            expected = points()
            self.assertTrue(expected)
            config(True)
            self.drain(self.session())
            self.assertEqual(points(), expected)

    def test_forest_patches_and_partial_cancel(self):
        self.assertIsNotNone(args.assets, 'Pro qualification requires --assets')
        import xml.etree.ElementTree as ET
        from math import ceil
        root = ET.Element('osm', version='0.6')
        for index, x in enumerate((0., .0007)):
            refs = []
            for n, (dy, dx) in enumerate(((0, 0), (0, .0005), (.0005, .0005), (.0005, 0))):
                identity = str(index*10+n+1)
                refs.append(identity)
                ET.SubElement(root, 'node', id=identity, lat=str(52+dy), lon=str(13+x+dx))
            way = ET.SubElement(root, 'way', id=str(100+index))
            for ref in refs+[refs[0]]:
                ET.SubElement(way, 'nd', ref=ref)
            ET.SubElement(way, 'tag', k='landuse', v='forest')
        path = args.artifacts/'forest.osm'
        ET.ElementTree(root).write(path)
        settings = reset(args)
        settings.dataType, settings.mode = 'osm', '3Drealistic'
        settings.osmSource, settings.osmFilepath = 'file', str(path)
        settings.assetsDir = str(args.assets)
        settings.buildings, settings.forests = False, True
        settings.loadMissingMembers = False
        settings.progressiveImport = True
        settings.singleObject = True
        bpy.context.preferences.addons[PACKAGE].preferences.enableExperimentalFeatures = False
        session = self.session()
        while session.committed_units < 1 and not session.job.done:
            session.job.step()
        self.assertIsNone(session.job.error, session.job.error)
        terrain = runtime.app.terrain.terrain
        self.assertEqual(len(terrain.particle_systems), 1)
        initial_count = terrain.particle_systems[0].settings.count
        self.assertGreater(initial_count, 0)
        self.drain(session)
        self.assertEqual(session.committed_units, 2)
        self.assertEqual(len(terrain.particle_systems), 1)
        self.assertGreater(terrain.particle_systems[0].settings.count, initial_count)
        self.assertTrue(settings.singleObject)
        bpy.context.view_layer.update()
        evaluated = terrain.evaluated_get(bpy.context.evaluated_depsgraph_get())
        self.assertGreater(len(evaluated.particle_systems[0].particles), 0)
        # Start another import on that existing terrain; cancel retains only
        # its completed patch and does not change the previous system.
        settings.terrainObject = terrain.name
        previous_count = terrain.particle_systems[0].settings.count
        session = self.session()
        while not session.committed_units and not session.job.done:
            session.job.step()
        self.assertIsNone(session.job.error, session.job.error)
        retained_count = terrain.particle_systems[-1].settings.count
        # Fail after changing cumulative density: slot/count rollback is local.
        with patch.object(session, '_terrain_key', side_effect=RuntimeError('Injected forest publication failure')):
            session.job.step()
        self.assertTrue(session.job.done)
        session.finish(silent=True)
        self.assertEqual(terrain.particle_systems[0].settings.count, previous_count)
        self.assertEqual(terrain.particle_systems[-1].settings.count, retained_count)

    def test_individual_trees_once_with_colliding_asset_name(self):
        self.assertIsNotNone(args.assets, 'Pro qualification requires --assets')
        path = args.artifacts/'trees.osm'
        path.write_text('''<osm version="0.6">
<node id="1" lat="52.0003" lon="13.0003"><tag k="natural" v="tree"/></node>
<node id="2" lat="52.0006" lon="13.0006"><tag k="natural" v="tree"/></node></osm>''')
        settings = reset(args)
        settings.mode, settings.dataType = '3Drealistic', 'osm'
        settings.osmSource, settings.osmFilepath = 'file', str(path)
        settings.assetsDir = str(args.assets)
        settings.forests, settings.buildings = True, False
        settings.loadMissingMembers = False
        settings.progressiveImport = True
        bpy.context.preferences.addons[PACKAGE].preferences.enableExperimentalFeatures = False
        collision = bpy.data.collections.new('trees')
        self.drain(self.session())
        trees = [obj for obj in bpy.context.scene.objects if obj.get('natural') == 'tree']
        self.assertEqual(len(trees), 2)
        self.assertTrue(all(obj.data is not None for obj in trees))
        self.assertEqual(len(collision.objects), 0)

    def test_realistic_and_uncached_export(self):
        self.assertIsNotNone(args.assets, 'Pro qualification requires --assets')
        if hasattr(runtime.app, 'preferableResult'):
            # Qualify the branch's Python volume/material path. Its legacy
            # footprint-GN groups require a separate matching binary asset.
            runtime.app.preferableResult = importlib.import_module(PACKAGE+'.defs').Result.VolumeMesh
        for export in (False, True):
            expected = None
            expected_uv = None
            if export:
                settings = self.prepare('rectangle', '3Drealistic')
                settings.assetsDir = str(args.assets)
                settings.importForExport = True
                prefs = bpy.context.preferences.addons[PACKAGE].preferences
                prefs.enableExperimentalFeatures = True
                baseline = args.artifacts/'blocking_export'
                baseline.mkdir(exist_ok=True)
                prefs.dataDir = str(baseline)
                self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
                expected = geometry(meshes())
                expected_uv = sorted(tuple(round(c, 6) for c in data.uv)
                                     for obj in meshes() for uv in obj.data.uv_layers for data in uv.data)
            settings = self.prepare('rectangle', '3Drealistic')
            settings.assetsDir = str(args.assets)
            settings.importForExport = export
            settings.progressiveImport = True
            prefs = bpy.context.preferences.addons[PACKAGE].preferences
            prefs.enableExperimentalFeatures = True
            folder = args.artifacts/('uncached_export' if export else 'realistic')
            folder.mkdir(exist_ok=True)
            prefs.dataDir = str(folder)
            session = self.session()
            self.drain(session)
            objects = meshes()
            self.assertTrue(objects)
            self.assertTrue(all(obj.data.materials for obj in objects))
            self.assertTrue(all(not obj.hide_viewport for obj in objects))
            if export:
                assert_geometry(self, expected, geometry(objects))
                actual_uv = sorted(tuple(round(c, 6) for c in data.uv)
                                   for obj in meshes() for uv in obj.data.uv_layers for data in uv.data)
                self.assertEqual(actual_uv, expected_uv)
                self.assertTrue(list(folder.rglob('*.png')) + list(folder.rglob('*.jpg')))
                self.assertFalse(list(folder.glob('texture/.blosm-bakes-*')))
                self.assertFalse([image for image in bpy.data.images if image.source == 'FILE' and image.size[0] == 1])


if PROFILE == 'base':
    # Edition selection is explicit, rather than marking a required case skipped.
    for method in ('test_forest_patches_and_partial_cancel', 'test_individual_trees_once_with_colliding_asset_name',
                   'test_realistic_and_uncached_export'):
        delattr(OsmProgressiveChecks, method)
run(args, OsmProgressiveChecks)
(args.artifacts/'osm_metrics.json').write_text(json.dumps(metrics, indent=2))
