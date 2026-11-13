"""Real OSM street producer parity, ownership and cancellation in isolated Blender.

Use --artifacts for every output and optionally repeat --osm with local real
extracts. No network or installed-addon lookup is used. --only selects test
method substrings for targeted qualification. Run only with the chat's factory
profile; the script deliberately removes all test-scene data between cases.
"""
import argparse
from contextlib import contextmanager
from dataclasses import replace
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import sys
import time
import unittest
from unittest.mock import patch
from xml.etree import ElementTree as ET

import bpy
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'tests'))


def enable():
    spec = importlib.util.spec_from_file_location('blosm', ROOT / '__init__.py',
                                                submodule_search_locations=[str(ROOT)])
    package = importlib.util.module_from_spec(spec)
    sys.modules['blosm'] = package
    spec.loader.exec_module(package)
    package.register()
    item = bpy.context.preferences.addons.new()
    item.module = package.__name__
    # Imports allocate historically process-global source IDs. Two controlled
    # parity runs must begin with identical counters, as fresh Blender runs do.
    # Preserve every emitted identity in the comparator instead of erasing IDs.
    for name in ('connectors', 'bundle', 'crosswalk', 'corner', 'intersection',
                 'dummy_node', 'pt_stop', 'section', 'side_lane', 'roadway',
                 'street', 'partial_intersection', 'sym_lane'):
        importlib.import_module('blosm.way.item.' + name)
    importlib.import_module('blosm.action.generate_streets')
    global SOURCE_IDS
    SOURCE_IDS = {value: value.ID for name, module in tuple(sys.modules.items())
                  if name.startswith('blosm.way') and module is not None
                  for value in vars(module).values()
                  if isinstance(value, type) and value.__module__ == name
                  and isinstance(vars(value).get('ID'), int)}
    return package


def synthetic_osm(path):
    """Three disconnected T networks: the scheduler must commit them separately."""
    osm = ET.Element('osm', version='0.6')
    ET.SubElement(osm, 'bounds', minlat='51.999', minlon='12.999', maxlat='52.004', maxlon='13.005')
    for component in range(3):
        base = component * 10
        x0, y0 = component * 120., 0.
        points = ((x0, y0), (x0 + 32., y0), (x0, y0 + 32.), (x0 - 32., y0))
        for number, (x, y) in enumerate(points):
            ET.SubElement(osm, 'node', id=str(base + number + 1),
                          lat=str(52. + y / 111320.), lon=str(13. + x / 68535.))
        for arm in range(3):
            way = ET.SubElement(osm, 'way', id=str(100 + base + arm))
            ET.SubElement(way, 'nd', ref=str(base + 1))
            ET.SubElement(way, 'nd', ref=str(base + arm + 2))
            ET.SubElement(way, 'tag', k='highway', v='residential')
            ET.SubElement(way, 'tag', k='name', v='Fixture %s arm %s' % (component, arm))
            ET.SubElement(way, 'tag', k='lanes', v='2')
    ET.ElementTree(osm).write(path, encoding='utf-8', xml_declaration=True)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def buffer(items, attribute, dtype, width=1):
    value = np.empty(len(items) * width, dtype=dtype)
    items.foreach_get(attribute, value)
    return hashlib.sha256(value.tobytes()).hexdigest()


def mesh_state(mesh):
    """Exact stored geometry, normals, shading, face materials and every UV map."""
    return dict(vertices=buffer(mesh.vertices, 'co', np.float32, 3),
        loops=buffer(mesh.loops, 'vertex_index', np.int32),
        sizes=buffer(mesh.polygons, 'loop_total', np.int32),
        normals=buffer(mesh.corner_normals, 'vector', np.float32, 3),
        smooth=buffer(mesh.polygons, 'use_smooth', bool),
        face_materials=buffer(mesh.polygons, 'material_index', np.int32),
        materials=[material.name if material else None for material in mesh.materials],
        uv={layer.name: buffer(layer.uv, 'vector', np.float32, 2) for layer in mesh.uv_layers},
        relief={attribute.name: buffer(attribute.data, 'value', np.float32)
                for attribute in mesh.attributes if attribute.name == 'st_terrain_prepared_z'},
        counts=(len(mesh.vertices), len(mesh.polygons)))


def controller_state(controller):
    source = json.loads(controller['st_py_source_data'])
    # A fresh import UUID is the only source field ignored. Construction
    # metadata, source provenance, section IDs and topology remain exact.
    source.pop('uuid', None)
    manifest = json.loads(controller.get('st_py_manifest', '{}'))
    manifest.pop('network_uuid', None)
    outputs = {}
    for obj in controller['st_py_outputs'].objects:
        key = obj['st_py_key']
        assert key not in outputs, key
        output = dict(type=obj.type, matrix=tuple(tuple(row) for row in obj.matrix_world),
                      owner=obj.get('st_py_owner'), metadata=json.loads(obj.get('st_py_metadata', '{}')))
        if obj.type == 'MESH':
            output['mesh'] = mesh_state(obj.data)
            if obj.modifiers:
                output['evaluated_mesh'] = mesh_state(obj.evaluated_get(bpy.context.evaluated_depsgraph_get()).data)
        else:
            output['asset'] = obj.get('st_py_asset')
        outputs[key] = output
    return dict(source=source, manifest=manifest, outputs=outputs,
                terrain_mode=controller.get('st_py_terrain_mode'))


def controllers(committed=False):
    return [obj for obj in bpy.context.scene.objects if obj.get('st_py_network_uuid')
            and (not committed or obj.get('st_py_manifest'))]


def scene_state():
    bpy.context.view_layer.update()
    result = {}
    for controller in controllers(committed=True):
        state = controller_state(controller)
        key = ','.join(str(section['id']) for section in state['source']['sections'])
        assert key not in result, key
        result[key] = state
    return result


def inventory():
    return {kind: {item.as_pointer() for item in getattr(bpy.data, kind)}
            for kind in ownership.KINDS}


@contextmanager
def offline():
    from contextlib import ExitStack
    with ExitStack() as stack:
        for name in ('socket.create_connection', 'socket.socket.connect',
                     'socket.socket.connect_ex', 'urllib.request.urlopen'):
            stack.enter_context(patch(name, side_effect=AssertionError('Unexpected network request')))
        yield


class StreetProgressiveChecks(unittest.TestCase):
    def tearDown(self):
        if driver.active is not None:
            driver.active.close(silent=True)
        session = getattr(runtime.app, 'osm_session', None)
        if session is not None and not session.finished:
            session.cancel()
            session.finish(silent=True)

    def prepare(self, progressive=True, path=None):
        scene_runtime.clear_runtime()
        if bpy.context.object and bpy.context.object.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')
        for obj in list(bpy.data.objects):
            bpy.data.objects.remove(obj, do_unlink=True)
        for collection in list(bpy.data.collections):
            bpy.data.collections.remove(collection)
        for text in list(bpy.data.texts):
            # Retained GN assets reference their persisted PADL documents.
            # Remove test/report Texts only after their consumers are gone.
            if text.users == 0:
                bpy.data.texts.remove(text)
        bpy.data.orphans_purge(do_recursive=True)
        scene = bpy.context.scene
        scene.unit_settings.scale_length = 1.
        for name in ('lat', 'lon'):
            if name in scene:
                del scene[name]
        scene.property_unset('blosm')
        settings = scene.blosm
        values = dict(mode='3Drealistic', dataType='osm', osmSource='file',
            osmFilepath=str(path or FIXTURE), buildings=False, forests=False,
            vegetation=False, aviation=False, highways=True, railways=False,
            water=False, loadMissingMembers=False, setupScript='',
            progressiveImport=progressive, singleObject=False, commandLineMode=False,
            terrainObject='', assetLibraryMode='APPEND', streetBackend='python',
            streetCellSize=.7, streetResolution=1, streetExaggeration=1., streetProfileRise=.08)
        for name, value in values.items():
            setattr(settings, name, value)
        prefs = bpy.context.preferences.addons['blosm'].preferences
        prefs.enableExperimentalFeatures = True
        prefs.dataDir = str(ARGS.artifacts / 'data')
        prefs.preferPerformanceDuringImport = False
        for cls, initial in SOURCE_IDS.items():
            cls.ID = initial
        random.seed(1729)
        return settings

    def new_session(self):
        session = runtime.OsmImport(bpy.context)
        self.addCleanup(lambda: None if session.finished else
                        (session.cancel(), session.finish(silent=True)))
        return session

    def advance(self, session, predicate):
        for _ in range(200000):
            if predicate() or session.job.done:
                break
            session.job.step(max_steps=1)
        self.assertIsNone(session.job.error, str(session.job.error))
        self.assertTrue(predicate(), 'Required checkpoint was not reached: ' + session.phase)

    def drain(self, session, *, successful=True):
        trace = []
        for _ in range(200000):
            if session.job.done:
                break
            start = time.perf_counter()
            session.job.step(max_steps=1)
            trace.append(dict(phase=session.phase, committed=session.committed_units,
                              seconds=time.perf_counter() - start))
        self.assertTrue(session.job.done, 'Progressive session did not terminate')
        result = session.finish(silent=True)
        METRICS.append(dict(test=self.id(), result=result, callbacks=trace))
        self.assertFalse(session.ownership.pending)
        self.assertIsNone(session.finish(silent=True), 'Cleanup ran twice')
        if successful:
            self.assertIsNone(session.job.error, str(session.job.error))
            self.assertEqual(result['outcome'], 'success', result)
            self.assertFalse(result['errors'])
        return result, trace

    def first_component(self, session):
        self.advance(session, lambda: session.committed_units == 1)
        self.assertFalse(session.job.done, 'First component had no visible UI boundary')
        self.assertEqual(len(controllers(committed=True)), 1)
        return controllers(committed=True)[0]

    def sentinel(self):
        mesh = bpy.data.meshes.new('Unrelated user mesh')
        mesh.from_pydata([(0., 0., 0.), (1., 0., 0.), (0., 1., 0.)], [], [(0, 1, 2)])
        obj = bpy.data.objects.new('Unrelated user object', mesh)
        bpy.context.scene.collection.objects.link(obj)
        text = bpy.data.texts.new('Unrelated user document')
        text.write('Preserve this user edit between import turns.')
        obj['user_document'] = text
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
        return obj, text

    def test_actual_osm_adapter_exact_parity(self):
        self.prepare(False)
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        expected = scene_state()
        self.assertEqual(len(expected), 3)
        self.prepare(True)
        session = self.new_session()
        selected, _ = self.sentinel()
        result, trace = self.drain(session)
        self.assertEqual(result['committed'], 3)
        self.assertEqual(expected, scene_state())
        self.assertEqual(bpy.context.view_layer.objects.active, selected)
        self.assertTrue(selected.select_get())
        self.assertIn(1, [row['committed'] for row in trace])
        self.assertIn(2, [row['committed'] for row in trace])
        self.assertGreater(sum(row['committed'] == 0 for row in trace), 10,
                           'Street generation provided no internal checkpoints')
        self.assertIsNotNone(result['first_result_seconds'])

    def test_cancel_before_first_step(self):
        self.prepare()
        self.sentinel()
        before = inventory()
        session = self.new_session()
        session.cancel()
        result = session.finish(silent=True)
        self.assertEqual(result['outcome'], 'cancelled')
        self.assertEqual(result['committed'], 0)
        self.assertFalse(session.ownership.pending)
        self.assertEqual(inventory(), before)
        self.assertIsNone(session.finish(silent=True))

    def test_cancel_hidden_component_keeps_user_datablocks(self):
        self.prepare()
        session = self.new_session()
        self.advance(session, lambda: bool(controllers()) and not session.committed_units)
        sentinel, note = self.sentinel()
        session.cancel()
        result = session.finish(silent=True)
        self.assertEqual(result['outcome'], 'cancelled')
        self.assertEqual(result['committed'], 0)
        self.assertEqual(controllers(), [])
        self.assertTrue(ownership.alive(sentinel))
        self.assertTrue(ownership.alive(sentinel.data))
        self.assertTrue(ownership.alive(note))
        self.assertFalse(session.ownership.pending)
        self.assertFalse(scene_runtime._ENGINES, 'Cancelled component retained numerical cache')

    def test_cancel_after_commit_preserves_exact_result_and_user_edit(self):
        self.prepare()
        session = self.new_session()
        self.first_component(session)
        sentinel, note = self.sentinel()
        expected, before = scene_state(), inventory()
        session.cancel()
        result = session.finish(silent=True)
        self.assertEqual(result['outcome'], 'cancelled')
        self.assertEqual(result['committed'], 1)
        self.assertEqual(scene_state(), expected)
        self.assertEqual(inventory(), before)
        self.assertEqual(sentinel['user_document'], note)
        self.assertFalse(session.ownership.pending)

    def test_cancel_staging_reclaims_pending_meshes(self):
        self.prepare()
        session = self.new_session()
        self.first_component(session)
        expected = scene_state()
        sentinel, note = self.sentinel()
        before = inventory()
        staged = lambda: [obj for obj in bpy.data.objects if obj.get('st_py_key')
                          and obj.get(ownership.PENDING)
                          and (obj.get('st_py_preview') or not obj.users_collection)]
        self.advance(session, lambda: bool(staged()))
        self.assertEqual(session.committed_units, 1)
        detached = [obj.as_pointer() for obj in staged()]
        session.cancel()
        result = session.finish(silent=True)
        self.assertEqual(result['outcome'], 'cancelled')
        self.assertEqual(scene_state(), expected)
        self.assertEqual(inventory(), before)
        self.assertFalse(set(detached) & inventory()['objects'])
        self.assertEqual(sentinel['user_document'], note)
        self.assertFalse(session.ownership.pending)

    def test_precommit_mesh_failure_preserves_prior_and_user_data(self):
        self.prepare()
        session = self.new_session()
        self.first_component(session)
        self.sentinel()
        expected, before = scene_state(), inventory()
        original = scene_runtime.SceneWriter._mesh
        allocated = []

        def fail_after_mesh(*args, **kwargs):
            mesh = original(*args, **kwargs)
            allocated.append(mesh.as_pointer())
            raise TypeError('Injected precommit mesh failure')

        with patch.object(scene_runtime.SceneWriter, '_mesh', side_effect=fail_after_mesh):
            result, _ = self.drain(session, successful=False)
        self.assertTrue(allocated, 'Failure did not reach actual mesh allocation')
        self.assertEqual(result['outcome'], 'partial', result)
        self.assertEqual(result['committed'], 1)
        self.assertTrue(any('Injected precommit' in error for error in result['errors']))
        self.assertEqual(scene_state(), expected)
        self.assertEqual(inventory(), before, 'Failed component leaked or removed a datablock')

    def test_staged_source_edit_preserves_last_accepted_network(self):
        from _python_street_fixtures import fixture_networks
        from blosm.way.street_network.python_backend import model
        self.prepare()
        network = next(network for network in fixture_networks(model)
                       if network.uuid == 'bent_two_intersections')
        controller = scene_runtime.NetworkScene.create(network)
        bpy.context.view_layer.update()
        expected, manifest = controller_state(controller), controller['st_py_manifest']
        source = next(obj for obj in controller['st_py_sources'].objects
                      if len(obj.data.vertices) == 3)
        source.data.vertices[1].co.y += .1
        source.data.update()
        steps = scene_runtime.NetworkScene.update_steps(controller)
        try:
            for _ in range(200000):
                next(steps)
                if any(obj.get('st_py_key') and obj.get(ownership.PENDING)
                       and not obj.users_collection for obj in bpy.data.objects):
                    break
            else:
                self.fail('No staged output checkpoint')
            sentinel, note = self.sentinel()
            source.data.vertices[1].co.y += .1
            source.data.update()
            with self.assertRaisesRegex(RuntimeError, 'sources changed'):
                for _ in steps:
                    pass
        finally:
            steps.close()
        self.assertEqual(controller['st_py_manifest'], manifest)
        self.assertEqual(controller_state(controller), expected)
        self.assertEqual(sentinel['user_document'], note)
        self.assertFalse(any(obj.get(ownership.PENDING) for obj in bpy.data.objects))
        self.assertEqual(scene_runtime._ENGINES[network.uuid].manifest(), json.loads(manifest))
        # A fresh update consumes the user's current source edit and succeeds.
        scene_runtime.NetworkScene.update(controller)
        self.assertNotEqual(controller_state(controller)['outputs'], expected['outputs'])

    def test_failure_after_engine_commit_restores_last_accepted_network(self):
        from _python_street_fixtures import fixture_networks
        from blosm.way.street_network.python_backend import model
        self.prepare()
        network = next(network for network in fixture_networks(model)
                       if network.uuid == 'bent_two_intersections')
        controller = scene_runtime.NetworkScene.create(network)
        bpy.context.view_layer.update()
        expected, manifest = controller_state(controller), controller['st_py_manifest']
        identities = {obj['st_py_key']: obj.as_pointer()
                      for obj in controller['st_py_outputs'].objects}
        source = next(obj for obj in controller['st_py_sources'].objects if len(obj.data.vertices) == 3)
        source.data.vertices[1].co.y += .2
        source.data.update()
        sentinel, note = self.sentinel()
        before = inventory()
        engine = scene_runtime._ENGINES[network.uuid]
        original, committed_revisions = engine.commit, []

        def fail_after_commit(update):
            result = original(update)
            committed_revisions.append(engine.updates.revision)
            raise RuntimeError('Injected failure after real engine commit')

        with patch.object(engine, 'commit', side_effect=fail_after_commit):
            with self.assertRaisesRegex(RuntimeError, 'after real engine commit'):
                for _ in scene_runtime.NetworkScene.update_steps(controller):
                    pass
        self.assertEqual(committed_revisions, [json.loads(manifest)['revision'] + 1])
        self.assertEqual(controller['st_py_manifest'], manifest)
        self.assertEqual(controller_state(controller), expected)
        self.assertEqual(engine.manifest(), json.loads(manifest))
        self.assertEqual(inventory(), before)
        self.assertEqual(identities, {obj['st_py_key']: obj.as_pointer()
                                     for obj in controller['st_py_outputs'].objects})
        self.assertEqual(sentinel['user_document'], note)
        scene_runtime.NetworkScene.update(controller)
        self.assertNotEqual(controller_state(controller)['outputs'], expected['outputs'])

    def test_live_collection_unlink_invalidates_staged_import(self):
        for target in ('component', 'outputs', 'street_root'):
            with self.subTest(collection=target):
                self.prepare()
                session = self.new_session()
                first = self.first_component(session)
                expected = controller_state(first)
                staged = lambda: [obj for obj in bpy.data.objects if obj.get('st_py_key')
                                  and obj.get(ownership.PENDING)
                                  and (obj.get('st_py_preview') or not obj.users_collection)]
                self.advance(session, lambda: bool(staged()))
                pending_controller = staged()[0]['st_py_controller']
                component = pending_controller.users_collection[0]
                output_collection = pending_controller['st_py_outputs']
                if target in ('component', 'street_root'):
                    parent = next(collection for collection in bpy.data.collections
                                  if component.name in collection.children)
                    if target == 'component':
                        parent.children.unlink(component)
                    else:
                        bpy.context.scene.collection.children.unlink(parent)
                else:
                    component.children.unlink(output_collection)
                self.assertTrue(ownership.alive(component))
                self.assertTrue(ownership.alive(output_collection))
                sentinel, note = self.sentinel()
                result, _ = self.drain(session, successful=False)
                self.assertEqual(result['outcome'], 'partial')
                self.assertTrue(result['errors'], 'Unlinked collection was accepted without an error')
                if ownership.alive(pending_controller):
                    self.assertFalse(pending_controller.get('st_py_manifest'),
                                     'Unlinked component advanced its accepted manifest')
                if ownership.alive(output_collection):
                    self.assertEqual(tuple(output_collection.objects), ())
                self.assertEqual(controller_state(first), expected)
                self.assertEqual(sentinel['user_document'], note)
                self.assertFalse(staged())
                self.assertFalse(any(obj.get(ownership.PENDING) for obj in bpy.data.objects))

    def test_carriageway_material_override_is_captured(self):
        settings = self.prepare()
        material = bpy.data.materials.new('Captured pavement override')
        alternate = bpy.data.materials.new('Changed UI pavement override')
        settings.streetCarriagewayMaterial = material
        session = self.new_session()
        settings.streetCarriagewayMaterial = alternate
        self.drain(session)
        pavement_objects = []
        for controller in controllers(committed=True):
            self.assertEqual(controller['st_py_pavement_material'], material)
            for obj in controller['st_py_outputs'].objects:
                metadata = json.loads(obj.get('st_py_metadata', '{}'))
                if obj.type == 'MESH' and metadata.get('terrain_material') == 'pavement':
                    pavement_objects.append(obj)
                    self.assertEqual(tuple(obj.data.materials), (material,))
        self.assertTrue(pavement_objects, 'No carriageway material assignment exercised')
        self.assertEqual(settings.streetCarriagewayMaterial, alternate)

    def test_deleted_material_override_invalidates_import(self):
        settings = self.prepare()
        material = bpy.data.materials.new('Deleted captured pavement')
        settings.streetCarriagewayMaterial = material
        session = self.new_session()
        sentinel, note = self.sentinel()
        bpy.data.materials.remove(material)
        result, _ = self.drain(session, successful=False)
        self.assertEqual(result['outcome'], 'failed')
        self.assertEqual(result['committed'], 0)
        self.assertTrue(any('deleted' in error.lower() for error in result['errors']))
        self.assertEqual(sentinel['user_document'], note)
        self.assertEqual(controllers(), [])

    def test_mixed_native_buildings_and_streets_exact_parity(self):
        fixture = ARGS.artifacts / 'streets_and_buildings.osm'
        tree = ET.parse(FIXTURE)
        osm = tree.getroot()
        for number, x0 in enumerate((50., 170.)):
            ids = []
            for corner, (x, y) in enumerate(((x0, 60.), (x0 + 10., 60.),
                                             (x0 + 10., 70.), (x0, 70.))):
                identity = str(1000 + number * 10 + corner)
                ET.SubElement(osm, 'node', id=identity,
                              lat=str(52. + y / 111320.), lon=str(13. + x / 68535.))
                ids.append(identity)
            way = ET.SubElement(osm, 'way', id=str(2000 + number))
            for identity in (*ids, ids[0]):
                ET.SubElement(way, 'nd', ref=identity)
            for key, value in (('building', 'yes'), ('building:levels', '3'), ('roof:shape', 'flat')):
                ET.SubElement(way, 'tag', k=key, v=value)
        tree.write(fixture, encoding='utf-8', xml_declaration=True)
        setup = ROOT / 'tests/fixtures/progressive_streets_buildings_setup.py'
        for backend in ('geometry_nodes', 'python'):
            with self.subTest(buildings=backend):
                states = []
                for progressive in (False, True):
                    settings = self.prepare(progressive, fixture)
                    settings.mode = '3Dsimple'
                    settings.buildings = True
                    settings.buildingBackend = backend
                    settings.buildingSeed = 113
                    settings.setupScript = str(setup)
                    session = self.new_session()
                    outcome, _ = self.drain(session)
                    self.assertEqual(outcome['committed'], 5 if progressive else 4)
                    self.assertEqual(len(controllers(committed=True)), 3)
                    results = runtime.app.padlBuildingResults
                    self.assertEqual(len(results), 2)
                    self.assertTrue(all(not row.get('error') for result in results for row in result.statuses))
                    bpy.context.view_layer.update()
                    graph = bpy.context.evaluated_depsgraph_get()
                    buildings = []
                    for result in results:
                        obj = result.output.evaluated_get(graph)
                        self.assertGreater(len(obj.data.polygons), 0)
                        buildings.append(dict(matrix=tuple(tuple(row) for row in obj.matrix_world),
                                              mesh=mesh_state(obj.data)))
                    states.append(dict(streets=scene_state(), buildings=buildings))
                    self.assertFalse(any(obj.get(ownership.PENDING) for obj in bpy.data.objects))
                self.assertEqual(states[0], states[1])

    def test_scene_units_change_invalidates_pending_import(self):
        self.prepare()
        session = self.new_session()
        self.first_component(session)
        self.sentinel()
        before = inventory()
        bpy.context.scene.unit_settings.scale_length = .01
        result, _ = self.drain(session, successful=False)
        self.assertEqual(result['outcome'], 'partial')
        self.assertTrue(any('units changed' in error for error in result['errors']))
        self.assertEqual(inventory(), before)

    def test_scene_switch_pauses_without_scene_mutation(self):
        self.prepare()
        session = self.new_session()
        self.first_component(session)
        controller = driver.Controller(session, bpy.context.window)
        driver.active = controller
        original = bpy.context.window.scene
        other = bpy.data.scenes.new('Other user scene')
        try:
            bpy.context.window.scene = other
            before, revision = inventory(), session.revision
            for _ in range(3):
                controller.tick()
            self.assertTrue(controller.paused)
            self.assertEqual(inventory(), before)
            self.assertEqual(session.revision, revision)
            bpy.context.window.scene = original
            controller.tick()
            self.assertFalse(controller.paused)
            controller.close(silent=True)
            self.assertIsNone(driver.active)
            self.assertTrue(session.finished)
            self.assertEqual(runtime.last_result['outcome'], 'cancelled')
            self.assertFalse(session.ownership.pending)
        finally:
            bpy.context.window.scene = original
            if driver.active is controller:
                controller.close(silent=True)
            bpy.data.scenes.remove(other)

    def test_target_collection_deletion_stops_pending_import(self):
        self.prepare()
        session = self.new_session()
        self.first_component(session)
        sentinel, note = self.sentinel()
        expected = scene_state()
        bpy.data.collections.remove(session.destination)
        result, _ = self.drain(session, successful=False)
        self.assertEqual(result['outcome'], 'partial')
        self.assertTrue(any('collection' in error for error in result['errors']))
        self.assertTrue(ownership.alive(sentinel))
        self.assertEqual(sentinel['user_document'], note)
        self.assertEqual(scene_state(), expected)

    def test_real_extracts_exact_parity(self):
        for path in ARGS.osm:
            with self.subTest(osm=str(path)):
                self.prepare(False, path)
                self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
                expected = scene_state()
                self.assertTrue(expected, 'Real extract produced no accepted street network')
                self.prepare(True, path)
                result, _ = self.drain(self.new_session())
                self.assertEqual(scene_state(), expected)
                self.assertEqual(result['committed'], len(expected))
                METRICS.append(dict(osm=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                    exact_geometry_sha256=digest(expected), components=len(expected)))

    def test_fixture_matrix_exact_parity(self):
        from _python_street_fixtures import (fixture_networks, transition_fixture_networks,
            attached_overlap_fixture_networks, source_support_fixture_networks,
            compound_route_fixture_networks)
        from blosm.way.street_network.python_backend import model
        fixtures = fixture_networks(model)
        if ARGS.extended:
            fixtures += (transition_fixture_networks(model) + attached_overlap_fixture_networks(model)
                         + source_support_fixture_networks(model) + compound_route_fixture_networks(model))
        else:
            fixtures = tuple(network for network in fixtures if network.uuid in
                             ('T', 'unequal_widths', 'bent_two_intersections'))
        known_rejections = {
            ('seamed_approach', 0):
                'physical_intersection:10: Corner between Section 3 end 0 and Section 1 end 0, '
                'radius 2 m, available approaches 40.000/2.000 m: '
                'Requested-radius physical-curb contact did not converge',
        }
        for network in fixtures:
            for interpolation in (0, 3):
                with self.subTest(network=network.uuid, interpolation=interpolation):
                    self.prepare()
                    network = replace(network, settings=replace(network.settings, resolution=1),
                        sections=tuple(replace(section, interpolation=interpolation)
                                       for section in network.sections))
                    controller = scene_runtime.NetworkScene.create(network, generate=False)
                    rejection = None
                    try:
                        scene_runtime.NetworkScene.update(controller)
                    except ValueError as exc:
                        self.assertIn((network.uuid, interpolation), known_rejections,
                                      'An otherwise valid fixture unexpectedly failed')
                        self.assertEqual(str(exc), known_rejections[(network.uuid, interpolation)])
                        rejection = type(exc), str(exc)
                        self.assertFalse(controller.get('st_py_manifest'))
                        self.assertEqual(tuple(controller['st_py_outputs'].objects), ())
                        self.assertFalse(any(obj.get(ownership.PENDING) for obj in bpy.data.objects))
                    bpy.context.view_layer.update()
                    expected = controller_state(controller)
                    self.prepare()
                    controller = scene_runtime.NetworkScene.create(network, generate=False)
                    phases = []
                    steps = scene_runtime.NetworkScene.update_steps(controller, on_progress=phases.append)
                    callbacks = 0
                    if rejection is not None:
                        with self.assertRaises(rejection[0]) as failure:
                            for _ in steps:
                                callbacks += 1
                        self.assertEqual(str(failure.exception), rejection[1])
                        self.assertFalse(controller.get('st_py_manifest'))
                        self.assertEqual(tuple(controller['st_py_outputs'].objects), ())
                        self.assertFalse(any(obj.get(ownership.PENDING) for obj in bpy.data.objects))
                        self.assertNotIn(network.uuid, scene_runtime._ENGINES)
                        METRICS.append(dict(fixture=network.uuid, interpolation=interpolation,
                            callbacks=callbacks, phases=phases, rejection_parity=True,
                            error_type=rejection[0].__name__, message=rejection[1]))
                        continue
                    while True:
                        try:
                            next(steps)
                            callbacks += 1
                        except StopIteration:
                            break
                    bpy.context.view_layer.update()
                    self.assertEqual(controller_state(controller), expected)
                    self.assertGreater(callbacks, 3)
                    self.assertTrue(any('Staging' in phase for phase in phases))
                    self.assertFalse(any(obj.get(ownership.PENDING) for obj in bpy.data.objects))
                    METRICS.append(dict(fixture=network.uuid, interpolation=interpolation,
                                        callbacks=callbacks, phases=phases, exact_sha256=digest(expected)))

    def test_fixture_terrain_modes_exact_parity(self):
        from _python_street_fixtures import fixture_networks
        from check_street_terrain import make_terrain
        from blosm.way.street_network.python_backend import model
        network = next(network for network in fixture_networks(model)
                       if network.uuid == 'bent_two_intersections')
        for mode in ('PROFILE', 'DIRECT'):
            for scale in (1., .1):
                with self.subTest(terrain=mode, scale=scale):
                    states, timings = [], []
                    for progressive in (False, True):
                        self.prepare()
                        bpy.context.scene.unit_settings.scale_length = scale
                        controller = scene_runtime.NetworkScene.create(network, generate=False)
                        terrain = make_terrain('Progressive terrain', extent=150., scale=scale)
                        controller.blosm_street_terrain.target = terrain
                        controller.blosm_street_terrain.mode = mode
                        before = mesh_state(terrain.data)
                        started = time.perf_counter()
                        if progressive:
                            for _ in scene_runtime.NetworkScene.update_steps(controller):
                                pass
                        else:
                            scene_runtime.NetworkScene.update(controller)
                        bpy.context.view_layer.update()
                        timings.append(time.perf_counter() - started)
                        self.assertEqual(mesh_state(terrain.data), before, 'Imported streets changed terrain')
                        states.append(controller_state(controller))
                    self.assertEqual(states[0], states[1])
                    METRICS.append(dict(terrain=mode, scale=scale, seconds=timings,
                                        exact_sha256=digest(states[0])))


class ReportResult(unittest.TextTestResult):
    records = []

    def startTest(self, test):
        super().startTest(test)
        self.started = time.monotonic()

    def record(self, test, status, detail=None):
        self.records.append(dict(id=test.id(), status=status, detail=detail,
                                 seconds=time.monotonic() - self.started))

    def addSuccess(self, test):
        super().addSuccess(test)
        self.record(test, 'passed')

    def addFailure(self, test, err):
        super().addFailure(test, err)
        self.record(test, 'failed', self._exc_info_to_string(err, test))

    def addError(self, test, err):
        super().addError(test, err)
        self.record(test, 'error', self._exc_info_to_string(err, test))

    def addSubTest(self, test, subtest, err):
        super().addSubTest(test, subtest, err)
        if err is not None:
            self.record(subtest, 'failed' if issubclass(err[0], test.failureException) else 'error',
                        self._exc_info_to_string(err, test))


def main():
    global ARGS, FIXTURE, METRICS, runtime, driver, scene_runtime, ownership
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifacts', type=Path, required=True)
    parser.add_argument('--osm', type=Path, action='append', default=[])
    parser.add_argument('--only', action='append', default=[])
    parser.add_argument('--extended', action='store_true',
                        help='Run all authored junction, Bundle, boundary and internal-route fixtures')
    ARGS = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    ARGS.artifacts.mkdir(parents=True, exist_ok=True)
    (ARGS.artifacts / 'data').mkdir(exist_ok=True)
    FIXTURE = ARGS.artifacts / 'three_disconnected_networks.osm'
    synthetic_osm(FIXTURE)
    METRICS = []
    enable()
    from blosm.app import osm_import as runtime
    from blosm.progressive import blender as driver, ownership
    from blosm.way.street_network.python_backend import scene as scene_runtime
    names = unittest.defaultTestLoader.getTestCaseNames(StreetProgressiveChecks)
    if ARGS.only:
        names = [name for name in names if any(fragment in name for fragment in ARGS.only)]
        assert names, 'No selected test methods'
    if not ARGS.osm:
        names = [name for name in names if name != 'test_real_extracts_exact_parity']
    suite = unittest.TestSuite(StreetProgressiveChecks(name) for name in names)
    with offline():
        result = unittest.TextTestRunner(verbosity=2, resultclass=ReportResult).run(suite)
    report = dict(passed=result.wasSuccessful(), blender=bpy.app.version_string,
                  checkout=str(ROOT), cases=result.records, metrics=METRICS,
                  real_extracts=[str(path) for path in ARGS.osm])
    (ARGS.artifacts / 'checks.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('STREET_PROGRESSIVE_CHECKS', json.dumps(dict(passed=report['passed'], tests=result.testsRun)), flush=True)
    if not result.wasSuccessful():
        raise RuntimeError('Street progressive checks failed; see checks.json')


if __name__ == '__main__':
    main()
