"""Qualify visible early street previews in a factory-profile Blender process.

Reuses the real OSM import, exact geometry and ownership helpers from
check_street_progressive. All generated fixtures and reports use --artifacts.
All fixtures are local and never access the network.
"""
import argparse
from dataclasses import replace
import importlib
import json
from pathlib import Path
import sys
import time
import unittest
from unittest.mock import patch
from xml.etree import ElementTree as ET

import bpy

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_street_progressive as support


def connected_osm(path):
    """Two joined T junctions: every emitted street belongs to one component."""
    osm = ET.Element('osm', version='0.6')
    ET.SubElement(osm, 'bounds', minlat='51.999', minlon='12.999',
                  maxlat='52.002', maxlon='13.002')
    points = ((0., 0.), (20., 2.), (40., 0.), (-30., 0.),
              (0., 30.), (70., 0.), (40., -30.))
    for number, (x, y) in enumerate(points, 1):
        ET.SubElement(osm, 'node', id=str(number), lat=str(52. + y / 111320.),
                      lon=str(13. + x / 68535.))
    for number, refs in enumerate(((1, 2, 3), (1, 4), (1, 5), (3, 6), (3, 7)), 101):
        way = ET.SubElement(osm, 'way', id=str(number))
        for ref in refs:
            ET.SubElement(way, 'nd', ref=str(ref))
        for key, value in (('highway', 'residential'), ('lanes', '2'),
                           ('name', 'Connected preview fixture ' + str(number))):
            ET.SubElement(way, 'tag', k=key, v=value)
    ET.ElementTree(osm).write(path, encoding='utf-8', xml_declaration=True)


def previews(controller=None, *, visible=False):
    layer = bpy.context.view_layer
    return [obj for obj in bpy.data.objects
            if obj.get('st_py_preview') and obj.type == 'MESH' and len(obj.data.polygons)
            and (controller is None or obj.get('st_py_controller') == controller)
            and (not visible or obj.name in layer.objects and obj.visible_get(view_layer=layer))]


def pending():
    return [(kind, item) for kind in support.ownership.KINDS
            for item in getattr(bpy.data, kind) if item.get(support.ownership.PENDING)]


def fixture_network():
    from _python_street_fixtures import fixture_networks
    from blosm.way.street_network.python_backend import model
    return next(network for network in fixture_networks(model)
                if network.uuid == 'bent_two_intersections')


class StreetPreviewChecks(support.StreetProgressiveChecks):
    def observer(self, controller, phases):
        records, identities = [], {}
        started = time.perf_counter()

        def callback(count):
            bpy.context.view_layer.update()
            visible = previews(controller, visible=True)
            records.append(dict(count=count, visible=len(visible),
                seconds=time.perf_counter() - started, phase=phases[-1] if phases else None))
            if not count:
                return
            self.assertFalse(controller.get('st_py_manifest'), 'Preview accepted a canonical manifest')
            self.assertEqual(len(controller['st_py_outputs'].objects), 0,
                             'Preview entered the canonical output collection')
            self.assertEqual(len(visible), count, 'Reported preview is not actually visible')
            for obj in visible:
                self.assertTrue(obj.get(support.ownership.PENDING))
                self.assertTrue(obj.data.get(support.ownership.PENDING))
                self.assertTrue(obj.users_collection, 'Preview mesh is still detached')
                self.assertTrue(all(collection != controller['st_py_outputs']
                    and collection.get('st_py_preview') and collection.get(support.ownership.PENDING)
                    for collection in obj.users_collection))
                identities.setdefault(obj['st_py_key'], (obj.as_pointer(), obj.data.as_pointer()))
            if len([record for record in records if record['count']]) == 1:
                self.assertNotIn('Street geometry prepared', phases)
                self.assertNotIn('Validating street overlaps', phases)
                self.assertNotIn('Checking stored street overlaps (complete component)', phases)
        return callback, records, identities

    def assert_no_previews(self):
        self.assertFalse(previews())
        self.assertFalse(any(collection.get('st_py_preview') for collection in bpy.data.collections))
        self.assertFalse(pending(), 'Preview cleanup left pending datablocks')

    def test_connected_initial_preview_exact_parity_and_reuse(self):
        from check_street_terrain import make_terrain
        network = fixture_network()
        # Both junctions share Section 1 and cover every authored Section.
        ends = [{end.section_id for end in junction.ends} for junction in network.intersections]
        self.assertTrue(set.intersection(*ends))
        self.assertEqual(set.union(*ends), {section.id for section in network.sections})
        for mode in (None, 'PROFILE', 'DIRECT'):
            with self.subTest(terrain=mode or 'NONE'):
                states = []
                for progressive in (False, True):
                    self.prepare()
                    controller = support.scene_runtime.NetworkScene.create(network, generate=False)
                    terrain = None
                    if mode:
                        terrain = make_terrain('Preview parity terrain', extent=150.)
                        controller.blosm_street_terrain.target = terrain
                        controller.blosm_street_terrain.mode = mode
                        terrain_state = support.mesh_state(terrain.data)
                    if progressive:
                        phases = []
                        callback, records, identities = self.observer(controller, phases)
                        for _ in support.scene_runtime.NetworkScene.update_steps(
                                controller, on_progress=phases.append, on_preview=callback):
                            pass
                        positive = [record for record in records if record['count']]
                        self.assertGreater(len({record['count'] for record in positive}), 1,
                                           'One connected component had no progressive preview growth')
                        self.assertEqual(records[-1]['count'], 0)
                        outputs = {obj['st_py_key']: obj for obj in controller['st_py_outputs'].objects}
                        for key, pair in identities.items():
                            self.assertIn(key, outputs)
                            self.assertEqual((outputs[key].as_pointer(), outputs[key].data.as_pointer()), pair,
                                             'Final publication replaced a preview object or mesh: ' + key)
                        self.assert_no_previews()
                        support.METRICS.append(dict(test=self.id(), terrain=mode,
                            preview_records=records, reused_preview_objects=len(identities)))
                    else:
                        support.scene_runtime.NetworkScene.update(controller)
                    bpy.context.view_layer.update()
                    if terrain:
                        self.assertEqual(support.mesh_state(terrain.data), terrain_state)
                    states.append(support.controller_state(controller))
                self.assertEqual(states[0], states[1], 'Preview changed the exact final street result')
                support.METRICS[-1]['exact_geometry_sha256'] = support.digest(states[0])

    def test_initial_preview_requires_opt_in_and_existing_update_never_previews(self):
        self.prepare()
        controller = support.scene_runtime.NetworkScene.create(fixture_network(), generate=False)
        for _ in support.scene_runtime.NetworkScene.update_steps(controller):
            self.assertFalse(previews())
        expected = support.controller_state(controller)
        observed = []
        for _ in support.scene_runtime.NetworkScene.update_steps(controller, force=True,
                                                               on_preview=observed.append):
            self.assertFalse(previews(), 'An existing accepted network exposed unvalidated previews')
        self.assertFalse(any(observed), observed)
        # A forced regeneration increments the revision; geometry remains exact.
        self.assertEqual(support.controller_state(controller)['outputs'], expected['outputs'])
        self.assert_no_previews()

    def test_close_during_first_preview_removes_only_pending_ids(self):
        self.prepare()
        controller = support.scene_runtime.NetworkScene.create(fixture_network(), generate=False)
        sentinel, note = self.sentinel()
        before = support.inventory()
        phases = []
        callback, records, _ = self.observer(controller, phases)
        steps = support.scene_runtime.NetworkScene.update_steps(controller,
            on_progress=phases.append, on_preview=callback)
        try:
            for _ in steps:
                if previews(controller, visible=True):
                    break
            self.assertTrue(previews(controller, visible=True))
        finally:
            steps.close()
        self.assertEqual(records[-1]['count'], 0)
        self.assertEqual(support.inventory(), before)
        self.assertEqual(sentinel['user_document'], note)
        self.assertFalse(controller.get('st_py_manifest'))
        self.assertEqual(len(controller['st_py_outputs'].objects), 0)
        self.assert_no_previews()
        self.assertNotIn(fixture_network().uuid, support.scene_runtime._ENGINES)

    def test_source_edit_during_early_preview_invalidates_new_network(self):
        self.prepare()
        controller = support.scene_runtime.NetworkScene.create(fixture_network(), generate=False)
        sentinel, note = self.sentinel()
        before = support.inventory()
        phases = []
        callback, records, _ = self.observer(controller, phases)
        steps = support.scene_runtime.NetworkScene.update_steps(controller,
            on_progress=phases.append, on_preview=callback)
        try:
            for _ in steps:
                if previews(controller, visible=True):
                    break
            self.assertTrue(previews(controller, visible=True))
            source = next(obj for obj in controller['st_py_sources'].objects if len(obj.data.vertices) == 3)
            source.data.vertices[1].co.y += .1
            source.data.update()
            with self.assertRaisesRegex(RuntimeError, 'sources changed'):
                for _ in steps:
                    pass
            self.assertAlmostEqual(source.data.vertices[1].co.y, 2.1, places=5)
        finally:
            steps.close()
        self.assertEqual(records[-1]['count'], 0)
        self.assertEqual(support.inventory(), before)
        self.assertEqual(sentinel['user_document'], note)
        self.assertFalse(controller.get('st_py_manifest'))
        self.assert_no_previews()

    def test_early_preview_deleted_unlinked_moved_or_hidden_aborts_cleanly(self):
        cases = ('delete_object', 'delete_mesh', 'unlink_object', 'unlink_collection',
                 'delete_collection', 'move_object', 'hide_object', 'add_modifier')
        for mutation in cases:
            with self.subTest(mutation=mutation):
                self.prepare()
                network = fixture_network()
                controller = support.scene_runtime.NetworkScene.create(network, generate=False)
                sentinel, note = self.sentinel()
                before, counts = support.inventory(), []
                steps = support.scene_runtime.NetworkScene.update_steps(controller, on_preview=counts.append)
                try:
                    for _ in steps:
                        if counts and counts[-1] > 0:
                            break
                    objects = previews(controller)
                    self.assertTrue(objects, 'No early preview checkpoint reached')
                    self.assertFalse(controller.get('st_py_manifest'))
                    obj = objects[0]
                    collection = obj.users_collection[0]
                    if mutation == 'delete_object':
                        bpy.data.objects.remove(obj, do_unlink=True)
                    elif mutation == 'delete_mesh':
                        bpy.data.meshes.remove(obj.data, do_unlink=True)
                    elif mutation == 'unlink_object':
                        collection.objects.unlink(obj)
                    elif mutation == 'unlink_collection':
                        controller.users_collection[0].children.unlink(collection)
                    elif mutation == 'delete_collection':
                        bpy.data.collections.remove(collection, do_unlink=True)
                    elif mutation == 'move_object':
                        obj.location.x += .25
                    elif mutation == 'hide_object':
                        obj.hide_viewport = True
                    elif mutation == 'add_modifier':
                        obj.modifiers.new('User change to pending surface', 'BEVEL')
                    with self.assertRaisesRegex(support.scene_runtime.StreetInputChanged,
                                                'Pending street (preview|output)'):
                        for _ in steps:
                            pass
                finally:
                    steps.close()
                self.assertEqual(counts[-1], 0)
                self.assertEqual(support.inventory(), before, mutation)
                self.assertEqual(sentinel['user_document'], note)
                self.assertEqual(bpy.context.view_layer.objects.active, sentinel)
                self.assertFalse(controller.get('st_py_manifest'))
                self.assertEqual(len(controller['st_py_outputs'].objects), 0)
                self.assertNotIn(network.uuid, support.scene_runtime._ENGINES)
                self.assert_no_previews()

    def test_post_refinement_preview_buffers_cannot_bypass_final_validation(self):
        from blosm.way.street_network.python_backend.terrain.blender import RAW_Z
        cases = {'coordinates': 'geometry or materials changed',
                 'topology': 'geometry or materials changed',
                 'materials': 'geometry or materials changed',
                 'uv': 'UV coordinates changed',
                 'raw_z': 'terrain coordinates changed'}
        for mutation, message in cases.items():
            with self.subTest(mutation=mutation):
                self.prepare()
                controller = support.scene_runtime.NetworkScene.create(fixture_network(), generate=False)
                sentinel, note = self.sentinel()
                user_material = bpy.data.materials.new('User material retained across preview failure')
                before, counts, phases = support.inventory(), [], []
                steps = support.scene_runtime.NetworkScene.update_steps(controller,
                    on_progress=phases.append, on_preview=counts.append)
                objects = []
                try:
                    for _ in steps:
                        objects = previews(controller)
                        if objects and all(obj.get('st_py_metadata') and obj.data.attributes.get(RAW_Z)
                                           for obj in objects):
                            break
                    self.assertTrue(objects, 'No refined preview checkpoint reached')
                    self.assertTrue(all(obj.data.attributes.get(RAW_Z) for obj in objects))
                    self.assertTrue(any(phase.startswith('Staging street object ') for phase in phases))
                    self.assertNotIn('Committing street component', phases,
                                     'Mutation was too late to exercise a real UI checkpoint')
                    self.assertFalse(controller.get('st_py_manifest'))
                    self.assertEqual(len(controller['st_py_outputs'].objects), 0)
                    obj = objects[0]
                    mesh = obj.data
                    if mutation == 'coordinates':
                        mesh.vertices[0].co.z += .02
                        mesh.update()
                    elif mutation == 'topology':
                        values = support.np.empty(len(mesh.loops), dtype=support.np.int32)
                        mesh.loops.foreach_get('vertex_index', values)
                        values[0] = (int(values[0]) + 1) % len(mesh.vertices)
                        mesh.loops.foreach_set('vertex_index', values)
                        mesh.update()
                    elif mutation == 'materials':
                        mesh.materials.append(user_material)
                    elif mutation == 'uv':
                        mesh.uv_layers[0].uv[0].vector.x += .125
                    elif mutation == 'raw_z':
                        mesh.attributes[RAW_Z].data[0].value += .02
                    with self.assertRaisesRegex(support.scene_runtime.StreetInputChanged, message):
                        for _ in steps:
                            pass
                finally:
                    steps.close()
                self.assertEqual(counts[-1], 0)
                self.assertEqual(support.inventory(), before, mutation)
                self.assertEqual(sentinel['user_document'], note)
                self.assertTrue(support.ownership.alive(user_material))
                self.assertFalse(controller.get('st_py_manifest'))
                self.assertEqual(len(controller['st_py_outputs'].objects), 0)
                self.assertFalse(support.ownership.alive(obj))
                self.assertFalse(support.ownership.alive(mesh))
                self.assert_no_previews()

    def test_failure_after_preview_promotion_restores_pending_state_and_prior_result(self):
        preview_runtime = importlib.import_module('blosm.way.street_network.python_backend.previews')
        self.prepare()
        accepted = support.scene_runtime.NetworkScene.create(fixture_network())
        expected = support.controller_state(accepted)
        network = replace(fixture_network(), uuid='pending_promotion_failure')
        controller = support.scene_runtime.NetworkScene.create(network, generate=False)
        sentinel, note = self.sentinel()
        before, counts, promoted = support.inventory(), [], []
        original = preview_runtime.PendingPreviews.promote

        def fail_after_promotion(preview):
            original(preview)
            self.assertIsNone(preview.collection)
            self.assertTrue(controller.get('st_py_manifest'))
            self.assertTrue(controller['st_py_outputs'].objects)
            self.assertFalse(any(obj.get('st_py_preview') for obj in preview.objects.values()))
            promoted.extend((obj, obj.data) for obj in preview.objects.values())
            raise RuntimeError('Injected failure after preview promotion')

        with patch.object(preview_runtime.PendingPreviews, 'promote', fail_after_promotion):
            with self.assertRaisesRegex(RuntimeError, 'Injected failure after preview promotion'):
                for _ in support.scene_runtime.NetworkScene.update_steps(controller, on_preview=counts.append):
                    pass
        self.assertTrue(promoted, 'Failure hook did not reach actual preview promotion')
        self.assertEqual(counts[-1], 0)
        self.assertEqual(support.controller_state(accepted), expected)
        self.assertEqual(support.inventory(), before)
        self.assertEqual(sentinel['user_document'], note)
        self.assertFalse(controller.get('st_py_manifest'))
        self.assertEqual(len(controller['st_py_outputs'].objects), 0)
        self.assertNotIn(network.uuid, support.scene_runtime._ENGINES)
        self.assertFalse(any(support.ownership.alive(item) for pair in promoted for item in pair))
        self.assert_no_previews()

    def test_production_connected_import_previews_before_commit_and_exact_parity(self):
        path = support.ARGS.artifacts / 'connected_network.osm'
        connected_osm(path)
        self.prepare(False, path=path)
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        expected = support.scene_state()
        self.assertEqual(len(expected), 1, 'Fixture unexpectedly split into disconnected components')
        self.prepare(True, path=path)
        session = self.new_session()
        records = []
        began = time.perf_counter()
        for _ in range(200000):
            if session.job.done:
                break
            session.job.step(max_steps=1)
            bpy.context.view_layer.update()
            count = len(previews(visible=True))
            if count and (not records or records[-1]['count'] != count):
                self.assertEqual(session.committed_units, 0)
                self.assertEqual(support.controllers(committed=True), [])
                records.append(dict(count=count, phase=session.phase, seconds=time.perf_counter() - began))
        self.assertTrue(session.job.done)
        self.assertIsNone(session.job.error, str(session.job.error))
        result = session.finish(silent=True)
        self.assertEqual(result['outcome'], 'success', result)
        self.assertGreater(len(records), 1, 'Single-component production import never grew visible previews')
        self.assertEqual(result['committed'], 1)
        self.assertIsNotNone(result['first_preview_seconds'])
        self.assertLess(result['first_preview_seconds'], result['first_result_seconds'])
        self.assertGreater(result['preview_batches'], 1)
        self.assertEqual(session.preview_count, 0)
        self.assertEqual(support.scene_state(), expected)
        self.assert_no_previews()
        support.METRICS.append(dict(test=self.id(), result=result, previews=records,
                                    exact_geometry_sha256=support.digest(expected)))

    def test_production_cancel_early_preview_preserves_prior_component_and_user_edit(self):
        self.prepare()
        session = self.new_session()
        self.first_component(session)
        sentinel, note = self.sentinel()
        expected, before = support.scene_state(), support.inventory()
        self.advance(session, lambda: bool(previews(visible=True)))
        self.assertEqual(session.committed_units, 1)
        preview_pointers = {obj.as_pointer() for obj in previews()}
        session.cancel()
        result = session.finish(silent=True)
        self.assertEqual(result['outcome'], 'cancelled', result)
        self.assertEqual(result['committed'], 1)
        self.assertEqual(support.scene_state(), expected)
        self.assertEqual(support.inventory(), before)
        self.assertFalse(preview_pointers & support.inventory()['objects'])
        self.assertEqual(sentinel['user_document'], note)
        self.assertEqual(bpy.context.view_layer.objects.active, sentinel)
        self.assertTrue(sentinel.select_get())
        self.assertFalse(session.ownership.pending)
        self.assert_no_previews()

    def test_late_validation_failure_rolls_back_preview_preserves_prior_result(self):
        backend = importlib.import_module('blosm.way.street_network.python_backend.backend')
        self.prepare()
        session = self.new_session()
        self.first_component(session)
        sentinel, note = self.sentinel()
        expected, before = support.scene_state(), support.inventory()
        seen = []

        def reject(*args, **kwargs):
            self.assertTrue(previews(visible=True), 'Failure did not occur after real visible preview')
            seen.append(dict(objects=tuple(previews()), meshes=tuple(obj.data for obj in previews())))
            raise ValueError('Injected late complete-component overlap rejection')
            yield  # Keep the real numerical validation generator contract.

        with patch.object(backend, 'validate_overlaps_steps', reject):
            result, _ = self.drain(session, successful=False)
        self.assertTrue(seen)
        self.assertEqual(result['outcome'], 'partial', result)
        self.assertEqual(result['committed'], 1)
        self.assertTrue(any('Injected late complete-component' in error for error in result['errors']))
        self.assertEqual(support.scene_state(), expected)
        # Recoverable component failures retain only their authored controls and
        # empty canonical collections for diagnostics. They must not retain any
        # provisional product, nor accept a manifest for the failed component.
        failed = [controller for controller in support.controllers()
                  if not controller.get('st_py_manifest')]
        self.assertEqual(len(failed), 2)
        allowed = {kind: set() for kind in support.ownership.KINDS}
        for controller in failed:
            self.assertIn('Injected late complete-component', controller['st_py_error'])
            self.assertEqual(len(controller['st_py_outputs'].objects), 0)
            allowed['objects'].add(controller.as_pointer())
            allowed['collections'].update(collection.as_pointer() for collection in controller.users_collection)
            for role in ('sources', 'outputs'):
                allowed['collections'].add(controller['st_py_' + role].as_pointer())
            for source in controller['st_py_sources'].objects:
                allowed['objects'].add(source.as_pointer())
                allowed['meshes'].add(source.data.as_pointer())
        actual = support.inventory()
        for kind, pointers in before.items():
            self.assertEqual(actual[kind], pointers | allowed[kind],
                             'Failed diagnostics leaked or removed ' + kind)
        for record in seen:
            self.assertFalse(any(support.ownership.alive(obj) for obj in record['objects']))
            self.assertFalse(any(support.ownership.alive(mesh) for mesh in record['meshes']))
        self.assertEqual(sentinel['user_document'], note)
        self.assert_no_previews()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifacts', type=Path, required=True)
    parser.add_argument('--only', action='append', default=[])
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.artifacts.mkdir(parents=True, exist_ok=True)
    (args.artifacts / 'data').mkdir(exist_ok=True)
    support.ARGS = args
    support.FIXTURE = args.artifacts / 'three_disconnected_networks.osm'
    support.synthetic_osm(support.FIXTURE)
    support.METRICS = []
    support.enable()
    from blosm.app import osm_import as runtime
    from blosm.progressive import blender as driver, ownership
    from blosm.way.street_network.python_backend import scene as scene_runtime
    support.runtime, support.driver = runtime, driver
    support.ownership, support.scene_runtime = ownership, scene_runtime
    names = sorted(name for name in StreetPreviewChecks.__dict__ if name.startswith('test_'))
    if args.only:
        names = [name for name in names if any(fragment in name for fragment in args.only)]
        assert names, 'No selected preview test methods'
    support.ReportResult.records = []
    suite = unittest.TestSuite(StreetPreviewChecks(name) for name in names)
    with support.offline():
        result = unittest.TextTestRunner(verbosity=2, resultclass=support.ReportResult).run(suite)
    report = dict(passed=result.wasSuccessful(), blender=bpy.app.version_string,
                  checkout=str(support.ROOT), cases=result.records, metrics=support.METRICS)
    (args.artifacts / 'preview_checks.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('STREET_PREVIEW_CHECKS', json.dumps(dict(passed=report['passed'], tests=result.testsRun)), flush=True)
    if not result.wasSuccessful():
        raise RuntimeError('Street preview checks failed; see preview_checks.json')


if __name__ == '__main__':
    main()
