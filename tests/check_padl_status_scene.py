"""Private progressive GN status graph: parity, reuse and failure cleanup."""
from pathlib import Path
import importlib
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bpy
from release_support import arguments, enable, reset, run, PACKAGE

args = arguments()
enable(args)
backends = importlib.import_module(PACKAGE + '.building.padl_integration.backends')
gn = importlib.import_module(PACKAGE + '.building.padl_integration.gn_runtime')
ownership = importlib.import_module(PACKAGE + '.progressive.ownership')
Footprint = importlib.import_module(PACKAGE + '.padl.geometry').Footprint
PadlError = importlib.import_module(PACKAGE + '.padl.errors').PadlError


class StatusSceneChecks(unittest.TestCase):
    def setUp(self):
        reset(args)
        self.service = backends.BuildingGenerationService(backend='geometry_nodes')

    def generate(self, identifier=1, *, terrain=None):
        x = (identifier - 1) * 20
        footprint = Footprint('status-scene:' + str(identifier), identifier,
                              ((x, 0), (x + 12, 0), (x + 12, 8), (x, 8)), (),
                              {'building_use': 'house', 'roof_type': 'flat',
                               'max_elevation': 14.0, 'roof_height': 4.0,
                               'levels_wall_count': 3})
        return self.service.generate([footprint], collection=bpy.context.collection,
                                     terrain=terrain, seed=71)

    def resources(self, reader):
        return (reader.scene, reader.owner, reader.empty_mesh, reader.wrapper)

    def assert_live_private(self, reader, main_scene):
        for item in self.resources(reader):
            self.assertTrue(ownership.alive(item))
            self.assertTrue(item.get(ownership.TRANSIENT))
        self.assertEqual(bpy.context.scene, main_scene)
        self.assertNotEqual(reader.scene, main_scene)
        self.assertNotIn(reader.owner.name, main_scene.objects)
        self.assertIn(reader.owner.name, reader.scene.objects)
        self.assertEqual(reader.owner.data, reader.empty_mesh)

    def assert_released(self, resources):
        self.assertTrue(all(not ownership.alive(item) for item in resources))
        self.assertFalse(any(item.get(ownership.TRANSIENT)
                             for kind in ('scenes', 'objects', 'meshes', 'node_groups')
                             for item in getattr(bpy.data, kind)))

    def test_reuse_status_parity_and_borrowed_mesh_release(self):
        main_scene = bpy.context.scene
        with self.service.batch(progressive=True):
            reader = self.service.backend._batch_reader
            resources = self.resources(reader)
            pointers = tuple(item.as_pointer() for item in resources)
            tree = self.service.backend._batch_tree
            results = []
            for identifier in (1, 2):
                result = self.generate(identifier)
                results.append(result)
                self.assert_live_private(reader, main_scene)
                self.assertEqual(tuple(item.as_pointer() for item in self.resources(reader)), pointers)
                self.assertEqual(result.owner.data.users, 1)
                self.assertTrue(result.statuses)
                self.assertTrue(all(not row['error'] for row in result.statuses))
                expected = gn.read_status(result.owner, result.owner.modifiers['PADL Buildings'])
                self.assertEqual(result.statuses, expected)
                self.assertNotIn(ownership.TRANSIENT, result.owner.data)
        self.assert_released(resources)
        self.assertTrue(ownership.alive(tree))
        for result in results:
            self.assertTrue(ownership.alive(result.owner))
            self.assertTrue(ownership.alive(result.owner.data))

    def test_failed_read_releases_borrowed_mesh_and_batch_helpers(self):
        with self.assertRaisesRegex(RuntimeError, 'injected status extraction'):
            with self.service.batch(progressive=True):
                reader = self.service.backend._batch_reader
                resources = self.resources(reader)
                result = self.generate()
                try:
                    with patch.object(reader, '_statuses', side_effect=RuntimeError('injected status extraction')):
                        reader.read(result.owner, result.owner.modifiers['PADL Buildings'])
                finally:
                    self.assertEqual(reader.owner.data, reader.empty_mesh)
                    self.assertEqual(result.owner.data.users, 1)
        self.assert_released(resources)
        self.assertTrue(ownership.alive(result.owner))
        self.assertTrue(ownership.alive(result.owner.modifiers['PADL Buildings'].node_group))

    def test_deleted_private_resources_fail_cleanly(self):
        for kind, attribute in (('scenes', 'scene'), ('objects', 'owner'),
                                ('meshes', 'empty_mesh'), ('node_groups', 'wrapper')):
            with self.subTest(resource=kind):
                with self.service.batch(progressive=True):
                    reader = self.service.backend._batch_reader
                    resources = self.resources(reader)
                    result = self.generate()
                    getattr(bpy.data, kind).remove(getattr(reader, attribute), do_unlink=True)
                    with self.assertRaisesRegex(PadlError, 'status evaluator was deleted or changed'):
                        reader.read(result.owner, result.owner.modifiers['PADL Buildings'])
                self.assert_released(resources)
                self.assertTrue(ownership.alive(result.owner))
                self.assertTrue(ownership.alive(result.owner.data))

    def test_deleted_private_modifier_fails_cleanly(self):
        with self.service.batch(progressive=True):
            reader = self.service.backend._batch_reader
            resources = self.resources(reader)
            result = self.generate()
            reader.owner.modifiers.remove(reader.modifier)
            with self.assertRaisesRegex(PadlError, 'status evaluator was deleted or changed'):
                reader.read(result.owner, result.owner.modifiers['PADL Buildings'])
        self.assert_released(resources)

    def test_terrain_dependencies_remain_canonical(self):
        mesh = bpy.data.meshes.new('Status test terrain')
        mesh.from_pydata([(-10, -10, 0), (40, -10, 0), (40, 30, 0), (-10, 30, 0)], [], [(0, 1, 2, 3)])
        terrain = bpy.data.objects.new('Status test terrain', mesh)
        bpy.context.collection.objects.link(terrain)
        terrain.location.z = 36
        bpy.context.view_layer.update()
        with self.service.batch(progressive=True):
            reader = self.service.backend._batch_reader
            resources = self.resources(reader)
            # The second footprint has no terrain coverage; its diagnostic
            # must match the legacy reader as well as the successful case.
            for identifier in (1, 5):
                result = self.generate(identifier, terrain=terrain)
                modifier = result.owner.modifiers['PADL Buildings']
                self.assertEqual(result.statuses, gn.read_status(result.owner, modifier))
                self.assertEqual(bool(result.statuses[0]['error']), identifier == 5)
                self.assertEqual(reader.owner.data, reader.empty_mesh)
            helpers = [obj for obj in reader.scene.objects if obj not in (reader.owner, terrain)]
            self.assertEqual(len(helpers), 1)
            self.assertNotIn(ownership.TRANSIENT, terrain)
            self.assertNotIn(ownership.TRANSIENT, helpers[0])
        self.assert_released(resources)
        self.assertTrue(ownership.alive(terrain))
        self.assertTrue(ownership.alive(mesh))
        self.assertTrue(ownership.alive(helpers[0]))
        self.assertIn(helpers[0].name, bpy.context.scene.objects)

    def test_private_graph_notifications_do_not_invalidate_assets(self):
        with self.service.batch(progressive=True):
            reader = self.service.backend._batch_reader
            guard = self.service.backend._batch_guard
            root = self.service.backend._batch_tree
            update = SimpleNamespace(id=root)
            guard._updated(reader.scene, SimpleNamespace(scene=reader.scene, updates=[update]))
            self.assertFalse(guard.dirty)
            guard._updated(bpy.context.scene, SimpleNamespace(scene=bpy.context.scene, updates=[update]))
            self.assertTrue(guard.dirty)
            guard.dirty = False
            other_scene = bpy.data.scenes.new('Status test user scene')
            try:
                guard._updated(other_scene, SimpleNamespace(scene=other_scene, updates=[update]))
                self.assertTrue(guard.dirty, 'Shared generator edits in another user scene were ignored')
            finally:
                bpy.data.scenes.remove(other_scene)

    def test_animated_terrain_uses_import_scene_frame_and_fps(self):
        source = bpy.context.scene
        original_time = (source.frame_current, source.frame_subframe,
                         source.render.fps, source.render.fps_base,
                         source.render.frame_map_old, source.render.frame_map_new)
        mesh = bpy.data.meshes.new('Time-dependent terrain')
        mesh.from_pydata([(-1.5, -1, 0), (13, -1, 0), (13, 9, 0), (-1.5, 9, 0)], [], [(0, 1, 2, 3)])
        terrain = bpy.data.objects.new('Time-dependent terrain', mesh)
        bpy.context.collection.objects.link(terrain)
        tree = bpy.data.node_groups.new('Status test terrain clock', 'GeometryNodeTree')
        tree.interface.new_socket(name='Geometry', in_out='INPUT', socket_type='NodeSocketGeometry')
        tree.interface.new_socket(name='Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
        incoming = tree.nodes.new('NodeGroupInput')
        outgoing = tree.nodes.new('NodeGroupOutput')
        clock = tree.nodes.new('GeometryNodeInputSceneTime')
        vector = tree.nodes.new('ShaderNodeCombineXYZ')
        transform = tree.nodes.new('GeometryNodeTransform')
        tree.links.new(incoming.outputs['Geometry'], transform.inputs['Geometry'])
        tree.links.new(clock.outputs['Seconds'], vector.inputs['X'])
        tree.links.new(vector.outputs['Vector'], transform.inputs['Translation'])
        tree.links.new(transform.outputs['Geometry'], outgoing.inputs['Geometry'])
        terrain.modifiers.new('Animated terrain', 'NODES').node_group = tree
        try:
            with self.service.batch(progressive=True):
                reader = self.service.backend._batch_reader
                resources = self.resources(reader)
                # Translation by seconds alternates covered/uncovered without
                # changing the source mesh or modifier inputs. Exercise FPS
                # alone, frame alone, and FPS base with a fractional frame.
                for frame, fps, fps_base, error in ((12, 6, 1.0, True),
                                                   (12, 24, 1.0, False),
                                                   (240, 24, 1.0, True),
                                                   (12, 24, 3.0, True)):
                    with self.subTest(frame=frame, fps=fps, fps_base=fps_base):
                        source.render.fps = fps
                        source.render.fps_base = fps_base
                        source.frame_set(frame, subframe=0.25)
                        result = self.generate(terrain=terrain)
                        self.assertEqual(bool(result.statuses[0]['error']), error)
                        self.assertEqual(result.statuses, gn.read_status(
                            result.owner, result.owner.modifiers['PADL Buildings']))
                        self.assertEqual(reader.scene.frame_current, frame)
                        self.assertEqual(reader.scene.frame_subframe, 0.25)
                        self.assertEqual(reader.scene.render.fps, fps)
                        self.assertEqual(reader.scene.render.fps_base, fps_base)
                        self.assertEqual(bpy.context.scene, source)
                source.render.fps = 10
                source.render.fps_base = 1.0
                remapped_errors = set()
                for old, new in ((100, 200), (200, 100)):
                    source.render.frame_map_old = old
                    source.render.frame_map_new = new
                    source.frame_set(12, subframe=0.25)
                    result = self.generate(terrain=terrain)
                    remapped_errors.add(bool(result.statuses[0]['error']))
                    self.assertEqual(result.statuses, gn.read_status(
                        result.owner, result.owner.modifiers['PADL Buildings']))
                    self.assertEqual(reader.scene.frame_current_final, source.frame_current_final)
                self.assertEqual(remapped_errors, {False, True}, 'Time remapping did not affect terrain coverage')
            self.assert_released(resources)
        finally:
            frame, subframe, fps, fps_base, old, new = original_time
            source.render.fps = fps
            source.render.fps_base = fps_base
            source.render.frame_map_old = old
            source.render.frame_map_new = new
            source.frame_set(frame, subframe=subframe)

    def test_subdivision_terrain_uses_import_scene_simplification(self):
        source = bpy.context.scene
        original = source.render.use_simplify, source.render.simplify_subdivision
        mesh = bpy.data.meshes.new('Curved subdivision terrain')
        mesh.from_pydata([(-1, -1, 0), (13, -1, 2), (13, 9, 0), (-1, 9, 2)], [], [(0, 1, 2, 3)])
        terrain = bpy.data.objects.new('Curved subdivision terrain', mesh)
        bpy.context.collection.objects.link(terrain)
        modifier = terrain.modifiers.new('Terrain subdivision', 'SUBSURF')
        modifier.levels = 2
        try:
            with self.service.batch(progressive=True):
                reader = self.service.backend._batch_reader
                resources = self.resources(reader)
                # With subdivision enabled the open boundary rounds inward,
                # leaving the footprint's corners uncovered. A viewport cap
                # of zero restores the original, fully covering quad.
                for enabled, cap, error in ((False, 0, True), (True, 0, False),
                                             (True, 2, True), (False, 0, True)):
                    with self.subTest(enabled=enabled, cap=cap):
                        source.render.use_simplify = enabled
                        source.render.simplify_subdivision = cap
                        bpy.context.view_layer.update()
                        result = self.generate(terrain=terrain)
                        self.assertEqual(bool(result.statuses[0]['error']), error)
                        self.assertEqual(result.statuses, gn.read_status(
                            result.owner, result.owner.modifiers['PADL Buildings']))
                        self.assertEqual(reader.scene.render.use_simplify, enabled)
                        self.assertEqual(reader.scene.render.simplify_subdivision, cap)
                        self.assertEqual(bpy.context.scene, source)
            self.assert_released(resources)
        finally:
            source.render.use_simplify, source.render.simplify_subdivision = original
            bpy.context.view_layer.update()


run(args, StatusSceneChecks)
