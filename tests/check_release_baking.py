"""Actual CPU diffuse/emission bakes, grouping/replacement and failure recovery."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
import json
from math import radians
import unittest
from unittest.mock import patch
import bpy
import numpy as np
from release_support import arguments, enable, reset, geometry, assert_geometry, meshes, REPO, PACKAGE, run

args = arguments()
enable(args)
baking = importlib.import_module(PACKAGE + '.util.blender_extra.baking')
baking.ensure_registered()


def material(channel, name='Known color'):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    nodes.clear()
    shader = nodes.new('ShaderNodeEmission' if channel == 'EMIT' else 'ShaderNodeBsdfDiffuse')
    shader.inputs['Color'].default_value = (0.8, 0.1, 0.2, 1)
    output = nodes.new('ShaderNodeOutputMaterial')
    mat.node_tree.links.new(shader.outputs[0], output.inputs['Surface'])
    return mat


def state():
    scene = bpy.context.scene
    bake = scene.render.bake
    return (scene.render.engine, scene.cycles.bake_type, bake.use_selected_to_active,
            bake.use_pass_direct, bake.use_pass_indirect, bake.use_pass_color,
            bake.margin, round(bake.cage_extrusion, 5))


class BakingChecks(unittest.TestCase):
    def check_tile_pipeline(self, joined):
        from gltf_test_support import fixture, make_glb, make_b3dm, no_tile_writes_or_operator
        manager = importlib.import_module(PACKAGE + '.threed_tiles.manager')
        settings = reset(args)
        settings.dataType = '3d-tiles'
        settings.threedTilesSource = 'custom'
        settings.threedTilesUrl = 'https://fixture.invalid/root.json'
        settings.progressiveImport = False
        settings.concurrentDownloads = False
        settings.join3dTilesObjects = joined
        settings.minLon, settings.maxLon = 24.93, 24.95
        settings.minLat, settings.maxLat = 60.16, 60.18
        doc, binary, _ = fixture()
        center = manager.BaseManager.fromGeographic(60.17, 24.94, 0)
        volume = {'region': [radians(24.93), radians(60.16), radians(24.95), radians(60.18), -100, 100]}
        root = {'asset': {'version': '1.0'}, 'geometricError': 0, 'root': {
            'geometricError': 0, 'boundingVolume': volume, 'children': [
                {'geometricError': 0, 'boundingVolume': volume, 'content': {'uri': name}}
                for name in ('one.b3dm', 'two.b3dm')]}}
        payloads = {'root.json': json.dumps(root).encode()}
        for i, name in enumerate(('one.b3dm', 'two.b3dm')):
            doc['nodes'][0]['translation'] = [i*5, 0, 0]
            payloads[name] = make_b3dm(make_glb(doc, binary), list(center))
        with patch.object(manager.BaseManager, 'download',
                          lambda self, uri: payloads[uri.rsplit('/', 1)[-1]]), no_tile_writes_or_operator():
            self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        objects = meshes()
        self.assertEqual(len(objects), 1 if joined else 2)
        expected = geometry(objects)
        bpy.ops.object.select_all(action='DESELECT')
        for obj in objects:
            obj.select_set(True)
        bpy.context.view_layer.objects.active = objects[0]
        settings.bake_image_size = 128
        settings.bake_replace_active = True
        bpy.context.scene.cycles.device = 'CPU'
        bpy.context.scene.cycles.samples = 1
        bpy.context.scene.render.threads_mode = 'FIXED'
        bpy.context.scene.render.threads = 2
        self.assertEqual(bpy.ops.blosm.bake_to_single_material(), {'FINISHED'})
        target = bpy.context.object
        assert_geometry(self, expected, geometry([target]))
        image = next(n.image for n in target.data.materials[0].node_tree.nodes if n.type == 'TEX_IMAGE')
        self.assertIsNotNone(image.packed_file)
        pixels = np.array(image.pixels[:]).reshape(-1, 4)
        self.assertGreater(np.mean(pixels[:, 0] > 0.5), 0.05, 'Red texture region lost')
        self.assertGreater(np.mean(pixels[:, 1] > 0.5), 0.05, 'Green texture region lost')
        output = args.artifacts / f'import_bake_export_joined_{joined}.glb'
        self.assertEqual(bpy.ops.export_scene.gltf(filepath=str(output), use_selection=True, export_format='GLB'),
                         {'FINISHED'})
        before = set(bpy.data.objects)
        self.assertEqual(bpy.ops.import_scene.gltf(filepath=str(output)), {'FINISHED'})
        imported = [o for o in meshes() if o not in before]
        assert_geometry(self, expected, geometry(imported))
        self.assertTrue(all(o.data.materials and o.data.uv_layers for o in imported))

    def test_tiles_import_bake_export_joined(self):
        self.check_tile_pipeline(True)

    def test_tiles_import_bake_export_separate(self):
        self.check_tile_pipeline(False)

    def sources(self, channel, joined):
        settings = reset(args)
        settings.dataType = '3d-tiles'
        settings.bake_image_size = 128
        settings.bake_margin_px = 2
        settings.bake_cage_extrusion = 0.02
        bpy.context.scene.cycles.device = 'CPU'
        bpy.context.scene.cycles.samples = 1
        bpy.context.scene.render.threads_mode = 'FIXED'
        bpy.context.scene.render.threads = 2
        objects = []
        for x in (0, 3):
            bpy.ops.mesh.primitive_plane_add(size=2, location=(x, 0, 1))
            obj = bpy.context.object
            obj.data.materials.append(material(channel))
            objects.append(obj)
        bpy.ops.object.select_all(action='DESELECT')
        for obj in objects:
            obj.select_set(True)
        bpy.context.view_layer.objects.active = objects[0]
        if joined:
            bpy.ops.object.join()
            objects = [bpy.context.object]
        return settings, objects

    def check_bake(self, channel, joined, replace):
        settings, sources = self.sources(channel, joined)
        settings.bake_replace_active = replace
        expected = geometry(sources)
        source_meshes = [o.data for o in sources]
        before = state()
        self.assertEqual(bpy.ops.blosm.bake_to_single_material(), {'FINISHED'})
        self.assertEqual(state(), before, 'Baking changed scene render settings')
        target = bpy.context.object
        assert_geometry(self, expected, geometry([target]))
        self.assertEqual(len(target.data.materials), 1)
        self.assertTrue(target.data.uv_layers)
        for loop in target.data.uv_layers.active.data:
            self.assertTrue(all(-0.001 <= v <= 1.001 for v in loop.uv))
        images = [n.image for n in target.data.materials[0].node_tree.nodes if n.type == 'TEX_IMAGE']
        self.assertEqual(len(images), 1)
        image = images[0]
        self.assertEqual(tuple(image.size), (128, 128))
        self.assertIsNotNone(image.packed_file)
        pixels = np.array(image.pixels[:]).reshape(-1, 4)
        rgb = pixels[:, :3]
        if image.colorspace_settings.name == 'sRGB':
            rgb = np.where(rgb <= 0.04045, rgb/12.92, ((rgb + 0.055)/1.055)**2.4)
        error = np.max(np.abs(rgb - (0.8, 0.1, 0.2)), axis=1)
        self.assertGreater(np.mean(error < 0.06), 0.5,
                           f'Baked pixels lost the known surface color: min={pixels.min(axis=0)}, '
                           f'max={pixels.max(axis=0)}, median={np.median(pixels, axis=0)}, '
                           f'colorspace={image.colorspace_settings.name}')
        if replace and joined:
            self.assertEqual(target, sources[0])
            self.assertNotEqual(target.data, source_meshes[0])
        else:
            self.assertNotIn(target, sources)
            self.assertEqual([o.data for o in sources], source_meshes)
            assert_geometry(self, expected, geometry(sources))
        self.assertEqual(len([o for o in bpy.context.scene.objects if o.type == 'MESH']),
                         len(sources) + (0 if replace and joined else 1))
        if channel == 'EMIT' and joined and replace:
            path = args.artifacts / 'packed_bake.blend'
            name = image.name
            bpy.ops.wm.save_as_mainfile(filepath=str(path))
            bpy.ops.wm.open_mainfile(filepath=str(path), load_ui=False)
            self.assertIsNotNone(bpy.data.images[name].packed_file)

    def test_channel_selection(self):
        selector = baking.TilesAutoSelector()
        diffuse, emit = material('DIFFUSE'), material('EMIT')
        self.assertEqual(selector.choose_channel([diffuse]), baking.BakeChannel.DIFFUSE)
        self.assertEqual(selector.choose_channel([emit]), baking.BakeChannel.EMIT)
        self.assertIsNone(selector.choose_channel([diffuse, emit]))
        self.assertIsNone(selector.choose_channel([]))

    def test_bake_failure_restores_settings_and_objects(self):
        _, sources = self.sources('DIFFUSE', False)
        before = state()
        objects = set(bpy.data.objects.keys())
        expected = geometry(sources)
        real_object_ops = bpy.ops.object
        class ObjectOps:
            def __getattr__(self, name):
                return getattr(real_object_ops, name)
            def bake(self, **kwargs):
                raise RuntimeError('Injected bake failure')
        with patch.object(bpy.ops, 'object', ObjectOps()):
            with self.assertRaisesRegex(RuntimeError, 'Injected bake failure'):
                bpy.ops.blosm.bake_to_single_material()
        self.assertEqual(before, state())
        self.assertEqual(set(bpy.data.objects.keys()), objects)
        assert_geometry(self, expected, geometry(sources))


for channel in ('DIFFUSE', 'EMIT'):
    for joined in (False, True):
        for replace in (False, True):
            def check(self, channel=channel, joined=joined, replace=replace):
                self.check_bake(channel, joined, replace)
            setattr(BakingChecks, f'test_{channel}_joined_{joined}_replace_{replace}', check)

run(args, BakingChecks)
