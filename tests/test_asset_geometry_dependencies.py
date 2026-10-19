"""Native GN object dependencies loaded through the feature-independent service.

Run in a disposable background Blender. A construction pattern is library data,
not a visible scene instance. Check embedded and retained loading, shared users,
modifier inputs, and persistence without registering a street update handler.
"""
import os
from pathlib import Path
import sys
import tempfile
import unittest

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from blosm.app.asset_library import AssetCatalog, AssetRef, BlenderAssetLoader, LibrarySpec
from blosm.app.library_store import file_digest


class AssetGeometryDependencyTests(unittest.TestCase):
    def test_uninstanced_pattern_append_and_retained_reopen(self):
        with tempfile.TemporaryDirectory(prefix='blosm_pattern_library_',
                                         dir=os.environ.get('BLOSM_TEST_TMP')) as directory:
            root = Path(directory)
            bpy.ops.wm.read_factory_settings(use_empty=True)
            mesh = bpy.data.meshes.new('Shared pattern mesh')
            mesh.from_pydata([(0, 0, 0), (2, 0, 0), (2, 1, 0), (0, 1, 0)], [], [(0, 1, 2, 3)])
            marker = mesh.attributes.new('pattern_marker', 'INT', 'POINT')
            marker.data.foreach_set('value', [17, 18, 19, 20])
            pattern = bpy.data.objects.new('Shared pattern object', mesh)
            tree = bpy.data.node_groups.new('Pattern consumer', 'GeometryNodeTree')
            tree.interface.new_socket(name='Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
            height = tree.interface.new_socket(name='Height', in_out='INPUT', socket_type='NodeSocketFloat')
            inp = tree.nodes.new('NodeGroupInput')
            out = tree.nodes.new('NodeGroupOutput')
            info = tree.nodes.new('GeometryNodeObjectInfo')
            info.inputs['Object'].default_value = pattern
            info.transform_space = 'ORIGINAL'
            xyz = tree.nodes.new('ShaderNodeCombineXYZ')
            move = tree.nodes.new('GeometryNodeSetPosition')
            tree.links.new(inp.outputs['Height'], xyz.inputs['Z'])
            tree.links.new(info.outputs['Geometry'], move.inputs['Geometry'])
            tree.links.new(xyz.outputs[0], move.inputs['Offset'])
            tree.links.new(move.outputs[0], out.inputs['Geometry'])
            height_id = height.identifier
            library = root/'patterns.blend'
            bpy.data.libraries.write(str(library), {tree}, fake_user=True)
            bpy.ops.wm.read_factory_settings(use_empty=True)
            catalog = AssetCatalog([LibrarySpec('patterns', library.name,
                sha256=file_digest(library), self_contained=True)])
            ref = AssetRef('patterns', 'node_groups', 'Pattern consumer')
            for linked in (False, True):
                with self.subTest(linked=linked):
                    loader = BlenderAssetLoader(root, catalog, linked_root=root/'retained')
                    loaded = loader.load(ref, link=linked)
                    dependency = next(n for n in loaded.nodes if n.type == 'OBJECT_INFO').inputs['Object'].default_value
                    self.assertFalse(dependency.users_scene)
                    self.assertEqual(bool(loaded.library), linked)
                    self.assertEqual(bool(dependency.library), linked)
                    receivers = []
                    for index in range(2):
                        receiver = bpy.data.objects.new('Receiver '+str(index), bpy.data.meshes.new('Empty source'))
                        bpy.context.scene.collection.objects.link(receiver)
                        modifier = receiver.modifiers.new('Pattern', 'NODES')
                        modifier.node_group = loader.load(ref, link=linked)
                        self.assertEqual(modifier.node_group, loaded)
                        getattr(modifier.properties.inputs, height_id).value = .25*(index+1)
                        receivers.append(receiver)
                    def check():
                        bpy.context.view_layer.update()
                        for index, receiver in enumerate(receivers):
                            result = receiver.evaluated_get(bpy.context.evaluated_depsgraph_get()).data
                            self.assertEqual(len(result.polygons), 1)
                            self.assertEqual(len(result.vertices), 4)
                            self.assertEqual([round(v.co.z, 6) for v in result.vertices], [.25*(index+1)]*4)
                            self.assertEqual([p.value for p in result.attributes['pattern_marker'].data], [17, 18, 19, 20])
                    check()
                    self.assertFalse(dependency.users_scene)
                    scene = root/('linked.blend' if linked else 'embedded.blend')
                    bpy.ops.wm.save_as_mainfile(filepath=str(scene))
                    bpy.ops.wm.open_mainfile(filepath=str(scene), use_scripts=False, load_ui=False)
                    receivers = [bpy.data.objects['Receiver '+str(i)] for i in range(2)]
                    check()
                    loaded = receivers[0].modifiers[0].node_group
                    self.assertEqual(receivers[1].modifiers[0].node_group, loaded)
                    dependency = next(n for n in loaded.nodes if n.type == 'OBJECT_INFO').inputs['Object'].default_value
                    self.assertFalse(dependency.users_scene)
                    fresh = BlenderAssetLoader(root, catalog, linked_root=root/'retained')
                    self.assertEqual(fresh.load(ref, link=linked), loaded)
                    bpy.ops.wm.read_factory_settings(use_empty=True)


if __name__ == '__main__':
    unittest.main(argv=[sys.argv[0]])
