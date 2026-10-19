"""Native saved-file recovery through the actual Blosm operator.

Run in a disposable background Blender; reopening the test scene is intentional.
"""
from pathlib import Path
import os
import sys
import tempfile
import unittest

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from blosm.app.asset_library import AssetCatalog, AssetRef, BlenderAssetLoader, LibrarySpec
from blosm.app.library_store import file_digest
from blosm.gui.asset_libraries import classes, BLOSM_PT_AssetLibraries


class RetainedLibraryRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='blosm_native_recovery_',
                                                     dir=os.environ.get('BLOSM_TEST_TMP'))
        self.root = Path(self.temporary.name).resolve()
        if os.environ.get('BLOSM_TEST_TMP'):
            assert self.root.is_relative_to(Path(os.environ['BLOSM_TEST_TMP']).resolve())
        self.addCleanup(self.temporary.cleanup)
        from blosm import gui
        self.registered = []
        for cls in (gui.BlosmDefaultLevelsEntry, gui.BlosmProperties, *classes):
            if not cls.is_registered:
                bpy.utils.register_class(cls); self.registered.append(cls)
        self.added_pointer = not hasattr(bpy.types.Scene, 'blosm')
        if self.added_pointer:
            bpy.types.Scene.blosm = bpy.props.PointerProperty(type=gui.BlosmProperties)
        self.addCleanup(self.unregister)

    def unregister(self):
        if self.added_pointer:
            del bpy.types.Scene.blosm
        for cls in reversed(self.registered): bpy.utils.unregister_class(cls)

    def fixture(self):
        tree = bpy.data.node_groups.new('Recovery GN', 'GeometryNodeTree')
        tree.interface.new_socket(name='Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
        cube = tree.nodes.new('GeometryNodeMeshCube'); cube.inputs['Size'].default_value = (2., 3., 4.)
        out = tree.nodes.new('NodeGroupOutput'); tree.links.new(cube.outputs['Mesh'], out.inputs['Geometry'])
        tree['fixture_revision'] = 'exact saved revision'
        source = self.root/'backup.blend'
        bpy.data.libraries.write(str(source), {tree})
        bpy.data.node_groups.remove(tree)
        self.digest = file_digest(source)
        spec = LibrarySpec('test', 'backup.blend', sha256=self.digest, self_contained=True)
        loader = BlenderAssetLoader(self.root, AssetCatalog((spec,)),
                                    linked_root=self.root/'retained', policy='RETAINED')
        tree = loader.load(AssetRef('test', 'node_groups', 'Recovery GN'))
        obj = bpy.data.objects.new('Recovery output', bpy.data.meshes.new('Recovery input'))
        bpy.context.scene.collection.objects.link(obj)
        obj['editable_scene_marker'] = 17
        obj.modifiers.new('Live GN', 'NODES').node_group = tree
        self.assertEqual(self.faces(), 6)
        self.project = self.root/'project.blend'
        bpy.ops.wm.save_as_mainfile(filepath=str(self.project), compress=True)
        target = loader.linked_store.path(self.digest)
        self.assertTrue(target.resolve().is_relative_to(self.root))
        target.rename(self.root/'removed-library.blend')
        bpy.ops.wm.open_mainfile(filepath=str(self.project), load_ui=False, use_scripts=False)
        return source, target

    @staticmethod
    def library():
        return bpy.data.objects['Recovery output'].modifiers[0].node_group.library

    @staticmethod
    def faces():
        bpy.context.view_layer.update()
        obj = bpy.data.objects['Recovery output'].evaluated_get(bpy.context.evaluated_depsgraph_get())
        mesh = obj.to_mesh()
        try: return len(mesh.polygons)
        finally: obj.to_mesh_clear()

    def test_missing_saved_library_recovers_then_relocates_without_revision_change(self):
        source, target = self.fixture()
        self.assertTrue(BLOSM_PT_AssetLibraries.poll(bpy.context))
        self.assertTrue(bpy.data.objects['Recovery output'].modifiers[0].node_group.is_missing)
        library_name = self.library().name
        wrong = self.root/'wrong.blend'; wrong.write_bytes(b'Different revision')
        with self.assertRaisesRegex(RuntimeError, 'fingerprint mismatch'):
            bpy.ops.blosm.restore_retained_library(library_name=library_name,
                expected_digest=self.digest, filepath=str(wrong))
        self.assertFalse(target.exists())
        result = bpy.ops.blosm.restore_retained_library(library_name=library_name,
            expected_digest=self.digest, filepath=str(source))
        self.assertEqual(result, {'FINISHED'})
        self.assertEqual(self.faces(), 6)
        self.assertEqual(bpy.data.objects['Recovery output']['editable_scene_marker'], 17)
        self.assertEqual(self.library().filepath, str(target))
        self.assertEqual(file_digest(target), self.digest)

        # Reopen the original saved project with the file absent again, then
        # exercise the explicit destination-directory option through real RNA.
        target.unlink()
        bpy.ops.wm.open_mainfile(filepath=str(self.project), load_ui=False, use_scripts=False)
        bpy.context.scene.blosm.assetLibraryDirectory = str(self.root/'relocated')
        result = bpy.ops.blosm.restore_retained_library(library_name=self.library().name,
            expected_digest=self.digest, filepath=str(source), use_current_directory=True)
        self.assertEqual(result, {'FINISHED'})
        moved = self.root/'relocated'/self.digest/'library.blend'
        self.assertEqual(Path(self.library().filepath), moved)
        self.assertFalse(target.exists())
        self.assertEqual(self.faces(), 6)
        recovered = self.root/'recovered.blend'
        bpy.ops.wm.save_as_mainfile(filepath=str(recovered), compress=True)
        source.unlink()
        bpy.ops.wm.open_mainfile(filepath=str(recovered), load_ui=False, use_scripts=False)
        self.assertEqual(self.faces(), 6)
        self.assertFalse(bpy.data.objects['Recovery output'].modifiers[0].node_group.is_missing)
        self.assertEqual(bpy.data.objects['Recovery output'].modifiers[0].node_group['fixture_revision'],
                         'exact saved revision')
        self.assertEqual(Path(self.library().filepath), moved)


if __name__ == '__main__':
    if not bpy.app.background: raise RuntimeError('Use a disposable background Blender')
    if not unittest.TextTestRunner(verbosity=2).run(
            unittest.defaultTestLoader.loadTestsFromTestCase(RetainedLibraryRecoveryTests)).wasSuccessful():
        raise RuntimeError('Native library recovery failed')
