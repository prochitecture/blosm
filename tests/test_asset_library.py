"""Run in Blender background with factory startup and --python-exit-code 1."""
from pathlib import Path
from dataclasses import replace
import os
import shutil
import sys
import tempfile
import unittest
from types import SimpleNamespace
import bpy

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from blosm.app.asset_library import (AssetCatalog, AssetError, AssetRef,
                                    AssetUnavailable, BlenderAssetLoader, LibrarySpec)
from blosm.app.library_store import file_digest
from blosm.way.dimension_policy import StreetDimensionPolicy


class AssetLibraryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix='blosm_asset_tests_',
                                             dir=os.environ.get('BLOSM_TEST_TMP'))
        cls.root = Path(cls.tmp.name)
        for filename, value in [('base.blend', .2), ('external.blend', .7)]:
            mat = bpy.data.materials.new('Shared name')
            mat.diffuse_color = (value, .1, .1, 1.)
            tree = bpy.data.node_groups.new('Example GN', 'GeometryNodeTree')
            shared = bpy.data.node_groups.new('Shared dependency GN', 'GeometryNodeTree')
            optional = bpy.data.node_groups.new('Optional GN', 'GeometryNodeTree')
            tree.nodes.new('GeometryNodeGroup').node_tree = shared
            optional.nodes.new('GeometryNodeGroup').node_tree = shared
            mesh = bpy.data.meshes.new('Asset mesh')
            obj = bpy.data.objects.new('Asset object', mesh)
            collection = bpy.data.collections.new('Asset collection')
            collection.objects.link(obj)
            bpy.data.libraries.write(str(cls.root/filename), {mat, tree, optional, collection})
            for data, block in [(bpy.data.collections, collection), (bpy.data.objects, obj),
                                (bpy.data.meshes, mesh), (bpy.data.node_groups, tree),
                                (bpy.data.node_groups, optional), (bpy.data.node_groups, shared),
                                (bpy.data.materials, mat)]:
                data.remove(block)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def setUp(self):
        self.catalog = AssetCatalog([
            LibrarySpec('base', 'base.blend'),
            LibrarySpec('external', str(self.root/'external.blend'), source='external'),
            LibrarySpec('pro', 'base.blend', edition='PRO')])
        self.loader = BlenderAssetLoader(self.root, self.catalog)

    def test_source_identity_and_persistent_reuse(self):
        unrelated = bpy.data.materials.new('Shared name')
        a = self.loader.load(AssetRef('base', 'materials', 'Shared name'))
        b = self.loader.load(AssetRef('external', 'materials', 'Shared name'))
        self.assertNotEqual(a, unrelated)
        self.assertNotEqual(a, b)
        self.assertAlmostEqual(a.diffuse_color[0], .2)
        self.assertAlmostEqual(b.diffuse_color[0], .7)
        fresh = BlenderAssetLoader(self.root, self.catalog)
        self.assertEqual(fresh.load(AssetRef('base','materials','Shared name')), a)

    def test_multiple_datablock_kinds_without_scene_instancing(self):
        refs = [AssetRef('base', k, n) for k, n in [('node_groups', 'Example GN'),
                ('objects', 'Asset object'), ('collections', 'Asset collection')]]
        loaded = self.loader.load_many(refs)
        self.assertEqual(len(loaded), 3)
        self.assertFalse(loaded[refs[1]].users_scene)

    def test_pro_availability_is_explicit(self):
        ref = AssetRef('pro', 'materials', 'Shared name')
        with self.assertRaisesRegex(AssetUnavailable, 'Pro'):
            self.loader.load(ref)
        pro = BlenderAssetLoader(self.root, self.catalog, is_pro=True)
        self.assertIsNotNone(pro.load(ref))

    def test_missing_batch_preflight(self):
        count = len(bpy.data.materials)
        with self.assertRaises(AssetUnavailable):
            self.loader.load_many((AssetRef('base','materials','Shared name'),
                                   AssetRef('base','materials','Missing')))
        self.assertEqual(count, len(bpy.data.materials))

    def test_optional_catalog_query_does_not_load_datablocks(self):
        before = {kind:len(getattr(bpy.data,kind)) for kind in
                  ('objects','meshes','materials','collections','node_groups')}
        self.assertTrue(self.loader.contains(AssetRef('base','node_groups','Example GN')))
        self.assertFalse(self.loader.contains(AssetRef('base','node_groups','Missing optional feature')))
        self.assertEqual(before,{kind:len(getattr(bpy.data,kind)) for kind in before})
        with self.assertRaisesRegex(AssetUnavailable,'Pro'):
            self.loader.contains(AssetRef('pro','node_groups','Example GN'))
        with self.assertRaises(AssetError):
            self.loader.contains(AssetRef('base','not_datablocks','Example GN'))
        self.catalog.register(replace(self.catalog.get('base'),sha256='0'*64),replace=True)
        with self.assertRaisesRegex(AssetError,'fingerprint mismatch'):
            self.loader.contains(AssetRef('base','node_groups','Example GN'))

    def test_path_and_missing_library(self):
        self.catalog.register(LibrarySpec('escape', '../outside.blend'))
        self.catalog.register(LibrarySpec('missing', 'missing.blend'))
        with self.assertRaises(AssetError):
            self.loader.resolve('escape')
        with self.assertRaises(AssetUnavailable):
            self.loader.resolve('missing')

    def test_link_and_append_are_distinct(self):
        ref = AssetRef('external', 'materials', 'Shared name')
        linked = self.loader.load(ref, link=True)
        appended = self.loader.load(ref)
        self.assertIsNotNone(linked.library)
        self.assertIsNone(appended.library)
        self.assertEqual(self.loader.load(ref, link=True), linked)

    def test_managed_linking_requires_certified_pinned_library(self):
        loader = BlenderAssetLoader(self.root, self.catalog, linked_root=self.root/'retained')
        ref = AssetRef('base', 'materials', 'Shared name')
        with self.assertRaisesRegex(AssetError, 'certified self-contained'):
            loader.load(ref, link=True)
        pinned = replace(self.catalog.get('base'), sha256=file_digest(self.root/'base.blend'))
        self.catalog.register(pinned, replace=True)
        with self.assertRaisesRegex(AssetError, 'certified self-contained'):
            loader.load(ref, link=True)
        self.catalog.register(replace(pinned, self_contained=True), replace=True)
        linked = loader.load(ref, link=True)
        expected = loader.linked_store.path(pinned.sha256)
        self.assertEqual(Path(bpy.path.abspath(linked.library.filepath)).resolve(), expected)
        self.assertEqual(loader.load(ref, link=True), linked)
        self.assertIsNone(loader.load(ref).library)

    def test_retained_pro_library_does_not_bypass_edition_check(self):
        pinned = replace(self.catalog.get('pro'), sha256=file_digest(self.root/'base.blend'),
                         self_contained=True)
        self.catalog.register(pinned, replace=True)
        root = self.root/'retained_pro'
        pro = BlenderAssetLoader(self.root, self.catalog, is_pro=True, linked_root=root)
        ref = AssetRef('pro', 'materials', 'Shared name')
        pro.load(ref, link=True)
        base = BlenderAssetLoader(self.root, self.catalog, linked_root=root)
        with self.assertRaisesRegex(AssetUnavailable, 'requires Blosm Pro'):
            base.load(ref, link=True)

    def test_pinned_mismatch_prevents_loading_or_reusing_ids(self):
        # Preflight must run even if a previous call cached the same ID.
        ref = AssetRef('external', 'materials', 'Shared name')
        self.loader.load(ref)
        self.catalog.register(replace(self.catalog.get('external'), sha256='0'*64), replace=True)
        count = len(bpy.data.materials)
        with self.assertRaisesRegex(AssetError, 'fingerprint mismatch'):
            self.loader.load(ref)
        self.assertEqual(len(bpy.data.materials), count)

    def test_pinned_append_identity_includes_bytes(self):
        source = self.root/'updated.blend'
        shutil.copyfile(self.root/'base.blend', source)
        entry = LibrarySpec('updated', str(source), source='external', sha256=file_digest(source))
        self.catalog.register(entry)
        ref = AssetRef('updated', 'materials', 'Shared name')
        first = self.loader.load(ref)
        shutil.copyfile(self.root/'external.blend', source)
        self.catalog.register(replace(entry, sha256=file_digest(source)), replace=True)
        second = self.loader.load(ref)
        self.assertNotEqual(first, second)
        self.assertAlmostEqual(first.diffuse_color[0], .2)
        self.assertAlmostEqual(second.diffuse_color[0], .7)

    def test_pinned_link_cannot_reuse_a_mutable_source_path(self):
        self.catalog.register(replace(self.catalog.get('base'), sha256=file_digest(self.root/'base.blend')),
                              replace=True)
        with self.assertRaisesRegex(AssetError, 'requires a retained library root'):
            self.loader.load(AssetRef('base', 'materials', 'Shared name'), link=True)

    def test_dimension_style_aliases_and_validation(self):
        class Street:
            def __init__(self, values): self.values = values
            def getStyleBlockAttr(self, name): return self.values.get(name)
        policy = StreetDimensionPolicy()
        for key in ('totalLanes', 'totalNumLanes'):
            values = policy.resolve(Street({key: 3, 'laneWidth': 2.8}), 'primary', False)
            self.assertEqual(values, {'nrLanes': 3, 'laneWidth': 2.8})
        self.assertEqual(policy.resolve(Street({}), 'service', False)['nrLanes'], 1)
        with self.assertRaises(ValueError):
            policy.resolve(Street({'totalLanes':0}), 'primary', False)

    def test_app_policy_mixed_batch_and_late_optional_dependency_reuse(self):
        self.catalog.register(replace(self.catalog.get('base'), sha256=file_digest(self.root/'base.blend'),
                                      self_contained=True), replace=True)
        app = SimpleNamespace(basePath=str(self.root), isPro=False, assetLibraryMode='AUTO',
                              assetLibraryDirectory=str(self.root/'app-retained'))
        loader = BlenderAssetLoader.for_app(app, self.catalog)
        # The real app uses basePath/assets. Fixture sources are external to it.
        self.catalog.register(replace(self.catalog.get('base'), path=str(self.root/'base.blend'),
                                      source='external'), replace=True)
        main = AssetRef('base', 'node_groups', 'Example GN')
        other = AssetRef('external', 'materials', 'Shared name')
        values = loader.load_many((main, other))
        self.assertIsNotNone(values[main].library)
        self.assertIsNone(values[other].library)
        optional = loader.load(AssetRef('base', 'node_groups', 'Optional GN'))
        first_child = next(n.node_tree for n in values[main].nodes if n.type == 'GROUP')
        second_child = next(n.node_tree for n in optional.nodes if n.type == 'GROUP')
        self.assertEqual(first_child, second_child)
        self.assertEqual(optional.library, values[main].library)
        self.assertEqual(sum(t.library == optional.library and t.name == 'Shared dependency GN'
                             for t in bpy.data.node_groups), 1)
        embedded = loader.load(main, link=False)
        self.assertIsNone(embedded.library)
        self.assertNotEqual(embedded, values[main])

    def test_automatic_mixed_batch_preflight_prevents_partial_load(self):
        self.catalog.register(replace(self.catalog.get('base'), sha256=file_digest(self.root/'base.blend'),
                                      self_contained=True), replace=True)
        loader = BlenderAssetLoader(self.root, self.catalog, policy='AUTO', linked_root=self.root/'mixed-retained')
        before = (len(bpy.data.node_groups), len(bpy.data.materials))
        with self.assertRaises(AssetUnavailable):
            loader.load_many((AssetRef('base','node_groups','Example GN'),
                              AssetRef('external','materials','Missing')))
        self.assertEqual(before, (len(bpy.data.node_groups), len(bpy.data.materials)))

    def test_registered_import_settings_reach_the_app_loader(self):
        from blosm import gui
        from blosm.app.blender import BlenderApp
        registered = []
        added_pointer = not hasattr(bpy.types.Scene, 'blosm')
        saved_settings = None
        try:
            for cls in (gui.BlosmDefaultLevelsEntry, gui.BlosmProperties):
                if not cls.is_registered:
                    bpy.utils.register_class(cls); registered.append(cls)
            if added_pointer:
                bpy.types.Scene.blosm = bpy.props.PointerProperty(type=gui.BlosmProperties)
            settings = bpy.context.scene.blosm
            saved_settings = settings.assetLibraryMode, settings.assetLibraryDirectory
            self.assertEqual(settings.assetLibraryMode, 'AUTO')
            settings.assetLibraryDirectory = str(self.root/'import-retained')
            app = SimpleNamespace(basePath=str(self.root))
            BlenderApp.setAttributes(app, bpy.context)
            loader = BlenderAssetLoader.for_app(app, self.catalog)
            self.assertEqual(loader.policy.mode, 'AUTO')
            self.assertEqual(loader.linked_store.root, self.root/'import-retained')
            settings.assetLibraryMode = 'APPEND'
            BlenderApp.setAttributes(app, bpy.context)
            self.assertIsNone(BlenderAssetLoader.for_app(app, self.catalog).linked_store)
        finally:
            if saved_settings is not None:
                settings.assetLibraryMode, settings.assetLibraryDirectory = saved_settings
            if added_pointer and hasattr(bpy.types.Scene, 'blosm'):
                del bpy.types.Scene.blosm
            for cls in reversed(registered): bpy.utils.unregister_class(cls)


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(AssetLibraryTests)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from test_asset_loading_policy import AssetLoadingPolicyTests
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(AssetLoadingPolicyTests))
    from test_asset_recovery import AssetRecoveryTests
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(AssetRecoveryTests))
    # Last: this check intentionally reopens a small saved fixture scene.
    from test_asset_library_recovery import RetainedLibraryRecoveryTests
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(RetainedLibraryRecoveryTests))
    if not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful():
        raise RuntimeError('Asset library regression failed')
