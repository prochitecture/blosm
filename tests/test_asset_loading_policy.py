"""Storage selection and retained-path preflight, without a Blender process."""
import importlib
import os
from pathlib import Path
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch

# Import only the feature-independent app modules, not the Blender addon entry.
package = ModuleType('blosm_asset_policy_tests')
package.__path__ = [str(Path(__file__).resolve().parents[1]/'app')]
sys.modules[package.__name__] = package
assets = importlib.import_module(package.__name__+'.asset_library')
policy = importlib.import_module(package.__name__+'.asset_policy')


class AssetLoadingPolicyTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='blosm_policy_', dir=os.environ.get('BLOSM_TEST_TMP'))
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root/'base.blend'
        self.source.write_bytes(b'Fixture bytes for storage preflight; no Blender data parsed')
        self.spec = assets.LibrarySpec('base', 'base.blend')
        self.certified = assets.LibrarySpec('certified', 'base.blend',
            sha256=assets.file_digest(self.source), self_contained=True)
        self.catalog = assets.AssetCatalog((self.spec, self.certified))
        self.store = self.root/'retained'

    def loader(self, mode='AUTO'):
        return assets.BlenderAssetLoader(self.root, self.catalog, linked_root=self.store, policy=mode)

    def test_automatic_retains_only_certified_revisions(self):
        loader = self.loader()
        self.assertEqual(loader.resolve('base')[1], self.source)
        self.assertFalse(self.store.exists())
        retained = loader.resolve('certified')[1]
        self.assertEqual(retained.read_bytes(), self.source.read_bytes())
        self.assertNotEqual(retained, self.source)
        self.source.unlink()
        self.assertEqual(loader.resolve('certified')[1], retained)

    def test_required_retention_does_not_silently_embed(self):
        with self.assertRaisesRegex(assets.AssetError, 'certified self-contained'):
            self.loader('RETAINED').resolve('base')
        self.assertFalse(self.store.exists())

    def test_automatic_never_falls_back_after_store_corruption(self):
        loader = self.loader()
        retained = loader.resolve('certified')[1]
        retained.write_bytes(b'damaged')
        with self.assertRaisesRegex(assets.AssetError, 'damaged or was replaced'):
            loader.resolve('certified')

    def test_explicit_override_and_legacy_constructor(self):
        loader = self.loader()
        self.assertTrue(loader.uses_link('certified'))
        self.assertFalse(loader.uses_link('certified', False))
        self.assertEqual(loader.resolve('certified', link=False)[1], self.source)
        legacy = assets.BlenderAssetLoader(self.root, self.catalog)
        self.assertFalse(legacy.uses_link('certified'))
        with self.assertRaisesRegex(assets.AssetError, 'durable library directory'):
            assets.BlenderAssetLoader(self.root, self.catalog, policy='AUTO')
        with self.assertRaises(TypeError):
            loader.uses_link('base', 1)
        with self.assertRaises(ValueError):
            policy.AssetLoadPolicy('unknown')

    def fake_blender(self, saved=False):
        return SimpleNamespace(
            data=SimpleNamespace(filepath=str(self.root/'project.blend') if saved else ''),
            path=SimpleNamespace(abspath=lambda value: str(self.root/value[2:]) if value.startswith('//') else value),
            utils=SimpleNamespace(resource_path=lambda kind: str(self.root/'blender_user')))

    def test_app_default_resolves_user_storage_without_creating_it(self):
        with patch.dict(sys.modules, bpy=self.fake_blender()):
            loader = assets.BlenderAssetLoader.for_app(SimpleNamespace(basePath=str(self.root)), self.catalog)
        self.assertEqual(loader.policy.mode, 'AUTO')
        self.assertEqual(loader.linked_store.root, self.root/'blender_user'/'datafiles'/'blosm'/'libraries')
        self.assertFalse(loader.linked_store.root.exists())

    def test_app_scene_relative_store_requires_saved_project(self):
        app = SimpleNamespace(basePath=str(self.root), assetLibraryMode='AUTO', assetLibraryDirectory='//libraries')
        with patch.dict(sys.modules, bpy=self.fake_blender()):
            with self.assertRaisesRegex(assets.AssetError, 'Save the scene'):
                assets.BlenderAssetLoader.for_app(app, self.catalog)
        with patch.dict(sys.modules, bpy=self.fake_blender(saved=True)):
            loader = assets.BlenderAssetLoader.for_app(app, self.catalog)
        self.assertEqual(loader.linked_store.root, self.root/'libraries')
        app.assetLibraryDirectory = 'ambiguous-relative'
        with patch.dict(sys.modules, bpy=self.fake_blender(saved=True)):
            with self.assertRaisesRegex(assets.AssetError, 'absolute or saved-project-relative'):
                assets.BlenderAssetLoader.for_app(app, self.catalog)

    def test_embedded_app_does_not_require_a_library_directory(self):
        app = SimpleNamespace(basePath=str(self.root), assetLibraryMode='APPEND', assetLibraryDirectory='//unused')
        with patch.dict(sys.modules, bpy=self.fake_blender()):
            loader = assets.BlenderAssetLoader.for_app(app, self.catalog)
        self.assertIsNone(loader.linked_store)
        self.assertFalse(loader.uses_link('certified'))


if __name__ == '__main__':
    unittest.main()
