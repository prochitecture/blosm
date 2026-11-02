"""Native regression for shared street materials and Python-only installation.

Use factory-startup background Blender. --asset-root can select a freshly
material library before installation. --output retains a small Python
scene to reopen independently with --reopen (no editing scripts required).
"""
import argparse
from dataclasses import replace
from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import bpy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
from blosm.app.asset_library import BlenderAssetLoader, AssetUnavailable
from blosm.app.library_store import file_digest
from blosm.way.street_network import assets
from blosm.way.street_network.design import StreetDesign
from blosm.way.street_network.python_backend import model
from blosm.way.street_network.python_backend.renderer import PythonStreetRenderer
from blosm.way.street_network.python_backend.scene import NetworkScene, clear_runtime

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--asset-root', type=Path, default=ROOT/'assets')
parser.add_argument('--output', type=Path)
parser.add_argument('--reopen', type=Path)
parser.add_argument('--temp-root', type=Path, help='Writable parent for temporary test libraries')
OPTIONS = parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])


def temporary_library_directory():
    root = OPTIONS.temp_root or (OPTIONS.output.parent if OPTIONS.output else None)
    if root:
        root.mkdir(parents=True, exist_ok=True)
    return TemporaryDirectory(dir=root)


class StreetMaterialLibraryTests(unittest.TestCase):
    def test_material_library_contains_no_geometry(self):
        with bpy.data.libraries.load(str(OPTIONS.asset_root/assets.MATERIAL_LIBRARY.path)) as (data, _):
            self.assertEqual(set(data.materials), {assets.PAVEMENT.name, assets.GRASS.name,
                                                 assets.DIAGNOSTICS.name, assets.DIAGNOSTIC_RED.name})
            for kind in ('node_groups', 'objects', 'meshes', 'texts', 'images'):
                self.assertEqual(getattr(data, kind), [], kind)

    def test_python_mesh_and_update_need_only_the_material_library(self):
        initial_groups = set(bpy.data.node_groups)
        with temporary_library_directory() as directory:
            root = Path(directory)
            (root/'assets').mkdir()
            shutil.copy2(OPTIONS.asset_root/assets.MATERIAL_LIBRARY.path, root/'assets'/assets.MATERIAL_LIBRARY.path)
            app = SimpleNamespace(basePath=str(root), assetLibraryMode='APPEND')
            renderer = PythonStreetRenderer(app, design=StreetDesign())
            renderer.prepare()
            with self.assertRaises(AssetUnavailable):
                renderer.loader.resolve('streets.native.v8')
            blocks = renderer.loader.load_many((assets.PAVEMENT, assets.GRASS))
            materials = {'pavement': blocks[assets.PAVEMENT], 'grass': blocks[assets.GRASS]}
            sections = (model.Section(1, 1, ((0., 0., 0.), (12., 0., 0.), (24., 0., 0.))),
                        model.Section(2, 2, ((0., 10., 0.), (12., 10., 0.), (24., 10., 0.))))
            bundle = model.Bundle(7, (model.BundleMember((model.Traversal(2),)),
                                      model.BundleMember((model.Traversal(1),))))
            controller = NetworkScene.create(model.Network('shared_material_check', sections, bundles=(bundle,)),
                name='Shared material Python test', materials=materials)
            def verify():
                outputs = [o for o in controller['st_py_outputs'].objects if o.type == 'MESH']
                self.assertTrue(outputs)
                used = {m for obj in outputs for m in obj.data.materials}
                self.assertEqual(used, set(materials.values()))
                self.assertTrue(all(obj.data.polygons and obj.data.materials for obj in outputs))
                self.assertFalse(any(m.type == 'NODES' for obj in outputs for m in obj.modifiers))
            verify()
            source = next(iter(controller['st_py_sources'].objects))
            source.data.vertices[1].co.y += .1
            source.data.update()
            NetworkScene.update(controller)
            verify()
            self.assertEqual(initial_groups, set(bpy.data.node_groups))
            if OPTIONS.output:
                OPTIONS.output.parent.mkdir(parents=True, exist_ok=True)
                bpy.context.preferences.filepaths.file_preview_type = 'NONE'
                bpy.ops.wm.save_as_mainfile(filepath=str(OPTIONS.output.resolve()))
            clear_runtime()

    def test_external_material_selection_is_independent_of_solver(self):
        path = str((OPTIONS.asset_root/assets.MATERIAL_LIBRARY.path).resolve())
        loader = BlenderAssetLoader(ROOT/'assets', assets.material_catalog(path))
        self.assertEqual(loader.resolve(assets.MATERIAL_LIBRARY.key)[1], Path(path))
        material = loader.load(assets.PAVEMENT)
        self.assertEqual(material.name.split('.')[0], assets.PAVEMENT.name)
        self.assertEqual(material['blosm_material_revision'], assets.MATERIAL_LIBRARY.revision)

    def test_retained_materials_share_one_identity(self):
        with temporary_library_directory() as directory:
            spec = replace(assets.MATERIAL_LIBRARY,
                sha256=file_digest(OPTIONS.asset_root/assets.MATERIAL_LIBRARY.path), self_contained=True)
            catalog = assets.material_catalog()
            catalog.register(spec, replace=True)
            loader = BlenderAssetLoader(OPTIONS.asset_root, catalog, policy='RETAINED', linked_root=directory)
            loaded = loader.load_many((assets.PAVEMENT, assets.GRASS, assets.DIAGNOSTICS))
            self.assertIs(loader.load(assets.PAVEMENT), loaded[assets.PAVEMENT])
            self.assertIsNotNone(loaded[assets.PAVEMENT].library)
            for material in set(loaded.values()):
                bpy.data.materials.remove(material)
            for library in list(bpy.data.libraries):
                if str(Path(directory).resolve()) in library.filepath:
                    bpy.data.libraries.remove(library)


if __name__ == '__main__':
    if OPTIONS.reopen:
        bpy.ops.wm.open_mainfile(filepath=str(OPTIONS.reopen.resolve()), use_scripts=False, load_ui=False)
        outputs = [o for o in bpy.data.objects if o.get('st_py_key') and o.type == 'MESH']
        assert outputs and all(o.data.polygons and o.data.materials for o in outputs)
        assert {m.name.split('.')[0] for o in outputs for m in o.data.materials} == {assets.PAVEMENT.name, assets.GRASS.name}
        assert not bpy.data.node_groups and not bpy.data.libraries
        print('SHARED STREET MATERIALS REOPEN PASSED', flush=True)
    else:
        unittest.main(argv=[sys.argv[0]], verbosity=2)
