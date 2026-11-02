"""Native addon registration, retired backend migration and provider extension.

Run in factory-startup Blender with --python-exit-code 1.
"""
from pathlib import Path
import sys
import unittest

import bpy
from _bpy_restrict_state import RestrictBlend

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
import blosm
from blosm import gui
from blosm.way.street_network.backends import StreetRenderBackend, StreetBackendDefinition, registry
from blosm.way.renderer_streets import StreetRenderer


class PythonStreetBackendTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.observers = tuple(bpy.app.handlers.depsgraph_update_post)
        with RestrictBlend():
            blosm.register()

    @classmethod
    def tearDownClass(cls):
        blosm.unregister()

    def test_python_install_has_no_experimental_street_gn_backend(self):
        self.assertFalse((ROOT/'way/gn').exists())
        self.assertFalse((ROOT/'way/item_renderer').exists())
        self.assertFalse((ROOT/'way/asset_store.py').exists())
        self.assertTrue((ROOT/'app/asset_store.py').is_file())
        for name in ('prochitecture_streets_v8.blend', 'prochitecture_streets_original_spans.blend'):
            self.assertFalse((ROOT/'assets'/name).exists(), name)
        self.assertTrue((ROOT/'assets/prochitecture_street_materials.blend').is_file())
        self.assertFalse(any(name == 'blosm.way.gn' or name.startswith('blosm.way.gn.') for name in sys.modules))
        self.assertEqual(tuple(bpy.app.handlers.depsgraph_update_post), self.observers)

    def test_gui_keeps_python_materials_and_building_gn(self):
        props = gui.BlosmProperties.bl_rna.properties
        modes = {key:number for key, label, description, number in gui.getStreetBackends(None, None)}
        self.assertEqual(modes, {'python':0})
        self.assertEqual(bpy.context.scene.blosm.streetBackend, 'python')
        for retired in ('streetConstruction', 'streetAssetFile'):
            self.assertNotIn(retired, props)
        self.assertIn('streetMaterialFile', props)
        self.assertIn('geometry_nodes', {item.identifier for item in props['buildingBackend'].enum_items})
        from blosm.building.padl_integration.backends import BuildingGenerationService
        from blosm.building.padl_integration.gn_runtime import read_status
        self.assertTrue(callable(BuildingGenerationService))
        self.assertTrue(callable(read_status))

    def test_saved_legacy_and_old_gn_selections_use_python(self):
        legacy = bpy.data.scenes.new('Saved Legacy street selection')
        retired = bpy.data.scenes.new('Saved retired street selection')
        try:
            legacy.blosm['streetBackend'] = 2
            retired.blosm['streetBackend'] = 1
            gui._onFileLoaded(None)
            self.assertEqual(legacy.blosm.streetBackend, 'python')
            self.assertEqual(retired.blosm.streetBackend, 'python')
            self.assertIn(gui._onFileLoaded, bpy.app.handlers.load_post)
        finally:
            bpy.data.scenes.remove(legacy)
            bpy.data.scenes.remove(retired)

    def test_new_backend_appears_in_real_rna_and_is_used_without_setup_branches(self):
        calls = []
        class FutureGeometryNodesBackend(StreetRenderBackend):
            def __init__(self, app, **options):
                calls.append(('construct', app))
            def prepare(self):
                calls.append('prepare')
            def render(self, manager, data):
                calls.append(('render', manager, data))
                return 'future context'
            def finalize(self):
                calls.append('finalize')
            def cleanup(self):
                calls.append('cleanup')
        scene = bpy.data.scenes.new('Future backend extension')
        registry.register(StreetBackendDefinition('geometry_nodes', 'Future GN', 'Test implementation',
                                                 3, FutureGeometryNodesBackend))
        try:
            self.assertEqual(len(gui.getStreetBackends(None, None)), 2)
            self.assertEqual(calls, [])
            scene.blosm.streetBackend = 'geometry_nodes'
            gui._onFileLoaded(None)
            self.assertEqual(scene.blosm.streetBackend, 'geometry_nodes')
            self.assertEqual(scene.blosm.get('streetBackend'), 3)
            renderer = StreetRenderer(scene.blosm)
            renderer.prepare()
            self.assertEqual(renderer.render('manager', 'osm'), 'future context')
            renderer.finalize()
            renderer.cleanup()
            self.assertEqual(calls, [('construct', scene.blosm), 'prepare',
                                     ('render', 'manager', 'osm'), 'finalize', 'cleanup'])
        finally:
            bpy.data.scenes.remove(scene)
            registry.unregister('geometry_nodes')


if __name__ == '__main__':
    unittest.main(argv=[sys.argv[0]], verbosity=2)
