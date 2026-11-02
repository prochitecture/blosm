"""Native addon registration after retiring the experimental street GN backend.

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
        for name in ('prochitecture_streets_v8.blend', 'prochitecture_streets_original_spans.blend'):
            self.assertFalse((ROOT/'assets'/name).exists(), name)
        self.assertTrue((ROOT/'assets/prochitecture_street_materials.blend').is_file())
        self.assertFalse(any(name == 'blosm.way.gn' or name.startswith('blosm.way.gn.') for name in sys.modules))
        self.assertEqual(tuple(bpy.app.handlers.depsgraph_update_post), self.observers)

    def test_gui_keeps_python_legacy_materials_and_building_gn(self):
        props = gui.BlosmProperties.bl_rna.properties
        modes = {item.identifier:item.value for item in props['streetBackend'].enum_items}
        self.assertEqual(modes, {'python':0, 'legacy':2})
        for retired in ('streetConstruction', 'streetAssetFile'):
            self.assertNotIn(retired, props)
        self.assertIn('streetMaterialFile', props)
        self.assertIn('geometry_nodes', {item.identifier for item in props['buildingBackend'].enum_items})
        from blosm.building.padl_integration.backends import BuildingGenerationService
        from blosm.building.padl_integration.gn_runtime import read_status
        self.assertTrue(callable(BuildingGenerationService))
        self.assertTrue(callable(read_status))

    def test_saved_legacy_selection_is_preserved_and_retired_selection_uses_python(self):
        legacy = bpy.data.scenes.new('Saved Legacy street selection')
        retired = bpy.data.scenes.new('Saved retired street selection')
        try:
            legacy.blosm['streetBackend'] = 2
            retired.blosm['streetBackend'] = 1
            gui._onFileLoaded(None)
            self.assertEqual(legacy.blosm.streetBackend, 'legacy')
            self.assertEqual(retired.blosm.streetBackend, 'python')
            self.assertIn(gui._onFileLoaded, bpy.app.handlers.load_post)
        finally:
            bpy.data.scenes.remove(legacy)
            bpy.data.scenes.remove(retired)


if __name__ == '__main__':
    unittest.main(argv=[sys.argv[0]], verbosity=2)
