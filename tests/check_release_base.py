"""Base-only negative contracts, exercised against the actual installed ZIP."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
from types import SimpleNamespace, MethodType
import unittest
from unittest.mock import patch
import bpy
from release_profiles import PRO_MODULES, PRO_SETUP, PRO_OPERATORS, inspect_candidate
from release_support import (arguments, enable, reset, PACKAGE, PROFILE, ADDON_ROOT,
                             catalog, configure_osm, geometry, meshes, assert_geometry, run)

args = arguments()
addon = enable(args)
gui = importlib.import_module(PACKAGE + '.gui')


class Layout:
    """Record the actual panel draw callbacks without needing screen coordinates."""
    def __init__(self):
        self.properties, self.operators = [], []
    def prop(self, data, name, **kwargs):
        self.properties.append(name)
    def prop_search(self, data, name, *args, **kwargs):
        self.properties.append(name)
    def operator(self, name, **kwargs):
        self.operators.append(name)
        return SimpleNamespace()
    def __getattr__(self, name):
        if name in ('box', 'row', 'column', 'split'):
            return lambda *args, **kwargs: self
        if name in ('label', 'separator', 'template_list'):
            return lambda *args, **kwargs: None
        raise AttributeError(name)


class BaseChecks(unittest.TestCase):
    def setUp(self):
        self.assertEqual(PROFILE, 'base', 'Base-specific tests require --profile base')
        self.assertFalse(addon.blenderApp.app.isPro)
        reset(args)

    def test_pro_modules_and_setup_scripts_absent(self):
        _, names = inspect_candidate(args.candidate_zip)
        for module in PRO_MODULES:
            prefix = module.replace('.', '/')
            self.assertFalse([n for n in names if n == prefix + '.py' or n.startswith(prefix + '/')], module)
            self.assertFalse((ADDON_ROOT / prefix).exists(), module)
            with self.assertRaises(ModuleNotFoundError) as error:
                importlib.import_module(PACKAGE + '.' + module)
            self.assertEqual(error.exception.name, PACKAGE + '.' + module)
        for name in PRO_SETUP:
            self.assertNotIn(f'setup/{name}.py', names)
            with self.assertRaises(ModuleNotFoundError):
                importlib.import_module(PACKAGE + '.setup.' + name)

    def test_common_capabilities_and_realistic_mode_rejected(self):
        app = addon.blenderApp.app
        self.assertFalse(app.has('mode3dRealistic'))
        self.assertTrue(app.has('mode3d'))
        self.assertTrue(app.has('overlay'))
        settings = bpy.context.scene.blosm
        modes = {entry.identifier for entry in settings.bl_rna.properties['mode'].enum_items}
        self.assertEqual(modes, {'2D', '3Dsimple'})
        self.assertEqual({entry[0] for entry in gui.getDataTypes()},
                         {'osm', 'terrain', 'overlay', 'gpx', '3d-tiles'})
        with self.assertRaises(TypeError):
            settings.mode = '3Drealistic'

    def test_pro_operators_unregistered_and_not_callable(self):
        for name in PRO_OPERATORS:
            operator = getattr(bpy.ops.blosm, name)
            with self.assertRaises(KeyError):
                operator.get_rna_type()
            with self.assertRaises(AttributeError):
                operator()
        self.assertIsNotNone(bpy.ops.blosm.replace_materials.get_rna_type())

    def test_pro_controls_absent_from_draw_callbacks(self):
        settings = bpy.context.scene.blosm
        prefs = bpy.context.preferences.addons[PACKAGE].preferences
        prefs.enableExperimentalFeatures = True  # Stale preferences must not reveal Pro controls.
        settings.importForExport = True
        for mode in ('2D', '3Dsimple'):
            settings.mode = mode
            layout = Layout()
            panel = SimpleNamespace(layout=layout)
            for name in ('_drawBuildingSettings', '_drawDefaultLevels'):
                setattr(panel, name, MethodType(getattr(gui.BLOSM_PT_Settings, name), panel))
            gui.BLOSM_PT_Settings.drawOsm(panel, bpy.context)
            self.assertFalse({'importForExport', 'assetPackage', 'litWindows', 'treeDensity'}.intersection(layout.properties))
            self.assertIn('singleObject', layout.properties)
        settings.dataType = '3d-tiles'
        bpy.ops.mesh.primitive_cube_add()
        layout = Layout()
        gui.BLOSM_PT_Tools.draw(SimpleNamespace(layout=layout), bpy.context)
        self.assertNotIn('blosm.bake_to_single_material', layout.operators)
        self.assertFalse(any(p.startswith('bake_') for p in layout.properties))
        self.assertIn('blosm.replace_materials', layout.operators)
        class PreferencesView:
            def __getattr__(self, name):
                return getattr(prefs, name)
        view = PreferencesView()
        view.layout = Layout()
        addon.BlosmPreferences.draw(view, bpy.context)
        self.assertFalse({'assetsDir', 'enableExperimentalFeatures', 'screenType'}.intersection(view.layout.properties))

    def test_google_highest_lod_blocked_before_download(self):
        settings = bpy.context.scene.blosm
        settings.dataType = '3d-tiles'
        settings.threedTilesSource = 'google'
        settings.lodOf3dTiles = 'lod6'
        settings.progressiveImport = settings.concurrentDownloads = False
        entry = settings.bl_rna.properties['lodOf3dTiles'].enum_items['lod6']
        self.assertIn('pro', entry.name.lower())  # Upsell choice may remain visible.
        manager = importlib.import_module(PACKAGE + '.threed_tiles.manager')
        with patch.object(bpy.types.WindowManager, 'popup_menu') as popup, \
                patch.object(manager.BaseManager, 'download') as download:
            self.assertEqual(bpy.ops.blosm.import_data(), {'CANCELLED'})
            download.assert_not_called()
            popup.assert_called_once()
        self.assertFalse(meshes())

    def test_stale_export_flags_do_not_enable_pro_renderer(self):
        entry = next(e for e in catalog() if e['id'] == 'rectangle')
        snapshots = []
        for export in (False, True):
            settings = configure_osm(args, entry)
            settings.importForExport = export
            bpy.context.preferences.addons[PACKAGE].preferences.enableExperimentalFeatures = export
            self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
            self.assertFalse(addon.blenderApp.app.enableExperimentalFeatures)
            snapshots.append(geometry(meshes()))
        assert_geometry(self, *snapshots)


run(args, BaseChecks)
