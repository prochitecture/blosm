"""Operator errors must reject invalid settings without adding scene objects."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import unittest
import importlib
from types import SimpleNamespace
from unittest.mock import patch
import bpy
import addon_utils
from release_support import arguments, enable, REPO, PACKAGE, catalog, configure_osm, reset, run

args = arguments()
enable(args)


class ValidationChecks(unittest.TestCase):
    def test_osm_requires_existing_data_folder_in_every_mode(self):
        prefs = bpy.context.preferences.addons[PACKAGE].preferences
        old = prefs.dataDir
        addon = importlib.import_module(PACKAGE)
        app = addon.blenderApp.app
        not_folder = args.artifacts / 'data-is-a-file.txt'
        not_folder.write_text('not a directory', encoding='utf-8')
        modes = ['2D', '3Dsimple'] + (['3Drealistic'] if app.isPro else [])
        try:
            for mode in modes:
                for source in ('file', 'server'):
                    for invalid in ('', '   ', str(args.artifacts / 'missing-data-folder'), str(not_folder)):
                        with self.subTest(mode=mode, source=source, directory=invalid):
                            settings = configure_osm(args, next(e for e in catalog() if e['id'] == 'rectangle'), mode=mode)
                            settings.osmSource = source
                            prefs.dataDir = invalid
                            with patch.object(app, 'downloadOsmFile') as download, patch.object(
                                    bpy.utils, 'extension_path_user', side_effect=AssertionError('Unexpected storage fallback')):
                                self.rejected('directory')
                                download.assert_not_called()
            with self.assertRaisesRegex(Exception, 'directory'):
                app.setDataDir(SimpleNamespace(preferences=SimpleNamespace(addons={})),
                               str(REPO), PACKAGE, required=True)
        finally:
            prefs.dataDir = old

    def test_osm_download_uses_configured_folder(self):
        addon = importlib.import_module(PACKAGE)
        prefs = bpy.context.preferences.addons[PACKAGE].preferences
        old = prefs.dataDir
        folder = args.artifacts / 'chosen-data-folder'
        folder.mkdir(exist_ok=True)
        settings = configure_osm(args, next(e for e in catalog() if e['id'] == 'rectangle'), mode='2D')
        fixture = Path(settings.osmFilepath).read_bytes()
        settings.osmSource = 'server'
        destinations = []
        def download(directory, *bounds):
            destination = Path(directory) / 'download.osm'
            destinations.append(destination)
            destination.write_bytes(fixture)
            addon.blenderApp.app.osmFilepath = str(destination)
        try:
            prefs.dataDir = str(folder.resolve())
            with patch.object(addon.blenderApp.app, 'downloadOsmFile', side_effect=download):
                self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
            self.assertEqual(destinations, [folder / 'osm' / 'download.osm'])
        finally:
            prefs.dataDir = old

    def rejected(self, message):
        before = set(bpy.data.objects.keys())
        try:
            result = bpy.ops.blosm.import_data()
        except RuntimeError as exc:
            self.assertRegex(str(exc).lower(), message.lower())
        else:
            self.assertEqual(result, {'CANCELLED'})
        self.assertEqual(set(bpy.data.objects.keys()), before)

    def test_missing_osm(self):
        settings = configure_osm(args, catalog()[0])
        settings.osmFilepath = str(args.artifacts / 'does_not_exist.osm')
        self.rejected('exist|no such file')

    def test_missing_gpx(self):
        settings = reset(args)
        settings.dataType = 'gpx'
        settings.gpxFilepath = str(args.artifacts / 'does_not_exist.gpx')
        self.rejected('valid gpx')

    def test_missing_setup_script(self):
        settings = configure_osm(args, catalog()[0])
        settings.setupScript = str(args.artifacts / 'does_not_exist.py')
        self.rejected('exist')

    def test_missing_data_directory(self):
        configure_osm(args, catalog()[0])
        prefs = bpy.context.preferences.addons[PACKAGE].preferences
        old = prefs.dataDir
        try:
            prefs.dataDir = str(args.artifacts / 'does_not_exist')
            self.rejected('directory')
        finally:
            prefs.dataDir = old

    def test_disable_reenable(self):
        reset(args)
        for _ in range(2):
            addon_utils.disable(PACKAGE, default_set=True)
            self.assertFalse(hasattr(bpy.types.Scene, 'blosm'))
            enable(args)
            self.assertTrue(hasattr(bpy.types.Scene, 'blosm'))
        configure_osm(args, next(e for e in catalog() if e['id'] == 'rectangle'))
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})


run(args, ValidationChecks)
