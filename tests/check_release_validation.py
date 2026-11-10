"""Operator errors must reject invalid settings without adding scene objects."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import unittest
import bpy
import addon_utils
from release_support import arguments, enable, REPO, PACKAGE, catalog, configure_osm, reset, run

args = arguments()
enable(args)


class ValidationChecks(unittest.TestCase):
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
