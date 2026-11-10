"""Common checks on the extension installed by run_candidate_suite.py."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import tomllib
import unittest
import bpy
import addon_utils
from release_profiles import inspect_candidate
from release_support import (arguments, enable, ADDON_ROOT, PACKAGE, PROFILE,
                             catalog, configure_osm, run)

args = arguments()
addon = enable(args)


class PackageChecks(unittest.TestCase):
    def test_installed_identity_and_manifest(self):
        self.assertIn(PROFILE, ('base', 'pro'), 'Use a build profile with a candidate ZIP')
        self.assertTrue(ADDON_ROOT.is_relative_to(args.artifacts.resolve()))
        self.assertTrue(PACKAGE.startswith('bl_ext.'))
        manifest = tomllib.loads((ADDON_ROOT / 'blender_manifest.toml').read_text(encoding='utf-8'))
        self.assertEqual(manifest['id'], 'blosm')
        self.assertEqual(addon.BlosmPreferences.bl_idname, PACKAGE)
        self.assertEqual(addon.blenderApp.app.version, tuple(map(int, manifest['version'].split('.'))))
        self.assertGreaterEqual(bpy.app.version, tuple(map(int, manifest['blender_version_min'].split('.'))))

    def test_package_excludes_development_inputs(self):
        _, names = inspect_candidate(args.candidate_zip)
        self.assertFalse([n for n in names if n.startswith(('tests/', 'docs/', 'tools/', 'script/'))])
        self.assertFalse([n for n in names if '__pycache__' in n or n.endswith(('_access_token.txt', 'api_key.txt'))])
        self.assertEqual({n for n in names if n.endswith('.blend')},
                         {'assets/base.blend', 'assets/roofs.blend', 'assets/way_profiles.blend'})

    def test_osm_after_disable_reenable(self):
        addon_utils.disable(PACKAGE, default_set=True)
        self.assertFalse(hasattr(bpy.types.Scene, 'blosm'))
        enable(args)
        configure_osm(args, next(e for e in catalog() if e['id'] == 'rectangle'))
        self.assertEqual(bpy.ops.blosm.import_data(), {'FINISHED'})
        self.assertTrue([o for o in bpy.context.scene.objects if o.type == 'MESH' and o.data.polygons])


run(args, PackageChecks)
