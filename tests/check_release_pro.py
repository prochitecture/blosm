"""Pro identity checks; functional Pro coverage is reused from the general suites."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
import unittest
import bpy
from release_profiles import PRO_MODULES, PRO_SETUP, PRO_OPERATORS, inspect_candidate
from release_support import arguments, enable, PACKAGE, PROFILE, run

args = arguments()
addon = enable(args)


class ProChecks(unittest.TestCase):
    def test_pro_package_and_runtime_capabilities(self):
        self.assertEqual(PROFILE, 'pro')
        self.assertTrue(addon.blenderApp.app.isPro)
        self.assertTrue(addon.blenderApp.app.has('mode3dRealistic'))
        _, names = inspect_candidate(args.candidate_zip)
        for module in PRO_MODULES:
            self.assertIn(module.replace('.', '/') + '/__init__.py', names)
        for name in PRO_SETUP:
            self.assertIn(f'setup/{name}.py', names)
        for module in ('setup.premium_default', 'setup.realistic_dev', 'pml', 'util.blender_extra.baking'):
            self.assertIsNotNone(importlib.import_module(PACKAGE + '.' + module))
        settings = bpy.context.scene.blosm
        self.assertIn('3Drealistic', settings.bl_rna.properties['mode'].enum_items)
        settings.mode = '3Drealistic'
        self.assertNotIn('[pro]', settings.bl_rna.properties['lodOf3dTiles'].enum_items['lod6'].name.lower())

    def test_pro_operators_registered(self):
        for name in PRO_OPERATORS:
            self.assertIsNotNone(getattr(bpy.ops.blosm, name).get_rna_type(), name)


run(args, ProChecks)
