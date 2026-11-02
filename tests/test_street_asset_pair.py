"""Load both street GN libraries and shared materials through the asset service.

Run in a fresh background Blender. Optional candidate paths/revisions exercise
the exact pair before installation; with no overrides, use the real bundled
declarations. Loading must preserve source identity, reuse and scene isolation.
"""
import argparse
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from blosm.app.asset_library import BlenderAssetLoader
from blosm.app.library_store import file_digest
from blosm.way.street_network import assets
from blosm.way.street_network.patch_programs import OriginalPatchCatalog, capabilities
from blosm.way.gn.editable_input.patch_program_selection import KINDS


OPTIONS = SimpleNamespace(v8_library=None, v8_revision=None,
                          original_library=None, original_revision=None)


def inventory():
    return {kind: len(getattr(bpy.data, kind)) for kind in
            ('objects', 'meshes', 'curves', 'collections', 'node_groups',
             'materials', 'texts', 'libraries')}


def candidate(spec, path, revision):
    if bool(path) != bool(revision):
        raise ValueError('A candidate library path and revision must be supplied together')
    if path is None:
        return spec
    path = Path(path).resolve()
    return replace(spec, path=str(path), source='external', revision=revision,
                   sha256=file_digest(path))


class StreetAssetPairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not bpy.app.background:
            raise RuntimeError('Run this regression in a fresh background Blender')
        cls.previous = assets.LIBRARY, assets.ORIGINAL_LIBRARY
        assets.LIBRARY = candidate(assets.LIBRARY, OPTIONS.v8_library, OPTIONS.v8_revision)
        assets.ORIGINAL_LIBRARY = candidate(assets.ORIGINAL_LIBRARY,
                                             OPTIONS.original_library, OPTIONS.original_revision)
        cls.loader = BlenderAssetLoader(Path(__file__).resolve().parents[1]/'assets', assets.catalog())

    @classmethod
    def tearDownClass(cls):
        assets.LIBRARY, assets.ORIGINAL_LIBRARY = cls.previous

    def test_pair_preserves_identity_interfaces_catalog_and_reuse(self):
        scene_objects = {obj.as_pointer() for obj in bpy.context.scene.objects}
        references = (assets.SOLVER, assets.PAVEMENT, assets.GRASS, assets.V8_DIAGNOSTICS,
                      assets.ORIGINAL_SOLVER, assets.ORIGINAL_PATCH_PROVIDER,
                      assets.ORIGINAL_STRIP_PROVIDER, assets.ORIGINAL_FAILURE_DISPLAY,
                      assets.ORIGINAL_PAVEMENT, assets.ORIGINAL_GRASS, assets.ORIGINAL_DIAGNOSTICS)
        loaded = self.loader.load_many(references)
        self.assertEqual(set(loaded), set(references))
        self.assertEqual(scene_objects, {obj.as_pointer() for obj in bpy.context.scene.objects})
        self.assertTrue(all(block.library is None for block in loaded.values()))
        # All geometry backends share the same independent material assets.
        self.assertEqual(assets.PAVEMENT.name, assets.ORIGINAL_PAVEMENT.name)
        self.assertEqual(assets.GRASS.name, assets.ORIGINAL_GRASS.name)
        self.assertEqual(loaded[assets.PAVEMENT], loaded[assets.ORIGINAL_PAVEMENT])
        self.assertEqual(loaded[assets.GRASS], loaded[assets.ORIGINAL_GRASS])
        self.assertEqual(assets.PAVEMENT.library, assets.MATERIAL_LIBRARY.key)
        self.assertNotEqual(assets.PAVEMENT.library, assets.SOLVER.library)
        before_validation = inventory()
        assets.validate_solver(loaded[assets.SOLVER])
        root = loaded[assets.ORIGINAL_SOLVER]
        patch = loaded[assets.ORIGINAL_PATCH_PROVIDER]
        strip = loaded[assets.ORIGINAL_STRIP_PROVIDER]
        assets.validate_original_solver(root, cached=True, cached_strips=True, separators=True)
        assets.validate_original_patch_provider(patch)
        assets.validate_original_strip_provider(strip)
        assets.validate_original_failure_display(loaded[assets.ORIGINAL_FAILURE_DISPLAY])
        assets.validate_original_internal_ends(root, patch)
        assets.validate_original_shared_owners(root, patch)
        assets.validate_original_seams(root, patch, strip)
        assets.validate_original_category_connections(root, patch, strip)
        assets.validate_original_quad_energy(root, patch)
        assets.validate_original_failure_support(root, strip)
        self.assertEqual(before_validation, inventory(), 'Validation generated Blender data')
        self.assertEqual(loaded[assets.SOLVER].get('blosm_street_asset_version'), assets.LIBRARY.revision)
        for ref in (assets.ORIGINAL_SOLVER, assets.ORIGINAL_PATCH_PROVIDER,
                    assets.ORIGINAL_STRIP_PROVIDER, assets.ORIGINAL_FAILURE_DISPLAY):
            self.assertEqual(loaded[ref].get('blosm_street_asset_version'), assets.ORIGINAL_LIBRARY.revision)

        descriptor = OriginalPatchCatalog.load_optional(self.loader)
        self.assertEqual(descriptor is not None, bool(root.get('blosm_finite_owner_programs')))
        program_refs = {}
        programs = {}
        if descriptor is not None:
            descriptor.validate_solver(root)
            program_refs = descriptor.references(set(KINDS))
            programs = self.loader.load_many(program_refs.values())
            self.assertEqual(set(program_refs), set(KINDS))
            for ref, tree in programs.items():
                self.assertEqual(capabilities(tree), dict(descriptor.patch))
                assets.validate_original_patch_provider(tree)
                self.assertTrue(tree.get(BlenderAssetLoader.provenance_key))

        stable_counts = inventory()
        repeated = self.loader.load_many(references)
        fresh = BlenderAssetLoader(Path(__file__).resolve().parents[1]/'assets', assets.catalog())
        reloaded = fresh.load_many(references)
        for ref in references:
            self.assertEqual(loaded[ref], repeated[ref])
            self.assertEqual(loaded[ref], reloaded[ref])
            self.assertTrue(loaded[ref].get(BlenderAssetLoader.provenance_key))
        if descriptor is not None:
            same_descriptor = OriginalPatchCatalog.load_optional(fresh)
            self.assertEqual(dict(descriptor.programs), dict(same_descriptor.programs))
            repeated_programs = fresh.load_many(program_refs.values())
            for ref in programs:
                self.assertEqual(programs[ref], repeated_programs[ref])
        self.assertEqual(stable_counts, inventory(), 'A fresh loader duplicated existing asset data')
        self.assertEqual(scene_objects, {obj.as_pointer() for obj in bpy.context.scene.objects})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--v8-library', type=Path)
    parser.add_argument('--v8-revision')
    parser.add_argument('--original-library', type=Path)
    parser.add_argument('--original-revision')
    OPTIONS = parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    unittest.main(argv=[sys.argv[0]], verbosity=2)
