"""Parsed-program reuse/invalidation; Python-only tests plus Blender integration.

Run directly with Python -B, or Blender --background --factory-startup
--python-exit-code 1 --python this_file. Fixtures stay under workspace/artifacts.
"""
import copy
import importlib
import json
import os
from pathlib import Path
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
PREFIX = '_blosm_program_cache_tests'
for suffix, path in (('', ROOT), ('.app', ROOT / 'app'),
                     ('.building', ROOT / 'building'),
                     ('.building.padl_integration', ROOT / 'building/padl_integration')):
    module = ModuleType(PREFIX + suffix)
    module.__path__ = [str(path)]
    module.__file__ = str(path / '__init__.py')
    sys.modules[PREFIX + suffix] = module
assets = importlib.import_module(PREFIX + '.building.padl_integration.assets')
parser_module = importlib.import_module(PREFIX + '.padl.parser')
evaluation = importlib.import_module(PREFIX + '.padl.evaluation')
PadlParser, PadlError = parser_module.PadlParser, assets.PadlError
try:
    import bpy
except ImportError:
    bpy = None


class ProgramCacheTests(unittest.TestCase):
    def setUp(self):
        artifact_root = (ROOT.parent / 'artifacts/padl_program_cache/unit').resolve()
        artifact_root.mkdir(parents=True, exist_ok=True)
        self.tmp = tempfile.TemporaryDirectory(prefix='library-', dir=artifact_root)
        self.root = Path(self.tmp.name).resolve()
        assert self.root.is_relative_to(artifact_root)
        self.addCleanup(self.tmp.cleanup)
        self.document = {'schema_version': '0.1', 'definitions': {'parameters': {
            'wall_height': {'type': 'uniform', 'min': 5, 'max': 10}}},
            'operations': [{'type': 'generate_volume', 'params': {'roof': {'type': 'flat'},
                'wall_min_height': {'type': 'parameter', 'scope': 'footprint', 'name': 'wall_height'}}}]}
        self.path = self.root / 'building.json'
        self.path.write_text(json.dumps(self.document), encoding='utf-8')
        self.fingerprint = PadlParser().parse(self.document).fingerprint
        self.record = {'document': self.path.name, 'program_sha256': self.fingerprint}
        catalog = {'schema': 'blosm.buildings.v1', 'programs': {'simple_3d': self.record}}
        (self.root / 'prochitecture_buildings.runtime.json').write_text(json.dumps(catalog), encoding='utf-8')
        self.library = assets.BuildingProgramLibrary(self.root)
        loader = patch.object(PadlParser, 'load', autospec=True, side_effect=PadlParser.load)
        self.loader = loader.start()
        self.addCleanup(loader.stop)

    def program(self):
        return self.library.program('simple_3d')

    def test_all_objects_share_one_parse_within_a_library(self):
        program = self.program()
        for _ in range(100):
            self.assertIs(self.program(), program)
        self.assertEqual(self.loader.call_count, 1)
        self.assertIsNot(assets.BuildingProgramLibrary(self.root).program('simple_3d'), program)
        self.assertEqual(self.loader.call_count, 2)

    def test_each_role_keeps_its_own_document_and_diagnostics(self):
        second = self.root / 'realistic.json'
        second.write_text(self.path.read_text(encoding='utf-8'), encoding='utf-8')
        self.library.catalog['programs']['realistic'] = {**self.record, 'document': second.name}
        first, other = self.program(), self.library.program('realistic')
        self.assertEqual(other.fingerprint, first.fingerprint)
        self.assertNotEqual(other.source.filename, first.source.filename)
        self.assertIs(self.library.program('realistic'), other)
        self.assertIs(self.program(), first)
        self.assertEqual(self.loader.call_count, 2)

    def test_unchanged_semantics_with_new_text_refreshes_source_map(self):
        first = self.program()
        self.path.write_text('\n' + self.path.read_text(encoding='utf-8'), encoding='utf-8')
        other = self.program()
        self.assertIsNot(first, other)
        self.assertEqual(first.fingerprint, other.fingerprint)
        self.assertEqual(other.source.at(()).line, first.source.at(()).line + 1)

    def test_same_size_edit_still_revalidates_the_fingerprint(self):
        self.program()
        before = self.path.stat()
        content = self.path.read_text(encoding='utf-8').replace('"min": 5', '"min": 6')
        self.path.write_text(content, encoding='utf-8')
        os.utime(self.path, ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000_000))
        self.assertEqual(self.path.stat().st_size, before.st_size)
        with self.assertRaises(PadlError) as caught:
            self.program()
        self.assertEqual(caught.exception.code, 'program_mismatch')
        self.assertEqual(self.loader.call_count, 2)

    def test_atomic_replacement_with_preserved_size_and_mtime_is_not_stale(self):
        self.program()
        before = self.path.stat()
        replacement = self.root / 'replacement.json'
        replacement.write_text(self.path.read_text(encoding='utf-8').replace('"min": 5', '"min": 6'), encoding='utf-8')
        os.utime(replacement, ns=(before.st_atime_ns, before.st_mtime_ns))
        replacement.replace(self.path)
        self.assertEqual(self.path.stat().st_mtime_ns, before.st_mtime_ns)
        with self.assertRaises(PadlError) as caught:
            self.program()
        self.assertEqual(caught.exception.code, 'program_mismatch')

    def test_catalog_revision_is_checked_on_a_hit(self):
        self.program()
        self.library.catalog['programs']['simple_3d']['program_sha256'] = '0' * 64
        with self.assertRaises(PadlError) as caught:
            self.program()
        self.assertEqual(caught.exception.code, 'program_mismatch')

    def test_publishing_a_new_revision_returns_the_new_program(self):
        first = self.program()
        changed = copy.deepcopy(self.document)
        changed['definitions']['parameters']['wall_height']['min'] = 6
        self.path.write_text(json.dumps(changed), encoding='utf-8')
        fingerprint = PadlParser().parse(changed).fingerprint
        self.library.catalog['programs']['simple_3d']['program_sha256'] = fingerprint
        other = self.program()
        self.assertIsNot(first, other)
        self.assertEqual(other.fingerprint, fingerprint)
        self.assertIs(self.program(), other)

    def test_failed_parse_does_not_poison_the_cache(self):
        self.program()
        self.path.write_text('{', encoding='utf-8')
        with self.assertRaises(PadlError):
            self.program()
        self.path.write_text(json.dumps(self.document), encoding='utf-8')
        self.assertEqual(self.program().fingerprint, self.fingerprint)
        self.assertEqual(self.loader.call_count, 3)

    def test_role_access_checks_are_not_bypassed_by_cache(self):
        self.library.app = SimpleNamespace(isPro=True)
        self.library.catalog['programs']['simple_3d']['edition'] = 'PRO'
        self.program()
        self.library.app.isPro = False
        with self.assertRaises(PadlError) as caught:
            self.program()
        self.assertEqual(caught.exception.code, 'unavailable_program')
        self.library.app.isPro = True
        self.library.catalog['programs']['simple_3d']['document'] = '../outside.json'
        with self.assertRaises(PadlError) as caught:
            self.program()
        self.assertEqual(caught.exception.code, 'incompatible_asset')

    def test_explicit_clear_for_one_or_all_roles(self):
        first = self.program()
        self.library.catalog['programs']['other'] = dict(self.record)
        other = self.library.program('other')
        self.library.clear_program_cache('simple_3d')
        self.assertIsNot(self.program(), first)
        self.assertIs(self.library.program('other'), other)
        self.library.clear_program_cache()
        self.assertIsNot(self.library.program('other'), other)

    def test_cached_structure_retains_no_inputs_or_random_samples(self):
        program = self.program()
        def run(identity, inputs, seed):
            return evaluation.Evaluator(self.program(), seed=seed).evaluate(
                evaluation.FootprintContext(identity, inputs))
        self.assertEqual(run('a', {}, 17), run('a', {}, 17))
        self.assertNotEqual(run('a', {}, 17).samples, run('b', {}, 17).samples)
        self.assertNotEqual(run('a', {}, 17).samples, run('a', {}, 19).samples)
        for height in (6, 9, 7):
            self.assertEqual(run('a', {'wall_height': height}, 17).volumes[0].wall_min_height, height)
        self.assertIs(self.program(), program)
        self.assertEqual(self.loader.call_count, 1)


@unittest.skipUnless(bpy is not None, 'Requires Blender')
class BackendIntegrationTests(unittest.TestCase):
    def test_python_and_gn_reuse_one_file_parse_across_generation_calls(self):
        backend_module = importlib.import_module(PREFIX + '.building.padl_integration.backends')
        geometry = importlib.import_module(PREFIX + '.padl.geometry')
        for backend in ('python', 'geometry_nodes'):
            with self.subTest(backend=backend), patch.object(
                    PadlParser, 'load', autospec=True, side_effect=PadlParser.load) as loader:
                service = backend_module.BuildingGenerationService(backend=backend)
                program = service.library.program(service.role)  # Import-manager prevalidation.
                first = geometry.Footprint('cache:1', 1, ((0, 0), (12, 0), (12, 8), (0, 8)),
                    attributes={'building_use': 'house', 'roof_type': 'flat', 'levels_wall_count': 2})
                second = geometry.Footprint('cache:2', 2, ((20, 0), (32, 0), (32, 8), (20, 8)),
                    attributes={'building_use': 'office', 'roof_type': 'skillion', 'levels_wall_count': 3})
                with service.batch():
                    for group in ((first, second), (first,), (second,)):
                        result = service.generate(group, collection=bpy.context.collection, seed=37)
                        self.assertTrue(result.statuses)
                        self.assertTrue(all(not row['error'] for row in result.statuses), result.statuses)
                self.assertIs(service.library.program(service.role), program)
                self.assertEqual(loader.call_count, 1)


if __name__ == '__main__':
    unittest.main(argv=[sys.argv[0]], verbosity=2)
