"""A partial aggregate surface must not qualify as a complete imported network."""
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest


spec = importlib.util.spec_from_file_location('_aggregate_surface_coverage_tests',
    Path(__file__).resolve().parents[1] / 'way/street_network/aggregate_coverage.py')
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
Coverage = module.AggregateSurfaceCoverage


class AggregateCoverageTests(unittest.TestCase):
    def setUp(self):
        def table(field, values):
            return SimpleNamespace(attributes={field: SimpleNamespace(
                data=[SimpleNamespace(value=v) for v in values])})
        plan = SimpleNamespace(ports=SimpleNamespace(sites=[
            SimpleNamespace(junction_id=70), SimpleNamespace(junction_id=90)]))
        tables = {
            'Corridors': table('st_cs_street_id', [10, 20]),
            'External Approaches': table('st_ca_street_id', [30, 10]),
            'Separators': table('st_sb_bundle_id', [800, 800, 900]),
        }
        self.coverage = Coverage.prepare(plan, tables)
        self.complete = {(0, 0), (0, 1), (1, 0), (1, 1),
                         (3, 0), (3, 1), (4, 800), (4, 900)}

    def test_saved_ownership_roundtrip_and_face_repetitions(self):
        restored = Coverage.from_object({
            'blosm_aggregate_surface_owners': self.coverage.encode()})
        self.assertEqual(restored.owners, self.coverage.owners)
        self.assertTrue(restored.evaluate(list(self.complete) * 3)['complete'])

    def test_missing_junction_or_street_reports_original_identity(self):
        missing = self.coverage.evaluate(self.complete - {(0, 1), (3, 1)})
        self.assertFalse(missing['complete'])
        self.assertEqual(missing['missing_junctions'], [90])
        self.assertEqual(missing['missing_streets'], [10])
        # Street 10's corridor exists, but its separate external approach is missing.
        self.assertEqual(missing['missing_owners'], [(0, 1), (3, 1)])

    def test_missing_separator_cannot_be_hidden_by_complete_pavement(self):
        missing = self.coverage.evaluate(self.complete - {(4, 900)})
        self.assertFalse(missing['complete'])
        self.assertEqual(missing['missing_bundles'], [900])
        self.assertEqual(missing['missing_streets'], [])

    def test_wrong_owner_or_kind_and_empty_output_are_incomplete(self):
        wrong = self.coverage.evaluate((self.complete - {(1, 0)}) | {(3, 10)})
        self.assertFalse(wrong['complete'])
        self.assertEqual(wrong['unexpected_owners'], [(3, 10)])
        self.assertEqual(wrong['missing_streets'], [10])
        self.assertFalse(self.coverage.evaluate([])['complete'])


if __name__ == '__main__':
    unittest.main()
