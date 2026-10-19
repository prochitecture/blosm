"""Pure-Python evidence bookkeeping checks; no Blender geometry is substituted."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest

path = Path(__file__).resolve().parents[1]/'way/street_network/coverage.py'
spec = importlib.util.spec_from_file_location('street_coverage_under_test',path)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
CoverageReport = module.CoverageReport


class JunctionSelectionTests(unittest.TestCase):
    def snapshot(self):
        return SimpleNamespace(intersections={i:dict(id=i,surface_status='GENERATED') for i in (2,7,1001)},
                               streets={},bundles={})

    def test_explicit_physical_and_logical_sites_are_pending_not_generated(self):
        snapshot = self.snapshot()
        CoverageReport.initialize_junctions(snapshot,{2,1001})
        self.assertEqual([snapshot.intersections[i]['surface_status'] for i in (2,7,1001)],
                         ['UNSUPPORTED','DELIBERATELY_EXCLUDED','UNSUPPORTED'])

    def test_full_input_does_not_exclude_or_reuse_previous_surface_success(self):
        snapshot = self.snapshot()
        CoverageReport.initialize_junctions(snapshot)
        self.assertEqual({r['surface_status'] for r in snapshot.intersections.values()},{'UNSUPPORTED'})

    def test_empty_explicit_selection_differs_from_full_input(self):
        snapshot = self.snapshot()
        CoverageReport.initialize_junctions(snapshot,set())
        self.assertEqual({r['surface_status'] for r in snapshot.intersections.values()},
                         {'DELIBERATELY_EXCLUDED'})

    def test_report_keeps_excluded_failed_and_generated_counts_distinct(self):
        snapshot = self.snapshot()
        CoverageReport.initialize_junctions(snapshot,{2,1001})
        snapshot.intersections[2]['surface_status'] = 'GENERATED'
        snapshot.intersections[1001]['surface_status'] = 'FAILED'
        result = CoverageReport().build(snapshot)
        self.assertEqual(result['junction_counts'],
                         {'GENERATED':1,'DELIBERATELY_EXCLUDED':1,'FAILED':1})


if __name__=='__main__':
    unittest.main()
