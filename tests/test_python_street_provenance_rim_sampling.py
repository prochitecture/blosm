"""Derived sampling spends only documented remaining source-curve budgets."""
from copy import copy
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_one_sided_paired_junctions import fixture
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.one_sided_paired_junctions import OneSidedPairedSpecification, certify
from tested_python_streets.provenance_rim_sampling import BudgetedSourceCurb, resample, grading_retry, grading_retry_steps
from tested_python_streets.partition import QuadChartGradingError
from tested_python_streets.one_sided_charts import OneSidedQuadChart
from tested_python_streets.sampling import AnalyticalBoundarySampler


def profiles(group):
    sources = {source.region.owner: source.layout for source in group.footprints}
    result = {}
    for loop in group.loops:
        data = BoundaryConditions().build(group, loop, sources)
        for port in data.ports:
            approach = sources[port.owner].approaches[port.approach_index]
            result[id(approach)] = loop.xy[list(port.vertices)].copy()
    return result


class ProvenanceRimSamplingTests(unittest.TestCase):
    def test_mouths_real_holes_and_original_source_references_survive(self):
        for dip, holes in ((.1, 0), (2.2, 1)):
            context, site, layout, visits, opposite = fixture(middle_dip=dip, near_bend=0.)
            original = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0).group
            before = profiles(original)
            source = context.network.to_dict()
            sampled = resample(original, context)
            after = profiles(sampled)
            self.assertEqual(set(before), set(after))
            for key in before:
                np.testing.assert_array_equal(after[key], before[key])
            self.assertEqual(sum(loop.hole for loop in sampled.loops), holes)
            self.assertEqual(certify(sampled)['holes'], holes)
            self.assertEqual(context.network.to_dict(), source)
            self.assertEqual(sampled.footprints[1:], original.footprints)
            ledger = sampled.physical_layout.provenance_rim_sampling
            self.assertLessEqual(max(ledger['local_curve_combined_bounds']), context.network.settings.curve_tolerance)
            self.assertGreater(certify(sampled)['boundary_band'], 0.)
            AnalyticalBoundarySampler().build(sampled)

    def test_original_circle_error_is_not_misreported_as_the_strip_chord_error(self):
        context, site, layout, visits, opposite = fixture()
        original = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0).group
        sampled = resample(original, context)
        ledger = sampled.physical_layout.provenance_rim_sampling
        self.assertGreater(max(ledger['original_local_chord_bounds']), ledger['original_strip_chord_error'])
        self.assertLessEqual(max(ledger['local_curve_combined_bounds']), ledger['authored_curve_budget'])

    def test_exhausted_budget_and_repeated_spending_reject_without_moving_input(self):
        context, site, layout, visits, opposite = fixture()
        original = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0).group
        before = original.loops[0].xy.copy()
        sampled = resample(original, context)
        with self.assertRaisesRegex(ValueError, 'twice'):
            resample(sampled, context)
        context.network = replace(context.network, settings=replace(context.network.settings, curve_tolerance=1e-7))
        with self.assertRaisesRegex(ValueError, 'remaining authored'):
            resample(original, context)
        np.testing.assert_array_equal(original.loops[0].xy, before)
        with self.assertRaises(ValueError):
            BudgetedSourceCurb(np.array(((0., 0.), (1., 0.))), (), (), float('nan'))

    def test_retry_requires_typed_grading_and_does_not_hide_source_rejection(self):
        context, site, layout, visits, opposite = fixture()
        group = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0).group
        calls = []
        class Chart:
            def build(self, value, context):
                calls.append(value)
                if len(calls) == 1:
                    raise QuadChartGradingError('positive atlas exhausted its grading budget')
                return SimpleNamespace(metadata={}, group=value)
        result = grading_retry(Chart, group, context)
        self.assertEqual(len(calls), 2)
        self.assertTrue(result.metadata['provenance_rim_sampling'])
        self.assertFalse(getattr(calls[0].physical_layout, 'provenance_rim_sampling', None))
        class InvalidSource:
            def build(self, value, context):
                raise ValueError('actual original source is not contained')
        with self.assertRaisesRegex(ValueError, 'not contained'):
            grading_retry(InvalidSource, group, context)

    def test_complete_opt_in_source_chart_keeps_the_unchanged_quality_contract(self):
        context, site, layout, visits, opposite = fixture()
        group = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0).group
        result = OneSidedQuadChart(context.patterns).build(group, context)
        self.assertLessEqual(result.metadata['stored_quality']['maximum_edge_ratio'], 8.)
        self.assertGreaterEqual(result.metadata['stored_quality']['minimum_corner_sine'], .03)

    def test_cancel_during_chart_build_closes_the_active_candidate(self):
        context, site, layout, visits, opposite = fixture()
        group = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0).group
        events = []

        class Chart:
            def build_steps(self, value, context):
                try:
                    yield 'chart started'
                    events.append('finished')
                finally:
                    events.append('closed')

        steps = grading_retry_steps(Chart, group, context)
        self.assertTrue(any(event == 'chart started' for event in steps))
        steps.close()
        self.assertEqual(events, ['closed'])


if __name__ == '__main__':
    unittest.main()
