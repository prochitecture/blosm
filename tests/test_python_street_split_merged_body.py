"""A cell-budget fallback splits computation without moving actual fronts."""
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from test_python_street_common_fibre_bodies import fixture
from tested_python_streets.bundle_merged_body import MergedRunBody, MergedBodyCellBudgetError
from tested_python_streets.split_merged_body import BoundedMergedRunBody, BoundedGappedRunBody


class BoundedMergedBodyTests(unittest.TestCase):
    def test_explicit_separated_computation_welds_pavement_and_whole_grass_fibre(self):
        for separated in (True, 'mixed'):
            frame, field, approach, cuts, context, bundle = fixture(separated,
                angle=.713, translation=(1200., -1800.))
            before = {key: points.copy() for key, points in cuts.items()}
            build = MergedRunBody.build_between
            def exhausted(body, *args, **kwargs):
                if isinstance(body, BoundedGappedRunBody):
                    raise MergedBodyCellBudgetError(9000, 8000)
                return build(body, *args, **kwargs)
            with patch.object(MergedRunBody, 'build_between', exhausted):
                bounded = BoundedGappedRunBody()
                pavement, grass, report = bounded.build_between(context, bundle,
                    approach, cuts, namespace=('gapped', separated))
            self.assertEqual(report['child_cell_budgets'], [8000, 8000])
            self.assertEqual(report['complete_leaf_limit'], 16000)
            self.assertIsNotNone(grass)
            if separated is True:
                self.assertFalse(report['computational_fibre_merged'])
                self.assertGreater(report['exact_shared_grass_edges'], 0)
            lookup = dict(zip(bounded.pavement_keys, pavement.xyz))
            for key in cuts:
                np.testing.assert_array_equal(cuts[key], before[key])
                np.testing.assert_array_equal(np.array([lookup[k] for k in bounded.cut_keys[key]]), cuts[key])
            for row in report['stored_checks']:
                self.assertLessEqual(row['stored_quality']['maximum_edge_ratio'], 8.)
                self.assertEqual(row['drainage']['maximum_fill'], 0.)

    def test_separated_child_grass_partitions_must_match_exactly(self):
        frame, field, approach, cuts, context, bundle = fixture(True)
        build = MergedRunBody.build_between
        def changed(body, *args, **kwargs):
            if isinstance(body, BoundedGappedRunBody):
                raise MergedBodyCellBudgetError(9000, 8000)
            result = build(body, *args, **kwargs)
            if kwargs['namespace'][-1] == 1:
                points, names = body.open_grass_cuts[0]
                points = points.copy();points[len(points)//2, 0] += 1e-8
                body.open_grass_cuts[0] = points, names
            return result
        with patch.object(MergedRunBody, 'build_between', changed), self.assertRaisesRegex(
                ValueError, 'identical complete grass fibre'):
            BoundedGappedRunBody().build_between(context, bundle, approach, cuts, namespace='bad grass')

    def test_default_body_is_unchanged_and_budget_fallback_welds_one_complete_crown(self):
        frame, field, approach, cuts, context, bundle = fixture(angle=.713, translation=(1200., -1800.))
        original = MergedRunBody()
        expected = original.build_between(context, bundle, approach, cuts, namespace='same')
        bounded = BoundedMergedRunBody()
        result = bounded.build_between(context, bundle, approach, cuts, namespace='same')
        np.testing.assert_array_equal(expected[0].xyz, result[0].xyz)
        np.testing.assert_array_equal(expected[0].quads, result[0].quads)
        self.assertEqual(expected[2], result[2])
        source = {key: points.copy() for key, points in cuts.items()}
        build = MergedRunBody.build_between
        def exhausted(body, *args, **kwargs):
            if isinstance(body, BoundedMergedRunBody): raise MergedBodyCellBudgetError(9000, 8000)
            return build(body, *args, **kwargs)
        with patch.object(MergedRunBody, 'build_between', exhausted):
            bounded = BoundedMergedRunBody()
            pavement, grass, report = bounded.build_between(context, bundle, approach, cuts, namespace='split')
        self.assertIsNone(grass)
        self.assertEqual(report['child_cell_budgets'], [8000, 8000])
        self.assertEqual(report['complete_leaf_limit'], 16000)
        self.assertGreater(report['exact_shared_crown_edges'], 0)
        self.assertFalse(report['physical_source_changed'])
        self.assertFalse(report['source_intervals_changed'])
        by_key = dict(zip(bounded.pavement_keys, pavement.xyz))
        for key, points in source.items():
            np.testing.assert_array_equal(cuts[key], points)
            np.testing.assert_array_equal(bounded.cuts[key], points)
            np.testing.assert_array_equal(np.asarray([by_key[name] for name in bounded.cut_keys[key]]), points)
        self.assertLessEqual(report['stored_checks'][0]['stored_quality']['maximum_edge_ratio'], 8.)
        self.assertGreaterEqual(report['stored_checks'][0]['stored_quality']['minimum_corner_sine'], .03)
        self.assertEqual(report['stored_checks'][0]['drainage']['maximum_fill'], 0.)

    def test_actual_separated_fronts_and_nonbudget_errors_do_not_enter_the_new_family(self):
        frame, field, approach, cuts, context, bundle = fixture('mixed')
        build = MergedRunBody.build_between
        def exhausted(body, *args, **kwargs):
            if isinstance(body, BoundedMergedRunBody): raise MergedBodyCellBudgetError(9000, 8000)
            return build(body, *args, **kwargs)
        with patch.object(MergedRunBody, 'build_between', exhausted), self.assertRaisesRegex(ValueError, 'two actual merged'):
            BoundedMergedRunBody().build_between(context, bundle, approach, cuts, namespace='invalid')
        broken = {key: points.copy() for key, points in cuts.items()}
        broken[0, 0] = broken[0, 0][1:]
        with self.assertRaisesRegex(ValueError, 'curb or crown'):
            BoundedMergedRunBody().build_between(context, bundle, approach, broken, namespace='invalid')


if __name__ == '__main__': unittest.main()
