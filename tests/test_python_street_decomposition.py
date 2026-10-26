import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.compounds import CompoundGroup, RimSegment, SampledFootprint
from tested_python_streets.decomposition import VisibilitySplitter
from tested_python_streets.partition import CompoundPartitioner
from tested_python_streets.quad_charts import visibility_center
from tested_python_streets.regions import Region, RegionUnion, signed_area
from tested_python_streets.bundles import CenterlineRange, certify_disjoint
from tested_python_streets.paths import Path


U = np.array(((0., 0.), (4., 0.), (4., 4.), (3., 4.), (3., 1.), (1., 1.), (1., 4.), (0., 4.)))


class PythonVisibilityDecompositionTests(unittest.TestCase):
    def test_indented_region_splits_into_visible_blocks_and_retains_boundary(self):
        owner = model.Key('intersection', 1)
        footprint = SampledFootprint(Region(owner, U), tuple(RimSegment('curb', i, 0., 1.) for i in range(len(U))))
        group = CompoundGroup((owner,), (), (footprint,), RegionUnion().build((footprint.region,)))
        result = CompoundPartitioner().build(group, None)
        self.assertEqual(len(result.blocks), 2)
        self.assertEqual(len(result.cuts), 1)
        self.assertIsNone(result.cuts[0].section_id, 'Computational cuts must not invent a source Section identity')
        self.assertAlmostEqual(sum(block.area for block in result.blocks), signed_area(U))
        for block in result.blocks:
            self.assertGreaterEqual(len(block.vertices), 4)
            self.assertTrue(np.all(np.isfinite(visibility_center(result.xy[list(block.vertices)]))))
        for point in U:
            self.assertTrue(np.any(np.all(result.xy == point, axis=1)))
        sampled = result.sample(.7)
        self.assertTrue(all(len(rim) % 2 == 0 for rim in sampled.rims))

    def test_translation_and_warm_runs_keep_the_same_computational_cuts(self):
        splitter = VisibilitySplitter()
        expected = splitter.split(U)
        self.assertEqual(expected, splitter.split(U))
        self.assertEqual(expected, splitter.split(U + (1_000_000., -2_000_000.)))
        for first, last in expected:
            self.assertTrue(splitter.visible(U, first, last))
            self.assertGreater(splitter.cut_angle(U, first, last), np.radians(5.))

    def test_visible_polygon_needs_no_cut_and_invalid_budgets_are_rejected(self):
        self.assertEqual(VisibilitySplitter().split(np.array(((0., 0.), (3., 0.), (3., 2.), (0., 2.)))), ())
        for settings in ({'max_cuts': 0}, {'max_candidates': 0}, {'epsilon': float('nan')}):
            with self.assertRaises(ValueError):
                VisibilitySplitter(**settings)
        with self.assertRaises(ValueError):
            VisibilitySplitter().split(U[::-1])

    def test_centerline_range_matches_the_original_cubic_and_metric_endpoints(self):
        path = Path(model.Section(1, 1, ((0., 0., 0.), (20., 2., 0.), (40., 0., 0.))))
        first, last = 5., path.length - 7.
        cropped = CenterlineRange.build(path, first, last)
        self.assertEqual(len(cropped.curves), 2)
        np.testing.assert_allclose(cropped.curves[0][0], path.frame(first)[0], atol=1e-10)
        np.testing.assert_allclose(cropped.curves[-1][-1], path.frame(last)[0], atol=1e-10)
        np.testing.assert_allclose(cropped.curves[0][-1], path.points[1], atol=1e-10)
        np.testing.assert_allclose(cropped.curves[1][0], path.points[1], atol=1e-10)
        a, b = float(path.parameter(first)), path.knots[1]
        t = np.linspace(0., 1., 101)
        c = cropped.curves[0]
        actual = ((1. - t) ** 3)[:, None] * c[0] + (3. * t * (1. - t) ** 2)[:, None] * c[1]
        actual += (3. * t ** 2 * (1. - t))[:, None] * c[2] + (t ** 3)[:, None] * c[3]
        np.testing.assert_allclose(actual, path.jet(a + t * (b - a))[0], atol=1e-10)

    def test_declared_junction_trims_remove_terminal_contact_but_not_visible_crossings(self):
        a = Path(model.Section(1, 1, ((0., 0., 0.), (-10., 30., 0.))))
        b = Path(model.Section(2, 2, ((0., 0., 0.), (10., 30., 0.))))
        with self.assertRaisesRegex(ValueError, 'cross or touch'):
            certify_disjoint((a,), (b,))
        self.assertGreater(certify_disjoint((CenterlineRange.build(a, 5., a.length),),
                                           (CenterlineRange.build(b, 5., b.length),)), 0)
        c = Path(model.Section(3, 3, ((-5., 15., 0.), (5., 15., 0.))))
        d = Path(model.Section(4, 4, ((0., 10., 0.), (0., 20., 0.))))
        with self.assertRaisesRegex(ValueError, 'cross or touch'):
            certify_disjoint((CenterlineRange.build(c, 1., c.length - 1.),),
                             (CenterlineRange.build(d, 1., d.length - 1.),))


if __name__ == '__main__':
    unittest.main()
