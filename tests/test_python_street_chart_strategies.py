"""General paired-rim eligibility, authoritative boundaries and safe fallback."""
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.chart_strategies import FlatPairedRimStrategy
from tested_python_streets.separator_surfaces import FixedRimSeparator
from tested_python_streets.patterns import PatternLibrary


def ribbon(length=80., width=4.):
    corners = np.array(((0., 0.), (length, 0.), (length, width), (0., width)))
    counts = (80, 8, 80, 8)
    xy = np.concatenate([np.linspace(a, b, n + 1)[:-1]
        for a, b, n in zip(corners, np.roll(corners, -1, axis=0), counts)])
    return np.column_stack((xy, np.full(len(xy), .08)))


class PairedRimStrategies(unittest.TestCase):
    def test_flat_ribbon_eligibility_uses_geometry_and_not_orientation_or_ids(self):
        strategy = FlatPairedRimStrategy(PatternLibrary(), spacing=.7, max_cells=50000)
        boundary = ribbon()
        for angle, shift, scale in ((0., (0., 0.), 1.), (.83, (1000., -200.), .5),
                                    (2.5, (-500., 1e6), 2.)):
            c, s = np.cos(angle), np.sin(angle)
            candidate = boundary.copy()
            candidate[:, :2] = scale * candidate[:, :2] @ np.array(((c, s), (-s, c))) + shift
            self.assertTrue(strategy.classify(candidate, datum=0., height=.08).eligible)
            # Perturb a source-side sample without changing the physical plateau.
            candidate[20, :2] += scale * np.array((-.001 * s, .001 * c))
            self.assertTrue(strategy.classify(candidate, datum=0., height=.08).eligible)
        nonflat = boundary.copy()
        nonflat[0, 2] -= 1e-6
        self.assertEqual(strategy.classify(nonflat, datum=0., height=.08).reason, 'nonconstant plateau')
        self.assertFalse(strategy.classify(ribbon(length=10.), datum=0., height=.08).eligible)
        nonflat[:, 2] += 1e12
        nonflat[0, 2] = 1e12 + .08 - .001
        self.assertFalse(strategy.classify(nonflat, datum=1e12, height=.08).eligible)

    def test_selected_ribbon_keeps_exact_rim_flat_heights_and_existing_quality_limits(self):
        boundary = ribbon()
        saved = boundary.copy()
        result = FixedRimSeparator(PatternLibrary()).build(boundary, datum=0., height=.08, exaggeration=1.)
        self.assertEqual(result.metadata['selected_chart'], 'paired-rim')
        np.testing.assert_array_equal(result.mesh.xyz[:len(boundary)], saved)
        np.testing.assert_array_equal(result.mesh.xyz[:, 2], np.full(len(result.mesh.xyz), .08))
        np.testing.assert_array_equal(boundary, saved)
        self.assertLessEqual(result.mesh.validate()['maximum_edge_ratio'], 8.)

    def test_strategy_rejection_retains_the_original_generator(self):
        boundary = ribbon(length=50., width=3.)
        reference = FixedRimSeparator(PatternLibrary(), paired_rim=False).build(
            boundary, datum=0., height=.08, exaggeration=1.)
        with patch.object(FlatPairedRimStrategy, 'build', side_effect=ValueError('test strategy rejection')):
            result = FixedRimSeparator(PatternLibrary()).build(boundary, datum=0., height=.08, exaggeration=1.)
        np.testing.assert_array_equal(result.mesh.xyz, reference.mesh.xyz)
        np.testing.assert_array_equal(result.mesh.quads, reference.mesh.quads)
        self.assertIn('paired-rim: test strategy rejection', result.metadata['failed_charts'])
        self.assertEqual(result.metadata['selected_chart'], reference.metadata['selected_chart'])


if __name__ == '__main__':
    unittest.main()
