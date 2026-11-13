"""Analytical boundary rules retain jets, envelopes and the fitting fallback."""
from dataclasses import replace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from _python_surface_performance_reference import ReferenceBoundaryJetField
from test_python_street_surfaces import boundary_data, LinearRule
from tested_python_streets.conditions import CurbEdge, PlateauEdge
from tested_python_streets.surface_fields import BoundaryJetField


class BadCoefficients(LinearRule):
    def polynomial(self):
        return np.full((6, 3), 10.)


class BoundaryPolynomialTests(unittest.TestCase):
    def test_analytical_coefficients_match_sampled_heights_and_gradients(self):
        random = np.random.default_rng(8732)
        t = np.r_[0., 1., np.linspace(0., 1., 401), random.uniform(0., 1., 200)]
        for datum in (0., 12.5, -302.75):
            rules = [PlateauEdge(datum + .08)]
            for _ in range(30):
                length = float(random.uniform(.05, 100.))
                first, last = sorted(random.uniform(0., length, 2))
                angle = float(random.uniform(-np.pi, np.pi))
                rules.append(CurbEdge(datum, datum + float(random.uniform(-.08, .08)),
                    float(random.uniform(0., .1)), float(random.uniform(0., .1)),
                    length, first, last, (np.cos(angle), np.sin(angle))))
            for rule in rules:
                height, gradient = rule.sample(t)
                expected = np.column_stack((height, gradient))
                actual = np.polynomial.polynomial.polyval(t, rule.polynomial()).T
                np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=0.)

    def test_field_preserves_original_fitted_surface_and_slope_checks(self):
        boundary = boundary_data()
        random = np.random.default_rng(382)
        points = np.vstack((random.uniform(0., 1., (1001, 2)), boundary.xy))
        for datum in (0., 12.5):
            for constant in (True, False):
                rules = (PlateauEdge(datum + .04),) * 4 if constant else tuple(
                    CurbEdge(datum + .01, datum + .07, .025, .04, 4., i, i + 1,
                             tuple(tangent))
                    for i, tangent in enumerate(((1., 0.), (0., 1.), (-1., 0.), (0., -1.))))
                value = replace(boundary, rules=rules, datum=datum, ceiling=datum + .1)
                old = ReferenceBoundaryJetField((value,))
                new = BoundaryJetField((value,))
                np.testing.assert_allclose(new(points), old(points), atol=2e-12, rtol=0.)
                self.assertAlmostEqual(new.maximum_boundary_slope, old.maximum_boundary_slope, places=14)
                self.assertLessEqual(new.metadata['maximum_sampled_fit_error'], 1e-8)
                with patch.object(np.polynomial.polynomial, 'polyfit', side_effect=AssertionError('unexpected fit')):
                    BoundaryJetField((value,))

    def test_unqualified_coefficients_use_the_original_residual_checked_fit(self):
        value = boundary_data()
        value = replace(value, rules=tuple(BadCoefficients(r.first, r.last) for r in value.rules))
        old, new = ReferenceBoundaryJetField((value,)), BoundaryJetField((value,))
        np.testing.assert_array_equal(new.coefficients, old.coefficients)
        self.assertEqual(new.maximum_boundary_slope, old.maximum_boundary_slope)
        np.testing.assert_array_equal(new(np.array(((.25, .5), (.8, .3)))), old(np.array(((.25, .5), (.8, .3)))))


if __name__ == '__main__':
    unittest.main()
