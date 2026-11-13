"""Output equivalence for faster pattern selection and polygon clipping."""
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from _python_surface_performance_reference import (ReferencePatternLibrary, ReferenceSlopeProjection,
    ReferenceQuadEnergy, ReferenceQuadGradient, visibility_center as original_center)
from tested_python_streets.mesh import Mesh
from tested_python_streets.slope_projection import SlopeProjection
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.quad_charts import visibility_center
from tested_python_streets import quad_charts
from tested_python_streets.quad_workspace import QuadEnergy, QuadGradient


class SurfacePerformanceTests(unittest.TestCase):
    def test_reused_gradient_buffers_keep_bits_and_descent_steps(self):
        random = np.random.default_rng(71940)
        for rows, columns in ((2, 3), (12, 17), (70, 80)):
            pattern = PatternLibrary().grid(rows, columns)
            xy = pattern.uv * (rows, columns)
            xy += random.uniform(-.2, .2, xy.shape)
            mesh = Mesh(np.column_stack((xy, np.zeros(len(xy)))), pattern.quads)
            for power in (1, 2):
                energy = QuadEnergy(pattern.quads, power)
                reference_energy = ReferenceQuadEnergy(pattern.quads, power)
                self.assertTrue(np.isfinite(energy.evaluate(xy)))
                old = ReferenceQuadGradient(pattern.quads, len(xy) + 7)
                new = QuadGradient(pattern.quads, len(xy) + 7)
                for amount in (0., .01, -.02):
                    energy.evaluate(xy + amount * np.sin(xy))
                    self.assertEqual(energy.value,
                        reference_energy.evaluate(xy + amount * np.sin(xy)))
                    np.testing.assert_array_equal(energy.numerator.view(np.uint64),
                        reference_energy.numerator.view(np.uint64))
                    expected = old.evaluate(energy).copy()
                    actual = new.evaluate(energy).copy()
                    np.testing.assert_array_equal(actual.view(np.uint64), expected.view(np.uint64))
                if rows <= 12:
                    solver = quad_charts.FixedRimQuality(iterations=12, target_ratio=1.05, power=power)
                    actual, report = solver.improve(mesh)
                    with patch.object(quad_charts, 'QuadGradient', ReferenceQuadGradient), \
                            patch.object(quad_charts, 'QuadEnergy', ReferenceQuadEnergy):
                        expected, original = solver.improve(mesh)
                    np.testing.assert_array_equal(actual.xyz.view(np.uint64), expected.xyz.view(np.uint64))
                    self.assertEqual(report, original)

    def test_prefetched_slope_coefficients_and_low_memory_fallback_are_exact(self):
        grid = PatternLibrary().grid(5, 4)
        fixed = np.any((grid.uv == 0.) | (grid.uv == 1.), axis=1)
        random = np.random.default_rng(7891)
        for datum in (0., 12.5):
            xyz = np.column_stack((grid.uv * (4., 3.), datum + random.uniform(0., .08, len(grid.uv))))
            xyz[fixed, 2] = datum
            mesh = Mesh(xyz, grid.quads)
            base_bytes = len(mesh.quads) * 2048 + len(mesh.xyz) * 96
            for budget in (base_bytes, 96 * 1024 * 1024):
                args = dict(datum=datum, ceiling=datum + .08, maximum_slope=.035)
                expected, old_report = ReferenceSlopeProjection(memory_bytes=budget).solve(mesh, fixed, **args)
                actual, new_report = SlopeProjection(memory_bytes=budget).solve(mesh, fixed, **args)
                np.testing.assert_array_equal(actual.xyz, expected.xyz)
                np.testing.assert_array_equal(actual.quads, expected.quads)
                self.assertEqual(new_report, old_report)

    def test_priority_queue_preserves_pattern_coordinates_faces_and_boundary_order(self):
        random = np.random.default_rng(9102)
        cases = [(1, 1, 1, 1), (1, 3, 23, 9), (54, 40, 4, 2), (81, 8, 15, 26), (241, 8, 31, 16)]
        cases += [tuple(map(int, 2 * random.integers(0, 20, 4) + 1)) for _ in range(60)]
        for counts in cases:
            with self.subTest(counts=counts):
                expected = ReferencePatternLibrary().quad(*counts)
                actual = PatternLibrary().quad(*counts)
                np.testing.assert_array_equal(actual.uv, expected.uv)
                np.testing.assert_array_equal(actual.quads, expected.quads)
                for side in range(4):
                    np.testing.assert_array_equal(PatternLibrary.boundary(actual, side),
                                                  PatternLibrary.boundary(expected, side))
                self.assertFalse(actual.uv.flags.writeable)
                self.assertFalse(actual.quads.flags.writeable)

    def test_eviction_zero_budget_and_repeated_patterns_preserve_output(self):
        for budget in (0, 1024, 8 * 1024 * 1024):
            cache = PatternLibrary(budget)
            for counts in ((35, 4, 7, 20), (3, 19, 21, 7), (35, 4, 7, 20)):
                actual = cache.quad(*counts)
                expected = ReferencePatternLibrary().quad(*counts)
                np.testing.assert_array_equal(actual.uv, expected.uv)
                np.testing.assert_array_equal(actual.quads, expected.quads)
                self.assertLessEqual(cache.bytes, budget)

    def test_kernel_keeps_original_centroid_and_rejections(self):
        random = np.random.default_rng(712)
        polygons = [np.array(((0.,0.), (4.,0.), (4.,4.), (3.,4.), (3.,1.),
                             (1.,1.), (1.,4.), (0.,4.))),
                    np.array(((0.,0.), (2.,0.), (4.,0.), (4.,2.), (0.,2.))),
                    np.array(((0.,0.), (4.,0.), (4.,1e-11), (0.,1e-11)))]
        for n in (4, 7, 16, 60, 180):
            for _ in range(10):
                angles = np.arange(n) * 2*np.pi/n
                radii = random.uniform(2., 6., n)
                polygons.append(radii[:,None] * np.column_stack((np.cos(angles), np.sin(angles))))
        for polygon in polygons:
            for offset in ((0.,0.), (1e6,-2e6)):
                for epsilon in (1e-10, 1e-7):
                    points = polygon + offset
                    with self.subTest(vertices=len(points), offset=offset, epsilon=epsilon):
                        try:
                            expected = original_center(points, epsilon)
                        except ValueError as error:
                            with self.assertRaises(ValueError) as caught:
                                visibility_center(points, epsilon)
                            self.assertEqual(str(caught.exception), str(error))
                        else:
                            np.testing.assert_array_equal(visibility_center(points, epsilon), expected)


if __name__ == '__main__':
    unittest.main()
