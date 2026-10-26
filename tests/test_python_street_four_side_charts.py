import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.four_side_charts import FourSidedQuadChart
from tested_python_streets.mesh import cross
from tested_python_streets.model import Key
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.quad_charts import FixedRimQuality
from tested_python_streets.validation import validate_overlaps


def ribbon():
    first, last = -.45, .45
    def arc(radius, a, b, count):
        angle = np.linspace(a, b, count + 1)[:-1]
        return radius * np.column_stack((np.sin(angle), np.cos(angle)))
    def profile(angle, a, b):
        return np.linspace(a, b, 9)[:-1, None] * np.array((np.sin(angle), np.cos(angle)))
    return np.vstack((arc(20., first, last, 64), profile(last, 20., 26.),
                      arc(26., last, first, 32), profile(first, 26., 20.)))


class PythonFourSidedChartTests(unittest.TestCase):
    def test_streamed_corner_candidates_match_the_scalar_oracle_in_every_chunk(self):
        from _python_four_side_reference import corners
        from tested_python_streets import four_side_charts
        angle = np.linspace(0., 2. * np.pi, 24, endpoint=False)
        radius = 5. + .15 * np.cos(3. * angle)
        curved = radius[:, None] * np.column_stack((np.cos(angle), np.sin(angle)))
        for boundary in (ribbon(), curved, curved + (1234., -5678.)):
            expected = corners(boundary)
            for chunk in (1, 17, 1024):
                with self.subTest(vertices=len(boundary), chunk=chunk):
                    with patch.object(four_side_charts, '_CORNER_CHUNK', chunk):
                        self.assertEqual(FourSidedQuadChart.corners(boundary), expected)

    def test_curved_strip_retains_every_independent_rim_sample_and_is_a_quad_disk(self):
        boundary = ribbon()
        result = FourSidedQuadChart(PatternLibrary()).build(boundary)
        mesh = result.mesh
        np.testing.assert_array_equal(mesh.xyz[mesh.seams[0], :2], boundary)
        self.assertLessEqual(result.metadata['quality']['maximum_edge_ratio'], 8.)
        self.assertGreaterEqual(result.metadata['quality']['minimum_corner_sine'], .03)
        self.assertEqual(mesh.quads.shape[1], 4)
        self.assertLess(len(mesh.quads), 2000)
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        self.assertEqual(len(mesh.xyz) - len(unique) + len(mesh.quads), 1)
        self.assertTrue(set(counts).issubset({1, 2}))
        self.assertEqual(np.sum(counts == 1), len(boundary))
        p = mesh.xyz[mesh.quads, :2]
        area = .5 * cross(p, np.roll(p, -1, axis=1)).sum()
        expected = .5 * cross(boundary, np.roll(boundary, -1, axis=0)).sum()
        self.assertAlmostEqual(area, expected, places=9)
        validate_overlaps({Key('chart', 0): mesh})

    def test_cache_translation_and_repeated_builds_keep_the_same_topology(self):
        boundary, patterns = ribbon(), PatternLibrary()
        chart = FourSidedQuadChart(patterns)
        first, second = chart.build(boundary), chart.build(boundary)
        np.testing.assert_array_equal(first.mesh.xyz, second.mesh.xyz)
        np.testing.assert_array_equal(first.mesh.quads, second.mesh.quads)
        offset = np.array((1_000_000., -2_000_000.))
        moved = chart.build(boundary + offset)
        np.testing.assert_array_equal(first.mesh.quads, moved.mesh.quads)
        np.testing.assert_allclose(first.mesh.xyz[:, :2], moved.mesh.xyz[:, :2] - offset, atol=1e-7, rtol=0.)
        self.assertGreater(patterns.hits, 0)

    def test_one_quad_needs_no_internal_points_or_relaxation(self):
        boundary = np.array(((0., 0.), (4., 0.), (4., 2.), (0., 2.)))
        result = FourSidedQuadChart(PatternLibrary()).build(boundary)
        self.assertEqual(result.mesh.xyz.shape, (4, 3))
        self.assertEqual(result.mesh.quads.shape, (1, 4))
        np.testing.assert_array_equal(result.mesh.xyz[:, :2], boundary)

    def test_squared_energy_repairs_a_nearly_flat_internal_corner_with_fixed_rim(self):
        from tested_python_streets.mesh import Mesh
        pattern = PatternLibrary().grid(2, 2)
        points = np.column_stack((2. * pattern.uv, np.zeros(len(pattern.uv))))
        points[4, :2] = .51, .51
        original = Mesh(points.copy(), pattern.quads)
        result, report = FixedRimQuality(power=2, target_ratio=1.6, iterations=100).improve(original)
        boundary = np.r_[0:4, 5:9]
        np.testing.assert_array_equal(result.xyz[boundary], points[boundary])
        np.testing.assert_array_equal(original.xyz, points)
        self.assertGreater(report['steps'], 0)
        self.assertGreater(result.validate()['minimum_corner'], original.validate()['minimum_corner'])
        self.assertLess(result.validate()['maximum_edge_ratio'], original.validate()['maximum_edge_ratio'])

    def test_invalid_boundaries_and_budgets_are_explicit(self):
        for settings in ({'max_cells': 0}, {'max_candidates': 0}, {'quality_iterations': -1}, {'maximum_ratio': float('nan')}):
            with self.assertRaisesRegex(ValueError, 'budgets'):
                FourSidedQuadChart(PatternLibrary(), **settings)
        with self.assertRaisesRegex(ValueError, 'cell budget'):
            FourSidedQuadChart(PatternLibrary(), max_cells=1).build(ribbon())
        with self.assertRaisesRegex(ValueError, 'even'):
            FourSidedQuadChart(PatternLibrary()).build(ribbon()[:-1])
        with self.assertRaises(ValueError):
            FourSidedQuadChart(PatternLibrary()).build(np.array(((0., 0.), (2., 2.), (0., 2.), (2., 0.))))
        with self.assertRaises(ValueError):
            FixedRimQuality(power=3)


if __name__ == '__main__':
    unittest.main()
