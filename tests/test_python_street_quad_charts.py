import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import fixture_networks
from test_python_street_compounds import footprints
from tested_python_streets import model
from tested_python_streets.compounds import sample_footprint
from tested_python_streets.mesh import cross
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.quad_charts import LayeredQuadChart, RadialTransition, QuadChartQualityError
from tested_python_streets.validation import validate_overlaps


class PythonDirectQuadChartTests(unittest.TestCase):
    def test_three_to_one_template_is_positive_for_varied_concave_radial_fronts(self):
        patterns = PatternLibrary()
        transition = RadialTransition(patterns)
        random = np.random.default_rng(881)
        for _ in range(120):
            angles = np.r_[0., np.cumsum(random.uniform(.03, .4, 3))]
            directions = np.column_stack((np.cos(angles), np.sin(angles)))
            outer = random.uniform(3., 8., 4)[:, None] * directions
            inner = random.uniform(.3, 1., 2)[:, None] * directions[[0, 3]]
            candidates = transition.candidates(outer, inner)
            self.assertTrue(candidates)
            for ratio, points in candidates:
                np.testing.assert_array_equal(points[transition.fine], outer)
                np.testing.assert_array_equal(points[transition.coarse], inner)
                p = points[transition.quads]
                e = np.roll(p, -1, axis=1) - p
                self.assertGreater(float(cross(e, np.roll(e, -1, axis=1)).min()), 1e-10)
                self.assertEqual(transition.quads.shape, (4, 4))
        self.assertIs(transition.pattern, patterns.strip(1, 3, 1))

    def test_six_families_preserve_rim_and_pass_positive_quad_disk_checks(self):
        for network in fixture_networks(model)[:6]:
            for owner, layout in footprints(network).items():
                with self.subTest(network=network.uuid, site=owner):
                    polygon = sample_footprint(owner, layout).region.xy
                    result = LayeredQuadChart(PatternLibrary()).build(polygon)
                    mesh = result.mesh
                    np.testing.assert_array_equal(mesh.xyz[mesh.seams[0], :2], polygon)
                    self.assertLessEqual(result.metadata['quality']['maximum_edge_ratio'], 8.)
                    self.assertGreater(result.metadata['quality']['minimum_corner_sine'], .03)
                    edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
                    unique, counts = np.unique(edges, axis=0, return_counts=True)
                    self.assertEqual(len(mesh.xyz) - len(unique) + len(mesh.quads), 1)
                    self.assertEqual(int(np.sum(counts == 1)), len(polygon))
                    self.assertTrue(set(counts).issubset({1, 2}))
                    validate_overlaps({model.Key('chart', 1): mesh})

    def test_cache_reuse_and_translation_preserve_topology_and_dimensions(self):
        network = fixture_networks(model)[0]
        owner, layout = next(iter(footprints(network).items()))
        polygon = sample_footprint(owner, layout).region.xy
        patterns = PatternLibrary()
        chart = LayeredQuadChart(patterns)
        first, second = chart.build(polygon), chart.build(polygon)
        np.testing.assert_array_equal(first.mesh.xyz, second.mesh.xyz)
        np.testing.assert_array_equal(first.mesh.quads, second.mesh.quads)
        offset = np.array((1_000_000., -2_000_000.))
        translated = chart.build(polygon + offset)
        np.testing.assert_array_equal(first.mesh.quads, translated.mesh.quads)
        np.testing.assert_allclose(first.mesh.xyz[:, :2], translated.mesh.xyz[:, :2] - offset, atol=1e-7, rtol=0.)
        self.assertGreater(patterns.hits, 0)

    def test_quality_failure_retains_the_positive_chart_for_an_assembled_retry(self):
        network = fixture_networks(model)[0]
        owner, layout = next(iter(footprints(network).items()))
        polygon = sample_footprint(owner, layout).region.xy
        chart = LayeredQuadChart(PatternLibrary(), maximum_ratio=1.01, quality_iterations=0)
        with self.assertRaises(QuadChartQualityError) as caught:
            chart.build(polygon)
        result = caught.exception.result
        self.assertGreater(result.mesh.validate()['minimum_corner'], 1e-10)
        self.assertGreater(result.metadata['quality']['maximum_edge_ratio'], 1.01)
        np.testing.assert_array_equal(result.mesh.xyz[result.mesh.seams[0], :2], polygon)
        validate_overlaps({model.Key('positive_retry', 1): result.mesh})

    def test_unsupported_kernel_parity_and_resource_requests_are_diagnosed(self):
        chart = LayeredQuadChart(PatternLibrary())
        u_shape = np.array(((0., 0.), (4., 0.), (4., 4.), (3., 4.), (3., 1.), (1., 1.), (1., 4.), (0., 4.)))
        with self.assertRaisesRegex(ValueError, 'visibility'):
            chart.build(u_shape)
        angles = np.arange(5) * 2. * np.pi / 5
        with self.assertRaisesRegex(ValueError, 'even'):
            chart.build(np.column_stack((np.cos(angles), np.sin(angles))))
        network = fixture_networks(model)[0]
        owner, layout = next(iter(footprints(network).items()))
        rim = sample_footprint(owner, layout).region.xy
        with self.assertRaisesRegex(ValueError, 'cell budget'):
            LayeredQuadChart(PatternLibrary(), max_cells=4).build(rim)
        with self.assertRaisesRegex(ValueError, 'layer budget'):
            LayeredQuadChart(PatternLibrary(), max_layers=1).build(rim)


if __name__ == '__main__':
    unittest.main()
