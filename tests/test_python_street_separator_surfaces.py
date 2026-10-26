"""Alternative median charts preserve the complete boundary and height envelope."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.separator_surfaces import FixedRimSeparator, PartitionedSeparatorChart
from tested_python_streets.mesh import cross


def rectangle():
    points = []
    for a, b in (((0., 0.), (10., 0.)), ((10., 0.), (10., 6.)),
                 ((10., 6.), (0., 6.)), ((0., 6.), (0., 0.))):
        points.extend(np.linspace(a, b, 9)[:-1])
    xy = np.asarray(points)
    z = .8 * np.clip(xy[:, 0] / 3., 0., 1.)
    return np.column_stack((xy, z))


class PythonFixedRimSeparators(unittest.TestCase):
    def test_rotated_nonvisible_rim_retains_exact_source_endpoints_through_partitioning(self):
        corners = np.array(((0., 0.), (4., 0.), (4., 4.), (3., 4.),
                            (3., 1.), (1., 1.), (1., 4.), (0., 4.)))
        boundary = np.concatenate([np.linspace(a, b, 9)[:-1]
            for a, b in zip(corners, np.roll(corners, -1, axis=0))])
        angle = .43
        rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        boundary = boundary @ rotation.T + (.013540461395702724, 1.7)
        saved = boundary.copy()
        result = PartitionedSeparatorChart(PatternLibrary()).build(boundary)
        mesh = result.mesh
        np.testing.assert_array_equal(mesh.xyz[:len(boundary), :2], saved)
        np.testing.assert_array_equal(boundary, saved)
        self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
        restoration = result.metadata['source_boundary_restoration']
        self.assertTrue(restoration['source_identity_bijection'])
        self.assertTrue(restoration['no_new_physical_samples'])
        self.assertGreater(restoration['changed_vertices'], 0)
        faces = mesh.xyz[mesh.quads, :2]
        self.assertAlmostEqual(.5 * np.sum(cross(faces, np.roll(faces, -1, axis=1))), 10., places=8)
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        expected = {tuple(sorted((i, (i + 1) % len(boundary)))) for i in range(len(boundary))}
        self.assertEqual({tuple(edge) for edge in unique[counts == 1]}, expected)
        self.assertTrue(np.all(counts <= 2))

    def test_nonvisible_u_shaped_median_uses_shared_quad_blocks_without_changing_its_boundary(self):
        corners = np.array(((0., 0.), (4., 0.), (4., 4.), (3., 4.), (3., 1.), (1., 1.), (1., 4.), (0., 4.)))
        xy = np.concatenate([np.linspace(a, b, 9)[:-1] for a, b in zip(corners, np.roll(corners, -1, axis=0))])
        boundary = np.column_stack((xy, np.full(len(xy), .8)))
        result = FixedRimSeparator(PatternLibrary()).build(boundary, datum=0., height=.8, exaggeration=10.)
        self.assertEqual(result.metadata['selected_chart'], 'partitioned')
        mesh = result.mesh
        np.testing.assert_array_equal(mesh.xyz[:len(boundary)], boundary)
        np.testing.assert_array_equal(mesh.xyz[:, 2], np.full(len(mesh.xyz), .8))
        self.assertGreater(mesh.validate()['minimum_corner'], 0.)
        faces = mesh.xyz[mesh.quads, :2]
        self.assertAlmostEqual(.5 * np.sum(cross(faces, np.roll(faces, -1, axis=1))), 10., places=8)
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, count = np.unique(edges, axis=0, return_counts=True)
        expected = {tuple(sorted((i, (i + 1) % len(boundary)))) for i in range(len(boundary))}
        self.assertEqual({tuple(edge) for edge in unique[count == 1]}, expected)
        self.assertTrue(np.all(count <= 2))

    def test_chart_keeps_all_boundary_edges_and_bounds_the_plateau(self):
        boundary = rectangle()
        saved = boundary.copy()
        result = FixedRimSeparator(PatternLibrary()).build(boundary, datum=0., height=.8, exaggeration=10.)
        mesh = result.mesh
        np.testing.assert_array_equal(boundary, saved)
        np.testing.assert_array_equal(mesh.xyz[:len(boundary)], boundary)
        self.assertGreater(mesh.validate()['minimum_corner'], 0.)
        self.assertGreaterEqual(mesh.xyz[:, 2].min(), 0.)
        self.assertLessEqual(mesh.xyz[:, 2].max(), .8)
        xy = mesh.xyz[mesh.quads, :2]
        self.assertAlmostEqual(.5 * np.sum(cross(xy, np.roll(xy, -1, axis=1))), 60., places=8)
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, count = np.unique(edges, axis=0, return_counts=True)
        expected = {tuple(sorted((i, (i + 1) % len(boundary)))) for i in range(len(boundary))}
        self.assertEqual({tuple(edge) for edge in unique[count == 1]}, expected)
        self.assertTrue(np.all(count <= 2))

    def test_common_datum_and_zero_exaggeration_are_supported(self):
        boundary = rectangle()
        raised = boundary.copy()
        raised[:, 2] += 12.5
        first = FixedRimSeparator(PatternLibrary()).build(boundary, datum=0., height=.8, exaggeration=10.).mesh
        second = FixedRimSeparator(PatternLibrary()).build(raised, datum=12.5, height=.8, exaggeration=10.).mesh
        np.testing.assert_allclose(second.xyz[:, 2] - 12.5, first.xyz[:, 2], atol=1e-12)
        boundary[:, 2] = 12.5
        flat = FixedRimSeparator(PatternLibrary()).build(boundary, datum=12.5, height=0., exaggeration=0.)
        np.testing.assert_array_equal(flat.mesh.xyz[:, 2], np.full(len(flat.mesh.xyz), 12.5))

    def test_incompatible_rims_and_exhausted_budgets_are_diagnostic(self):
        for boundary in (rectangle()[:-1], rectangle()[::-1]):
            with self.assertRaisesRegex(ValueError, 'even counterclockwise'):
                FixedRimSeparator(PatternLibrary()).build(boundary, datum=0., height=.8, exaggeration=10.)
        with self.assertRaisesRegex(ValueError, 'No accepted.*cell budget'):
            FixedRimSeparator(PatternLibrary(), max_cells=4).build(rectangle(), datum=0., height=.8, exaggeration=10.)


if __name__ == '__main__':
    unittest.main()
