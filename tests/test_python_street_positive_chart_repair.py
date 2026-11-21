"""Bounded topology repair keeps authoritative vertex IDs and profile aliases."""
import unittest
import numpy as np
import _python_street_package
from test_python_street_finite_scope_charts import boundary
from tested_python_streets.four_side_charts import FourSidedQuadChart, quality
from tested_python_streets.mesh import Mesh, publication_coordinates
from tested_python_streets.partition import BoundaryRef, CompoundQuadChart
from tested_python_streets.patterns import PatternLibrary


class PositiveChartRepairTests(unittest.TestCase):
    def test_progressive_repair_can_stop_before_changing_the_input(self):
        patterns = PatternLibrary()
        mesh, _ = FourSidedQuadChart(patterns).map(boundary(), (0, 38, 50, 76))
        original = mesh.xyz.copy()
        faces = mesh.quads.copy()
        steps = CompoundQuadChart(patterns)._repair_positive_candidate_steps(mesh)
        next(steps)
        steps.close()
        np.testing.assert_array_equal(mesh.xyz, original)
        np.testing.assert_array_equal(mesh.quads, faces)

    def test_existing_ids_preserve_seams_boundary_references_and_profile_aliases(self):
        for offset in ((0., 0., 0.), (1200., -1800., 0.)):
            patterns = PatternLibrary()
            mesh, _ = FourSidedQuadChart(patterns).map(boundary(), (0, 38, 50, 76))
            mesh.xyz += offset
            original = mesh.xyz.copy()
            rim = mesh.seams[0].copy()
            references = {int(node): BoundaryRef(0, i, 0.) for i, node in enumerate(rim)}
            aliases = {(99, 0, i/len(rim)): int(node) for i, node in enumerate(rim)}
            result, report = CompoundQuadChart(patterns)._repair_positive_candidate(mesh)
            self.assertTrue(report['original_vertex_ids_preserved'])
            self.assertGreater(len(result.quads), len(mesh.quads))
            np.testing.assert_array_equal(result.seams[0], rim)
            np.testing.assert_array_equal(result.xyz[rim], original[rim])
            self.assertTrue(set(range(len(original))).issubset(set(result.quads.ravel())))
            for node in (*references, *aliases.values()):
                np.testing.assert_array_equal(result.xyz[node], original[node])
            stored = quality(result.with_positions(publication_coordinates(result.xyz, 1., physical=True)))
            self.assertLessEqual(stored['maximum_edge_ratio'], 8.)
            self.assertGreaterEqual(stored['minimum_corner_sine'], .03)
            np.testing.assert_array_equal(mesh.xyz, original)

    def test_a_collapsed_stored_boundary_remains_unqualified(self):
        xyz = np.array(((1e8, 0., 0.), (1e8+1e-4, 0., 0.),
                        (1e8+1e-4, 1e-4, 0.), (1e8, 1e-4, 0.)))
        mesh = Mesh(xyz.copy(), np.array(((0, 1, 2, 3),)), seams={0: np.arange(4)})
        mesh.validate()
        with self.assertRaises(ValueError):
            CompoundQuadChart(PatternLibrary())._repair_positive_candidate(mesh)
        np.testing.assert_array_equal(mesh.xyz, xyz)


if __name__ == '__main__':
    unittest.main()
