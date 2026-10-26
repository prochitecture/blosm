import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh, cross
from tested_python_streets.model import Key
from tested_python_streets.quad_charts import FixedRimQuality
from tested_python_streets.quad_rotations import QuadEdgeRotations
from tested_python_streets.validation import validate_overlaps


def pair():
    # A six-point patch extracted from an indented physical junction. Its old
    # interior edge creates an 8.45-ratio cell; another all-quad split fits it.
    xy = np.array(((0., 0.), (-1.40264232, .07436300),
                   (-.50666791, .17756867), (-.97282109, .18626701),
                   (.05257885, -.15767405), (-1.22692156, -.09843152)))
    return Mesh(np.column_stack((xy, np.zeros(6))),
                np.array(((1, 0, 2, 3), (0, 1, 5, 4))), 'grass',
                {7: np.array((1, 5, 4, 0, 2, 3))})


def boundary_edges(mesh):
    edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    return unique[counts == 1]


class PythonQuadRotationTests(unittest.TestCase):
    def test_two_quad_rotation_repairs_grading_without_changing_the_region(self):
        original = pair()
        source_xyz, source_faces = original.xyz.copy(), original.quads.copy()
        result, report = QuadEdgeRotations().improve(original, FixedRimQuality(iterations=100))
        self.assertGreater(original.validate()['maximum_edge_ratio'], 8.)
        self.assertLess(result.validate()['maximum_edge_ratio'], 7.)
        self.assertEqual(result.quads.shape, (2, 4))
        np.testing.assert_array_equal(result.xyz, source_xyz)
        np.testing.assert_array_equal(original.quads, source_faces)
        np.testing.assert_array_equal(boundary_edges(result), boundary_edges(original))
        np.testing.assert_array_equal(result.seams[7], original.seams[7])
        self.assertEqual(result.material, 'grass')
        def area(mesh):
            p = mesh.xyz[mesh.quads, :2]
            return .5 * cross(p, np.roll(p, -1, axis=1)).sum()
        self.assertAlmostEqual(area(result), area(original), places=12)
        validate_overlaps({Key('chart', 0): result})
        self.assertEqual(report['history'][0]['rotations'], 1)

    def test_translation_and_repeated_calls_keep_the_same_topology(self):
        original = pair()
        solver = QuadEdgeRotations()
        first, _ = solver.improve(original, FixedRimQuality())
        second, _ = solver.improve(original, FixedRimQuality())
        offset = np.array((1_000_000., -2_000_000., 0.))
        moved, _ = solver.improve(Mesh(original.xyz + offset, original.quads), FixedRimQuality())
        np.testing.assert_array_equal(first.quads, second.quads)
        np.testing.assert_array_equal(first.xyz, second.xyz)
        np.testing.assert_array_equal(first.quads, moved.quads)
        np.testing.assert_allclose(first.xyz, moved.xyz - offset, atol=1e-9, rtol=0.)

    def test_already_qualified_mesh_and_disabled_repair_preserve_buffers(self):
        original = pair()
        disabled, report = QuadEdgeRotations(max_passes=0).improve(original, FixedRimQuality())
        self.assertIs(disabled, original)
        self.assertEqual(report['passes'], 0)
        repaired, _ = QuadEdgeRotations().improve(original, FixedRimQuality())
        unchanged, report = QuadEdgeRotations().improve(repaired, FixedRimQuality())
        self.assertIs(unchanged, repaired)
        self.assertEqual(report['passes'], 0)

    def test_nonmanifold_input_and_invalid_or_exceeded_budgets_are_rejected(self):
        for kwargs in ({'max_passes': -1}, {'max_candidates': 0}, {'max_cells': 0}, {'maximum_ratio': float('nan')}):
            with self.assertRaisesRegex(ValueError, 'budgets'):
                QuadEdgeRotations(**kwargs)
        with self.assertRaisesRegex(ValueError, 'cell budget'):
            QuadEdgeRotations(max_cells=1).improve(pair(), FixedRimQuality())
        original = pair()
        invalid = Mesh(original.xyz, np.vstack((original.quads, original.quads[:1])))
        with self.assertRaisesRegex(ValueError, 'manifold'):
            QuadEdgeRotations().improve(invalid, FixedRimQuality())


if __name__ == '__main__':
    unittest.main()
