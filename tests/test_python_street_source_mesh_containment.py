"""A body cannot bridge a source gap or hide a whole separator in one quad."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.source_mesh_containment import certify_mesh


def strips():
    x = np.array((0., 1., 5., 6.))
    inner = np.array((-1., 1., 1., -1.))
    return [np.c_[x, np.full(4, 4.)], np.c_[x, inner],
            np.c_[x, -inner], np.c_[x, np.full(4, -4.)]]


def quad(xy):
    return Mesh(np.c_[xy, np.zeros(4)], np.array(((0, 1, 2, 3),)))


class SourceMeshContainmentTests(unittest.TestCase):
    def test_valid_source_quad_and_original_error_allocation_survive_rigid_motion(self):
        curves = strips()
        points = np.array(((2., 2.), (3., 2.), (3., 3.), (2., 3.)))
        angle = .713
        rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        for transform, offset in ((np.eye(2), np.zeros(2)), (rotation, np.array((1800., -1200.)))):
            originals = [curve.copy() for curve in curves]
            result = certify_mesh(quad(points@transform.T+offset),
                [curve@transform.T+offset for curve in curves],
                source_error=.0001, mesh_error=.0002, maximum_error=.0005)
            self.assertEqual(result['holes'], 1)
            self.assertEqual(result['complete_quads'], 1)
            self.assertAlmostEqual(result['maximum_combined_error'], .0003)
            for a, b in zip(curves, originals):np.testing.assert_array_equal(a, b)

    def test_whole_hole_inside_quad_is_rejected_even_when_no_quad_edge_crosses_it(self):
        mesh = quad(np.array(((0., -4.), (6., -4.), (6., 4.), (0., 4.))))
        with self.assertRaisesRegex(ValueError, 'encloses an original source boundary'):
            certify_mesh(mesh, strips(), source_error=.0001, mesh_error=.0002, maximum_error=.0005)

    def test_original_finite_caps_and_real_source_gaps_cannot_be_extended_or_filled(self):
        for points in (np.array(((-.1, 2.), (.5, 2.), (.5, 3.), (-.1, 3.))),
                       np.array(((2., -.5), (3., -.5), (3., .5), (2., .5)))):
            with self.subTest(points=points.tolist()), self.assertRaisesRegex(ValueError, 'finite original|real original'):
                certify_mesh(quad(points), strips(), source_error=.0001, mesh_error=.0002, maximum_error=.0005)
        with self.assertRaisesRegex(ValueError, 'containment budget'):
            certify_mesh(quad(np.array(((2., 2.), (3., 2.), (3., 3.), (2., 3.)))), strips(),
                source_error=.0004, mesh_error=.0004, maximum_error=.0005)


if __name__ == '__main__':unittest.main()
