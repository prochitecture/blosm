"""Singular fixed-rim parity cells need a different direct-quad topology."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.positive_quad_patches import positive_patches


def fixture(angle=0., shift=(0., 0.)):
    # The boundary has one collinear and one reflex vertex. Neither existing
    # two-quad diagonal gives a positive chart; three quads can connect both
    # constrained vertices to one certified interior point.
    xy = np.array(((0., 0.), (1., 0.), (2., 0.), (1., .3), (1., 1.), (0., 1.)))
    rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
    return Mesh(np.c_[xy @ rotation.T+shift, np.zeros(len(xy))],
                np.array(((0, 1, 2, 3), (0, 3, 4, 5)), dtype=np.int32), seams={0: np.arange(6)})


def boundary_edges(mesh):
    edge = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    unique, counts = np.unique(edge, axis=0, return_counts=True)
    return unique[counts == 1], len(unique), counts


class PositiveQuadPatchTests(unittest.TestCase):
    def test_direct_three_quad_patch_keeps_the_whole_fixed_rim(self):
        for angle, shift in ((0., (0., 0.)), (.713, (1200., -1800.))):
            with self.subTest(angle=angle):
                original = fixture(angle, shift)
                before = original.xyz.copy()
                candidates = positive_patches(original)
                self.assertTrue(candidates)
                for mesh, report in candidates:
                    self.assertEqual(len(mesh.quads), 3)
                    self.assertEqual(len(mesh.xyz), 7)
                    self.assertEqual(report['added_interior_vertices'], 1)
                    self.assertFalse(report['triangles_created'])
                    self.assertFalse(report['shape_qualified'])
                    mesh.validate()
                    mesh.with_positions(mesh.xyz.astype(np.float32)).validate()
                    np.testing.assert_array_equal(mesh.xyz[:6], before)
                    np.testing.assert_array_equal(mesh.seams[0], np.arange(6))
                    exterior, edge_count, counts = boundary_edges(mesh)
                    np.testing.assert_array_equal(exterior, boundary_edges(original)[0])
                    self.assertLessEqual(counts.max(), 2)
                    self.assertEqual(len(mesh.xyz)-edge_count+len(mesh.quads), 1)
                np.testing.assert_array_equal(original.xyz, before)

    def test_budget_and_authored_uvs_cannot_be_silently_changed(self):
        mesh = fixture()
        with self.assertRaisesRegex(ValueError, 'No bounded positive'):
            positive_patches(mesh, max_cells=2)
        charted = mesh.with_uv_layer('Authored', np.zeros((2, 4, 2)))
        with self.assertRaisesRegex(ValueError, 'uncharted planar'):
            positive_patches(charted)
        with self.assertRaisesRegex(ValueError, 'bounded uncharted'):
            positive_patches(mesh, max_candidates=5)
        for invalid in ((True,),(-1,),(2,),(0,0),(0,1,0,1,0),[0]):
            with self.assertRaisesRegex(ValueError,'distinct existing faces'):
                positive_patches(mesh,constrained_faces=invalid)

    def test_positive_double_quad_can_need_an_opt_in_stored_topology_repair(self):
        mesh = fixture(shift=(1200., 1800.))
        mesh.xyz[1, 1] -= 1e-7
        mesh.validate()
        with self.assertRaisesRegex(ValueError, 'collapsed, inverted or concave'):
            mesh.with_positions(mesh.xyz.astype(np.float32)).validate()
        with self.assertRaisesRegex(ValueError, 'one bounded singular patch'):
            positive_patches(mesh)
        candidates = positive_patches(mesh, include_stored_failures=True)
        self.assertTrue(candidates)
        for candidate, report in candidates:
            self.assertTrue(report['admitted_stored_failures'])
            np.testing.assert_array_equal(candidate.xyz[:len(mesh.xyz)], mesh.xyz)
            np.testing.assert_array_equal(boundary_edges(candidate)[0], boundary_edges(mesh)[0])
            candidate.validate()
            candidate.with_positions(candidate.xyz.astype(np.float32)).validate()


if __name__ == '__main__':
    unittest.main()
