"""A computational inversion can move only its interior one-ring vertex."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.quad_seed_untangling import positive_seeds


class QuadSeedUntanglingTests(unittest.TestCase):
    def test_all_incident_quads_and_rim_survive_rotation_and_translation(self):
        for angle, offset in ((0., (0., 0.)), (.71, (1200., -1800.))):
            points = np.array([(x, y) for y in (0., 1., 2.) for x in (0., 1., 2.)])
            points[4] = (-.1, .7)
            rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
            xyz = np.c_[points @ rotation.T+offset, np.zeros(len(points))]
            faces = np.array(((0, 1, 4, 3), (1, 2, 5, 4), (3, 4, 7, 6), (4, 5, 8, 7)))
            mesh = Mesh(xyz, faces, seams={0: np.array((0, 1, 2, 5, 8, 7, 6, 3))})
            uv = np.arange(len(faces)*8,dtype=float).reshape(len(faces),4,2)
            mesh = mesh.with_uv_layer('UVMap',uv)
            with self.assertRaisesRegex(ValueError, 'collapsed, inverted or concave'): mesh.validate()
            candidates = positive_seeds(mesh)
            self.assertTrue(candidates)
            for result, proof in candidates:
                result.validate(); Mesh(result.xyz.astype(np.float32), result.quads).validate()
                self.assertEqual(proof['vertex'], 4)
                self.assertEqual(proof['incident_faces'], [0, 1, 2, 3])
                self.assertLessEqual(proof['displacement'], proof['displacement_bound'])
                np.testing.assert_array_equal(result.xyz[mesh.seams[0]], xyz[mesh.seams[0]])
                np.testing.assert_array_equal(result.quads, faces)
                np.testing.assert_array_equal(result.uv_layers['UVMap'], uv)
                np.testing.assert_array_equal(mesh.xyz, xyz)
                self.assertTrue(proof['exact_physical_rim'] and proof['complete_stored_mesh_positive'])
                self.assertFalse(proof['shape_qualified'])

    def test_an_inverted_fixed_boundary_cannot_be_moved_and_budgets_are_strict(self):
        mesh = Mesh(np.array(((0., 0., 0.), (1., 0., 0.), (.2, -.1, 0.), (0., 1., 0.))), np.array(((0, 1, 2, 3),)))
        with self.assertRaisesRegex(ValueError, 'No bounded one-vertex'): positive_seeds(mesh)
        for args in ({'max_bad_faces': 5}, {'max_candidates': 5}, {'max_cells': 50001}):
            with self.assertRaisesRegex(ValueError, 'bounded.*budgets'): positive_seeds(mesh, **args)


if __name__ == '__main__': unittest.main()
