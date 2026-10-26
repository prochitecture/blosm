"""Independent derivatives, buffer isolation and fixed-rim descent invariants."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.quad_charts import FixedRimQuality
from tested_python_streets.quad_workspace import QuadEnergy, QuadGradient


def distorted_mesh(hole=False):
    pattern = PatternLibrary().grid(4, 4)
    xy = pattern.uv * 4.
    xy[12] += (.8, .1)
    faces = pattern.quads
    if hole:
        faces = np.delete(faces, 5, axis=0)
    z = .03 * xy[:, 0] + .07 * xy[:, 1]
    return Mesh(np.column_stack((xy, z)), faces, 'test_paving', {7: np.array((0, 1, 2))})


def reference_energy(xy, faces, power):
    # Independent scalar corner formula; no production workspace/gradient calls.
    total = 0.
    for face in faces:
        for i in range(4):
            a = xy[face[(i + 1) % 4]] - xy[face[i]]
            b = xy[face[(i + 2) % 4]] - xy[face[(i + 1) % 4]]
            area = a[0] * b[1] - a[1] * b[0]
            shape = (np.dot(a, a) + np.dot(b, b)) / area
            total += shape if power == 1 else shape * shape / 4.
    return total


class PythonQuadWorkspaceTests(unittest.TestCase):
    def test_gradient_accumulation_keeps_the_exact_legacy_stream_and_unused_vertices(self):
        mesh = distorted_mesh(hole=True)
        random = np.random.default_rng(4902)
        order = random.permutation(len(mesh.quads))
        vertices = len(mesh.xyz) + 7
        for power in (1, 2):
            energy = QuadEnergy(mesh.quads[order], power)
            energy.evaluate(mesh.xyz[:, :2])
            scatter = QuadGradient(mesh.quads[order], vertices)
            actual = scatter.evaluate(energy).copy()
            expected = np.zeros_like(actual)
            np.add.at(expected, scatter.targets, scatter.contributions.reshape(-1, 2))
            np.testing.assert_array_equal(actual.view(np.uint64), expected.view(np.uint64))
            np.testing.assert_array_equal(actual[-7:], 0.)
            # Reusing the buffer cannot retain values on unreferenced vertices.
            scatter.gradient[-7:] = 99.
            np.testing.assert_array_equal(scatter.evaluate(energy).view(np.uint64), expected.view(np.uint64))

    def test_energy_and_directional_derivatives_match_independent_corner_formula(self):
        mesh = distorted_mesh()
        xy = mesh.xyz[:, :2].copy()
        random = np.random.default_rng(3038)
        for power in (1, 2):
            energy = QuadEnergy(mesh.quads, power)
            self.assertAlmostEqual(energy.evaluate(xy), reference_energy(xy, mesh.quads, power), places=10)
            gradient = QuadGradient(mesh.quads, len(xy)).evaluate(energy).copy()
            for _ in range(5):
                direction = random.normal(size=xy.shape)
                direction /= np.linalg.norm(direction)
                epsilon = 1e-6
                numeric = (reference_energy(xy + epsilon * direction, mesh.quads, power)
                           - reference_energy(xy - epsilon * direction, mesh.quads, power)) / (2. * epsilon)
                self.assertAlmostEqual(float(np.sum(gradient * direction)), numeric,
                                       delta=max(1., abs(numeric)) * 1e-7)

    def test_rejected_trial_and_repeated_gradient_leave_accepted_state_unchanged(self):
        mesh = distorted_mesh()
        xy = mesh.xyz[:, :2].copy()
        for power in (1, 2):
            current, trial = QuadEnergy(mesh.quads, power), QuadEnergy(mesh.quads, power)
            value = current.evaluate(xy)
            gradient = QuadGradient(mesh.quads, len(xy))
            expected = gradient.evaluate(current).copy()
            edge, jac, numerator = current.edge.copy(), current.jac.copy(), current.numerator.copy()
            self.assertEqual(trial.evaluate(xy * (-1., 1.)), float('inf'))
            np.testing.assert_array_equal(current.edge, edge)
            np.testing.assert_array_equal(current.jac, jac)
            np.testing.assert_array_equal(current.numerator, numerator)
            np.testing.assert_array_equal(gradient.evaluate(current), expected)
            self.assertEqual(current.evaluate(xy), value)
            np.testing.assert_array_equal(gradient.evaluate(current), expected)
            self.assertEqual(trial.evaluate(xy), value)

    def test_descent_preserves_exterior_and_hole_rims_heights_topology_and_metadata(self):
        for hole in (False, True):
            mesh = distorted_mesh(hole)
            before = mesh.xyz.copy()
            edges = np.sort(np.concatenate((mesh.quads[:, :2], mesh.quads[:, 1:3],
                mesh.quads[:, 2:4], mesh.quads[:, (3, 0)])), axis=1)
            edges, counts = np.unique(edges, axis=0, return_counts=True)
            fixed = np.unique(edges[counts == 1])
            for power in (1, 2):
                result, report = FixedRimQuality(iterations=8, target_ratio=1.1, power=power).improve(mesh)
                self.assertGreater(report['steps'], 0)
                self.assertLessEqual(report['steps'], 8)
                self.assertLess(report['energy'], reference_energy(before[:, :2], mesh.quads, power))
                self.assertEqual(report['boundary_error'], 0.)
                self.assertGreater(result.validate()['minimum_corner'], 1e-10)
                np.testing.assert_array_equal(result.xyz[fixed], before[fixed])
                np.testing.assert_array_equal(result.xyz[:, 2], before[:, 2])
                np.testing.assert_array_equal(result.quads, mesh.quads)
                np.testing.assert_array_equal(mesh.xyz, before)
                np.testing.assert_array_equal(result.seams[7], mesh.seams[7])
                self.assertEqual(result.material, mesh.material)

    def test_iteration_zero_and_already_regular_mesh_remain_exact(self):
        mesh = distorted_mesh()
        result, report = FixedRimQuality(iterations=0).improve(mesh)
        np.testing.assert_array_equal(result.xyz, mesh.xyz)
        self.assertEqual(report['steps'], 0)
        pattern = PatternLibrary().grid(3, 3)
        regular = Mesh(np.column_stack((pattern.uv, np.zeros(len(pattern.uv)))), pattern.quads)
        result, report = FixedRimQuality().improve(regular)
        np.testing.assert_array_equal(result.xyz, regular.xyz)
        self.assertEqual(report['steps'], 0)

    def test_nonpositive_initial_mesh_is_rejected_before_descent(self):
        mesh = distorted_mesh()
        mesh.xyz[:, 0] *= -1.
        for iterations in (0, 10):
            with self.assertRaisesRegex(ValueError, 'positive initial mesh'):
                FixedRimQuality(iterations=iterations).improve(mesh)

    def test_rejections_recheck_current_coordinates_and_connectivity(self):
        xy = np.array(((0., 0.), (1., 0.), (1., 1.), (0., 1.),
                       (3., 0.), (4., 0.), (4., 1.), (3., 1.)))
        for power in (1, 2):
            faces = np.array(((0, 1, 2, 3), (4, 5, 6, 7)), dtype=np.int32)
            energy = QuadEnergy(faces, power)
            self.assertAlmostEqual(energy.evaluate(xy), reference_energy(xy, faces, power))
            bad = xy.copy()
            bad[2] = (1., -.5)
            for _ in range(3):
                self.assertEqual(energy.evaluate(bad), float('inf'))
            # The invalid vertices cease to belong to the face. A rejection
            # remembered from another topology must not reject this valid input.
            faces[0] = faces[1]
            self.assertAlmostEqual(energy.evaluate(bad), reference_energy(bad, faces, power))
            np.testing.assert_array_equal(energy.edge, np.roll(bad[faces], -1, axis=1) - bad[faces])
            # Restore both topology and coordinates, then move the bad corner
            # to the other face; the previous witness is no longer sufficient.
            faces[0] = (0, 1, 2, 3)
            self.assertAlmostEqual(energy.evaluate(xy), reference_energy(xy, faces, power))
            other = xy.copy()
            other[6] = (4., -.5)
            self.assertEqual(energy.evaluate(other), float('inf'))
            self.assertAlmostEqual(energy.evaluate(xy), reference_energy(xy, faces, power))
            self.assertTrue(np.all(np.isfinite(QuadGradient(faces, len(xy)).evaluate(energy))))

    def test_rejection_threshold_is_unchanged_for_neighboring_float_values(self):
        faces = np.array(((0, 1, 2, 3),), dtype=np.int32)
        for power in (1, 2):
            energy = QuadEnergy(faces, power)
            for height in (-1., 0., 1e-11, np.nextafter(1e-10, 0.), 1e-10,
                           np.nextafter(1e-10, 1.), 1e-9, 1., -1., 2.):
                xy = np.array(((0., 0.), (1., 0.), (1., height), (0., height)))
                value = energy.evaluate(xy)
                if height <= 1e-10:
                    self.assertEqual(value, float('inf'))
                else:
                    expected = reference_energy(xy, faces, power)
                    self.assertAlmostEqual(value, expected, delta=expected * 1e-14)

    def test_reused_energy_matches_independent_formula_for_moving_meshes(self):
        pattern = PatternLibrary().grid(4, 3)
        random = np.random.default_rng(5401)
        for power in (1, 2):
            for shift in ((0., 0.), (1500., -6000.), (1e9, -1e9)):
                energy = QuadEnergy(pattern.quads, power)
                for amount in (.7, .7, .1, .7, .1):
                    xy = pattern.uv * (4., 3.) + shift
                    xy += random.uniform(-amount, amount, xy.shape)
                    # Classify validity independently of the workspace, using
                    # each ordered triple's signed area before scalar energy.
                    areas = []
                    for face in pattern.quads:
                        for corner in range(4):
                            p, q, r = (xy[face[(corner + i) % 4]] for i in range(3))
                            areas.append((q[0]-p[0])*(r[1]-q[1]) - (q[1]-p[1])*(r[0]-q[0]))
                    value = energy.evaluate(xy)
                    if min(areas) <= 1e-10:
                        self.assertEqual(value, float('inf'))
                    else:
                        expected = reference_energy(xy, pattern.quads, power)
                        self.assertAlmostEqual(value, expected, delta=expected * 1e-13)


if __name__ == '__main__':
    unittest.main()
