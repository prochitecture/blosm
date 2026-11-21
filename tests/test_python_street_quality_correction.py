"""Independent gradient, feasibility, safety and budget checks for planar repair."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh, cross
from tested_python_streets.model import Key
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.quality_correction import QuadShapeCorrection
from tested_python_streets.validation import validate_overlaps


def distorted():
    pattern = PatternLibrary().grid(4, 4)
    xy = pattern.uv * 4.
    center = int(np.argmin(np.linalg.norm(xy - (2., 2.), axis=1)))
    xy[center, 0] = 2.95
    return Mesh(np.column_stack((xy, np.zeros(len(xy)))), pattern.quads)


class PythonQualityCorrectionTests(unittest.TestCase):
    def test_independent_vertex_steps_repair_dense_bands_without_moving_the_seam(self):
        from tested_python_streets.local_quad_grading import grade_dense_steps
        from tested_python_streets.updates import drain_steps
        pattern = PatternLibrary().grid(26,26)
        xy = np.rint(26.*pattern.uv)
        rim = np.flatnonzero(np.any((xy==0.) | (xy==26.), axis=1))
        movable = np.ones(len(xy), dtype=bool)
        movable[rim] = False
        xy[movable & (xy[:,0]%2 == 1),0] += .93
        mesh = Mesh(np.column_stack((xy + (1932.,-721.),np.zeros(len(xy)))),
                    pattern.quads, seams={0:rim})
        before = mesh.xyz.copy()
        attempt = grade_dense_steps(mesh,2000)
        next(attempt)
        attempt.close()
        np.testing.assert_array_equal(mesh.xyz,before)
        result, report = drain_steps(grade_dense_steps(mesh,2000))
        self.assertGreater(report['initial_defects'],256)
        self.assertLessEqual(report['stored_quality']['maximum_edge_ratio'],8.)
        self.assertGreaterEqual(report['stored_quality']['minimum_corner_sine'],.03)
        np.testing.assert_array_equal(result.xyz[rim],before[rim])
        np.testing.assert_array_equal(mesh.xyz,before)

    def test_native_stencil_grade_derivative_includes_profile_height_changes(self):
        from tested_python_streets.native_stencil_alignment import _grade_penalty
        xy=np.array(((0.,0.),(1.7,.1),(1.2,1.),(-.1,.8)))
        z=np.array((0.,.18,.27,.02))
        slopes=np.array(((.01,.02),(.04,-.03),(.05,.01),(-.02,.03)))
        faces=np.array(((0,1,2,3),))
        value,gradient=_grade_penalty(xy,z,slopes,faces,.10)
        self.assertGreater(value,0.)
        for i,axis in np.ndindex(xy.shape):
            delta=np.zeros_like(xy)
            delta[i,axis]=1e-6
            dz=np.sum(slopes*delta,axis=1)
            hi,_=_grade_penalty(xy+delta,z+dz,slopes,faces,.10)
            lo,_=_grade_penalty(xy-delta,z-dz,slopes,faces,.10)
            self.assertAlmostEqual(gradient[i,axis],(hi-lo)/(2e-6),places=6)

    def test_narrow_grass_disk_keeps_the_fixed_tip_and_pavement_requires_its_own_angle_limit(self):
        from tested_python_streets.local_quad_grading import grade
        mesh = Mesh(np.array(((0., 0., 0.), (1., 0., 0.), (1.5, .003, 0.), (.5, .003, 0.))),
                    np.array(((0, 1, 2, 3),)), 'grass', {0: np.arange(4)})
        disk, details = grade(mesh, 1, narrow_grass_disk=True)
        np.testing.assert_array_equal(disk.xyz, mesh.xyz)
        self.assertLessEqual(details['stored_quality']['maximum_edge_ratio'], 8.)
        self.assertGreater(details['stored_quality']['minimum_corner'], 0.)
        mesh.material = 'pavement'
        with self.assertRaises(ValueError):
            grade(mesh, 1, narrow_grass_disk=True)
        with self.assertRaises(ValueError):
            grade(mesh, 1)

    def test_local_grading_derivative_matches_a_directional_difference(self):
        from tested_python_streets.local_quad_grading import _objective
        xy = np.array(((0., 0.), (1., 0.), (1.2, .006), (.43, .05)))
        faces = np.array(((0, 1, 2, 3),))
        free = np.array((1, 2))
        scales = np.array((.25, .48))
        offsets = np.zeros(4)
        direction = np.array((.23, -.06, .15, .1))
        def evaluate(x):
            return _objective(xy, faces, free, scales, x, .06)
        _, gradient = evaluate(offsets)
        step = 1e-6
        difference = (evaluate(offsets+step*direction)[0]-evaluate(offsets-step*direction)[0])/(2*step)
        self.assertAlmostEqual(difference, gradient @ direction, places=6)

    def test_local_collar_separates_incompatible_fixed_edge_lengths(self):
        from tested_python_streets.local_quad_grading import grade
        from tested_python_streets.four_side_charts import quality
        xy = np.array(((1240., 860.), (1241., 860.), (1240.55, 861.), (1240.45, 861.)))
        mesh = Mesh(np.column_stack((xy, np.zeros(4))), np.array([[0, 1, 2, 3]]),
                    seams={0: np.arange(4)})
        repaired, report = grade(mesh, 64)
        np.testing.assert_array_equal(repaired.xyz[:4], mesh.xyz)
        np.testing.assert_array_equal(repaired.seams[0], mesh.seams[0])
        self.assertEqual(len(repaired.quads), 5)
        self.assertEqual(report['inserted_quads'], 4)
        self.assertLessEqual(quality(repaired)['maximum_edge_ratio'], 8.)
        validate_overlaps({Key('chart', 1): repaired})
        with self.assertRaises(ValueError):
            grade(mesh, 4)

    def test_joint_loss_gradient_matches_independent_directional_differences(self):
        xy = np.array(((0., 0.), (10., 0.), (10.2, .002), (-.1, .003)))
        faces = np.array(((0, 1, 2, 3),))
        solver = QuadShapeCorrection()
        value, report, gradient = solver.evaluate(xy, faces, gradient=True)
        losses = solver.face_losses(np.stack((xy,xy[::-1])))
        self.assertAlmostEqual(losses[0],value)
        self.assertTrue(np.isinf(losses[1]))
        self.assertLess(report['minimum_corner_sine'], .05)
        for direction in (np.array(((.2, .03), (-.1, .01), (.07, -.02), (-.05, -.01))),
                          np.array(((0., .02), (0., -.01), (0., .03), (0., -.02)))):
            epsilon = 1e-7
            a = solver.evaluate(xy + epsilon * direction, faces)[0]
            b = solver.evaluate(xy - epsilon * direction, faces)[0]
            numeric = (a - b) / (2. * epsilon)
            self.assertAlmostEqual(float(np.sum(gradient * direction)) / numeric, 1., places=5)

    def test_distorted_grid_becomes_feasible_without_boundary_or_connectivity_changes(self):
        mesh = distorted()
        original = mesh.xyz.copy()
        result, report = QuadShapeCorrection().improve(mesh)
        self.assertTrue(report['qualified'])
        self.assertLess(result.validate()['maximum_edge_ratio'], 8.)
        self.assertTrue(report['exact_fixed_rim'])
        np.testing.assert_array_equal(mesh.xyz, original)
        np.testing.assert_array_equal(result.quads, mesh.quads)
        rim = np.any((original[:, :2] == 0.) | (original[:, :2] == 4.), axis=1)
        np.testing.assert_array_equal(result.xyz[rim], original[rim])
        validate_overlaps({Key('chart', 1): result}, cell_size=1.)
        area = lambda p: .5 * np.sum(cross(p, np.roll(p, -1, axis=1)))
        self.assertAlmostEqual(area(result.xyz[result.quads, :2]), 16., places=10)

    def test_already_qualified_input_keeps_identity_and_exact_buffers(self):
        pattern = PatternLibrary().grid(2, 2)
        mesh = Mesh(np.column_stack((pattern.uv, np.zeros(len(pattern.uv)))), pattern.quads)
        result, report = QuadShapeCorrection().improve(mesh)
        self.assertIs(result, mesh)
        self.assertEqual(report['stop'], 'already_qualified')
        self.assertEqual(report['iterations'], 0)

    def test_semantic_interior_guides_stay_fixed_while_other_vertices_repair_the_chart(self):
        mesh = distorted()
        guide = int(np.argmin(np.linalg.norm(mesh.xyz[:, :2] - (2.95, 2.), axis=1)))
        result, report = QuadShapeCorrection().improve(mesh, fixed_vertices=np.array((guide,)))
        self.assertTrue(report['qualified'])
        self.assertEqual(report['fixed_guides'], 1)
        np.testing.assert_array_equal(result.xyz[guide], mesh.xyz[guide])
        self.assertLessEqual(result.validate()['maximum_edge_ratio'], 8.)
        for invalid in ([-1], [len(mesh.xyz)], [1.5], [[1]], [True]):
            with self.assertRaisesRegex(ValueError, 'guide vertex'):
                QuadShapeCorrection().improve(mesh, fixed_vertices=invalid)

    def test_hole_rim_and_material_are_fixed_along_with_the_exterior(self):
        pattern = PatternLibrary().grid(6, 6)
        xy = pattern.uv * 6.
        centers = xy[pattern.quads].mean(axis=1)
        faces = pattern.quads[~np.all(abs(centers - (2.5, 2.5)) < 1e-8, axis=1)]
        hole = np.flatnonzero(np.all((xy >= (2., 2.)) & (xy <= (3., 3.)), axis=1))
        displaced = int(np.argmin(np.linalg.norm(xy - (4., 4.), axis=1)))
        xy[displaced, 0] = 4.95
        mesh = Mesh(np.column_stack((xy, np.zeros(len(xy)))), faces, 'grass', {7: hole})
        before = mesh.xyz.copy()
        result, report = QuadShapeCorrection().improve(mesh)
        self.assertTrue(report['qualified'])
        np.testing.assert_array_equal(result.xyz[hole], before[hole])
        np.testing.assert_array_equal(result.seams[7], hole)
        self.assertEqual(result.material, 'grass')
        self.assertEqual(len(result.quads), 35)
        p = result.xyz[result.quads, :2]
        self.assertAlmostEqual(.5 * np.sum(cross(p, np.roll(p, -1, axis=1))), 35., places=10)
        validate_overlaps({Key('ring', 1): result}, cell_size=1.)

    def test_fixed_infeasible_rectangle_is_reported_without_weakening_the_limit(self):
        mesh = Mesh(np.array(((0., 0., 0.), (100., 0., 0.), (100., 1., 0.), (0., 1., 0.))),
                    np.array(((0, 1, 2, 3),)))
        result, report = QuadShapeCorrection().improve(mesh)
        self.assertFalse(report['qualified'])
        np.testing.assert_array_equal(result.xyz, mesh.xyz)
        self.assertEqual(result.validate()['maximum_edge_ratio'], 100.)

    def test_invalid_settings_budget_and_nonplanar_input_are_rejected(self):
        for kw in (dict(max_iterations=0), dict(max_cells=0), dict(target_ratio=8.),
                   dict(minimum_sine=0.), dict(target_sine=2.), dict(target_ratio=float('nan'))):
            with self.assertRaises(ValueError):
                QuadShapeCorrection(**kw)
        with self.assertRaisesRegex(ValueError, 'cell budget'):
            QuadShapeCorrection(max_cells=1).improve(distorted())
        mesh = distorted();mesh.xyz[12, 2] = .1
        with self.assertRaisesRegex(ValueError, 'planar'):
            QuadShapeCorrection().improve(mesh)

    def test_iteration_budget_is_finite_and_original_input_is_never_mutated(self):
        mesh = distorted();before = mesh.xyz.copy()
        result, report = QuadShapeCorrection(max_iterations=1).improve(mesh)
        self.assertEqual(report['iterations'], 1)
        self.assertLessEqual(report['line_search_evaluations'], 24)
        self.assertGreater(result.validate()['minimum_corner'], 0.)
        np.testing.assert_array_equal(mesh.xyz, before)


if __name__ == '__main__':
    unittest.main()
