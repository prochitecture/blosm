import math
import unittest
import numpy as np
import _python_street_package
from test_python_street_compounds import footprint_context
from tested_python_streets import model
from tested_python_streets.continuations import WidthContinuation
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.mesh import cross


def through_junction(widths=(8.1, 5., 5.4), angle=180.3):
    sections = tuple(model.Section(i, i, ((0., 0., 0.),
        (40. * math.cos(math.radians(a)), 40. * math.sin(math.radians(a)), 0.)), width=w)
        for i, (a, w) in enumerate(zip((0., 90., angle), widths)))
    return model.Network('width_continuation', sections,
        (model.Intersection(1, tuple(model.End(s.id, 0) for s in sections)),))


class PythonWidthContinuationTests(unittest.TestCase):
    def test_endpoint_tangents_and_curvature_bound_over_complete_curve(self):
        curve = WidthContinuation((-5., 3.), (5., 4.35), (1., 0.), (1., .01), 2.)
        np.testing.assert_allclose(curve([0., 1.]), ((-5., 3.), (5., 4.35)), atol=1e-14)
        velocity = curve.derivative([0., 1.])
        np.testing.assert_allclose(velocity / np.linalg.norm(velocity, axis=1)[:, None],
            ((1., 0.), np.array((1., .01)) / np.linalg.norm((1., .01))), atol=1e-14)
        t = np.linspace(0., 1., 10001)
        # Independently differentiate the sampled velocity, rather than use the
        # Bernstein numerator used by the production certificate.
        derivative = curve.derivative(t)
        second = np.gradient(derivative, t, axis=0, edge_order=2)
        curvature = abs(cross(derivative, second)) / np.linalg.norm(derivative, axis=1) ** 3
        self.assertLessEqual(float(curvature.max()), curve.curvature_bound + 1e-7)
        self.assertLessEqual(curve.curvature_bound, .5)
        np.testing.assert_allclose(curve.second_derivative[[0, -1]], 0., atol=1e-12)

    def test_sampling_respects_five_mm_error_and_physical_edge_length(self):
        curve = WidthContinuation((-5., 3.), (5., 4.35), (1., 0.), (1., 0.), 2.)
        intervals = curve.intervals(.7, .005)
        t = np.linspace(0., 1., intervals + 1)
        points = curve(t)
        dense = np.linspace(0., 1., 10001)
        chords = np.column_stack([np.interp(dense, t, points[:, i]) for i in (0, 1)])
        self.assertLessEqual(float(np.linalg.norm(curve(dense) - chords, axis=1).max()), .005)
        self.assertLessEqual(float(np.linalg.norm(np.diff(points, axis=0), axis=1).max()), .7)

    def test_local_curvature_certificate_keeps_the_same_valid_polynomial(self):
        a = math.radians(60.)
        arguments = ((-10., 0.), (10., 0.), (math.cos(a), math.sin(a)),
                     (math.cos(a), -math.sin(a)))
        reference = WidthContinuation(*arguments, 1e-6)
        self.assertGreater(reference.curvature_bound, .5)
        curve = WidthContinuation(*arguments, 2., refine_curvature=True)
        np.testing.assert_array_equal(curve.controls, reference.controls)
        self.assertGreater(curve.curvature_intervals, 1)
        self.assertLessEqual(curve.curvature_intervals, 128)
        self.assertLessEqual(curve.curvature_bound, .5)
        t = np.linspace(0., 1., 20001)
        velocity = curve.derivative(t)
        acceleration = np.gradient(velocity, t, axis=0, edge_order=2)
        curvature = abs(cross(velocity, acceleration)) / np.linalg.norm(velocity, axis=1) ** 3
        self.assertLessEqual(float(curvature.max()), curve.curvature_bound + 1e-7)
        with self.assertRaisesRegex(ValueError, 'minimum radius'):
            WidthContinuation(*arguments, 10., refine_curvature=True)  # Actual radius is below 10 m.

    def test_local_certificate_is_rigid_transform_invariant_and_bounded(self):
        a = math.radians(60.)
        points = np.array(((-10., 0.), (10., 0.)))
        tangents = np.array(((math.cos(a), math.sin(a)), (math.cos(a), -math.sin(a))))
        angle = .713
        rotation = np.array(((math.cos(angle), -math.sin(angle)), (math.sin(angle), math.cos(angle))))
        moved = points @ rotation.T + (1234., -5678.)
        directed = tangents @ rotation.T
        curve = WidthContinuation(*moved, *directed, 2., refine_curvature=True)
        reference = WidthContinuation(*points, *tangents, 2., refine_curvature=True)
        # Rotation and evaluation at kilometre coordinates differ by a few
        # absolute-coordinate ulps, without a geometric tolerance change.
        roundoff = 8. * np.spacing(abs(moved).max())
        np.testing.assert_allclose(curve(np.linspace(0., 1., 101)),
            reference(np.linspace(0., 1., 101)) @ rotation.T + (1234., -5678.), atol=roundoff, rtol=0.)
        self.assertLessEqual(curve.curvature_bound, .5)
        self.assertLessEqual(curve.curvature_intervals, 128)

    def test_near_straight_width_step_keeps_turning_circles_and_full_mouths(self):
        network = through_junction()
        layout = IntersectionProvider.footprint(footprint_context(network), network.intersections[0])
        self.assertEqual(sum(b.continuation is not None for b in layout.boundaries), 1)
        self.assertEqual(sum(c.center is not None for c in layout.corners), 2)
        for boundary in layout.boundaries:
            if boundary.continuation is not None:
                self.assertLessEqual(boundary.continuation.curvature_bound, .5)
            else:
                angle = boundary.corner.start_angle + np.linspace(0., 1., 31) * boundary.corner.sweep
                points = boundary.corner.center + boundary.corner.radius * np.column_stack((np.cos(angle), np.sin(angle)))
                np.testing.assert_allclose(np.linalg.norm(points - boundary.corner.center, axis=1), 2., atol=1e-12)
        self.assertEqual(set(layout.mouths), {(s.id, 0) for s in network.sections})

    def test_equal_width_near_straight_approaches_keep_circular_contact_branch(self):
        network = through_junction((6., 6., 6.))
        layout = IntersectionProvider.footprint(footprint_context(network), network.intersections[0])
        self.assertTrue(all(not c.continuation and c.center is not None for c in layout.corners))

    def test_insufficient_reach_is_diagnosed_instead_of_reducing_radius(self):
        with self.assertRaisesRegex(ValueError, 'minimum radius'):
            WidthContinuation((0., 0.), (1., 1.), (1., 0.), (1., 0.), 2.)
        with self.assertRaisesRegex(ValueError, 'turn back'):
            WidthContinuation((0., 0.), (10., 1.), (1., 0.), (-1., 0.), 2.)

    def test_invalid_dimensions_tangents_and_sampling_are_rejected(self):
        for first, tangent, radius in (((0., 0., 0.), (1., 0.), 2.),
                ((0., 0.), (0., 0.), 2.), ((float('nan'), 0.), (1., 0.), 2.),
                ((0., 0.), (1., 0.), float('nan')), ((0., 0.), (1., 0.), 0.)):
            with self.subTest(first=first, tangent=tangent, radius=radius), self.assertRaises(ValueError):
                WidthContinuation(first, (10., 1.), tangent, (1., 0.), radius)
        curve = WidthContinuation((0., 0.), (10., 1.), (1., 0.), (1., 0.), 2.)
        for spacing, tolerance in ((0., .005), (.7, -1.), (float('inf'), .005)):
            with self.assertRaises(ValueError):
                curve.intervals(spacing, tolerance)


if __name__ == '__main__':
    unittest.main()
