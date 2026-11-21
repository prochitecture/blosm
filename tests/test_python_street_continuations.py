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
    def test_outward_offset_preserves_endpoint_jets_and_certifies_curvature(self):
        args = ((0., 0.), (2., .04), (1., 0.), (1., .08), 2.)
        before = WidthContinuation(*args, refine_curvature=True)
        after = WidthContinuation(*args, refine_curvature=True, outward_offset=.025)
        t = np.array((0., 1.))
        np.testing.assert_array_equal(after(t), before(t))
        np.testing.assert_allclose(after.derivative(t), before.derivative(t), rtol=0., atol=1e-14)
        np.testing.assert_allclose(after.second_derivative[[0, -1]], 0., atol=1e-13)
        self.assertAlmostEqual(np.linalg.norm(after(.5)-before(.5)), .025)
        self.assertLessEqual(after.curvature_bound, .5)
        with self.assertRaises(ValueError):
            WidthContinuation(*args, refine_curvature=True, outward_offset=.8)

    def test_real_near_straight_width_change_can_use_a_certified_style_continuation(self):
        # Original Bergamo source geometry. Its 4.9-degree through bend and
        # 7-to-3.5m width change have no circular contact in the style range.
        # They still admit a finite smooth curb at the original 2m minimum.
        origin = (-671.226806640625, -760.1535034179688, 0.)
        sections = (
            model.Section(1, 1, (origin, (-613.0177001953125, -754.5712890625, 0.),
                (-597.9947509765625, -753.0810546875, 0.)), width=7.),
            model.Section(2, 2, (origin, (-684.250244140625, -760.2857055664062, 0.),
                (-699.218994140625, -763.4566650390625, 0.),
                (-707.994873046875, -766.027099609375, 0.),
                (-712.1651000976562, -767.7855224609375, 0.),
                (-715.8685913085938, -770.2897338867188, 0.),
                (-718.2806396484375, -773.2060546875, 0.)), width=3.5),
            model.Section(3, 3, ((-696.20361328125, -792.1217651367188, 0.),
                (-695.0670166015625, -785.7098388671875, 0.),
                (-692.8023681640625, -779.1756591796875, 0.),
                (-690.1954956054688, -773.543212890625, 0.),
                (-686.74853515625, -769.0796508789062, 0.),
                (-680.58642578125, -764.7388916015625, 0.), origin), width=3.5))
        junction = model.Intersection(1, (model.End(1, 0), model.End(2, 0), model.End(3, 1)),
                                      radius=2., radius_modes=('STYLE',)*3)
        network = model.Network('near_straight_style_width', sections, (junction,))
        before = network.to_dict()
        from tested_python_streets.backend import PythonMeshBackend
        from tested_python_streets.corner_radii import fitted_radius_report
        update = PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(network.to_dict(), before)
        layout = update.products[junction.layout_key]
        boundary, = [b for b in layout.boundaries if b.continuation is not None]
        self.assertEqual(boundary.corner.fit_policy, 'STYLE_WIDTH_CONTINUATION')
        self.assertEqual(boundary.corner.radius, 2.)
        self.assertLessEqual(boundary.continuation.curvature_bound, .5)
        self.assertIn('STYLE_WIDTH_CONTINUATION',
                      [r['fit_policy'] for r in fitted_radius_report(junction.key, layout)])
        for product in update.products.values():
            if getattr(product, 'mesh', None) is not None:
                self.assertLessEqual(product.mesh.validate()['maximum_edge_ratio'], 8.)
        strict = model.Network('strict_width', sections,
            (model.Intersection(1, junction.ends, radius=2.),))
        with self.assertRaisesRegex(ValueError, 'Requested-radius'):
            PythonMeshBackend(strict.settings).prepare(strict)

    def test_short_real_width_change_refines_certificate_without_moving_mouths(self):
        # Bergamo service junction: the through road changes from 7 to 3.5 m
        # beside a 2 m side street. Exact original controls and fitted style
        # radii demonstrate a conservative-certificate failure, not a reason
        # to shrink the prescribed 2 m back-curb radius or move the mouths.
        origin = (-850.3736572265625, -824.6301879882812, 0.)
        sections = (
            model.Section(1, 1, (origin, (-840.5625, -819.7556762695312, 0.)), width=7.),
            model.Section(2, 2, ((-859.437744140625, -828.6809692382812, 0.), origin), width=3.5),
            model.Section(3, 3, ((-846.4586181640625, -869.1251220703125, 0.),
                (-848.5791625976562, -844.4564208984375, 0.),
                (-850.1173706054688, -828.0365600585938, 0.), origin), width=2.))
        junction = model.Intersection(1, (model.End(1, 0), model.End(2, 1), model.End(3, 1)),
                                      radii=(2., 2., 1.))
        network = model.Network('short_width_continuation', sections, (junction,))
        layout = IntersectionProvider.footprint(footprint_context(network), junction)
        np.testing.assert_allclose(layout.outward_stations,
            (.5869040880496814, 4.692148672141053, 5.118109650396801), atol=1e-7, rtol=0.)
        boundary, = [b for b in layout.boundaries if b.continuation is not None]
        curve = boundary.continuation
        self.assertEqual(boundary.corner.radius, 2.)
        self.assertGreater(curve.curvature_intervals, 1)
        self.assertLessEqual(curve.curvature_bound, .5)
        reference = WidthContinuation(boundary.first, boundary.last,
            curve.first_derivative[0], curve.first_derivative[-1], 1e-6)
        np.testing.assert_allclose(curve.controls, reference.controls, atol=1e-12, rtol=0.)
        from tested_python_streets.backend import PythonMeshBackend
        update = PythonMeshBackend(network.settings).prepare(network)
        from tested_python_streets.ownership import PLAN
        prepared=update.products[PLAN].by_site[junction.key].physical_layout
        original_planes={(a.source.id,a.end):s for a,s in zip(layout.approaches,layout.outward_stations)}
        prepared_planes={(a.source.id,a.end):s for a,s in zip(prepared.approaches,prepared.outward_stations)}
        # The old narrow plane is already feasible. A coarse dimensional
        # reach must not move it from 4.692148672... to 4.725 m.
        self.assertEqual(prepared_planes[2,1],original_planes[2,1])
        dimensional,=[c.dimensional_corner for c in prepared.corners if hasattr(c,'dimensional_corner')]
        self.assertTrue(dimensional.report['fixed_original_contact_candidate'])
        self.assertEqual(dimensional.report['minimum_radius'],2.)
        self.assertLessEqual(prepared.source_cycle_full_width_proof['maximum_error'],.002)
        mesh = update.products[model.Key('junction_mesh', 1)]
        self.assertGreater(len(mesh.mesh.quads), 0)
        self.assertLess(mesh.metadata['quality']['maximum_edge_ratio'], 8.)
        for section in sections:
            self.assertGreater(len(update.products[model.Key('road', section.id)].mesh.quads), 0)

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
