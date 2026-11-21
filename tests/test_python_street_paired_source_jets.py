"""One-sided derivatives follow actual incoming paired curb chord domains."""
from copy import copy
from dataclasses import replace
from functools import lru_cache
import unittest
import numpy as np
import _python_street_package
from test_python_street_merged_entrances import network, updater, context_for
from tested_python_streets.bundle_angular_curves import AngularCoordinates, RayCurb
from tested_python_streets.bundle_merged_entrances import PairedAngularEntrance
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.finite_terminal_entrances import LinearCrown, LinearFrame
from tested_python_streets.intersections import Approach
from tested_python_streets.model import End, Key
from tested_python_streets.ownership import surface_key
from tested_python_streets.paired_source_jet_domains import PairedSourceJetDomain
from tested_python_streets.source_offset_rules import SourceOffsetRule
from tested_python_streets.surface_fields import BoundaryJetField


@lru_cache(maxsize=1)
def fixture():
    source = network(3)
    engine = updater(source)
    update = engine.prepare(source)
    context = context_for(source, update, engine.patterns)
    product = context[surface_key(Key('intersection', 500))]
    approach = next(a for a in product.group.physical_layout.approaches
                    if type(a) is PairedAngularEntrance)
    return context, product, approach


def two_chord_curve():
    frame = object.__new__(AngularCoordinates)
    frame.center, frame.scale, frame.angle, frame.direction = np.zeros(2), 10., 0., 1.
    angles, radii = np.array((-.2, 0., .2)), np.array((10., 11., 10.5))
    return RayCurb(frame, radii[:, None]*np.c_[np.cos(angles), np.sin(angles)])


class PairedSourceJetTests(unittest.TestCase):
    def test_complete_source_condition_field_fits_with_unchanged_tolerance(self):
        context, product, _ = fixture()
        original = context.network.to_dict()
        curves = tuple((curve, curve.points.copy()) for arm in product.group.physical_layout.approaches
                       if type(arm) is PairedAngularEntrance for curve in arm.field.curves)
        data, _, _ = CompoundSurfaceBuilder(context.patterns, source_heights=True).boundary_conditions(product.group, context)
        field = BoundaryJetField(data)
        self.assertLessEqual(field.metadata['maximum_sampled_fit_error'], 1e-8)
        self.assertLess(field.metadata['declared_polynomial_pieces']['prediction_branches'], 20000)
        rules = [rule for boundary in data for rule in boundary.rules if type(rule) is SourceOffsetRule]
        self.assertTrue(any(rule.paired_domain is not None and rule.polynomial_breaks() for rule in rules))
        self.assertEqual(context.network.to_dict(), original)
        for curve, points in curves:
            np.testing.assert_array_equal(curve.points, points)

    def test_actual_chord_endpoint_has_the_correct_one_sided_derivative(self):
        curve = two_chord_curve()
        knot = float(curve.stations[1])
        _, left = curve.jet(knot, interval=(float(curve.stations[0]), knot))
        _, right = curve.jet(knot)
        step = 1e-5
        values = curve.frame.local(curve(np.array((knot-step, knot, knot+step))))[:, 1]
        self.assertAlmostEqual(float(left), float((values[1]-values[0])/step), places=5)
        self.assertAlmostEqual(float(right), float((values[2]-values[1])/step), places=5)
        self.assertGreater(abs(float(left-right)), .1)

    def test_undeclared_chord_crossing_and_outside_queries_are_refused(self):
        curve = two_chord_curve()
        with self.assertRaisesRegex(ValueError, 'undeclared'):
            curve.jet(np.array((-.5, .5)), interval=(-1., 1.))
        with self.assertRaisesRegex(ValueError, 'undeclared'):
            curve.jet(-.1, interval=(0., 1.))
        with self.assertRaises(ValueError):
            curve.jet(0., interval=(float('nan'), 1.))

    def test_one_sided_jet_is_covariant_under_rigid_world_translation(self):
        curve = two_chord_curve()
        angle = .73
        rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        offset = np.array((20000., -30000.))
        frame = copy(curve.frame)
        frame.center = curve.frame.center@rotation.T+offset
        frame.angle += angle
        moved = RayCurb(frame, curve.points@rotation.T+offset)
        # Each declared endpoint uses its own exact source lattice. Its
        # round-trip Cartesian position need not select the same default side.
        for index in (0, 1):
            end = index+1
            _, expected = curve.jet(curve.stations[end], interval=tuple(curve.stations[index:end+1]))
            _, actual = moved.jet(moved.stations[end], interval=tuple(moved.stations[index:end+1]))
            self.assertAlmostEqual(float(expected), float(actual), places=8)

    def test_borrowed_native_path_and_changed_shared_profile_are_refused(self):
        context, _, approach = fixture()
        borrowed = copy(approach)
        borrowed.paths = (copy(approach.paths[0]), approach.paths[1])
        with self.assertRaisesRegex(ValueError, 'original Path'):
            PairedSourceJetDomain.create(context, borrowed)
        changed_end = copy(approach)
        first = approach.original_ends[0]
        changed_end.original_ends = (End(first.section_id, 1-first.end), approach.original_ends[1])
        with self.assertRaisesRegex(ValueError, 'original Path'):
            PairedSourceJetDomain.create(context, changed_end)
        from tested_python_streets.updates import BuildContext
        changed = BuildContext(context.network, dict(context.products), context.patterns)
        key = Key('profile', approach.visits[0].section_id)
        changed.products[key] = replace(changed[key], rise=changed[key].rise+.01)
        with self.assertRaisesRegex(ValueError, 'shared profile'):
            PairedSourceJetDomain.create(changed, approach)

    def test_borrowed_field_frame_is_refused(self):
        context, _, approach = fixture()
        changed = copy(approach)
        changed.field = copy(approach.field)
        changed.field.frame = copy(approach.frame_chart)
        with self.assertRaisesRegex(ValueError, 'incoming field'):
            PairedSourceJetDomain.create(context, changed)

    def test_station_domain_is_finite_and_cannot_skip_a_chord(self):
        context, _, approach = fixture()
        domain = PairedSourceJetDomain.create(context, approach)
        knots = domain.outward_knots()
        first, last = map(float, knots[2:4])
        result = domain.gradient(np.array((first, last)), 1., (first, last))
        self.assertTrue(np.isfinite(result).all())
        with self.assertRaisesRegex(ValueError, 'station piece'):
            domain.gradient(np.array((first-.01,)), 1., (first, last))
        with self.assertRaisesRegex(ValueError, 'finite original'):
            domain.gradient(np.array((0.,)), 1., (-1., 0.))

    def test_ordinary_and_linear_sources_keep_their_existing_interfaces(self):
        context, _, _ = fixture()
        ordinary = Approach(context[Key('path', 300)], 0, context[Key('profile', 300)])
        self.assertIsNone(PairedSourceJetDomain.create(context, ordinary))
        frame = LinearFrame((0., 0.), (1., 0.))
        curves = tuple(np.array(((-2., y), (2., y))) for y in (3., -1., 1., -3.))
        field = LinearCrown(frame, curves, datum=0., rise=.08)
        points = np.array(((0., 3.), (0., -3.), (0., 0.)))
        self.assertTrue(np.isfinite(field.gradient(points)).all())
        np.testing.assert_allclose(field.gradient(points[:2])[:, 1], (-.08*2./3., .08*2./3.))


if __name__ == '__main__':
    unittest.main()
