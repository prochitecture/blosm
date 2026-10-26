"""Independent differential checks for the actual dynamic Bundle crown."""
import unittest
import numpy as np
import _python_street_package
from _python_angular_fixtures import angular_network
from tested_python_streets import model
from tested_python_streets.bundle_angular_curves import AngularCrown
from tested_python_streets.bundle_angular_internal import shared_frame
from tested_python_streets.bundle_angular_ports import ray_jet
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import UpdateCoordinator, BuildContext


def field_for(source):
    patterns = PatternLibrary(source.settings.pattern_cache_bytes)
    prepared = UpdateCoordinator((PathProvider(), ProfileProvider()), patterns).prepare(source)
    context = BuildContext(source, prepared.products, patterns)
    bundle = source.bundles[0]
    frame = shared_frame(context, bundle)
    curves = []
    for member in bundle.members:
        visit, = member.path
        path = context[model.Key('path', visit.section_id)]
        samples = path.samples(0., path.length, source.settings.spacing,
                               tolerance=.1*source.settings.curve_tolerance)
        for side in (-1., 1.):
            q = side*(1 if visit.forward else -1)
            xy = path.point(samples, path.section.offset+q*.5*path.section.width)
            curves.append(xy if visit.forward else xy[::-1])
    field = AngularCrown(frame, curves, datum=source.sections[0].datum, rise=.08)
    knots = np.unique(np.concatenate([curve.stations for curve in field.curves]))
    common = knots[(knots > max(curve.stations[0] for curve in field.curves)) &
                   (knots < min(curve.stations[-1] for curve in field.curves))]
    # Midpoints avoid one-sided Exact chord knots in central-difference checks.
    stations = .5*(common[:-1]+common[1:])
    return field, stations[::max(1, len(stations)//19)]


class PythonAngularGradientTests(unittest.TestCase):
    def test_dynamic_world_gradient_matches_independent_height_differences(self):
        for mode in (0, 3):
            for variable in (False, True):
                source = angular_network(model, 180., variable, mode, widths=(5.5, 6.5),
                                         offsets=(.2, -.1), datum=1.7)
                before = source.to_dict()
                field, station = field_for(source)
                ao, ai, bi, bo = field.edges(station)
                high, low = np.maximum(ao, bi), np.minimum(ai, bo)
                first, last = field.peaks(station)
                # Interior points on both crown halves, plus plateau points
                # where a real separator remains. No clipped curb derivative.
                y = np.stack((.8*high+.2*first, .3*high+.7*first,
                              .5*(first+last), .7*last+.3*low, .2*last+.8*low), axis=-1)
                xy = field.frame.world(np.stack(np.broadcast_arrays(station[:, None], y), axis=-1))
                epsilon = 1e-5
                expected = np.stack([(field.values(xy+axis)-field.values(xy-axis))/(2.*epsilon)
                    for axis in (np.array((epsilon, 0.)), np.array((0., epsilon)))], axis=-1)
                with self.subTest(mode=mode, variable=variable):
                    np.testing.assert_allclose(field.gradient(xy), expected, rtol=0., atol=4e-6)
                    self.assertGreater(float(abs(field.gradient(xy)).max()), 1e-3)
                    self.assertEqual(source.to_dict(), before)

    def test_crests_and_plateaus_have_zero_gradient(self):
        for mode in (0, 3):
            field, station = field_for(angular_network(model, 180., True, mode))
            first, last = field.peaks(station)
            for y in (first, last, .5*(first+last)):
                xy = field.frame.world(np.stack((station, y), axis=-1))
                np.testing.assert_allclose(field.values(xy), field.ceiling, rtol=0., atol=1e-12)
                np.testing.assert_allclose(field.gradient(xy), 0., rtol=0., atol=1e-12)

    def test_one_sided_curb_gradients_keep_zero_height_along_actual_curbs(self):
        for mode in (0, 3):
            field, station = field_for(angular_network(model, 180., True, mode, widths=(5.5, 6.5)))
            for curb in (field.curves[0], field.curves[3]):
                xy = curb(station)
                step = 1e-5
                tangent = (curb(station+step)-curb(station-step))/(2.*step)
                gradient = field.gradient(xy)
                np.testing.assert_allclose(field.values(xy), field.datum, rtol=0., atol=1e-12)
                np.testing.assert_allclose(np.sum(gradient*tangent, axis=-1), 0., rtol=0., atol=1e-8)
                self.assertTrue(np.all(np.linalg.norm(gradient, axis=-1) > 1e-3))
                # The new shared chord jet retains the old port calculation.
                for actual, reference in zip(curb.jet(station), ray_jet(curb, station)):
                    np.testing.assert_array_equal(actual, reference)

    def test_rigid_transform_datum_and_reversed_storage_preserve_physical_gradients(self):
        angle = .73
        matrix = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        translation = np.array((1000., -500.))
        for mode in (0, 3):
            source = angular_network(model, 180., True, mode, widths=(5.5, 6.5))
            field, station = field_for(source)
            ao, ai, bi, bo = field.edges(station)
            y = .7*np.maximum(ao, bi)+.3*np.minimum(ai, bo)
            xy = field.frame.world(np.stack((station, y), axis=-1))
            expected = field.gradient(xy)
            transformed, _ = field_for(angular_network(model, 180., True, mode, widths=(5.5, 6.5),
                reversed_sources=(17, 29), datum=37., rotation=angle, translation=translation))
            actual = transformed.gradient(xy@matrix.T+translation)
            np.testing.assert_allclose(actual, expected@matrix.T, rtol=0., atol=2e-9)


if __name__ == '__main__':
    unittest.main()
