"""Boundary conditions must retain both derivatives of a moving crown."""
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.conditions import ProfileEdge, BoundaryConditions
from tested_python_streets.intersections import Approach
from tested_python_streets.paths import Path, right
from tested_python_streets.profiles import Crown
from tested_python_streets import model
from test_python_street_conditions import inputs
from _python_street_fixtures import fixture_networks


class MovingCrownApproach(Approach):
    """Analytic surface z=h[1-(y-a*x)^2/b^2] on a straight XY path.

    A mouth at x=0 has exactly the ordinary crown's values, but its incoming
    longitudinal derivative is nonzero away from the crest. This distinguishes
    gradient forwarding from a value-only seam or zero-longitudinal shortcut.
    """
    crest_rate = .15

    def physical_gradient(self, station, q):
        q = np.asarray(q)
        across = 4.*self.profile.rise*q/self.width
        return np.stack((-self.crest_rate*across, across), axis=-1)


class PythonProfileGradientTests(unittest.TestCase):
    def ordinary(self, mode=3, end=0):
        section = model.Section(1, 1, ((0., 0., 0.), (20., 2., 0.), (40., -1., 0.)),
                                width=5.5, offset=.4, interpolation=mode)
        return Approach(Path(section), end, Crown(.12, 37., 1))

    def test_ordinary_gradient_preserves_previous_arithmetic_exactly(self):
        for mode in (0, 3):
            for end in (0, 1):
                approach = self.ordinary(mode, end)
                for station in (-.3, 0., .55*approach.length, approach.length):
                    for q in (np.asarray(.2), np.linspace(-1., 1., 17), np.array(((-1., 0.), (.2, 1.)))):
                        normal = approach.frame(station)[2]
                        velocity = approach.jet(station, approach.offset+.5*approach.width*q)[1]
                        covector = right(velocity)
                        denominator = np.sum(covector*(.5*approach.width*normal), axis=-1)
                        previous = approach.profile.derivative(q)[..., None]*covector/denominator[..., None]
                        actual = approach.physical_gradient(station, q)
                        self.assertEqual(actual.dtype, previous.dtype)
                        self.assertEqual(actual.shape, previous.shape)
                        self.assertEqual(actual.tobytes(), previous.tobytes())

    def test_profile_edge_forwards_nonzero_longitudinal_gradient(self):
        section = model.Section(1, 1, ((0., 0., 0.), (20., 0., 0.)), width=6.)
        approach = MovingCrownApproach(Path(section), 0, Crown(.08, 1.7))
        for t in (np.asarray(.2), np.linspace(0., 1., 13), np.array(((0., .3), (.7, 1.)))):
            height, gradient = ProfileEdge(approach, 0., 0., 1.).sample(t)
            xy = approach.point(0., approach.width*(t-.5))
            # The path's right coordinate is world -Y. Use the independent
            # analytic world surface and its derivatives at the same points.
            y = xy[..., 1]
            h = approach.profile.rise
            expected_z = approach.profile.datum+h*(1.-(y/(.5*approach.width))**2)
            expected_y = -2.*h*y/(.5*approach.width)**2
            expected = np.stack((-approach.crest_rate*expected_y, expected_y), axis=-1)
            np.testing.assert_allclose(height, expected_z, rtol=0., atol=1e-14)
            np.testing.assert_allclose(gradient, expected, rtol=0., atol=1e-14)
            if np.ndim(t) == 0:
                self.assertNotEqual(float(gradient[0]), 0.)

    def test_compound_vertices_and_profile_rules_use_the_same_dynamic_gradient(self):
        network = fixture_networks(model)[0]
        group, layouts = inputs(network)
        marker = np.array((.012, -.023))
        ordinary = Approach.physical_gradient
        def supplied(approach, station, q):
            q = np.asarray(q)
            return ordinary(approach, station, q)+(1.-q*q)[..., None]*marker
        with patch.object(Approach, 'physical_gradient', supplied):
            data = BoundaryConditions().build(group, group.loops[0], layouts)
            self.assertEqual(len(data.ports), 3)
            for port in data.ports:
                approach = layouts[port.owner].approaches[port.approach_index]
                station = layouts[port.owner].outward_stations[port.approach_index]
                q = 2.*np.asarray(port.outward_u)-1.
                np.testing.assert_array_equal(data.gradients[list(port.vertices)], supplied(approach, station, q))
                self.assertGreater(float(abs(data.gradients[list(port.vertices)]-ordinary(approach, station, q)).max()), 0.)
            for rule in data.rules:
                if isinstance(rule, ProfileEdge):
                    t = np.linspace(0., 1., 7)
                    _, gradient = rule.sample(t)
                    q = 2.*(rule.first_u+t*(rule.last_u-rule.first_u))-1.
                    np.testing.assert_array_equal(gradient, supplied(rule.approach, rule.station, q))

    def test_regular_gradient_keeps_singular_frame_diagnostic(self):
        approach = self.ordinary()
        with patch.object(approach, 'jet', return_value=(np.zeros(2), np.zeros((3, 2)))):
            with self.assertRaisesRegex(ValueError, 'singular tangent frame'):
                approach.physical_gradient(0., np.array((-1., 0., 1.)))


if __name__ == '__main__':
    unittest.main()
