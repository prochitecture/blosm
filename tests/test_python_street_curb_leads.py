"""Width-change curbs may use longer finite leads without changing profiles."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.model import Section, Network, Intersection, End, Key


def width_step(length=60.):
    return Network('short-width-step', (
        Section(1, 1, ((0., 0., 0.), (0., length, 0.)), width=2.),
        Section(2, 2, ((0., 0., 0.), (-length, 0., 0.)), width=2.),
        Section(3, 3, ((0., 0., 0.), (.05*length, -length, 0.)), width=7.)),
        (Intersection(1, tuple(End(i, 0) for i in (1, 2, 3)), 2.),))


class CurbLeadTests(unittest.TestCase):
    def test_width_step_keeps_source_profiles_and_certifies_the_whole_back_curb(self):
        network = width_step()
        original = network.to_dict()
        update = PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(network.to_dict(), original)
        layout = update.products[Key('mouths', 1)]
        boundary, = [b for b in layout.boundaries if b.continuation is not None]
        self.assertEqual(boundary.corner.radius, 2.)
        self.assertLessEqual(boundary.continuation.curvature_bound, .5)
        # A dense independent derivative check also confirms the same radius.
        t = np.linspace(0., 1., 10001)
        velocity = boundary.continuation.derivative(t)
        acceleration = np.gradient(velocity, t, axis=0, edge_order=2)
        numerator = abs(velocity[:, 0]*acceleration[:, 1] - velocity[:, 1]*acceleration[:, 0])
        self.assertLessEqual(float(np.max(numerator/np.linalg.norm(velocity, axis=1)**3)), .5+1e-6)
        mesh = update.products[Key('junction_mesh', 1)].mesh
        stored = mesh.with_positions(mesh.xyz.astype(np.float32).astype(float))
        self.assertLessEqual(stored.validate()['maximum_edge_ratio'], 8.)
        for approach, station in zip(layout.approaches, layout.outward_stations):
            self.assertGreater(station, 0.)
            self.assertLess(station, approach.length)
            q = np.linspace(-1., 1., 2*layout.radial+1)
            expected = np.column_stack((approach.point(station, approach.offset+.5*approach.width*q),
                approach.profile.value(q)*network.settings.exaggeration))
            error = np.linalg.norm(mesh.xyz[:, None] - expected[None], axis=-1).min(axis=0)
            self.assertLess(float(error.max()), 1e-8)
        for corner in layout.corners:
            if corner.center is None:
                continue
            self.assertEqual(corner.radius, 2.)
            for a, s, offset in ((corner.first, corner.stations[0], -.5*corner.first.width),
                                 (corner.second, corner.stations[1], .5*corner.second.width)):
                self.assertAlmostEqual(float(np.linalg.norm(a.point(s, offset)-corner.center)), 2., places=7)

    def test_finite_occupied_source_end_is_never_extrapolated_for_a_width_step(self):
        from tested_python_streets.model import Connection
        network = width_step()
        # The short broad arm meets a separately declared second source road;
        # give that road a scope end so the combined route is also finite.
        short = replace(network.sections[2], points=((0., 0., 0.), (.015, -.3, 0.)))
        tail = Section(4, 4, ((.015, -.3, 0.), (.02, -.4, 0.)), width=7., scope_ends=(1,))
        network = replace(network, sections=(*network.sections[:2], short, tail),
            connections=(Connection(1, (End(3, 1), End(4, 0))),))
        original = network.to_dict()
        with self.assertRaises(ValueError):
            PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(network.to_dict(), original)


if __name__ == '__main__':
    unittest.main()
