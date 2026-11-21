"""Returning OSM Sections keep two incidences and one path/profile resource."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section, Network, Intersection, Connection, End, Key
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.updates import DependencyIndex
from tested_python_streets.paths import Path


def returning_network(interpolation=3):
    return Network('returning', (
        Section(1, 1, ((0, 0, 0), (30, 0, 0), (30, 30, 0), (0, 30, 0), (0, 0, 0)),
                width=3., interpolation=interpolation),
        Section(2, 2, ((0, 0, 0), (-30, -30, 0)), width=3., interpolation=interpolation)),
        (Intersection(1, (End(1, 0), End(1, 1), End(2, 0)), 2.),))


def beijing_returning_network():
    # OSM ways 724548765 (service loop) and 724548762 (attached service).
    # Original lon/lat controls mapped into a local metre plane, then rotated
    # to the actual shared OSM node; no geometric control is inserted/dropped.
    points = np.array(((0., 0., 0.), (-3.074232, -1.870176, 0.),
        (-4.525953, -3.45092, 0.), (-6.41319, -7.825796, 0.),
        (-6.42173, -11.421432, 0.), (-5.12372, -14.928012, 0.),
        (-2.775349, -17.64422, 0.), (1.400484, -19.70364, 0.),
        (6.541283, -19.603452, 0.), (9.641134, -18.156292, 0.),
        (11.682082, -16.30838, 0.), (13.407068, -13.425192, 0.),
        (14.124389, -9.773896, 0.), (13.55224, -6.23392, 0.),
        (11.793096, -3.094696, 0.), (8.180873, -.267168, 0.),
        (3.936725, .712448, 0.)))
    origin = points[-1].copy()
    loop = np.concatenate((points[-1:], points)) - origin
    spur = np.array((origin, (3.595144, 16.25272, 0.))) - origin
    return Network('beijing_returning_service', (
        Section(1, 1, tuple(map(tuple, loop)), width=2.3),
        Section(2, 2, tuple(map(tuple, spur)), width=2.3)),
        (Intersection(1, (End(1, 0), End(1, 1), End(2, 0)), 1.),))


class ReturningSectionTests(unittest.TestCase):
    def test_task_dependencies_are_unique_without_losing_endpoint_incidence(self):
        network = returning_network()
        engine = PythonMeshBackend(network.settings)
        tasks = tuple(task for provider in engine.providers for task in provider.tasks(network))
        DependencyIndex(tasks)
        self.assertEqual(len(network.intersections[0].ends), 3)
        keys = {task.key for task in tasks}
        self.assertIn(Key('approach', 2), keys)
        self.assertIn(Key('approach', 3), keys)
        self.assertEqual(sum(task.key == Key('profile', 1) for task in tasks), 1)

    def test_returning_road_is_trimmed_and_matches_both_junction_mouths(self):
        for network in (returning_network(), returning_network(0), beijing_returning_network()):
            with self.subTest(network=network.uuid, mode=network.sections[0].interpolation):
                before = network.to_dict()
                update = PythonMeshBackend(network.settings).prepare(network)
                self.assertEqual(network.to_dict(), before)
                path = update.products[Key('path', 1)]
                road = update.products[Key('road', 1)]
                junction = update.products[Key('junction_mesh', 1)]
                self.assertGreater(road.metadata['start'], 0.)
                self.assertLess(road.metadata['end'], path.length)
                self.assertEqual(len(road.metadata['junction_spans']), 2)
                for seam in road.mesh.seams.values():
                    distances = np.linalg.norm(junction.mesh.xyz[:, None, :] - road.mesh.xyz[seam][None], axis=2)
                    self.assertLess(float(distances.min(axis=0).max()), 1e-8)
                for part in (road, junction):
                    quality = part.mesh.validate()
                    self.assertLessEqual(quality['maximum_edge_ratio'], 8.)
                    self.assertGreaterEqual(float(part.mesh.xyz[:, 2].min()), -1e-9)
                    self.assertLessEqual(float(part.mesh.xyz[:, 2].max()), .08 * network.settings.exaggeration + 1e-9)

    def test_loop_closure_requires_explicit_incidence(self):
        source = returning_network().sections[0]
        with self.assertRaisesRegex(ValueError, 'explicit start/end incidence'):
            Path(source)

    def test_self_connection_uses_both_endpoint_bend_halves(self):
        source = returning_network().sections[0]
        network = Network('explicit_loop_seam', (source,),
            connections=(Connection(1, (End(1, 0), End(1, 1))),))
        update = PythonMeshBackend(network.settings).prepare(network)
        path = update.products[Key('path', 1)]
        for first, last in zip(path.frame(0.), path.frame(path.length)):
            np.testing.assert_allclose(first, last, rtol=0., atol=1e-10)
        self.assertEqual(set(path.endpoint_bends), {0, 1})

    def test_edit_of_returning_road_matches_a_fresh_build(self):
        network = returning_network()
        engine = PythonMeshBackend(network.settings)
        engine.commit(engine.prepare(network))
        source = network.sections[0]
        points = list(source.points)
        points[2] = (32., 31., 0.)
        edited = replace(network, sections=(replace(source, points=tuple(points)), network.sections[1]))
        incremental = engine.prepare(edited)
        fresh = PythonMeshBackend(network.settings).prepare(edited)
        for key, part in incremental.products.items():
            if hasattr(part, 'mesh'):
                np.testing.assert_array_equal(part.mesh.xyz, fresh.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, fresh.products[key].mesh.quads)

    def test_returning_bulb_cannot_use_the_same_curb_point_as_both_circle_contacts(self):
        from tested_python_streets.providers import PathProvider
        from tested_python_streets.profiles import ProfileProvider
        from tested_python_streets.patterns import PatternLibrary
        from tested_python_streets.updates import UpdateCoordinator
        from tested_python_streets.intersections import Approach, _corner, _local_corner
        from tested_python_streets.corner_radii import CornerContactError
        # Original Bergamo turning bulb: the two loop ends are one physical
        # curb traversed in opposite directions. The old root had s0+s1=L.
        points = ((-463.0650329589844, -691.7221069335938),
            (-476.7802734375, -686.4111328125), (-478.12615966796875, -686.0325317382812),
            (-479.5187072753906, -686.0435791015625), (-480.88018798828125, -686.4442138671875),
            (-482.0783386230469, -687.2011108398438), (-483.0275573730469, -688.2697143554688),
            (-483.6422424316406, -689.5609741210938), (-483.8757629394531, -691.1305541992188),
            (-483.6113586425781, -692.7001342773438), (-482.88018798828125, -694.1028442382812),
            (-481.8611145019531, -695.13818359375), (-480.6008605957031, -695.839599609375),
            (-479.2005310058594, -696.1514282226562), (-477.76904296875, -696.0624389648438),
            (-476.430908203125, -695.5838623046875), (-472.5719909667969, -693.5470581054688),
            (-470.42474365234375, -692.8347778320312), (-467.9429626464844, -692.144775390625),
            (-463.0650329589844, -691.7221069335938))
        loop = Section(1, 1, tuple((*p, 0.) for p in points), width=7.)
        spur = Section(2, 2, ((*points[0], 0.), (-410., -706., 0.)), width=7.)
        network = Network('bulb_coincident_root', (loop, spur),
            (Intersection(1, (End(1, 0), End(1, 1), End(2, 0))),))
        products = UpdateCoordinator((PathProvider(), ProfileProvider()), PatternLibrary()).prepare(network).products
        path, profile = products[Key('path', 1)], products[Key('profile', 1)]
        first, second = Approach(path, 0, profile), Approach(path, 1, profile)
        with self.assertRaisesRegex(CornerContactError, 'contacts coincide'):
            _local_corner(first, second, .75)
        corner = _corner(first, second, .75)
        self.assertLess(sum(corner.stations), path.length - 5.)
        self.assertLess(corner.sweep, -.1)
        contacts = (first.point(corner.stations[0], -3.5),
                    second.point(corner.stations[1], 3.5))
        self.assertGreater(np.linalg.norm(contacts[0] - contacts[1]), .1)
        for approach, station, contact, offset in zip((first, second), corner.stations, contacts, (-3.5, 3.5)):
            vector = contact - corner.center
            self.assertAlmostEqual(float(np.linalg.norm(vector)), .75, places=7)
            self.assertAlmostEqual(float(vector @ approach.jet(station, offset)[1]), 0., places=6)


if __name__ == '__main__':
    unittest.main()
