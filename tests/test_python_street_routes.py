from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.providers import RoadSpan
from tested_python_streets.routes import route_key
from _python_street_fixtures import seamed_approach_network


def seamed_network(*, reverse_first=False, reverse_next=False, curved=False):
    return seamed_approach_network(model, reverse_first=reverse_first, reverse_next=reverse_next, curved=curved)


class PythonApproachRouteTests(unittest.TestCase):
    def test_trim_crosses_section_seam_without_changing_source_identity(self):
        for reverse_first in (False, True):
            for reverse_next in (False, True):
                with self.subTest(first=reverse_first, next=reverse_next):
                    network = seamed_network(reverse_first=reverse_first, reverse_next=reverse_next)
                    update = PythonMeshBackend(network.settings).prepare(network)
                    route = update.products[route_key(model.End(1, int(reverse_first)))]
                    self.assertEqual([p.section.id for p in route.pieces], [1, 4])
                    self.assertAlmostEqual(route.length, 40.)
                    p, t, n = route.frame(np.array((0., 1., 2., 5., 40.)))
                    np.testing.assert_allclose(p[:, 0], 0., atol=1e-12)
                    np.testing.assert_allclose(p[:, 1], (0., 1., 2., 5., 40.), atol=1e-12)
                    np.testing.assert_allclose(t, np.tile((0., 1.), (5, 1)), atol=1e-12)
                    np.testing.assert_allclose(n, np.tile((1., 0.), (5, 1)), atol=1e-12)
                    consumed = update.products[model.Key('road', 1)]
                    self.assertIsInstance(consumed, RoadSpan)
                    self.assertEqual(consumed.metadata['consumed_by'], ['intersection:10'])
                    road = update.products[model.Key('road', 4)]
                    seam = road.mesh.xyz[road.mesh.seams[1 if reverse_next else 0]]
                    np.testing.assert_allclose(seam[:, 1], 5.25, atol=1e-10)
                    junction = update.products[model.Key('junction_mesh', 10)].mesh
                    distance = np.linalg.norm(junction.xyz[:, None] - seam[None], axis=-1)
                    self.assertLess(float(distance.min(axis=0).max()), 1e-8)
                    self.assertEqual(network.sections[0].id, 1)
                    self.assertEqual(network.sections[3].id, 4)

    def test_continuation_edit_updates_intersection_and_matches_clean_rebuild(self):
        network = seamed_network(curved=True)
        engine = PythonMeshBackend(network.settings)
        initial = engine.prepare(network)
        engine.commit(initial)
        points = ((0., 2., 0.), (1.2, 20., 0.), (0., 40., 0.))
        moved = replace(network.sections[3], points=points)
        edited = replace(network, sections=(*network.sections[:3], moved, network.sections[4]))
        update = engine.prepare(edited)
        full = PythonMeshBackend(network.settings).prepare(edited)
        self.assertIn(model.Key('junction_mesh', 10), update.rebuilt)
        self.assertNotIn(model.Key('road', 5), update.rebuilt)
        self.assertNotEqual(initial.products[model.Key('placement', 81)], update.products[model.Key('placement', 81)])
        for key, part in update.products.items():
            if hasattr(part, 'mesh'):
                np.testing.assert_array_equal(part.mesh.xyz, full.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, full.products[key].mesh.quads)

    def test_incompatible_datum_stops_route_and_preserves_last_good_state(self):
        network = seamed_network()
        engine = PythonMeshBackend(network.settings)
        engine.commit(engine.prepare(network))
        manifest = engine.manifest()
        changed = replace(network.sections[3], datum=1.)
        with self.assertRaises(ValueError):
            engine.prepare(replace(network, sections=(*network.sections[:3], changed, network.sections[4])))
        self.assertEqual(manifest, engine.manifest())

    def test_real_intersection_stops_route_instead_of_becoming_a_seam(self):
        network = seamed_network()
        extra = (model.Section(6, 60, ((0., 2., 0.), (-40., 2., 0.))),)
        junction = model.Intersection(11, (model.End(1, 1), model.End(4, 0), model.End(6, 0)))
        network = replace(network, sections=(*network.sections, *extra),
                          connections=(), intersections=(*network.intersections, junction))
        with self.assertRaises(ValueError):
            PythonMeshBackend(network.settings).prepare(network)


if __name__ == '__main__':
    unittest.main()
