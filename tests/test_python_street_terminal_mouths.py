"""Keep a tiny free-end remnant in its junction without removing source data."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import terminal_fixture_networks
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.providers import PathProvider, RoadSpan
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext, UpdateCoordinator, FeatureProvider
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.routes import RouteTopology, Route, route_key
from tested_python_streets.ownership import footprint_key


def footprint(network):
    patterns = PatternLibrary()
    update = UpdateCoordinator((PathProvider(), ProfileProvider()), patterns).prepare(network)
    context = BuildContext(network, update.products, patterns)
    site = network.intersections[0]
    topology = RouteTopology(network)
    for end in site.ends:
        context.products[route_key(end)] = Route.build(context, topology.trace(end))
    return IntersectionProvider.footprint(context, site)


class PythonTerminalMouthTests(unittest.TestCase):
    def test_unresolved_corner_identifies_the_sources_radius_and_available_lengths(self):
        network, = terminal_fixture_networks(model)
        short = replace(network, sections=(replace(network.sections[0], points=((0., 0., 0.), (1., 0., 0.))),
                                           *network.sections[1:]), furniture=())
        with self.assertRaisesRegex(ValueError, r'Section 1 end 0.*Section 2 end 0, radius 2 m, available approaches 1\.000/30\.000 m'):
            PythonMeshBackend(short.settings).prepare(short)

    def test_tiny_terminal_remnant_keeps_its_full_profile_and_original_circular_curbs(self):
        network, = terminal_fixture_networks(model)
        update = PythonMeshBackend(network.settings).prepare(network)
        span = update.products[model.Key('road', 1)]
        self.assertIsInstance(span, RoadSpan)
        self.assertEqual(span.metadata['consumed_by'], ['intersection:10'])
        current = update.products[footprint_key(model.Key('intersection', 10))]
        long = replace(network, sections=(replace(network.sections[0], points=((0., 0., 0.), (30., 0., 0.))), *network.sections[1:]))
        reference = footprint(long)
        for a, b in zip(current.corners, reference.corners):
            self.assertEqual(a.radius, b.radius)
            np.testing.assert_allclose(a.stations, b.stations, atol=1e-8)
            if a.center is not None:
                np.testing.assert_allclose(a.center, b.center, atol=1e-8)
                self.assertAlmostEqual(a.start_angle, b.start_angle, places=8)
                self.assertAlmostEqual(a.sweep, b.sweep, places=8)
        index = next(i for i, a in enumerate(current.approaches) if a.source.id == 1)
        a, distance = current.approaches[index], current.outward_stations[index]
        self.assertAlmostEqual(distance, a.length, places=10)
        q = np.linspace(-1., 1., 2 * current.radial + 1)
        expected = np.column_stack((a.point(distance, .5 * a.width * q),
                                   a.profile.value(q) * network.settings.exaggeration))
        mesh = update.products[model.Key('junction_mesh', 10)].mesh
        self.assertLess(np.linalg.norm(mesh.xyz[:, None] - expected[None], axis=-1).min(axis=0).max(), 1e-8)
        self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)

    def test_extension_restores_a_regular_road_and_incremental_results_match_full(self):
        network, = terminal_fixture_networks(model)
        backend = PythonMeshBackend(network.settings)
        initial = backend.prepare(network)
        backend.commit(initial)
        longer = replace(network, sections=(replace(network.sections[0], points=((0., 0., 0.), (5.168, 0., 0.))), *network.sections[1:]))
        changed = backend.prepare(longer)
        road = changed.products[model.Key('road', 1)]
        self.assertIsNotNone(road.mesh)
        self.assertLess(road.mesh.validate()['maximum_edge_ratio'], 8.)
        fresh = PythonMeshBackend(longer.settings).prepare(longer)
        for key in (model.Key('road', 1), model.Key('junction_mesh', 10)):
            np.testing.assert_array_equal(changed.products[key].mesh.xyz, fresh.products[key].mesh.xyz)
            np.testing.assert_array_equal(changed.products[key].mesh.quads, fresh.products[key].mesh.quads)
        self.assertNotEqual(initial.products[model.Key('placement', 85)], changed.products[model.Key('placement', 85)])
        anchor = changed.products[model.Key('placement', 85)]
        q = 1. - 1. / network.sections[0].width
        self.assertAlmostEqual(anchor.position[2], network.settings.exaggeration * .08 * (1. - q * q), places=10)
        backend.commit(changed)
        restored = backend.prepare(network)
        self.assertIsInstance(restored.products[model.Key('road', 1)], RoadSpan)
        np.testing.assert_array_equal(initial.products[model.Key('junction_mesh', 10)].mesh.xyz,
                                      restored.products[model.Key('junction_mesh', 10)].mesh.xyz)

    def test_dimensional_continuation_keeps_its_incidence_when_the_junction_absorbs_it(self):
        network, = terminal_fixture_networks(model)
        end = network.sections[0].points[-1]
        continued = replace(network, sections=(*network.sections,
            model.Section(4, 1, (end, (30., 0., 0.)), width=7.)),
            connections=(model.Connection(99, (model.End(1, 1), model.End(4, 0))),))
        layout = footprint(continued)
        index = next(i for i, a in enumerate(layout.approaches) if a.source.id == 1)
        approach=layout.approaches[index]
        self.assertGreater(layout.outward_stations[index],end[0])
        self.assertLess(layout.outward_stations[index],approach.length-1.)
        self.assertEqual([p.section.id for p in approach.path.pieces],[1,4])
        self.assertEqual(approach.width,7.)
        self.assertEqual(layout.absorbed_transitions,(model.Key('connection',99),))
        update=PythonMeshBackend(continued.settings).prepare(continued)
        self.assertIsNone(update.products[model.Key('seam_mesh',99)])
        self.assertIsInstance(update.products[model.Key('road',1)],RoadSpan)
        self.assertEqual(update.products[model.Key('seam_layout',99)].absorbed_by,model.Key('intersection',10))

    def test_compatible_seam_absorbs_the_remnant_and_keeps_the_next_section_complete(self):
        network, = terminal_fixture_networks(model)
        end = network.sections[0].points[-1]
        continued = replace(network, sections=(*network.sections,
            model.Section(4, 1, (end, (30., 0., 0.)), width=5.4)),
            connections=(model.Connection(99, (model.End(1, 1), model.End(4, 0))),))
        update = PythonMeshBackend(continued.settings).prepare(continued)
        span = update.products[model.Key('road', 1)]
        self.assertIsInstance(span, RoadSpan)
        self.assertEqual(span.metadata['consumed_by'], ['intersection:10'])
        road = update.products[model.Key('road', 4)]
        self.assertEqual(road.metadata['start'], 0.)
        cap = road.mesh.xyz[road.mesh.seams[0]]
        junction = update.products[model.Key('junction_mesh', 10)].mesh
        self.assertLess(np.linalg.norm(junction.xyz[:, None] - cap[None], axis=-1).min(axis=0).max(), 1e-8)
        self.assertLess(road.mesh.validate()['maximum_edge_ratio'], 8.)

    def test_new_remote_incidence_invalidates_the_former_free_cap(self):
        network, = terminal_fixture_networks(model)
        class RootFootprint(FeatureProvider):
            def tasks(self, source):
                keys = {route_key(e) for e in source.intersections[0].ends}
                keys.add(footprint_key(source.intersections[0].key))
                yield from (t for t in IntersectionProvider().tasks(source) if t.key in keys)
        coordinator = UpdateCoordinator((PathProvider(), ProfileProvider(), RootFootprint()), PatternLibrary())
        initial = coordinator.prepare(network)
        coordinator.commit(initial)
        end = network.sections[0].points[-1]
        connected = replace(network, sections=(*network.sections,
            model.Section(4, 4, (end, (end[0], 30., 0.)), width=5.4),
            model.Section(5, 5, (end, (end[0], -30., 0.)), width=5.4)),
            intersections=(*network.intersections,
                model.Intersection(20, (model.End(1, 1), model.End(4, 0), model.End(5, 0)), radius=2.)))
        update = coordinator.prepare(connected)
        self.assertIn(route_key(model.End(1, 0)), update.rebuilt)
        layout = update.products[footprint_key(network.intersections[0].key)]
        index = next(i for i, a in enumerate(layout.approaches) if a.source.id == 1)
        self.assertFalse(layout.approaches[index].path.terminal_free)
        self.assertAlmostEqual(layout.outward_stations[index], 4.95, places=8)
        coordinator.commit(update)
        remote_radius = replace(connected, intersections=(connected.intersections[0],
                                                        replace(connected.intersections[1], radius=3.)))
        local = coordinator.prepare(remote_radius)
        self.assertNotIn(route_key(model.End(1, 0)), local.rebuilt)
        self.assertNotIn(footprint_key(network.intersections[0].key), local.rebuilt)

    def test_inset_anchor_follows_the_profile_on_either_side_and_outside_offsets_keep_edge_height(self):
        source = model.Section(1, 1, ((0., 0., 0.), (20., 0., 0.)), width=6., datum=2., offset=.3)
        placements = tuple(model.Placement(i, 1, 10., side=side, offset=offset, height=.7)
                           for i, (side, offset) in enumerate((('LEFT', -1.), ('RIGHT', -1.), ('LEFT', 1.), ('RIGHT', 1.))))
        network = model.Network('inset_anchor', (source,), furniture=placements)
        update = PythonMeshBackend(network.settings).prepare(network)
        for i in range(4):
            anchor = update.products[model.Key('placement', i)]
            expected = 2.7 + (10. * .08 * (1. - (2./3.) ** 2) if i < 2 else 0.)
            self.assertAlmostEqual(anchor.position[2], expected, places=10)


if __name__ == '__main__':
    unittest.main()
