"""Dimensional seam geometry, source ownership and incremental regeneration."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section, Connection, End, Network, Settings, Key, Placement, fingerprint
from tested_python_streets.paths import Path
from tested_python_streets.intersections import Approach
from tested_python_streets.profiles import Crown
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.transitions import SectionTransition, TransitionProduct, required
from tested_python_streets.backend import PythonMeshBackend


def sources(*, reverse=False, offset=0., rises=(.08, .08), length=40.):
    a = Section(1, 7, ((-length, 0., 0.), (0., 0., 0.)), width=5.4, rise=rises[0])
    b = Section(2, 7, ((0., 0., 0.), (length, 0., 0.)), width=8.1, offset=offset, rise=rises[1])
    if reverse:
        b = replace(b, points=b.points[::-1], offset=-offset)
    seam = Connection(6, (End(1, 1), End(2, int(reverse))))
    return Network('dimensional_seam', (a, b), connections=(seam,))


def construct(network):
    context = SimpleNamespace(network=network, patterns=PatternLibrary())
    site = network.connections[0]
    approaches = tuple(Approach(Path(network.section_index[e.section_id]), e.end,
        Crown(network.section_index[e.section_id].rise, network.section_index[e.section_id].datum)) for e in site.ends)
    plan = SectionTransition.solve(context, site, approaches)
    return context, plan, TransitionProduct.build(context, plan)


class PythonSectionTransitionPrimitiveTests(unittest.TestCase):
    def test_width_transition_keeps_both_complete_mouth_profiles_and_zero_outer_curbs(self):
        network = sources()
        before = fingerprint(network)
        context, plan, product = construct(network)
        mesh = product.pavement.mesh
        self.assertTrue(required(network, network.connections[0]))
        self.assertEqual(fingerprint(network), before)
        self.assertLessEqual(plan.certification['physical_slopes']['maximum'], plan.certification['maximum_slope'] + 1e-8)
        self.assertLessEqual(plan.certification['maximum_entrance_longitudinal_slope'], .03 + 1e-8)
        for i, approach in enumerate(plan.approaches):
            q = np.linspace(-1., 1., plan.across + 1)
            outward_q = -q if i == 0 else q
            row = mesh.xyz[mesh.seams[i]]
            np.testing.assert_allclose(row[:, :2], approach.point(plan.distances[i], approach.offset + .5 * approach.width * outward_q), atol=1e-12)
            np.testing.assert_allclose(row[:, 2], 10. * approach.profile.value(outward_q), atol=1e-12)
        for side in (-1., 1.):
            t = np.linspace(0., 1., plan.rows + 1)
            self.assertEqual(np.max(abs(plan.surface(t, side)[1])), 0.)
        self.assertGreater(plan.certification['minimum_jacobian_bound'], 0.)
        self.assertLess(mesh.validate()['maximum_edge_ratio'], 8.)

    def test_reversing_a_source_does_not_change_the_physical_surface(self):
        _, first, product = construct(sources(offset=.7, rises=(.06, .1)))
        _, reversed_plan, reversed_product = construct(sources(reverse=True, offset=.7, rises=(.06, .1)))
        np.testing.assert_allclose(product.pavement.mesh.xyz, reversed_product.pavement.mesh.xyz, atol=1e-12)
        np.testing.assert_array_equal(product.pavement.mesh.quads, reversed_product.pavement.mesh.quads)
        for t in (.0, .2, .5, .9, 1.):
            xy, z = first.surface(np.array((t,)), np.array((0.,)))
            self.assertAlmostEqual(float(product.height(xy)[0]), 10. * float(z[0]), delta=.002)
        bounds = product.pavement.mesh.bounds
        product.release_geometry()
        self.assertIsNone(product.pavement.mesh)
        self.assertEqual(product.pavement.metadata['bounds'], bounds)
        product.release_geometry()
        self.assertTrue(np.isfinite(product.height(first.surface(np.array((.5,)), np.array((0.,)))[0])).all())

    def test_compatible_seam_needs_no_transition_and_short_or_conflicting_data_is_rejected(self):
        network = sources()
        compatible = replace(network, sections=(network.sections[0], replace(network.sections[1], width=5.4)))
        self.assertFalse(required(compatible, compatible.connections[0]))
        with self.assertRaisesRegex(ValueError, 'No bounded dimensional'):
            construct(sources(length=.1, offset=10.))
        conflicting = replace(network, sections=(network.sections[0], replace(network.sections[1], datum=1.)))
        with self.assertRaisesRegex(ValueError, 'common minimum'):
            construct(conflicting)


class PythonSectionTransitionIntegrationTests(unittest.TestCase):
    def test_exhausted_clearance_is_owned_by_the_junction_without_changing_sources(self):
        from _python_street_fixtures import clearance_fixture_networks
        from tested_python_streets import model
        network = clearance_fixture_networks(model)[0]
        backend = PythonMeshBackend(network.settings)
        backend.commit(backend.prepare(network))
        sections = tuple(replace(s, points=((0., 0., 0.), (2.5, 0., 0.), (5., 0., 0.))) if s.id == 1
                         else replace(s, points=((5., 0., 0.), *s.points[1:])) if s.id == 2 else s for s in network.sections)
        changed = replace(network, sections=sections, furniture=())
        snapshot = fingerprint(changed)
        update = backend.prepare(changed)
        self.assertIsNone(update.products[Key('seam_mesh',6)])
        self.assertEqual(update.products[Key('seam_layout',6)].absorbed_by,Key('intersection',1))
        self.assertIsNone(getattr(update.products[Key('road',1)],'mesh',None))
        self.assertLess(update.products[Key('junction_mesh',1)].mesh.validate()['maximum_edge_ratio'],8.)
        self.assertEqual(fingerprint(changed), snapshot)

    def test_clearance_follows_a_mouth_beyond_a_compatible_seam_in_both_source_directions(self):
        from _python_street_fixtures import clearance_fixture_networks
        from tested_python_streets import model
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                network = clearance_fixture_networks(model)[0]
                short = replace(network.sections[0], points=((0., 0., 0.), (.5, 0., 0.), (1., 0., 0.)))
                continuation = replace(network.sections[0], id=6,
                    points=((1., 0., 0.), (3.75, 0., 0.), (7.5, 0., 0.)))
                if reverse:
                    continuation = replace(continuation, points=continuation.points[::-1])
                network = replace(network, sections=(short, *network.sections[1:], continuation),
                    connections=(Connection(6, (End(6, int(not reverse)), End(2, 0))),
                                 Connection(7, (End(1, 1), End(6, int(reverse))))),
                    furniture=(replace(network.furniture[0], section_id=6,
                                       station=.7 if reverse else 5.8, side='LEFT' if reverse else 'RIGHT'),))
                result = PythonMeshBackend(network.settings).prepare(network)
                self.assertIn(Key('physical_intersection', 1), result.dependencies[Key('seam_layout', 6)])
                self.assertIsNone(getattr(result.products[Key('road', 1)], 'mesh', None))
                road = result.products[Key('road', 6)]
                self.assertAlmostEqual(road.metadata['end'] - road.metadata['start'], .25, places=9)
                self.assertTrue(np.all(np.isfinite(result.products[Key('placement', 81)].position)))

    def test_junction_clearance_preserves_a_complete_profile_strip_and_declares_its_dependency(self):
        from _python_street_fixtures import clearance_fixture_networks
        from tested_python_streets import model
        network = clearance_fixture_networks(model)[0]
        result = PythonMeshBackend(network.settings).prepare(network)
        self.assertIn(Key('physical_intersection', 1), result.dependencies[Key('seam_layout', 6)])
        road = result.products[Key('road', 1)]
        self.assertAlmostEqual(road.metadata['end'] - road.metadata['start'], .25, places=9)
        plan = result.products[Key('seam_surface', 6)].plan
        self.assertEqual(plan.rows % 2, 0)
        for approach in plan.approaches:
            rim = plan.source_half_rim(approach.source.id, 1., network.settings.exaggeration)
            self.assertEqual((len(rim) - 1) % 2, 0)
        constraint = next(c for c in plan.certification['clearance_constraints'] if c['section'] == 1)
        self.assertAlmostEqual(constraint['clearance'], .25)
        profile = result.products[Key('profile', 1)]
        for row in road.mesh.seams.values():
            np.testing.assert_allclose(road.mesh.xyz[row, 2], 10. * profile.value(np.linspace(-1., 1., len(row))), atol=1e-12)

    def test_editing_only_the_neighboring_junction_rebuilds_the_transition_and_matches_a_clean_build(self):
        from _python_street_fixtures import clearance_fixture_networks
        from tested_python_streets import model
        network = clearance_fixture_networks(model)[0]
        backend = PythonMeshBackend(network.settings)
        first = backend.prepare(network)
        backend.commit(first)
        changed = replace(network, sections=tuple(replace(s, width=s.width + .4) if s.id == 3 else s for s in network.sections))
        edited = backend.prepare(changed)
        self.assertIn(Key('seam_surface', 6), edited.rebuilt)
        self.assertIs(first.products[Key('path', 1)], edited.products[Key('path', 1)])
        self.assertIs(first.products[Key('road', 5)], edited.products[Key('road', 5)])
        self.assertFalse(np.array_equal(first.products[Key('seam_mesh', 6)].mesh.xyz,
                                       edited.products[Key('seam_mesh', 6)].mesh.xyz))
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        for key, part in edited.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, fresh.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, fresh.products[key].mesh.quads)
        self.assertEqual(edited.products[Key('placement', 81)], fresh.products[Key('placement', 81)])

    def test_curved_seam_and_rigid_transforms_keep_the_same_physical_geometry(self):
        from _python_street_fixtures import transition_fixture_networks
        from tested_python_streets import model
        network = transition_fixture_networks(model)[2]
        original = PythonMeshBackend(network.settings).prepare(network).products[Key('seam_mesh', 6)].mesh
        angle = .731
        rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        offset = np.array((1234., -5678.))
        sections = []
        for source in network.sections:
            points = np.asarray(source.points).copy()
            points[:, :2] = points[:, :2] @ rotation.T + offset
            sections.append(replace(source, points=tuple(map(tuple, points))))
        moved_network = replace(network, sections=tuple(sections))
        moved = PythonMeshBackend(network.settings).prepare(moved_network).products[Key('seam_mesh', 6)].mesh
        # Angular sorting can swap the two ends. Compare the same surface in
        # source order rather than depending on an arbitrary first direction.
        actual = (moved.xyz[:, :2] - offset) @ rotation
        expected = original.xyz[:, :2]
        if np.linalg.norm(actual[0] - expected[0]) > 1e-7:
            actual = actual.reshape(-1, 17, 2)[::-1, ::-1].reshape(-1, 2)
        np.testing.assert_allclose(actual, expected, atol=2e-7, rtol=0.)

    def test_default_backend_trims_roads_and_publishes_the_seam_with_exact_interfaces(self):
        network = sources()
        result = PythonMeshBackend(network.settings).prepare(network)
        transition = result.products[Key('seam_surface', 6)]
        self.assertIs(result.products[Key('seam_mesh', 6)], transition.pavement)
        for i, approach in enumerate(transition.plan.approaches):
            road = result.products[Key('road', approach.source.id)]
            row = road.mesh.xyz[road.mesh.seams[approach.end]]
            expected = transition.pavement.mesh.xyz[transition.pavement.mesh.seams[i]]
            same_direction = (i == 0 and approach.end == 1) or (i == 1 and approach.end == 0)
            np.testing.assert_allclose(row, expected if same_direction else expected[::-1], atol=1e-12)
            self.assertEqual(len(road.metadata['transition_spans']), 1)

    def test_anchor_inside_a_transition_uses_its_surface_after_release(self):
        network = replace(sources(rises=(.06, .1)), furniture=(Placement(81, 1, 39., offset=-1.),))
        backend = PythonMeshBackend(network.settings)
        initial = backend.prepare(network)
        owner = initial.products[Key('seam_surface', 6)]
        anchor = initial.products[Key('placement', 81)]
        self.assertAlmostEqual(anchor.position[2], owner.height(np.array((anchor.position[:2],)))[0], places=12)
        backend.commit(initial)
        backend.release_published_meshes()
        changed = replace(network, furniture=(replace(network.furniture[0], station=39.2),))
        edited = backend.prepare(changed)
        self.assertEqual(edited.rebuilt, (Key('placement', 81),))
        self.assertNotEqual(anchor.position, edited.products[Key('placement', 81)].position)
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        self.assertEqual(edited.products[Key('placement', 81)], fresh.products[Key('placement', 81)])

    def test_changing_to_compatible_dimensions_retires_only_the_transition_area(self):
        network = sources()
        backend = PythonMeshBackend(network.settings)
        backend.commit(backend.prepare(network))
        changed = replace(network, sections=(network.sections[0], replace(network.sections[1], width=5.4)))
        result = backend.prepare(changed)
        self.assertIn(Key('seam_mesh', 6), result.removed)
        self.assertNotIn(Key('seam_surface', 6), result.products)
        self.assertEqual({k.kind for k in result.products}, {'path', 'profile', 'road', 'seam_layout'})
        self.assertAlmostEqual(result.products[Key('road', 1)].metadata['end'], 40.)

    def test_bundle_rims_include_the_owned_transition_and_no_separator_gap(self):
        from _python_street_fixtures import transition_fixture_networks
        from tested_python_streets import model
        network = transition_fixture_networks(model)[1]
        result = PythonMeshBackend(network.settings).prepare(network)
        self.assertIn(Key('separator', 7), result.products)
        self.assertTrue(result.products[Key('separator', 7)].mesh.validate()['minimum_corner'] > 0.)
        self.assertEqual(result.products[Key('bundle_plan', 7)].mode, 'separated')
        from tested_python_streets.bundles import Member
        context = type('Context', (), {'network': network, '__getitem__': lambda self, key: result.products[key]})()
        first = Member.build(context, network.bundles[0].members[0], 0)
        # Centerline checks retain the whole explicit member route, including
        # the two source ranges whose pavement belongs to the transition.
        self.assertEqual(len(first.path), 2)
        np.testing.assert_allclose(first.path[0].curves[-1][-1], (0., 8.), atol=1e-12)
        np.testing.assert_allclose(first.path[1].curves[0][0], (0., 8.), atol=1e-12)
        self.assertTrue(np.all(np.diff(first.rim[:, 0]) > 0.))

    def test_connection_and_intersection_identity_spaces_do_not_alias_transition_products(self):
        from tested_python_streets.model import Intersection
        network = sources()
        shifted = tuple(replace(s, id=s.id + 2, street_id=8,
            points=tuple((x, y + 30., z) for x, y, z in s.points)) for s in network.sections)
        network = replace(network, sections=(*network.sections, *shifted),
            intersections=(Intersection(6, (End(3, 1), End(4, 0))),))
        result = PythonMeshBackend(network.settings).prepare(network)
        self.assertIn(Key('seam_surface', 6), result.products)
        self.assertIn(Key('profile_transition', 6), result.products)
        self.assertIsNot(result.products[Key('seam_surface', 6)], result.products[Key('profile_transition', 6)])
        a, b = result.products[Key('seam_mesh', 6)].mesh, result.products[Key('junction_mesh', 6)].mesh
        np.testing.assert_allclose(b.xyz - a.xyz, np.tile((0., 30., 0.), (len(a.xyz), 1)), atol=1e-12)


if __name__ == '__main__':
    unittest.main()
