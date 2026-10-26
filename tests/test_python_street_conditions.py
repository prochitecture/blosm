from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import (fixture_networks, compound_cross_network, short_junction_pair,
                                    seamed_approach_network, compound_seam_fixture)
from test_python_street_compounds import footprints
from tested_python_streets import model
from tested_python_streets.compounds import CompoundGroup, CompoundPlanner, sample_footprint
from tested_python_streets.conditions import BoundaryConditions, CurbEdge
from tested_python_streets.profiles import ProfileProvider, memberships
from tested_python_streets.regions import RegionUnion
from tested_python_streets.backend import PythonMeshBackend


def inputs(network):
    layouts = footprints(network)
    groups = CompoundPlanner().build(layouts)
    if not groups:
        owner = next(iter(layouts))
        sampled = sample_footprint(owner, layouts[owner])
        groups = (CompoundGroup((owner,), (), (sampled,), RegionUnion().build((sampled.region,))),)
    return groups[0], layouts


class PythonBoundaryConditionTests(unittest.TestCase):
    def assert_ports(self, network, data, layouts):
        ownership = memberships(network)
        for port in data.ports:
            section = network.section_index[port.section_id]
            expected = ProfileProvider.build(section, ownership.get(section.id), network.section_index)
            u = np.asarray(port.source_u)
            z = expected.datum + expected.value(2. * u - 1.)
            np.testing.assert_allclose(data.heights[list(port.vertices)], z, rtol=0., atol=1e-12)
            self.assertAlmostEqual(port.outward_u[0], 1.)
            self.assertAlmostEqual(port.outward_u[-1], 0.)
            approach = layouts[port.owner].approaches[port.approach_index]
            station = layouts[port.owner].outward_stations[port.approach_index]
            q = 2. * np.asarray(port.outward_u) - 1.
            velocity = approach.jet(station, approach.offset + .5 * approach.width * q)[1]
            gradient = data.gradients[list(port.vertices)]
            np.testing.assert_allclose(np.sum(gradient * velocity, axis=1), 0., rtol=0., atol=1e-9)

    def test_full_profiles_and_zero_curbs_remain_independent_of_mesh_layout(self):
        network = fixture_networks(model)[0]
        group, layouts = inputs(network)
        data = BoundaryConditions().build(group, group.loops[0], layouts)
        self.assertEqual(len(data.ports), 3)
        self.assert_ports(network, data, layouts)
        for rule in data.rules:
            if isinstance(rule, CurbEdge):
                z, gradient = rule.sample(np.linspace(0., 1., 11))
                np.testing.assert_array_equal(z, 0.)
                self.assertTrue(np.all(np.isfinite(gradient)))
        self.assertFalse(data.heights.flags.writeable)
        self.assertFalse(data.gradients.flags.writeable)

    def test_compound_cross_has_eight_external_ports_and_a_flat_separator_rim(self):
        network = compound_cross_network(model)
        group, layouts = inputs(network)
        data = [BoundaryConditions().build(group, loop, layouts) for loop in group.loops]
        outer, hole = data
        self.assertEqual(len(outer.ports), 8)
        self.assertEqual(hole.ports, ())
        internal = {claim.section_id for claim in group.claims}
        self.assertTrue(internal.isdisjoint({port.section_id for port in outer.ports}))
        np.testing.assert_array_equal(hole.heights, .08)
        np.testing.assert_array_equal(hole.gradients, 0.)
        self.assert_ports(network, outer, layouts)
        for rule in outer.rules:
            if isinstance(rule, CurbEdge):
                z, _ = rule.sample(np.linspace(0., 1., 11))
                self.assertTrue(np.all((abs(z) < 1e-12) | (abs(z - .08) < 1e-12)))

    def test_consumed_short_link_adds_no_internal_profile_constraint(self):
        network = short_junction_pair(model)
        group, layouts = inputs(network)
        data = BoundaryConditions().build(group, group.loops[0], layouts)
        self.assertEqual(len(data.ports), 4)
        self.assertNotIn(1, {p.section_id for p in data.ports})
        self.assert_ports(network, data, layouts)

    def test_seamed_approach_port_uses_the_actual_continuation_section(self):
        for reverse in (False, True):
            network = seamed_approach_network(model, reverse_next=reverse, curved=True)
            group, layouts = inputs(network)
            data = BoundaryConditions().build(group, group.loops[0], layouts)
            port = next(p for p in data.ports if p.section_id == 4)
            self.assertEqual(port.end, int(reverse))
            self.assertNotIn(1, {p.section_id for p in data.ports})
            self.assert_ports(network, data, layouts)

    def test_unequal_rises_obey_the_highest_crown_envelope_and_edge_continuity(self):
        network = fixture_networks(model)[4]
        group, layouts = inputs(network)
        data = BoundaryConditions().build(group, group.loops[0], layouts)
        self.assert_ports(network, data, layouts)
        self.assertAlmostEqual(data.ceiling, .08)
        for i in range(len(data.rules)):
            _, z, g = data.edge_samples(i, np.linspace(0., 1., 17))
            self.assertGreaterEqual(float(z.min()), data.datum - 1e-10)
            self.assertLessEqual(float(z.max()), data.ceiling + 1e-10)
            self.assertAlmostEqual(float(z[0]), data.heights[i])
            self.assertAlmostEqual(float(z[-1]), data.heights[(i + 1) % len(data.heights)])
            self.assertTrue(np.all(np.isfinite(g)))

    def test_unrelated_layouts_cannot_change_datum_ceiling_or_boundary_conditions(self):
        for network in (fixture_networks(model)[0], compound_cross_network(model)):
            group, layouts = inputs(network)
            expected = [BoundaryConditions().build(group, loop, layouts) for loop in group.loops]
            for datum, rise in ((0., 2.), (37., .08)):
                remote = fixture_networks(model)[0]
                remote = replace(remote, sections=tuple(replace(s, datum=datum, rise=rise) for s in remote.sections))
                unrelated = next(iter(footprints(remote).values()))
                expanded = {model.Key('intersection', 999): unrelated, **layouts}
                for loop, before in zip(group.loops, expected):
                    with self.subTest(network=network.uuid, datum=datum, rise=rise, hole=loop.hole):
                        after = BoundaryConditions().build(group, loop, expanded)
                        self.assertEqual((before.datum, before.ceiling), (after.datum, after.ceiling))
                        np.testing.assert_array_equal(before.heights, after.heights)
                        np.testing.assert_array_equal(before.gradients, after.gradients)

    def test_consumed_internal_crowns_cannot_set_the_external_height_ceiling(self):
        base = short_junction_pair(model)
        reference = PythonMeshBackend(base.settings).prepare(base).products[model.Key('junction_mesh', 1)].mesh
        for rise in (0., .04, .18, .5):
            n = replace(base, sections=tuple(replace(s, rise=rise) if s.id == 1 else s for s in base.sections))
            with self.subTest(rise=rise):
                group, layouts = inputs(n)
                data = BoundaryConditions().build(group, group.loops[0], layouts)
                self.assertEqual(data.ceiling, .08)
                self.assertNotIn(1, {port.section_id for port in data.ports})
                result = PythonMeshBackend(n.settings).prepare(n)
                mesh = result.products[model.Key('junction_mesh', 1)].mesh
                self.assertLessEqual(float(mesh.xyz[:, 2].max()), .08 * n.settings.exaggeration + 1e-8)
                np.testing.assert_array_equal(reference.xyz, mesh.xyz)
                np.testing.assert_array_equal(reference.quads, mesh.quads)

    def test_consumed_crown_projection_does_not_relax_the_common_minimum_input_contract(self):
        for base, consumed in ((short_junction_pair(model), {1}), (compound_seam_fixture(model), {1, 6})):
            n = replace(base, sections=tuple(replace(s, datum=1.) if s.id in consumed else s for s in base.sections))
            with self.subTest(network=base.uuid):
                with self.assertRaisesRegex(ValueError, 'common minimum'):
                    PythonMeshBackend(n.settings).prepare(n)

    def test_partial_profile_metadata_is_diagnosed_instead_of_reinterpreted(self):
        network = short_junction_pair(model)
        group, layouts = inputs(network)
        loop = group.loops[0]
        source = {p.region.owner: p for p in group.footprints}
        index = next(i for i, span in enumerate(loop.spans)
            if source[span.owner].segments[span.segment].kind == 'mouth'
            and source[span.owner].segments[span.segment].first == 1.)
        span = replace(loop.spans[index], first=.1)
        clipped = replace(loop, spans=(*loop.spans[:index], span, *loop.spans[index + 1:]))
        with self.assertRaisesRegex(ValueError, 'only part'):
            BoundaryConditions().build(group, clipped, layouts)


if __name__ == '__main__':
    unittest.main()
