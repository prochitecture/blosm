"""Once-owned compound inputs from shared physical connector rings."""
from copy import deepcopy
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest

from test_physical_junctions import snapshot, sections, module as physical, contract as physical_contract
from test_bounded_bundle_network import inputs

contract = NS(**vars(physical_contract), STREET_ITEM=inputs.STREET_ITEM,
              BUNDLE_ITEM=inputs.BUNDLE_ITEM, validate=inputs.validate)

spec = importlib.util.spec_from_file_location('tested_physical_clusters',
    Path(__file__).resolve().parents[1]/'way/street_network/physical_clusters.py')
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module; spec.loader.exec_module(module)


def fixture(reverse=False):
    result = snapshot(reverse)
    # This fixture's three physical sites have known CCW incidence orders.
    # The origin recorder tests use arbitrary rings; those are not an embedding.
    for site, order in zip(result.junction_origins['physical'], ((1, 0, 2), (2, 0, 1, 3), (2, 1, 0))):
        site['ring'] = tuple(site['ring'][i] for i in order)
    # The second Bundle points south, so its authored lateral order reverses.
    for member in result.bundles[2]['topology']['members']:
        member['head_order'] = 1-member['head_order']
        member['tail_order'] = 1-member['tail_order']
    connector = result.junction_origins['logical'][1]['junction']['ring'][0]
    connector['endpoints'] = connector['endpoints'][::-1]
    return result


class PhysicalClusterTests(unittest.TestCase):
    def setUp(self):
        self.physical = physical.PreparedPhysicalJunctions(contract, sections, {'residential': {'radius': 2.}})
        self.reader = module.PreparedPhysicalClusters(self.physical)

    def test_shared_component_has_one_owner_and_preserves_complete_sources(self):
        for reverse in (False, True):
            source = fixture(reverse); before = deepcopy(source)
            packet = self.reader.prepare(source, 100)
            self.assertEqual(packet, self.reader.prepare(source, 200))
            self.assertEqual(source, before)
            self.assertEqual(packet.owner_id, 100)
            self.assertEqual(packet.physical_junction_ids, (10, 20, 30))
            self.assertEqual(packet.logical_junction_ids, (100, 200))
            self.assertEqual([s.source_id for s in packet.internal_streets], [5, 6])
            self.assertEqual(len(packet.junction.connectors), 4)
            self.assertEqual(sorted(s.source_id for s in packet.junction.streets), [1, 2, 3, 4, 7, 8])
            for street in (*packet.junction.streets, *packet.internal_streets):
                self.assertEqual(street.sections, sections(source, street.source_id))
            self.assertEqual(len(packet.external_ring), 6)
            all_incidents = list(packet.external_ring)+[x for _, ends in packet.internal_incidents for x in ends]
            self.assertEqual(len(all_incidents), 10)
            self.assertEqual(len(set(all_incidents)), 10)
            self.assertEqual([p.separator for p in packet.junction.profiles].count(True), 2)
            self.assertTrue(all(not p.half_profile for p in packet.junction.profiles))
            self.assertTrue(all(c.radius == 2. for c in packet.junction.connectors))

    def test_traversal_preserves_cyclic_order_under_rotated_input_rings(self):
        original = self.reader.prepare(fixture(), 100)
        source = fixture()
        for site in source.junction_origins['physical']:
            site['ring'] = site['ring'][1:]+site['ring'][:1]
        candidate = self.reader.prepare(source, 200)
        self.assertEqual(original, candidate)

    def test_curved_internal_street_is_not_replaced_with_a_chord(self):
        source = fixture()
        source.streets[5]['items'][0]['points'] = ((-20., 0.), (-12., 4.), (-5., 3.), (0., 0.))
        packet = self.reader.prepare(source, 100)
        self.assertEqual(packet.internal_streets[0].sections[0].centerline,
                         source.streets[5]['items'][0]['points'])

    def test_no_angle_or_coordinate_pairing_for_interleaved_members(self):
        source = fixture()
        # Keep all IDs, endpoints and positions; change only the recorded ring.
        site = source.junction_origins['physical'][1]
        site['ring'] = site['ring'][::-1]
        with self.assertRaisesRegex(ValueError, 'interleaved'):
            self.reader.prepare(source, 100)

    def test_missing_ambiguous_or_damaged_logical_provenance_rejects(self):
        for damage in ('missing', 'duplicate', 'duplicate_internal', 'missing_internal', 'terminal', 'endpoint_set'):
            source = fixture(); logical = source.junction_origins['logical'][0]
            endpoints = list(logical['endpoints'])
            if damage == 'missing':
                logical['endpoints'] = tuple(endpoints[1:])
            elif damage == 'duplicate':
                logical['endpoints'] = (*endpoints, endpoints[0])
            elif damage == 'duplicate_internal':
                internal = next(e for e in endpoints if e['endpoint']['street_id'] == 5)
                logical['endpoints'] = (*endpoints, internal)
            elif damage == 'missing_internal':
                internal = next(e for e in endpoints if e['endpoint']['street_id'] == 5)
                logical['endpoints'] = tuple(e for e in endpoints if e is not internal)
            elif damage == 'terminal':
                endpoints[0]['endpoint']['terminal_id'] = -10
            else:
                logical['junction']['ring'][0]['endpoints'] = logical['junction']['ring'][0]['endpoints'][:1]
            with self.subTest(damage=damage), self.assertRaises(ValueError):
                self.reader.prepare(source, 100)

    def test_disconnected_or_cyclic_internal_graph_does_not_fill_an_island(self):
        source = fixture(); packet = self.physical.prepare(source, 100)
        # Feed an already transported disconnected component to the graph stage.
        reader = module.PreparedPhysicalClusters(NS(contract=contract, rise=.08, cell_size=.7,
            prepare=lambda *args: physical.PhysicalJunctionInputs(packet.logical_junction_ids,
                (packet.junctions[0], packet.junctions[2]), packet.endpoint_bundles)))
        with self.assertRaisesRegex(ValueError, 'connected tree'):
            reader.prepare(source, 100)
        # Add a third internal edge around the three sites. It encloses an
        # island; a single compound disk must not swallow it.
        street = contract.StreetInput(9, (contract.SectionInput(90, ((-20., 0.), (20., 0.)), 6.),))
        sites = list(packet.junctions)
        for i, leaving in ((0, True), (2, False)):
            site = sites[i]
            connector = contract.ConnectorInput(900+i, (contract.MemberInput(len(site.streets), leaving),), 0)
            sites[i] = replace(site, streets=(*site.streets, street), connectors=(*site.connectors, connector))
        cyclic = replace(packet, junctions=tuple(sites),
            endpoint_bundles=(*packet.endpoint_bundles, (10, 900, -1), (30, 902, -1)))
        reader = module.PreparedPhysicalClusters(NS(contract=contract, rise=.08, cell_size=.7,
                                                   prepare=lambda *args: cyclic))
        with self.assertRaisesRegex(ValueError, 'cycles retain islands'):
            reader.prepare(source, 100)

    def test_different_widths_and_declared_corner_radii_survive_pairing(self):
        source = fixture()
        source.streets[1]['items'][0].update(width=7.5, category='wide')
        source.streets[2]['items'][0].update(width=4.5, category='narrow')
        reader = module.PreparedPhysicalClusters(physical.PreparedPhysicalJunctions(contract, sections,
            {'residential': {'radius': 2.}, 'wide': {'radius': 5.}, 'narrow': {'radius': 3.}}))
        packet = reader.prepare(source, 100)
        pair = next(c for c in packet.junction.connectors if c.bundle_source_id == 1)
        self.assertEqual([packet.junction.streets[m.street_index].sections[0].width for m in pair.members], [7.5, 4.5])
        self.assertEqual(pair.radius, 5.)  # Original last exterior incidence is Street 1.

    def test_internal_bundle_member_is_not_absorbed_as_ordinary_pavement(self):
        source = fixture(); original = self.physical.prepare(source, 100)
        roles = tuple((pid, cid, 1 if cid == 101 else bid) for pid, cid, bid in original.endpoint_bundles)
        altered = replace(original, endpoint_bundles=roles)
        reader = module.PreparedPhysicalClusters(NS(contract=contract, rise=.08, cell_size=.7,
                                                   prepare=lambda *args: altered))
        with self.assertRaisesRegex(ValueError, 'separator ownership'):
            reader.prepare(source, 100)

    def test_unsupported_original_source_is_not_dropped(self):
        source = fixture(); source.streets[5]['status'] = 'UNSUPPORTED'
        with self.assertRaisesRegex(ValueError, 'unsupported original Street'):
            self.reader.prepare(source, 100)

    def test_corridor_profiles_are_derived_only_inside_the_compound_region(self):
        for reverse in (False, True):
            source = fixture(reverse); before = deepcopy(source)
            reader = module.PreparedPhysicalCorridors(self.reader)
            result = reader.prepare(source, 100)
            self.assertEqual(source, before)
            self.assertEqual(result, reader.prepare(source, 200))
            original = {j.source_id: j for j in self.physical.prepare(source, 100).junctions}
            internal = dict(result.internal_profiles)
            self.assertEqual(set(internal), {5, 6})
            for junction in result.junctions:
                baseline = original[junction.source_id]
                self.assertEqual(junction.streets, baseline.streets)
                self.assertEqual(junction.connectors, baseline.connectors)
                for connector in junction.connectors:
                    member, = connector.members
                    street = junction.streets[member.street_index]
                    profile = junction.profiles[connector.profile_index]
                    if street.source_id in internal:
                        self.assertEqual(profile.half_profile, internal[street.source_id]*(1 if member.leaving else -1))
                    else:
                        self.assertEqual(profile, baseline.profiles[connector.profile_index])
                signs = [junction.profiles[c.profile_index].half_profile for c in junction.connectors]
                self.assertTrue(all((a == -1) == (b == 1) for a, b in zip(signs, signs[1:]+signs[:1])))

    def test_incompatible_external_median_roles_are_not_overwritten(self):
        source = fixture()
        for member in source.bundles[2]['topology']['members']:
            member['head_order'] = 1-member['head_order']
            member['tail_order'] = 1-member['tail_order']
        before = deepcopy(source)
        with self.assertRaisesRegex(ValueError, 'incompatible median'):
            module.PreparedPhysicalCorridors(self.reader).prepare(source, 100)
        self.assertEqual(source, before)


if __name__ == '__main__':
    unittest.main()
