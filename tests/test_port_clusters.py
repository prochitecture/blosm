"""Logical divided junctions retain original incidences without inferring holes."""
from copy import deepcopy
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest

from test_physical_ports import section, m as ports

spec = importlib.util.spec_from_file_location('tested_port_clusters',
    Path(__file__).resolve().parents[1]/'way/street_network/port_clusters.py')
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)


def fixture():
    locations = ((0., 0.), (20., 0.), (20., 20.), (0., 20.))
    streets, bundles, physical, provenance = {}, {}, [], []

    def source(sid, a, b, start, end, bundle=None):
        streets[sid] = dict(id=sid, status='SUPPORTED', endpoints=(start, end),
            bundle_id=bundle, items=[section(sid*10, (a, b))])

    def record(sid, end):
        return dict(street_id=sid, terminal_id=sid*10, terminal_kind='Section', end=end)

    for i, xy in enumerate(locations):
        nxt = (i+1) % 4
        source(i+1, xy, locations[nxt], i+10, nxt+10, 200+i % 2)
        for k, offset in ((0, (-30., -12.)), (1, (12., -30.))):
            sid = 100+2*i+k
            bid = 300+(i if k == 0 else (i-1) % 4)
            source(sid, xy, (xy[0]+offset[0], xy[1]+offset[1]), i+10, None, bid)
        ends = ((100+2*i, 0), (i+1, 0), ((i-1) % 4+1, 1), (101+2*i, 0))
        physical.append(dict(junction_id=i+10, location=xy, ring=[dict(
            item_kind='Street', connector_id=1000+4*i+j, endpoints=[record(sid, end)])
            for j, (sid, end) in enumerate(ends)]))
    logical_ring = []
    for i in range(4):
        sids = (100+2*i, 101+2*((i+1) % 4))
        members = [dict(street_id=sid, head_order=k, tail_order=k, forward=True, complete=True)
                   for k, sid in enumerate(sids)]
        bundles[300+i] = dict(id=300+i, status='SUPPORTED',
                              topology=dict(status='MAPPED', members=members))
        ends = [record(sid, 0) for sid in sids]
        logical_ring.append(dict(connector_id=2000+i, item_id=300+i, item_kind='Bundle',
                                 leaving=True, endpoints=ends))
        for end in ends:
            sid = end['street_id']; site = (sid-100)//2
            provenance.append(dict(logical_connector_id=2000+i, endpoint=end,
                physical_junction_id=site+10, physical_connector_id=1000+4*site+(0 if sid % 2 == 0 else 3),
                status='EXACT'))
    for i in range(2):
        bundles[200+i] = dict(id=200+i, status='SUPPORTED', source_role='RETAINED_PHYSICAL',
            topology=dict(status='MAPPED', members=[dict(street_id=sid, head_order=k,
                tail_order=k, forward=True, complete=True) for k, sid in enumerate((i+1, i+3))]))
    return NS(streets=streets, bundles=bundles, junction_origins=dict(physical=physical,
        components=[dict(logical_junction_ids=[5000], physical_junction_ids=list(range(10, 14)),
                         shared_physical_junction_ids=[], status='MAPPED', unresolved=[])],
        logical=[dict(junction=dict(junction_id=5000, ring=logical_ring), endpoints=provenance)]))


def retain_internal_connectors(snapshot, street_id=1):
    """Model a real logical self-return: both physical ends remain in its ring."""
    row = snapshot.junction_origins['logical'][0]
    for site in snapshot.junction_origins['physical']:
        for connector in site['ring']:
            endpoint = connector['endpoints'][0]
            if endpoint['street_id'] != street_id:
                continue
            cid = 3000+endpoint['end']
            row['junction']['ring'].append(dict(connector_id=cid, item_kind='Street',
                item_id=street_id, leaving=endpoint['end'] == 0, endpoints=[deepcopy(endpoint)]))
            row['endpoints'].append(dict(logical_connector_id=cid, endpoint=deepcopy(endpoint),
                physical_junction_id=site['junction_id'],
                physical_connector_id=connector['connector_id'], status='EXACT'))


class LogicalPortClusterTests(unittest.TestCase):
    def setUp(self):
        self.snapshot = fixture()
        self.reader = m.PreparedLogicalPortClusters()
        self.plan = self.plan_for(self.snapshot)

    @staticmethod
    def plan_for(snapshot):
        return ports.PreparedPhysicalPorts(('residential',)).prepare(snapshot, physical_ids=range(10, 14))

    def result(self):
        return self.reader.prepare(self.snapshot, self.plan, 5000)

    def test_original_cycle_keeps_members_and_boundary_provenance(self):
        before = deepcopy(self.snapshot)
        result = self.result()
        self.assertEqual(self.snapshot, before)
        self.assertEqual(result.physical_junction_ids, (10, 11, 12, 13))
        self.assertEqual([s.identity for s in result.internal_spans], [(i, i*10, i*10) for i in range(1, 5)])
        self.assertEqual([g.item_id for g in result.groups], [300, 301, 302, 303])
        self.assertEqual([p.endpoint.street_id for p in result.external_ring], [100, 103, 102, 105, 104, 107, 106, 101])
        self.assertEqual([len(b.traversed_spans) for b in result.boundaries], [1, 0]*4)
        for at, boundary in enumerate(result.boundaries):
            self.assertEqual(boundary.start_port, result.external_ring[at])
            self.assertEqual(boundary.end_port, result.external_ring[(at+1) % 8])
            self.assertEqual(len(boundary.original_corners), 1+len(boundary.traversed_spans))
        self.assertEqual(result.source_cycle_rank, 1)
        # The result is a source-topology contract, not a decision about grass.
        self.assertFalse(hasattr(result, 'hole_count'))

    def test_storage_order_and_ring_rotations_do_not_change_contract(self):
        expected = self.result()
        sites = tuple(replace(s, ports=s.ports[1:]+s.ports[:1]) for s in reversed(self.plan.sites))
        self.plan = replace(self.plan, sites=sites, spans=tuple(reversed(self.plan.spans)))
        self.snapshot.junction_origins['components'][0]['physical_junction_ids'].reverse()
        self.assertEqual(self.result(), expected)

    def test_geometry_width_and_object_grouping_do_not_reassign_members(self):
        expected = self.result()
        for sid, street in self.snapshot.streets.items():
            street['object'] = 'All Streets' if sid % 2 else 'Separate Street'
            item = street['items'][0]; a, b = item['points']
            item['points'] = (a, ((a[0]+b[0])*.5+.2, (a[1]+b[1])*.5), b)
            item['width'] += .4
        # Rebuilding private source-control bindings changes only those indices.
        self.plan = self.plan_for(self.snapshot)
        result = self.result()
        self.assertEqual(result.external_ring, expected.external_ring)
        self.assertEqual(result.boundaries, expected.boundaries)
        self.assertEqual(result.groups, expected.groups)
        self.assertEqual([s.identity for s in result.internal_spans], [s.identity for s in expected.internal_spans])

    def test_remote_surface_support_is_distinct_from_mapped_topology(self):
        expected = self.result()
        self.snapshot.bundles[300].update(status='UNSUPPORTED', reason='Remote native patch pending')
        self.assertEqual(self.result(), expected)
        self.assertEqual(self.snapshot.bundles[300]['status'], 'UNSUPPORTED')

    def test_logical_mouth_uses_chain_terminals_not_all_remote_members(self):
        expected = self.result()
        topology = self.snapshot.bundles[300]['topology']
        first, second = topology['members']
        topology['chains'] = ((first['street_id'], 900), (second['street_id'],))
        topology['members'].append(dict(first, street_id=900))
        topology['required_features'] = ('LONGITUDINAL_CHAINS',)
        # A mapped remote continuation is not a third mouth at this junction.
        self.assertEqual(self.result(), expected)
        first['forward'] = False
        with self.assertRaisesRegex(ValueError, 'wrong original Street end'):
            self.result()

    def test_partial_selection_and_missing_incidence_reject(self):
        for mutation in ('site', 'span', 'port'):
            with self.subTest(mutation=mutation):
                plan = self.plan
                if mutation == 'site':
                    plan = replace(plan, sites=plan.sites[1:])
                elif mutation == 'span':
                    plan = replace(plan, spans=plan.spans[1:])
                else:
                    plan = replace(plan, sites=(replace(plan.sites[0], ports=plan.sites[0].ports[1:]),)+plan.sites[1:])
                with self.assertRaises(ValueError): self.reader.prepare(self.snapshot, plan, 5000)

    def test_shared_multiple_or_unmapped_logical_owners_reject(self):
        for changes in (dict(shared_physical_junction_ids=[10]), dict(logical_junction_ids=[5000, 5001]),
                        dict(status='UNRESOLVED'), dict(unresolved=['missing original site'])):
            snapshot = deepcopy(self.snapshot)
            snapshot.junction_origins['components'][0].update(changes)
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.reader.prepare(snapshot, self.plan, 5000)

    def test_duplicate_components_or_original_rings_reject(self):
        for key in ('components', 'physical', 'logical'):
            snapshot = deepcopy(self.snapshot)
            snapshot.junction_origins[key].append(deepcopy(snapshot.junction_origins[key][0]))
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.reader.prepare(snapshot, self.plan, 5000)

    def test_wrong_physical_or_logical_cyclic_order_reject(self):
        site = self.plan.sites[0]
        self.plan = replace(self.plan, sites=(replace(site, ports=tuple(reversed(site.ports))),)+self.plan.sites[1:])
        with self.assertRaisesRegex(ValueError, 'original cyclic order'): self.result()
        self.plan = self.plan_for(self.snapshot)
        self.snapshot.junction_origins['logical'][0]['junction']['ring'].reverse()
        with self.assertRaisesRegex(ValueError, 'group order'): self.result()

    def test_bad_logical_endpoint_provenance_rejects(self):
        for mutation in ('missing', 'duplicate', 'wrong_port', 'wrong_section', 'not_exact'):
            snapshot = deepcopy(self.snapshot)
            rows = snapshot.junction_origins['logical'][0]['endpoints']
            if mutation == 'missing': rows.pop()
            elif mutation == 'duplicate': rows.append(deepcopy(rows[0]))
            elif mutation == 'wrong_port': rows[0]['physical_connector_id'] = rows[1]['physical_connector_id']
            elif mutation == 'wrong_section': rows[0]['endpoint']['terminal_id'] = 999
            else: rows[0]['status'] = 'UNRESOLVED'
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                self.reader.prepare(snapshot, self.plan, 5000)

    def test_incomplete_or_wrong_bundle_roles_reject(self):
        for mutation in ('unmapped', 'missing', 'wrong_street', 'incomplete', 'wrong_owner', 'order'):
            snapshot = deepcopy(self.snapshot)
            topology = snapshot.bundles[300]['topology']; member = topology['members'][0]
            if mutation == 'unmapped': topology['status'] = 'UNRESOLVED'
            elif mutation == 'missing': topology['members'].pop()
            elif mutation == 'wrong_street': member['street_id'] = 999
            elif mutation == 'incomplete': member['complete'] = False
            elif mutation == 'wrong_owner': snapshot.streets[100]['bundle_id'] = 999
            else: member['head_order'] = 1
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                self.reader.prepare(snapshot, self.plan, 5000)

    def test_duplicate_span_or_endpoint_owner_reject(self):
        plans = [replace(self.plan, spans=self.plan.spans+(self.plan.spans[0],))]
        span = self.plan.spans[1]
        plans.append(replace(self.plan, spans=(replace(self.plan.spans[0], end=span.end),)+self.plan.spans[1:]))
        for plan in plans:
            with self.assertRaises(ValueError): self.reader.prepare(self.snapshot, plan, 5000)

    def test_invalid_logical_ids_are_not_coerced(self):
        for value in (True, 5000., -1, '5000'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.reader.prepare(self.snapshot, self.plan, value)

    def test_internal_logical_connectors_keep_both_ends_and_source_span(self):
        expected = self.result()
        retain_internal_connectors(self.snapshot)
        before = deepcopy(self.snapshot)
        with self.assertRaisesRegex(ValueError, 'external port'):
            self.result()  # Existing assets must not silently admit the extension.
        result = self.reader.prepare(self.snapshot, self.plan, 5000, allow_internal_connectors=True)
        self.assertEqual(self.snapshot, before)
        for field in ('external_ring', 'boundaries', 'groups', 'internal_spans', 'source_cycle_rank'):
            self.assertEqual(getattr(result, field), getattr(expected, field))
        self.assertEqual(result.footprint_spans, ((1, 10, 10),))
        self.assertEqual({(g.item_id, g.ports[0].endpoint.end) for g in result.internal_groups}, {(1, 0), (1, 1)})
        self.assertEqual({g.connector_id for g in result.internal_groups}, {3000, 3001})

    def test_incomplete_internal_connector_pair_rejects(self):
        retain_internal_connectors(self.snapshot)
        self.snapshot.junction_origins['logical'][0]['endpoints'].pop()
        with self.assertRaisesRegex(ValueError, 'both original terminal'):
            self.reader.prepare(self.snapshot, self.plan, 5000, allow_internal_connectors=True)

    def test_internal_ring_and_provenance_must_still_agree(self):
        retain_internal_connectors(self.snapshot)
        self.snapshot.junction_origins['logical'][0]['junction']['ring'].pop()
        with self.assertRaisesRegex(ValueError, 'omits an original port'):
            self.reader.prepare(self.snapshot, self.plan, 5000, allow_internal_connectors=True)

    def test_internal_connector_direction_is_not_reinterpreted(self):
        retain_internal_connectors(self.snapshot)
        connector = self.snapshot.junction_origins['logical'][0]['junction']['ring'][-1]
        connector['leaving'] = not connector['leaving']
        with self.assertRaisesRegex(ValueError, 'direction disagrees'):
            self.reader.prepare(self.snapshot, self.plan, 5000, allow_internal_connectors=True)


if __name__ == '__main__': unittest.main()
