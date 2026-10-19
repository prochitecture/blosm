"""Connected discovery must retain original ownership, complete rings and gaps."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import sys
import unittest
from test_physical_ports import fixture, mixed_fixture, seam_fixture, section, m as ports
from test_port_network import m as network, compounds

spec = importlib.util.spec_from_file_location('tested_port_connections',
    Path(__file__).resolve().parents[1] / 'way/street_network/port_connections.py')
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)


def add_link(snapshot, sid, start, end, points):
    snapshot.streets[sid] = dict(id=sid, items=[section(sid * 10, points)],
        status='SUPPORTED', endpoints=(start, end), bundle_id=None)
    sites = {s['junction_id']: s for s in snapshot.junction_origins['physical']}
    for at, jid in enumerate((start, end)):
        sites[jid]['ring'].append(dict(item_kind='Street', connector_id=sid * 10 + at,
            endpoints=[dict(street_id=sid, terminal_id=sid * 10, terminal_kind='Section', end=at)]))


class PortConnectionTests(unittest.TestCase):
    def setUp(self):
        self.reader = ports.PreparedPhysicalPorts(('residential',))
        self.discovery = m.PreparedPortConnections(self.reader)

    def test_capable_discovery_crosses_recorded_seam_without_merging_streets(self):
        data=seam_fixture();before=deepcopy(data)
        discovery=m.PreparedPortConnections(ports.PreparedPhysicalPorts(('residential',),seam_connections=True))
        result=discovery.expand(data,physical_ids=(1001,))
        self.assertEqual(result.seam_ids,(2000,));self.assertFalse(result.frontiers)
        self.assertEqual(data,before)
        routes=network.PreparedPortNetwork(compounds.PreparedPortCompound()).prepare(result.plan)
        self.assertEqual({r.span.street_id:r.kind for r in routes.routes if r.span.street_id in (4,8)},
                         {4:'CORRIDOR',8:'APPROACH'})
        self.assertEqual(result,discovery.expand(data,physical_ids=(1001,),seam_ids=()))
        reverse=discovery.expand(data,seam_ids=(2000,))
        self.assertEqual(result.plan,reverse.plan)

    def test_retained_seam_is_explicit_frontier_until_native_connection_exists(self):
        data=fixture();location=(-60.,-20.)
        data.streets[8]=dict(id=8,items=[section(80,(location,(-80.,-20.)))],
                            status='SUPPORTED',endpoints=(None,None),bundle_id=None)
        data.junction_origins['physical'].append(dict(junction_id=2000,location=location,ring=[]))
        data.junction_origins['seams']=[dict(junction_id=2000,location=location,endpoints=[
            dict(street_id=4,terminal_id=40,terminal_kind='Section',end=1),
            dict(street_id=8,terminal_id=80,terminal_kind='Section',end=0)])]
        result=self.discovery.expand(data,physical_ids=(1001,))
        issue,=result.frontiers
        self.assertEqual((issue.owner_kind,issue.junction_id,issue.span_ids),('SEAM',2000,((4,40,40),)))
        self.assertIn('native connection provider',issue.reason)
        routes=network.PreparedPortNetwork(compounds.PreparedPortCompound()).prepare(result.plan)
        route=next(r for r in routes.routes if r.span.street_id==4)
        self.assertEqual(route.kind,'REMOTE_CONNECTION')
        self.assertEqual(route.span.end.owner_kind,'SEAM')
        self.assertEqual(len(routes.unresolved),1)

    def test_one_interior_seed_reaches_both_original_ends_and_preserves_sources(self):
        snapshot = fixture()
        before = deepcopy(snapshot)
        result = self.discovery.expand(snapshot, interior_ids=(21,))
        self.assertEqual(snapshot, before)
        self.assertEqual(result.physical_ids, (1001, 1002))
        self.assertEqual(result.interior_ids, (21, 22))
        expected = self.reader.prepare(snapshot, physical_ids=(1001, 1002), interior_ids=(21, 22))
        self.assertEqual(result.plan, expected)
        self.assertFalse(result.frontiers)
        self.assertEqual(result.components, ((21, 22, 1001, 1002),))
        self.assertEqual(result.cycle_rank, 0)
        routes = network.PreparedPortNetwork(compounds.PreparedPortCompound()).prepare(result.plan)
        self.assertFalse(routes.unresolved)
        self.assertEqual([r.span.identity for r in routes.routes if r.kind == 'CORRIDOR'],
                         [(1, 10, 10), (1, 11, 12), (1, 13, 13)])

    def test_seed_order_does_not_change_output(self):
        a = self.discovery.expand(fixture(), physical_ids=(1002, 1001), interior_ids=(22, 21))
        b = self.discovery.expand(fixture(), physical_ids=(1001, 1002), interior_ids=(21, 22))
        self.assertEqual(a, b)
        self.assertFalse(a.added_sites)

    def test_deferred_physical_incidence_does_not_stop_road_discovery(self):
        source = mixed_fixture()
        before = deepcopy(source)
        result = self.discovery.expand(source, interior_ids=(22,))
        reference = self.discovery.expand(fixture(), interior_ids=(22,))
        self.assertEqual(source, before)
        self.assertEqual(result.physical_ids, reference.physical_ids)
        self.assertEqual(result.interior_ids, reference.interior_ids)
        self.assertEqual(result.components, reference.components)
        self.assertEqual(result.plan.sites, reference.plan.sites)
        self.assertEqual(result.plan.spans, reference.plan.spans)
        self.assertEqual(result.plan.source_streets, reference.plan.source_streets)
        self.assertFalse(result.frontiers)
        self.assertEqual(result.plan.deferred,
            (ports.DeferredPort(1001, ports.SectionEnd(99, 990, 0), 999, 'footway'),))

    def test_failed_ring_is_complete_rejection_with_other_connections_retained(self):
        snapshot = fixture()
        snapshot.streets[7]['status'] = 'UNSUPPORTED'
        result = self.discovery.expand(snapshot, interior_ids=(21,))
        self.assertEqual(result.physical_ids, (1001,))
        self.assertEqual(result.interior_ids, (21, 22))
        issue, = result.frontiers
        self.assertEqual((issue.owner_kind, issue.junction_id, issue.span_ids),
                         ('PHYSICAL', 1002, ((1, 13, 13),)))
        self.assertIn('Street 7', issue.reason)
        self.assertNotIn(6, result.plan.source_streets)
        self.assertNotIn(7, result.plan.source_streets)
        route = next(r for r in network.PreparedPortNetwork(compounds.PreparedPortCompound()).prepare(result.plan).routes
                     if r.span.identity == (1, 13, 13))
        self.assertEqual(route.kind, 'REMOTE_CONNECTION')
        self.assertEqual(route.span.end.junction_id, 1002)

    def test_whole_street_does_not_jump_across_a_rejected_interior_site(self):
        snapshot = fixture()
        snapshot.streets[3]['status'] = 'UNSUPPORTED'
        result = self.discovery.expand(snapshot, physical_ids=(1001,))
        self.assertEqual(result.physical_ids, (1001,))
        self.assertEqual(result.interior_ids, (21,))
        self.assertEqual(result.frontiers[0].junction_id, 22)
        # Section 13 is still present on the whole source, but its far site is
        # not reachable without passing through the failed original site 22.
        self.assertIn((1, 13, 13), [s.identity for s in result.plan.spans])
        self.assertNotIn(1002, [s.junction_id for s in result.added_sites])

    def test_frontier_preserves_the_source_rejection_reason(self):
        snapshot = fixture()
        snapshot.streets[7].update(status='UNSUPPORTED', reason='Grade/elevation treatment is not yet implemented')
        before = deepcopy(snapshot)
        result = self.discovery.expand(snapshot, interior_ids=(21,))
        issue, = result.frontiers
        self.assertEqual(issue.reason, 'Port source Street 7 is unsupported: Grade/elevation treatment is not yet implemented')
        self.assertEqual(snapshot, before)
        self.assertEqual(issue.junction_id, 1002)
        self.assertNotIn(1002, result.physical_ids)

    def test_missing_source_is_distinguished_from_unsupported_source(self):
        snapshot = fixture()
        del snapshot.streets[7]
        result = self.discovery.expand(snapshot, interior_ids=(21,))
        issue, = result.frontiers
        self.assertEqual(issue.reason, 'Port source Street 7 is missing from the original inventory')
        self.assertEqual(issue.span_ids, ((1, 13, 13),))

    def test_missing_original_record_stays_remote_logical(self):
        snapshot = fixture()
        snapshot.junction_origins['physical'].pop()
        result = self.discovery.expand(snapshot, interior_ids=(21,))
        issue, = result.frontiers
        self.assertEqual((issue.owner_kind, issue.junction_id), ('REMOTE_LOGICAL', 1002))
        self.assertIn('No exact', issue.reason)
        self.assertEqual(result.cycle_rank, 0)

    def test_cycle_closes_without_duplicate_source_spans(self):
        snapshot = fixture()
        add_link(snapshot, 8, 1001, 1002, ((-40., 0.), (0., 50.), (80., 0.)))
        result = self.discovery.expand(snapshot, interior_ids=(21,))
        self.assertEqual(result.cycle_rank, 1)
        self.assertEqual(len(result.components), 1)
        self.assertEqual(len(result.plan.spans), 10)
        self.assertEqual(len({s.identity for s in result.plan.spans}), 10)
        self.assertEqual(len(result.added_sites), 3)

    def test_parallel_spans_and_same_site_loop_count_separately(self):
        snapshot = fixture()
        add_link(snapshot, 8, 1001, 1002, ((-40., 0.), (0., 50.), (80., 0.)))
        add_link(snapshot, 9, 1001, 1002, ((-40., 0.), (0., -50.), (80., 0.)))
        add_link(snapshot, 10, 1002, 1002, ((80., 0.), (110., 20.), (110., -20.), (80., 0.)))
        result = self.discovery.expand(snapshot, interior_ids=(21,))
        self.assertEqual(result.cycle_rank, 3)
        self.assertFalse(result.frontiers)

    def test_issue_lists_all_reachable_spans_to_failed_owner(self):
        snapshot = fixture()
        add_link(snapshot, 8, 1001, 1002, ((-40., 0.), (0., 50.), (80., 0.)))
        snapshot.streets[7]['status'] = 'UNSUPPORTED'
        result = self.discovery.expand(snapshot, interior_ids=(21,))
        issue, = result.frontiers
        self.assertEqual(issue.span_ids, ((1, 13, 13), (8, 80, 80)))

    def test_disconnected_explicit_seeds_report_two_components(self):
        snapshot = fixture()
        ring = []
        for sid, p in ((90, (220., 0.)), (91, (190., 20.)), (92, (190., -20.))):
            snapshot.streets[sid] = dict(id=sid, items=[section(sid * 10, ((200., 0.), p))],
                status='SUPPORTED', endpoints=(2001, None), bundle_id=None)
            ring.append(dict(item_kind='Street', connector_id=sid, endpoints=[
                dict(street_id=sid, terminal_id=sid * 10, terminal_kind='Section', end=0)]))
        snapshot.junction_origins['physical'].append(dict(junction_id=2001, location=(200., 0.), ring=ring))
        result = self.discovery.expand(snapshot, physical_ids=(1001, 2001))
        self.assertEqual(result.components, ((21, 22, 1001, 1002), (2001,)))
        self.assertEqual(result.cycle_rank, 0)

    def test_reversal_and_coordinate_edits_preserve_the_connection_graph(self):
        snapshot = fixture()
        before = self.discovery.expand(snapshot, interior_ids=(21,))
        source = snapshot.streets[1]
        source['items'].reverse()
        source['endpoints'] = source['endpoints'][::-1]
        for item in source['items']:
            if item['kind'] == 'Section':
                item['points'] = tuple(reversed(item['points']))
            elif item.get('through'):
                item['through'] = tuple((sid, 1-end) for sid, end in reversed(item['through']))
        for site in snapshot.junction_origins['physical']:
            site['ring'][0]['endpoints'][0]['end'] = 1 - site['ring'][0]['endpoints'][0]['end']
        after = self.discovery.expand(snapshot, interior_ids=(21,))
        self.assertEqual(before.components, after.components)
        self.assertEqual(before.cycle_rank, after.cycle_rank)
        self.assertEqual(before.plan.source_streets, after.plan.source_streets)
        self.assertFalse(after.frontiers)
        source['items'][0]['points'] = ((80., 0.), (60., -1.), (40., 0.))
        self.assertEqual(after, self.discovery.expand(snapshot, interior_ids=(21,)))

    def test_invalid_explicit_seeds_are_not_silently_skipped(self):
        for kwargs in ({}, dict(physical_ids=(999,)), dict(interior_ids=(21, 21)),
                       dict(interior_ids=(True,)), dict(physical_ids=(21,), interior_ids=(21,))):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.discovery.expand(fixture(), **kwargs)

    def test_conflicts_between_locally_valid_rings_reject_the_union(self):
        snapshot = fixture()
        # Street 2 is absent from physical ring 1002's local source set. Its
        # Section ID must nevertheless remain unique in the final component.
        snapshot.streets[6]['items'][0]['id'] = 20
        snapshot.junction_origins['physical'][1]['ring'][1]['endpoints'][0]['terminal_id'] = 20
        with self.assertRaisesRegex(ValueError, 'multiple source-span owners'):
            self.discovery.expand(snapshot, interior_ids=(21,))

    def test_each_discovered_ring_is_prepared_once_before_union_validation(self):
        from unittest.mock import patch
        with patch.object(self.reader, 'prepare', wraps=self.reader.prepare) as prepare:
            self.discovery.expand(fixture(), interior_ids=(21,))
        # Seed, three new rings, final union: traversal does not re-prepare a
        # growing entire component after every newly discovered connection.
        self.assertEqual(prepare.call_count, 5)


if __name__ == '__main__':
    unittest.main()
