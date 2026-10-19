"""A lateral Bundle chain can cross Street IDs, never undeclared incidences."""
from copy import deepcopy
import unittest

from test_original_span_separator_chains import fixture as legacy_fixture, m


def fixture():
    snapshot, component = legacy_fixture()
    bundle = snapshot.bundles[7]; topology = bundle['topology']
    # Split the first lateral side at its existing physical junction. Preserve
    # the actual junction and Sections; give each original Street its own range.
    first, second, _ = component.network.routes
    first.span.source_control_count = 3
    second.span.street_id = 3
    second.span.first_control = 0; second.span.last_control = 2
    second.span.source_control_count = 3
    for end in (second.span.start, second.span.end, second.start_port, second.end_port):
        end.endpoint = m.SectionEnd(3, 19, end.endpoint.end)
    snapshot.streets[3] = deepcopy(snapshot.streets[1])
    snapshot.streets[1]['items'] = snapshot.streets[1]['items'][:1]
    snapshot.streets[3]['items'] = snapshot.streets[3]['items'][1:]
    for street in snapshot.streets.values(): street['bundle_id'] = 7
    original = topology['members'][0]
    next_member = deepcopy(original); next_member['street_id'] = 3
    next_member['spans'] = next_member['spans'][1:]; next_member['spans'][0]['ordinal'] = 0
    original['spans'] = original['spans'][:1]
    topology['members'].append(next_member)
    topology['chains'] = ((1, 3), (2,))
    topology['chain_joins'] = [dict(junction_id=13, before=[1, 10, 1], after=[3, 19, 0])]
    topology['junctions'] = []; topology['required_features'] = ['LONGITUDINAL_CHAINS']
    component.ports.source_streets = (1, 2, 3)
    return snapshot, component


class LongitudinalSeparatorTests(unittest.TestCase):
    def candidate(self, snapshot, component, free=False):
        planner = m.ChainedSeparatorPlanner(('residential',), longitudinal_chains=True)
        method = planner.free_chain_candidate if free else planner.candidate
        return method(snapshot.bundles[7], snapshot, component)

    def test_complete_closed_chain_keeps_original_ids_and_sources(self):
        snapshot, component = fixture(); before = deepcopy((snapshot, component))
        item = self.candidate(snapshot, component)
        self.assertEqual(item.member_chains, ((1, 3), (2,)))
        self.assertEqual(item.member_span_ids, (((1, 10, 10), (3, 19, 19)), ((2, 20, 20),)))
        self.assertEqual(item.tail_ends[0], m.SectionEnd(3, 19, 1))
        self.assertEqual(item.intermediate_sites, (((13,),), ()))
        self.assertEqual((snapshot, component), before)
        prepared = m.ChainedSeparatorPlanner(('residential',), longitudinal_chains=True).prepare(snapshot, component)
        self.assertEqual(prepared.chained, (item,)); self.assertFalse(prepared.pending)

    def test_opt_in_does_not_change_legacy_results(self):
        snapshot, component = legacy_fixture()
        self.assertEqual(m.ChainedSeparatorPlanner(('residential',)).prepare(snapshot, component),
            m.ChainedSeparatorPlanner(('residential',), longitudinal_chains=True).prepare(snapshot, component))

    def test_actual_topology_mapper_and_reversed_authored_streets(self):
        from types import SimpleNamespace as S
        from test_bundle_topology import chained, PreparedBundleTopology
        for reverse in (False, True):
            bundle, streets, origins = chained(reverse)
            bundle['topology'] = PreparedBundleTopology(('residential',)).read(bundle, streets, origins=origins).report()
            routes = []
            for sid, head, tail in ((1, 11, 8), (3, 8, 12), (2, 11, 12)):
                def end(at, jid):
                    return S(endpoint=m.SectionEnd(sid, sid*10, at), junction_id=jid, owner_kind='PHYSICAL')
                sites = (tail, head) if reverse else (head, tail)
                start, finish = (end(at, jid) for at, jid in enumerate(sites))
                span = S(street_id=sid, section_ids=(sid*10,), first_control=0, last_control=1,
                         source_control_count=2, in_scope=True, start=start, end=finish)
                routes.append(S(kind='CORRIDOR', span=span, start_port=deepcopy(start), end_port=deepcopy(finish)))
            snapshot = S(streets=streets, bundles={7: bundle})
            component = S(ports=S(source_streets=(1, 2, 3)), network=S(ordinary_site_ids=(11, 8, 12),
                           compounds=(), logical_clusters=(), routes=routes))
            item = self.candidate(snapshot, component)
            self.assertEqual(item.member_chains, ((1, 3), (2,)))
            self.assertEqual(item.head_ends[0], m.SectionEnd(1, 10, int(reverse)))
            self.assertEqual(item.tail_ends[0], m.SectionEnd(3, 30, int(not reverse)))
            self.assertEqual(item.intermediate_sites, (((8,),), ()))

    def test_free_end_at_either_chain_terminal(self):
        for head in (True, False):
            with self.subTest(head=head):
                snapshot, component = fixture(); routes = component.network.routes
                for route, at in ((routes[0], 0), (routes[2], 1)) if head else ((routes[1], 1), (routes[2], 0)):
                    end = (route.span.start, route.span.end)[at]
                    end.owner_kind = 'FREE'; end.junction_id = None; route.kind = 'APPROACH'
                    setattr(route, 'end_port' if at else 'start_port', None)
                component.network.ordinary_site_ids = (12, 13) if head else (11, 13)
                item = self.candidate(snapshot, component, True)
                self.assertEqual(item.attached_at_head, not head)
                self.assertEqual(item.free_ends[0].street_id, 1 if head else 3)
                result = m.ChainedSeparatorPlanner(('residential',), longitudinal_chains=True).prepare(
                    snapshot, component, free_ends=True, free_chains=True)
                self.assertEqual(result.free_chained, (item,)); self.assertFalse(result.pending)

    def test_every_original_identity_and_chain_join_is_required(self):
        for change in ('wrong_join', 'missing_join', 'extra_join', 'duplicate_join', 'wrong_owner',
                       'missing_source', 'wrong_role', 'missing_route', 'gap', 'wrong_extent',
                       'partial', 'shared_chain', 'extra_member', 'unknown_feature'):
            with self.subTest(change=change):
                snapshot, component = fixture(); topology = snapshot.bundles[7]['topology']
                joins = topology['chain_joins']; member = topology['members'][-1]
                if change == 'wrong_join': joins[0]['after'][0] = 1
                if change == 'missing_join': joins.clear()
                if change == 'extra_join': joins.append(dict(junction_id=99, before=[3, 19, 1], after=[2, 20, 0]))
                if change == 'duplicate_join': joins.append(deepcopy(joins[0]))
                if change == 'wrong_owner': snapshot.streets[3]['bundle_id'] = 8
                if change == 'missing_source': component.ports.source_streets = (1, 2)
                if change == 'wrong_role': member['tail_order'] = 1
                if change == 'missing_route': component.network.routes.pop(1)
                if change == 'gap': component.network.routes[1].span.first_control = 1
                if change == 'wrong_extent': member['spans'][0]['start']['bundle_end'] = 'HEAD'
                if change == 'partial': member['complete'] = False
                if change == 'shared_chain': topology['chains'] = ((1, 3), (2, 3))
                if change == 'extra_member': topology['members'].append(deepcopy(member))
                if change == 'unknown_feature': topology['required_features'].append('UNKNOWN')
                with self.assertRaises(ValueError): self.candidate(snapshot, component)


if __name__ == '__main__': unittest.main()
