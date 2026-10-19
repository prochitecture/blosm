"""A crossing continuation must keep complete, unambiguous original ownership."""
from copy import deepcopy
import importlib
from types import SimpleNamespace as S
import unittest

from test_original_span_separator_regions import fixture as region_fixture, m as region_model

m = importlib.import_module(region_model.__package__+'.original_span_separator_continuations')


def fixture():
    snapshot, component = region_fixture()
    rest = component.network.routes[1]
    rest.kind = 'APPROACH'; rest.span.end.owner_kind = 'FREE'
    rest.span.end.junction_id = None; rest.end_port = None
    extra = deepcopy(component.network.routes[-1]); extra.kind = 'APPROACH'
    extra.span.street_id = 4; extra.span.section_ids = (40,)
    extra.span.start = S(endpoint=m.SectionEnd(4,40,0), junction_id=22, owner_kind='PHYSICAL')
    extra.span.end = S(endpoint=m.SectionEnd(4,40,1), junction_id=None, owner_kind='FREE')
    extra.start_port = deepcopy(extra.span.start); extra.end_port = None
    component.network.routes.append(extra)
    snapshot.streets[4] = dict(bundle_id=None, items=[dict(kind='Section',id=40,category='residential')])
    component.ports.source_streets = (1,2,3,4)
    component.ports.sites = (S(junction_id=22, ports=[deepcopy(component.network.routes[2].start_port),
                             deepcopy(component.network.routes[3].end_port), deepcopy(extra.start_port)]),)
    return snapshot, component


class CrossingContinuationTests(unittest.TestCase):
    def setUp(self):
        self.snapshot,self.component = fixture()
        self.planner = m.CrossingFreeRegionPlanner(('residential',))

    def candidate(self):
        return self.planner.continuation_candidate(self.snapshot.bundles[7],self.snapshot,self.component)

    def test_preparation_revalidates_original_prefix_and_reports_remaining_failures(self):
        from dataclasses import replace
        prefix=self.planner.candidate(self.snapshot.bundles[7],self.snapshot,self.component)
        result=self.planner.prepare_continuations(self.snapshot,self.component,prefixes=[prefix])
        self.assertEqual(result.regions,(self.candidate(),));self.assertEqual(result.pending,())
        stale=replace(prefix,crossing_span_id=(99,99,99))
        result=self.planner.prepare_continuations(self.snapshot,self.component,prefixes=[stale])
        self.assertFalse(result.regions);self.assertEqual(result.pending[0]['bundle_id'],7)
        with self.assertRaises(ValueError):
            self.planner.prepare_continuations(self.snapshot,self.component,prefixes=[prefix,prefix])
        self.component.network.routes[-1].span.in_scope=False
        result=self.planner.prepare_continuations(self.snapshot,self.component,prefixes=[prefix])
        self.assertFalse(result.regions);self.assertTrue(result.pending)

    def test_explicit_derived_relationship_does_not_rewrite_original_bundle(self):
        before = deepcopy((self.snapshot,self.component)); row = self.candidate()
        self.assertEqual(row.original_member_street_ids,(1,2))
        self.assertEqual(row.boundary_street_ids,(1,4))
        self.assertEqual(row.boundary_span_ids,((1,19,19),(4,40,40)))
        self.assertEqual(row.attached_sites,((13,),(22,)))
        self.assertEqual(row.attached_ends,(m.SectionEnd(1,19,0),m.SectionEnd(4,40,0)))
        self.assertEqual(row.free_ends,(m.SectionEnd(1,19,1),m.SectionEnd(4,40,1)))
        self.assertEqual((self.snapshot,self.component),before)
        self.assertFalse(hasattr(row,'generated') or hasattr(row,'geometry_accepted'))

    def test_reversed_authored_continuation(self):
        route = self.component.network.routes[-1]
        route.span.start = S(endpoint=m.SectionEnd(4,40,0),junction_id=None,owner_kind='FREE')
        route.span.end = S(endpoint=m.SectionEnd(4,40,1),junction_id=22,owner_kind='PHYSICAL')
        route.start_port = None; route.end_port = deepcopy(route.span.end)
        self.component.ports.sites[0].ports[-1] = deepcopy(route.end_port)
        row = self.candidate()
        self.assertEqual(row.attached_ends[1],m.SectionEnd(4,40,1))
        self.assertEqual(row.free_ends[1],m.SectionEnd(4,40,0))

    def test_original_incidence_is_complete_and_unambiguous(self):
        for case in ('missing','extra','duplicate','wrong_junction','wrong_member','wrong_port'):
            self.setUp(); ports = self.component.ports.sites[0].ports
            if case == 'missing': ports.pop()
            if case == 'extra': ports.append(S(endpoint=m.SectionEnd(8,80,0),junction_id=22))
            if case == 'duplicate': ports.append(deepcopy(ports[-1]))
            if case == 'wrong_junction': ports[-1].junction_id = 99
            if case == 'wrong_member': ports[0].endpoint = m.SectionEnd(2,99,0)
            if case == 'wrong_port': ports[-1].endpoint = m.SectionEnd(4,40,1)
            with self.subTest(case=case),self.assertRaises(ValueError): self.candidate()

    def test_no_invented_free_end_or_omitted_remainder(self):
        for case in ('rest_attached','rest_unresolved','next_attached','next_missing','next_partial','next_inactive'):
            self.setUp(); rest=self.component.network.routes[1]; other=self.component.network.routes[-1]
            if case == 'rest_attached':
                rest.kind='CORRIDOR';rest.span.end.owner_kind='PHYSICAL';rest.span.end.junction_id=12
                rest.end_port=deepcopy(rest.span.end)
            if case == 'rest_unresolved': rest.kind='UNRESOLVED'
            if case == 'next_attached': other.span.end.junction_id=99
            if case == 'next_missing': self.component.network.routes.pop()
            if case == 'next_partial': other.span.last_control-=1
            if case == 'next_inactive': other.span.in_scope=False
            with self.subTest(case=case),self.assertRaises(ValueError): self.candidate()

    def test_competing_bundle_membership_and_categories_reject(self):
        for case in ('owner','member','owned_street','category','absent_source'):
            self.setUp(); source=self.snapshot.streets[4]
            if case == 'owner': source['bundle_id']=88
            if case == 'member': self.snapshot.bundles[88]=dict(topology=dict(members=[dict(street_id=4)]))
            if case == 'owned_street': self.snapshot.bundles[88]=dict(owned_streets=[4])
            if case == 'category': source['items'][0]['category']='footway'
            if case == 'absent_source': self.component.ports.source_streets=(1,2,3)
            with self.subTest(case=case),self.assertRaises(ValueError): self.candidate()

    def test_input_order_and_unrelated_identifiers_do_not_select_the_branch(self):
        expected=self.candidate()
        self.component.network.routes.reverse();self.component.ports.sites[0].ports.reverse()
        self.assertEqual(self.candidate(),expected)
        # The identifiers are deliberately independent of the motivating map.
        self.snapshot.bundles[7000]=self.snapshot.bundles.pop(7);self.snapshot.bundles[7000]['id']=7000
        actual=self.planner.continuation_candidate(self.snapshot.bundles[7000],self.snapshot,self.component)
        self.assertEqual(actual.bundle_id,7000)
        self.assertEqual(actual.boundary_street_ids,(1,4))


if __name__ == '__main__': unittest.main()
