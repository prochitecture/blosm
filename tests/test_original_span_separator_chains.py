"""A chained separator declaration must account for every original route join."""
from copy import deepcopy
import importlib
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace as S
import unittest

package = ModuleType('_separator_chain_tests')
package.__path__ = [str(Path(__file__).resolve().parents[1]/'way/street_network')]
sys.modules[package.__name__] = package
m = importlib.import_module(package.__name__+'.original_span_separator_chains')


def fixture():
    def route(sid, sec, first, last, a, b):
        def end(at, site):
            return S(endpoint=m.SectionEnd(sid, sec, at), junction_id=site, owner_kind='PHYSICAL')
        start, tail = end(0, a), end(1, b)
        return S(kind='CORRIDOR', span=S(street_id=sid, section_ids=(sec,), first_control=first,
            last_control=last, source_control_count=5, in_scope=True, start=start, end=tail),
            start_port=deepcopy(start), end_port=deepcopy(tail))
    routes = [route(1,10,0,2,11,13), route(1,19,2,4,13,12), route(2,20,0,4,12,11)]
    def extent(ordinal, sections, start, end):
        return dict(ordinal=ordinal, sections=[dict(start=dict(section_id=s), end=dict(section_id=s)) for s in sections],
            start=dict(bundle_end=start if isinstance(start,str) else None,
                       junction_id=start if isinstance(start,int) else None),
            end=dict(bundle_end=end if isinstance(end,str) else None,
                     junction_id=end if isinstance(end,int) else None))
    members = [dict(street_id=1, head_order=0, tail_order=0, forward=True, complete=True,
                    spans=[extent(0,[10],'HEAD',13),extent(1,[19],13,'TAIL')]),
               dict(street_id=2, head_order=1, tail_order=1, forward=False, complete=True,
                    spans=[extent(0,[20],'HEAD','TAIL')])]
    bundle = dict(id=7, topology=dict(status='MAPPED',members=members, interior_streets=[],
        junctions=[dict(street_id=1,junction_id=13,through=[[10,1],[19,0]])],required_features=['INTERNAL_JUNCTIONS']))
    snapshot = S(streets={sid:dict(items=[dict(kind='Section',id=sec,category='residential') for sec in sections])
                         for sid,sections in ((1,(10,19)),(2,(20,)))},bundles={7:bundle})
    component = S(ports=S(source_streets=(1,2)),network=S(ordinary_site_ids=(11,12,13),
                  compounds=(),logical_clusters=(),routes=routes))
    return snapshot, component


class ChainedSeparatorTests(unittest.TestCase):
    def setUp(self):
        self.snapshot,self.component=fixture();self.planner=m.ChainedSeparatorPlanner(('residential',))
        self.bundle=self.snapshot.bundles[7]

    def candidate(self):return self.planner.candidate(self.bundle,self.snapshot,self.component)

    def test_complete_original_chains_keep_direction_and_ids(self):
        before=deepcopy((self.snapshot,self.component));actual=self.candidate()
        self.assertEqual(actual.member_span_ids,(((1,10,10),(1,19,19)),((2,20,20),)))
        self.assertEqual(actual.head_ends,(m.SectionEnd(1,10,0),m.SectionEnd(2,20,1)))
        self.assertEqual(actual.tail_ends,(m.SectionEnd(1,19,1),m.SectionEnd(2,20,0)))
        self.assertEqual(actual.intermediate_sites,(((13,),),()))
        self.assertEqual((actual.head_sites,actual.tail_sites),((11,),(12,)))
        self.assertEqual((self.snapshot,self.component),before)

    def test_reversed_multicorridor_member_uses_bundle_order(self):
        for member in self.bundle['topology']['members']:
            member['forward']=not member['forward'];member['spans'].reverse()
            for i,extent in enumerate(member['spans']):
                extent['ordinal']=i;extent['start'],extent['end']=extent['end'],extent['start']
                extent['sections'].reverse()
                for key in ('start','end'):
                    label=extent[key]['bundle_end']
                    if label:extent[key]['bundle_end']='TAIL' if label=='HEAD' else 'HEAD'
        actual=self.candidate()
        self.assertEqual(actual.member_span_ids[0],((1,19,19),(1,10,10)))
        self.assertEqual(actual.head_sites,(12,));self.assertEqual(actual.tail_sites,(11,))

    def test_invalid_route_coverage_is_rejected(self):
        for change in ('gap','overlap','repeated_section','missing_section','inactive','unresolved','wrong_source_count'):
            with self.subTest(change=change):
                self.setUp();first,second=self.component.network.routes[:2]
                if change=='gap':second.span.first_control=3
                if change=='overlap':second.span.first_control=1
                if change=='repeated_section':second.span.section_ids=(10,)
                if change=='missing_section':self.component.network.routes.pop(1)
                if change=='inactive':second.span.in_scope=False
                if change=='unresolved':second.kind='REMOTE_CONNECTION'
                if change=='wrong_source_count':second.span.source_control_count=6
                with self.assertRaises(ValueError):self.candidate()

    def test_every_original_join_and_port_is_required(self):
        for change in ('missing_join','extra_join','duplicate_join','wrong_through','wrong_site','missing_port','wrong_port'):
            with self.subTest(change=change):
                self.setUp();joins=self.bundle['topology']['junctions'];route=self.component.network.routes[0]
                if change=='missing_join':joins.clear()
                if change=='extra_join':joins.append(dict(street_id=2,junction_id=99,through=[]))
                if change=='duplicate_join':joins.append(deepcopy(joins[0]))
                if change=='wrong_through':joins[0]['through'][0][1]=0
                if change=='wrong_site':route.span.end.junction_id=99
                if change=='missing_port':route.end_port=None
                if change=='wrong_port':route.end_port.endpoint=m.SectionEnd(1,99,1)
                with self.assertRaises(ValueError):self.candidate()

    def test_original_extent_mapping_must_agree_with_routes(self):
        for change in ('section','order','join','label','count','partial'):
            with self.subTest(change=change):
                self.setUp();member=self.bundle['topology']['members'][0];extents=member['spans']
                if change=='section':extents[0]['sections'][0]['start']['section_id']=99
                if change=='order':extents[0]['ordinal']=1
                if change=='join':extents[0]['end']['junction_id']=99
                if change=='label':extents[0]['start']['bundle_end']='TAIL'
                if change=='count':extents.pop()
                if change=='partial':member['complete']=False
                with self.assertRaises(ValueError):self.candidate()

    def test_intermediate_owner_cannot_be_a_terminal_owner(self):
        self.component.network.ordinary_site_ids=(12,)
        self.component.network.compounds=(S(site_ids=(11,13)),)
        with self.assertRaisesRegex(ValueError,'multiple separator regions'):self.candidate()

    def test_active_interiors_and_unknown_features_remain_unsupported(self):
        self.bundle['topology']['interior_streets']=[3]
        self.snapshot.streets[3]=dict(items=[dict(kind='Section',id=30,category='residential')])
        with self.assertRaisesRegex(ValueError,'Active Bundle interior'):self.candidate()
        self.snapshot.streets[3]['items'][0]['category']='footway'
        self.assertEqual(self.candidate().bundle_id,7)
        self.bundle['topology']['required_features'].append('UNRECOGNIZED')
        with self.assertRaises(ValueError):self.candidate()

    def test_longitudinal_street_chains_are_not_mistaken_for_one_street_members(self):
        self.bundle['topology']['chains'] = ((1, 3), (2,))
        self.bundle['topology']['members'].append(dict(self.bundle['topology']['members'][0], street_id=3))
        actual = self.planner.prepare(self.snapshot, self.component)
        self.assertFalse(actual.chained or actual.closed or actual.free_ended or actual.selection_ended)
        self.assertEqual(len(actual.pending), 1)
        self.assertIn('multiple original Streets per lateral chain', actual.pending[0]['reason'])

    def test_opt_in_preparation_preserves_basic_dispatch_and_reports_deferred_paths(self):
        self.bundle['topology']['interior_streets']=[3]
        self.snapshot.streets[3]=dict(items=[dict(kind='Section',id=30,category='footway')])
        before=deepcopy((self.snapshot,self.component))
        basic=m.OriginalSpanSeparatorPlanner(('residential',)).prepare(self.snapshot,self.component)
        self.assertEqual(len(basic.pending),1)
        actual=self.planner.prepare(self.snapshot,self.component)
        self.assertEqual(actual.chained,(self.candidate(),))
        self.assertFalse(actual.closed or actual.free_ended or actual.pending)
        self.assertEqual(actual.deferred_interiors[0]['streets'][0]['street_id'],3)
        self.assertEqual((self.snapshot,self.component),before)

    def test_prepare_keeps_existing_simple_closed_branch(self):
        first,second=self.component.network.routes[:2]
        first.span.section_ids=(10,19);first.span.last_control=4
        first.span.end=second.span.end;first.end_port=second.end_port
        self.component.network.routes.pop(1)
        self.bundle['topology']['junctions']=[]
        self.bundle['topology']['required_features']=[]
        self.bundle['topology']['members'][0]['spans']=[{}]
        basic=m.OriginalSpanSeparatorPlanner(('residential',)).prepare(self.snapshot,self.component)
        actual=self.planner.prepare(self.snapshot,self.component)
        self.assertEqual(len(basic.closed),1)
        self.assertEqual(actual.closed,basic.closed)
        self.assertEqual(actual.pending,basic.pending)
        self.assertFalse(actual.chained)

    def test_prepare_retains_explicit_failure_and_ignores_unrelated_component(self):
        self.component.network.routes[1].span.first_control=3
        actual=self.planner.prepare(self.snapshot,self.component)
        self.assertFalse(actual.chained)
        self.assertEqual(actual.pending[0]['bundle_id'],7)
        self.assertEqual(actual.pending[0]['status'],'PENDING')
        self.assertIn('; chained:',actual.pending[0]['reason'])
        self.component.ports.source_streets=(99,)
        unrelated=self.planner.prepare(self.snapshot,self.component)
        self.assertFalse(unrelated.closed or unrelated.free_ended or unrelated.chained or unrelated.pending)


class FreeChainedSeparatorTests(unittest.TestCase):
    def setUp(self):
        self.snapshot,self.component=fixture()
        self.planner=m.ChainedSeparatorPlanner(('residential',))
        self.bundle=self.snapshot.bundles[7]
        self.make_free(head=True)

    def make_free(self, *, head):
        self.snapshot,self.component=fixture();self.bundle=self.snapshot.bundles[7]
        routes=self.component.network.routes
        for route,at in ((routes[0],0),(routes[2],1)) if head else ((routes[1],1),(routes[2],0)):
            end=(route.span.start,route.span.end)[at]
            end.owner_kind='FREE';end.junction_id=None
            setattr(route,'end_port' if at else 'start_port',None)
            route.kind='APPROACH'
        self.component.network.ordinary_site_ids=(12,13) if head else (11,13)

    def candidate(self):
        return self.planner.free_chain_candidate(self.bundle,self.snapshot,self.component)

    def test_free_head_keeps_complete_authored_chains_and_original_ends(self):
        before=deepcopy((self.snapshot,self.component));actual=self.candidate()
        self.assertIsInstance(actual,m.FreeChainedBundleSeparator)
        self.assertFalse(actual.attached_at_head)
        self.assertEqual((actual.head_sites,actual.tail_sites),((),(12,)))
        self.assertEqual(actual.attached_sites,(12,))
        self.assertEqual(actual.free_ends,(m.SectionEnd(1,10,0),m.SectionEnd(2,20,1)))
        self.assertEqual(actual.member_span_ids,(((1,10,10),(1,19,19)),((2,20,20),)))
        self.assertEqual(actual.intermediate_sites,(((13,),),()))
        self.assertEqual((self.snapshot,self.component),before)

    def test_free_tail_is_supported_without_reversing_authored_streets(self):
        self.make_free(head=False);actual=self.candidate()
        self.assertTrue(actual.attached_at_head)
        self.assertEqual((actual.head_sites,actual.tail_sites),((11,),()))
        self.assertEqual(actual.free_ends,(m.SectionEnd(1,19,1),m.SectionEnd(2,20,0)))

    def test_free_chain_dispatch_requires_separate_capability(self):
        previous=self.planner.prepare(self.snapshot,self.component,free_ends=True)
        self.assertEqual(len(previous.pending),1)
        self.assertFalse(hasattr(previous,'free_chained'))
        result=self.planner.prepare(self.snapshot,self.component,free_ends=True,free_chains=True)
        self.assertEqual(result.free_chained,(self.candidate(),))
        self.assertFalse(result.pending or result.chained or result.free_ended or result.closed)
        with self.assertRaisesRegex(ValueError,'free-end capability'):
            self.planner.prepare(self.snapshot,self.component,free_chains=True)

    def test_existing_closed_chain_declarations_are_exact_with_free_capability(self):
        self.snapshot,self.component=fixture()
        previous=self.planner.prepare(self.snapshot,self.component,free_ends=True)
        current=self.planner.prepare(self.snapshot,self.component,free_ends=True,free_chains=True)
        for key in ('closed','free_ended','chained','pending','deferred_interiors'):
            self.assertEqual(getattr(current,key),getattr(previous,key))
        self.assertFalse(current.free_chained)

    def test_free_end_requires_explicit_identity_and_absence_of_port(self):
        for change in ('site','port','endpoint','unresolved','corridor','both_free','interior_free'):
            with self.subTest(change=change):
                self.setUp();first,second=self.component.network.routes[:2]
                if change=='site':first.span.start.junction_id=99
                if change=='port':first.start_port=deepcopy(first.span.start)
                if change=='endpoint':first.span.start.endpoint=m.SectionEnd(1,19,0)
                if change=='unresolved':first.span.start.owner_kind='CATEGORY'
                if change=='corridor':first.kind='CORRIDOR'
                if change=='both_free':
                    second.span.end.owner_kind='FREE';second.span.end.junction_id=None;second.end_port=None
                    second.kind='APPROACH'
                if change=='interior_free':
                    first.span.end.owner_kind='FREE';first.span.end.junction_id=None;first.end_port=None
                with self.assertRaises(ValueError):self.candidate()

    def test_nonterminal_approach_and_route_gaps_are_rejected(self):
        for change in ('extra_approach','gap','missing_route','missing_join','missing_owner','wrong_extent'):
            with self.subTest(change=change):
                self.setUp();second=self.component.network.routes[1]
                if change=='extra_approach':second.kind='APPROACH'
                if change=='gap':second.span.first_control=3
                if change=='missing_route':self.component.network.routes.pop(1)
                if change=='missing_join':self.bundle['topology']['junctions'].clear()
                if change=='missing_owner':self.component.network.ordinary_site_ids=(12,)
                if change=='wrong_extent':self.bundle['topology']['members'][0]['spans'][0]['end']['junction_id']=99
                with self.assertRaises(ValueError):self.candidate()

    def test_members_cannot_put_their_free_ends_on_opposite_bundle_ends(self):
        route=self.component.network.routes[2]
        route.span.start.owner_kind='FREE';route.span.start.junction_id=None;route.start_port=None
        route.span.end.owner_kind='PHYSICAL';route.span.end.junction_id=12
        route.end_port=deepcopy(route.span.end)
        with self.assertRaisesRegex(ValueError,'free-end placement'):self.candidate()

    def test_active_interior_is_not_waived_by_free_chain_branch(self):
        self.bundle['topology']['interior_streets']=[3]
        self.snapshot.streets[3]=dict(items=[dict(kind='Section',id=30,category='residential')])
        result=self.planner.prepare(self.snapshot,self.component,free_ends=True,free_chains=True)
        self.assertFalse(result.free_chained)
        self.assertIn('Active Bundle interior',result.pending[0]['reason'])
        self.snapshot.streets[3]['items'][0]['category']='footway'
        result=self.planner.prepare(self.snapshot,self.component,free_ends=True,free_chains=True)
        self.assertEqual(result.free_chained,(self.candidate(),))
        self.assertEqual(result.deferred_interiors[0]['streets'][0]['street_id'],3)

    def test_shared_intermediate_and_attached_owner_stays_pending(self):
        self.component.network.ordinary_site_ids=()
        self.component.network.compounds=(S(site_ids=(12,13)),)
        with self.assertRaisesRegex(ValueError,'multiple separator regions'):self.candidate()


if __name__=='__main__':unittest.main()
