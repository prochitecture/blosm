"""Crossing regions must preserve the entire original inventory and its gaps."""
from copy import deepcopy
import importlib
from types import SimpleNamespace as S
import unittest

from test_original_span_separator_chains import fixture as chain_fixture, m as chain_model

m=importlib.import_module(chain_model.__package__+'.original_span_separator_regions')


def fixture():
    snapshot,component=chain_fixture()
    second=component.network.routes[2]
    second.span.start.junction_id=22;second.start_port.junction_id=22
    cross=deepcopy(component.network.routes[0])
    cross.span.street_id=3;cross.span.section_ids=(30,);cross.span.source_control_count=3
    for at,site in enumerate((13,22)):
        end=S(endpoint=m.SectionEnd(3,30,at),junction_id=site,owner_kind='PHYSICAL')
        if at==0:cross.span.start=end;cross.start_port=deepcopy(end)
        else:cross.span.end=end;cross.end_port=deepcopy(end)
    component.network.routes.append(cross)
    component.network.ordinary_site_ids=(11,12,13,22)
    component.ports.source_streets=(1,2,3)
    snapshot.streets[3]=dict(items=[dict(kind='Section',id=30,category='residential')])
    bundle=snapshot.bundles[7]
    bundle['topology']['interior_streets']=[3]
    bundle['topology']['required_features'].append('INTERIOR_STREETS')
    return snapshot,component


class CrossingRegionTests(unittest.TestCase):
    def setUp(self):
        self.snapshot,self.component=fixture()
        self.planner=m.CrossingSeparatorRegionPlanner(('residential',))
        self.bundle=self.snapshot.bundles[7]

    def candidate(self):return self.planner.candidate(self.bundle,self.snapshot,self.component)

    def test_exact_original_prefix_and_remaining_extent_are_both_reported(self):
        before=deepcopy((self.snapshot,self.component));row=self.candidate()
        self.assertEqual(row.member_span_ids,(((1,10,10),),((2,20,20),)))
        self.assertEqual(row.remaining_member_span_ids,(((1,19,19),),()))
        self.assertEqual(row.head_ends,(m.SectionEnd(1,10,0),m.SectionEnd(2,20,1)))
        self.assertEqual(row.tail_ends,(m.SectionEnd(1,10,1),m.SectionEnd(2,20,0)))
        self.assertEqual(row.head_sites,(11,));self.assertEqual(row.tail_sites,((13,),(22,)))
        self.assertEqual(row.crossing_span_id,(3,30,30))
        self.assertEqual(row.crossing_ends,(m.SectionEnd(3,30,0),m.SectionEnd(3,30,1)))
        self.assertEqual((self.snapshot,self.component),before)

    def test_common_attached_cap_can_be_at_either_bundle_end(self):
        expected=self.candidate()
        for member in self.bundle['topology']['members']:
            member['forward']=not member['forward'];member['spans'].reverse()
            for i,extent in enumerate(member['spans']):
                extent['ordinal']=i;extent['start'],extent['end']=extent['end'],extent['start']
                extent['sections'].reverse()
                for key in ('start','end'):
                    label=extent[key]['bundle_end']
                    if label:extent[key]['bundle_end']='TAIL' if label=='HEAD' else 'HEAD'
        self.assertEqual(self.candidate(),expected)

    def test_leftover_member_inventory_is_validated_not_discarded(self):
        for case in ('missing','gap','wrong_source_count','wrong_section','inactive','invalid_free_end'):
            with self.subTest(case=case):
                self.setUp();rest=self.component.network.routes[1]
                if case=='missing':self.component.network.routes.pop(1)
                if case=='gap':rest.span.first_control+=1
                if case=='wrong_source_count':rest.span.source_control_count+=1
                if case=='wrong_section':rest.span.section_ids=(88,)
                if case=='inactive':rest.span.in_scope=False
                if case=='invalid_free_end':rest.span.start.owner_kind='FREE'
                with self.assertRaises(ValueError):self.candidate()

    def test_crossing_requires_complete_corridor_and_both_original_endpoints(self):
        for case in ('missing','split','not_corridor','incomplete','inactive','wrong_section','missing_port','wrong_junction','duplicate_end'):
            with self.subTest(case=case):
                self.setUp();cross=self.component.network.routes[-1]
                if case=='missing':self.component.network.routes.pop()
                if case=='split':self.component.network.routes.append(deepcopy(cross))
                if case=='not_corridor':cross.kind='APPROACH'
                if case=='incomplete':cross.span.last_control-=1
                if case=='inactive':cross.span.in_scope=False
                if case=='wrong_section':cross.span.section_ids=(99,)
                if case=='missing_port':cross.end_port=None
                if case=='wrong_junction':cross.span.end.junction_id=999;cross.end_port.junction_id=999
                if case=='duplicate_end':cross.span.end.junction_id=13;cross.end_port.junction_id=13
                with self.assertRaises(ValueError):self.candidate()

    def test_scope_never_drops_an_active_interior(self):
        self.bundle['topology']['interior_streets'].append(4)
        self.snapshot.streets[4]=dict(items=[dict(kind='Section',id=40,category='footway')])
        row=self.candidate();self.assertEqual(row.deferred_interiors[0]['street_id'],4)
        self.snapshot.streets[4]['items'][0]['category']='residential'
        with self.assertRaisesRegex(ValueError,'one complete active'):self.candidate()
        self.snapshot.streets[4]['items'][0]['category']='footway'
        fake=deepcopy(self.component.network.routes[-1]);fake.span.street_id=4
        self.component.network.routes.append(fake)
        with self.assertRaisesRegex(ValueError,'one complete active'):self.candidate()

    def test_declared_extents_and_throughs_cannot_go_stale(self):
        for case in ('missing_join','wrong_through','wrong_extent','duplicate_join','unknown_feature'):
            with self.subTest(case=case):
                self.setUp();t=self.bundle['topology']
                if case=='missing_join':t['junctions'].clear()
                if case=='wrong_through':t['junctions'][0]['through'][0][0]=99
                if case=='wrong_extent':t['members'][0]['spans'][0]['end']['junction_id']=99
                if case=='duplicate_join':t['junctions'].append(deepcopy(t['junctions'][0]))
                if case=='unknown_feature':t['required_features'].append('UNKNOWN')
                with self.assertRaises(ValueError):self.candidate()

    def test_prepare_returns_regions_without_claiming_complete_bundle_coverage(self):
        result=self.planner.prepare(self.snapshot,self.component,bundle_ids=(7,))
        self.assertEqual(result.regions,(self.candidate(),));self.assertFalse(result.pending)
        self.assertFalse(hasattr(result,'chained') or hasattr(result,'closed'))
        self.component.network.routes[-1].kind='UNRESOLVED'
        rejected=self.planner.prepare(self.snapshot,self.component,bundle_ids=(7,))
        self.assertFalse(rejected.regions);self.assertEqual(rejected.pending[0]['bundle_id'],7)


if __name__=='__main__':unittest.main()
