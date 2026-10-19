"""Semantic half-profile orientation and explicit span/attachment regressions."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import sys
import unittest

name='tested_bundle_junctions'
spec=importlib.util.spec_from_file_location(name,Path(__file__).resolve().parents[1]/'way/street_network/bundle_junctions.py')
module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)


def sample(forward=True, order=0):
    def span(index, sid):
        return dict(ordinal=index,start=dict(junction_id=8 if index else None),
                    end=dict(junction_id=None if index else 8),
                    sections=[dict(start=dict(section_id=sid),end=dict(section_id=sid))])
    a=dict(street_id=1,forward=forward,head_order=order,tail_order=order,
           spans=[span(0,10 if forward else 11),span(1,11 if forward else 10)])
    b=dict(street_id=2,forward=True,head_order=1-order,tail_order=1-order,spans=[])
    inner_right=(order==0)==forward
    junction=dict(junction_id=8,street_id=1,through=((10,1),(11,0)),attachments=[
        dict(street_id=3,section_id=30,end=0,carriageway=True,side='LEFT' if inner_right else 'RIGHT')])
    return dict(status='MAPPED',bundle_id=7,members=[a,b],junctions=[junction]),junction


def span_chain(forward=True):
    t,j=sample(forward)
    second=deepcopy(j);second.update(junction_id=9,through=((12,1),(13,0)))
    second['attachments'][0].update(street_id=4,section_id=31)
    t['junctions'].append(second)
    chains=([10],[11,12],[13,14]) if forward else ([14,13],[12,11],[10])
    boundaries=(None,8,9,None) if forward else (None,9,8,None)
    member=t['members'][0];member['complete']=True;member['spans']=[]
    for i,ids in enumerate(chains):
        member['spans'].append(dict(ordinal=i,
            start=dict(junction_id=boundaries[i],bundle_end='HEAD' if i==0 else None),
            end=dict(junction_id=boundaries[i+1],bundle_end='TAIL' if i==2 else None),
            sections=[dict(start=dict(section_id=sid),end=dict(section_id=sid)) for sid in ids]))
    return t


class JunctionPlanTests(unittest.TestCase):
    def setUp(self):
        self.planner=module.BundleJunctionPlanner()

    def reject(self, topology, junction, code):
        with self.assertRaises(module.BundleJunctionError) as cm:self.planner.plan(topology,junction)
        self.assertEqual(cm.exception.code,code)

    def test_either_member_and_authored_direction(self):
        for forward in (True,False):
            for order in (0,1):
                with self.subTest(forward=forward,order=order):
                    t,j=sample(forward,order);before=deepcopy((t,j));p=self.planner.plan(t,j)
                    self.assertEqual((t,j),before)
                    self.assertEqual(p.inner_u,float((order==0)==forward))
                    sign=1 if p.inner_u else -1
                    self.assertEqual([(r.section_id,r.end,r.half_profile) for r in p.through],[(10,1,-sign),(11,0,sign)])
                    self.assertEqual((p.before_span,p.after_span),(0,1))
                    self.assertEqual(p.partner_street_id,2)

    def test_street_reversal_flips_local_side_not_physical_crown(self):
        a,j=sample();pa=self.planner.plan(a,j)
        b,k=sample(False);pb=self.planner.plan(b,k)
        self.assertEqual(pa.inner_u,1.-pb.inner_u)
        self.assertEqual(pa.through[0].half_profile,pb.through[1].half_profile)

    def test_inner_attachment_requires_opening_provider(self):
        t,j=sample();j['attachments'][0]['side']='RIGHT';self.reject(t,j,'INNER_ATTACHMENT')

    def test_missing_side_not_inferred(self):
        t,j=sample();j['attachments'][0]['side']=None;self.reject(t,j,'MISSING_SIDE')

    def test_unmapped_extents(self):
        t,j=sample();t['status']='UNRESOLVED';self.reject(t,j,'UNMAPPED_EXTENTS')

    def test_multiple_members_need_another_provider(self):
        t,j=sample();t['members'].append(deepcopy(t['members'][1]));self.reject(t,j,'MEMBER_COUNT')

    def test_lateral_order_must_persist(self):
        t,j=sample();t['members'][0]['tail_order']=1;self.reject(t,j,'MEMBER_ORDER')

    def test_member_reference(self):
        t,j=sample();j['street_id']=3;self.reject(t,j,'MEMBER_REFERENCE')

    def test_junction_at_extent_end_needs_another_patch(self):
        t,j=sample();t['members'][0]['spans'].pop();self.reject(t,j,'SPAN_INCIDENCE')

    def test_distant_spans_cannot_be_joined(self):
        t,j=sample();t['members'][0]['spans'][1]['ordinal']=3;self.reject(t,j,'SPAN_INCIDENCE')

    def test_through_refs_must_match_spans(self):
        t,j=sample();j['through']=((10,1),(12,0));self.reject(t,j,'THROUGH_INCIDENCE')

    def test_same_end_twice_rejected(self):
        t,j=sample();j['through']=((10,0),(11,0));self.reject(t,j,'THROUGH_INCIDENCE')

    def test_opposite_ends_of_the_correct_sections_are_not_incident(self):
        t,j=sample();j['through']=((10,0),(11,1));self.reject(t,j,'THROUGH_INCIDENCE')

    def test_deferred_attachment_does_not_open_median(self):
        t,j=sample();j['attachments'].append(dict(street_id=4,section_id=40,end=1,carriageway=False,side='RIGHT'))
        self.assertEqual(self.planner.plan(t,j).attachments,((3,30,0),))

    def test_only_deferred_attachments(self):
        t,j=sample();j['attachments'][0]['carriageway']=False;self.reject(t,j,'NO_CARRIAGEWAY_ATTACHMENT')

    def test_several_outer_branches_retain_semantic_ids(self):
        t,j=sample();j['attachments'].append(dict(street_id=4,section_id=40,end=1,carriageway=True,side='LEFT'))
        self.assertEqual(self.planner.plan(t,j).attachments,((3,30,0),(4,40,1)))

    def test_report_is_capability_plan_not_generated_surface(self):
        t,j=sample();r=self.planner.review(t)
        self.assertEqual(r[0]['status'],'PLANNED')
        j['attachments'][0]['side']='RIGHT';r=self.planner.review(t)
        self.assertEqual(r[0]['code'],'INNER_ATTACHMENT')

    def test_shared_section_has_two_ends_but_one_record(self):
        t, j = sample()
        second = deepcopy(j); second.update(junction_id=9, through=((11,1),(12,0)))
        second['attachments'][0].update(street_id=4, section_id=31)
        spans = t['members'][0]['spans']; spans[1]['end']['junction_id'] = 9
        spans.append(dict(ordinal=2, start=dict(junction_id=9), end=dict(junction_id=None),
                          sections=[dict(start=dict(section_id=12),end=dict(section_id=12))]))
        t['junctions'].append(second); original=deepcopy(t)
        networks=self.planner.incident_networks(t)
        self.assertEqual(t,original)
        self.assertEqual(len(networks),1)
        network=networks[0]
        self.assertEqual([p.junction_id for p in network.junctions],[8,9])
        self.assertEqual([s.section_id for s in network.sections],[10,11,12,30,31])
        shared=next(s for s in network.sections if s.section_id==11)
        self.assertEqual((shared.start_junction,shared.end_junction),(8,9))
        self.assertTrue(shared.half_profile)
        self.assertEqual(network.report()['status'],'PLANNED_INCIDENT_NETWORK')
        # Input list order is not an identifier or a generation order.
        t['junctions'].reverse()
        self.assertEqual(self.planner.incident_networks(t),networks)

    def test_unbound_end_is_not_reported_as_free(self):
        t,j=sample(); network=self.planner.incident_networks(t)[0]
        before=next(s for s in network.sections if s.section_id==10)
        self.assertEqual((before.start_junction,before.end_junction),(None,8))

    def test_unsupported_patch_remains_in_review(self):
        t,j=sample();j['attachments'][0]['side']='RIGHT'
        self.assertEqual(self.planner.incident_networks(t),())
        self.assertEqual(self.planner.review(t)[0]['code'],'INNER_ATTACHMENT')

    def test_duplicate_patch_is_rejected(self):
        t,j=sample();t['junctions'].append(deepcopy(j))
        with self.assertRaises(module.BundleJunctionError) as cm:self.planner.incident_networks(t)
        self.assertEqual(cm.exception.code,'DUPLICATE_JUNCTION')

    def test_section_cannot_be_both_branch_and_through_member(self):
        t,j=sample();j['attachments'][0]['section_id']=11
        with self.assertRaises(module.BundleJunctionError) as cm:self.planner.incident_networks(t)
        self.assertEqual(cm.exception.code,'SECTION_PROFILE_CONFLICT')

    def test_duplicate_attachment_end_is_rejected(self):
        t,j=sample();j['attachments'].append(deepcopy(j['attachments'][0]))
        with self.assertRaises(module.BundleJunctionError) as cm:self.planner.incident_networks(t)
        self.assertEqual(cm.exception.code,'SECTION_END_OWNER')

    def test_declared_multi_section_span_connects_junctions_once(self):
        t=span_chain();before=deepcopy(t)
        isolated=self.planner.incident_networks(t)
        self.assertEqual(len(isolated),2)
        network,=self.planner.incident_networks(t,include_member_spans=True)
        self.assertEqual(t,before)
        self.assertEqual([s.section_id for s in network.sections],[10,11,12,13,14,30,31])
        self.assertEqual([(s.before_section,s.after_section) for s in network.seams],[(11,12),(13,14)])
        self.assertEqual(network.member_spans,((1,0),(1,1),(1,2)))
        self.assertEqual(network.report()['status'],'PLANNED_MEMBER_NETWORK')

    def test_reverse_bundle_traversal_retains_authored_section_seams(self):
        a,=self.planner.incident_networks(span_chain(),include_member_spans=True)
        b,=self.planner.incident_networks(span_chain(False),include_member_spans=True)
        self.assertEqual(a.seams,b.seams)
        self.assertEqual([s.section_id for s in a.sections],[s.section_id for s in b.sections])
        self.assertEqual(next(s.inner_u for s in a.sections if s.half_profile),
                         1.-next(s.inner_u for s in b.sections if s.half_profile))

    def test_partial_terminal_extent_is_not_promoted_to_complete_section(self):
        t=span_chain();t['members'][0]['complete']=False
        network,=self.planner.incident_networks(t,include_member_spans=True)
        self.assertEqual(network.member_spans,((1,1),))
        self.assertNotIn(14,[s.section_id for s in network.sections])

    def test_unsupported_patch_is_not_bridged_by_a_member_span(self):
        t=span_chain();t['junctions'][1]['attachments'][0]['side']='RIGHT'
        network,=self.planner.incident_networks(t,include_member_spans=True)
        self.assertEqual([j.junction_id for j in network.junctions],[8])
        self.assertEqual(network.member_spans,((1,0),))
        self.assertEqual(network.seams,())

    def test_repeated_span_section_is_an_error(self):
        t=span_chain();t['members'][0]['spans'][1]['sections'].insert(1,
            deepcopy(t['members'][0]['spans'][1]['sections'][0]))
        with self.assertRaises(module.BundleJunctionError) as cm:
            self.planner.incident_networks(t,include_member_spans=True)
        self.assertEqual(cm.exception.code,'SPAN_SECTION_RANGE')

    def test_span_seam_cannot_replace_an_attached_end(self):
        t=span_chain();span=t['members'][0]['spans'][1]
        span['sections'].insert(1,dict(start=dict(section_id=10),end=dict(section_id=10)))
        with self.assertRaises(module.BundleJunctionError) as cm:
            self.planner.incident_networks(t,include_member_spans=True)
        self.assertEqual(cm.exception.code,'SEAM_AT_JUNCTION')


if __name__=='__main__':unittest.main(argv=[__file__])
