"""Complete Bundle planning is semantic and does not construct surfaces."""
from copy import deepcopy
from dataclasses import dataclass, replace
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ROOT = Path(__file__).resolve().parents[1]
plans = load('tested_bounded_plans', ROOT/'way/street_network/bundle_junctions.py')
topology = load('tested_bounded_topology', ROOT/'way/street_network/bundle_topology.py')
inputs = load('tested_bounded_inputs', ROOT/'way/gn/kernel/inputs.py')
bounded_inputs = load('tested_bounded_adapter', ROOT/'way/street_network/bounded_inputs.py')


@dataclass(frozen=True)
class HalfProfile(inputs.ProfileInput):
    half_profile: int = 0


contract = NS(**{key:getattr(inputs,key) for key in ('SectionInput','StreetInput','MemberInput',
    'ConnectorInput','PreparedJunction')},ProfileInput=HalfProfile)


def sample(reverse=False):
    def section(sid, points):
        return dict(id=sid, kind='Section', points=tuple((x,y,0.) for x,y in points),
            width=6., offset=0., category='residential', tags={})
    left = [section(10, [(-6,0),(-6,25),(-6,50)]),
            dict(id=8,kind='Intersection',location=(-6,50),surface_treatment='NATIVE_JUNCTION',
                through=((10,1),(11,0)),attachments=[dict(street_id=3,section_id=30,end=0,
                connector_id=17,category='residential',side='LEFT')]),
            section(11,[(-6,50),(-6,75),(-6,100)])]
    right = [(6,0),(6,50),(6,100)]
    streets = {
        1:dict(id=1,items=left,endpoints=(100,101),bundle_id=7,status='SUPPORTED'),
        2:dict(id=2,items=[section(20,right[::-1] if reverse else right)],
            endpoints=(101,100) if reverse else (100,101),bundle_id=7,status='SUPPORTED'),
        3:dict(id=3,items=[section(30,[(-6,50),(-40,50)])],endpoints=(8,None),bundle_id=None,status='SUPPORTED')}
    for sid, y, x in ((4,0,-40),(5,0,40),(6,100,-40),(7,100,40)):
        streets[sid] = dict(id=sid,items=[section(sid*10,[(0,y),(x,y)])],
            endpoints=(100 if y==0 else 101,None),bundle_id=None,status='SUPPORTED')
    bundle = dict(id=7,head=(1,2),tail=(1,2),head_locations=((-6,0),(6,0)),
        tail_locations=((-6,100),(6,100)),owned_streets=(1,2))
    bundle['topology'] = topology.PreparedBundleTopology(('residential',)).read(bundle,streets).report()
    junctions = {}
    for jid,y in ((100,0),(101,100)):
        local = []
        for sid in (1,2,4 if y==0 else 6,5 if y==0 else 7):
            ss = [s for s in streets[sid]['items'] if s['kind']=='Section']
            local.append(inputs.StreetInput(sid,tuple(inputs.SectionInput(s['id'],s['points'],s['width'],
                chain_index=i) for i,s in enumerate(ss))))
        paired = inputs.ConnectorInput(9,(inputs.MemberInput(0,y==0,0,0,0,0),
            inputs.MemberInput(1,(y==0)!=reverse,1,0,1,1)),
            0,item_kind=inputs.BUNDLE_ITEM,leaving=y==0,bundle_source_id=7)
        west=inputs.ConnectorInput(10,(inputs.MemberInput(2,True),),1)
        east=inputs.ConnectorInput(11,(inputs.MemberInput(3,True),),1)
        connectors = (east,paired,west) if y==0 else (east,west,paired)
        junctions[jid] = inputs.PreparedJunction(jid,(0,y),tuple(local),connectors,
            (inputs.ProfileInput(.08,separator=True),inputs.ProfileInput(.08)))
    local=[]
    for street_id,section_id in ((1,11),(3,30),(1,10)):
        s=next(s for s in streets[street_id]['items'] if s['id']==section_id)
        local.append(inputs.StreetInput(street_id,(inputs.SectionInput(section_id,s['points'],s['width']),)))
    junctions[8]=inputs.PreparedJunction(8,(-6,50),tuple(local),tuple(
        inputs.ConnectorInput(100+i,(inputs.MemberInput(i,i!=2),),i,radius=2.+i)
        for i in range(3)),tuple(inputs.ProfileInput(.08) for i in range(3)))
    return NS(streets=streets,bundles={7:bundle},junctions=junctions,intersections={})


def coupled_sample(reverse=False):
    """Two complete Bundles share the same middle four-arm junction."""
    result=sample();second=sample(reverse)
    sid=lambda value:6 if value==4 else 7 if value==5 else value+100
    section=lambda value:60 if value==40 else 70 if value==50 else value+1000
    jid=lambda value:None if value is None else 101 if value==100 else value+1000
    point=lambda p:(p[0],p[1]+100.,*p[2:])
    new_streets={}
    for old,source in second.streets.items():
        if old in (4,5):continue
        source=deepcopy(source);source['id']=sid(old)
        source['bundle_id']=8 if source['bundle_id']==7 else None
        source['endpoints']=tuple(jid(j) for j in source['endpoints'])
        for item in source['items']:
            if item['kind']=='Section':
                item['id']=section(item['id']);item['points']=tuple(point(p) for p in item['points'])
            else:
                item['id']=jid(item['id']);item['location']=point(item['location'])
                item['through']=tuple((section(s),end) for s,end in item['through'])
                for a in item['attachments']:
                    a.update(street_id=sid(a['street_id']),section_id=section(a['section_id']))
        new_streets[source['id']]=source
    result.streets.update(new_streets)
    bundle=second.bundles[7]
    bundle=dict(bundle,id=8,head=tuple(sid(s) for s in bundle['head']),tail=tuple(sid(s) for s in bundle['tail']),
        head_locations=tuple(point(p) for p in bundle['head_locations']),
        tail_locations=tuple(point(p) for p in bundle['tail_locations']),owned_streets=(101,102))
    bundle['topology']=topology.PreparedBundleTopology(('residential',)).read(bundle,result.streets).report()
    result.bundles[8]=bundle
    for old,j in second.junctions.items():
        streets=tuple(replace(s,source_id=sid(s.source_id),sections=tuple(replace(ss,
            source_id=section(ss.source_id),centerline=tuple(point(p) for p in ss.centerline))
            for ss in s.sections)) for s in j.streets)
        connectors=tuple(replace(c,bundle_source_id=8 if c.bundle_source_id==7 else -1) for c in j.connectors)
        j=replace(j,source_id=jid(old),location=point(j.location),streets=streets,connectors=connectors)
        if old==100:
            previous=result.junctions[101]
            paired=next(c for c in j.connectors if len(c.members)==2)
            paired=replace(paired,members=tuple(replace(m,street_index=m.street_index+4) for m in paired.members))
            east,west,south=previous.connectors
            j=replace(previous,streets=previous.streets+streets[:2],connectors=(east,paired,west,south))
        result.junctions[j.source_id]=j
    return result


def free_sample(reverse=False,sides=(1,),side_junction=True):
    s=sample(reverse)
    for side in sides:
        jid=100+side;s.junctions.pop(jid,None);s.intersections.pop(jid,None)
        for street in s.streets.values():
            street['endpoints']=tuple(None if value==jid else value for value in street['endpoints'])
        for sid in ((4,5) if side==0 else (6,7)):s.streets.pop(sid,None)
    if not side_junction:
        s.junctions.pop(8,None);s.intersections.pop(8,None);s.streets.pop(3,None)
        s.streets[1]['items']=[i for i in s.streets[1]['items'] if i['kind']=='Section']
    bundle=s.bundles[7]
    bundle['topology']=topology.PreparedBundleTopology(('residential',)).read(bundle,s.streets).report()
    return s


def staggered_free_sample(reverse=False):
    """Independent curved controls and two oblique separator caps."""
    s=free_sample(reverse,(0,1),False)
    points=((6.,2.,0.),(7.,30.,0.),(7.,68.,0.),(6.,98.,0.))
    s.streets[2]['items'][0]['points']=points[::-1] if reverse else points
    bundle=s.bundles[7]
    bundle.update(head_locations=((-6.,0.),(6.,2.)),tail_locations=((-6.,100.),(6.,98.)))
    bundle['topology']=topology.PreparedBundleTopology(('residential',)).read(bundle,s.streets).report()
    return s


class FreeEndPlanTests(unittest.TestCase):
    def test_original_transition_kinds_have_separate_provenance_namespaces(self):
        s=free_sample(sides=(0,1),side_junction=False)
        left=s.streets[1]['items'];left[-1]['width']=5.
        left.insert(1,dict(kind='SymLane',id=7,points=[],surface_treatment='WIDTH_TRANSITION',
            transition=dict(before_section=10,after_section=11,source_kind='SymLane',source_id=7,
                            minimum_half_length=1.5,maximum_edge_slope=.3)))
        a=s.streets[2]['items'][0];a['points']=((6.,0.,0.),(6.,25.,0.),(6.,50.,0.))
        b=dict(a,id=21,points=((6.,50.,0.),(6.,75.,0.),(6.,100.,0.)),width=5.)
        crossing=dict(kind='Intersection',id=7,points=[],location=(6.,50.),through=((20,1),(21,0)),
            attachments=[dict(category='footway',street_id=99,section_id=99,connector_id=99,end=0,side='RIGHT')],
            surface_treatment='DEFERRED_ATTACHMENT',transition=dict(before_section=20,after_section=21,
            source_kind='Intersection',source_id=7,minimum_half_length=1.5,maximum_edge_slope=.3))
        s.streets[2]['items']=[a,crossing,b]
        bundle=s.bundles[7];bundle['topology']=topology.PreparedBundleTopology(('residential',)).read(bundle,s.streets).report()
        plan=plans.BoundedBundlePlanner(allow_free_ends=True).plan(s,7)
        data=bounded_inputs.PreparedBoundedInputs(contract).prepare(s,plan)
        self.assertEqual({r.transition_id for r in data.seams},{0,1})
        self.assertEqual({(r.source_kind,r.source_id) for r in data.seams},{('SymLane',7),('Intersection',7)})

    def test_staggered_free_ends_preserve_every_original_control(self):
        for reverse in (False,True):
            s=staggered_free_sample(reverse);before=deepcopy(s)
            plan=plans.CoupledBundlePlanner(allow_free_ends=True).plan(s,7)
            data=bounded_inputs.PreparedBoundedInputs(contract).prepare(s,plan)
            self.assertEqual(s,before)
            for street in s.streets.values():
                for section in street['items']:
                    if section['kind']!='Section':continue
                    expected=tuple((p[0]-data.origin[0],p[1]-data.origin[1],p[2]) for p in section['points'])
                    self.assertEqual(data.sections[section['id']].centerline,expected)

    def test_explicit_head_tail_and_both_free_in_each_authoring_direction(self):
        for reverse in (False,True):
            for sides in ((0,),(1,),(0,1)):
                s=free_sample(reverse,sides);before=deepcopy(s)
                plan=plans.BoundedBundlePlanner(allow_free_ends=True).plan(s,7)
                self.assertEqual(s,before)
                self.assertEqual([e.side for e in plan.free_ends],list(sides))
                for end in plan.free_ends:
                    expected={(10 if end.side==0 else 11,end.side),(20,1-end.side if reverse else end.side)}
                    self.assertEqual(set(end.member_ends),expected)
                data=bounded_inputs.PreparedBoundedInputs(contract).prepare(s,plan)
                self.assertEqual(data.free_ends,plan.free_ends)
                self.assertEqual(set(data.section_profiles),set(data.sections))

    def test_no_junction_bundle_has_complete_sources_profiles_and_seam(self):
        s=free_sample(True,(0,1),False)
        plan=plans.CoupledBundlePlanner(allow_free_ends=True).plan(s,7)
        self.assertFalse(plan.terminal_junctions or plan.junctions)
        self.assertEqual(plan.seams,(plans.IncidentSeam(1,10,11),))
        data=bounded_inputs.PreparedBoundedInputs(contract,profile_rise=.12,datum=2.).prepare(s,plan)
        self.assertFalse(data.junctions)
        self.assertEqual(set(data.sections),{10,11,20})
        self.assertTrue(all(p[2:]==(.12,2.) for p in data.section_profiles.values()))

    def test_missing_or_rejected_terminal_never_becomes_free(self):
        for declared in ({'status':'UNSUPPORTED','order':3},{'status':'ENDPOINT','order':2},{}):
            s=sample();del s.junctions[101];s.intersections[101]=declared
            with self.assertRaises(plans.BundleJunctionError) as cm:
                plans.BoundedBundlePlanner(allow_free_ends=True).plan(s,7)
            self.assertEqual(cm.exception.code,'MISSING_TERMINAL')

    def test_single_declared_dead_end_may_supply_free_member(self):
        s=free_sample(sides=(1,));s.streets[1]['endpoints']=(100,777)
        s.intersections[777]=dict(status='ENDPOINT',order=1)
        self.assertEqual(len(plans.BoundedBundlePlanner(allow_free_ends=True).plan(s,7).free_ends),1)

    def test_one_free_member_does_not_hide_unresolved_opposite_member(self):
        s=free_sample(sides=(1,));s.streets[2]['endpoints']=(100,777)
        with self.assertRaises(plans.BundleJunctionError):
            plans.BoundedBundlePlanner(allow_free_ends=True).plan(s,7)


class BoundedPlanTests(unittest.TestCase):
    def setUp(self):
        self.planner = plans.BoundedBundlePlanner()

    def reject(self, snapshot, code):
        with self.assertRaises(plans.BundleJunctionError) as cm:
            self.planner.plan(snapshot,7)
        self.assertEqual(cm.exception.code,code)

    def test_full_members_and_three_junctions_are_once_owned(self):
        for reverse in (False,True):
            s=sample(reverse);before=deepcopy(s);p=self.planner.plan(s,7)
            self.assertEqual(s,before)
            self.assertEqual(p.terminal_junctions,(100,101))
            self.assertEqual(len(p.sections),8)
            self.assertEqual(p.seams,())
            rows={r.section_id:r for r in p.sections}
            self.assertEqual((rows[10].start_junction,rows[10].end_junction),(100,8))
            self.assertEqual((rows[11].start_junction,rows[11].end_junction),(8,101))
            self.assertEqual((rows[20].start_junction,rows[20].end_junction),
                (101,100) if reverse else (100,101))
            self.assertEqual(rows[20].inner_u,float(reverse))
            self.assertEqual(p.report()['status'],'PLANNED_BOUNDED_NETWORK')

    def test_input_table_orders_are_not_identity(self):
        s=sample();expected=self.planner.plan(s,7)
        s.junctions=dict(reversed(list(s.junctions.items())))
        s.bundles[7]['topology']['members']=s.bundles[7]['topology']['members'][::-1]
        s.streets=dict(reversed(list(s.streets.items())))
        self.assertEqual(self.planner.plan(s,7),expected)

    def test_partial_member(self):
        s=sample();s.bundles[7]['topology']['members'][0]['complete']=False
        self.reject(s,'PARTIAL_EXTENTS')

    def test_interior_member_requires_opening(self):
        s=sample();s.bundles[7]['topology']['interior_streets']=[42]
        self.reject(s,'INTERIOR_STREETS')

    def test_explicit_feature_selection_retains_deferred_interior_streets(self):
        s=sample();s.bundles[7]['topology']['interior_streets']=[42]
        s.streets[42]=dict(id=42,items=[dict(kind='Section',id=420,category='footway')])
        plan=plans.BoundedBundlePlanner(('residential',)).plan(s,7)
        self.assertEqual(plan.deferred_streets,(42,))
        self.assertNotIn(420,[r.section_id for r in plan.sections])
        s.streets[42]['items'][0]['category']='residential'
        with self.assertRaises(plans.BundleJunctionError) as cm:
            plans.BoundedBundlePlanner(('residential',)).plan(s,7)
        self.assertEqual(cm.exception.code,'INTERIOR_STREETS')

    def test_inner_branch_requires_opening(self):
        s=sample();s.bundles[7]['topology']['junctions'][0]['attachments'][0]['side']='RIGHT'
        self.reject(s,'INNER_ATTACHMENT')

    def test_incomplete_member_span(self):
        s=sample();member=s.bundles[7]['topology']['members'][0];member['spans']=member['spans'][:-1]
        self.reject(s,'SPAN_INCIDENCE')

    def test_missing_terminal(self):
        s=sample();del s.junctions[101];self.reject(s,'MISSING_TERMINAL')

    def test_duplicate_terminal_owner(self):
        s=sample();s.junctions[102]=replace(s.junctions[100],source_id=102)
        self.reject(s,'TERMINAL_OWNER')

    def test_mixed_terminal_ends(self):
        s=sample();j=s.junctions[100];c=j.connectors[1]
        c=replace(c,members=(c.members[0],replace(c.members[1],leaving=False)))
        s.junctions[100]=replace(j,connectors=(j.connectors[0],c,j.connectors[2]))
        self.reject(s,'TERMINAL_INCIDENCE')

    def test_wrong_terminal_section(self):
        s=sample();j=s.junctions[100];street=j.streets[0]
        s.junctions[100]=replace(j,streets=(replace(street,sections=street.sections[::-1]),*j.streets[1:]))
        self.reject(s,'TERMINAL_SECTION')

    def test_external_junction_is_not_treated_as_free(self):
        s=sample();s.streets[3]['endpoints']=(8,777)
        s.intersections[777]=dict(status='SUPPORTED',order=3)
        self.reject(s,'EXTERNAL_JUNCTION')

    def test_declared_free_endpoint_is_allowed(self):
        s=sample();s.streets[3]['endpoints']=(8,777)
        s.intersections[777]=dict(status='ENDPOINT',order=1)
        self.assertEqual(len(self.planner.plan(s,7).sections),8)

    def test_multisection_branch_preserves_seam(self):
        s=sample();first=s.streets[3]['items'][0]
        second=dict(first,id=31,points=((-40,50,0),(-60,50,0)))
        s.streets[3]['items'].append(second)
        p=self.planner.plan(s,7)
        self.assertEqual(p.seams,(plans.IncidentSeam(3,30,31),))
        self.assertEqual(len(p.sections),9)

    def test_another_bundle_in_branch(self):
        s=sample();s.streets[3]['bundle_id']=99;self.reject(s,'OTHER_BUNDLE')

    def test_unsupported_branch_is_not_omitted(self):
        s=sample();s.streets[3]['status']='UNSUPPORTED';self.reject(s,'STREET_SCOPE')

    def test_failed_plan_reports_reason_without_generated_claim(self):
        s=sample();s.junctions.clear();r=self.planner.review(s,7)
        self.assertEqual(r['status'],'UNSUPPORTED');self.assertEqual(r['code'],'MISSING_TERMINAL')

    def test_prepared_contract_preserves_sources_radii_and_authored_roles(self):
        for reverse in (False,True):
            s=sample(reverse);saved=deepcopy(s);plan=self.planner.plan(s,7)
            data=bounded_inputs.PreparedBoundedInputs(contract).prepare(s,plan)
            self.assertEqual(s,saved)
            self.assertEqual(len(data.sections),8)
            for j in data.junctions:
                self.assertEqual([c.source_id for c in j.connectors],[c.source_id for c in s.junctions[j.source_id].connectors])
                self.assertEqual([c.radius for c in j.connectors],[c.radius for c in s.junctions[j.source_id].connectors])
                self.assertTrue(all(len(x.sections)==1 and len(x.sections[0].centerline)==2 for x in j.streets))
            side=next(j for j in data.junctions if j.source_id==8)
            self.assertEqual([p.half_profile for p in side.profiles],[1,0,-1])
            self.assertEqual(data.section_profiles[20],(True,float(reverse),.08,0.))
            self.assertEqual(data.sections[20].centerline,s.streets[2]['items'][0]['points'])

    def test_local_frame_and_provenance_preserve_every_control(self):
        s=sample();dx,dy=200000.,500000.
        for street in s.streets.values():
            for item in street['items']:
                if item['kind']=='Section':
                    item['points']=tuple((p[0]+dx,p[1]+dy,p[2]) for p in item['points'])
                    item['tags']={'surface':'asphalt'}
        # Plan identity is independent of the translated coordinates. Prepared
        # local terminal controls are read from the original source Sections.
        s.junctions={jid:replace(j,location=(j.location[0]+dx,j.location[1]+dy)) for jid,j in s.junctions.items()}
        data=bounded_inputs.PreparedBoundedInputs(contract).prepare(s,self.planner.plan(s,7))
        self.assertEqual(data.origin,(dx,dy))
        self.assertEqual(data.sections[10].centerline,((-6,0,0),(-6,25,0),(-6,50,0)))
        self.assertEqual(data.sections[10].tags,(('surface','asphalt'),))

    def test_inconsistent_terminal_profile_is_not_silently_changed(self):
        s=sample();j=s.junctions[101]
        s.junctions[101]=replace(j,profiles=(replace(j.profiles[0],rise=.09),j.profiles[1]))
        with self.assertRaisesRegex(ValueError,'one profile rise'):
            bounded_inputs.PreparedBoundedInputs(contract).prepare(s,self.planner.plan(s,7))


class CoupledPlanTests(unittest.TestCase):
    def test_shared_physical_site_connects_unqualified_logical_terminals(self):
        s=coupled_sample();s.junctions.clear()
        s.junction_origins=dict(logical=[dict(junction=dict(junction_id=90,
            ring=[dict(item_id=7,item_kind='Bundle')])),dict(junction=dict(junction_id=91,
            ring=[dict(item_id=8,item_kind='Bundle')]))],components=[dict(
                logical_junction_ids=(90,91),physical_junction_ids=(80,81,82),
                shared_physical_junction_ids=(81,))])
        before=deepcopy(s);planner=plans.CoupledBundlePlanner()
        self.assertEqual(planner.component(s,7),(7,8))
        self.assertEqual(planner.component(s,8),(7,8))
        row,=planner.review(s)
        self.assertEqual(row['code'],'SHARED_PHYSICAL_JUNCTION')
        self.assertEqual(row['bundle_ids'],(7,8))
        self.assertEqual(s,before)

    def test_origin_bundle_absent_from_inventory_is_not_omitted(self):
        s=sample();s.junction_origins=dict(logical=[dict(junction=dict(junction_id=90,
            ring=[dict(item_id=999,item_kind='Bundle')]))])
        with self.assertRaises(plans.BundleJunctionError) as cm:
            plans.CoupledBundlePlanner.component(s,7)
        self.assertEqual(cm.exception.code,'MISSING_BUNDLE')

    def test_both_seeds_produce_one_once_owned_component(self):
        for reverse in (False,True):
            s=coupled_sample(reverse);saved=deepcopy(s);planner=plans.CoupledBundlePlanner()
            a,b=planner.plan(s,7),planner.plan(s,8)
            self.assertEqual(a,b);self.assertEqual(s,saved)
            self.assertEqual(a.bundle_ids,(7,8));self.assertEqual(len(a.sections),14)
            self.assertEqual(a.terminal_junctions,(100,101,1101))
            self.assertEqual(len(a.junctions),2)
            self.assertEqual(dict(a.section_bundles),{10:7,11:7,20:7,1010:8,1011:8,1020:8})
            review=planner.review(s)
            self.assertEqual(len(review),1);self.assertEqual(review[0]['status'],'PLANNED_COUPLED_NETWORK')

    def test_shared_junction_is_prepared_once_with_all_six_members(self):
        s=coupled_sample();plan=plans.CoupledBundlePlanner().plan(s,7)
        data=bounded_inputs.PreparedBoundedInputs(contract).prepare(s,plan)
        self.assertEqual(len(data.sections),14);self.assertEqual(len(data.junctions),5)
        shared=next(j for j in data.junctions if j.source_id==101)
        self.assertEqual(len(shared.connectors),4);self.assertEqual(len(shared.streets),6)
        self.assertEqual({c.bundle_source_id for c in shared.connectors if len(c.members)==2},{7,8})
        self.assertEqual(data.section_bundles,dict(plan.section_bundles))

    def test_rejected_neighbor_is_not_omitted_from_component(self):
        s=coupled_sample();s.bundles[8]['topology']['members'][0]['complete']=False
        rows=plans.CoupledBundlePlanner().review(s)
        self.assertEqual(len(rows),1);self.assertEqual(rows[0]['bundle_ids'],(7,8))
        self.assertEqual(rows[0]['code'],'PARTIAL_EXTENTS')

    def test_missing_terminal_remains_explicit(self):
        s=coupled_sample();del s.junctions[1101]
        rows=plans.CoupledBundlePlanner().review(s)
        self.assertEqual(len(rows),1);self.assertEqual(rows[0]['code'],'MISSING_TERMINAL')

    def test_table_permutation_keeps_component_identity(self):
        s=coupled_sample();planner=plans.CoupledBundlePlanner();expected=planner.plan(s,7)
        s.bundles=dict(reversed(list(s.bundles.items())))
        s.junctions=dict(reversed(list(s.junctions.items())))
        self.assertEqual(planner.plan(s,8),expected)


if __name__=='__main__':unittest.main(argv=[__file__])
