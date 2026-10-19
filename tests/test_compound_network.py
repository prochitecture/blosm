"""Complete component ownership retains finite internal source Streets."""
from copy import deepcopy
from dataclasses import replace
import importlib.util
import math
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module

base=load('compound_base_cases',ROOT/'tests/test_bounded_bundle_network.py')
interior=load('compound_source_records',ROOT/'way/street_network/junction_interiors.py')
compound=load('compound_network_tested',ROOT/'way/street_network/compound_network.py')


def sample(reverse=False,multiple=False):
    s=base.sample(reverse)
    points=((-6.,0.,0.),(0.,0.,0.),(6.,0.,0.))
    if reverse:points=points[::-1]
    sections=(base.inputs.SectionInput(300,points,4.,tags=(('surface','asphalt'),)),)
    if multiple:
        sections=(replace(sections[0],centerline=points[:2]),
                  replace(sections[0],source_id=301,centerline=points[1:],chain_index=1))
    s.streets[30]=dict(id=30,bundle_id=None,endpoints=(100,100),status='SUPPORTED',items=[dict(
        id=ss.source_id,kind='Section',points=ss.centerline,width=ss.width,offset=ss.offset,
        category=ss.category,tags=dict(ss.tags)) for ss in sections])
    j=s.junctions.pop(100)
    link=interior.JunctionInternalLink(100,30,(90,91),tuple(ss.source_id for ss in sections),7,0,
                                     (1,0) if reverse else (0,1))
    packet=interior.PreparedCompoundJunction(j,(base.inputs.StreetInput(30,sections),),(link,),
        (11,90,9,91,10),())
    return s,packet


def branch_sample(reverse=False,multiple=False,*,nested=False):
    """An ordinary complete Street reaches the Bundle through its own side junctions."""
    s,p=sample(reverse,multiple)
    def section(sid,points,width=6.):
        return dict(id=sid,kind='Section',points=tuple((*q,0.) for q in points),
                    width=width,offset=0.,category='residential',tags={})
    def event(jid,point,before,after,street,sid):
        return dict(id=jid,kind='Intersection',location=point,surface_treatment='NATIVE_JUNCTION',
            through=((before,1),(after,0)),attachments=[dict(street_id=street,section_id=sid,
            end=0,connector_id=jid+10000,category='residential',side='LEFT')])
    s.streets[4]['items']=[section(40,[(0,0),(-22,0),(-45,0)]),
        event(900,(-45,0),40,41,80,800),section(41,[(-45,0),(-67,0),(-90,0)])]
    s.streets[80]=dict(id=80,bundle_id=None,endpoints=(900,None),status='SUPPORTED',
        items=[section(800,[(-45,0),(-45,-36),(-45,-80)],4.)])
    if nested:
        s.streets[80]['items']=[section(800,[(-45,0),(-45,-18),(-45,-36)],4.),
            event(910,(-45,-36),800,801,81,810),section(801,[(-45,-36),(-45,-55),(-45,-80)],4.)]
        s.streets[81]=dict(id=81,bundle_id=None,endpoints=(910,None),status='SUPPORTED',
            items=[section(810,[(-45,-36),(-65,-36),(-85,-36)],5.)])
    def input_section(item):
        return base.inputs.SectionInput(item['id'],item['points'],item['width'],item['offset'],item['category'])
    streets=list(p.junction.streets)
    slot=next(i for i,street in enumerate(streets) if street.source_id==4)
    streets[slot]=base.inputs.StreetInput(4,tuple(input_section(i) for i in s.streets[4]['items'] if i['kind']=='Section'))
    p=replace(p,junction=replace(p.junction,streets=tuple(streets)))
    for owner in (4,80):
        for item in s.streets[owner]['items']:
            if item['kind']!='Intersection':continue
            refs=[(owner,sid,end) for sid,end in item['through']]
            refs.extend((a['street_id'],a['section_id'],a['end']) for a in item['attachments'])
            arms=[]
            for street,sid,end in refs:
                sec=next(v for v in s.streets[street]['items'] if v['id']==sid)
                q=sec['points'];a,b=q[:2] if end==0 else q[-1:-3:-1]
                arms.append((math.atan2(b[1]-a[1],b[0]-a[0])%math.tau,street,sec,end))
            arms.sort();local=[];connectors=[]
            for i,(_angle,street,sec,end) in enumerate(arms):
                local.append(base.inputs.StreetInput(street,(input_section(sec),)))
                connectors.append(base.inputs.ConnectorInput(item['id']*10+i,(base.inputs.MemberInput(i,end==0),),i,radius=3.))
            s.junctions[item['id']]=base.inputs.PreparedJunction(item['id'],item['location'],tuple(local),
                tuple(connectors),tuple(base.inputs.ProfileInput(.08) for _ in arms))
    return s,p


class CompoundComponentTests(unittest.TestCase):
    def planner(self):
        return compound.PreparedCompoundComponentPlanner(base.plans.CoupledBundlePlanner())

    def test_both_directions_preserve_snapshot_and_own_every_source_once(self):
        for reverse in (False,True):
            for multiple in (False,True):
                s,p=sample(reverse,multiple);before=deepcopy(s)
                plan=self.planner().plan(s,7,(p,))
                self.assertEqual(s,before);self.assertNotIn(100,s.junctions)
                swept={r.section_id for r in plan.network.sections}
                internal=set(plan.report()['owned_internal_sections'])
                expected={ss['id'] for street in s.streets.values() for ss in street['items'] if ss['kind']=='Section'}
                self.assertFalse(swept&internal);self.assertEqual(swept|internal,expected)
                self.assertEqual(plan.report()['status'],'PLANNED_COMPOUND_COMPONENT')

    def test_input_localization_preserves_original_chain_and_provenance(self):
        s,p=sample(True,True)
        # Force the common frame away from zero to verify one translation.
        origin=(37.,-21.)
        for street in s.streets.values():
            for item in street['items']:
                if item['kind']=='Section':
                    item['points']=tuple((q[0]+origin[0],q[1]+origin[1],q[2]) for q in item['points'])
        def moved(j):
            return replace(j,location=(j.location[0]+origin[0],j.location[1]+origin[1]),
                streets=tuple(replace(street,sections=tuple(replace(ss,centerline=tuple(
                    (q[0]+origin[0],q[1]+origin[1],q[2]) for q in ss.centerline)) for ss in street.sections))
                    for street in j.streets))
        s.junctions={jid:moved(j) for jid,j in s.junctions.items()}
        internal_streets=tuple(replace(street,sections=tuple(replace(ss,centerline=tuple(
            (q[0]+origin[0],q[1]+origin[1],q[2]) for q in ss.centerline)) for ss in street.sections))
            for street in p.internal_streets)
        p=replace(p,junction=moved(p.junction),internal_streets=internal_streets)
        before=deepcopy(s);plan=self.planner().plan(s,7,(p,))
        data=compound.PreparedCompoundComponentInputs(base.bounded_inputs.PreparedBoundedInputs(base.contract)).prepare(s,plan)
        self.assertEqual(s,before);self.assertEqual(tuple(data.network.origin),origin)
        localized=data.junction_interiors[0]
        self.assertIs(localized.junction,next(j for j in data.network.junctions if j.source_id==100))
        for original,new in zip(p.internal_streets[0].sections,localized.internal_streets[0].sections):
            self.assertEqual(new.centerline,tuple((q[0]-origin[0],q[1]-origin[1],q[2]) for q in original.centerline))
            self.assertEqual((new.source_id,new.tags,new.width,new.offset,new.chain_index),
                             (original.source_id,original.tags,original.width,original.offset,original.chain_index))
        self.assertNotIn(300,data.network.sections)

    def test_incomplete_component_is_rejected_instead_of_cutting_a_branch(self):
        s,p=sample();s.streets[4]['endpoints']=(100,900)
        with self.assertRaisesRegex(ValueError,'unowned junction'):
            self.planner().plan(s,7,(p,))


    def test_duplicate_foreign_and_changed_sources_rejected(self):
        s,p=sample()
        changes=[replace(p,original_connector_ids=(11,90,9,90,10)),
                 replace(p,original_connector_ids=(11,90,91,9,999)),
                 replace(p,original_connector_ids=(10,90,9,91,11)),
                 replace(p,links=(replace(p.links[0],junction_id=101),)),
                 replace(p,links=(replace(p.links[0],bundle_side=1),)),
                 replace(p,links=(replace(p.links[0],endpoint_members=(1,0)),)),
                 replace(p,links=(replace(p.links[0],section_ids=(301,)),)),
                 replace(p,internal_streets=(replace(p.internal_streets[0],sections=(
                     replace(p.internal_streets[0].sections[0],width=4.1),)),))]
        for changed in changes:
            with self.subTest(changed=changed):
                with self.assertRaises(ValueError):self.planner().plan(s,7,(changed,))
        with self.assertRaises(ValueError):self.planner().plan(s,7,(p,p))

    def test_internal_section_cannot_also_be_owned_by_a_sweep(self):
        s,p=sample()
        # Duplicated Section identity in another original Street must not pass
        # simply because its input is otherwise accepted by the base planner.
        s.streets[4]['items'][0]['id']=300
        j=p.junction;street=j.streets[2]
        streets=list(j.streets);streets[2]=replace(street,sections=(replace(street.sections[0],source_id=300),))
        p=replace(p,junction=replace(j,streets=tuple(streets)))
        with self.assertRaisesRegex(ValueError,'swept and absorbed'):
            self.planner().plan(s,7,(p,))

class OrdinaryBranchComponentTests(unittest.TestCase):
    def planner(self):
        return compound.PreparedCompoundComponentPlanner(base.plans.CoupledBundlePlanner(
            ('residential',),include_branch_junctions=True))

    def test_complete_recursive_branches_both_directions_and_source_preservation(self):
        for reverse in (False,True):
            for nested in (False,True):
                s,p=branch_sample(reverse,True,nested=nested);before=deepcopy(s)
                plan=self.planner().plan(s,7,(p,))
                self.assertEqual(s,before)
                self.assertEqual(plan.network.ordinary_junctions,(900,910) if nested else (900,))
                expected={i['id'] for street in s.streets.values() for i in street['items'] if i['kind']=='Section'}
                swept={v.section_id for v in plan.network.sections};internal=set(plan.report()['owned_internal_sections'])
                self.assertEqual(swept|internal,expected);self.assertFalse(swept&internal)
                owners={v.section_id:v for v in plan.network.sections}
                self.assertEqual((owners[40].start_junction,owners[40].end_junction),(100,900))
                self.assertEqual(owners[41].start_junction,900)
                if nested:self.assertEqual((owners[800].start_junction,owners[800].end_junction),(900,910))
                data=compound.PreparedCompoundComponentInputs(base.bounded_inputs.PreparedBoundedInputs(base.contract)).prepare(s,plan)
                self.assertEqual({j.source_id for j in data.network.junctions},{8,100,101,900}|({910} if nested else set()))
                for sid,source in data.network.sections.items():
                    original=next(i for street in s.streets.values() for i in street['items'] if i['kind']=='Section' and i['id']==sid)
                    self.assertEqual(source.centerline,original['points'])
                    self.assertEqual(source.width,original['width'])
                self.assertFalse(any(seam.before_section==40 and seam.after_section==41 for seam in data.network.seams))

    def test_previous_policy_remains_explicitly_unchanged(self):
        s,p=branch_sample()
        with self.assertRaisesRegex(ValueError,'requires another junction'):
            compound.PreparedCompoundComponentPlanner(base.plans.CoupledBundlePlanner()).plan(s,7,(p,))

    def test_missing_branch_contract_and_incomplete_incidence_rejected(self):
        for change in ('missing','omit_attachment','duplicate','wrong_through'):
            s,p=branch_sample()
            if change=='missing':s.junctions.pop(900)
            elif change=='omit_attachment':
                j=s.junctions[900];s.junctions[900]=replace(j,connectors=j.connectors[:2])
            elif change=='duplicate':
                j=s.junctions[900];s.junctions[900]=replace(j,connectors=j.connectors+(j.connectors[0],))
            else:s.streets[4]['items'][1]['through']=((40,0),(41,1))
            with self.subTest(change=change):
                with self.assertRaises(ValueError):self.planner().plan(s,7,(p,))

    def test_deferred_feature_is_retained_and_not_inferred_as_a_carriageway(self):
        s,p=branch_sample();item=s.streets[4]['items'][1]
        item['attachments'].append(dict(street_id=999,section_id=9990,end=0,category='cycleway'))
        before=deepcopy(s);plan=self.planner().plan(s,7,(p,))
        self.assertEqual(s,before);self.assertNotIn(999,{v.street_id for v in plan.network.sections})

    def test_external_endpoint_and_foreign_bundle_not_silently_dropped(self):
        for change in ('external','bundle'):
            s,p=branch_sample(nested=True)
            if change=='external':s.streets[81]['endpoints']=(910,999)
            else:s.streets[81]['bundle_id']=99
            with self.subTest(change=change):
                with self.assertRaises(ValueError):self.planner().plan(s,7,(p,))



if __name__=='__main__':unittest.main()
