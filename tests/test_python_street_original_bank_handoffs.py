"""Exact original-circle exits cannot acquire another nearly coincident bank."""
from copy import copy
from dataclasses import replace
from types import SimpleNamespace
import json
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,Intersection,Network,Key,End
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.routes import Piece,Route
from tested_python_streets.intersections import Approach
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.corner_patch import corner_patch
from tested_python_streets.compound_cycles import _port_footprint
from tested_python_streets.compounds import CompoundGroup,rim_label,SampledFootprint
from tested_python_streets.regions import RegionUnion,Region
from tested_python_streets.updates import BuildContext
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.sampling import AnalyticalContactError
from tested_python_streets.original_bank_handoffs import _ContactSampler,_candidate,sample
from tested_python_streets.source_bank_reconciliation import SourceBankDiscovery,signature,SOURCE_BANK_POLICY


def fixture(*,reverse=False,angle=0.,shift=(0.,0.),mode='STYLE'):
    # A small native three-arm regression, independent of imported feature IDs.
    # The near-straight sector's true circle sagitta is below union welding.
    controls=(((1344.2884521484375,688.14404296875),(1353.258544921875,642.4270629882812),(1353.0902099609375,629.6141357421875)),
              ((1341.866455078125,701.8914794921875),(1342.8089599609375,695.691162109375),(1344.2884521484375,688.14404296875)),
              ((1344.2884521484375,688.14404296875),(1324.2884521484375,688.14404296875)))
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    sections=[];products={};arms=[];grids={}
    for sid,points in enumerate(controls,1):
        xy=np.asarray(points)@rotation.T+shift
        if reverse:xy=xy[::-1]
        section=Section(sid,sid,tuple((x,y,0.) for x,y in xy),width=7.,interpolation=3)
        path,profile=Path(section,.005),Crown(.08,0.)
        products[Key('path',sid)]=path;products[Key('profile',sid)]=profile;sections.append(section)
        forward=(sid!=2)!=reverse
        arms.append(Approach(Route((Piece(path,forward),)),0 if forward else 1,profile))
    order=np.argsort([math.atan2(*arm.frame(0.)[1][::-1])%math.tau for arm in arms])
    corners=tuple(CornerRadiusPolicy.solve(arms[i],arms[j],2.,mode) for i,j in zip(order,np.roll(order,-1)))
    for arm in arms:
        path=products[Key('path',arm.source.id)];cuts=[]
        for corner in corners:
            for original,contact in zip((corner.first,corner.second),corner.stations):
                if original is arm:cuts.append(contact if arm.end==0 else path.length-contact)
        grids[arm.source.id]=np.unique(np.r_[path.samples(0.,path.length,.35,tolerance=.00025),cuts])
    outward=(10.,arms[1].length,arms[2].length)
    ports=tuple(_port_footprint(Key('cycle_port' if i==0 else 'cycle_strip',i+1),arm,cut,8,.7,.00025,grids)
                for i,(arm,cut) in enumerate(zip(arms,outward)))
    patches=[]
    for index,corner in enumerate(corners):
        source=corner_patch(Key('cycle_corner',index+1),corner,.7,.00025)
        source.layout.original_junction=Key('intersection',1);patches.append(source)
    site=Intersection(1,tuple(End(arm.source.id,arm.end) for arm in arms),radius=2.,radius_modes=(mode,)*3)
    network=Network('near-straight-original-handoff',tuple(sections),(site,))
    context=BuildContext(network,products,PatternLibrary())
    physical=SimpleNamespace(approaches=tuple(arms),outward_stations=outward,corners=corners,mouths={})
    sources=(*patches,*ports)
    group=CompoundGroup((site.key,),(),sources,RegionUnion().build(source.region for source in sources),physical)
    return SimpleNamespace(context=context,layouts={}),group,grids


def failure(group):
    sampler=_ContactSampler(tolerance=.005,spacing=.7,max_samples=20000)
    try:sampler.build(group)
    except AnalyticalContactError:return sampler.failed_contact
    raise AssertionError('The native regression did not exercise a failed analytical handoff')


class OriginalBankHandoffTests(unittest.TestCase):
    def test_normal_sampling_repairs_only_owned_provenance_and_retains_strict_circle(self):
        for reverse,mode in ((False,'STYLE'),(True,'STRICT')):
            with self.subTest(reverse=reverse,mode=mode):
                planner,group,grids=fixture(reverse=reverse,mode=mode)
                before=planner.context.network.to_dict()
                points=tuple(source.region.xy.tobytes() for source in group.footprints)
                result=sample(planner,group,grids)
                record,=result.physical_layout.source_bank_handoff_canonicalizations
                self.assertEqual(record['whole_original_arc'],(0.,1.))
                self.assertLessEqual(record['whole_run_reference_distance'],RegionUnion().epsilon)
                self.assertEqual(record['construction']['radius'],2.)
                self.assertEqual(record['construction']['mode'],mode)
                self.assertEqual(record['third_source_candidates'],0)
                self.assertTrue(record['source_exit_exact'] and record['other_end_exact'])
                self.assertEqual(record['future_strip_guard'],'PASS_WITHOUT_MOUTH_ADVANCE')
                self.assertLessEqual(record['fresh_width_proof']['maximum_error'],.002)
                self.assertEqual(SourceBankDiscovery(planner.context,result).discover(),())
                self.assertIs(result.footprints,group.footprints)
                self.assertIs(result.physical_layout.corners,group.physical_layout.corners)
                self.assertIs(result.physical_layout.mouths,group.physical_layout.mouths)
                self.assertEqual(tuple(source.region.xy.tobytes() for source in group.footprints),points)
                self.assertEqual(planner.context.network.to_dict(),before)

    def test_changed_profiles_paths_or_radius_contract_do_not_gain_a_handoff(self):
        planner,group,_=fixture();first,last,labels=failure(group)
        products=dict(planner.context.products)
        products[Key('profile',1)]=Crown(.09,0.)
        changed=BuildContext(planner.context.network,products,planner.context.patterns)
        with self.assertRaisesRegex(ValueError,'profile'):_candidate(changed,group,first,last,labels)
        products=dict(planner.context.products)
        products[Key('path',1)]=Path(planner.context[Key('path',1)].section,.005)
        with self.assertRaisesRegex(ValueError,'Path identity'):
            _candidate(BuildContext(planner.context.network,products,planner.context.patterns),group,first,last,labels)
        site,=planner.context.network.intersections
        network=replace(planner.context.network,intersections=(replace(site,radius=3.),))
        with self.assertRaisesRegex(ValueError,'radius'):
            _candidate(BuildContext(network,planner.context.products,planner.context.patterns),group,first,last,labels)

    def test_a_genuinely_partial_original_arc_is_not_completed_by_proximity(self):
        planner,group,_=fixture();first,last,labels=failure(group)
        clipped=replace(first,last=.9)
        sources={source.region.owner:source for source in group.footprints}
        with self.assertRaisesRegex(ValueError,'complete original circle'):
            _candidate(planner.context,group,clipped,last,(rim_label(clipped,sources),labels[1]))

    def test_third_source_coverage_and_a_run_outside_the_weld_band_are_refused(self):
        planner,group,_=fixture();first,last,labels=failure(group)
        sources={source.region.owner:source for source in group.footprints}
        exit=sources[first.owner].layout.boundaries[labels[0].index](1.)
        xy=exit+np.array(((.0002,0.),(.0004,0.),(.0003,.0002)))
        third=SampledFootprint(Region(Key('unrelated',1),xy),(),None)
        # Choose a third source inside the finite intended exit-edge box.
        direction=sources[first.owner].layout.corners[0].second.frame(0.)[1]
        origin=exit+.01*direction
        xy=origin+np.array(((-.001,-.001),(.001,-.001),(0.,.001)))
        third=replace(third,region=Region(third.region.owner,xy))
        with self.assertRaisesRegex(ValueError,'third source'):
            _candidate(planner.context,replace(group,footprints=(*group.footprints,third)),first,last,labels)
        loops=list(group.loops);loop=loops[0];points=loop.xy.copy()
        start=next(i for i,span in enumerate(loop.spans) if span is last)
        points[start+1]+=np.array((.001,.001))
        loops[0]=replace(loop,xy=points)
        with self.assertRaisesRegex(ValueError,'welding band'):
            _candidate(planner.context,replace(group,loops=tuple(loops)),first,last,labels)

    def test_handoff_only_signature_keeps_local_construction_and_ignores_proof_counters(self):
        planner,group,grids=fixture();result=sample(planner,group,grids)
        self.assertFalse(getattr(result.physical_layout,'source_bank_reconciliation',None))
        value=signature(result)
        self.assertEqual(value[0],SOURCE_BANK_POLICY);self.assertEqual(value[1],())
        before=json.dumps(value,sort_keys=True)
        records=result.physical_layout.source_bank_handoff_canonicalizations
        result.physical_layout.source_bank_handoff_canonicalizations=tuple(dict(record,diagnostic_counter=19) for record in records)
        self.assertEqual(json.dumps(signature(result),sort_keys=True),before)


if __name__=='__main__':unittest.main()
