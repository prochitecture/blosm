"""Finite paired terminal cuts preserve actual strips and source aliases."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.updates import BuildContext
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.bundle_visible_runs import BundleRun,MemberRun
from tested_python_streets.finite_terminal_entrances import FinitePairedEntrance
from tested_python_streets.cycle_terminal_ports import (TerminalPair,footprint,minimum_station,
    certify,discover,prepare_requested,original_section_gap)
from tested_python_streets.compounds import CompoundGroup,rim_label
from tested_python_streets.regions import RegionUnion,BoundaryLoop


def fixture(mode=3,reverse=False,reflection=False,angle=0.,translation=(0.,0.),separation=8.,divergence=0.):
    rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
    matrix=rotation@np.diag((1.,-1. if reflection else 1.))
    source=[];visits=[]
    for sid,y in ((1,.5*separation),(2,-.5*separation)):
        points=np.array(((0.,y),(30.,y+.5),(60.,y+1.),(100.,y+1.)))
        points[:,1]+=np.sign(y)*divergence*points[:,0]
        points=points@matrix.T+translation
        if reverse:points=points[::-1]
        source.append(model.Section(sid,sid,tuple((x,y,0.) for x,y in points),width=10.8,interpolation=mode))
        visits.append(model.Traversal(sid,not reverse))
    if reflection:visits=visits[::-1]
    bundle=model.Bundle(10,tuple(model.BundleMember((v,)) for v in visits))
    sides=[]
    for i,v in enumerate(visits):
        p=source[v.section_id-1]
        root=p.points[-1 if reverse else 0]
        remote=tuple(np.asarray(root)+np.array((0.,-30.,0.)))
        sides.append(model.Section(100+i,100+i,(root,remote)))
    sites=tuple(model.Intersection(200+i,(model.End(v.section_id,int(reverse)),model.End(100+i,0))) for i,v in enumerate(visits))
    n=model.Network('finite-terminal-test',tuple(source+sides),intersections=sites,bundles=(bundle,))
    products={model.Key('path',s.id):Path(s) for s in n.sections}
    for role,v in enumerate(visits):
        products[model.Key('profile',v.section_id)]=Crown(.08,0.,(1 if role==0 else -1)*(1 if v.forward else -1))
    ctx=BuildContext(n,products,PatternLibrary())
    tokens=(('bundle_end',0),('bundle_end',1))
    members=tuple(MemberRun((v,),(model.End(v.section_id,int(not v.forward)),model.End(v.section_id,int(v.forward))),tokens) for v in visits)
    run=BundleRun(members,tokens)
    originals=tuple(Approach(products[model.Key('path',v.section_id)],int(reverse),products[model.Key('profile',v.section_id)]) for v in visits)
    return ctx,TerminalPair(bundle.id,run,0,originals),sites


class TerminalPortsTests(unittest.TestCase):
    def test_required_fibre_can_keep_two_real_profiles_after_original_strips_diverge(self):
        for mode in (0,3):
            for reverse,reflection,angle,translation in ((False,False,0.,(0.,0.)),
                    (True,True,.713,(1800.,-1200.))):
                with self.subTest(mode=mode,reverse=reverse):
                    c,pair,_=fixture(mode,reverse,reflection,angle,translation,divergence=.1)
                    before=c.network.to_dict()
                    a=pair.build(c,{})
                    self.assertLess(a.gap(minimum_station(a)),0.)
                    stations={model.End(original.source.id,original.end):22. for original in pair.originals}
                    chosen,at,gap=prepare_requested(pair,c,{},stations)
                    self.assertIsNone(chosen)
                    self.assertGreater(gap.gap,gap.allowance+gap.margin)
                    self.assertEqual({row['section'] for row in gap.sources},{1,2})
                    for row in gap.sources:
                        path=c[model.Key('path',row['section'])]
                        for curb in row['curbs']:
                            np.testing.assert_array_equal(path.point(curb['original_station'],curb['offset']),curb['point'])
                            self.assertLess(curb['inversion_residual'],gap.margin)
                    self.assertEqual(c.network.to_dict(),before)
                    self.assertTrue(gap.report()['original_separate_profiles_retained'])

    def test_merged_success_stays_exact_and_sampled_field_cannot_certify_separation(self):
        c,pair,_=fixture()
        stations={model.End(original.source.id,original.end):6. for original in pair.originals}
        a,at,gap=prepare_requested(pair,c,{},stations)
        self.assertIsNone(gap)
        old=pair.build(c,{})
        from tested_python_streets.cycle_terminal_ports import requested_station
        self.assertEqual(at,requested_station(old,stations))
        expected=footprint(model.Key('terminal',10),old,at,4)
        actual=footprint(model.Key('terminal',10),a,at,4)
        np.testing.assert_array_equal(actual.region.xy,expected.region.xy)
        self.assertEqual(actual.segments,expected.segments)
        # The analytical original curbs, not a proposal field's apparent gap,
        # decide whether two separated source profiles may be retained.
        a.gap=lambda station:100.
        self.assertIsNone(original_section_gap(a,at))
        with self.assertRaisesRegex(ValueError,'finite original source support'):
            original_section_gap(a,a.length+1.)

    def test_original_profiles_and_complete_strip_union_survive_storage_and_rigid_changes(self):
        for mode in (0,3):
            for reverse,reflection,angle,translation in ((False,False,0.,(0.,0.)),(True,True,.713,(1800.,-1200.))):
                with self.subTest(mode=mode,reverse=reverse):
                    c,pair,sites=fixture(mode,reverse,reflection,angle,translation)
                    before=c.network.to_dict();a=pair.build(c,{})
                    station=minimum_station(a)+6.
                    f=footprint(model.Key('terminal',10),a,station,4)
                    loops=RegionUnion().build((f.region,));g=CompoundGroup(tuple(s.key for s in sites),(),(f,),loops)
                    self.assertEqual(len(certify(g)),1)
                    self.assertTrue(any(s.kind=='internal' for s in f.segments))
                    self.assertEqual(sum(s.first-s.last for s in f.segments if s.kind=='mouth'),1.)
                    self.assertLessEqual(a.source_curve_tolerance,.00025)
                    np.testing.assert_allclose(a.profile.value(np.array((-1.,1.))),0.,atol=1e-14)
                    aliases=a.source_mouths(station)
                    self.assertEqual({x.endpoint for x in aliases},set(pair.endpoints))
                    u=a.mouth_parameters(4,1.)
                    for alias in aliases:
                        path=c[model.Key('path',alias.endpoint.section_id)]
                        self.assertAlmostEqual(float(a.frame_chart.local(path.frame(alias.station)[0])[0]),float(a.chart_station(station)),places=8)
                        selected=u[alias.contains(u)]
                        p=c[model.Key('profile',alias.endpoint.section_id)].oriented(alias.endpoint.end)
                        np.testing.assert_allclose(p.value(2.*alias.outward_parameter(selected)-1.),a.profile.value(2.*selected-1.),atol=1e-12)
                    self.assertEqual(c.network.to_dict(),before)

    def test_dynamic_crown_gradient_is_the_actual_transverse_field_derivative(self):
        c,pair,_=fixture();a=pair.build(c,{})
        at=minimum_station(a)+5.;a.profile=a.cross_section(at)
        xy=a.point(at,a.width*np.array((-.3,-.15,.15,.3)))
        actual=a.field.gradient(xy)
        for axis in np.eye(2):
            h=1e-5
            finite=(a.field.values(xy+h*axis)-a.field.values(xy-h*axis))/(2*h)
            np.testing.assert_allclose(actual@axis,finite,atol=1e-8,rtol=1e-5)

    def test_only_both_original_terminal_owners_and_overlapping_finite_profiles_qualify(self):
        c,pair,sites=fixture();planner=SimpleNamespace(context=c)
        self.assertEqual(len(discover(planner,tuple(s.key for s in sites),set(),pair.originals)),1)
        self.assertFalse(discover(planner,(sites[0].key,),set(),pair.originals))
        self.assertFalse(discover(planner,tuple(s.key for s in sites),{1},pair.originals))
        self.assertFalse(discover(planner,tuple(s.key for s in sites),set(),pair.originals[:1]))
        for field in ('scope_ends','selection_ends'):
            changed=replace(c.network,sections=tuple(replace(s,**{field:(0,)}) if s.id==1 else s for s in c.network.sections),
                intersections=tuple(s for s in c.network.intersections if model.End(1,0) not in s.ends))
            modified=BuildContext(changed,c.products,c.patterns)
            self.assertFalse(discover(SimpleNamespace(context=modified),tuple(s.key for s in sites),set(),pair.originals))
        c2,pair2,sites2=fixture(separation=14.)
        self.assertFalse(discover(SimpleNamespace(context=c2),tuple(s.key for s in sites2),set(),pair2.originals))
        a=pair.build(c,{})
        for at in (minimum_station(a)-.01,a.length,a.length+1.):
            with self.assertRaises(ValueError):footprint(model.Key('terminal',10),a,at,4)
        for error in (0.,float('nan'),.001):
            with self.assertRaises(ValueError):FinitePairedEntrance(c,pair.run,0,source_tolerance=error)
        for station in (float('nan'),-1.,c[model.Key('path',1)].length+1.):
            with self.assertRaises(ValueError):FinitePairedEntrance(c,pair.run,0,source_tolerance=.00025,source_parameters={1:(station,)})
        wrong=replace(pair.run,members=(replace(pair.run.members[0],ends=(model.End(1,1),model.End(1,1))),pair.run.members[1]))
        with self.assertRaisesRegex(ValueError,'incidence'):FinitePairedEntrance(c,wrong,0,source_tolerance=.00025)

    def test_continuous_proof_rejects_a_hole_wholly_inside_a_source_prefix(self):
        c,pair,sites=fixture();a=pair.build(c,{})
        f=footprint(model.Key('terminal',10),a,minimum_station(a)+8.,4)
        loops=RegionUnion().build((f.region,))
        point=np.mean(f.layout.terminal_source_polygons[0],axis=0)
        hole=point+np.array(((-.02,-.02),(-.02,.02),(.02,.02),(.02,-.02)))
        g=CompoundGroup(tuple(s.key for s in sites),(),(f,),(*loops,BoundaryLoop(hole,())))
        with self.assertRaisesRegex(ValueError,'source strip'):certify(g)


if __name__=='__main__':unittest.main()
