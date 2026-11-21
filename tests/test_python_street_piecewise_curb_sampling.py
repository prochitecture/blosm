"""Separate lead/arc bounds retain physical curves and complete profiles."""
from dataclasses import replace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.compounds import CompoundGroup, rim_label, sample_footprint
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.intersections import Approach, BoundaryCurve, JunctionLayout, Mouth, _corner
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.paths import Path
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import Crown
from tested_python_streets.regions import RegionUnion
from tested_python_streets.sampling import AnalyticalBoundarySampler, CertifiedCurbSampler


def approaches(mode=0, *, reverse=False, angle=0., translation=(0., 0.), straight=False):
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    shapes=(((0.,0.),(120.,0.)) if straight else ((0.,0.),(30.,0.),(60.,.6),(120.,.6)),
            ((0.,0.),(0.,120.)), ((0.,0.),(-120.,0.)))
    result=[]
    for sid,points in enumerate(shapes,1):
        xy=np.asarray(points)@rotation.T+translation
        if reverse:xy=xy[::-1]
        source=model.Section(sid,sid,tuple((x,y,0.) for x,y in xy),width=2.,interpolation=mode)
        result.append(Approach(Path(source),int(reverse),Crown(.08,0.)))
    return tuple(result)


def ordinary_group(*, reverse=False, angle=0., translation=(0., 0.)):
    incoming=approaches(reverse=reverse,angle=angle,translation=translation)
    corners=tuple(_corner(a,incoming[(i+1)%3],1.) for i,a in enumerate(incoming))
    stations=(94.7,12.,15.)
    curves=tuple(BoundaryCurve(c,(stations[i],stations[(i+1)%3]),.005) for i,c in enumerate(corners))
    mouths={(a.source.id,a.end):Mouth(s if a.end==0 else a.length-s,8)
            for a,s in zip(incoming,stations)}
    layout=JunctionLayout(incoming,corners,curves,stations,mouths,1024,4,
                          incoming[0].origin.copy(),[0.]*3,mesh_strategy='boundary_chart')
    owner=model.Key('intersection',100)
    footprint=sample_footprint(owner,layout)
    loops=RegionUnion().build((footprint.region,))
    return CompoundGroup((owner,),(),(footprint,),loops),layout


class PiecewiseCurbSamplingTests(unittest.TestCase):
    def check_curve_geometry(self,curve,parameters,tolerance,spacing):
        points=curve(parameters)
        delta=np.diff(points,axis=0)
        lengths=np.linalg.norm(delta,axis=1)
        self.assertTrue(np.all(lengths>1e-7))
        self.assertLessEqual(float(lengths.max()),spacing+1e-8)
        q=np.linspace(0.,1.,33)
        actual=curve(parameters[:-1,None]+q*(np.diff(parameters)[:,None]))
        offset=actual-points[:-1,None,:]
        projection=np.clip(np.sum(offset*delta[:,None,:],axis=-1)/lengths[:,None]**2,0.,1.)
        distance=np.linalg.norm(offset-projection[...,None]*delta[:,None,:],axis=-1)
        self.assertLessEqual(float(distance.max()),tolerance+1e-8)

    def test_skewed_leads_keep_tight_arc_accuracy_and_all_exact_source_knots(self):
        for mode,reverse,angle,offset in ((0,False,0.,(0.,0.)),(0,True,.713,(1200.,-1800.)),
                                        (3,False,.23,(-350.,700.)),(3,True,.713,(1200.,-1800.))):
            with self.subTest(mode=mode,reverse=reverse):
                a,b,_=approaches(mode,reverse=reverse,angle=angle,translation=offset)
                old=(a.source.points,b.source.points)
                curve=BoundaryCurve(_corner(a,b,1.),(94.7,12.),.005)
                sampler=CertifiedCurbSampler(.005,.7,16384,piecewise=True)
                parameters=sampler.sample(curve,0.,1.)
                legacy=CertifiedCurbSampler(.005,.7,16384).sample(curve,0.,1.)
                self.assertLess(len(parameters),.6*len(legacy))
                self.assertEqual((len(parameters)-1)%2,0)
                np.testing.assert_array_equal(parameters[[0,-1]],(0.,1.))
                for contact in curve.lengths[1:3]/curve.length:
                    self.assertTrue(np.any(parameters==contact))
                self.check_curve_geometry(curve,parameters,.005,.7)
                if mode==0:
                    for lead,base in zip(curve.leads,(0.,curve.lengths[2])):
                        approach,source_parameters,cumulative,distance=lead
                        for knot in sampler.lead_metrics(approach,source_parameters,cumulative,
                                                        distance,base,curve.length)[2]:
                            self.assertTrue(np.any(parameters==knot))
                self.assertEqual(len(sampler.last_piecewise_bounds),3)
                self.assertTrue(all(piece['chord_error_bound']<=.005+1e-12
                                    for piece in sampler.last_piecewise_bounds))
                self.assertEqual(old,(a.source.points,b.source.points))
                clipped=sampler.sample(curve,.013,.987)
                np.testing.assert_array_equal(clipped[[0,-1]],(.013,.987))
                self.check_curve_geometry(curve,clipped,.005,.7)

    def test_unused_distant_bend_does_not_densify_a_straight_retained_lead(self):
        for reverse in (False,True):
            shapes=(((0.,0.),(40.,0.),(55.,0.),(55.,8.),(65.,8.)),((0.,0.),(0.,70.)))
            incoming=[]
            for sid,points in enumerate(shapes,1):
                xy=np.asarray(points)[::-1] if reverse else np.asarray(points)
                section=model.Section(sid,sid,tuple((x,y,0.) for x,y in xy),width=2.,interpolation=3)
                incoming.append(Approach(Path(section),int(reverse),Crown(.08,0.)))
            a,b=incoming
            curve=BoundaryCurve(_corner(a,b,1.),(30.,10.),.005)
            sampler=CertifiedCurbSampler(.005,.7,16384,piecewise=True)
            lead=curve.leads[0]
            legacy=sampler.lead_metrics(*lead,0.,curve.length)
            local=sampler.lead_metrics(*lead,0.,curve.length,local_curvature=True)
            self.assertGreater(legacy[0],.1)
            self.assertLess(local[0],1e-8)
            parameters=sampler.sample(curve,0.,1.)
            self.check_curve_geometry(curve,parameters,.005,.7)
            self.assertLess(len(parameters),200)

    def test_ordinary_boundary_keeps_every_complete_mouth_and_semantic_curb_endpoint(self):
        for reverse,angle,offset in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            with self.subTest(reverse=reverse):
                group,layout=ordinary_group(reverse=reverse,angle=angle,translation=offset)
                original=group.loops[0].xy.copy()
                sampler=AnalyticalBoundarySampler(piecewise=True)
                sampled=sampler.build(group)
                repeated=sampler.build(group)
                self.assertIsNone(sampled.physical_layout)
                self.assertEqual(sampled.members,group.members);self.assertEqual(sampled.claims,())
                self.assertEqual(len(sampled.loops),1)
                self.assertFalse(sampled.loops[0].hole)
                np.testing.assert_array_equal(sampled.loops[0].xy,repeated.loops[0].xy)
                np.testing.assert_array_equal(group.loops[0].xy,original)
                layouts={group.members[0]:layout}
                before=BoundaryConditions().build(group,group.loops[0],layouts)
                after=BoundaryConditions().build(sampled,sampled.loops[0],layouts)
                self.assertEqual(len(before.ports),3);self.assertEqual(len(after.ports),3)
                for port in before.ports:
                    other=next(p for p in after.ports if (p.section_id,p.end)==(port.section_id,port.end))
                    self.assertEqual(port.station,other.station);self.assertEqual(port.source_u,other.source_u)
                    for field in ('xy','heights','gradients'):
                        np.testing.assert_array_equal(getattr(before,field)[list(port.vertices)],
                                                      getattr(after,field)[list(other.vertices)])
                inventory={s.region.owner:s for s in sampled.footprints}
                labels=[rim_label(span,inventory) for span in sampled.loops[0].spans]
                for index,curve in enumerate(layout.boundaries):
                    own=[label for label in labels if label.kind=='curb' and label.index==index]
                    self.assertEqual(own[0].first,0.);self.assertEqual(own[-1].last,1.)
                    if len(curve.lengths)==4:
                        for contact in curve.lengths[1:3]/curve.length:
                            self.assertTrue(any(label.first==contact or label.last==contact for label in own))
                self.assertTrue(sampler.curve_reports)

    def test_default_sampling_retains_the_original_whole_curve_grid(self):
        a,b,_=approaches(straight=True)
        curve=BoundaryCurve(_corner(a,b,1.),(94.7,12.),.005)
        default=CertifiedCurbSampler(.005,.7,16384).sample(curve,0.,1.)
        explicit=CertifiedCurbSampler(.005,.7,16384,piecewise=False).sample(curve,0.,1.)
        # This straight-lead fixture's legacy whole-curve bound is the 1 m
        # arc's .2 m chord step, applied to the entire 104.270796 m curve.
        np.testing.assert_array_equal(default,np.linspace(0.,1.,523))
        np.testing.assert_array_equal(default,explicit)
        group,_=ordinary_group()
        first=AnalyticalBoundarySampler().build(group)
        second=AnalyticalBoundarySampler(piecewise=False).build(group)
        np.testing.assert_array_equal(first.loops[0].xy,second.loops[0].xy)
        self.assertEqual(first.loops[0].spans,second.loops[0].spans)

    def test_piecewise_allocation_and_recursive_nonordinary_retries_remain_bounded(self):
        a,b,_=approaches();curve=BoundaryCurve(_corner(a,b,1.),(94.7,12.),.005)
        for budget in (4,140):
            with self.subTest(curve_budget=budget),self.assertRaisesRegex(ValueError,'sample budget'):
                CertifiedCurbSampler(.005,.7,budget,piecewise=True).sample(curve,0.,1.)
        group,_=ordinary_group()
        for budget in (32,200):
            with self.subTest(boundary_budget=budget),self.assertRaisesRegex(ValueError,'sample budget'):
                AnalyticalBoundarySampler(max_samples=budget,piecewise=True).build(group)
        patterns=PatternLibrary()
        invalid=(dict(piecewise_sampling=True),dict(four_sided=True))
        for flags in invalid:
            with self.assertRaisesRegex(ValueError,'unchanged local junction'):
                CompoundQuadChart(patterns,**flags)._piecewise_ordinary(group,None,'old bounded chart failed')
        alternatives=(replace(group,members=(*group.members,model.Key('intersection',101))),
                      replace(group,claims=(object(),)),replace(group,physical_layout=group.footprints[0].layout),
                      replace(group,loops=(*group.loops,replace(group.loops[0],xy=group.loops[0].xy[::-1]))))
        for changed in alternatives:
            with self.assertRaisesRegex(ValueError,'unchanged local junction'):
                CompoundQuadChart(patterns)._piecewise_ordinary(changed,None,'old bounded chart failed')


if __name__=='__main__':unittest.main()
