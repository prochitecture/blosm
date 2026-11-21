"""A bounded complete-route proof may request another real physical site."""
from dataclasses import replace
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.intersections import Approach
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.routes import Route,RouteTopology
from tested_python_streets.updates import UpdateCoordinator,BuildContext
from tested_python_streets.overlapping_approaches import persistent_overlap,overlap_certificate
from tested_python_streets.cycle_overlap_expansion import required


def fixture(*,near_miss=False,angle=0.,shift=(0.,0.)):
    origin=np.array((-1249.9444580078125,-685.5804443359375))
    controls=(np.array((tuple(origin),(-1252.1441650390625,-675.4276733398438),
        (-1251.6607666015625,-670.6632690429688),(-1249.6832275390625,-662.7377319335938),
        (-1247.052734375,-658.1073608398438))),
        np.array((tuple(origin),(-1249.08740234375,-678.7122192382812),(-1248.059326171875,-673.0685424804688),
            (-1246.97705078125,-668.615966796875),(-1245.4512939453125,-664.1857299804688),
            (-1242.9532470703125,-660.0117797851562),(-1241.0933837890625,-657.57421875),
            (-1238.1754150390625,-654.3799438476562),(-1234.603759765625,-650.82958984375),
            (-1225.73388671875,-645.4879760742188))))
    if near_miss:controls[0][1:,0]-=.2
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    points=[(c-origin)@rotation.T+shift for c in controls]
    section=lambda i,p:model.Section(i,i,tuple((*v,0.) for v in p),width=3.5)
    sections=[section(i,p) for i,p in enumerate(points,1)]
    a,b,c=points[0][0],points[0][-1],points[1][-1]
    for i,node,delta in ((3,a,(-20.,0.)),(4,b,(-20.,20.)),(5,b,(20.,20.)),
                         (6,c,(-20.,20.)),(7,c,(20.,20.))):
        sections.append(section(i,(node,node+np.asarray(delta))))
    sites=(model.Intersection(10,(model.End(1,0),model.End(2,0),model.End(3,0))),
        model.Intersection(11,(model.End(1,1),model.End(4,0),model.End(5,0))),
        model.Intersection(12,(model.End(2,1),model.End(6,0),model.End(7,0))))
    return model.Network('finite-overlap-source-expansion',tuple(sections),sites)


def planner(n):
    patterns=PatternLibrary();products=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(n).products
    context=BuildContext(n,products,patterns);topology=RouteTopology(n)
    approaches=tuple(Approach(Route.build(context,topology.trace(model.End(i,0))),0,
        context[model.Key('profile',i)]) for i in (1,2))
    return SimpleNamespace(context=context,topology=topology,sites={j.key:j for j in n.intersections}),approaches


class CycleOverlapExpansionTests(unittest.TestCase):
    def test_refined_continuous_proof_requests_only_the_real_external_terminal(self):
        for angle,shift in ((0.,(0.,0.)),(.7,(1200.,-800.))):
            n=fixture(angle=angle,shift=shift);before=n.to_dict();p,(inside,outside)=planner(n)
            self.assertFalse(persistent_overlap(inside,outside,clearance=1.))
            proof=overlap_certificate(inside,outside,clearance=1.,refine=True)
            self.assertIsNotNone(proof)
            self.assertLess(proof['maximum_distance_bound'],proof['permitted_distance'])
            self.assertLessEqual(proof['attempts'],3)
            self.assertLessEqual(max(proof['source_samples']),4096)
            self.assertFalse(proof['endpoint_extrapolation'])
            internal={model.End(1,0),model.End(1,1)};members=(model.Key('intersection',10),model.Key('intersection',11))
            result=required(p,members[0],outside,inside,10.,'STYLE',members,internal,{1})
            self.assertIsNotNone(result);terminal,record=result
            self.assertEqual(terminal,model.Key('intersection',12))
            self.assertEqual(record['original_internal_sections'],[1])
            self.assertEqual(record['original_external_sections'],[2])
            self.assertEqual(record['minimum_style_radius'],.5)
            self.assertEqual(before,n.to_dict())

    def test_a_near_miss_has_no_complete_source_corridor_certificate(self):
        n=fixture(near_miss=True);p,(inside,outside)=planner(n)
        self.assertFalse(persistent_overlap(inside,outside,clearance=1.,refine=True))
        self.assertIsNone(required(p,model.Key('intersection',10),outside,inside,10.,'STYLE',
            (model.Key('intersection',10),model.Key('intersection',11)),{model.End(1,0),model.End(1,1)},{1}))

    def test_strict_missing_original_incidence_or_unknown_terminal_cannot_expand(self):
        p,(inside,outside)=planner(fixture());site=model.Key('intersection',10)
        members=(site,model.Key('intersection',11));internal={model.End(1,0),model.End(1,1)}
        self.assertIsNone(required(p,site,outside,inside,10.,'STRICT',members,internal,{1}))
        self.assertIsNone(required(p,site,outside,inside,10.,'STYLE',members,{model.End(1,0)},{1}))
        self.assertIsNone(required(p,site,outside,inside,10.,'STYLE',members,internal,set()))
        p.sites.pop(model.Key('intersection',12))
        self.assertIsNone(required(p,site,outside,inside,10.,'STYLE',members,internal,{1}))

    def test_scope_selection_and_explicit_allocation_bounds_are_retained(self):
        for name in ('scope_ends','selection_ends'):
            n=fixture();n=replace(n,sections=tuple(replace(s,**{name:(1,)}) if s.id==2 else s for s in n.sections),
                intersections=tuple(j for j in n.intersections if j.id!=12))
            p,(inside,outside)=planner(n)
            self.assertIsNone(required(p,model.Key('intersection',10),outside,inside,10.,'STYLE',
                (model.Key('intersection',10),model.Key('intersection',11)),{model.End(1,0),model.End(1,1)},{1}))
        p,(inside,outside)=planner(fixture())
        with self.assertRaises(ValueError):overlap_certificate(inside,outside,clearance=-1.,refine=True)
        with self.assertRaises(ValueError):overlap_certificate(inside,outside,clearance=1.,refine=1)


if __name__=='__main__':unittest.main()
