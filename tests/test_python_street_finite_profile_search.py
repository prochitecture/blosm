"""Complete finite original profiles advance without moving other mouths."""
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,End,Intersection,Network,Settings,Key
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.routes import Route,Piece,RouteTopology
from tested_python_streets.intersections import Approach,Mouth
from tested_python_streets.regions import Region
from tested_python_streets.updates import BuildContext
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.finite_profile_search import advance,FiniteProfileUnavailable
from tested_python_streets.compound_preparation import CompoundExpansionRequired


def fixture(*,reverse=False,angle=0.,blocked_until=80.5,scope=False,diverging=False):
    rot=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    def xy(points):return np.asarray(points)@rot.T+(1200.,-1800.)
    controls=[((0.,0.),(100.,0.)),((0.,0.),(100.,20.) if diverging else (0.,20.)),((0.,0.),(0.,-20.)),
              ((100.,0.),(120.,10.)),((100.,0.),(120.,-10.))]
    sections=[]
    for sid,points in enumerate(controls,1):
        points=xy(points)
        if sid==1 and reverse:points=points[::-1]
        sections.append(Section(sid,sid,tuple((x,y,0.) for x,y in points),width=4.,interpolation=0,
            scope_ends=(0 if reverse else 1,) if sid==1 and scope else ()))
    root=End(1,int(reverse));far=End(1,1-int(reverse))
    own=Intersection(1,(root,End(2,0),End(3,0)))
    remote=Intersection(2,(far,End(4,0),End(5,0)))
    n=Network('finite-profile',tuple(sections),(own,) if scope else (own,remote),settings=Settings())
    products={Key('path',s.id):Path(s) for s in sections}
    context=BuildContext(n,products,PatternLibrary())
    a=Approach(Route((Piece(products[Key('path',1)],not reverse),)),int(reverse),Crown(.08,0.))
    b=Approach(Route((Piece(products[Key('path',2)],True),)),0,Crown(.08,0.))
    def source(owner,points):return SimpleNamespace(region=Region(Key('test_source',owner),xy(points)))
    outward=[16.,3.] if diverging else [35.,7.]
    ports=tuple(SimpleNamespace(region=Region(Key('test_source',i),
        np.array((arm.point(0.,2.),arm.point(station,2.),arm.point(station,-2.),arm.point(0.,-2.)))))
        for i,(arm,station) in enumerate(zip((a,b),outward),1))
    # Neither endpoint nor the middle of the transverse mouth intersects this
    # real finite half-metre strip. Its whole-profile crossing must be found.
    blocker=source(3,((25.,.5),(blocked_until,.5),(blocked_until,1.),(25.,1.)))
    remote_layout=SimpleNamespace(mouths={(1,far.end):Mouth(5. if reverse else 95.,8)})
    planner=SimpleNamespace(context=context,topology=RouteTopology(n),layouts={own.key:SimpleNamespace(mouths={}),remote.key:remote_layout})
    stations={sid:np.linspace(0.,products[Key('path',sid)].length,101) for sid in (1,2)}
    return planner,(a,b),outward,ports,ports if diverging else (*ports,blocker),stations,(own.key,),remote.key


class FiniteProfileSearchTests(unittest.TestCase):
    def test_coupled_rebuilt_clips_use_only_full_finite_proof_cells_and_preserve_original_controls(self):
        from tested_python_streets.regions import RegionUnion
        from tested_python_streets.finite_containment import contains_segment
        for reverse,angle in ((False,0.),(True,.713)):
            planner,a,outward,ports,sources,stations,members,remote=fixture(reverse=reverse,angle=angle,diverging=True)
            before=planner.context.network.to_dict();ports=list(ports);current=list(outward)
            # The short second mouth is covered by the first actual finite
            # prefix. Only that proven failed mouth advances initially.
            first,report=advance(planner,a,current,(1,),ports,sources,stations,members,full_source_exclusion=True)
            self.assertEqual(first[0],current[0]);self.assertGreater(first[1],20.)
            self.assertTrue(report[0]['full_finite_source_exclusion'])
            def region(index,station):
                arm=a[index]
                return SimpleNamespace(region=Region(ports[index].region.owner,
                    np.array((arm.point(0.,2.),arm.point(station,2.),arm.point(station,-2.),arm.point(0.,-2.)))))
            ports[1]=region(1,first[1]);current=first
            # Its newly owned true prefix now clips the old first profile.
            # Full finite exclusion resolves that actual new clip once.
            result,report=advance(planner,a,current,(0,),ports,ports,stations,members,full_source_exclusion=True)
            self.assertEqual(result[1],current[1]);self.assertGreater(result[0],20.)
            ports[0]=region(0,result[0]);loops=RegionUnion().build(s.region for s in ports)
            self.assertEqual(len(loops),1)
            for arm,station in zip(a,result):
                self.assertTrue(contains_segment(*arm.point(station,np.array((-2.,2.))),loops[0].xy))
            self.assertEqual(planner.context.network.to_dict(),before)
            with self.assertRaisesRegex(ValueError,'candidate budget'):
                advance(planner,a,outward,(1,),ports,ports,stations,members,
                        max_candidates=2,full_source_exclusion=True)

    def test_long_available_original_lead_and_continuous_profile_keep_other_mouth_exact(self):
        for reverse,angle in ((False,0.),(True,.713)):
            planner,a,outward,ports,sources,stations,members,remote=fixture(reverse=reverse,angle=angle)
            before=planner.context.network.to_dict()
            result,report=advance(planner,a,outward,(0,),ports,sources,stations,members)
            self.assertAlmostEqual(result[0],81.,places=10)
            self.assertEqual(result[1],7.);self.assertEqual(outward,[35.,7.])
            self.assertGreater(report[0]['candidates_tested'],32)
            self.assertEqual(report[0]['complete_transverse_profile'],'continuous_outside_all_other_source_operands')
            self.assertFalse(report[0]['endpoint_extrapolation'])
            self.assertEqual(planner.context.network.to_dict(),before)

    def test_finite_scope_and_candidate_limits_cannot_be_relabelled_as_free(self):
        for scope in (False,True):
            planner,a,outward,ports,sources,stations,members,remote=fixture(scope=scope,blocked_until=101.)
            with self.assertRaisesRegex(FiniteProfileUnavailable,'finite boundary'):
                advance(planner,a,outward,(0,),ports,sources,stations,members)
            if scope:
                a[0].path.terminal_free=True
                with self.assertRaisesRegex(ValueError,'scope.*as free'):
                    advance(planner,a,outward,(0,),ports,sources,stations,members)
        planner,a,outward,ports,sources,stations,members,remote=fixture()
        with self.assertRaisesRegex(ValueError,'candidate budget'):
            advance(planner,a,outward,(0,),ports,sources,stations,members,max_candidates=2)

    def test_unknown_neighbor_is_rejected_and_known_opposite_claim_requires_closure(self):
        planner,a,outward,ports,sources,stations,members,remote=fixture()
        del planner.layouts[remote]
        with self.assertRaisesRegex(ValueError,'available original neighbor'):
            advance(planner,a,outward,(0,),ports,sources,stations,members)
        planner,a,outward,ports,sources,stations,members,remote=fixture(blocked_until=96.5)
        with self.assertRaises(CompoundExpansionRequired) as caught:
            advance(planner,a,outward,(0,),ports,sources,stations,members)
        self.assertEqual(caught.exception.neighbors,(remote,))


if __name__=='__main__':unittest.main()
