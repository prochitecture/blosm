"""A complete original profile must also clear remote owner protrusions."""
import math
import unittest
from types import SimpleNamespace
import numpy as np
import _python_street_package
from test_python_street_finite_profile_search import fixture as base_fixture
from tested_python_streets.source_mouth_clearance import advance_past_overlap
from tested_python_streets.compound_preparation import CompoundExpansionRequired
from tested_python_streets.model import Key
from tested_python_streets.regions import Region


def fixture(**kwargs):
    planner,arms,outward,ports,sources,stations,members,remote=base_fixture(**kwargs)
    angle=kwargs.get('angle',0.);last=kwargs.get('blocked_until',80.5)
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    def source(owner,points):
        return SimpleNamespace(region=Region(Key('test_source',owner),
            np.asarray(points,float)@rotation.T+(1200.,-1800.)))
    # The remote strip joins the current source pool through this separate
    # finite region above the road. Its narrow crossing is past the exposed
    # full mouth at 35 m, so a mouth-only check cannot find the overlap.
    blocker=source(3,((45.,.5),(last,.5),(last,1.),(45.,1.)))
    connector=source(4,((0.,5.),(last-.5,5.),(last-.5,.5),(last+.5,.5),
                        (last+.5,6.),(0.,6.)))
    return planner,arms,outward,ports,(*ports,blocker,connector),stations,members,remote


class SourceMouthClearanceTests(unittest.TestCase):
    def test_complete_profile_advances_past_remote_whole_width_overlap(self):
        for reverse,angle in ((False,0.),(True,.713)):
            with self.subTest(reverse=reverse,angle=angle):
                planner,arms,outward,ports,sources,stations,members,remote=fixture(reverse=reverse,angle=angle)
                before=planner.context.network.to_dict();paths=tuple(p.path for a in arms for p in a.path.pieces)
                # This profile itself is clear. The owner's narrow original
                # source strip crosses the finite road again much farther on.
                result,reports=advance_past_overlap(planner,arms,outward,ports,sources,stations,members)
                self.assertAlmostEqual(result[0],82.,places=9);self.assertEqual(result[1],7.)
                self.assertEqual(outward,[35.,7.]);self.assertEqual(len(reports),1)
                self.assertGreater(reports[0]['last_overlapping_cell_end'],80.5)
                self.assertEqual(reports[0]['overlap_operands'],['test_source:3','test_source:4'])
                from tested_python_streets.regions import RegionUnion,_classify
                loops=RegionUnion().build(s.region for s in sources)
                self.assertEqual(sum(not loop.hole for loop in loops),1)
                mouth=arms[0].point(outward[0],np.linspace(-2.,2.,17))
                _,boundary=_classify(mouth,next(loop.xy for loop in loops if not loop.hole),1e-7)
                self.assertTrue(np.all(boundary>=0))
                a=arms[0]
                grown=SimpleNamespace(region=Region(ports[0].region.owner,np.array((
                    a.point(0.,2.),a.point(result[0],2.),a.point(result[0],-2.),a.point(0.,-2.)))))
                rebuilt=RegionUnion().build(s.region for s in (grown,*sources[1:]))
                _,boundary=_classify(a.point(result[0],np.linspace(-2.,2.,17)),
                    next(loop.xy for loop in rebuilt if not loop.hole),1e-7)
                self.assertTrue(np.all(boundary>=0))
                from tested_python_streets.source_strip_containment import certify
                cut=a.path.trims(result[0])[1,int(reverse)]
                proof=certify(planner.context,SimpleNamespace(loops=rebuilt),(1,),
                    intervals={1:(cut,a.length) if reverse else (0.,cut)})
                self.assertEqual(proof['maximum_error'],.002)
                self.assertFalse(reports[0]['endpoint_extrapolation'])
                self.assertEqual(planner.context.network.to_dict(),before)
                self.assertTrue(all(before is after for before,after in zip(paths,
                    (p.path for a in arms for p in a.path.pieces))))

    def test_contained_operand_and_boundary_only_contact_are_distinguished(self):
        planner,arms,outward,ports,sources,stations,members,remote=fixture()
        def source(points):
            return SimpleNamespace(region=Region(Key('clearance',1),np.asarray(points,float)+(1200.,-1800.)))
        contained=source(((40.4,-.1),(40.6,-.1),(40.6,.1),(40.4,.1)))
        result,reports=advance_past_overlap(planner,arms,outward,ports,(*ports,contained),stations,members)
        self.assertEqual(result,[42.,7.]);self.assertEqual(len(reports),1)
        for depth in (0.,.0001):
            touching=source(((25.,2.-depth),(80.5,2.-depth),(80.5,3.),(25.,3.)))
            for budget in (0,4096):
                result,reports=advance_past_overlap(planner,arms,outward,ports,(*ports,touching),stations,members,
                                                  max_candidates=budget)
                self.assertEqual(result,outward);self.assertEqual(reports,[])

    def test_existing_neighbor_scope_and_finite_budgets_remain_mandatory(self):
        for blocked in (96.5,101.):
            planner,arms,outward,ports,sources,stations,members,remote=fixture(blocked_until=blocked)
            planner.sites={key:None for key in planner.layouts}
            with self.subTest(blocked=blocked),self.assertRaises(CompoundExpansionRequired) as caught:
                advance_past_overlap(planner,arms,outward,ports,sources,stations,members)
            self.assertEqual(caught.exception.neighbors,(remote,))
            if blocked>100.:
                proof=caught.exception.source_mouth_clearance
                self.assertEqual(proof['original_neighbor'],remote.text())
                self.assertEqual(proof['overlap_operands'],['test_source:3'])
                self.assertFalse(proof['proof_only_future_operands'])
                self.assertFalse(proof['endpoint_extrapolation'])
                self.assertIn(proof,planner.source_mouth_expansions.values())
        planner,arms,outward,ports,sources,stations,members,remote=fixture(scope=True,blocked_until=101.)
        with self.assertRaisesRegex(ValueError,'finite boundary'):
            advance_past_overlap(planner,arms,outward,ports,sources,stations,members)
        planner,arms,outward,ports,sources,stations,members,remote=fixture()
        for kwargs,reason in ((dict(maximum_cells=2),'cell budget'),
                              (dict(maximum_pairs=1),'pair budget'),
                              (dict(max_candidates=0),'candidate budget')):
            with self.subTest(kwargs=kwargs),self.assertRaisesRegex(ValueError,reason):
                advance_past_overlap(planner,arms,outward,ports,sources,stations,members,**kwargs)


if __name__=='__main__':unittest.main()
