"""Complete source width, explicit holes and continuous boundary-band proofs."""
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,Network,End,Intersection,Settings
from tested_python_streets.paths import Path
from tested_python_streets.regions import BoundaryLoop
from tested_python_streets.source_strip_containment import certify,certify_intervals
from tested_python_streets.finite_containment import contains_segment


class Context:
    def __init__(self,network,path):self.network,self.path=network,path
    def __getitem__(self,key):return self.path


def context(points=((0.,0.,0.),(10.,0.,0.)),mode=3):
    section=Section(1,1,points,width=2.,interpolation=mode)
    closed=points[0]==points[-1]
    sites=(Intersection(1,(End(1,0),End(1,1))),) if closed else ()
    return Context(Network('width-proof',(section,),sites,settings=Settings(exaggeration=1.)),
        Path(section,returning=closed))


def group(outer,holes=()):
    return SimpleNamespace(loops=(BoundaryLoop(np.array(outer,float),()),
        *(BoundaryLoop(np.array(hole,float),()) for hole in holes)))


class SourceStripContainmentTests(unittest.TestCase):
    def test_remote_bend_does_not_exhaust_budget_for_a_straight_owned_interval(self):
        c = context(((0,0,0), (2,.15,0), (4,.5,0), (6,.5,0), (8,.5,0), (160,.5,0)))
        owner = group(((19.,-1.), (61.,-1.), (61.,2.), (19.,2.)))
        proof = certify_intervals(c, owner, {1:((20., 60.),)}, maximum_cells=85)
        self.assertLessEqual(proof['total_cells'], 85)
        self.assertEqual(proof['sections'][0]['maximum_combined_error'], .002)
        clipped = group(((19.,-.45), (61.,-.45), (61.,2.), (19.,2.)))
        with self.assertRaisesRegex(ValueError, 'clips the complete source strip'):
            certify_intervals(c, clipped, {1:((20., 60.),)}, maximum_cells=85)

    def test_explicit_original_station_interval_keeps_a_complete_partial_prefix(self):
        g=group(((4.99,-1.1),(10.1,-1.1),(10.1,1.1),(4.99,1.1)))
        for mode in (0,3):
            c=context(mode=mode)
            with self.assertRaises(ValueError):certify(c,g,{1})
            last=c.path.length;first=.5*last
            proof=certify(c,g,{1},intervals={1:(first,last)})
            self.assertEqual(proof['sections'][0]['interval'],[first,last])
            reverse=context(((10.,0.,0.),(0.,0.,0.)),mode)
            certify(reverse,g,{1},intervals={1:(0.,.5*reverse.path.length)})
            for invalid in ((-1.,5.),(5.,11.),(5.,5.)):
                with self.assertRaisesRegex(ValueError,'finite Section'):certify(c,g,{1},intervals={1:invalid})
            with self.assertRaisesRegex(ValueError,'declared original'):certify(c,g,{1},intervals={2:(0.,5.)})

    def test_complete_width_and_adjacent_hole_are_retained(self):
        g=group(((-1.,-2.),(11.,-2.),(11.,2.),(-1.,2.)),
            (((4.,1.2),(4.,1.8),(6.,1.8),(6.,1.2)),))
        proof=certify(context(),g,{1})
        self.assertEqual(proof['holes'],1)
        self.assertLessEqual(proof['maximum_error'],.002)

    def test_midspan_outer_width_loss_is_not_centerline_ownership(self):
        g=group(((-1.,-2.),(11.,-2.),(11.,2.),(5.135,2.),(5.135,.5),
            (5.125,.5),(5.125,2.),(-1.,2.)))
        with self.assertRaisesRegex(ValueError,'clips the complete source strip'):
            certify(context(),g,{1})

    def test_a_tiny_hole_inside_the_strip_cannot_hide_on_a_sampling_cross_section(self):
        g=group(((-1.,-2.),(11.,-2.),(11.,2.),(-1.,2.)),
            (((4.9999,-.0001),(4.9999,.0001),(5.0001,.0001),(5.0001,-.0001)),))
        with self.assertRaisesRegex(ValueError,'hole.*complete source strip'):
            certify(context(),g,{1})

    def test_an_outer_bay_with_a_submillimetre_neck_cannot_hide_inside_a_source_cell(self):
        g=group(((-1.,-2.),(11.,-2.),(11.,2.),(5.2501,2.),(5.2501,.8),
            (5.4,.8),(5.4,.2),(5.1,.2),(5.1,.8),(5.2499,.8),(5.2499,2.),(-1.,2.)))
        with self.assertRaisesRegex(ValueError,'exterior enters'):
            certify(context(),g,{1})

    def test_returning_exact_source_strip_keeps_its_real_central_hole(self):
        c=context(((0.,0.,0.),(10.,0.,0.),(10.,10.,0.),(0.,10.,0.),(0.,0.,0.)),0)
        g=group(((-2.,-2.),(12.,-2.),(12.,12.),(-2.,12.)),
            (((2.,2.),(2.,8.),(8.,8.),(8.,2.)),))
        certify(c,g,{1})

    def test_returning_prefixes_do_not_claim_the_unowned_middle(self):
        points=((0.,0.,0.),(10.,0.,0.),(10.,10.,0.),(0.,10.,0.),(0.,0.,0.))
        g=group(((-2.,-2.),(3.,-2.),(3.,3.),(-2.,3.)))
        for controls in (points,points[::-1]):
            c=context(controls,0);original=c.path;length=original.length
            proof=certify_intervals(c,g,{1:((0.,2.),(length-2.,length))})
            self.assertEqual([row['interval'] for row in proof['sections']],[[0.,2.],[length-2.,length]])
            self.assertEqual(proof['total_cells'],2)
            self.assertIs(c.path,original)
            with self.assertRaisesRegex(ValueError,'complete source strip'):
                certify_intervals(c,g,{1:((0.,length),)})
            with self.assertRaisesRegex(ValueError,'finite .*budget'):
                certify_intervals(c,g,{1:((0.,2.),(length-2.,length))},maximum_cells=1)
            full=group(((-2.,-2.),(12.,-2.),(12.,12.),(-2.,12.)))
            with self.assertRaisesRegex(ValueError,'finite cell budget'):
                certify_intervals(c,full,{1:((0.,12.),(length-12.,length))},maximum_cells=3)

    def test_interval_merge_and_whole_claim_keep_the_original_single_call_proof(self):
        c=context(mode=0);g=group(((-1.,-2.),(11.,-2.),(11.,2.),(-1.,2.)))
        for intervals,expected in ((((0.,2.),(2.,3.),(1.,2.5)),(0.,3.)),
                                   (((0.,2.),(8.,10.),(0.,10.)),(0.,10.))):
            actual=certify_intervals(c,g,{1:intervals})
            self.assertEqual(actual,certify(c,g,{1},intervals={1:expected}))
        proof=certify_intervals(c,g,{1:((0.,2.),(2.000000001,3.))})
        self.assertEqual(len(proof['sections']),2)  # A numerical gap is not ownership.
        for invalid in ((-1.,2.),(2.,11.),(2.,2.)):
            with self.assertRaisesRegex(ValueError,'finite Section'):
                certify_intervals(c,g,{1:(invalid,)})

    def test_boundary_band_is_continuous_and_not_a_near_midpoint_excuse(self):
        polygon=np.array(((-1.,-2.),(11.,-2.),(11.,2.),(8.,2.),(8.,-1.),
            (5.0001,-1.),(5.0001,-.0001),(4.9999,-.0001),(4.9999,-1.),(2.,-1.),(2.,2.),(-1.,2.)))
        self.assertFalse(contains_segment((0.,0.),(10.,0.),polygon,tolerance=.002))
        rectangle=np.array(((-1.,.001),(11.,.001),(11.,2.),(-1.,2.)))
        self.assertTrue(contains_segment((0.,0.),(10.,0.),rectangle,tolerance=.002))


if __name__=='__main__':unittest.main()
