"""Persistent approach overlap is a bounded geometric planning condition."""
import unittest
from dataclasses import replace
import _python_street_package
from tested_python_streets import model
from tested_python_streets.intersections import Approach
from tested_python_streets.overlapping_approaches import persistent_overlap
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.routes import Piece, Route
from tested_python_streets.compound_preparation import DeferredCompoundPlanner


def approach(sid,points,*,width=5.,mode=3):
    section=model.Section(sid,sid,tuple((*p,0.) for p in points),width=width,interpolation=mode)
    return Approach(Route((Piece(Path(section),True),)),0,Crown(.08,0.))


class PersistentApproachTests(unittest.TestCase):
    def test_complete_shorter_route_stays_inside_combined_width(self):
        a=approach(1,((0.,0.),(20.,.4),(40.,1.)))
        b=approach(2,((0.,0.),(25.,2.),(60.,2.)))
        self.assertTrue(persistent_overlap(a,b))
        self.assertFalse(persistent_overlap(b,a),'A long tail cannot inherit an unfinished overlap claim')

    def test_interior_excursion_is_not_accepted_from_only_coincident_endpoints(self):
        a=approach(1,((0.,0.),(30.,12.),(60.,0.)))
        b=approach(2,((0.,0.),(60.,0.)))
        self.assertFalse(persistent_overlap(a,b))

    def test_brief_gap_can_be_too_small_for_the_required_fillet(self):
        a=approach(1,((0.,0.),(40.,0.)))
        b=approach(2,((0.,0.),(20.,5.1),(50.,5.1)))
        self.assertFalse(persistent_overlap(a,b))
        self.assertTrue(persistent_overlap(a,b,clearance=1.))
        with self.assertRaises(ValueError):persistent_overlap(a,b,clearance=-1.)

    def test_distinct_origins_exact_paths_and_unbounded_work_are_not_candidates(self):
        a=approach(1,((0.,0.),(20.,0.)))
        self.assertFalse(persistent_overlap(a,approach(2,((0.,.01),(25.,.01)))))
        self.assertFalse(persistent_overlap(a,approach(2,((0.,0.),(25.,0.)),mode=0)))
        self.assertFalse(persistent_overlap(a,approach(2,((0.,0.),(2000.,0.)))))

    def test_next_hop_requires_a_real_terminal_and_an_owned_original_chain(self):
        from types import SimpleNamespace
        from tested_python_streets.patterns import PatternLibrary
        from tested_python_streets.providers import PathProvider
        from tested_python_streets.profiles import ProfileProvider
        from tested_python_streets.updates import UpdateCoordinator,BuildContext
        points=(((0.,0.),(20.,0.)),((20.,0.),(60.,1.)),((0.,0.),(50.,3.)),
                ((0.,0.),(-20.,-10.)),((20.,0.),(20.,-20.)),
                ((50.,3.),(70.,20.)),((50.,3.),(70.,-20.)))
        sections=tuple(model.Section(i,i,tuple((*p,0.) for p in ps),width=5.) for i,ps in enumerate(points,1))
        sites=(model.Intersection(1,(model.End(1,0),model.End(3,0),model.End(4,0))),
               model.Intersection(2,(model.End(1,1),model.End(2,0),model.End(5,0))),
               model.Intersection(3,(model.End(3,1),model.End(6,0),model.End(7,0))))
        network=model.Network('next-source-hop',sections,sites)
        patterns=PatternLibrary()
        products=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(network).products
        context=BuildContext(network,products,patterns)
        def arm(sid):
            return Approach(Route((Piece(products[model.Key('path',sid)],True),)),0,products[model.Key('profile',sid)])
        planner=DeferredCompoundPlanner(context,sites,{})
        self.assertEqual(planner.next_source_hops((sites[0].key,sites[1].key),[arm(i) for i in (3,4,2,5)]),
                         {sites[2].key})
        # The same geometric proximity cannot invent a terminal junction.
        context.network=replace(network,intersections=sites[:2])
        planner=DeferredCompoundPlanner(context,sites[:2],{})
        self.assertEqual(planner.next_source_hops((sites[0].key,sites[1].key),[arm(i) for i in (3,4,2,5)]),set())


if __name__=='__main__':unittest.main()
