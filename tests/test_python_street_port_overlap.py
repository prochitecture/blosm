"""Crossed complete mouths require their actual common source neighbor."""
from types import SimpleNamespace
import unittest
import _python_street_package
from tested_python_streets.model import Section,End,Intersection,Network,Key
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider,Crown
from tested_python_streets.updates import UpdateCoordinator,BuildContext
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.routes import Route,RouteTopology
from tested_python_streets.intersections import Approach
from tested_python_streets.compound_port_overlap import shared_port_neighbors


def fixture(shared=True,half=0):
    def section(sid,a,b,width=2.):
        return Section(sid,sid,((*a,0.),(*b,0.)),width=width)
    sections=[section(1,(0.,0.),(10.,10.),6.),section(2,(0.,5.),(10.,10.),6.),
        section(3,(10.,10.),(30.,10.)),section(4,(-20.,0.),(0.,0.)),
        section(5,(0.,0.),(0.,-20.)),section(6,(-20.,5.),(0.,5.)),
        section(7,(0.,5.),(0.,25.))]
    sites=[Intersection(1,(End(1,0),End(4,1),End(5,0))),
           Intersection(2,(End(2,0),End(6,1),End(7,0)))]
    if shared:
        sites.append(Intersection(3,(End(1,1),End(2,1),End(3,0))))
    else:
        # Distinct source junction identities at the same XY must stay distinct.
        sections.extend((section(8,(10.,10.),(30.,11.)),section(9,(10.,10.),(20.,30.)),
                         section(10,(10.,10.),(20.,-10.))))
        sites.extend((Intersection(3,(End(1,1),End(3,0),End(10,0))),
                      Intersection(4,(End(2,1),End(8,0),End(9,0)))))
    network=Network('shared_mouth_neighbor',tuple(sections),tuple(sites))
    patterns=PatternLibrary()
    products=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(network).products
    context=BuildContext(network,products,patterns)
    approaches=[Approach(Route.build(context,RouteTopology(network).trace(End(sid,0))),0,
        Crown(.08,0.,half if sid==1 else -half)) for sid in (1,2)]
    layout=SimpleNamespace(approaches=approaches,radial=4,
        outward_stations=(approaches[0].length-1.,approaches[1].length-1.2))
    planner=SimpleNamespace(topology=RouteTopology(network),sites={j.key:j for j in sites})
    return network,planner,layout


class PortOverlapTests(unittest.TestCase):
    def test_complete_full_or_half_mouths_expand_only_to_their_real_shared_node(self):
        for half in (0,1):
            with self.subTest(half=half):
                n,planner,layout=fixture(half=half)
                before=n.to_dict()
                self.assertEqual(shared_port_neighbors(planner,{Key('intersection',1),Key('intersection',2)},layout),
                    {Key('intersection',3)})
                self.assertEqual(n.to_dict(),before)

    def test_equal_coordinates_do_not_invent_a_common_original_junction(self):
        _,planner,layout=fixture(shared=False)
        self.assertEqual(shared_port_neighbors(planner,{Key('intersection',1),Key('intersection',2)},layout),set())

    def test_separate_mouths_or_an_already_owned_neighbor_do_not_expand(self):
        _,planner,layout=fixture()
        members={Key('intersection',1),Key('intersection',2)}
        self.assertFalse(shared_port_neighbors(planner,members|{Key('intersection',3)},layout))
        layout.outward_stations=(0.,0.)
        self.assertFalse(shared_port_neighbors(planner,members,layout))


if __name__=='__main__':unittest.main()
