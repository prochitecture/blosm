"""Original-incidence contact supports must never become invented road seams."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,End,Intersection,Network,Key
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.updates import UpdateCoordinator,BuildContext
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.routes import Route,Piece
from tested_python_streets.intersections import Approach
from tested_python_streets.compound_contact_routes import extend,regular_contact,outward_station


def fixture(*,reverse=False,duplicate=False,width=6.,scope=False,unrelated=False):
    sections=[Section(1,1,((0.,0.,0.),(10.,0.,0.))),
        Section(2,2,((10.,0.,0.),(30.,3.,0.)),width=width),
        Section(3,3,((10.,0.,0.),(10.,20.,0.))),
        Section(4,4,((30.,3.,0.),(50.,3.,0.)))]
    second_end=End(2,0)
    if reverse:
        sections[1]=replace(sections[1],points=sections[1].points[::-1]);second_end=End(2,1)
    if scope:sections[0]=replace(sections[0],scope_ends=(1,))
    ends=[End(1,1),second_end,End(3,0)]
    if duplicate:
        sections.append(Section(5,5,((10.,0.,0.),(30.,-3.,0.))));ends.append(End(5,0))
    if unrelated:
        ends.remove(second_end)
        sections.append(Section(6,6,((10.,0.,0.),(-10.,0.,0.)),width=2.));ends.append(End(6,0))
    if scope:ends.remove(End(1,1))
    sites=[Intersection(1,tuple(ends)),Intersection(2,(End(2,0 if reverse else 1),End(4,0)))]
    n=Network('contact_route',tuple(sections),tuple(sites));patterns=PatternLibrary()
    products=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(n).products
    context=BuildContext(n,products,patterns)
    a=Approach(Route((Piece(context[Key('path',1)],True),)),0,context[Key('profile',1)])
    return n,context,a


class ContactRouteTests(unittest.TestCase):
    def test_cross_street_ids_and_bent_frames_keep_original_geometry_and_identity(self):
        for reverse in (False,True):
            with self.subTest(reverse=reverse):
                n,c,a=fixture(reverse=reverse);before=n.to_dict()
                b=extend(c,a,{Key('intersection',1),Key('intersection',2)},{1})
                self.assertIs(b.source,a.source);self.assertEqual(b.end,a.end)
                self.assertEqual([p.section.id for p in b.path.pieces],[1,2])
                self.assertEqual(b.path.contact_extension_through,(Key('intersection',1),))
                self.assertGreater(b.length,a.length)
                at=a.length+.7
                expected=b.path.pieces[-1].frame(at-a.length)
                for actual,wanted in zip(b.frame(at),expected):np.testing.assert_array_equal(actual,wanted)
                self.assertFalse(regular_contact(b,a.length));self.assertTrue(regular_contact(b,at))
                self.assertEqual(n.to_dict(),before)

    def test_ambiguous_changed_width_scope_or_unrelated_positions_do_not_extend(self):
        for options in (dict(duplicate=True),dict(width=5.),dict(scope=True),dict(unrelated=True)):
            with self.subTest(options=options):
                _,c,a=fixture(**options)
                self.assertIs(extend(c,a,{Key('intersection',1)},{1}),a)
        _,c,a=fixture()
        self.assertIs(extend(c,a,set(),{1}),a)
        self.assertIs(extend(c,a,{Key('intersection',1)},set()),a)

    def test_original_contact_aliases_move_the_exterior_mouth_in_both_directions(self):
        for reverse in (False,True):
            n,c,a=fixture(reverse=reverse)
            b=extend(c,a,{Key('intersection',1)},{1})
            trims=b.path.trims(a.length+2.5)
            exterior=Approach(Route((b.path.pieces[-1],)),1 if reverse else 0,c[Key('profile',2)])
            aliases={key:(value,) for key,value in trims.items()}
            self.assertAlmostEqual(outward_station(exterior,aliases),2.5)
            self.assertEqual(outward_station(exterior,{}),0.)


if __name__=='__main__':unittest.main()
