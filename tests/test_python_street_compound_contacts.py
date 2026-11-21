"""A curved near-through source may need the other local contact normal side."""
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,End,Intersection,Network,Key
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.updates import UpdateCoordinator,BuildContext
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.routes import Route,RouteTopology
from tested_python_streets.intersections import Approach
from tested_python_streets.corner_radii import CornerRadiusPolicy,CornerContactError
from tested_python_streets.corner_source_support import complete_source_roots
from tested_python_streets.compound_corner_contacts import retry
from tested_python_streets.corner_patch import corner_patch


def fixture(angle=0.,shift=(0.,0.)):
    origin=np.array((1610.6834716796875,1025.2938232421875))
    c,s=math.cos(angle),math.sin(angle);rotation=np.array(((c,-s),(s,c)))
    def transform(points):
        xy=(np.asarray(points)-origin)@rotation.T+shift
        return tuple((*p,0.) for p in xy)
    sections=(Section(1,1,transform(((1730.9678955078125,982.6461791992188),
        (1702.082763671875,979.109375),(1682.1099853515625,978.4583129882812),
        (1671.656005859375,981.2273559570312),(1657.2423095703125,987.9026489257812),
        (1643.3570556640625,995.9918823242188),tuple(origin))),width=5.4),
        Section(2,2,transform((tuple(origin),(1579.986083984375,1053.0267333984375))),width=5.4),
        Section(3,3,transform(((1662.5262451171875,977.1617431640625),
            (1622.7115478515625,1010.40234375),tuple(origin))),width=5.))
    site=Intersection(1,(End(1,1),End(2,0),End(3,1)))
    n=Network('curved_contact',sections,(site,));patterns=PatternLibrary()
    products=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(n).products
    context=BuildContext(n,products,patterns);topology=RouteTopology(n)
    approaches=tuple(Approach(Route.build(context,topology.trace(end)),end.end,
        context[Key('profile',end.section_id)]) for end in site.ends)
    return n,context,approaches


class CompoundContactTests(unittest.TestCase):
    def test_original_near_through_pair_keeps_its_complete_facing_sector(self):
        controls=(((0.,0.),(2.0301513671875,-1.837158203125),(4.37158203125,-3.25140380859375),
            (8.307861328125,-4.49896240234375),(11.2095947265625,-4.67767333984375),
            (14.095947265625,-4.2774658203125),(17.877197265625,-2.7642822265625)),
            ((-3.7713623046875,8.561279296875),(-3.2117919921875,5.622314453125),
                (-2.05322265625,2.861328125),(0.,0.)),
            ((0.,0.),(14.9854736328125,-30.23736572265625)))
        # Local contact unit test. Full original-network PLAN and strip-union
        # qualification are performed independently on the captured import.
        for angle,shift in ((0.,(0.,0.)),(.7,(1200.,-800.))):
            rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
            sections=tuple(Section(i,i,tuple((*p,0.) for p in np.asarray(points)@rotation.T+shift),
                width=5. if i<3 else 3.5) for i,points in enumerate(controls,1))
            site=Intersection(1,(End(1,0),End(2,1),End(3,0)))
            n=Network('facing-root-contact',sections,(site,));before=n.to_dict();patterns=PatternLibrary()
            products=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(n).products
            context=BuildContext(n,products,patterns);topology=RouteTopology(n)
            support=tuple(Approach(Route.build(context,topology.trace(end)),end.end,
                context[Key('profile',end.section_id)]) for end in site.ends)
            a,b=support[:2]
            full=lambda candidate:complete_source_roots(candidate,.7,.00025,support_approaches=support)
            with self.assertRaises(CornerContactError):CornerRadiusPolicy.solve(a,b,3.,'STYLE',accept=full)
            corner=retry(context,a,b,3.,'STYLE',{site.key},{1},support)
            self.assertIsNotNone(corner);self.assertEqual(corner.radius,3.)
            self.assertEqual(corner.fit_policy,'SOURCE_CORNER_FACING_ROOTS')
            self.assertTrue(complete_source_roots(corner,.7,.00025,support_approaches=support,
                audit_all_contacts=True,audit_facing_halves=True))
            corner_patch(Key('proof',1),corner,.7,.00025)
            self.assertEqual(n.to_dict(),before)
            self.assertIsNone(retry(context,a,b,3.,'STRICT',{site.key},{1},support))

    def test_secondary_side_keeps_original_curbs_and_requested_circle(self):
        for angle,shift in ((0.,(0.,0.)),(.7,(1500.,-800.))):
            with self.subTest(angle=angle):
                n,c,support=fixture(angle,shift);a,b=support[:2];before=n.to_dict()
                accept=lambda corner:complete_source_roots(corner,.7,.00025,support_approaches=support)
                with self.assertRaises(CornerContactError):
                    CornerRadiusPolicy.solve(a,b,2.,'STYLE',accept=accept)
                corner=retry(c,a,b,2.,'STYLE',{Key('intersection',1)},{1},support)
                self.assertIsNotNone(corner);self.assertEqual(corner.radius,2.)
                self.assertEqual(corner.fit_policy,'SOURCE_CURVED_CONTACT_SIDE')
                self.assertGreater(corner.stations[0],0.);self.assertLess(corner.stations[1],0.)
                self.assertTrue(complete_source_roots(corner,.7,.00025,support_approaches=support,
                    audit_all_contacts=True,audit_facing_halves=True))
                patch=corner_patch(Key('proof',0),corner,.7,.00025)
                self.assertGreater(len(patch.region.xy),3)
                self.assertEqual(n.to_dict(),before)

    def test_explicit_radius_does_not_enable_an_alternative_branch(self):
        _,c,support=fixture();a,b=support[:2]
        self.assertIsNone(retry(c,a,b,2.,'STRICT',{Key('intersection',1)},{1},support))


if __name__=='__main__':unittest.main()
