"""Displaced source mouths retain incidence order and local median contacts."""
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.compound_preparation import DeferredCompoundPlanner
from tested_python_streets.corner_radii import CornerContactError, CornerRadiusPolicy
from tested_python_streets.intersections import Approach, _corner
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.routes import RouteTopology


def median_pair():
    # Slight convergence at the mouth, followed by ordinary curved streets.
    # A local median cap needs just over a semicircle; an ordinary fillet
    # must retain its existing <=180 degree sector constraint.
    first = model.Section(1, 1, ((0.,8.,0.),(-13.5,7.,0.),(-60.,7.,0.),
                                 (-71.,3.5,0.)), width=5.4)
    second = model.Section(2, 2, ((-48.,0.,0.),(-13.,0.,0.),(0.,0.,0.)), width=5.4)
    return (Approach(Path(first), 0, Crown(.08,0.,-1)),
            Approach(Path(second), 1, Crown(.08,0.,-1)))


class DisplacedCompoundTests(unittest.TestCase):
    def test_incidence_walk_does_not_swap_nearly_parallel_external_mouths(self):
        sections = tuple(model.Section(i,i,((0.,0.,0.),(20.,0.,0.))) for i in range(1,6))
        sites = (model.Intersection(1,tuple(model.End(i,0) for i in (1,2,3))),
                 model.Intersection(2,(model.End(1,1),model.End(4,0),model.End(5,0))))
        network = model.Network('rotation',sections,sites)
        def arm(sid,end,angle):
            return SimpleNamespace(source=network.section_index[sid],end=end,
                direction=np.array((math.cos(angle),math.sin(angle))))
        layouts = {sites[0].key:SimpleNamespace(approaches=(arm(1,0,1.),arm(2,0,3.),arm(3,0,5.))),
                   sites[1].key:SimpleNamespace(approaches=(arm(1,1,4.),arm(4,0,3.1),arm(5,0,.5)))}
        planner=DeferredCompoundPlanner(SimpleNamespace(network=network),sites,layouts)
        order=planner.exterior_order(tuple(j.key for j in sites),{model.End(1,0),model.End(1,1)})
        self.assertEqual(order,tuple(model.End(i,0) for i in (2,3,5,4)))
        self.assertNotEqual(order,tuple(model.End(i,0) for i in (5,2,4,3)))

    def test_paired_median_circle_is_local_tangent_and_not_an_ordinary_corner(self):
        first,second=median_pair()
        with self.assertRaises(CornerContactError):
            CornerRadiusPolicy.solve(first,second,2.,'STYLE')
        cap=CornerRadiusPolicy.solve(first,second,2.,'STYLE',displaced=True,median_cap=True)
        self.assertTrue(cap.median_cap)
        self.assertGreater(abs(cap.sweep),math.pi)
        self.assertLessEqual(abs(cap.sweep),math.radians(210.))
        self.assertGreaterEqual(min(cap.stations),0.)
        self.assertLess(max(cap.stations),10.)
        for a,s,distance in zip((first,second),cap.stations,(-first.width*.5,second.width*.5)):
            contact,tangent=a.jet(s,distance)
            radial=contact-cap.center
            self.assertAlmostEqual(float(np.linalg.norm(radial)),cap.radius,delta=1e-7)
            self.assertAlmostEqual(float(radial@tangent)/np.linalg.norm(tangent),0.,delta=1e-7)
        # The same exact contact cannot silently change the ordinary sector.
        with self.assertRaises(CornerContactError):
            _corner(first,second,cap.radius,side_override=1.)

    def test_explicit_median_radius_remains_strict_and_needs_displaced_contract(self):
        first,second=median_pair()
        with self.assertRaises(CornerContactError):
            CornerRadiusPolicy.solve(first,second,2.,'STRICT',displaced=True,median_cap=True)
        with self.assertRaisesRegex(ValueError,'explicit displaced'):
            CornerRadiusPolicy.solve(first,second,2.,'STYLE',median_cap=True)

    def test_internal_source_cycle_cannot_silently_become_a_filled_disk(self):
        sections=(model.Section(1,1,((0.,0.,0.),(10.,0.,0.))),
                  model.Section(2,2,((10.,0.,0.),(5.,10.,0.))),
                  model.Section(3,3,((5.,10.,0.),(0.,0.,0.))),
                  model.Section(4,4,((0.,0.,0.),(-20.,-10.,0.))),
                  model.Section(5,5,((10.,0.,0.),(30.,-10.,0.))),
                  model.Section(6,6,((5.,10.,0.),(5.,30.,0.))))
        sites=(model.Intersection(1,(model.End(1,0),model.End(3,1),model.End(4,0))),
               model.Intersection(2,(model.End(1,1),model.End(2,0),model.End(5,0))),
               model.Intersection(3,(model.End(2,1),model.End(3,0),model.End(6,0))))
        network=model.Network('internal-cycle',sections,sites)
        layouts={j.key:SimpleNamespace(approaches=tuple(Approach(Path(network.section_index[e.section_id]),
            e.end,Crown(.08,0.)) for e in j.ends)) for j in sites}
        planner=DeferredCompoundPlanner(SimpleNamespace(network=network),sites,layouts)
        with self.assertRaisesRegex(ValueError,'complete source-incidence cycle'):
            planner.exterior_order(tuple(j.key for j in sites),
                {model.End(sid,end) for sid in (1,2,3) for end in (0,1)})
        order,chains,cycles=planner.exterior_order(tuple(j.key for j in sites),
            {model.End(sid,end) for sid in (1,2,3) for end in (0,1)},
            include_chains=True,include_cycles=True)
        self.assertEqual(set(order),{model.End(sid,0) for sid in (4,5,6)})
        self.assertEqual(set(chains),set(zip(order,order[1:]+order[:1])))
        self.assertEqual(len(cycles),1)
        self.assertEqual(len(cycles[0]),3)
        self.assertEqual({end.section_id for pair in cycles[0] for end in pair},{1,2,3})


if __name__=='__main__':unittest.main()
