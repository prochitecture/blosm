"""Compound planning must repair partial ports even without a deferred site."""
import unittest
from dataclasses import replace
import numpy as np
import _python_street_package
from _python_street_ordinary_port_fixture import fixture,returning_fixture
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.updates import UpdateCoordinator,BuildContext
from tested_python_streets.source_support import ClippedApproachSupport
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.compounds import CompoundPlanner
from tested_python_streets.compound_preparation import DeferredCompoundPlanner,DeferredCompoundFootprint
from tested_python_streets.compound_ports import CompoundPortPlanner
from tested_python_streets.sampling import AnalyticalBoundarySampler
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.compound_port_fragments import normalize


class OrdinaryPortTests(unittest.TestCase):
    def test_microscopic_internal_caps_refine_curbs_without_discarding_external_profiles(self):
        n=fixture();before=n.to_dict();patterns=PatternLibrary()
        products=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(n).products
        context=BuildContext(n,products,patterns);layouts={}
        for site in n.intersections:
            ClippedApproachSupport._approaches(context,site)
            layouts[site.key]=IntersectionProvider.footprint(context,site)
        group=CompoundPlanner().build(layouts)[0]
        self.assertTrue(CompoundPortPlanner.incomplete_ports(group))
        normalized=normalize(context,group)
        self.assertIsNot(normalized,group)
        self.assertEqual(normalized.claims,group.claims)
        self.assertFalse(CompoundPortPlanner.incomplete_ports(normalized))
        report=normalized.physical_layout.internal_profile_fragments
        self.assertEqual(len(report['fragments']),2)
        self.assertTrue(all(r['fragment_length']<.001 for r in report['fragments']))
        self.assertEqual(report['complete_source_width']['holes'],1)
        self.assertEqual(len(normalized.physical_layout.approaches),4)
        canonical=normalized.physical_layout.source_curb_owner_canonicalization
        self.assertTrue(canonical)
        self.assertTrue(all(row['endpoint_error']<=1e-7 for row in canonical))
        inventory={f.region.owner:f.layout for f in normalized.footprints}
        for loop in normalized.loops:
            boundary=BoundaryConditions().build(normalized,loop,inventory)
            for port in boundary.ports:
                layout=inventory[port.owner];a=layout.approaches[port.approach_index]
                expected=a.point(layout.outward_stations[port.approach_index],
                    a.offset+a.width*(np.asarray(port.outward_u)-.5))
                np.testing.assert_allclose(boundary.xy[list(port.vertices)],expected,atol=1e-10,rtol=0.)
        self.assertIs(normalize(context,normalized),normalized)
        # Without an actual internal ownership claim, the same clipped profile
        # is external and must never be discarded merely because it is small.
        external=replace(group,claims=())
        self.assertIs(normalize(context,external),external)
        self.assertEqual(n.to_dict(),before)

    def test_returning_source_cannot_be_replaced_by_a_shorter_displaced_circle(self):
        n=returning_fixture();before=n.to_dict();patterns=PatternLibrary()
        products=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(n).products
        context=BuildContext(n,products,patterns);layouts={}
        for site in n.intersections:
            ClippedApproachSupport._approaches(context,site)
            layouts[site.key]=IntersectionProvider.footprint(context,site)
        ordinary=CompoundPlanner().build(layouts)
        self.assertTrue(any(CompoundPortPlanner.incomplete_ports(group) for group in ordinary))
        resolved=DeferredCompoundPlanner(context,n.intersections,layouts).build(ordinary)
        self.assertEqual(len(resolved),1);group=resolved[0]
        self.assertFalse(CompoundPortPlanner.incomplete_ports(group))
        self.assertTrue(any(getattr(curve,'pieces',()) for curve in group.physical_layout.boundaries))
        consumed={claim.section_id for claim in group.claims}
        self.assertEqual(len(consumed),1)
        CompoundPortPlanner(context,n.intersections,layouts).check_support(consumed,group.loops[0])
        sampled=AnalyticalBoundarySampler().build(group)
        inventory={f.region.owner:f.layout for f in sampled.footprints}
        boundaries=[BoundaryConditions().build(sampled,loop,inventory) for loop in sampled.loops]
        self.assertEqual(sum(len(b.ports) for b in boundaries),4)
        self.assertEqual(n.to_dict(),before)

    def test_valid_individual_footprints_do_not_authorize_partial_compound_profiles(self):
        n=fixture();before=n.to_dict();patterns=PatternLibrary()
        products=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(n).products
        context=BuildContext(n,products,patterns);layouts={}
        for site in n.intersections:
            ClippedApproachSupport._approaches(context,site)
            layouts[site.key]=IntersectionProvider.footprint(context,site)
        self.assertFalse(any(isinstance(layout,DeferredCompoundFootprint) for layout in layouts.values()))
        ordinary=CompoundPlanner().build(layouts)
        self.assertTrue(any(CompoundPortPlanner.incomplete_ports(group) for group in ordinary))
        resolved=DeferredCompoundPlanner(context,n.intersections,layouts).build(ordinary)
        self.assertEqual(len(resolved),1);group=resolved[0]
        self.assertIsNotNone(group.physical_layout)
        self.assertFalse(CompoundPortPlanner.incomplete_ports(group))
        self.assertEqual(sum(loop.hole for loop in group.loops),1)
        sampled=AnalyticalBoundarySampler().build(group)
        inventory={f.region.owner:f.layout for f in sampled.footprints}
        boundaries=[BoundaryConditions().build(sampled,loop,inventory) for loop in sampled.loops]
        self.assertEqual(sum(len(b.ports) for b in boundaries),4)
        self.assertTrue(all(abs(port.outward_u[0]-1.)<1e-8 and abs(port.outward_u[-1])<1e-8
            for boundary in boundaries for port in boundary.ports))
        self.assertEqual(n.to_dict(),before)


if __name__=='__main__':unittest.main()
