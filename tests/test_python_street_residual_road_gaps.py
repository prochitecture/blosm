"""Tiny generated corner gaps become pavement without moving other boundaries."""
from dataclasses import replace
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.corner_patch import corner_patch
from tested_python_streets.compounds import CompoundGroup
from tested_python_streets.regions import RegionUnion, signed_area
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.residual_road_gaps import reconcile, signature, ResidualRoadGapFill


def fixture(*, distance=1.9, reverse=False, angle=0., shift=(0.,0.), mode='STYLE', profile_half=0):
    sections=[];sites=[];sources=[];products={};radius=2.
    for i in range(3):
        direction=angle+i*math.tau/3.
        beta=direction-math.pi/4.
        rotation=np.array(((math.cos(beta),-math.sin(beta)),(math.sin(beta),math.cos(beta))))
        center=distance*np.array((math.cos(direction),math.sin(direction)))+shift
        origin=center+np.array((-(radius+1.),-(radius+1.)))@rotation.T
        arms=[]
        for j,delta in enumerate(((30.,0.),(0.,30.))):
            sid=2*i+j+1;points=np.array((origin,origin+np.asarray(delta)@rotation.T))
            xyz=np.c_[points,np.zeros(2)]
            if reverse:xyz=xyz[::-1]
            section=model.Section(sid,sid,tuple(map(tuple,xyz)),width=2.,interpolation=0)
            path=Path(section,.005);profile=Crown(.08,0.,profile_half)
            sections.append(section);products[model.Key('path',sid)]=path;products[model.Key('profile',sid)]=profile
            arms.append(Approach(path,int(reverse),profile))
        corner=CornerRadiusPolicy.solve(*arms,radius,mode)
        source=corner_patch(model.Key('cycle_corner',i),corner,.2,.00025)
        source.layout.original_junction=model.Key('intersection',i)
        sources.append(source)
        sites.append(model.Intersection(i,tuple(model.End(a.source.id,a.end) for a in arms),
                                       radius=radius,radius_modes=(mode,mode)))
    network=model.Network('three-original-corner-gaps',tuple(sections),tuple(sites))
    context=BuildContext(network,products,PatternLibrary())
    physical=SimpleNamespace(approaches=(),outward_stations=(),mouths={},
        corners=tuple(source.layout.corners[0] for source in sources))
    group=CompoundGroup(tuple(site.key for site in sites),(),tuple(sources),
        RegionUnion().build(source.region for source in sources),physical)
    return SimpleNamespace(context=context),group


class ResidualRoadGapTests(unittest.TestCase):
    def test_actual_original_corners_fill_only_the_tiny_generated_hole(self):
        planner,group=fixture();before=planner.context.network.to_dict()
        self.assertEqual(len(group.loops),2)
        result=reconcile(planner,group)
        self.assertEqual(len(result.loops),1)
        self.assertIs(result.loops[0],group.loops[0])
        self.assertEqual(result.claims,group.claims);self.assertIs(result.members,group.members)
        self.assertIs(result.physical_layout.corners,group.physical_layout.corners)
        self.assertIs(result.physical_layout.mouths,group.physical_layout.mouths)
        self.assertEqual(result.footprints[:-1],group.footprints)
        addition=result.footprints[-1]
        self.assertIsInstance(addition.layout.generated_junction_area,ResidualRoadGapFill)
        np.testing.assert_array_equal(addition.region.xy,group.loops[1].xy[::-1])
        independent=RegionUnion().build(source.region for source in result.area_footprints)
        self.assertEqual(len(independent),1)
        self.assertAlmostEqual(abs(signed_area(independent[0].xy)),abs(signed_area(result.loops[0].xy)),places=9)
        self.assertIs(reconcile(planner,result),result)
        self.assertEqual(planner.context.network.to_dict(),before)
        self.assertIsNotNone(signature(result))

    def test_rotation_translation_and_storage_direction_preserve_eligibility(self):
        areas=[]
        for reverse,angle,shift in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            planner,group=fixture(reverse=reverse,angle=angle,shift=shift)
            result=reconcile(planner,group)
            self.assertEqual(len(result.loops),1)
            areas.append(result.physical_layout.residual_road_gap_reconciliation['added_pavement_area'])
        self.assertAlmostEqual(*areas,places=9)

    def test_large_holes_and_explicit_corner_radii_keep_original_geometry(self):
        for parameters in (dict(distance=1.3),dict(mode='STRICT')):
            planner,group=fixture(**parameters)
            self.assertGreater(len(group.loops),1)
            self.assertIs(reconcile(planner,group),group)

    def test_native_bank_labels_are_never_treated_as_generated_corner_gaps(self):
        planner,group=fixture()
        changed=replace(group,footprints=tuple(replace(source,
            region=replace(source.region,owner=model.Key('cycle_strip',source.region.owner.id)))
            for source in group.footprints))
        # Re-union supplies matching real strip-owner spans. Shape alone is
        # never sufficient to erase that native-road boundary.
        changed=replace(changed,loops=RegionUnion().build(source.region for source in changed.footprints))
        self.assertIs(reconcile(planner,changed),changed)

    def test_bundle_members_and_true_inner_plateaus_keep_their_holes(self):
        planner,group=fixture(profile_half=1)
        self.assertIs(reconcile(planner,group),group)
        planner,group=fixture()
        bundle=model.Bundle(10,tuple(model.BundleMember((model.Traversal(sid),)) for sid in (1,3)))
        planner.context.network=replace(planner.context.network,bundles=(bundle,))
        self.assertIs(reconcile(planner,group),group)

    def test_work_allocation_is_explicit_and_bounded(self):
        planner,group=fixture()
        for count in (0,17,True):
            with self.assertRaisesRegex(ValueError,'finite work allowance'):
                reconcile(planner,group,maximum_fills=count)


if __name__=='__main__':unittest.main()
