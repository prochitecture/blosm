"""A typed acute root retry owns only its actual exposed source half."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,End,Key
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.compound_cycles import _footprint
from tested_python_streets.regions import BoundaryLoop,SourceSpan
from tested_python_streets.cycle_dimensional_retry import _exposed_facing_cap,eligible,required
from tested_python_streets.cycle_root_support import prepare_required
from tested_python_streets.model import Intersection,Settings
from tested_python_streets.regions import Region


def fixture():
    first=Approach(Path(Section(1,1,((0.,0.,0.),(40.,0.,0.)),width=3.5)),0,Crown(.08,0.))
    angle=np.deg2rad(30.)
    second=Approach(Path(Section(2,2,((0.,0.,0.),(40.*np.cos(angle),40.*np.sin(angle),0.)),width=7.)),0,Crown(.08,0.,-1))
    corner=CornerRadiusPolicy.solve(first,second,2.,'STYLE')
    xy=np.vstack((second.point(np.array((0.,40.)),3.5),second.point(np.array((40.,0.)),-3.5)))
    source=_footprint(Key('source_strip',2),xy,internal_edges=(3,))
    source.layout.original_root_caps={3:End(2,0)}
    spans=tuple(SourceSpan(source.region.owner,i,.999,1.) if i==3 else SourceSpan(source.region.owner,i,0.,1.) for i in range(4))
    return corner,source,BoundaryLoop(source.region.xy,spans)


class CycleRootHalfTests(unittest.TestCase):
    def test_half_profile_transition_requires_its_actual_zero_height_facing_root(self):
        a=Approach(Path(Section(1,1,((0.,0.,0.),(40.,0.,0.)),width=3.5)),0,Crown(.08,0.))
        b=Approach(Path(Section(2,2,((0.,0.,0.),(-40.,0.,0.)),width=7.)),0,Crown(.08,0.,-1))
        co=CornerRadiusPolicy.solve(a,b,2.,'STYLE');self.assertTrue(eligible(co))
        self.assertFalse(eligible(replace(co,first=b,second=a)))
        self.assertFalse(eligible(replace(co,radius_mode='STRICT')))
        xy=np.vstack((b.point(np.array((0.,40.)),3.5),b.point(np.array((40.,0.)),-3.5)))
        source=_footprint(Key('source_strip',2),xy,internal_edges=(3,))
        source.layout.original_root_caps={3:End(2,0)}
        spans=tuple(SourceSpan(source.region.owner,i,.8,1.) if i==3 else SourceSpan(source.region.owner,i,0.,1.) for i in range(4))
        loop=BoundaryLoop(source.region.xy,spans);site=Intersection(1,(End(1,0),End(2,0)))
        planner=SimpleNamespace(sites={site.key:site})
        pair=(End(1,0),End(2,0));candidates={pair:(site.key,co)}
        self.assertEqual(required(planner,(loop,),(source,),candidates,()),pair)
        self.assertIsNone(required(planner,(loop,),(source,),candidates,(pair,)))
        moved=replace(source,region=Region(source.region.owner,source.region.xy+.01*b.direction))
        self.assertIsNone(required(planner,(loop,),(moved,),candidates,()))

    def test_exposed_right_root_belongs_to_the_second_facing_corner_only(self):
        corner,source,loop=fixture()
        self.assertTrue(_exposed_facing_cap((loop,),(source,),corner))
        opposite=replace(corner,first=corner.second,second=corner.first)
        self.assertFalse(_exposed_facing_cap((loop,),(source,),opposite))

    def test_coordinate_agreement_cannot_replace_the_original_root_incidence(self):
        corner,source,loop=fixture()
        source.layout.original_root_caps={3:End(99,0)}
        self.assertFalse(_exposed_facing_cap((loop,),(source,),corner))
        source.layout.original_root_caps={}
        self.assertFalse(_exposed_facing_cap((loop,),(source,),corner))

    def test_a_root_tag_must_also_match_the_actual_declared_cross_section_geometry(self):
        corner,source,loop=fixture()
        third=Approach(Path(Section(3,3,((0.,0.,0.),(0.,-40.,0.)),width=7.)),0,Crown(.08,0.))
        site=Intersection(100,(End(1,0),End(2,0),End(3,0)))
        planner=SimpleNamespace(sites={site.key:site},layouts={site.key:SimpleNamespace(approaches=(corner.first,corner.second,third))})
        moved=replace(source,region=Region(source.region.owner,source.region.xy+.01*corner.second.direction))
        with self.assertRaisesRegex(ValueError,'actual original root cross-section'):
            prepare_required(planner,(loop,),(moved,),(),(),set(),(site.key,),Settings())


if __name__=='__main__':unittest.main()
