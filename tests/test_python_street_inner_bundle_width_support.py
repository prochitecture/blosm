"""Only an exposed plateau root can add an original member width transition."""
from dataclasses import replace
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import (Section, End, Key, Intersection, Network,
    Bundle, BundleMember, Traversal, Settings)
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import UpdateCoordinator, BuildContext
from tested_python_streets.routes import Route, RouteTopology
from tested_python_streets.intersections import Approach
from tested_python_streets.compound_cycles import _footprint
from tested_python_streets.regions import BoundaryLoop, SourceSpan, Region
from tested_python_streets.inner_bundle_width_support import required


def fixture(reverse=False, angle=0., shift=(0., 0.)):
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    rows=((1,(-40.,0.),(0.,0.),6.),(2,(0.,0.),(40.,0.),4.),
          (3,(-40.,-12.),(40.,-12.),4.),(4,(0.,0.),(0.,30.),4.),
          (5,(0.,0.),(0.,-30.),4.))
    sections=[]
    for sid,first,last,width in rows:
        points=np.array((first,last))@rotation.T+shift
        if reverse and sid in (1,2):points=points[::-1]
        sections.append(Section(sid,1 if sid in (1,2) else sid,
            tuple((x,y,0.) for x,y in points),width=width))
    ends=(End(1,0 if reverse else 1),End(2,1 if reverse else 0),End(4,0),End(5,0))
    site=Intersection(10,ends,radius_modes=('STYLE',)*4)
    bundle=Bundle(7,(BundleMember(tuple(Traversal(sid,not reverse) for sid in (1,2))),
                     BundleMember((Traversal(3,True),))))
    network=Network('original_inner_width_root',tuple(sections),(site,),(bundle,))
    patterns=PatternLibrary()
    products=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(network).products
    context=BuildContext(network,products,patterns)
    approaches=tuple(Approach(Route.build(context,RouteTopology(network).trace(end)),end.end,
        products[Key('profile',end.section_id)]) for end in ends)
    a=approaches[0]
    xy=np.vstack((a.point(np.array((0.,30.)),.5*a.width),
                  a.point(np.array((30.,0.)),-.5*a.width)))
    source=_footprint(Key('original_port',1),xy,internal_edges=(3,))
    source.layout.original_root_caps={3:ends[0]}
    spans=tuple(SourceSpan(source.region.owner,i,0.,.1 if i==3 else 1.) for i in range(4))
    loop=BoundaryLoop(source.region.xy[::-1],spans)
    planner=SimpleNamespace(context=context,layouts={site.key:SimpleNamespace(approaches=approaches)})
    return planner,source,loop,site,bundle


class InnerBundleWidthSupportTests(unittest.TestCase):
    def test_original_profiles_widths_radius_and_finite_endpoints_survive_transforms(self):
        for reverse,angle,shift in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            planner,source,loop,site,bundle=fixture(reverse,angle,shift)
            before=planner.context.network.to_dict()
            token,record=required(planner,(loop,),(source,),(),{4,5},(site.key,),Settings())
            self.assertEqual(token,(bundle.id,0,site.key))
            self.assertEqual(record.report['minimum_radius'],2.)
            self.assertLessEqual(record.report['curvature_bound'],.5*(1.+1e-8))
            for approach,station,side,t in zip((record.corner.first,record.corner.second),
                    record.stations,(-1.,1.),(0.,1.)):
                self.assertLess(station,approach.length)
                self.assertEqual(float(approach.profile.value(side)),.08)
                self.assertEqual(float(approach.profile.derivative(side)),0.)
                np.testing.assert_allclose(record.curve(t),approach.point(station,side*.5*approach.width),atol=1e-9)
            self.assertEqual(planner.context.network.to_dict(),before)
            self.assertIsNone(required(planner,(loop,),(source,),(token,),{4,5},(site.key,),Settings()))

    def test_outer_roots_strict_overrides_consumed_pairs_and_missing_ownership_do_not_qualify(self):
        planner,source,loop,site,bundle=fixture()
        strict=replace(site,radius_modes=('STRICT',)*4)
        for consumed,members in (({1,4,5},(site.key,)),({4,5},())):
            self.assertIsNone(required(planner,(loop,),(source,),(),consumed,members,Settings()))
        original=planner.context.network
        planner.context.network=replace(original,intersections=(strict,))
        self.assertIsNone(required(planner,(loop,),(source,),(),{4,5},(site.key,),Settings()))
        planner.context.network=original
        outer=replace(loop,spans=tuple(replace(span,first=.8,last=1.) if span.segment==3 else span for span in loop.spans))
        self.assertIsNone(required(planner,(outer,),(source,),(),{4,5},(site.key,),Settings()))
        self.assertIsNone(required(planner,(replace(loop,xy=loop.xy[::-1]),),(source,),(),{4,5},(site.key,),Settings()))
        invalid=replace(bundle,members=(replace(bundle.members[0],path=bundle.members[0].path[::-1]),bundle.members[1]))
        planner.context.network=replace(original,bundles=(invalid,))
        self.assertIsNone(required(planner,(loop,),(source,),(),{4,5},(site.key,),Settings()))

    def test_a_root_tag_requires_the_actual_original_cross_section(self):
        planner,source,loop,site,bundle=fixture()
        moved=replace(source,region=Region(source.region.owner,source.region.xy+np.array((.01,0.))))
        with self.assertRaisesRegex(ValueError,'actual original root cross-section'):
            required(planner,(loop,),(moved,),(),{4,5},(site.key,),Settings())

    def test_a_zero_length_internal_tag_is_not_an_exposed_hole_root(self):
        planner,source,loop,site,bundle=fixture()
        zero=replace(loop,spans=tuple(replace(span,last=span.first) if span.segment==3 else span
                                      for span in loop.spans))
        self.assertIsNone(required(planner,(zero,),(source,),(),{4,5},(site.key,),Settings()))


if __name__=='__main__':unittest.main()
