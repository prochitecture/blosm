"""Original ordinary mouth planes bound existing finite root-support curves."""
import math
from dataclasses import replace
from types import SimpleNamespace
from unittest import mock
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,End,Intersection,Network,Settings,Key
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import Crown
from tested_python_streets.providers import PathProvider
from tested_python_streets.routes import Route,RouteTopology,route_key
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.updates import BuildContext
from tested_python_streets.ordinary_mouth_clearance import correct
from tested_python_streets.ordinary_singleton_sources import build
from tested_python_streets.source_strip_containment import SourceStripContainmentError
from tested_python_streets.cyclic_root_support import prepare
from tested_python_streets.compound_preparation import DeferredCompoundPlanner
from tested_python_streets.compound_cycles import source_cycle_group
from tested_python_streets.compounds import rim_label
from tested_python_streets.continuations import ContinuationReachError


def fixture(reverse=False,angle=0.,shift=(0.,0.),mode='STYLE'):
    # Renumbered actual three-arm source controls. The seven-metre branch's
    # right root protrudes beyond its circular support beside two five-metre arms.
    controls=(
      ((-34.99005126953125,13.373291015625),(-20.8743896484375,6.85955810546875),(0.,0.)),
      ((-42.6630859375,-25.4429931640625),(-21.88720703125,-9.11468505859375),
       (-3.53118896484375,-1.46905517578125),(0.,0.)),
      ((0.,0.),(11.20745849609375,-.69134521484375),(23.51171875,-.302978515625),
       (37.84619140625,2.9683837890625),(49.63751220703125,7.3309326171875)))
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    sections=[]
    for sid,(points,width) in enumerate(zip(controls,(5.,7.,5.)),1):
        points=np.asarray(points)
        if reverse:points=points[::-1]
        points=points@rotation.T+shift
        sections.append(Section(sid,sid,tuple((x,y,0.) for x,y in points),width=width,interpolation=3))
    ends=(End(1,1),End(2,1),End(3,0))
    if reverse:ends=tuple(End(e.section_id,1-e.end) for e in ends)
    site=Intersection(1,ends,3.,(3.,2.,3.),(mode,)*3)
    n=Network('ordinary-root-bounds',tuple(sections),(site,),settings=Settings(exaggeration=1.))
    p={Key('profile',s.id):Crown(.08,0.) for s in sections};c=BuildContext(n,p,PatternLibrary())
    for s in sections:p[Key('path',s.id)]=PathProvider.build(c,s,())
    topology=RouteTopology(n)
    for end in ends:p[route_key(end)]=Route.build(c,topology.trace(end))
    raw=IntersectionProvider.footprint(c,site)
    return SimpleNamespace(context=c,topology=topology,layouts={site.key:raw}),site,raw


def arguments(raw,site):
    donor=next(c for c in raw.corners if (c.first.source.id,c.second.source.id)==(3,1))
    return donor,dict(site=site,original_approaches=raw.approaches,
        support_corners=tuple(c for c in raw.corners if c is not donor),
        required_roots=((next(a for a in raw.approaches if a.source.id==2),1),))


class OrdinaryRootBoundTests(unittest.TestCase):
    def test_strict_source_and_an_untyped_three_arm_request_do_not_enable_support(self):
        planner,site,raw=fixture(mode='STRICT')
        with self.assertRaises(SourceStripContainmentError) as caught:
            correct(planner,site.key,raw)
        with self.assertRaisesRegex(ValueError,'internal fillet preparation seam'):
            build(planner.context,site.key,raw,planner.layouts,(site,),caught.exception)
        self.assertTrue(all(c.radius_mode=='STRICT' for c in raw.corners))
        planner,site,raw=fixture()
        ordinary=DeferredCompoundPlanner(planner.context,(site,),planner.layouts)
        ordinary.ordinary_singleton_root_bounds=object()
        with self.assertRaisesRegex(ValueError,'typed original singleton context'):
            source_cycle_group(ordinary,(site.key,),set(),set(),{},set())

    def test_actual_three_arm_root_is_opt_in_and_keeps_reversed_transformed_mouths(self):
        for reverse,angle,shift in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            planner,site,raw=fixture(reverse,angle,shift);before=planner.context.network.to_dict()
            with self.assertRaises(SourceStripContainmentError) as caught:
                correct(planner,site.key,raw)
            ordinary=DeferredCompoundPlanner(planner.context,(site,),planner.layouts)
            with self.assertRaisesRegex(ValueError,'internal fillet preparation seam'):
                source_cycle_group(ordinary,(site.key,),set(),set(),{},set())
            g=build(planner.context,site.key,raw,planner.layouts,(site,),caught.exception)
            self.assertEqual(g.members,(site.key,));self.assertEqual(g.claims,())
            self.assertEqual(g.physical_layout.mouths,raw.mouths)
            self.assertEqual(g.physical_layout.outward_stations,raw.outward_stations)
            records=g.physical_layout.source_cycle_root_support_envelopes
            self.assertEqual(len(records),1);record=records[0]
            self.assertEqual(record['preferred_minimum_radius'],3.)
            self.assertLessEqual(record['curvature_bound'],1./3.)
            bounds=record['original_mouth_contact_bounds']
            self.assertLessEqual(bounds['anchored_searches'],2)
            self.assertTrue(bounds['both_original_lead_bounds_preserved'])
            source={s.region.owner:s for s in g.footprints}
            self.assertFalse(any(rim_label(span,source).kind=='internal' for loop in g.loops for span in loop.spans))
            self.assertEqual(planner.context.network.to_dict(),before)

    def test_default_none_parity_and_original_lead_bound_exhaustion(self):
        planner,site,raw=fixture();donor,kw=arguments(raw,site)
        settings=replace(planner.context.network.settings,curve_tolerance=.00025)
        default=prepare(Key('patch',1),donor,settings,**kw)
        explicit=prepare(Key('patch',1),donor,settings,maximum_stations=None,**kw)
        self.assertEqual(default.report,explicit.report)
        self.assertEqual(default.stations,explicit.stations)
        np.testing.assert_array_equal(default.source.region.xy,explicit.source.region.xy)
        np.testing.assert_array_equal(default.curve(np.linspace(0.,1.,101)),explicit.curve(np.linspace(0.,1.,101)))
        self.assertNotIn('original_mouth_contact_bounds',default.report)
        from tested_python_streets import cyclic_root_support
        search=cyclic_root_support._search_steps;anchors=[]
        def observed(*args,**kwargs):
            if kwargs.get('anchor') is not None:anchors.append(kwargs['anchor'])
            return search(*args,**kwargs)
        with mock.patch.object(cyclic_root_support,'_search_steps',side_effect=observed):
            with self.assertRaisesRegex(ContinuationReachError,'original mouth lead bounds'):
                prepare(Key('patch',1),donor,settings,maximum_stations=(.1,.1),**kw)
        self.assertEqual(anchors,[0,1])

    def test_unrelated_errors_strict_donors_and_wrong_source_bounds_are_not_recovered(self):
        planner,site,raw=fixture();donor,kw=arguments(raw,site)
        settings=replace(planner.context.network.settings,curve_tolerance=.00025)
        with self.assertRaisesRegex(ValueError,'STYLE'):
            prepare(Key('patch',1),replace(donor,radius_mode='STRICT'),settings,maximum_stations=(3.,3.),**kw)
        for bounds in ((float('nan'),3.),(donor.first.length+1.,3.),(3.,)):
            with self.assertRaisesRegex(ValueError,'finite original donor stations'):
                prepare(Key('patch',1),donor,settings,maximum_stations=bounds,**kw)
        baseline=prepare(Key('patch',1),donor,settings,**kw)
        with mock.patch('tested_python_streets.cyclic_root_support._search_steps',
                side_effect=(_python_street_package.prepared_steps(baseline),ValueError('allocation witness'))):
            with self.assertRaisesRegex(ValueError,'allocation witness'):
                prepare(Key('patch',1),donor,settings,maximum_stations=(.1,.1),**kw)


if __name__=='__main__':unittest.main()
