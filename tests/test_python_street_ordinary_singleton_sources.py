"""A real ordinary width step retains its original finite source prefixes."""
from copy import copy
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
from tested_python_streets.ordinary_mouth_clearance import correct,finalize
from tested_python_streets.ordinary_singleton_sources import build,OrdinarySourceSingletonUnavailable
from tested_python_streets.source_strip_containment import SourceStripContainmentError


def fixture(mode='STYLE'):
    controls=(
      ((0.,0.),(19.07122802734375,5.69775390625),(208.25439453125,15.97979736328125),
       (235.59246826171875,13.68487548828125),(262.40936279296875,8.86322021484375),
       (314.51092529296875,-4.052490234375),(327.2507019042969,-2.65045166015625)),
      ((0.,0.),(8.98382568359375,7.156982421875),(16.871337890625,17.94305419921875),
       (20.7220458984375,26.3250732421875),(22.75274658203125,34.707275390625),
       (22.77703857421875,44.6146240234375),(21.35479736328125,55.9693603515625),
       (16.8447265625,65.5321044921875),(11.2457275390625,74.0152587890625)),
      ((-49.63751220703125,-7.3309326171875),(-38.4300537109375,-8.02227783203125),
       (-26.12579345703125,-7.6339111328125),(-11.79132080078125,-4.362548828125),(0.,0.)))
    sections=tuple(Section(i,i,tuple((x,y,0.) for x,y in points),width=width,interpolation=3)
        for i,(points,width) in enumerate(zip(controls,(3.5,5.,5.)),1))
    site=Intersection(1,(End(1,0),End(2,0),End(3,1)),10.,(10.,3.,3.),(mode,)*3)
    n=Network('ordinary-singleton-source',sections,(site,),settings=Settings(exaggeration=1.))
    p={Key('profile',s.id):Crown(.08,0.) for s in sections};c=BuildContext(n,p,PatternLibrary())
    for s in sections:p[Key('path',s.id)]=PathProvider.build(c,s,())
    topology=RouteTopology(n)
    for end in site.ends:p[route_key(end)]=Route.build(c,topology.trace(end))
    raw=IntersectionProvider.footprint(c,site)
    planner=SimpleNamespace(context=c,topology=topology,layouts={site.key:raw})
    return planner,site,raw


class OrdinarySingletonSourceTests(unittest.TestCase):
    def test_normal_coordinator_keeps_the_real_single_owner_and_source_operands(self):
        from tested_python_streets.backend import PythonMeshBackend
        from tested_python_streets.ownership import PLAN,footprint_key,surface_key
        planner,site,raw=fixture();n=planner.context.network;before=n.to_dict()
        with self.assertRaises(SourceStripContainmentError):correct(planner,site.key,raw)
        backend=PythonMeshBackend(n.settings);update=backend.prepare(n)
        group=update.products[PLAN].by_site[site.key]
        original=update.products[footprint_key(site.key)]
        self.assertEqual(group.members,(site.key,));self.assertEqual(group.claims,())
        self.assertGreater(len(group.area_footprints),1)
        self.assertEqual(group.physical_layout.mouths,original.mouths)
        old={(a.source.id,a.end):(a,s) for a,s in zip(original.approaches,original.outward_stations)}
        for a,s in zip(group.physical_layout.approaches,group.physical_layout.outward_stations):
            prior,t=old[a.source.id,a.end]
            self.assertIs(a.path,prior.path);self.assertIs(a.profile,prior.profile);self.assertEqual(s,t)
        product=update.products[surface_key(site.key)]
        self.assertLessEqual(product.pavement.mesh.validate()['maximum_edge_ratio'],8.)
        record=product.pavement.metadata['ordinary_complete_mouth_clearance']['corrections'][0]
        self.assertEqual(record['correction'],'FINITE_ORIGINAL_SINGLE_JUNCTION_SOURCE_UNION')
        self.assertEqual(record['consumed_sections'],[])
        self.assertTrue(product.pavement.metadata['source_cycle_reconciliation']['exact_boundary'])
        self.assertEqual(n.to_dict(),before)
        backend.commit(update);self.assertEqual(backend.prepare(n).rebuilt,())

    def test_strict_radius_inventory_and_resource_failures_stay_fatal(self):
        planner,site,raw=fixture('STRICT')
        with self.assertRaises(SourceStripContainmentError) as captured:
            correct(planner,site.key,raw)
        failure=captured.exception
        g=build(planner.context,site.key,raw,planner.layouts,(site,),failure)
        self.assertTrue(all(c.radius_mode=='STRICT' for c in g.physical_layout.corners))
        invalid=copy(g.physical_layout)
        invalid.corners=(replace(invalid.corners[0],radius=invalid.corners[0].radius*.8),*invalid.corners[1:])
        with mock.patch('tested_python_streets.compound_cycles.source_cycle_group_steps',
                        return_value=_python_street_package.prepared_steps(replace(g,physical_layout=invalid))):
            with self.assertRaisesRegex(OrdinarySourceSingletonUnavailable,'STRICT'):
                build(planner.context,site.key,raw,planner.layouts,(site,),failure)
        invalid=copy(g.physical_layout)
        invalid.outward_stations=invalid.outward_stations[:-1]
        with mock.patch('tested_python_streets.compound_cycles.source_cycle_group_steps',
                        return_value=_python_street_package.prepared_steps(replace(g,physical_layout=invalid))):
            with self.assertRaisesRegex(OrdinarySourceSingletonUnavailable,'End station'):
                build(planner.context,site.key,raw,planner.layouts,(site,),failure)
        with mock.patch('tested_python_streets.compound_cycles.source_cycle_group_steps',side_effect=ValueError('allocation witness')):
            with self.assertRaisesRegex(ValueError,'allocation witness'):
                build(planner.context,site.key,raw,planner.layouts,(site,),failure)
        with self.assertRaisesRegex(OrdinarySourceSingletonUnavailable,'real complete-prefix'):
            build(planner.context,site.key,raw,planner.layouts,(site,),ValueError('unrelated witness'))
        with self.assertRaisesRegex(OrdinarySourceSingletonUnavailable,'failing original Section'):
            build(planner.context,site.key,raw,planner.layouts,(site,),
                  SourceStripContainmentError('unowned source witness',section_id=99))
        with mock.patch('tested_python_streets.compound_cycles.source_cycle_group_steps',
                        return_value=_python_street_package.prepared_steps(g)),\
                mock.patch('tested_python_streets.cycle_dimensional_retry.certify_layout',
                    side_effect=SourceStripContainmentError('fresh width proof failed',section_id=3)):
            with self.assertRaisesRegex(SourceStripContainmentError,'fresh width proof failed'):
                build(planner.context,site.key,raw,planner.layouts,(site,),failure)

    def test_an_already_valid_layout_never_tries_the_source_union(self):
        from test_python_street_ordinary_mouth_clearance import fixture as mouth_fixture
        planner,site,raw=mouth_fixture();good=correct(planner,site.key,raw)
        with mock.patch('tested_python_streets.ordinary_singleton_sources.build_steps',
                side_effect=AssertionError('A successful ordinary layout must keep its path')):
            self.assertEqual(finalize(planner.context,(site,),{site.key:good},()),())


if __name__=='__main__':unittest.main()
