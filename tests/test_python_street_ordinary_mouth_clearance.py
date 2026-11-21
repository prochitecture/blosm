"""An acute branch keeps its actual root width beyond another complete mouth."""
import math
from dataclasses import replace
from types import SimpleNamespace
from unittest import mock
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section, End, Intersection, Network, Settings, Key
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown, ProfileProvider
from tested_python_streets.routes import Route, RouteTopology, route_key
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.ordinary_mouth_clearance import correct, prepare, singleton_groups, signature
from tested_python_streets.source_strip_containment import SourceStripContainmentError
from tested_python_streets.cycle_dimensional_retry import certify_layout
from tested_python_streets.finite_profile_search import FiniteProfileUnavailable
from tested_python_streets.compounds import CompoundPlanner


def fixture(*, reverse=False, angle=0., shift=(0.,0.), finite_length=None):
    # Renumbered local controls of an authored acute three-arm road, retaining
    # the curved branch which protrudes beyond the wide through-road mouth.
    controls=(
      ((0.,0.),(-20.879272,35.693695),(-97.411133,166.534073),(-139.452026,250.323257),
       (-168.495728,324.023926),(-195.788208,399.472),(-218.401978,486.663208),(-229.52832,536.03627)),
      ((0.,0.),(-3.769897,13.47052),(-17.672852,41.748779),(-19.592407,50.009094),
       (-19.707153,58.61412),(-16.819092,67.641495),(-12.297119,77.625839),(-2.687988,86.050583),
       (9.308472,90.745605),(20.814209,92.446274),(33.820557,90.539871),(45.487915,84.704224),(58.66394,75.829239)),
      ((18.010254,-29.859863),(0.,0.)))
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    sources=[];ends=[]
    for sid,(points,width,offset,root) in enumerate(zip(controls,(7.,3.5,10.5),(0.,0.,1.75),(0,0,1)),1):
        points=np.asarray(points,float)
        if sid==3 and finite_length is not None:
            points[0]=points[0]/np.linalg.norm(points[0])*finite_length
        if reverse:points=points[::-1];root=1-root;offset=-offset
        points=points@rotation.T+shift
        sources.append(Section(sid,sid,tuple((x,y,0.) for x,y in points),width=width,offset=offset,
                               interpolation=3,scope_ends=(1-root,)))
        ends.append(End(sid,root))
    site=Intersection(1,tuple(ends),10.,(),('STYLE',)*3)
    network=Network('ordinary-mouth-source',tuple(sources),(site,),settings=Settings(exaggeration=1.))
    products={Key('profile',s.id):Crown(.08,0.) for s in sources}
    context=BuildContext(network,products,PatternLibrary())
    for s in sources:products[Key('path',s.id)]=PathProvider.build(context,s,())
    topology=RouteTopology(network)
    for end in ends:products[route_key(end)]=Route.build(context,topology.trace(end))
    layout=IntersectionProvider.footprint(context,site)
    planner=SimpleNamespace(context=context,topology=topology,layouts={site.key:layout})
    return planner,site,layout


class NativeBundleMouthTests(unittest.TestCase):
    def test_accepts_actual_oriented_half_crowns_and_rejects_substituted_profile(self):
        from test_python_street_bundle_terminal_expansion import fixture as bundle_fixture
        from tested_python_streets.ordinary_mouth_clearance import _eligible
        from tested_python_streets.profiles import memberships

        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                setup, _, _ = bundle_fixture(reverse)
                network = setup.network
                owners = memberships(network)
                products = {}
                context = BuildContext(network, products, PatternLibrary())
                topology = RouteTopology(network)
                for section in network.sections:
                    products[Key('path', section.id)] = Path(section)
                    products[Key('profile', section.id)] = ProfileProvider.build(
                        section, owners.get(section.id), network.section_index)
                site = network.intersections[0]
                for end in site.ends:
                    products[route_key(end)] = Route.build(context, topology.trace(end))
                layout = IntersectionProvider.footprint(context, site)
                half, = (arm for arm in layout.approaches if arm.profile.half)
                self.assertTrue(_eligible(context, layout))

                # A consistent pair of substituted profiles can pass the ordinary
                # oriented-profile equality test. The Bundle still owns its role.
                key = Key('profile', half.source.id)
                products[key] = replace(products[key], half=-products[key].half)
                half.profile = products[key].oriented(half.end)
                with self.assertRaisesRegex(ValueError, 'membership-resolved source profile'):
                    _eligible(context, layout)


def returning_fixture(*,reverse=False,angle=0.,shift=(0.,0.),mode='STYLE'):
    # Renumbered original smooth returning-road controls. Both real Ends meet
    # the same three-arm site; only their finite endpoint prefixes are owned.
    controls=(np.array(((0.,0.),(17.5615234375,1.29632568359375),
        (36.065673828125,-2.47198486328125),(48.9251708984375,-12.809814453125),
        (13.493896484375,-25.67755126953125),(0.,0.))),
        np.array(((-9.2698974609375,-3.79864501953125),(0.,0.))))
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    sections=[]
    for sid,points in enumerate(controls,1):
        if reverse:points=points[::-1]
        points=points@rotation.T+shift
        sections.append(Section(sid,sid,tuple((x,y,0.) for x,y in points),width=2.,interpolation=3))
    ends=(End(1,0),End(2,1),End(1,1))
    if reverse:ends=tuple(End(e.section_id,1-e.end) for e in ends)
    site=Intersection(1,ends,1.,(1.,1.,1.),(mode,)*3)
    network=Network('returning-prefixes',tuple(sections),(site,),settings=Settings(exaggeration=1.))
    products={Key('profile',s.id):Crown(.08,0.) for s in sections}
    context=BuildContext(network,products,PatternLibrary())
    for s in sections:products[Key('path',s.id)]=PathProvider.build(context,s,(),returning=s.id==1)
    topology=RouteTopology(network)
    for end in ends:products[route_key(end)]=Route.build(context,topology.trace(end))
    layout=IntersectionProvider.footprint(context,site)
    return SimpleNamespace(context=context,topology=topology,layouts={site.key:layout}),site,layout


class OrdinaryMouthClearanceTests(unittest.TestCase):
    def test_returning_original_prefixes_keep_radius_and_certified_curve_phase(self):
        from tested_python_streets.ordinary_mouth_clearance import _group
        from tested_python_streets.sampling import AnalyticalBoundarySampler
        for reverse,angle,shift,mode in ((False,0.,(0.,0.),'STYLE'),
                (True,.713,(1200.,-1800.),'STRICT')):
            planner,site,raw=returning_fixture(reverse=reverse,angle=angle,shift=shift,mode=mode)
            original=planner.context.network.to_dict()
            with self.assertRaises(SourceStripContainmentError):
                certify_layout(planner.context,_group(site.key,raw),())
            result=correct(planner,site.key,raw)
            self.assertIs(result.corners,raw.corners)
            self.assertIs(result.approaches,raw.approaches)
            self.assertEqual(result.outward_stations,raw.outward_stations)
            self.assertEqual(result.mouths,raw.mouths)
            self.assertTrue(result.ordinary_complete_mouth_clearance['certified_same_curve_sampling'])
            self.assertEqual(result.ordinary_complete_mouth_clearance['candidates_tested'],0)
            self.assertTrue(all(c.radius==1. and c.radius_mode==mode for c in result.corners))
            g=_group(site.key,result)
            sampled=AnalyticalBoundarySampler(tolerance=planner.context.network.settings.curve_tolerance,
                spacing=planner.context.network.settings.spacing).build(g)
            proof=certify_layout(planner.context,sampled,())
            ranges=[row['interval'] for row in proof['sections'] if row['section']==1]
            self.assertEqual(len(ranges),2)
            self.assertLess(ranges[0][1],ranges[1][0])
            self.assertEqual(ranges[0][0],0.)
            self.assertEqual(ranges[1][1],planner.context[Key('path',1)].length)
            with self.assertRaises(SourceStripContainmentError):
                certify_layout(planner.context,sampled,{1})
            self.assertEqual(planner.context.network.to_dict(),original)
            self.assertIs(correct(planner,site.key,result),result)

    def test_complete_widths_advance_without_changing_original_contacts_or_sources(self):
        for reverse,angle,shift in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            planner,site,raw=fixture(reverse=reverse,angle=angle,shift=shift)
            before=planner.context.network.to_dict();stations=raw.outward_stations;mouths=dict(raw.mouths)
            fixed=tuple((c.center.copy() if c.center is not None else None,c.stations,c.radius) for c in raw.corners)
            result=correct(planner,site.key,raw)
            self.assertIsNot(result,raw)
            self.assertEqual(raw.outward_stations,stations);self.assertEqual(raw.mouths,mouths)
            self.assertEqual(planner.context.network.to_dict(),before)
            self.assertIs(result.approaches,raw.approaches);self.assertIs(result.corners,raw.corners)
            for old,c in zip(fixed,result.corners):
                if old[0] is not None:np.testing.assert_array_equal(old[0],c.center)
                self.assertEqual(old[1:],(c.stations,c.radius));self.assertEqual(c.radius,10.)
            for a in result.approaches:
                self.assertEqual(a.profile,planner.context[Key('profile',a.source.id)].oriented(a.end))
                self.assertIs(a.path.pieces[0].path,planner.context[Key('path',a.source.id)])
            wide=next(i for i,a in enumerate(raw.approaches) if a.source.id==3)
            self.assertGreater(result.outward_stations[wide],raw.outward_stations[wide])
            self.assertLess(result.outward_stations[wide],1.)
            groups=singleton_groups({site.key:result},())
            proof=certify_layout(planner.context,groups[0],())
            self.assertEqual({v['section'] for v in proof['sections']},{1,2,3})
            self.assertEqual(proof['maximum_error'],.002)
            self.assertTrue(signature(groups[0]))
            # A fresh check preserves an already qualified ordinary layout exactly.
            self.assertIs(correct(planner,site.key,result),result)

    def test_noop_success_and_typed_unrelated_or_resource_failures(self):
        planner,site,raw=fixture()
        result=correct(planner,site.key,raw)
        with mock.patch('tested_python_streets.compound_cycles._canonical_samples',
                        side_effect=AssertionError('No grid needed for qualified input')):
            self.assertIs(correct(planner,site.key,result),result)
        with mock.patch('tested_python_streets.cycle_dimensional_retry.certify_layout',side_effect=ValueError('allocation witness')):
            with self.assertRaisesRegex(ValueError,'allocation witness'):correct(planner,site.key,raw)
        with mock.patch('tested_python_streets.ordinary_mouth_clearance.segment_on_polygon_side',return_value=True):
            with self.assertRaises(SourceStripContainmentError):correct(planner,site.key,raw)
        with self.assertRaisesRegex(ValueError,'bounded'):correct(planner,site.key,raw,maximum_passes=17)

    def test_finite_scope_exhaustion_and_wrong_path_ownership_remain_fatal(self):
        planner,site,raw=fixture(finite_length=.35)
        with self.assertRaises(FiniteProfileUnavailable):correct(planner,site.key,raw)
        planner,site,raw=fixture()
        old=planner.context[Key('path',1)]
        planner.context.products[Key('path',1)]=Path(old.section)
        with self.assertRaisesRegex(ValueError,'Path ownership'):correct(planner,site.key,raw)

    def test_effective_layout_and_singleton_owner_preserve_original_incidence(self):
        planner,site,raw=fixture()
        layouts={site.key:raw}
        effective,changed=prepare(planner.context,(site,),layouts)
        self.assertIs(layouts[site.key],raw)
        groups=CompoundPlanner().build(effective)
        groups=(*groups,*singleton_groups(changed,groups))
        self.assertEqual(len(groups),1);group=groups[0]
        self.assertEqual(group.members,(site.key,));self.assertEqual(group.claims,())
        self.assertEqual(set(group.physical_layout.mouths),set(raw.mouths))
        for a,b in zip(group.physical_layout.approaches,raw.approaches):self.assertIs(a,b)

    def test_normal_backend_publishes_complete_profiles_and_coupled_checks(self):
        from tested_python_streets.backend import PythonMeshBackend
        from tested_python_streets.ownership import PLAN, surface_key
        planner,site,raw=fixture()
        network=planner.context.network;before=network.to_dict()
        backend=PythonMeshBackend(network.settings)
        update=backend.prepare(network)
        group=update.products[PLAN].by_site[site.key]
        product=update.products[surface_key(site.key)]
        self.assertEqual(group.members,(site.key,));self.assertEqual(group.claims,())
        self.assertTrue(product.pavement.metadata['ordinary_complete_mouth_clearance'])
        self.assertTrue(product.pavement.metadata['source_cycle_reconciliation'])
        for approach,station in zip(group.physical_layout.approaches,group.physical_layout.outward_stations):
            sid,end=list(approach.path.trims(station))[-1]
            for u in (0.,.5,1.):
                loop,at,source_station=product.port_points[sid,end,u]
                expected=approach.point(station,approach.offset+approach.width*(u-.5))
                np.testing.assert_allclose(product.rims[loop][at,:2],expected,atol=1e-10,rtol=0.)
                self.assertAlmostEqual(product.rims[loop][at,2],approach.profile.value(2*u-1),places=12)
        self.assertEqual(network.to_dict(),before)
        backend.commit(update)
        self.assertEqual(backend.prepare(network).rebuilt,())

    def test_prepared_neighbor_inventory_keeps_real_owner_and_opposite_trims(self):
        from tested_python_streets.compound_preparation import DeferredCompoundFootprint
        from tested_python_streets.intersections import Mouth, Approach
        from tested_python_streets.compounds import CompoundGroup
        planner,site,raw=fixture()
        network=planner.context.network
        wide=network.section_index[3]
        location=np.asarray(wide.points[0])
        source4=Section(4,4,(tuple(location),tuple(location+(10.,0.,0.))),width=3.5)
        source5=Section(5,5,(tuple(location+(10.,0.,0.)),tuple(location+(20.,0.,0.))),width=3.5)
        source6=Section(6,6,(tuple(location),tuple(location+(0.,-10.,0.))),width=3.5)
        source7=Section(7,7,(tuple(location+(10.,0.,0.)),tuple(location+(10.,-10.,0.))),width=3.5)
        neighbor=Intersection(2,(End(3,0),End(4,0),End(6,0)))
        other=Intersection(3,(End(4,1),End(5,0),End(7,0)))
        changed=replace(network,sections=tuple(replace(s,scope_ends=()) if s.id==3 else s for s in network.sections)
            +(source4,source5,source6,source7),intersections=(site,neighbor,other))
        planner.context.network=changed
        planner.context.products[Key('path',3)]=PathProvider.build(planner.context,changed.section_index[3],())
        topology=RouteTopology(changed)
        planner.context.products[route_key(End(3,1))]=Route.build(planner.context,topology.trace(End(3,1)))
        raw=IntersectionProvider.footprint(planner.context,site)
        raw_neighbors={neighbor.key:SimpleNamespace(mouths={(3,0):Mouth(34.5,8)}),
                       other.key:SimpleNamespace(mouths={})}
        layouts={site.key:raw,**raw_neighbors}
        opposite=Approach(Route.build(planner.context,topology.trace(End(3,0))),0,
                          planner.context[Key('profile',3)])
        prepared=CompoundGroup((neighbor.key,other.key),(),(),(),
            SimpleNamespace(approaches=(opposite,),outward_stations=(1.,)))
        effective,corrected=prepare(planner.context,changed.intersections,layouts,
            targets={site.key},prepared_groups=(prepared,))
        self.assertIn(site.key,corrected)
        self.assertIs(effective[neighbor.key],raw_neighbors[neighbor.key])
        self.assertIs(effective[other.key],raw_neighbors[other.key])
        # The prepared owner's actual cut, rather than the stale raw cut,
        # now closes the available interval. It must request BOTH real members.
        prepared=replace(prepared,physical_layout=SimpleNamespace(approaches=(opposite,),outward_stations=(34.5,)))
        effective,corrected=prepare(planner.context,changed.intersections,layouts,
            targets={site.key},prepared_groups=(prepared,))
        self.assertNotIn(site.key,corrected)
        self.assertIsInstance(effective[site.key],DeferredCompoundFootprint)
        self.assertEqual(effective[site.key].neighbors,(neighbor.key,other.key))

    def test_real_source_tasks_execute_the_certified_late_scope_composer(self):
        from test_python_street_finite_scope_compounds import fixture as scope_fixture, OwnershipOnly, prepare as reference
        from tested_python_streets.ownership import OwnershipPlan, PLAN, footprint_key
        from tested_python_streets.finite_scope_ports import groups as scope_groups
        from tested_python_streets.ordinary_mouth_clearance import finalize
        from tested_python_streets.updates import UpdateCoordinator
        captured={}
        class LateScopeStage(OwnershipOnly):
            def tasks(self,network):
                for task in super().tasks(network):
                    if task.key!=PLAN:
                        yield task;continue
                    def build(context):
                        layouts={j.key:context[footprint_key(j.key)] for j in context.network.intersections}
                        existing=scope_groups(context,layouts)
                        captured['original']=existing[0]
                        captured['loop']=existing[0].loops[0].xy.copy()
                        final=finalize(context,context.network.intersections,layouts,existing)
                        return OwnershipPlan(final,{site:g for g in final for site in g.members})
                    yield replace(task,build=build)
        network=scope_fixture();before=network.to_dict()
        update=UpdateCoordinator((PathProvider(),ProfileProvider(),LateScopeStage()),PatternLibrary()).prepare(network)
        group=update.products[PLAN].by_site[Key('intersection',1)]
        expected=reference(network).products[PLAN].by_site[Key('intersection',1)]
        self.assertEqual(group.members,(Key('intersection',1),Key('intersection',2)))
        self.assertEqual(group.claims,expected.claims)
        np.testing.assert_array_equal(group.loops[0].xy,expected.loops[0].xy)
        np.testing.assert_array_equal(captured['loop'],captured['original'].loops[0].xy)
        self.assertTrue(group.physical_layout.finite_scope_source_proof)
        self.assertEqual(network.to_dict(),before)

    def test_real_source_task_rejects_unsupported_owner_requirement_explicitly(self):
        from test_python_street_finite_scope_compounds import fixture as scope_fixture, OwnershipOnly
        from tested_python_streets.ownership import OwnershipPlan, PLAN, footprint_key
        from tested_python_streets.compound_preparation import DeferredCompoundFootprint
        from tested_python_streets.ordinary_mouth_clearance import finalize,PreparedOwnerClearanceRequired
        from tested_python_streets.updates import UpdateCoordinator
        from tested_python_streets.profiles import ProfileProvider
        class RequiredClosureStage(OwnershipOnly):
            def tasks(self,network):
                for task in super().tasks(network):
                    if task.key!=PLAN:
                        yield task;continue
                    def build(context):
                        accepted=OwnershipPlan.build(context,context.network.intersections)
                        layouts={j.key:context[footprint_key(j.key)] for j in context.network.intersections}
                        # Fault injection at the typed source-requirement boundary:
                        # source6 really connects this site to the composed owner2.
                        # No pavement, profile, alias or accepted product is injected.
                        site=Key('intersection',3);raw=layouts[site]
                        layouts[site]=DeferredCompoundFootprint(site,raw.approaches,
                            (Key('intersection',1),Key('intersection',2)),'explicit connected owner requirement')
                        return finalize(context,context.network.intersections,layouts,accepted.groups)
                    yield replace(task,build=build)
        network=scope_fixture();before=network.to_dict()
        with self.assertRaises(ValueError) as caught:
            UpdateCoordinator((PathProvider(),ProfileProvider(),RequiredClosureStage()),PatternLibrary()).prepare(network)
        error=caught.exception.__cause__
        self.assertIsInstance(error,PreparedOwnerClearanceRequired)
        self.assertEqual(error.owner_members,(Key('intersection',1),Key('intersection',2)))
        self.assertEqual(error.family,'composed_finite_scope')
        self.assertEqual(network.to_dict(),before)


if __name__=='__main__':unittest.main()
