"""Authored curved-fork sources, with no saved scenes or diagnostic pickles."""
from dataclasses import replace
from types import SimpleNamespace
from copy import copy
from unittest import mock
import math
import unittest
import numpy as np
import _python_street_package

from tested_python_streets.model import Section,End,Intersection,Network,Settings,Key
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import Crown
from tested_python_streets.providers import PathProvider
from tested_python_streets.routes import Route,RouteTopology,route_key
from tested_python_streets.updates import BuildContext
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.ordinary_mouth_clearance import correct
from tested_python_streets.ordinary_singleton_sources import OriginalSingletonRootBounds,_ends
from tested_python_streets.source_strip_containment import SourceStripContainmentError
from tested_python_streets.compound_preparation import DeferredCompoundPlanner
from tested_python_streets import compound_cycles,cycle_root_support
from tested_python_streets.source_root_templates import prepare,OriginalRootReference,OriginalRootTemplateUnavailable,NativeRootEnvelope,construction_signature
from tested_python_streets.source_root_ports import OriginalRootPort,OriginalRootPortUnavailable
from tested_python_streets.source_root_witnesses import StagedRecipientUnavailable,admission
from tested_python_streets.generated_junction_curbs import is_generated_junction_area
from tested_python_streets.compounds import SampledFootprint,RimSegment
from tested_python_streets.regions import Region


def fixture(reverse=False,angle=0.,shift=(0.,0.)):
    controls=(
        ((70.5457763671875,19.67291259765625),(55.2774658203125,13.29608154296875),
         (41.47955322265625,10.02484130859375),(9.75360107421875,3.17138671875),(0.,0.)),
        ((0.,0.),(-25.846435546875,-9.3809814453125),(-29.564453125,-11.7850341796875)),
        ((-2.98968505859375,-24.47882080078125),(-7.14990234375,-17.8770751953125),
         (-7.85693359375,-11.54290771484375),(-5.997314453125,-5.7099609375),(0.,0.)))
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    sections=[]
    for sid,(xy,width) in enumerate(zip(controls,(7.,10.5,7.)),1):
        xy=np.asarray(xy)
        if reverse:xy=xy[::-1]
        xy=xy@rotation.T+shift
        sections.append(Section(sid,sid,tuple((x,y,0.) for x,y in xy),width=width,rise=.08,interpolation=3))
    ends=(End(1,1),End(2,0),End(3,1))
    if reverse:ends=tuple(End(e.section_id,1-e.end) for e in ends)
    site=Intersection(1,ends,2.,(2.,2.,2.),('STYLE',)*3)
    network=Network('curved-source-root',tuple(sections),(site,),settings=Settings(exaggeration=10.))
    products={Key('profile',s.id):Crown(.08,0.) for s in sections}
    context=BuildContext(network,products,PatternLibrary(),{},
        {Key(kind,s.id):f'{kind}:{s.id}:authored' for s in sections for kind in ('path','profile')})
    context.source_digests[site.key]='authored-junction-radius-2'
    for s in sections:products[Key('path',s.id)]=PathProvider.build(context,s,())
    topology=RouteTopology(network)
    for end in ends:products[route_key(end)]=Route.build(context,topology.trace(end))
    raw=IntersectionProvider.footprint(context,site);layouts={site.key:raw}
    try:correct(SimpleNamespace(context=context,topology=topology,layouts=layouts),site.key,raw)
    except SourceStripContainmentError as failure:original_failure=failure
    else:raise AssertionError('Authored fork must expose the wider original source root')
    planner=DeferredCompoundPlanner(context,(site,),layouts)
    planner.ordinary_singleton_root_bounds=OriginalSingletonRootBounds(site.key,raw,original_failure)
    pair=(ends[0],ends[1]);captured={}
    class Captured(RuntimeError):pass
    def stop(*args,**kwargs):
        captured['group']=kwargs['source_group'];raise Captured()
    try:
        with mock.patch.object(cycle_root_support,'prepare_required',stop):
            compound_cycles.source_cycle_group(planner,(site.key,),set(),set(),{},set(),retry_dimensional=(pair,))
    except Captured:pass
    pre=captured['group']
    donor=next(c for c in pre.physical_layout.corners if c.first.source.id==3 and c.second.source.id==1)
    source=next(s for s in pre.area_footprints if any(c is donor for c in getattr(s.layout,'corners',())))
    required=next(a for a in pre.physical_layout.approaches if a.source.id==2)
    settings=replace(network.settings,curve_tolerance=min(.001,.05*network.settings.curve_tolerance))
    bounds=tuple(planner.ordinary_singleton_root_bounds.limits()[End(a.source.id,a.end)] for a in (donor.first,donor.second))
    return SimpleNamespace(context=context,planner=planner,site=site,raw=raw,pre=pre,donor=donor,owner=source.region.owner,
        required=required,settings=settings,bounds=bounds,dimensional_pair=pair)


def arguments(f):
    return dict(site=f.site,original_approaches=f.raw.approaches,required_roots=((f.required,-1),),
        maximum_stations=f.bounds,original_bounds=f.planner.ordinary_singleton_root_bounds)


class SourceRootTemplateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.f=fixture()

    def test_reverse_route_coincidence_needs_a_retained_circle(self):
        f=fixture(reverse=True)
        reference=OriginalRootReference.build(f.context,f.pre,f.site,f.owner,f.donor,
            f.raw.approaches,f.bounds,f.planner.ordinary_singleton_root_bounds)
        port=next(p for p in reference.ports if p.approach.source.id==2)
        station=.25
        neighbor=float(np.nextafter(station,np.inf))
        self.assertEqual(port.approach.path.trims(station),port.approach.path.trims(neighbor))
        with self.assertRaisesRegex(OriginalRootTemplateUnavailable,'original retained circle'):
            reference.contact_grid(f.context,port,np.array([station,neighbor]),station)
        distinct=np.array([station,station+1e-6])
        grid,replaced=reference.contact_grid(f.context,port,distinct,station)
        np.testing.assert_array_equal(grid,distinct)
        self.assertEqual(replaced,())

    def test_authored_curved_fork_reverse_and_rigid_keep_real_planes_and_full_sources(self):
        for reverse,angle,shift in ((False,0.,(0.,0.)),(True,0.,(0.,0.)),(False,.713,(1200.,-1800.))):
            with self.subTest(reverse=reverse,angle=angle):
                f=self.f if not reverse and angle==0. else fixture(reverse,angle,shift)
                before=f.context.network.to_dict()
                record=prepare(f.context,f.pre,f.owner,f.donor,f.settings,**arguments(f))
                self.assertIs(type(record),NativeRootEnvelope)
                self.assertTrue(is_generated_junction_area(record.source))
                self.assertIs(record.source.layout.generated_junction_area,record)
                self.assertEqual(record.corner.radius,2.)
                self.assertLessEqual(record.report['curvature_bound'],.5*(1.+1e-8))
                self.assertLessEqual(record.report['bell_maximum_displacement'],1.75)
                self.assertEqual(record.report['reach_candidates'],1)
                self.assertEqual(len(record.staged_witnesses),1)
                record.reference.validate(f.context)
                token=((End(record.corner.first.source.id,record.corner.first.end),End(record.corner.second.source.id,record.corner.second.end)),
                    End(f.required.source.id,f.required.end))
                final=compound_cycles.source_cycle_group(f.planner,(f.site.key,),set(),set(),{},set(),
                    retry_dimensional=(f.dimensional_pair,),retry_support=((token,record),))
                generated=[source for source in final.area_footprints if is_generated_junction_area(source)]
                self.assertEqual(len(generated),1)
                final_root=next(c.root_support_envelope for c in final.physical_layout.corners if hasattr(c,'root_support_envelope'))
                self.assertIs(generated[0].layout.generated_junction_area,final_root)
                witness=record.staged_witnesses[0];a,s=_ends(f.raw)[witness.endpoint.section_id,witness.endpoint.end]
                b,t=_ends(final.physical_layout)[witness.endpoint.section_id,witness.endpoint.end]
                proof=admission(f.context,final,a,s,b,t)
                self.assertTrue(proof['root_support_did_not_advance_original_recipient'])
                self.assertEqual(t,witness.plane)
                self.assertTrue(final.physical_layout.source_cycle_full_width_proof)
                self.assertEqual(f.context.network.to_dict(),before)
                signature=construction_signature(f.context,final)
                self.assertIsNotNone(signature)
                final.physical_layout.generated_curb_prepared_neighbor_guard={'checked':999999}
                final.physical_layout.source_root_prepared_neighbor_guard={'checked':123456}
                self.assertEqual(signature,construction_signature(f.context,final))
                if not reverse and angle==0.:
                    changed=copy(final.physical_layout);stations=list(changed.outward_stations)
                    i=next(i for i,arm in enumerate(changed.approaches) if arm.source.id==3)
                    stations[i]+=.01;changed.outward_stations=tuple(stations)
                    with self.assertRaisesRegex(StagedRecipientUnavailable,'advanced'):
                        witness.validate_final(f.context,replace(final,physical_layout=changed),a,s,b,t)
                    root_corner=next(c for c in final.physical_layout.corners if hasattr(c,'root_support_envelope'))
                    actual=root_corner.root_support_envelope
                    try:
                        root_corner.root_support_envelope=replace(actual,contact_bounds=(actual.contact_bounds[0]+1.,actual.contact_bounds[1]))
                        with self.assertRaisesRegex(StagedRecipientUnavailable,'forged contact'):
                            witness.validate_final(f.context,final,a,s,b,t)
                        with self.assertRaisesRegex(OriginalRootTemplateUnavailable,'forged contact bounds'):
                            construction_signature(f.context,final)
                    finally:root_corner.root_support_envelope=actual

    def test_wrong_graph_unowned_support_and_original_bound_changes_are_refused(self):
        f=self.f
        with self.assertRaises(OriginalRootTemplateUnavailable):
            prepare(f.context,f.pre,f.owner,f.donor,f.settings,**dict(arguments(f),maximum_stations=(f.bounds[0]+.1,f.bounds[1])))
        with self.assertRaises(OriginalRootTemplateUnavailable):
            prepare(f.context,f.pre,Key('foreign',1),f.donor,f.settings,**arguments(f))
        other=replace(f.context.network,uuid='other-root-graph')
        other_context=BuildContext(other,f.context.products,f.context.patterns,{},f.context.source_digests)
        reference=OriginalRootReference.build(f.context,f.pre,f.site,f.owner,f.donor,f.raw.approaches,
            f.bounds,f.planner.ordinary_singleton_root_bounds)
        with self.assertRaises(OriginalRootTemplateUnavailable):reference.validate(other_context)
        changed_radius=replace(f.site,radius=3.,radii=(3.,3.,3.))
        changed_network=replace(f.context.network,intersections=(changed_radius,))
        with self.assertRaises(OriginalRootTemplateUnavailable):
            reference.validate(BuildContext(changed_network,f.context.products,f.context.patterns,{},f.context.source_digests))
        with self.assertRaises(OriginalRootTemplateUnavailable):
            prepare(f.context,f.pre,f.owner,f.donor,f.settings,**dict(arguments(f),site=replace(f.site,radius_modes=('STRICT',)*3)))
        foreign=SampledFootprint(Region(Key('foreign_area',1),np.array(((100.,100.),(101.,100.),(100.,101.)))),
            (RimSegment('curb',0,0.,1.),)*3,SimpleNamespace(approaches=(),corners=()))
        with self.assertRaisesRegex(OriginalRootTemplateUnavailable,'unowned physical'):
            prepare(f.context,replace(f.pre,footprints=(*f.pre.footprints,foreign)),f.owner,f.donor,f.settings,**arguments(f))
        previous=f.donor.radius_mode
        try:
            f.donor.radius_mode='STRICT'
            with self.assertRaisesRegex(OriginalRootTemplateUnavailable,'STYLE donor'):
                prepare(f.context,f.pre,f.owner,f.donor,f.settings,**arguments(f))
        finally:f.donor.radius_mode=previous

    def test_stale_profiles_original_planes_and_forged_witnesses_are_refused(self):
        f=self.f
        reference=OriginalRootReference.build(f.context,f.pre,f.site,f.owner,f.donor,f.raw.approaches,
            f.bounds,f.planner.ordinary_singleton_root_bounds)
        witness=reference.staged_witnesses[0]
        altered=dict(f.context.products);altered[Key('profile',1)]=Crown(.081,0.)
        with self.assertRaises(StagedRecipientUnavailable):
            witness.validate_pre(BuildContext(f.context.network,altered,f.context.patterns))
        with self.assertRaises(StagedRecipientUnavailable):replace(witness,previous=witness.previous+.01).validate_pre(f.context)
        with self.assertRaises(StagedRecipientUnavailable):replace(witness,plane=witness.plane+.01).validate_pre(f.context)
        wrong=replace(f.pre,members=(Key('intersection',999),))
        with self.assertRaises(StagedRecipientUnavailable):
            witness.validate_final(f.context,wrong,witness.approach,witness.previous,witness.approach,witness.plane)
        port=reference.ports[0]
        with self.assertRaises(OriginalRootPortUnavailable):
            OriginalRootPort.build(f.context,port.footprint,reference.ports[1].approach)


if __name__=='__main__':unittest.main()
