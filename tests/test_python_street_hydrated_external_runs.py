"""Actual merged and separated products bound finite multi-Section bodies."""
from dataclasses import replace
from functools import lru_cache
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_compounds import footprint_context
from tested_python_streets import model
from tested_python_streets.bundle_angular_curves import AngularCrown
from tested_python_streets.bundle_angular_routes import AngularMemberRoute
from tested_python_streets.bundle_merged_bridge import MergedProfileBridge
from tested_python_streets.bundle_visible_ports import ProfileSamplingPatterns
from tested_python_streets.bundle_visible_runs import BundleRun, MemberRun
from tested_python_streets.hydrated_external_runs import HydratedExternalRunWithCap
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.merged_run_inventory import SourceSpan, SpanEnd, inventory_frame, hydrate_evaluation_spans
from tested_python_streets.one_sided_paired_junctions import OneSidedPairedSpecification, DerivedSourcePortKey
from tested_python_streets.ownership import PLAN, CompoundProduct, OwnershipPlan, JunctionOwnership, footprint_key, surface_key
from tested_python_streets.prepared_merged_profiles import PreparedMergedMouth
from tested_python_streets.separator_ports import PORTS, SeparatorPortIndex
from tested_python_streets.providers import CarriagewayProvider
from tested_python_streets.routes import RouteTopology
from tested_python_streets.bundle_merged_body import MergedRunBody
from tested_python_streets.uv import apply_mapping, PlanarMapping


@lru_cache(maxsize=4)
def fixture(reverse=False, angle=0., offset=(0., 0.), short_terminal=False):
    rotate=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
    shapes=[]
    cut=89. if short_terminal else 60.
    for sid,first,last in ((10,-30.,0.),(1,0.,30.),(2,30.,cut),(3,cut,90.)):
        shapes.append((sid,((first,4.+.05*first),(last,4.+.05*last)),10.,sid))
    for sid,first,last in ((4,-30.,30.),(5,30.,60.),(6,60.,90.)):
        shapes.append((sid,((first,-4.-.05*first),(last,-4.-.05*last)),10.,sid))
    shapes.extend(((20,((0.,4.),(0.,30.)),2.,20),
        (30,((90.,-8.5),(90.,8.5)),3.,800),(31,((90.,8.5),(90.,30.)),3.,800),
        (32,((90.,-30.),(90.,-8.5)),3.,800)))
    sections=[]
    for sid,points,width,street in shapes:
        moved=np.asarray(points)@rotate.T+offset
        xyz=tuple((x,y,0.) for x,y in moved)
        sections.append(model.Section(sid,street,xyz[::-1] if reverse else xyz,width=width,interpolation=0))
    e=lambda sid,end:model.End(sid,1-end if reverse else end)
    sites=(model.Intersection(100,(e(10,1),e(1,0),e(20,0)),radius=2.),
           model.Intersection(201,(e(3,1),e(30,1),e(31,0)),radius=2.),
           model.Intersection(202,(e(6,1),e(32,1),e(30,0)),radius=2.))
    connections=tuple(model.Connection(900+i,(e(a,1),e(b,0))) for i,(a,b) in
                      enumerate(((1,2),(2,3),(4,5),(5,6))))
    v=lambda sid:model.Traversal(sid,not reverse)
    bundle=model.Bundle(50,(model.BundleMember(tuple(v(i) for i in (10,1,2,3))),
                            model.BundleMember(tuple(v(i) for i in (4,5,6)))))
    n=model.Network('hydrated-external-run',tuple(sections),sites,connections=connections,bundles=(bundle,),
                    settings=model.Settings(exaggeration=1.))
    c=footprint_context(n)
    for site in sites:
        c.products[footprint_key(site.key)]=IntersectionProvider.footprint(c,site)
    specification=OneSidedPairedSpecification.build(c,sites[0],c[footprint_key(sites[0].key)],
                                                   (v(10),v(1)),v(4),0)
    root=CompoundProduct.build(c,sites[0].key,specification.group)
    c.products[surface_key(sites[0].key)]=root
    c.products[sites[0].layout_key]=root.group.physical_layout
    c.products[PLAN]=OwnershipPlan((specification.group,),{sites[0].key:specification.group})
    for site in sites[1:]:
        c.products[site.layout_key]=IntersectionProvider.layout(c,site,c[footprint_key(site.key)])
        c.products[surface_key(site.key)]=JunctionOwnership.surface(c,site.key)
    c.task_digests={key:'accepted-explicit-source-test-product' for key in c.products}
    c.products[PORTS]=SeparatorPortIndex.build(c)
    index=next(i for i,a in enumerate(root.group.physical_layout.approaches)
               if hasattr(a,'near_role') and a.source.id==1)
    mouth=PreparedMergedMouth.from_product(root,index)
    boundaries=('bundle head','bundle tail')
    run=BundleRun(tuple(MemberRun(m.path,(e(m.path[0].section_id,0),e(m.path[-1].section_id,1)),boundaries)
                        for m in bundle.members),boundaries)
    frame=inventory_frame(c,bundle,SimpleNamespace(run=run,source_ranges=True))
    spans=[]
    for role,ids in enumerate(((1,2,3),(4,5,6))):
        alias=mouth.approach.source_mouths(mouth.station)[role]
        head=SpanEnd('MERGED' if alias.original_incidence else 'DERIVED_MERGED',
            alias.endpoint if alias.original_incidence else None,mouth.owner,owner=mouth.owner,
            source_station=alias.station,merged_entrance=mouth.approach,
            registry_key=None if alias.original_incidence else DerivedSourcePortKey(mouth.owner,ids[0],alias.endpoint.end))
        tail=SpanEnd('PHYSICAL',e(ids[-1],1),sites[role+1].key)
        visits=tuple(v(i) for i in ids)
        spans.append(SourceSpan(role,visits,AngularMemberRoute(c,model.BundleMember(visits),role),head,tail))
    spans=tuple(role[0] for role in hydrate_evaluation_spans(c,bundle,mouth.approach,(specification,),
        spans=tuple((span,) for span in spans),frame=frame))
    curves=tuple(curve for span in spans for curve in span.curves)
    field=AngularCrown(frame,curves,datum=mouth.approach.field.datum,rise=mouth.approach.field.rise)
    target=SimpleNamespace(frame=frame,field=field,source_sections=tuple(v.section_id for span in spans for v in span.visits),
                           source_visits=tuple(span.visits for span in spans),entrance=mouth.approach)
    bridge=MergedProfileBridge()
    connector,report=bridge.build(c,mouth.approach,mouth.station,mouth.xyz,target_inventory=target)
    return c,bundle,spans,mouth,frame,bridge,connector,report


class HydratedExternalRunTests(unittest.TestCase):
    def test_three_original_visits_per_role_keep_actual_profiles_and_external_cap(self):
        for reverse,angle,offset in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            with self.subTest(reverse=reverse):
                c,bundle,spans,mouth,frame,bridge,connector,_=fixture(reverse,angle,offset)
                source=c.network.to_dict();before=mouth.xyz.copy()
                body=HydratedExternalRunWithCap()
                pavement,grass,report=body.build(c,bundle,spans,mouth,bridge.body_cut,
                    frame=frame,namespace=('external-source-test',reverse))
                np.testing.assert_array_equal(pavement.xyz[body.merged_cut_ids],bridge.body_cut)
                np.testing.assert_array_equal(mouth.xyz,before)
                self.assertEqual(report['source_sections'],[1,2,3,4,5,6])
                self.assertTrue(all(len(span.visits)==3 for span in spans))
                self.assertTrue(report['legacy']['actual_external_rounded_cap'])
                self.assertFalse(report['runtime_discovery_enabled'])
                for (role,end),ids in body.external_profile_ids.items():
                    port=body.domain.ports[role,end]
                    np.testing.assert_array_equal(pavement.xyz[ids],port.mouth.physical_xyz)
                    self.assertIsInstance(body.domain.patterns_for(port),ProfileSamplingPatterns)
                for row in report['stored_checks']:
                    self.assertLessEqual(row['quality']['maximum_edge_ratio'],8.)
                    if row['material']!='grass':self.assertGreaterEqual(row['quality']['minimum_corner_sine'],.03)
                    self.assertLessEqual(row['maximum_slope'],.100001)
                    if row['material']!='grass':self.assertLessEqual(row['longitudinal_bound'],.03000001)
                    self.assertLessEqual(row['drainage']['maximum_fill'],1e-6)
                self.assertTrue(np.isfinite(apply_mapping(pavement,PlanarMapping()).uv_layers['UVMap']).all())
                self.assertTrue(np.isfinite(apply_mapping(grass,PlanarMapping()).uv_layers['UVMap']).all())
                self.assertEqual(c.network.to_dict(),source)
                # The final cap is an actual boundary of the two original
                # remote junctions and their connecting crossroad pavement.
                road=CarriagewayProvider.build(c,c.network.section_index[30],
                                               RouteTopology(c.network).attachments(c.network)[30]).mesh
                pv=(pavement,connector,road,*(c[surface_key(model.Key('intersection',sid))].pavement.mesh
                                            for sid in (100,201,202)))
                edges={tuple(sorted((tuple(mesh.xyz[a]),tuple(mesh.xyz[b]))))
                       for mesh in pv for a,b in MergedRunBody._boundary_vertices(mesh)}
                for a,b in MergedRunBody._boundary_vertices(grass):
                    self.assertIn(tuple(sorted((tuple(grass.xyz[a]),tuple(grass.xyz[b])))),edges)

    def test_unhydrated_profiles_changed_cut_and_missing_original_incidence_are_rejected(self):
        c,bundle,spans,mouth,frame,bridge,_,_=fixture()
        changed=replace(spans[0],tail=replace(spans[0].tail,port=None))
        with self.assertRaisesRegex(ValueError,'actual separated original sampled profiles'):
            HydratedExternalRunWithCap().build(c,bundle,(changed,spans[1]),mouth,bridge.body_cut,
                                               frame=frame,namespace='wrong')
        cut=bridge.body_cut.copy();cut[len(cut)//2,2]+=.001
        with self.assertRaisesRegex(ValueError,'target-frame merged crown'):
            HydratedExternalRunWithCap().build(c,bundle,spans,mouth,cut,frame=frame,namespace='wrong')
        rows=list(spans[0].source_intervals);rows[-1]=dict(rows[-1],last=rows[-1]['last']-.01)
        with self.assertRaisesRegex(ValueError,'complete physical profile'):
            HydratedExternalRunWithCap().build(c,bundle,(replace(spans[0],source_intervals=tuple(rows)),spans[1]),
                mouth,bridge.body_cut,frame=frame,namespace='wrong')
        original=c.network
        c.network=replace(original,connections=tuple(j for j in original.connections if j.id!=900))
        try:
            with self.assertRaisesRegex(ValueError,'consecutive source incidence'):
                HydratedExternalRunWithCap().build(c,bundle,spans,mouth,bridge.body_cut,frame=frame,namespace='wrong')
        finally:c.network=original


if __name__=='__main__':unittest.main()
