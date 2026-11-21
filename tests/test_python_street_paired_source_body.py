"""Two real three-arm source junctions bound one paired source pavement piece."""
from dataclasses import replace
from functools import lru_cache
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_compounds import footprint_context
from tested_python_streets import model
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.one_sided_paired_junctions import DerivedSourcePortKey
from tested_python_streets.paired_source_sampling import OneSidedPreparationRequest,canonical_reprepare
from tested_python_streets.ownership import CompoundProduct,surface_key
from tested_python_streets.prepared_merged_profiles import PreparedMergedMouth
from tested_python_streets.merged_run_inventory import SourceSpan,SpanEnd,inventory_frame,hydrate_evaluation_spans
from tested_python_streets.bundle_visible_runs import BundleRun,MemberRun
from tested_python_streets.bundle_angular_routes import AngularMemberRoute
from tested_python_streets.paired_source_body import PairedSourceBodyPiece
from tested_python_streets.uv import apply_mapping,PlanarMapping
from tested_python_streets.separator_ports import PORTS,SeparatorPortIndex


@lru_cache(maxsize=4)
def fixture(mode=0,reverse=False,rotation=0.,translation=(0.,0.),dip=0.):
    matrix=np.array(((np.cos(rotation),-np.sin(rotation)),(np.sin(rotation),np.cos(rotation))))
    shapes=((1,((-60.,4.),(-40.,4.),(-20.,4.)),10.),
            (2,((-20.,4.),(0.,4.),(20.,4.)),10.),
            (3,((20.,4.),(40.,4.),(60.,4.)),10.),
            (4,((-60.,-4.),(-20.,-4.),(0.,-4.-dip),(20.,-4.),(60.,-4.)),10.),
            (5,((-20.,4.),(-20.,30.)),2.),(6,((20.,4.),(20.,30.)),2.))
    sections=[]
    for sid,points,width in shapes:
        xy=np.asarray(points)@matrix.T+translation
        xyz=tuple((x,y,0.) for x,y in xy)
        sections.append(model.Section(sid,sid,xyz[::-1] if reverse else xyz,width=width,interpolation=mode))
    visits={sid:model.Traversal(sid,not reverse) for sid in (1,2,3,4)}
    sites=(model.Intersection(100,(model.End(1,int(not reverse)),model.End(2,int(reverse)),model.End(5,int(reverse))),radius=2.),
           model.Intersection(101,(model.End(2,int(not reverse)),model.End(3,int(reverse)),model.End(6,int(reverse))),radius=2.))
    bundle=model.Bundle(10,(model.BundleMember(tuple(visits[i] for i in (1,2,3))),model.BundleMember((visits[4],))))
    network=model.Network('two-original-one-sided-junctions',tuple(sections),sites,bundles=(bundle,),
        settings=model.Settings(exaggeration=1.))
    c=footprint_context(network)
    requests=tuple(OneSidedPreparationRequest(site,IntersectionProvider.footprint(c,site),
        tuple(visits[i] for i in near),visits[4],0) for site,near in zip(sites,((1,2),(2,3))))
    specifications,parameters,phase=canonical_reprepare(c,requests)
    products=[];mouths=[]
    for end,spec in enumerate(specifications):
        product=CompoundProduct.build(c,spec.owner,spec.group);products.append(product)
        c.products[surface_key(spec.owner)]=product
        site=next(site for site in sites if site.key==spec.owner)
        c.products[site.layout_key]=product.group.physical_layout
        indices=[i for i,a in enumerate(product.group.physical_layout.approaches)
                 if hasattr(a,'near_role') and a.source.id==2 and a.end_index==end]
        if len(indices)!=1:raise AssertionError('Fixture lost its actual inward-facing physical mouth')
        mouths.append(PreparedMergedMouth.from_product(product,indices[0]))
    tokens=(('bundle_end',0),('bundle_end',1))
    run=BundleRun(tuple(MemberRun(member.path,(model.End(member.path[0].section_id,int(not member.path[0].forward)),
        model.End(member.path[-1].section_id,int(member.path[-1].forward))),tokens) for member in bundle.members),tokens)
    frame=inventory_frame(c,bundle,SimpleNamespace(run=run,source_ranges=True))
    spans=[]
    for role,sid in ((0,2),(1,4)):
        path=c[model.Key('path',sid)];visit=visits[sid]
        station=np.asarray(parameters[sid]);q=-1. if visit.forward else 1.
        curves=tuple(path.point(station,path.section.offset+side*q*.5*path.section.width)
            for side in (1.,-1.))
        if not visit.forward:curves=tuple(points[::-1] for points in curves)
        aliases=[mouth.approach.source_mouths(mouth.station)[role] for mouth in mouths]
        ends=[]
        for mouth,alias in zip(mouths,aliases):
            key=None if alias.original_incidence else DerivedSourcePortKey(mouth.owner,sid,alias.endpoint.end)
            ends.append(SpanEnd('MERGED' if alias.original_incidence else 'DERIVED_MERGED',
                alias.endpoint if alias.original_incidence else None,mouth.owner,owner=mouth.owner,
                source_station=alias.station,merged_entrance=mouth.approach,registry_key=key))
        first,last=sorted(alias.station for alias in aliases)
        spans.append(SourceSpan(role,(visit,),AngularMemberRoute(c,model.BundleMember((visit,)),role),
            *ends,source_intervals=(dict(section=sid,forward=visit.forward,first=first,last=last),),curves=curves))
    c.task_digests={site.layout_key:'explicit-canonical-paired-source-fixture' for site in sites}
    c.task_digests.update({surface_key(site.key):'accepted-canonical-paired-source-fixture' for site in sites})
    c.products[PORTS]=SeparatorPortIndex.build(c)
    hydrated=hydrate_evaluation_spans(c,bundle,mouths[0].approach,specifications,
        spans=tuple((span,) for span in spans),frame=frame)
    return c,bundle,tuple(role[0] for role in hydrated),tuple(mouths),frame,tuple(products),phase


class PairedSourceBodyTests(unittest.TestCase):
    def test_two_actual_junction_products_connect_without_replacing_their_area_or_source_aliases(self):
        for mode,reverse,angle,offset in ((0,False,0.,(0.,0.)),(3,True,.713,(1200.,-1800.))):
            with self.subTest(mode=mode,reverse=reverse):
                c,bundle,spans,mouths,frame,products,phase=fixture(mode,reverse,angle,offset)
                before=c.network.to_dict();rims=[mouth.xyz.copy() for mouth in mouths]
                piece=PairedSourceBodyPiece()
                pavement,grass,report=piece.build(c,bundle,spans,mouths,frame=frame,namespace=('test piece',mode,reverse))
                for end,mouth in enumerate(mouths):
                    np.testing.assert_array_equal(mouth.xyz,rims[end])
                    np.testing.assert_array_equal(pavement.xyz[piece.front_ids[end]],mouth.xyz)
                    aliases=mouth.approach.source_mouths(mouth.station)
                    self.assertEqual([a.original_incidence for a in aliases],[True,False])
                    self.assertIn((mouth.owner,4,aliases[1].endpoint.end),products[end].derived_peaks)
                self.assertEqual(report['source_sections'],[2,4])
                self.assertTrue(phase['overlapping_station_arrays_exact'])
                self.assertFalse(report['junction_area_remeshed'])
                self.assertLess(report['total_quads'],16000)
                self.assertLessEqual(report['stored_checks'][0]['quality']['maximum_edge_ratio'],8.)
                self.assertGreaterEqual(report['stored_checks'][0]['quality']['minimum_corner_sine'],.03)
                self.assertTrue(np.isfinite(apply_mapping(pavement,PlanarMapping()).uv_layers['UVMap']).all())
                self.assertEqual(c.network.to_dict(),before)

    def test_wrong_original_interval_and_owner_scoped_alias_are_rejected(self):
        c,bundle,spans,mouths,frame,products,phase=fixture()
        wrong=replace(spans[1],head=replace(spans[1].head,registry_key=None))
        with self.assertRaisesRegex(ValueError,'derived alias registry'):
            PairedSourceBodyPiece().build(c,bundle,(spans[0],wrong),mouths,frame=frame,namespace='wrong')
        rows=(dict(spans[0].source_intervals[0],first=spans[0].source_intervals[0]['first']+.01),)
        with self.assertRaisesRegex(ValueError,'actual prepared mouth'):
            PairedSourceBodyPiece().build(c,bundle,(replace(spans[0],source_intervals=rows),spans[1]),
                mouths,frame=frame,namespace='wrong')
        with self.assertRaisesRegex(ValueError,'ordered actual head/tail'):
            PairedSourceBodyPiece().build(c,bundle,spans,mouths[::-1],frame=frame,namespace='wrong')
        changed=mouths[0].xyz.copy();changed[len(changed)//2,2]+=.001
        with self.assertRaisesRegex(ValueError,'completed physical junction rim'):
            PairedSourceBodyPiece().build(c,bundle,spans,(replace(mouths[0],xyz=changed),mouths[1]),frame=frame,namespace='wrong')
        false_original=replace(spans[1],head=replace(spans[1].head,endpoint=model.End(4,0)))
        with self.assertRaisesRegex(ValueError,'derived alias registry'):
            PairedSourceBodyPiece().build(c,bundle,(spans[0],false_original),mouths,frame=frame,namespace='wrong')
        with self.assertRaisesRegex(ValueError,'hydrated original visits'):
            PairedSourceBodyPiece().build(c,bundle,(replace(spans[0],curves=()),spans[1]),mouths,frame=frame,namespace='wrong')


if __name__=='__main__':unittest.main()
