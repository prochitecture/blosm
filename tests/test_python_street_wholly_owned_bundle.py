"""Whole source-owner banks meet true free source profiles without new Ends."""
from dataclasses import replace
from functools import lru_cache
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from test_python_street_source_cycles import triangle,prepare
from tested_python_streets import model
from tested_python_streets.bundle_axis import SeparatedPlan
from tested_python_streets.ownership import PLAN,OwnershipPlan,CompoundProduct,JunctionOwnership,surface_key,footprint_key
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.wholly_owned_bundle import WhollyOwnedBundleDomain,WhollyOwnedBundleSeparator,WhollyOwnedBundleUnavailable,_plane_anchor
from tested_python_streets.uv import apply_mapping,PlanarMapping
from tested_python_streets.bundles import BundleProvider
from tested_python_streets.bundle_compounds import UnequalBundleVisibilityError
from tested_python_streets.providers import CarriagewayProvider
from tested_python_streets.routes import RouteTopology


@lru_cache(maxsize=4)
def fixture(reverse=False,rotation=0.,translation=(0.,0.)):
    source=triangle(8.)
    sections=(*source.sections,model.Section(7,7,((3.,-10.,0.),(6.,-10.,0.),(9.,-10.,0.)),width=3.,interpolation=0))
    matrix=np.array(((np.cos(rotation),-np.sin(rotation)),(np.sin(rotation),np.cos(rotation))))
    transformed=[]
    for section in sections:
        xyz=np.asarray(section.points).copy();xyz[:,:2]=xyz[:,:2]@matrix.T+translation
        transformed.append(replace(section,points=tuple(map(tuple,xyz[::-1] if reverse else xyz))))
    sites=tuple(replace(site,ends=tuple(model.End(e.section_id,1-e.end) for e in site.ends)) for site in source.intersections) if reverse else source.intersections
    bundle=model.Bundle(9,(model.BundleMember((model.Traversal(1,not reverse),)),
                           model.BundleMember((model.Traversal(7,not reverse),))))
    source=replace(source,sections=tuple(transformed),intersections=sites,bundles=(bundle,))
    context,group=prepare(source)
    context.products[PLAN]=OwnershipPlan((group,),{key:group for key in group.members})
    product=CompoundProduct.build(context,group.members[0],group)
    context.products[surface_key(product.owner)]=product
    for site in source.intersections:
        context.products[footprint_key(site.key)]=IntersectionProvider.footprint(context,site)
        context.products[site.layout_key]=JunctionOwnership.layout(context,site)
    context.products[model.Key('bundle_plan',bundle.id)]=SeparatedPlan()
    return context,bundle,product


class WhollyOwnedBundleTests(unittest.TestCase):
    def test_free_plane_requires_one_crossing_in_both_bank_directions(self):
        points=np.array(((-1.,0.,.08),(.2,0.,.08),(.8,0.,.08),(1.2,0.,.08),(2.,0.,.08),
            (1.,1.,.08),(.8,1.,.08),(.2,1.,.08),(-1.,1.,.08)))
        for origin,end in (((0.,0.),0),((1.,0.),1)):
            with self.assertRaisesRegex(WhollyOwnedBundleUnavailable,'unique'):
                _plane_anchor(points,np.array(origin),np.array((1.,0.)),2.,end)
        for chain in (points[::-1],np.array(((0.,0.,.08),(0.,1.,.08),(1.,1.,.08)))):
            with self.assertRaises(WhollyOwnedBundleUnavailable):
                _plane_anchor(chain,np.array((0.,0.)),np.array((1.,0.)),2.,0)
        unique=points[:5]
        selected=_plane_anchor(unique,np.array((0.,0.)),np.array((1.,0.)),2.,0)
        self.assertEqual(selected[0],1)
        np.testing.assert_array_equal(unique[selected[0]],points[1])

    def test_absent_physical_layout_preserves_original_runtime_unavailability(self):
        context,bundle,product=fixture();previous=UnequalBundleVisibilityError('original unsupported visibility')
        changed=replace(product,group=replace(product.group,physical_layout=None))
        saved=dict(context.products)
        attachments=RouteTopology(context.network).attachments(context.network)
        for sid in bundle.sections:
            context.products[model.Key('road',sid)]=CarriagewayProvider.build(context,context.network.section_index[sid],attachments[sid])
        context.products[surface_key(product.owner)]=changed
        try:
            with self.assertRaisesRegex(WhollyOwnedBundleUnavailable,'physical layout'):
                WhollyOwnedBundleDomain.build(context,bundle)
            with patch('tested_python_streets.bundle_compounds.CompoundBundleBuilder.build',side_effect=previous):
                with self.assertRaises(UnequalBundleVisibilityError) as caught:BundleProvider.build(context,bundle)
                self.assertIs(caught.exception,previous)
        finally:
            context.products.clear();context.products.update(saved)

    def test_original_owned_bank_and_free_road_keep_every_whole_edge(self):
        for reverse,rotation,translation in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            with self.subTest(reverse=reverse):
                context,bundle,owner=fixture(reverse,rotation,translation)
                before=context.network.to_dict();original=owner.physical_rims[0].copy()
                domain=WhollyOwnedBundleDomain.build(context,bundle)
                part=WhollyOwnedBundleSeparator.build(context,bundle)
                self.assertEqual((domain.owned_section,domain.free_section),(1,7))
                self.assertEqual(part.metadata['source_aliases_created'],0)
                self.assertFalse(part.metadata['pavement_republished'])
                self.assertTrue(part.metadata['original_bank_whole_edges'])
                self.assertTrue(domain.proof['free_original_ends'])
                self.assertEqual(domain.proof['source_width']['maximum_error'],.002)
                for row in part.metadata['physical_and_stored_quality']:
                    self.assertLessEqual(row['maximum_edge_ratio'],8.)
                    self.assertGreater(row['minimum_corner'],0.)
                self.assertLessEqual(part.metadata['maximum_physical_slope'],.1+1e-6)
                self.assertEqual(part.metadata['drainage']['maximum_fill'],0.)
                self.assertTrue(np.isfinite(apply_mapping(part.mesh,PlanarMapping()).uv_layers['UVMap']).all())
                np.testing.assert_array_equal(owner.physical_rims[0],original)
                self.assertEqual(context.network.to_dict(),before)

    def test_scope_selection_and_missing_original_source_claim_stay_ineligible(self):
        context,bundle,product=fixture();network=context.network
        for flags in ({'scope_ends':(0,)},{'selection_ends':(1,)}):
            sections=tuple(replace(section,**flags) if section.id==7 else section for section in network.sections)
            context.network=replace(network,sections=sections)
            try:
                with self.assertRaisesRegex(WhollyOwnedBundleUnavailable,'unscoped free'):
                    WhollyOwnedBundleDomain.build(context,bundle)
            finally:context.network=network
        changed=replace(product,group=replace(product.group,claims=tuple(c for c in product.group.claims if c.section_id!=1)))
        context.products[surface_key(product.owner)]=changed
        try:
            with self.assertRaisesRegex(WhollyOwnedBundleUnavailable,'actual compound source claim'):
                WhollyOwnedBundleDomain.build(context,bundle)
        finally:context.products[surface_key(product.owner)]=product

    def test_runtime_dispatch_keeps_old_success_failure_and_eligible_fatal_error(self):
        context,bundle,product=fixture();attachments=RouteTopology(context.network).attachments(context.network)
        saved=dict(context.products)
        try:
            for sid in bundle.sections:
                context.products[model.Key('road',sid)]=CarriagewayProvider.build(context,context.network.section_index[sid],attachments[sid])
            with patch('tested_python_streets.bundle_compounds.CompoundBundleBuilder.build',return_value='old success'), \
                 patch('tested_python_streets.wholly_owned_bundle.WhollyOwnedBundleSeparator.build') as new:
                self.assertEqual(BundleProvider.build(context,bundle),'old success');new.assert_not_called()
            previous=UnequalBundleVisibilityError('original unsupported visibility')
            with patch('tested_python_streets.bundle_compounds.CompoundBundleBuilder.build',side_effect=previous), \
                 patch('tested_python_streets.wholly_owned_bundle.WhollyOwnedBundleSeparator.build',side_effect=WhollyOwnedBundleUnavailable('not eligible')):
                with self.assertRaises(UnequalBundleVisibilityError) as caught:BundleProvider.build(context,bundle)
                self.assertIs(caught.exception,previous)
            with patch('tested_python_streets.bundle_compounds.CompoundBundleBuilder.build',side_effect=previous), \
                 patch('tested_python_streets.wholly_owned_bundle.WhollyOwnedBundleSeparator.build',side_effect=ValueError('eligible source guard failed')):
                with self.assertRaisesRegex(ValueError,'eligible source guard failed'):BundleProvider.build(context,bundle)
        finally:
            context.products.clear();context.products.update(saved)


if __name__=='__main__':unittest.main()
