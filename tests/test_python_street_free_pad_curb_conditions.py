"""Whole-pad envelope conditions retain original profiles and source authority."""
from copy import copy
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_owned_free_pads import original_pad
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.ownership import surface_key
from tested_python_streets.updates import BuildContext
from tested_python_streets.free_pad_curb_conditions import validate,FreePadCurbUnavailable


def case():
    network=original_pad();backend=PythonMeshBackend(network.settings);update=backend.prepare(network)
    product=update.products[surface_key(network.intersections[0].key)]
    context=BuildContext(network,update.products,backend.patterns,update.task_fingerprints,update.source_fingerprints)
    group=product.group;sources={source.region.owner:source for source in group.footprints}
    source=next(iter(sources.values()));curve=source.layout.boundaries[0]
    return context,group,sources,curve,product


class FreePadCurbConditionTests(unittest.TestCase):
    def test_normal_product_zero_curbs_and_original_native_profiles_remain_exact(self):
        context,group,sources,curve,product=case();before=context.network.to_dict()
        domain=validate(context,group,sources,curve)
        self.assertEqual(domain.heights,(0.,0.))
        self.assertTrue(all(value>0. for value in domain.inward_slopes))
        for li,kinds in product.kinds.items():
            for i,kind in enumerate(kinds):
                if kind=='curb':self.assertAlmostEqual(product.physical_rims[li][i,2],0.,places=10)
        self.assertLessEqual(product.physical_pavement.slope_statistics()['maximum'],.1+1e-6)
        for (sid,end,u),(li,index,station) in product.port_points.items():
            profile=context[model.Key('profile',sid)].oriented(end)
            actual=product.physical_rims[li][index]
            self.assertAlmostEqual(actual[2],profile.datum+float(profile.value(2*u-1)),places=10)
        self.assertEqual(context.network.to_dict(),before)

    def test_borrowed_envelope_and_changed_complete_pad_claim_are_rejected(self):
        context,group,sources,curve,_=case()
        with self.assertRaisesRegex(FreePadCurbUnavailable,'borrowed'):
            validate(context,group,sources,copy(curve))
        layout=next(iter(sources.values())).layout;records=layout.finite_free_pads
        layout.finite_free_pads=(replace(records[0],length=records[0].length+.01),)
        try:
            with self.assertRaises(FreePadCurbUnavailable):validate(context,group,sources,curve)
        finally:layout.finite_free_pads=records

    def test_strict_style_and_new_bundle_membership_cannot_reuse_old_free_pad_proof(self):
        context,group,sources,curve,_=case();network=context.network
        node=network.intersections[0]
        strict=replace(network,intersections=(replace(node,radius_modes=('STRICT',)*3),))
        altered=BuildContext(strict,context.products,context.patterns)
        with self.assertRaises(FreePadCurbUnavailable):validate(altered,group,sources,curve)
        bundle=model.Bundle(9,tuple(model.BundleMember((model.Traversal(sid),)) for sid in (1,3)))
        altered=BuildContext(replace(network,bundles=(bundle,)),context.products,context.patterns)
        with self.assertRaises(FreePadCurbUnavailable):validate(altered,group,sources,curve)

    def test_equal_source_reread_preserves_the_same_native_resource_proof(self):
        context,group,sources,curve,_=case()
        reread=BuildContext(model.Network.from_dict(context.network.to_dict()),context.products,context.patterns)
        self.assertEqual(validate(context,group,sources,curve).heights,validate(reread,group,sources,curve).heights)


if __name__=='__main__':unittest.main()
