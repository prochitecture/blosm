"""A truly free opposite source remains outside the complete owned exterior."""
from dataclasses import replace
import unittest
import copy
import _python_street_package
from test_python_street_source_cycles import triangle,prepare
from test_python_street_wholly_owned_bundle import fixture
from tested_python_streets import model
from tested_python_streets.wholly_owned_bundle_preparation import certify,WhollyOwnedBundleUnavailable
from tested_python_streets.ownership import bundle_signature
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.updates import BuildContext,UpdateCoordinator


def network(width,points,free_width):
    source=triangle(width)
    free=model.Section(7,7,tuple((x,y,0.) for x,y in points),width=free_width,interpolation=0)
    bundle=model.Bundle(9,(model.BundleMember((model.Traversal(1),)),model.BundleMember((model.Traversal(7),))))
    return replace(source,sections=(*source.sections,free),bundles=(bundle,))


class WhollyOwnedPreparationTests(unittest.TestCase):
    def test_ordinary_group_without_physical_layout_is_typed_unavailable(self):
        context,bundle,product=fixture()
        with self.assertRaisesRegex(WhollyOwnedBundleUnavailable,'physical layout'):
            certify(context,replace(product.group,physical_layout=None),store=False)

    def test_normal_source_exterior_records_independent_full_width_counterpart(self):
        context,bundle,product=fixture()
        proof=product.group.physical_layout.wholly_owned_free_counterparts
        self.assertTrue(proof['opposite_source_remains_independent'])
        self.assertEqual(proof['physical_source_connections_created'],0)
        self.assertEqual(proof['acquired_opposite_owners'],0)
        self.assertEqual(proof['candidates'][0]['free_section'],7)
        self.assertTrue(proof['candidates'][0]['outer_region_including_holes_excluded'])
        self.assertGreater(proof['candidates'][0]['certified_positive_clearance'],0.)
        with self.assertRaisesRegex(ValueError,'cell allocation'):
            certify(context,product.group,maximum_cells=0)

    def test_counterpart_inside_real_hole_is_ineligible_even_without_pavement_overlap(self):
        # Width-two source triangle has a genuine grass hole. Being in that
        # hole does not make an unconnected road an exterior counterpart.
        with self.assertRaisesRegex(ValueError,'complete free counterpart ribbon'):
            prepare(network(2.,((5.,3.),(7.,3.)),.25))

    def test_source_approximation_clearance_and_original_centerlines_remain_mandatory(self):
        with self.assertRaisesRegex(ValueError,'complete free counterpart ribbon'):
            prepare(network(8.,((3.,-5.501),(9.,-5.501)),3.))
        with self.assertRaisesRegex(ValueError,'centerlines cross or touch'):
            prepare(network(8.,((6.,-10.),(6.,10.)),1.))

    def test_valid_free_width_edit_rechecks_certificate_and_invalidates_owner_signature(self):
        context,bundle,product=fixture()
        group=replace(product.group,physical_layout=copy.copy(product.group.physical_layout))
        old=group.physical_layout.wholly_owned_free_counterparts['candidates'][0]
        changed=replace(context.network,sections=tuple(replace(section,width=section.width+.1)
            if section.id==7 else section for section in context.network.sections))
        update=UpdateCoordinator((PathProvider(),ProfileProvider()),context.patterns).prepare(changed)
        fresh=BuildContext(changed,{**context.products,**update.products},context.patterns,
            update.task_fingerprints,update.source_fingerprints)
        before_context=BuildContext(context.network,context.products,context.patterns,
            dict(update.task_fingerprints,**{}),update.source_fingerprints)
        before_context.task_digests={**update.task_fingerprints,model.Key('path',7):'original free Path',
                                    model.Key('profile',7):'original free Profile'}
        before=bundle_signature(before_context,group)
        certify(fresh,group)
        after=bundle_signature(fresh,group)
        self.assertNotEqual(before,after)
        new=group.physical_layout.wholly_owned_free_counterparts['candidates'][0]
        self.assertNotEqual(old['free_source_fingerprint'],new['free_source_fingerprint'])
        self.assertEqual(new['free_source_width'],3.1)
        self.assertTrue(group.physical_layout.wholly_owned_free_counterparts['opposite_source_remains_independent'])


if __name__=='__main__':unittest.main()
