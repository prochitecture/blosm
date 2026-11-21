"""Computational Bundle normals use exact owned ranges, never new claims."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import _python_street_package
from test_python_street_paired_interval_atlas import fixture
from tested_python_streets.model import Key
from tested_python_streets.routes import Route,Piece
from tested_python_streets.intersections import Approach
from tested_python_streets.updates import BuildContext
from tested_python_streets.source_bank_reconciliation import owned_ranges
from tested_python_streets.owned_bundle_normal_atlas import (OwnedBundleNormalScope,
    OwnedBundleNormalCuts,eligible,signature)


def setup(reverse=False):
    context,group=fixture(reverse)
    group.physical_layout.outward_stations=(2.,2.)
    # This family does not require, publish or infer terminal member claims.
    group.physical_layout.partial_bundle_claims=()
    group.physical_layout.partial_bundle_width_proof={}
    ranges=owned_ranges(context,group)
    proof=dict(sections=[dict(section=sid,interval=interval) for sid,rows in ranges.items() for interval in rows],maximum_error=.002)
    group.physical_layout.source_cycle_full_width_proof=proof
    return context,group,proof


class OwnedBundleNormalTests(unittest.TestCase):
    def test_actual_roles_select_only_proven_owned_native_ranges_without_partial_claims(self):
        for reverse in (False,True):
            context,group,proof=setup(reverse);before=context.network.to_dict();claims=group.claims
            self.assertTrue(eligible(group,context));scope=OwnedBundleNormalScope.build(context,group,proof)
            self.assertIs(scope.bundle,context.network.bundles[0])
            for role,expected,paired in ((0,{2,3,9},{2,3}),(1,{5,6,9},{5,6})):
                policy=OwnedBundleNormalCuts(scope,role,source_corridors=True,hole_bridges=True)
                self.assertEqual(policy._source_inventory(group,context),(expected,paired))
                self.assertEqual(policy.source_ranges,owned_ranges(context,group))
                self.assertEqual(policy.cadence_multiplier,2)
                self.assertEqual(policy.maximum_unpaired_span_widths,8.)
                self.assertEqual(policy.maximum_cuts,128)
            self.assertIs(group.claims,claims);self.assertEqual(context.network.to_dict(),before)
            self.assertFalse(group.physical_layout.partial_bundle_claims)
            self.assertTrue({Key(kind,sid) for sid in range(1,7) for kind in ('path','profile')}<=set(scope.dependencies))

    def test_disjoint_native_intervals_cannot_be_replaced_by_their_hull_or_shorter_scope(self):
        context,group,proof=setup();path=context[Key('path',2)]
        first=Approach(Route((Piece(path,True),)),0,context[Key('profile',2)])
        group.physical_layout.approaches=(*group.physical_layout.approaches,first)
        group.physical_layout.outward_stations=(2.,2.,2.)
        ranges=owned_ranges(context,group)
        proof=dict(sections=[dict(section=sid,interval=interval) for sid,rows in ranges.items() for interval in rows],maximum_error=.002)
        scope=OwnedBundleNormalScope.build(context,group,proof)
        self.assertEqual(dict(scope.ranges)[2],((0.,2.),(path.length-2.,path.length)))
        for changed in (((0.,path.length),),((0.,1.),(path.length-2.,path.length))):
            rows=[row for row in proof['sections'] if row['section']!=2]
            rows.extend(dict(section=2,interval=interval) for interval in changed)
            with self.assertRaisesRegex(ValueError,'actual prepared native ownership'):
                OwnedBundleNormalScope.build(context,group,dict(proof,sections=rows))

    def test_missing_role_ambiguous_bundle_and_incomplete_members_do_not_admit_a_scope(self):
        context,group,proof=setup();bundle=context.network.bundles[0]
        no_role=dict(proof,sections=[row for row in proof['sections'] if row['section'] not in (5,6)])
        with self.assertRaises(ValueError):OwnedBundleNormalScope.build(context,group,no_role)
        changed=replace(bundle,members=(replace(bundle.members[0],complete=False),bundle.members[1]))
        incomplete=BuildContext(replace(context.network,bundles=(changed,)),context.products,context.patterns)
        self.assertFalse(eligible(group,incomplete))
        with self.assertRaisesRegex(ValueError,'complete two-role Bundle'):
            OwnedBundleNormalScope.build(incomplete,group,proof)
        other=replace(bundle,id=bundle.id+100)
        ambiguous=BuildContext(replace(context.network,bundles=(bundle,other)),context.products,context.patterns)
        self.assertFalse(eligible(group,ambiguous))
        for role in (-1,2,True):
            with self.assertRaises(ValueError):OwnedBundleNormalCuts(OwnedBundleNormalScope.build(context,group,proof),role)

    def test_missing_width_proof_and_replaced_original_bundle_are_refused(self):
        context,group,proof=setup();scope=OwnedBundleNormalScope.build(context,group,proof)
        group.physical_layout.source_cycle_full_width_proof={}
        self.assertFalse(eligible(group,context));self.assertIsNone(signature(context,group))
        self.assertEqual(OwnedBundleNormalCuts(scope,0)._source_inventory(group,context),({2,3,9},{2,3}))
        changed=replace(scope,bundle=replace(scope.bundle))
        with self.assertRaisesRegex(ValueError,'original Bundle identity'):
            OwnedBundleNormalCuts(changed,0)._source_inventory(group,context)
        with self.assertRaisesRegex(ValueError,'source-width proof'):
            OwnedBundleNormalScope.build(context,group,{})

    def test_optional_cache_signature_never_rejects_an_unneeded_recovery(self):
        context,group,proof=setup()
        projected=SimpleNamespace(claims=(),physical_layout=None)
        self.assertFalse(eligible(projected,context));self.assertIsNone(signature(context,projected))
        first,*rest=proof['sections']
        stale=dict(proof,sections=[dict(first,interval=(first['interval'][0]+.01,first['interval'][1])),*rest])
        group.physical_layout.source_cycle_full_width_proof=stale
        self.assertTrue(eligible(group,context))
        self.assertIsNone(signature(context,group))
        with self.assertRaisesRegex(ValueError,'actual prepared native ownership'):
            OwnedBundleNormalScope.build(context,group,stale)
        group.physical_layout.source_cycle_full_width_proof=dict(proof,sections=[dict(section=row['section']) for row in proof['sections']])
        self.assertIsNone(signature(context,group))


if __name__=='__main__':unittest.main()
