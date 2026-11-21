"""Side compounds retain both original half-profile through ports."""
from dataclasses import replace
import unittest
import _python_street_package
from test_python_street_partial_bundle_claims import sources
from tested_python_streets.model import End,Key
from tested_python_streets.routes import Route,RouteTopology
from tested_python_streets.intersections import Approach
from tested_python_streets.retained_bundle_through import prepare


class RetainedBundleThroughTests(unittest.TestCase):
    @staticmethod
    def pair(reverse=False):
        context,bundle,_,members=sources(reverse)
        ends=(End(1,0 if reverse else 1),End(2,1 if reverse else 0))
        approaches=[Approach(Route.build(context,RouteTopology(context.network).trace(end)),end.end,
                    context[Key('profile',end.section_id)]) for end in ends]
        return context,bundle,approaches,members,ends

    def test_complete_source_pair_preserves_body_controls_and_both_oriented_half_profiles(self):
        for reverse in (False,True):
            context,bundle,approaches,members,ends=self.pair(reverse)
            before=context.network.to_dict()
            claim,=prepare(context,bundle,{5},approaches,members)
            self.assertEqual(claim.cuts,ends);self.assertEqual(claim.sections,())
            self.assertEqual(claim.junction,Key('intersection',10))
            self.assertFalse(claim.report()['bundle_body_consumed'])
            self.assertEqual(context.network.to_dict(),before)

    def test_missing_unowned_reprofiled_or_consumed_half_ports_do_not_grant_a_through_pair(self):
        context,bundle,approaches,members,_=self.pair()
        self.assertEqual(prepare(context,bundle,{5},approaches[:1],members),())
        self.assertEqual(prepare(context,bundle,{5},approaches,members-{Key('intersection',10)}),())
        self.assertEqual(prepare(context,bundle,{2},approaches,members),())
        original=approaches[0].offset
        approaches[0].offset+=.25
        self.assertEqual(prepare(context,bundle,{5},approaches,members),())
        approaches[0].offset=original
        approaches[0].profile=replace(approaches[0].profile,half=0)
        self.assertEqual(prepare(context,bundle,{5},approaches,members),())

    def test_nonadjacent_original_member_traversals_are_never_paired_by_position(self):
        context,bundle,approaches,members,_=self.pair()
        member=replace(bundle.members[0],path=tuple(reversed(bundle.members[0].path)))
        changed=replace(bundle,members=(member,bundle.members[1]))
        self.assertEqual(prepare(context,changed,{5},approaches,members),())
        incomplete=replace(bundle,members=(replace(bundle.members[0],complete=False),bundle.members[1]))
        self.assertEqual(prepare(context,incomplete,{5},approaches,members),())

    def test_consumption_of_the_opposite_member_does_not_hide_an_original_through_pair(self):
        for reverse in (False,True):
            context,bundle,approaches,members,ends=self.pair(reverse)
            before=context.network.to_dict()
            claim,=prepare(context,bundle,{4,5},approaches,members)
            self.assertEqual(claim.cuts,ends)
            self.assertEqual(claim.sections,())
            self.assertEqual(claim.junction,Key('intersection',10))
            # This permission grants no ownership of the other member and
            # cannot expose either a consumed or a missing half profile.
            self.assertEqual(prepare(context,bundle,{2,4,5},approaches,members),())
            self.assertEqual(prepare(context,bundle,{4,5},approaches[:1],members),())
            self.assertEqual(context.network.to_dict(),before)


if __name__=='__main__':unittest.main()
