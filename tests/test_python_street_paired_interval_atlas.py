"""One original member supplies cut axes; both source intervals stay owned."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import _python_street_package
from test_python_street_partial_bundle_claims import paired_ends
from tested_python_streets.partial_bundle_claims import prepare
from tested_python_streets.paired_interval_atlas import eligible,original_intervals,PairedIntervalCuts
from tested_python_streets.physical_hole_cuts import PhysicalHoleCuts
from tested_python_streets.model import End
from tested_python_streets.model import Key
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.routes import Route,RouteTopology
from tested_python_streets.intersections import Approach
from tested_python_streets.updates import BuildContext,UpdateCoordinator


def fixture(reverse=False):
    context,bundle,approaches,members=paired_ends()
    if reverse:
        def end(e):return End(e.section_id,1-e.end) if e.section_id<=6 else e
        network=replace(context.network,
            sections=tuple(replace(s,points=s.points[::-1]) if s.id<=6 else s for s in context.network.sections),
            intersections=tuple(replace(j,ends=tuple(map(end,j.ends))) for j in context.network.intersections),
            connections=tuple(replace(j,ends=tuple(map(end,j.ends))) for j in context.network.connections),
            bundles=(replace(bundle,members=tuple(replace(member,path=tuple(replace(v,forward=not v.forward)
                for v in member.path)) for member in bundle.members)),))
        products=UpdateCoordinator((PathProvider(),ProfileProvider()),PatternLibrary()).prepare(network).products
        context=BuildContext(network,products,PatternLibrary());bundle=network.bundles[0]
        approaches=[Approach(Route.build(context,RouteTopology(network).trace(e)),e.end,
            context[Key('profile',e.section_id)]) for e in (End(2,0),End(5,0))]
    ids={3,6,9}
    intervals=prepare(context,bundle,ids,approaches,members)
    layout=SimpleNamespace(approaches=approaches,partial_bundle_claims=intervals,
        partial_bundle_width_proof=dict(sections=[dict(section=sid) for sid in ids],maximum_error=.002),
        source_cycle_proofs={})
    group=SimpleNamespace(physical_layout=layout,members=tuple(sorted(members)),
        claims=tuple(SimpleNamespace(section_id=sid) for sid in sorted(ids)),
        loops=(SimpleNamespace(hole=False),SimpleNamespace(hole=True)))
    return context,group


class PairedIntervalAtlasTests(unittest.TestCase):
    def test_only_computational_normal_selection_changes_for_both_original_roles(self):
        for reverse in (False,True):
            with self.subTest(reverse=reverse):
                context,group=fixture(reverse);before=context.network.to_dict();claims=group.claims
                self.assertTrue(eligible(group))
                intervals=original_intervals(context,group)
                self.assertEqual(PhysicalHoleCuts()._source_inventory(group,context),({3,6,9},set()))
                for role,expected in ((0,{3,9}),(1,{6,9})):
                    selected,paired=PairedIntervalCuts(intervals,role)._source_inventory(group,context)
                    self.assertEqual(selected,expected);self.assertEqual(paired,set())
                    self.assertIs(group.claims,claims)
                    self.assertEqual({claim.section_id for claim in group.claims},{3,6,9})
                self.assertEqual(context.network.to_dict(),before)

    def test_original_bundle_intervals_and_complete_half_profile_cuts_are_rechecked(self):
        context,group=fixture();layout=group.physical_layout;original=layout.partial_bundle_claims
        for changed,reason in ((replace(original[0],bundle=999),'original Bundle'),
                               (replace(original[0],sections=(2,)),'original source interval'),
                               (replace(original[0],cuts=(End(1,1),)),'half-profile cut')):
            layout.partial_bundle_claims=(changed,original[1])
            with self.subTest(changed=changed),self.assertRaisesRegex(ValueError,reason):
                original_intervals(context,group)
        layout.partial_bundle_claims=original
        for proof in ({},dict(sections=[dict(section=3)],maximum_error=.00201)):
            layout.partial_bundle_width_proof=proof
            with self.assertRaisesRegex(ValueError,'source-width proof'):
                original_intervals(context,group)

    def test_retained_empty_sides_and_whole_bundles_keep_the_existing_policy(self):
        context,group=fixture();layout=group.physical_layout;original=layout.partial_bundle_claims
        layout.partial_bundle_claims=(original[0],replace(original[1],sections=()))
        self.assertFalse(eligible(group))
        layout.partial_bundle_claims=tuple(replace(claim,bundle_end='WHOLE') for claim in original)
        self.assertFalse(eligible(group))
        with self.assertRaisesRegex(ValueError,'original member role'):
            PairedIntervalCuts(original,2)


if __name__=='__main__':unittest.main()
