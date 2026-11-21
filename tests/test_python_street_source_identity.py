"""Equal scene rereads preserve provenance; copied or stale native paths do not."""
from copy import copy
from dataclasses import replace
from types import SimpleNamespace
import unittest
import _python_street_package
from tested_python_streets.model import Key, Network
from tested_python_streets.updates import BuildContext
from tested_python_streets.source_identity import matches_current_section
from tested_python_streets.junction_bank_transitions import original_corner
from tested_python_streets.owned_bundle_normal_atlas import OwnedBundleNormalScope
from test_python_street_retained_arc_fillets import primitive_fixture
from test_python_street_owned_bundle_normal_atlas import setup


class SourceIdentityTests(unittest.TestCase):
    def test_reconstructed_equal_input_keeps_the_canonical_native_corner(self):
        context,group=primitive_fixture()
        rebuilt=Network.from_dict(context.network.to_dict())
        reread=BuildContext(rebuilt,context.products,context.patterns)
        source=group.footprints[0];corner=source.layout.corners[0]
        node=next(node for node in rebuilt.intersections if node.key==source.layout.original_junction)
        self.assertIsNot(corner.first.source,rebuilt.section_index[corner.first.source.id])
        self.assertTrue(original_corner(SimpleNamespace(context=reread),node,corner))

    def test_both_owned_roles_keep_exact_native_ranges_after_a_scene_reread(self):
        context,group,proof=setup()
        initial=OwnedBundleNormalScope.build(context,group,proof)
        rebuilt=Network.from_dict(context.network.to_dict())
        reread=BuildContext(rebuilt,context.products,context.patterns)
        current=OwnedBundleNormalScope.build(reread,group,proof)
        self.assertEqual(initial.ranges,current.ranges)
        self.assertEqual(initial.dependencies,current.dependencies)
        self.assertIs(current.bundle,rebuilt.bundles[0])

    def test_current_geometry_dimensions_and_source_identity_must_still_match(self):
        context,_=primitive_fixture();section=context.network.sections[0]
        for changed in (replace(section,width=section.width+.01),replace(section,datum=.01),
                        replace(section,street_id=section.street_id+100),
                        replace(section,points=((section.points[0][0]+.01,*section.points[0][1:]),*section.points[1:]))):
            network=replace(context.network,sections=(changed,*context.network.sections[1:]))
            stale=BuildContext(network,context.products,context.patterns)
            with self.subTest(changed=changed):
                self.assertFalse(matches_current_section(stale,section))
        self.assertFalse(matches_current_section(context,replace(section)))

    def test_borrowed_native_path_stays_rejected_even_with_equal_source_values(self):
        context,group=primitive_fixture();source=group.footprints[0];corner=copy(source.layout.corners[0])
        borrowed=copy(corner.first);borrowed.path=copy(borrowed.path);corner.first=borrowed
        node=next(node for node in context.network.intersections if node.key==source.layout.original_junction)
        self.assertFalse(original_corner(SimpleNamespace(context=context),node,corner))
        products=dict(context.products);products.pop(Key('path',borrowed.source.id))
        missing=BuildContext(context.network,products,context.patterns)
        self.assertFalse(matches_current_section(missing,borrowed.source))


if __name__=='__main__':unittest.main()
