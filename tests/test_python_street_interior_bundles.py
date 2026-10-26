"""Interior crossing inventory, original ports and independent profiles."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import interior_fixture_networks, side_junction_fixture_networks
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_regions import CrossingFreeContinuation, continuation_sections, interior_crossings
from tested_python_streets.bundles import BundleProvider, Cap, Member
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import Crown
from tested_python_streets.updates import BuildContext


def fixture():
    network = side_junction_fixture_networks(model)[1]
    return replace(network, bundles=(replace(network.bundles[0], interior_section_ids=(4,)),))


def continuation_fixture():
    return interior_fixture_networks(model)[0]


class PythonInteriorBundleTests(unittest.TestCase):
    def test_convex_but_badly_graded_cells_are_never_silently_accepted(self):
        # A long rectangular median with independently dense end-cap samples.
        # The old tensor product creates 48,000 convex but elongated quads.
        network = model.Network('graded_separator', (
            model.Section(1, 1, ((0., 12.7, 0.), (280., 12.7, 0.)), width=5.4),
            model.Section(2, 2, ((0., -2.7, 0.), (280., -2.7, 0.)), width=5.4)),
            bundles=(model.Bundle(7, (model.BundleMember((model.Traversal(1),)),
                                     model.BundleMember((model.Traversal(2),)))),))
        bundle, = network.bundles
        x = np.linspace(0., 280., 401)
        left = np.column_stack((x, np.full_like(x, 10.), np.full_like(x, .8)))
        right = np.column_stack((x, np.zeros_like(x), np.full_like(x, .8)))
        members = (Member((), left, (model.End(1, 0), model.End(1, 1))),
                   Member((), right, (model.End(2, 0), model.End(2, 1))))
        caps = [Cap.free(left[i], right[i], 120) for i in (0, -1)]
        context = BuildContext(network, {model.Key('profile', 1): Crown(.08, 0., 1),
            model.Key('profile', 2): Crown(.08, 0., -1)}, PatternLibrary())
        # This deliberately severe sample imbalance is not yet handled by the
        # available charts. Positive cells alone must not bypass their limits.
        with self.assertRaisesRegex(ValueError, 'No accepted fixed-rim'):
            BundleProvider.region(context, bundle, members, caps, 0)
    def test_explicit_interior_keeps_the_same_crossing_pavement_and_two_separators(self):
        before = side_junction_fixture_networks(model)[1]
        before = replace(before, bundles=(replace(before.bundles[0], interior_section_ids=()),))
        after = fixture()
        crossing, = interior_crossings(after, after.bundles[0])
        self.assertEqual(crossing.sections, (4,))
        self.assertEqual(crossing.sites, (model.Key('intersection', 1), model.Key('intersection', 2)))
        expected = PythonMeshBackend(before.settings).prepare(before)
        actual = PythonMeshBackend(after.settings).prepare(after)
        self.assertEqual(len(actual.products[model.Key('separator', 7)].metadata['runs']), 2)
        self.assertEqual(actual.products[model.Key('profile', 4)].half, 0)
        for key, part in expected.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, actual.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, actual.products[key].mesh.quads)

    def test_a_complete_multisection_crossing_preserves_authored_directions(self):
        network = fixture()
        sections = tuple(replace(s, points=(s.points[0], s.points[1])) if s.id == 4 else s
                         for s in network.sections)
        sections += (model.Section(6, 4, ((0., -8., 0.), (0., 0., 0.)), width=5.4),)
        junctions = tuple(replace(j, ends=tuple(model.End(6, 0) if e == model.End(4, 1) else e
                                              for e in j.ends)) for j in network.intersections)
        bundle = replace(network.bundles[0], interior_section_ids=(6, 4))
        network = replace(network, sections=sections, intersections=junctions, bundles=(bundle,),
            connections=(model.Connection(3, (model.End(4, 1), model.End(6, 1))),))
        crossing, = interior_crossings(network, bundle)
        self.assertEqual(crossing.path, (model.Traversal(4), model.Traversal(6, False)))
        incomplete = replace(bundle, interior_section_ids=(4,))
        with self.assertRaisesRegex(ValueError, 'every Section'):
            interior_crossings(network, incomplete)

    def test_disconnected_or_branching_interiors_are_not_silently_omitted(self):
        network = fixture()
        for case in ('free_end', 'wrong_member', 'isolated'):
            with self.subTest(case=case):
                changed = network
                if case == 'free_end':
                    changed = replace(network, intersections=network.intersections[:1])
                if case == 'wrong_member':
                    changed = replace(network, intersections=tuple(replace(j,
                        ends=(model.End(4, 1), model.End(2, 1), model.End(5, 1)))
                        if j.id == 2 else j for j in network.intersections))
                if case == 'isolated':
                    changed = replace(network, sections=(*network.sections,
                        model.Section(8, 8, ((-20., 0., 0.), (-10., 0., 0.)))),
                        bundles=(replace(network.bundles[0], interior_section_ids=(4, 8)),))
                with self.assertRaises(ValueError):
                    interior_crossings(changed, changed.bundles[0])

    def test_crossing_street_can_continue_outside_both_physical_member_ports(self):
        source = fixture()
        source = replace(source, sections=(*source.sections,
            model.Section(6, 4, ((0., 8., 0.), (0., 36., 0.)), width=5.4),
            model.Section(7, 4, ((0., -8., 0.), (0., -36., 0.)), width=5.4)),
            intersections=tuple(replace(site, ends=(*site.ends, model.End(6 if site.id == 1 else 7, 0)))
                                for site in source.intersections))
        original = source.to_dict()
        crossing, = interior_crossings(source, source.bundles[0])
        self.assertEqual(crossing.sections, (4,))
        self.assertEqual(crossing.street_id, 4)
        self.assertEqual(crossing.sites, (model.Key('intersection', 1), model.Key('intersection', 2)))
        # Only the persistent Street identity differs in this reference. All
        # original vertices, physical incidences, widths and profiles match.
        reference = replace(source, sections=tuple(replace(section, street_id=section.id)
            if section.id in (6, 7) else section for section in source.sections))
        expected = PythonMeshBackend(reference.settings).prepare(reference)
        actual = PythonMeshBackend(source.settings).prepare(source)
        for key, part in expected.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, actual.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, actual.products[key].mesh.quads)
        for sid in (4, 6, 7):
            self.assertEqual(actual.products[model.Key('profile', sid)].half, 0)
        self.assertEqual(source.to_dict(), original)

    def test_interior_inventory_cannot_include_a_same_street_exterior_branch(self):
        source = fixture()
        source = replace(source, sections=(*source.sections,
            model.Section(6, 4, ((0., 8., 0.), (0., 36., 0.)), width=5.4)),
            intersections=tuple(replace(site, ends=(*site.ends, model.End(6, 0)))
                                if site.id == 1 else site for site in source.intersections),
            bundles=(replace(source.bundles[0], interior_section_ids=(4, 6)),))
        with self.assertRaises(ValueError):
            interior_crossings(source, source.bundles[0])

    def test_crossing_to_free_region_retains_original_membership_and_full_profile(self):
        network = continuation_fixture()
        bundle = network.bundles[0]
        crossing, = interior_crossings(network, bundle)
        plan = CrossingFreeContinuation.build(network, bundle, crossing, 1, 1)
        self.assertEqual(continuation_sections(network, bundle), {5})
        self.assertEqual(plan.member.path, (model.Traversal(5),))
        result = PythonMeshBackend(network.settings).prepare(network)
        separator = result.products[model.Key('separator', 7)]
        self.assertEqual(len(separator.metadata['runs']), 2)
        self.assertEqual(separator.metadata['derived_continuations'], [plan.record()])
        self.assertNotIn(5, bundle.sections)
        self.assertEqual(result.products[model.Key('profile', 5)].half, 0)
        self.assertEqual(result.products[model.Key('profile', 5)].rise, .12)
        self.assertGreater(separator.mesh.validate()['minimum_corner'], 0.)
        self.assertLessEqual(separator.mesh.validate()['maximum_edge_ratio'], 8.)
        self.assertTrue(np.all(separator.mesh.xyz[:, 2] >= -1e-12))
        self.assertTrue(np.all(separator.mesh.xyz[:, 2] <= .8 + 1e-12))

    def test_continuation_at_the_other_end_and_authored_reverse_direction(self):
        network = continuation_fixture()
        bundle = network.bundles[0]
        reversed_members = tuple(model.BundleMember(tuple(model.Traversal(t.section_id, not t.forward)
            for t in reversed(member.path))) for member in reversed(bundle.members))
        network = replace(network, bundles=(replace(bundle, members=reversed_members),),
            sections=tuple(replace(s, points=s.points[::-1]) if s.id == 5 else s for s in network.sections),
            intersections=tuple(replace(j, ends=tuple(model.End(5, 1) if e == model.End(5, 0) else e
                for e in j.ends)) for j in network.intersections))
        result = PythonMeshBackend(network.settings).prepare(network)
        part = result.products[model.Key('separator', 7)]
        record, = part.metadata['derived_continuations']
        self.assertEqual(record['bundle_end'], 0)
        self.assertEqual(record['member_side'], 0)
        self.assertEqual(record['path'], [dict(section_id=5, forward=True)])
        self.assertEqual(len(part.metadata['runs']), 2)

    def test_full_profile_continuation_edits_update_both_regions_without_stale_outputs(self):
        network = continuation_fixture()
        backend = PythonMeshBackend(network.settings)
        before = backend.prepare(network)
        backend.commit(before)
        changed = replace(network, sections=tuple(replace(s, width=6.2) if s.id == 5 else s
                                                  for s in network.sections))
        actual = backend.prepare(changed)
        self.assertIn(model.Key('separator', 7), actual.rebuilt)
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        for key, part in actual.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, fresh.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, fresh.products[key].mesh.quads)
        backend.commit(actual)
        restored = backend.prepare(network)
        for key, part in before.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, restored.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, restored.products[key].mesh.quads)

    def test_ambiguous_or_attached_continuations_do_not_invent_a_free_boundary(self):
        network = continuation_fixture()
        crossing, = interior_crossings(network, network.bundles[0])
        for case in ('ambiguous', 'attached'):
            extra = model.Section(6, 6, ((0., -8., 0.), (0., -30., 0.))) if case == 'ambiguous' else (
                model.Section(6, 6, ((40., -8., 0.), (60., -8., 0.))))
            junctions = tuple(replace(j, ends=(*j.ends, model.End(6, 0))) if j.id == 2 else j
                              for j in network.intersections) if case == 'ambiguous' else network.intersections
            changed = replace(network, sections=(*network.sections, extra), intersections=junctions,
                connections=() if case == 'ambiguous' else (
                    model.Connection(3, (model.End(5, 1), model.End(6, 0))),))
            with self.subTest(case=case), self.assertRaises(ValueError):
                CrossingFreeContinuation.build(changed, changed.bundles[0], crossing, 1, 1)

        competing = model.Bundle(88, (model.BundleMember((model.Traversal(5),)),
                                     model.BundleMember((model.Traversal(99),))))
        changed = replace(network, bundles=(*network.bundles, competing))
        with self.assertRaisesRegex(ValueError, 'competing'):
            CrossingFreeContinuation.build(changed, changed.bundles[0], crossing, 1, 1)

    def test_a_later_section_of_the_derived_street_remains_in_the_update_closure(self):
        network = continuation_fixture()
        sections = tuple(replace(s, points=((0., -8., 0.), (20., -8., 0.))) if s.id == 5 else s
                         for s in network.sections)
        sections += (model.Section(6, 5, ((40., -8., 0.), (20., -8., 0.)), width=5.4, rise=.12),)
        network = replace(network, sections=sections, furniture=(),
            connections=(model.Connection(3, (model.End(5, 1), model.End(6, 1))),))
        self.assertEqual(continuation_sections(network, network.bundles[0]), {5, 6})
        backend = PythonMeshBackend(network.settings)
        original = backend.prepare(network)
        backend.commit(original)
        changed = replace(network, sections=tuple(replace(s, width=6.2) if s.id == 6 else s for s in network.sections))
        actual = backend.prepare(changed)
        self.assertIn(model.Key('separator', 7), actual.rebuilt)
        self.assertNotIn(model.Key('road', 99), actual.rebuilt)
        record, = actual.products[model.Key('separator', 7)].metadata['derived_continuations']
        self.assertEqual(record['path'], [dict(section_id=5, forward=True), dict(section_id=6, forward=False)])
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        for key, part in actual.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, fresh.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, fresh.products[key].mesh.quads)

    def test_translated_common_datum_and_zero_exaggeration_keep_the_same_regions(self):
        network = continuation_fixture()
        for exaggeration in (0., 10.):
            changed = replace(network, sections=tuple(replace(s, datum=7.3) for s in network.sections),
                settings=replace(network.settings, exaggeration=exaggeration))
            result = PythonMeshBackend(changed.settings).prepare(changed)
            part = result.products[model.Key('separator', 7)]
            self.assertEqual(len(part.metadata['runs']), 2)
            self.assertGreaterEqual(float(part.mesh.xyz[:, 2].min()), 7.3 - 1e-12)
            self.assertLessEqual(float(part.mesh.xyz[:, 2].max()), 7.3 + exaggeration * .08 + 1e-12)


if __name__ == '__main__':
    unittest.main()
