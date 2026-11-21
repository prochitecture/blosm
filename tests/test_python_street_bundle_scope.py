"""Ground Bundle extents retain original paired grade-transition incidence."""
from copy import deepcopy
from dataclasses import replace
from types import SimpleNamespace
import unittest

import _python_street_package
from tested_python_streets.adapter import PreparedPythonAdapter
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_scope import ground_bundle_scope
from tested_python_streets.model import Key, Settings, Section, BundleMember, Traversal, Network, Bundle, ground_extent_id


def fixture(deferred=(frozenset((2,)), frozenset((2,))), reverse=True):
    streets, members = {}, []
    for role in (0, 1):
        forward = not (role == 1 and reverse)
        items = []
        for index in range(4):
            a, b = ((20. * index, 5. - 10. * role), (20. * (index + 1), 5. - 10. * role))
            items.append(dict(kind='Section', id=role * 100 + index, category='residential',
                points=(a, b) if forward else (b, a), width=6., offset=0.,
                tags={'bridge': 'yes'} if index in deferred[role] else {}))
        streets[role] = dict(id=role, items=items if forward else items[::-1])
        members.append(dict(street_id=role, head_order=role, forward=forward, complete=True))
    return SimpleNamespace(streets=streets, bundles={7: dict(id=7, topology=dict(
        status='MAPPED', members=members, interior_streets=[]))}, junction_origins={'physical': []})


def convert(source):
    return PreparedPythonAdapter(('residential',),
        lambda tags: 'Bridge outside ground scope' if tags.get('bridge') else None,
        {'residential': {'radius': 2.}}).read(source, 'ground_bundle', Settings())


class GroundBundleScopeTests(unittest.TestCase):
    def test_split_member_keeps_the_unique_alongside_run_and_reports_ordinary_remainder(self):
        sections = {0: Section(0, 0, ((0., 5., 0.), (10., 5., 0.)), scope_ends=(1,)),
                    2: Section(2, 0, ((20., 5., 0.), (60., 5., 0.)), scope_ends=(0,)),
                    100: Section(100, 1, ((30., -5., 0.), (50., -5., 0.)))}
        original = (BundleMember(tuple(Traversal(sid) for sid in (0, 1, 2))),
                    BundleMember((Traversal(100),)))
        before = deepcopy(sections)
        members, record = ground_bundle_scope(7, original, sections, {1: 'bridge'})
        self.assertEqual(sections, before)
        self.assertEqual([[t.section_id for t in m.path] for m in members], [[2], [100]])
        self.assertEqual(record['unpaired_section_ids'], [0])
        self.assertEqual(record['deferred_section_ids'], [1])
        self.assertEqual(record['policy'], 'UNIQUE_FINITE_GROUND_RUN_PAIR')
        self.assertEqual(record['finite_source_selection']['overlap_length'], 20.)
        self.assertEqual(record['unpaired_ground_runs'][0]['boundaries'][0]['end'], 1)

        # Two plausible runs, a touching/crossing counterpart and a missing
        # scope incidence must remain errors rather than arbitrary selections.
        for points in (((0., -5., 0.), (60., -5., 0.)),
                       ((30., 5., 0.), (50., 5., 0.)),
                       ((30., -5., 0.), (50., 15., 0.))):
            invalid = dict(sections)
            invalid[100] = replace(sections[100], points=points)
            with self.assertRaisesRegex(ValueError, 'unique finite'):
                ground_bundle_scope(7, original, invalid, {1: 'bridge'})
        missing_cut = dict(sections)
        missing_cut[2] = replace(sections[2], scope_ends=())
        with self.assertRaisesRegex(ValueError, 'scope endpoint'):
            ground_bundle_scope(7, original, missing_cut, {1: 'bridge'})

    def test_corresponding_end_cuts_preserve_direction_controls_profiles_and_diagnostics(self):
        for deferred in (frozenset((0,)), frozenset((3,)), frozenset((0, 3))):
            for reverse in (False, True):
                with self.subTest(deferred=deferred, reverse=reverse):
                    source = fixture((deferred, deferred), reverse)
                    before = deepcopy(source)
                    adaptation = convert(source)
                    self.assertFalse(adaptation.blocked_sections)
                    self.assertEqual(source, before)
                    bundle, = adaptation.network.bundles
                    self.assertTrue(all(member.complete for member in bundle.members))
                    retained = tuple(i for i in range(4) if i not in deferred)
                    self.assertEqual([t.section_id for t in bundle.members[0].path], list(retained))
                    self.assertEqual([t.section_id for t in bundle.members[1].path], [100 + i for i in retained])
                    self.assertTrue(all(t.forward == (not reverse) for t in bundle.members[1].path))
                    record, = [r for r in adaptation.issues if r['kind'] == 'ScopeBundle']
                    self.assertEqual(record['original_section_ids'], [0, 1, 2, 3, 100, 101, 102, 103])
                    self.assertEqual(set(record['deferred_section_ids']),
                                     {role * 100 + i for role in (0, 1) for i in deferred})
                    self.assertTrue(record['complete_in_ground_scope'])
                    update = PythonMeshBackend(adaptation.network.settings).prepare(adaptation.network)
                    self.assertGreater(len(update.products[Key('separator', 7)].mesh.quads), 0)
                    for role, member in enumerate(bundle.members):
                        for traversal in member.path:
                            profile = update.products[Key('profile', traversal.section_id)]
                            expected = (1 if role == 0 else -1) * (1 if traversal.forward else -1)
                            self.assertEqual(profile.half, expected)
                            self.assertEqual(profile.rise, .08)

    def test_middle_holes_unequal_cuts_and_empty_member_do_not_invent_an_extent(self):
        for a, b in (((3,), (0,)), ((3,), ()), ((0, 1, 2, 3), (3,))):
            with self.subTest(first=a, second=b):
                result = convert(fixture((frozenset(a), frozenset(b))))
                self.assertFalse(result.network.bundles)
                self.assertTrue(result.blocked_sections)
                self.assertFalse(any(r['kind'] == 'ScopeBundle' for r in result.issues))

    def test_two_disjoint_ground_pairs_keep_all_roads_and_independent_half_profiles(self):
        source = fixture((frozenset((1,)), frozenset((1,))))
        before = deepcopy(source)
        result = convert(source)
        self.assertFalse(result.blocked_sections)
        self.assertEqual(source, before)
        bundles = result.network.bundles
        self.assertEqual([b.id for b in bundles], [ground_extent_id(7, 1), ground_extent_id(7, 2)])
        self.assertEqual([(b.source_bundle_id, b.scope_ordinal) for b in bundles], [(7, 1), (7, 2)])
        self.assertEqual([[[t.section_id for t in m.path] for m in b.members] for b in bundles],
                         [[[0], [100]], [[2, 3], [102, 103]]])
        self.assertEqual(set(result.network.section_index), {0, 2, 3, 100, 102, 103})
        self.assertEqual(Network.from_dict(result.network.to_dict()), result.network)
        update = PythonMeshBackend(result.network.settings).prepare(result.network)
        for bundle in bundles:
            part = update.products[Key('separator', bundle.id)]
            self.assertEqual(part.metadata['source_bundle_id'], 7)
            self.assertEqual(part.metadata['scope_ordinal'], bundle.scope_ordinal)
            for role, member in enumerate(bundle.members):
                for t in member.path:
                    self.assertEqual(update.products[Key('profile', t.section_id)].half,
                                     (1 if role == 0 else -1) * (1 if t.forward else -1))

    def test_mixed_grade_interior_retains_two_exterior_ground_approaches_and_scope_handoffs(self):
        source = fixture((frozenset((3,)), frozenset((3,))))
        points = (((10., 30.), (20., 25.)), ((20., 25.), (40., -25.)), ((40., -25.), (50., -30.)))
        source.streets[2] = dict(id=2, items=[dict(kind='Section', id=300+i, category='residential',
            width=6., offset=0., points=p, tags={'bridge':'yes'} if i == 1 else {}) for i,p in enumerate(points)])
        source.bundles[7]['topology']['interior_streets'] = [2]
        result = convert(source)
        self.assertFalse(result.blocked_sections)
        bundle, = result.network.bundles
        self.assertEqual(bundle.id, 7)
        self.assertEqual(bundle.interior_section_ids, ())
        self.assertTrue({300, 302} <= result.network.section_index.keys())
        self.assertEqual(result.network.section_index[300].scope_ends, (1,))
        self.assertEqual(result.network.section_index[302].scope_ends, (0,))
        record, = [r for r in result.issues if r['kind'] == 'ScopeBundle']
        roles = record['interior_source_roles']
        self.assertEqual([(r['section_ids'], r['role']) for r in roles],
            [([301], 'OUTSIDE_GROUND_SCOPE'), ([300], 'RETAINED_EXTERIOR_GROUND_RUN'),
             ([302], 'RETAINED_EXTERIOR_GROUND_RUN')])
        update = PythonMeshBackend(result.network.settings).prepare(result.network)
        self.assertEqual(update.products[Key('profile',300)].half, 0)
        self.assertGreater(len(update.products[Key('road',300)].mesh.quads), 0)
        self.assertGreater(len(update.products[Key('road',302)].mesh.quads), 0)
        # A disconnected ground piece in the finite median is not discarded.
        source.streets[2]['items'][0]['points'] = ((10., 0.), (20., 0.))
        invalid = convert(source)
        self.assertTrue(invalid.blocked_sections)
        self.assertFalse(invalid.network.bundles)
        self.assertTrue(any('intrudes' in r['reason'] for r in invalid.issues if r['kind'] == 'Bundle'))

    def test_derived_lineage_ranges_collisions_and_old_rows(self):
        result = convert(fixture((frozenset((3,)), frozenset((3,)))))
        data = result.network.to_dict()
        for row in data['bundles']:
            row.pop('source_bundle_id')
            row.pop('scope_ordinal')
        restored = Network.from_dict(data)
        self.assertEqual(restored, result.network)
        for source, ordinal in ((-1,1),(2**31,1),(7,0),(7,2**16),(7,True)):
            with self.assertRaisesRegex(ValueError,'range'):
                ground_extent_id(source,ordinal)
        with self.assertRaisesRegex(ValueError,'reserved runtime'):
            Bundle(7, restored.bundles[0].members, source_bundle_id=7, scope_ordinal=1)
        source = fixture((frozenset((1,)), frozenset((1,))))
        identity = ground_extent_id(7,1)
        source.bundles[identity] = dict(id=identity, topology={})
        result = convert(source)
        self.assertTrue(result.blocked_sections)
        self.assertFalse(result.network.bundles)
        self.assertTrue(any('collides' in r['reason'] for r in result.issues if r['kind'] == 'Bundle'))

    def test_invalid_dimensions_or_unsupported_category_are_not_grade_exclusions(self):
        for invalid in ('width', 'category'):
            source = fixture((frozenset((3,)), frozenset((3,))))
            item = source.streets[0]['items'][-1]
            item['tags'] = {}
            item[invalid] = -1. if invalid == 'width' else 'unknown_mode'
            result = convert(source)
            self.assertFalse(result.network.bundles)
            self.assertTrue(result.blocked_sections)

    def test_cut_requires_a_recorded_scope_end_and_complete_original_members(self):
        result = convert(fixture((frozenset((3,)), frozenset((3,)))))
        record, = [r for r in result.issues if r['kind'] == 'ScopeBundle']
        from tested_python_streets.model import BundleMember, Traversal
        original = tuple(BundleMember(tuple(Traversal(t['section'], t['forward'])
            for t in member['original'])) for member in record['members'])
        sections = result.network.section_index
        reasons = {sid: 'bridge' for sid in record['deferred_section_ids']}
        sections = {sid: replace(s, scope_ends=()) for sid, s in sections.items()}
        with self.assertRaisesRegex(ValueError, 'scope endpoint'):
            ground_bundle_scope(7, original, sections, reasons)
        with self.assertRaisesRegex(ValueError, 'complete original'):
            ground_bundle_scope(7, (replace(original[0], complete=False), original[1]),
                                result.network.section_index, reasons)


if __name__ == '__main__':
    unittest.main()
