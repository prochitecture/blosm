"""Separator dispatch must follow prepared endpoint identities, never proximity."""
from copy import deepcopy
import importlib
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace as S
import unittest

package = ModuleType('_original_separator_tests')
package.__path__ = [str(Path(__file__).resolve().parents[1]/'way/street_network')]
sys.modules[package.__name__] = package
m = importlib.import_module(package.__name__+'.original_span_separators')


def fixture():
    def route(sid, start, end):
        return S(kind='CORRIDOR', span=S(street_id=sid, section_ids=(sid*10,),
            in_scope=True,
            first_control=0, last_control=2, source_control_count=3,
            start=S(junction_id=start), end=S(junction_id=end)))
    members = [dict(street_id=1, forward=True, head_order=0, tail_order=0, complete=True, spans=[{}]),
               dict(street_id=2, forward=False, head_order=1, tail_order=1, complete=True, spans=[{}])]
    snapshot = S(streets={sid: dict(items=[dict(kind='Section', id=sid*10)]) for sid in (1, 2)},
                 bundles={7: dict(id=7, topology=dict(status='MAPPED', members=members,
                     interior_streets=[], junctions=[], required_features=[]))})
    component = S(ports=S(source_streets=(1, 2)), network=S(ordinary_site_ids=(11, 12),
        compounds=(), logical_clusters=(), routes=(route(1, 11, 12), route(2, 12, 11))))
    return snapshot, component


class OriginalSpanSeparatorTests(unittest.TestCase):
    def setUp(self):
        self.snapshot, self.component = fixture()
        self.planner = m.OriginalSpanSeparatorPlanner()

    def plan(self):
        return self.planner.prepare(self.snapshot, self.component)

    def test_reversed_member_preserves_declared_head_and_tail_without_mutation(self):
        before = deepcopy((self.snapshot, self.component))
        result = self.plan()
        self.assertEqual(result.closed, (m.ClosedBundleSeparator(7, (1, 2), (11,), (12,)),))
        self.assertFalse(result.pending)
        self.assertEqual((self.snapshot, self.component), before)

    def test_distinct_physical_ends_can_share_one_logical_surface_owner(self):
        n = self.component.network
        n.ordinary_site_ids = ()
        n.logical_clusters = (S(physical_junction_ids=(11, 13)), S(physical_junction_ids=(12, 14)))
        n.routes[1].span.start.junction_id = 14
        n.routes[1].span.end.junction_id = 13
        self.assertEqual(self.plan().closed[0].head_sites, (11, 13))

    def test_invalid_ownership_is_not_silently_a_pending_geometry_case(self):
        self.component.network.compounds = (S(site_ids=(11, 12)),)
        with self.assertRaisesRegex(ValueError, 'multiple surface owners'):
            self.plan()

    def test_partial_or_unresolved_member_never_becomes_a_closed_separator(self):
        for change in ('partial', 'unresolved', 'missing', 'free'):
            with self.subTest(change=change):
                self.setUp()
                span = self.component.network.routes[0].span
                if change == 'partial': span.first_control = 1
                if change == 'unresolved': self.component.network.routes[0].kind = 'REMOTE_CONNECTION'
                if change == 'missing': self.component.ports.source_streets = (2,)
                if change == 'free': span.start.junction_id = None
                result = self.plan()
                self.assertFalse(result.closed)
                self.assertEqual(len(result.pending), 1)

    def test_split_extent_and_order_change_remain_explicit_pending_cases(self):
        for change in ('split', 'order', 'interior', 'attachment', 'partial_extent'):
            with self.subTest(change=change):
                self.setUp(); t = self.snapshot.bundles[7]['topology']
                if change == 'split': t['members'][0]['spans'].append({})
                if change == 'order': t['members'][0]['tail_order'] = 1
                if change == 'interior': t['interior_streets'] = [3]
                if change == 'attachment': t['junctions'] = [dict(junction_id=99)]
                if change == 'partial_extent': t['members'][0]['complete'] = False
                result = self.plan()
                self.assertFalse(result.closed)
                self.assertEqual(result.pending[0]['status'], 'PENDING')

    def test_unrelated_component_does_not_claim_or_report_a_bundle(self):
        self.component.ports.source_streets = (99,)
        self.assertEqual(self.plan(), m.SeparatorPreparation((), ()))

    def test_coincident_owner_is_not_a_closed_two_terminal_space(self):
        self.component.network.ordinary_site_ids = ()
        self.component.network.compounds = (S(site_ids=(11, 12)),)
        self.assertFalse(self.plan().closed)
        self.assertIn('two distinct', self.plan().pending[0]['reason'])

    def interior(self, categories=('footway', 'cycleway')):
        self.snapshot.bundles[7]['topology'].update(interior_streets=[3], required_features=['INTERIOR_STREETS'])
        self.snapshot.streets[3] = dict(items=[dict(kind='Section', id=30+i, category=category)
                                             for i, category in enumerate(categories)])

    def test_explicitly_deferred_paths_retain_identities_and_do_not_change_closed_boundary(self):
        self.interior(); before = deepcopy(self.snapshot)
        result = m.OriginalSpanSeparatorPlanner(('residential',)).prepare(self.snapshot, self.component)
        self.assertEqual(result.closed, (m.ClosedBundleSeparator(7, (1, 2), (11,), (12,)),))
        self.assertEqual(result.deferred_interiors, (dict(bundle_id=7, streets=(
            dict(street_id=3, section_ids=(30, 31), categories=('cycleway', 'footway')),)),))
        self.assertFalse(result.pending)
        self.assertEqual(self.snapshot, before)

    def test_deferred_paths_need_scope_and_complete_original_evidence(self):
        for change in ('no_scope', 'missing', 'empty', 'unknown', 'duplicate', 'active_section', 'active_route'):
            with self.subTest(change=change):
                self.setUp(); self.interior()
                planner = m.OriginalSpanSeparatorPlanner(None if change == 'no_scope' else ('residential',))
                if change == 'missing': del self.snapshot.streets[3]
                if change == 'empty': self.snapshot.streets[3]['items'] = []
                if change == 'unknown': self.snapshot.streets[3]['items'][0]['category'] = None
                if change == 'duplicate': self.snapshot.bundles[7]['topology']['interior_streets'].append(3)
                if change == 'active_section': self.snapshot.streets[3]['items'][1]['category'] = 'residential'
                if change == 'active_route':
                    self.component.network.routes += (S(span=S(street_id=3, in_scope=True)),)
                result = planner.prepare(self.snapshot, self.component)
                self.assertFalse(result.closed)
                self.assertFalse(result.deferred_interiors)
                self.assertEqual(len(result.pending), 1)

    def test_enabling_an_interior_mode_requires_its_surface_to_be_resolved(self):
        self.interior(('cycleway',))
        result = m.OriginalSpanSeparatorPlanner(('residential', 'cycleway')).prepare(self.snapshot, self.component)
        self.assertFalse(result.closed)
        self.assertIn('Active Bundle interior', result.pending[0]['reason'])

    def test_deferred_paths_do_not_hide_an_attachment_or_unknown_feature(self):
        for feature in ('INTERNAL_JUNCTIONS', 'UNRECOGNIZED_FEATURE'):
            with self.subTest(feature=feature):
                self.setUp(); self.interior()
                self.snapshot.bundles[7]['topology']['required_features'].append(feature)
                result = m.OriginalSpanSeparatorPlanner(('residential',)).prepare(self.snapshot, self.component)
                self.assertFalse(result.closed)
                self.assertEqual(len(result.pending), 1)


class FreeEndSeparatorTests(unittest.TestCase):
    def setUp(self):
        self.snapshot, self.component = fixture()
        self.component.network.ordinary_site_ids = (11,)
        for route, attached_index in zip(self.component.network.routes, (0, 1)):
            route.kind = 'APPROACH'
            ends, ports = [], []
            for index in (0, 1):
                endpoint = m.SectionEnd(route.span.street_id, route.span.section_ids[0], index)
                attached = index == attached_index
                ends.append(S(endpoint=endpoint, owner_kind='PHYSICAL' if attached else 'FREE',
                              junction_id=11 if attached else None))
                ports.append(S(endpoint=endpoint, junction_id=11) if attached else None)
            route.span.start, route.span.end = ends
            route.start_port, route.end_port = ports
        self.planner = m.OriginalSpanSeparatorPlanner()

    def plan(self):
        return self.planner.prepare(self.snapshot, self.component, free_ends=True)

    def test_preparation_is_opt_in_until_native_open_boundary_support_is_installed(self):
        result = self.planner.prepare(self.snapshot, self.component)
        self.assertFalse(result.free_ended)
        self.assertFalse(result.closed)
        self.assertEqual(len(result.pending), 1)

    def test_reversed_street_keeps_its_actual_free_endpoint(self):
        before = deepcopy((self.snapshot, self.component))
        result = self.plan()
        self.assertEqual(result.free_ended, (m.FreeEndBundleSeparator(
            7, (1, 2), (11,), True, (m.SectionEnd(1, 10, 1), m.SectionEnd(2, 20, 0))),))
        self.assertFalse(result.closed)
        self.assertFalse(result.pending)
        self.assertEqual((self.snapshot, self.component), before)

    def test_attached_tail_does_not_reverse_authored_member_order(self):
        for route in self.component.network.routes:
            route.span.start, route.span.end = route.span.end, route.span.start
            for index, end in enumerate((route.span.start, route.span.end)):
                end.endpoint = m.SectionEnd(route.span.street_id, route.span.section_ids[0], index)
            route.start_port, route.end_port = (
                S(endpoint=end.endpoint, junction_id=end.junction_id) if end.junction_id else None
                for end in (route.span.start, route.span.end))
        result = self.plan().free_ended[0]
        self.assertFalse(result.attached_at_head)
        self.assertEqual(result.member_street_ids, (1, 2))
        self.assertEqual(result.free_ends, (m.SectionEnd(1, 10, 0), m.SectionEnd(2, 20, 1)))

    def test_multiple_sections_retain_the_last_original_section_id(self):
        route = self.component.network.routes[0]
        self.snapshot.streets[1]['items'].append(dict(kind='Section', id=19))
        route.span.section_ids = (10, 19)
        route.span.end.endpoint = m.SectionEnd(1, 19, 1)
        self.assertEqual(self.plan().free_ended[0].free_ends[0], m.SectionEnd(1, 19, 1))

    def test_nonfree_or_missing_port_cannot_masquerade_as_free(self):
        for change in ('category', 'remote', 'fake_free', 'missing_port', 'wrong_port', 'wrong_section', 'partial'):
            with self.subTest(change=change):
                self.setUp()
                route = self.component.network.routes[0]
                if change == 'category': route.span.end.owner_kind = 'CATEGORY_BOUNDARY'
                if change == 'remote': route.span.end.owner_kind = 'REMOTE_LOGICAL'
                if change == 'fake_free': route.span.end.junction_id = 99
                if change == 'missing_port': route.start_port = None
                if change == 'wrong_port': route.start_port.endpoint = m.SectionEnd(1, 10, 1)
                if change == 'wrong_section': route.span.end.endpoint = m.SectionEnd(1, 99, 1)
                if change == 'partial': route.span.first_control = 1
                result = self.plan()
                self.assertFalse(result.free_ended)
                self.assertEqual(len(result.pending), 1)

    def test_members_must_attach_at_the_same_bundle_end_and_owner(self):
        for change in ('end', 'owner'):
            with self.subTest(change=change):
                self.setUp()
                if change == 'end':
                    self.snapshot.bundles[7]['topology']['members'][1]['forward'] = True
                else:
                    self.component.network.ordinary_site_ids = (11, 12)
                    route = self.component.network.routes[1]
                    route.span.end.junction_id = route.end_port.junction_id = 12
                self.assertFalse(self.plan().free_ended)
                self.assertIn('same attached', self.plan().pending[0]['reason'])

    def test_compound_terminal_keeps_all_actual_surface_sites(self):
        self.component.network.ordinary_site_ids = ()
        self.component.network.compounds = (S(site_ids=(11, 12)),)
        route = self.component.network.routes[1]
        route.span.end.junction_id = route.end_port.junction_id = 12
        self.assertEqual(self.plan().free_ended[0].attached_sites, (11, 12))

    def test_closed_dispatch_is_unchanged_when_free_end_branch_is_enabled(self):
        snapshot, component = fixture()
        expected = self.planner.prepare(snapshot, component)
        self.assertEqual(self.planner.prepare(snapshot, component, free_ends=True), expected)


if __name__ == '__main__':
    unittest.main()
