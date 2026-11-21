"""Unsupported grade continuations are handoffs, not missing road junctions."""
from dataclasses import replace
import json
from types import SimpleNamespace
import unittest
import numpy as np

import _python_street_package
from tested_python_streets.adapter import PreparedPythonAdapter, components
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.model import Connection, End, Key, Network, Settings
from tested_python_streets.routes import RouteTopology, route_key


def section(identity, a, b, *, bridge=False):
    return dict(kind='Section', id=identity, points=(a, b), width=6., offset=0.,
                category='residential', tags={'bridge': 'yes'} if bridge else {})


def endpoint(identity, end):
    return dict(terminal_kind='Section', terminal_id=identity, end=end)


def fixture(kind='physical', reverse=False):
    road = section(10, (-40., 0.), (0., 0.))
    bridge = section(20, (0., 0.), (40., 0.), bridge=True)
    if reverse:
        road['points'] = road['points'][::-1]
        bridge['points'] = bridge['points'][::-1]
    ends = [endpoint(10, 0 if reverse else 1), endpoint(20, 1 if reverse else 0)]
    streets = {1: dict(id=1, items=[road]), 2: dict(id=2, items=[bridge])}
    origins = {'physical': [], 'seams': []}
    if kind == 'physical':
        origins['physical'].append(dict(junction_id=7, ring=[dict(endpoints=[end]) for end in ends]))
    elif kind == 'seam':
        origins['seams'].append(dict(junction_id=7, endpoints=ends))
    else:
        streets = {1: dict(id=1, items=[road, bridge])}
    return SimpleNamespace(streets=streets, bundles={}, junction_origins=origins)


def convert(source):
    return PreparedPythonAdapter(('residential',),
        lambda tags: 'Unsupported bridge grade' if tags.get('bridge') == 'yes' else None,
        {'residential': {'radius': 2.}}).read(source, 'scope', Settings())


class ScopeBoundaryTests(unittest.TestCase):
    def test_directed_source_handoffs_are_preserved_and_roundtrip(self):
        for kind in ('physical', 'seam'):
            for reverse in (False, True):
                with self.subTest(kind=kind, reverse=reverse):
                    result = convert(fixture(kind, reverse))
                    self.assertFalse(result.blocked_sections)
                    self.assertEqual(len(tuple(components(result.network))), 1)
                    road, = result.network.sections
                    self.assertEqual(road.scope_ends, (0 if reverse else 1,))
                    self.assertEqual(road.selection_ends, ())
                    self.assertEqual(result.network.intersections, ())
                    self.assertEqual(result.network.connections, ())
                    record, = [r for r in result.issues if r['kind'] == 'ScopeBoundary']
                    self.assertEqual((record['section'], record['deferred_section']), (10, 20))
                    self.assertEqual(record['location'], (0., 0., 0.))
                    self.assertEqual(Network.from_dict(json.loads(json.dumps(result.network.to_dict()))), result.network)

    def test_explicit_street_sequence_also_records_the_grade_handoff(self):
        result = convert(fixture('sequence'))
        self.assertFalse(result.blocked_sections)
        self.assertEqual(result.network.section_index[10].scope_ends, (1,))
        self.assertEqual([r['source_kind'] for r in result.issues if r['kind'] == 'ScopeBoundary'],
                         ['StreetSequence'])

    def test_grade_handoff_does_not_change_the_carriageway_geometry(self):
        network = convert(fixture()).network
        # The finite full-width cross-section is unchanged. The difference is
        # the explicit topology contract used by approaching junction builders.
        ordinary = replace(network, sections=tuple(replace(s, scope_ends=()) for s in network.sections))
        first = PythonMeshBackend(network.settings).prepare(network).products[Key('road', 10)].mesh
        second = PythonMeshBackend(ordinary.settings).prepare(ordinary).products[Key('road', 10)].mesh
        np.testing.assert_array_equal(first.xyz, second.xyz)
        np.testing.assert_array_equal(first.quads, second.quads)
        self.assertLessEqual(first.validate()['maximum_edge_ratio'], 8.)
        candidate = RouteTopology(network).trace(End(10, 0))
        self.assertTrue(candidate.terminal_scope)
        self.assertFalse(candidate.terminal_selection)

    def test_scope_end_cannot_be_owned_by_a_second_junction_or_selection(self):
        network = convert(fixture()).network
        with self.assertRaisesRegex(ValueError, 'scope boundary'):
            replace(network, connections=(Connection(1, (End(10, 0), End(10, 1))),))
        with self.assertRaisesRegex(ValueError, 'scope boundary'):
            replace(network.sections[0], selection_ends=(1,))

    def test_ground_junction_stays_buildable_when_a_distant_end_meets_a_bridge(self):
        source = fixture()
        source.streets[1]['items'][0]['points'] = ((0., 0.), (40., 0.))
        source.streets[2]['items'][0]['points'] = ((40., 0.), (80., 0.))
        source.streets[3] = dict(id=3, items=[section(30, (0., 0.), (0., 40.))])
        source.streets[4] = dict(id=4, items=[section(40, (0., 0.), (-40., 0.))])
        source.junction_origins['physical'].append(dict(junction_id=6,
            ring=[dict(endpoints=[endpoint(sid, 0)]) for sid in (10, 30, 40)]))
        result = convert(source)
        self.assertFalse(result.blocked_sections)
        network, = components(result.network)
        self.assertEqual({s.id for s in network.sections}, {10, 30, 40})
        update = PythonMeshBackend(network.settings).prepare(network)
        self.assertFalse(update.products[route_key(End(10, 0))].terminal_free)
        self.assertIn(Key('junction_mesh', 6), update.products)
        self.assertEqual(update.products[Key('road', 10)].mesh.xyz[:, 0].max(), 40.)

    def test_invalid_or_displaced_partner_is_not_treated_as_a_handoff(self):
        for invalid in ('position', 'dimensions'):
            source = fixture()
            bridge = source.streets[2]['items'][0]
            if invalid == 'position':
                bridge['points'] = ((1., 0.), (40., 0.))
            else:
                bridge['tags'] = {}
                bridge['width'] = -1.
            result = convert(source)
            self.assertIn(10, result.blocked_sections)
            self.assertEqual(result.network.section_index[10].scope_ends, ())

    def test_multiarm_ground_projection_retains_every_deferred_incidence(self):
        source = fixture()
        source.streets[3] = dict(id=3, items=[section(30, (0., 0.), (0., 40.))])
        source.junction_origins['physical'][0]['ring'].append(dict(endpoints=[endpoint(30, 0)]))
        result = convert(source)
        self.assertFalse(result.blocked_sections)
        self.assertFalse(any(r['kind'] == 'ScopeBoundary' for r in result.issues))
        seam, = result.network.connections
        self.assertEqual(seam.ends, (End(10, 1), End(30, 0)))
        record, = [r for r in result.issues if r['kind'] == 'ScopeJunction']
        self.assertEqual(record['available'], [dict(section=10, end=1), dict(section=30, end=0)])
        self.assertEqual(record['deferred'][0]['section'], 20)
        self.assertTrue(record['deferred'][0]['reason'])
        update = PythonMeshBackend(result.network.settings).prepare(result.network)
        self.assertIn(Key('road', 10), update.products)
        self.assertNotIn(Key('road', 20), update.products)

    def test_three_ground_arms_keep_a_rounded_junction_without_an_elevated_arm(self):
        source = fixture()
        for sid, point in ((30, (0., 40.)), (40, (0., -40.))):
            source.streets[sid] = dict(id=sid, items=[section(sid, (0., 0.), point)])
            source.junction_origins['physical'][0]['ring'].append(dict(endpoints=[endpoint(sid, 0)]))
        result = convert(source)
        self.assertFalse(result.blocked_sections)
        update = PythonMeshBackend(result.network.settings).prepare(result.network)
        mesh = update.products[Key('junction_mesh', 7)].mesh
        self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
        self.assertGreater(len(mesh.quads), 0)
        self.assertEqual({s.id for s in result.network.sections}, {10, 30, 40})

    def test_only_deferred_arms_produce_no_spurious_junction_error(self):
        source = fixture()
        source.streets[1]['items'][0]['tags']['bridge'] = 'yes'
        result = convert(source)
        self.assertFalse(result.network.sections)
        self.assertTrue(all(r['kind'] == 'Section' for r in result.issues))

    def test_one_ground_arm_with_two_elevated_arms_is_an_explicit_scope_end(self):
        source = fixture()
        source.streets[3] = dict(id=3, items=[section(30, (0., 0.), (0., 40.), bridge=True)])
        source.junction_origins['physical'][0]['ring'].append(dict(endpoints=[endpoint(30, 0)]))
        result = convert(source)
        self.assertFalse(result.blocked_sections)
        self.assertEqual(result.network.section_index[10].scope_ends, (1,))
        record, = [r for r in result.issues if r['kind'] == 'ScopeBoundary']
        self.assertEqual({p['section'] for p in record['deferred']}, {20, 30})


if __name__ == '__main__':
    unittest.main()
