import json
from types import SimpleNamespace
import unittest
import _python_street_package
from tested_python_streets.adapter import PreparedPythonAdapter, components
from tested_python_streets.model import Network, Settings


def section(identity, a, b):
    return dict(kind='Section', id=identity, points=(a, b), category='residential', width=6., offset=0., tags={})


def source():
    streets = {1: dict(id=1, items=[section(10, (0., 5.), (20., 5.)),
                                   section(11, (20., 5.), (40., 5.))]),
               2: dict(id=2, items=[section(20, (0., -5.), (20., -5.)),
                                   section(21, (20., -5.), (40., -5.))])}
    bundle = dict(id=7, topology=dict(status='MAPPED', interior_streets=[], members=[
        dict(street_id=1, head_order=0, forward=True, complete=True),
        dict(street_id=2, head_order=1, forward=False, complete=True)]))
    return SimpleNamespace(streets=streets, bundles={7: bundle}, junction_origins={'physical': []})


def adapter():
    return PreparedPythonAdapter(('residential',), lambda tags: 'Unsupported grade' if tags.get('bridge') else None,
                                 {'residential': {'radius': 2.}})


class PythonAdapterTests(unittest.TestCase):
    def test_unresolved_placement_is_reported_as_invalid_ground_input(self):
        original = source()
        original.bundles.clear()
        original.streets[1]['items'][0]['source_placement'] = dict(
            status='UNRESOLVED', reason='Variable placement needs an explicit dimensional transition')
        result = adapter().read(original, 'import', Settings())
        self.assertNotIn(10, {s.id for s in result.network.sections})
        self.assertTrue(any(i['id'] == 10 and 'placement' in i['reason'] for i in result.issues))

    def test_category_style_radii_are_preferred_and_an_override_remains_strict(self):
        original = source()
        original.bundles.clear()
        original.junction_origins['physical'] = [dict(junction_id=8, ring=[dict(endpoints=[
            dict(terminal_id=10, terminal_kind='Section', end=0),
            dict(terminal_id=20, terminal_kind='Section', end=0),
            dict(terminal_id=11, terminal_kind='Section', end=1)])])]
        network = adapter().read(original, 'import', Settings()).network
        site, = network.intersections
        self.assertTrue(all(site.corner_radius(end) == (2., 'STYLE') for end in site.ends))
        override = site.with_corner_radius(site.ends[0], 3.)
        self.assertEqual(override.corner_radius(site.ends[0]), (3., 'STRICT'))
        self.assertTrue(all(override.corner_radius(end) == (2., 'STYLE') for end in site.ends[1:]))

    def test_multi_section_bundle_roles_and_direction_are_retained(self):
        result = adapter().read(source(), 'import', Settings())
        self.assertEqual(result.issues, ())
        network = result.network
        bundle = network.bundles[0]
        self.assertEqual([t.section_id for t in bundle.members[0].path], [10, 11])
        self.assertEqual([t.section_id for t in bundle.members[1].path], [21, 20])
        self.assertTrue(all(not t.forward for t in bundle.members[1].path))
        self.assertEqual({s.street_id for s in network.sections}, {1, 2})
        self.assertEqual(len(network.connections), 2)
        self.assertEqual(len(tuple(components(network))), 1)
        reopened = Network.from_dict(json.loads(json.dumps(network.to_dict())))
        self.assertEqual(reopened, network)

    def test_unmapped_bundle_blocks_its_original_ground_ownership(self):
        for declaration in ('owned', 'terminals', 'street'):
            original = source()
            bundle = original.bundles[7]
            bundle['topology'] = dict(status='UNRESOLVED', reason='Original lateral order is unresolved')
            if declaration == 'owned': bundle['owned_streets'] = (1, 2)
            elif declaration == 'terminals': bundle.update(head=(1, 2), tail=(2, 1))
            else:
                for street in original.streets.values(): street['bundle_id'] = 7
            with self.subTest(declaration=declaration):
                result = adapter().read(original, 'import', Settings())
                self.assertEqual(result.blocked_sections, {10, 11, 20, 21})
                self.assertFalse(result.network.bundles)
                issue, = result.issues
                self.assertEqual(issue['code'], 'UNMAPPED_BUNDLE_MEMBERSHIP')
                self.assertEqual(issue['sections'], [10, 11, 20, 21])
                self.assertEqual(issue['reason'], bundle['topology']['reason'])

    def test_unmapped_bundle_entirely_outside_scope_needs_no_ground_block(self):
        original = source()
        original.bundles[7].update(topology=dict(status='UNRESOLVED'), owned_streets=(1, 2))
        for street in original.streets.values():
            for item in street['items']: item['tags']['bridge'] = 'yes'
        result = adapter().read(original, 'import', Settings())
        self.assertFalse(result.blocked_sections)
        self.assertFalse(result.network.sections)
        self.assertFalse(any(issue['kind'] == 'Bundle' for issue in result.issues))

    def test_regrouping_objects_does_not_change_logical_preparation(self):
        original = source()
        first = adapter().read(original, 'import', Settings()).network
        for street in original.streets.values():
            street['object_group'] = 'some other Blender object'
        second = adapter().read(original, 'import', Settings()).network
        self.assertEqual(first, second)

    def test_missing_required_ground_approach_blocks_its_neighbors(self):
        original = source()
        original.bundles.clear()
        original.streets[2]['items'][0]['tags'] = {'bridge': 'yes'}
        original.junction_origins['physical'] = [dict(junction_id=8, ring=[dict(endpoints=[
            dict(terminal_id=10, terminal_kind='Section', end=0),
            dict(terminal_id=20, terminal_kind='Section', end=0)])])]
        result = adapter().read(original, 'import', Settings())
        self.assertIn(10, result.blocked_sections)
        self.assertTrue(any(r['reason'] == 'Required source approaches absent' for r in result.issues))

    def test_absent_topology_does_not_connect_coincident_unrelated_sources(self):
        original = source()
        original.bundles.clear()
        original.streets[2]['items'][0]['points'] = ((0., 5.), (20., -5.))
        result = adapter().read(original, 'import', Settings())
        self.assertEqual(len(tuple(components(result.network))), 2)
        self.assertEqual(len(result.network.intersections), 0)


if __name__ == '__main__':
    unittest.main()
