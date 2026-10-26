"""Visible separator incidence across fully consumed terminal Sections."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import consumed_cap_fixture_networks
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.separator_ports import PORTS, SeparatorPortIndex, VisiblePort


class PythonSeparatorPortTests(unittest.TestCase):
    def test_cap_uses_the_visible_port_without_rewriting_the_original_junction(self):
        network, = consumed_cap_fixture_networks(model)
        result = PythonMeshBackend(network.settings).prepare(network)
        endpoint = model.End(2, 1)
        port = result.products[PORTS].resolve(endpoint)
        self.assertEqual(port.site.id, 2)
        self.assertNotIn(endpoint, port.site.ends)
        self.assertIn(model.End(3, 1), port.site.ends)
        self.assertIsNone(result.products[PORTS].resolve(model.End(3, 1)))
        layout = result.products[port.site.layout_key]
        self.assertEqual(layout.approaches[port.approach_index].source.id, 3)
        self.assertIsNone(getattr(result.products[model.Key('road', 3)], 'mesh', None))
        part = result.products[model.Key('separator', 7)]
        self.assertEqual(len(part.metadata['runs']), 1)
        self.assertEqual(part.metadata['runs'][0]['attached_caps'], [False, True])
        self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)

    def test_reversed_authored_section_resolves_the_opposite_directed_end(self):
        network, = consumed_cap_fixture_networks(model)
        bundle, = network.bundles
        changed = replace(network, sections=tuple(replace(s, points=s.points[::-1]) if s.id == 2 else s
            for s in network.sections), bundles=(replace(bundle, members=(bundle.members[0],
                model.BundleMember((model.Traversal(2, False), model.Traversal(3))))),),
            connections=(model.Connection(3, (model.End(2, 0), model.End(3, 0))),), furniture=())
        result = PythonMeshBackend(changed.settings).prepare(changed)
        self.assertEqual(result.products[PORTS].resolve(model.End(2, 0)).site.id, 2)
        self.assertIsNone(result.products[PORTS].resolve(model.End(2, 1)))
        self.assertGreater(result.products[model.Key('separator', 7)].mesh.validate()['minimum_corner'], 0.)

    def test_visible_port_width_update_matches_full_build_and_restores_exactly(self):
        network, = consumed_cap_fixture_networks(model)
        backend = PythonMeshBackend(network.settings)
        initial = backend.prepare(network)
        backend.commit(initial)
        changed = replace(network, sections=tuple(replace(s, width=6.2) if s.street_id == 2 else s
                                                  for s in network.sections))
        actual = backend.prepare(changed)
        self.assertIn(model.Key('separator', 7), actual.rebuilt)
        self.assertNotIn(model.Key('road', 99), actual.rebuilt)
        self.assertIn(model.Key('placement', 83), actual.rebuilt)
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        for key, part in actual.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, fresh.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, fresh.products[key].mesh.quads)
        backend.commit(actual)
        restored = backend.prepare(network)
        for key, part in initial.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, restored.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, restored.products[key].mesh.quads)

    def test_remote_junction_changes_the_index_without_rebuilding_this_separator(self):
        network, = consumed_cap_fixture_networks(model)
        sections = tuple(replace(s, points=(s.points[0], s.points[1])) if s.id == 99 else s
                         for s in network.sections)
        sections += (model.Section(100, 99, ((0., 40., 0.), (40., 40., 0.))),
                     model.Section(101, 101, ((0., 40., 0.), (0., 60., 0.))))
        network = replace(network, sections=sections, intersections=(*network.intersections,
            model.Intersection(99, (model.End(99, 1), model.End(100, 0), model.End(101, 0)))))
        backend = PythonMeshBackend(network.settings)
        initial = backend.prepare(network)
        backend.commit(initial)
        changed = replace(network, sections=tuple(replace(s, width=6.4) if s.id == 101 else s
                                                  for s in network.sections))
        result = backend.prepare(changed)
        self.assertIn(PORTS, result.rebuilt)
        self.assertNotIn(model.Key('separator', 7), result.rebuilt)
        self.assertIs(result.products[model.Key('separator', 7)], initial.products[model.Key('separator', 7)])

    def test_ambiguous_visible_incidence_is_diagnostic(self):
        end = model.End(2, 1)
        first = model.Intersection(1, (model.End(1, 0), model.End(2, 1)))
        second = model.Intersection(2, (model.End(3, 0), model.End(4, 1)))
        index = SeparatorPortIndex({end: [VisiblePort(first, 0), VisiblePort(second, 0)]}, {})
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            index.resolve(end)
        self.assertIsNone(index.resolve(model.End(99, 0)))

    def test_local_signature_keeps_canonical_order_without_scanning_remote_ports(self):
        stamps = {model.End(5, 0): [('intersection:2', 'b', 'surface-b')],
                  model.End(999, 1): [('intersection:99', 'remote', None)],
                  model.End(3, 1): [('intersection:1', 'a', 'surface-a')]}
        index = SeparatorPortIndex({}, stamps)
        self.assertEqual(index.signature((5, 3, 5)), [
            (3, 1, [('intersection:1', 'a', 'surface-a')]),
            (5, 0, [('intersection:2', 'b', 'surface-b')])])


if __name__ == '__main__':
    unittest.main()
