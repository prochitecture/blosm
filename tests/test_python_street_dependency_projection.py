"""Exact consumer inputs, including composite heights outside the road task."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import locality_fixture_networks
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend


class PythonDependencyProjections(unittest.TestCase):
    def test_far_potential_mouth_does_not_rebuild_the_unchanged_road_or_separator(self):
        network, = locality_fixture_networks(model)
        backend = PythonMeshBackend(network.settings)
        initial = backend.prepare(network)
        backend.commit(initial)
        changed = replace(network, sections=tuple(replace(s, width=5.8) if s.id == 8 else s for s in network.sections))
        update = backend.prepare(changed)
        self.assertIn(model.Key('compound_plan', 0), update.dependencies[model.Key('bundle_plan', 7)])
        self.assertNotIn(model.Key('bundle_plan', 7), update.rebuilt)
        for key in (model.Key('road', 4), model.Key('separator', 7), model.Key('placement', 84)):
            self.assertNotIn(key, update.rebuilt)
            self.assertIs(initial.products[key], update.products[key])
        self.assertIn(model.Key('junction_mesh', 2), update.rebuilt)
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        for key, part in update.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, fresh.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, fresh.products[key].mesh.quads)

    def test_nearby_actual_mouth_still_rebuilds_its_road_and_separator(self):
        network, = locality_fixture_networks(model)
        backend = PythonMeshBackend(network.settings)
        backend.commit(backend.prepare(network))
        changed = replace(network, sections=tuple(replace(s, width=5.8) if s.id == 4 else s for s in network.sections))
        update = backend.prepare(changed)
        self.assertIn(model.Key('road', 1), update.rebuilt)
        self.assertIn(model.Key('separator', 7), update.rebuilt)
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        for key in (model.Key('road', 1), model.Key('separator', 7)):
            np.testing.assert_array_equal(update.products[key].mesh.xyz, fresh.products[key].mesh.xyz)
            np.testing.assert_array_equal(update.products[key].mesh.quads, fresh.products[key].mesh.quads)

    @staticmethod
    def free_bundle(separation):
        return model.Network('free_projection', (
            model.Section(1, 1, ((0., separation, 0.), (20., separation, 0.)), width=5.4),
            model.Section(2, 2, ((0., -separation, 0.), (20., -separation, 0.)), width=5.4)),
            bundles=(model.Bundle(7, (model.BundleMember((model.Traversal(1),)),
                                      model.BundleMember((model.Traversal(2),)))),),
            furniture=(model.Placement(8, 2, 10., side='LEFT', offset=.1),))

    def test_separated_member_edit_keeps_the_other_road_and_anchor(self):
        network = self.free_bundle(8.)
        backend = PythonMeshBackend(network.settings)
        initial = backend.prepare(network)
        backend.commit(initial)
        changed = replace(network, sections=(replace(network.sections[0], width=6.2), network.sections[1]))
        update = backend.prepare(changed)
        self.assertNotIn(model.Key('road', 2), update.rebuilt)
        self.assertNotIn(model.Key('placement', 8), update.rebuilt)
        self.assertIn(model.Key('separator', 7), update.rebuilt)

    def test_union_anchor_depends_on_the_composite_height_when_its_own_road_is_unchanged(self):
        network = self.free_bundle(2.)
        backend = PythonMeshBackend(network.settings)
        initial = backend.prepare(network)
        backend.commit(initial)
        changed = replace(network, sections=(replace(network.sections[0], width=6.2), network.sections[1]))
        update = backend.prepare(changed)
        key = model.Key('placement', 8)
        self.assertNotIn(model.Key('road', 2), update.rebuilt)
        self.assertIn(key, update.rebuilt)
        self.assertEqual(initial.products[key].position[:2], update.products[key].position[:2])
        self.assertNotEqual(initial.products[key].position[2], update.products[key].position[2])
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        self.assertEqual(update.products[key], fresh.products[key])


if __name__ == '__main__':
    unittest.main()
