"""Median continuity through outside turnoffs and interruption by crossroads."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import side_junction_fixture_networks
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.ownership import CurbChainUnavailable, surface_key
from tested_python_streets.mesh import cross


def area(mesh):
    xy = mesh.xyz[mesh.quads, :2]
    return .5 * np.sum(cross(xy, np.roll(xy, -1, axis=1)))


def components(mesh):
    parent = np.arange(len(mesh.xyz))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for face in mesh.quads:
        first = root(int(face[0]))
        for i in face[1:]:
            parent[root(int(i))] = first
    return len({root(int(i)) for i in np.unique(mesh.quads)})


class PythonBundleCurbContinuations(unittest.TestCase):
    def test_curved_member_rims_continue_over_the_actual_rounded_junction_boundary(self):
        network = side_junction_fixture_networks(model)[0]
        sections = tuple(replace(s, points=(s.points[0], (s.points[1][0], 9.5, 0.), s.points[-1]))
                         if s.id in (1, 2) else s for s in network.sections)
        network = replace(network, sections=sections)
        result = PythonMeshBackend(network.settings).prepare(network)
        part = result.products[model.Key('separator', 7)]
        self.assertEqual(len(part.metadata['runs']), 1)
        self.assertEqual(components(part.mesh), 1)
        self.assertGreater(part.mesh.validate()['minimum_corner'], 0.)
        np.testing.assert_allclose(part.mesh.xyz[:, 2], .8, atol=1e-12)
        curve = result.products[surface_key(model.Key('intersection', 1))].cap((1, 1), (2, 0))
        self.assertGreater(float(np.ptp(curve[:, 1])), 1e-3)

    def test_the_outward_curb_rule_also_works_on_the_other_bundle_member(self):
        network = side_junction_fixture_networks(model)[0]
        mirrored = replace(network, sections=tuple(replace(s,
            points=tuple((x, -y, z) for x, y, z in s.points)) for s in network.sections),
            bundles=(replace(network.bundles[0], members=network.bundles[0].members[::-1]),))
        result = PythonMeshBackend(mirrored.settings).prepare(mirrored)
        part = result.products[model.Key('separator', 7)]
        self.assertEqual(len(part.metadata['runs']), 1)
        self.assertEqual(components(part.mesh), 1)
        self.assertAlmostEqual(area(part.mesh), 80. * 10.6, places=7)
        np.testing.assert_allclose(part.mesh.xyz[:, 2], .8, atol=1e-12)

    def test_outward_turnoff_preserves_the_complete_rectangular_median_and_its_plateau(self):
        network = side_junction_fixture_networks(model)[0]
        result = PythonMeshBackend(network.settings).prepare(network)
        part = result.products[model.Key('separator', 7)]
        self.assertEqual(len(part.metadata['runs']), 1)
        self.assertEqual(components(part.mesh), 1)
        self.assertAlmostEqual(area(part.mesh), 80. * (16. - 5.4), places=7)
        np.testing.assert_allclose(part.mesh.xyz[:, 2], .8, atol=1e-12)
        self.assertGreater(part.mesh.validate()['minimum_corner'], 0.)
        surface = result.products[surface_key(model.Key('intersection', 1))]
        curve = surface.cap((1, 1), (2, 0))
        np.testing.assert_allclose(curve[:, 1], 5.3, atol=1e-12)
        np.testing.assert_allclose(curve[:, 2], .8, atol=1e-12)

    def test_a_crossroad_interrupts_the_median_and_cannot_be_bridged_by_a_curb(self):
        network = side_junction_fixture_networks(model)[1]
        result = PythonMeshBackend(network.settings).prepare(network)
        part = result.products[model.Key('separator', 7)]
        self.assertEqual(len(part.metadata['runs']), 2)
        self.assertEqual(components(part.mesh), 2)
        self.assertTrue(np.all(abs(part.mesh.xyz[:, 0]) >= 2.7 - 1e-8))
        self.assertGreaterEqual(float(part.mesh.xyz[:, 2].min()), -1e-12)
        self.assertLessEqual(float(part.mesh.xyz[:, 2].max()), .8 + 1e-12)
        self.assertLess(area(part.mesh), 80. * (16. - 5.4) - 5.4 * (16. - 5.4))
        surface = result.products[surface_key(model.Key('intersection', 1))]
        with self.assertRaises(CurbChainUnavailable):
            surface.cap((1, 1), (2, 0))

    def test_turnoff_width_edit_matches_a_clean_rebuild_and_restores_exactly(self):
        network = side_junction_fixture_networks(model)[0]
        backend = PythonMeshBackend(network.settings)
        initial = backend.prepare(network)
        backend.commit(initial)
        changed = replace(network, sections=tuple(replace(s, width=s.width + .4) if s.id == 4 else s for s in network.sections))
        edited = backend.prepare(changed)
        self.assertIn(model.Key('separator', 7), edited.rebuilt)
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        for key, part in edited.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, fresh.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, fresh.products[key].mesh.quads)
        self.assertAlmostEqual(area(edited.products[model.Key('separator', 7)].mesh), 80. * 10.6, places=7)
        backend.commit(edited)
        restored = backend.prepare(network)
        for key, part in initial.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, restored.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, restored.products[key].mesh.quads)


if __name__ == '__main__':
    unittest.main()
