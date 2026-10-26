"""Higher valence and two-ended unions, including independent boundary checks."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import higher_valence_fixture_networks, two_merged_terminal_fixture
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_terminal_caps import SkewMergedCapPlan
from test_python_street_dimensional_bundles import meshes


def edge_data(mesh):
    edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    return np.unique(edges, axis=0, return_counts=True)


def distances(points, edges):
    a, delta = edges[:, 0], edges[:, 1]-edges[:, 0]
    lengths = np.sum(delta*delta, axis=1)
    result = []
    for start in range(0, len(points), 128):
        offset = points[start:start+128, None]-a
        t = np.clip(np.sum(offset*delta, axis=-1)/lengths, 0., 1.)
        result.extend(np.linalg.norm(offset-t[..., None]*delta, axis=-1).min(axis=1))
    return np.asarray(result)


class PythonExtendedFixtureTests(unittest.TestCase):
    def equal(self, first, second):
        self.assertEqual(first.keys(), second.keys())
        for key in first:
            np.testing.assert_array_equal(first[key][0], second[key][0])
            np.testing.assert_array_equal(first[key][1], second[key][1])

    def test_six_and_seven_arms_keep_whole_profiles_and_analytical_curbs(self):
        for original in higher_valence_fixture_networks(model):
            for mode in (0, 3):
                with self.subTest(network=original.uuid, mode=mode):
                    network = replace(original, sections=tuple(replace(s, interpolation=mode) for s in original.sections))
                    authored = network.to_dict()
                    result = PythonMeshBackend(network.settings).prepare(network)
                    self.assertEqual(network.to_dict(), authored)
                    mesh = result.products[model.Key('junction_mesh', 1)].mesh
                    unique, counts = edge_data(mesh)
                    self.assertEqual(len(mesh.xyz)-len(unique)+len(mesh.quads), 1)
                    self.assertLessEqual(counts.max(), 2)
                    self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
                    self.assertGreaterEqual(mesh.xyz[:, 2].min(), -1e-8)
                    self.assertLessEqual(mesh.xyz[:, 2].max(), max(s.rise for s in network.sections)*network.settings.exaggeration+1e-8)
                    rim = mesh.xyz[unique[counts == 1]]
                    for section in network.sections:
                        road = result.products[model.Key('road', section.id)].mesh
                        p = road.xyz[road.seams[0]]
                        for edge in np.stack((p[:-1], p[1:]), axis=1):
                            error = np.minimum(abs(rim-edge).max(axis=(1, 2)), abs(rim-edge[::-1]).max(axis=(1, 2)))
                            self.assertLess(float(error.min()), 1e-8)
                    boundary_ids = np.unique(unique[counts == 1])
                    for curve in result.products[model.Key('mouths', 1)].boundaries:
                        points = curve(np.linspace(0., 1., 2001))
                        self.assertLessEqual(float(distances(points, rim[:, :, :2]).max()), network.settings.curve_tolerance+1e-8)
                        reference = np.stack((points[:-1], points[1:]), axis=1)
                        at = boundary_ids[distances(mesh.xyz[boundary_ids, :2], reference) <= 1e-5]
                        self.assertGreaterEqual(len(at), 2)
                        self.assertLessEqual(float(abs(mesh.xyz[at, 2]).max()), 1e-8)

    def test_both_merged_caps_keep_annular_topology_profiles_and_furniture_queries(self):
        for mode in (0, 3):
            with self.subTest(mode=mode):
                original = two_merged_terminal_fixture(model)
                network = replace(original, sections=tuple(replace(s, interpolation=mode) for s in original.sections))
                result = PythonMeshBackend(network.settings).prepare(network)
                plan = result.products[model.Key('bundle_plan', 7)]
                self.assertIsInstance(plan, SkewMergedCapPlan)
                self.assertEqual(plan.merged_ends, (0, 1))
                mesh = result.products[model.Key('bundle_pavement', 7)].mesh
                self.assertLess(len(mesh.quads), 12000)
                self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
                unique, counts = edge_data(mesh)
                self.assertLessEqual(counts.max(), 2)
                self.assertEqual(len(mesh.xyz)-len(unique)+len(mesh.quads), 0)
                self.assertEqual(result.products[model.Key('separator', 7)].metadata['islands'], 1)
                for end, profile in plan.union_profiles.items():
                    cap = plan.cap_meshes[-1, end]
                    self.assertLess(float(np.linalg.norm(cap.xyz[:, None]-profile[None, :], axis=-1).min(axis=0).max()), 1e-8)
                    self.assertLessEqual(plan.union_cap_metadata[end]['maximum_body_join_longitudinal_slope'], .03+1e-8)
                for identity in (85, 86):
                    position = result.products[model.Key('placement', identity)].position
                    self.assertTrue(np.isfinite(position).all())
                    self.assertAlmostEqual(position[2], float(plan.height(np.asarray(position[:2]))), places=10)

    def test_two_end_local_update_matches_full_build_and_restores_anchors(self):
        for mode in (0, 3):
            with self.subTest(mode=mode):
                original = two_merged_terminal_fixture(model)
                network = replace(original, sections=tuple(replace(s, interpolation=mode) for s in original.sections))
                backend = PythonMeshBackend(network.settings)
                first = backend.prepare(network); backend.commit(first)
                saved = meshes(first)
                remote = first.products[model.Key('road', 99)]
                old_pose = first.products[model.Key('placement', 85)].position
                source = network.sections[0]
                points = list(source.points)
                points[1] = (points[1][0], points[1][1]+.02, points[1][2])
                changed = replace(network, sections=(replace(source, points=tuple(points), width=source.width+.05), *network.sections[1:]))
                local = backend.prepare(changed)
                self.equal(meshes(local), meshes(PythonMeshBackend(network.settings).prepare(changed)))
                self.assertIs(local.products[model.Key('road', 99)], remote)
                self.assertNotEqual(local.products[model.Key('placement', 85)].position, old_pose)
                backend.commit(local)
                restored = backend.prepare(network)
                self.equal(meshes(restored), saved)
                self.assertEqual(restored.products[model.Key('placement', 85)].position, old_pose)


if __name__ == '__main__':
    unittest.main(verbosity=2)
