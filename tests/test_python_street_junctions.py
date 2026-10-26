from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import fixture_networks
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend


class PythonJunctionTests(unittest.TestCase):
    def test_six_families_and_connected_bent_road_are_positive_quad_disks(self):
        for network in fixture_networks(model):
            with self.subTest(network=network.uuid):
                update = PythonMeshBackend(network.settings).prepare(network)
                for junction in network.intersections:
                    mesh = update.products[model.Key('junction_mesh', junction.id)].mesh
                    xyz, faces = mesh.xyz, mesh.quads
                    edges = np.sort(np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
                    unique, counts = np.unique(edges, axis=0, return_counts=True)
                    self.assertEqual(len(xyz) - len(unique) + len(faces), 1)
                    self.assertTrue(set(counts).issubset({1, 2}))
                    rim = unique[counts == 1]
                    degree = np.bincount(rim.ravel(), minlength=len(xyz))
                    self.assertTrue(np.all(degree[degree > 0] == 2))
                    self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
                    self.assertGreaterEqual(float(xyz[:, 2].min()), 0.)
                    self.assertLessEqual(float(xyz[:, 2].max()), .08 * network.settings.exaggeration + 1e-8)

    def test_round_curbs_have_zero_height_and_the_requested_radius(self):
        network = fixture_networks(model)[0]
        update = PythonMeshBackend(network.settings).prepare(network)
        layout = update.products[model.Key('mouths', 1)]
        for boundary in layout.boundaries:
            corner = boundary.corner
            if corner.center is None:
                continue
            distance = np.linspace(boundary.lengths[1], boundary.lengths[2], 20)
            xy = boundary(distance / boundary.length)
            np.testing.assert_allclose(np.linalg.norm(xy - corner.center, axis=1), network.intersections[0].radius, atol=1e-8)
        mesh = update.products[model.Key('junction_mesh', 1)].mesh
        # Band row zero is the physical curb. Locate it by analytical samples,
        # independently of the final welded mesh vertex numbering.
        for boundary in layout.boundaries:
            samples = boundary(np.linspace(0., 1., layout.rows + 1))
            for point in samples:
                index = np.argmin(np.linalg.norm(mesh.xyz[:, :2] - point, axis=1))
                self.assertLess(np.linalg.norm(mesh.xyz[index, :2] - point), 1e-8)
                self.assertAlmostEqual(mesh.xyz[index, 2], 0.)

    def test_full_entrance_profiles_match_road_meshes(self):
        network = fixture_networks(model)[4]
        update = PythonMeshBackend(network.settings).prepare(network)
        intersection = update.products[model.Key('junction_mesh', 1)].mesh
        for section in network.sections:
            road = update.products[model.Key('road', section.id)].mesh
            points = road.xyz[road.seams[0]]
            distances = np.linalg.norm(intersection.xyz[:, None, :] - points[None, :, :], axis=-1)
            self.assertLess(float(distances.min(axis=0).max()), 1e-8)

    def test_edit_updates_both_incident_junctions_and_preserves_identity(self):
        network = fixture_networks(model)[-1]
        engine = PythonMeshBackend(network.settings)
        engine.commit(engine.prepare(network))
        changed = replace(network, sections=(replace(network.sections[0],
            points=((0., 0., 0.), (20., 3., 0.), (40., 0., 0.))), *network.sections[1:]))
        incremental = engine.prepare(changed)
        full = PythonMeshBackend(changed.settings).prepare(changed)
        for jid in (1, 2):
            self.assertIn(model.Key('junction_mesh', jid), incremental.rebuilt)
        for key, part in incremental.products.items():
            if hasattr(part, 'mesh'):
                np.testing.assert_allclose(part.mesh.xyz, full.products[key].mesh.xyz, rtol=0., atol=0.)
                np.testing.assert_array_equal(part.mesh.quads, full.products[key].mesh.quads)

    def test_fixture_slopes_are_measured_without_inspection_exaggeration(self):
        for network in fixture_networks(model):
            with self.subTest(network=network.uuid):
                update = PythonMeshBackend(network.settings).prepare(network)
                for junction in network.intersections:
                    part = update.products[junction.mesh_key]
                    slopes = part.mesh.slope_statistics(network.settings.exaggeration)
                    self.assertLess(slopes['maximum'], .10)
                    # The same physical slope is recovered from the 1x mesh.
                    from tested_python_streets.mesh import Mesh
                    xyz = part.mesh.xyz.copy()
                    xyz[:, 2] /= network.settings.exaggeration
                    actual = Mesh(xyz, part.mesh.quads).slope_statistics()
                    self.assertAlmostEqual(slopes['maximum'], actual['maximum'])


if __name__ == '__main__':
    unittest.main()
