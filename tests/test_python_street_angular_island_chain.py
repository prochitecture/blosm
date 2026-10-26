"""Multiple physical islands, whole quad interfaces, source edits and anchors."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_angular_fixtures import angular_island_chain_network
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_angular_curves import polyline_distance
from tested_python_streets.mesh import Mesh
from test_python_street_angular_bundles import buffers


def edges(mesh):
    unique, count = np.unique(np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)),
                                  axis=-1).reshape(-1, 2), axis=1), axis=0, return_counts=True)
    if int(count.max()) > 2:
        raise AssertionError('Nonmanifold quad edge')
    return unique, unique[count == 1]


def components(pairs):
    parent = {int(index): int(index) for index in np.unique(pairs)}
    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index
    for first, last in pairs:
        parent[root(int(first))] = root(int(last))
    return len({root(index) for index in parent})


class AngularIslandChainTests(unittest.TestCase):
    @staticmethod
    def equal(first, second):
        if first.keys() != second.keys():
            raise AssertionError('Different output owners')
        for key in first:
            np.testing.assert_array_equal(first[key][0], second[key][0])
            np.testing.assert_array_equal(first[key][1], second[key][1])

    def qualify(self, source, update, islands=3):
        plan = update.products[model.Key('bundle_plan', 13)]
        self.assertEqual(type(plan).__name__, 'AngularBundlePlan')
        self.assertEqual(len(plan.islands), islands)
        pavement = update.products[model.Key('bundle_pavement', 13)]
        separator = update.products[model.Key('separator', 13)]
        self.assertEqual(separator.metadata['islands'], islands)
        self.assertLessEqual(pavement.metadata['angular_longitudinal_bound'], .03+1e-8)
        self.assertLessEqual(pavement.metadata['physical_slopes']['maximum'], plan.maximum_slope+1e-8)
        published = []
        for part in (pavement, separator):
            self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
            mesh = Mesh(part.mesh.xyz.astype(np.float32), part.mesh.quads)
            self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
            self.assertEqual(len(np.unique(mesh.quads)), len(mesh.xyz))
            self.assertGreaterEqual(float(mesh.xyz[:, 2].min()), plan.datum-1e-7)
            self.assertLessEqual(float(mesh.xyz[:, 2].max()),
                                 plan.datum+source.settings.exaggeration*plan.rise+1e-7)
            published.append(mesh)
        paving, grass = published
        all_paving, paving_rim = edges(paving)
        all_grass, grass_rim = edges(grass)
        self.assertEqual(components(all_paving), 1)
        self.assertEqual(components(all_grass), islands)
        self.assertEqual(components(paving_rim), islands+1)
        self.assertEqual(components(grass_rim), islands)
        for mesh, rim in ((paving, paving_rim), (grass, grass_rim)):
            degrees = np.bincount(rim.ravel(), minlength=len(mesh.xyz))
            self.assertTrue(np.all(degrees[np.unique(rim)] == 2))
        def segments(mesh, rim):
            return {tuple(sorted((tuple(mesh.xyz[a]), tuple(mesh.xyz[b])))) for a, b in rim}
        self.assertTrue(segments(grass, grass_rim) <= segments(paving, paving_rim))
        self.assertTrue(np.all(grass.xyz[:, 2] == grass.xyz[0, 2]))
        # Audit entire physical edge spans, independently of semantic seam keys.
        boundary = edges(separator.mesh)[1]
        xy = separator.mesh.xyz[boundary, :2]
        query = np.concatenate([(1.-t)*xy[:, 0]+t*xy[:, 1] for t in np.linspace(0., 1., 9)])
        distances = np.stack([polyline_distance(query, curve.points) for curve in plan.field.curves])
        self.assertLessEqual(float(distances.min(axis=0).max()), source.settings.curve_tolerance+1e-8)
        if source.settings.resolution == 1:
            self.assertLess(len(pavement.mesh.quads)+len(separator.mesh.quads), 10000)
        return plan

    def test_exact_blender_precision_three_islands_and_complete_interfaces(self):
        source = angular_island_chain_network(model, 0, float32=True)
        before = source.to_dict()
        self.qualify(source, PythonMeshBackend(source.settings).prepare(source))
        self.assertEqual(source.to_dict(), before)

    def test_smooth_blender_precision_three_islands_and_complete_interfaces(self):
        source = angular_island_chain_network(model, 3, float32=True)
        before = source.to_dict()
        self.qualify(source, PythonMeshBackend(source.settings).prepare(source))
        self.assertEqual(source.to_dict(), before)

    def test_reversed_source_storage_retains_meshes(self):
        for mode in (0, 3):
            source = angular_island_chain_network(model, mode, float32=True)
            reversed_source = angular_island_chain_network(model, mode, float32=True, reversed_sources=(17, 29))
            first = PythonMeshBackend(source.settings).prepare(source)
            second = PythonMeshBackend(source.settings).prepare(reversed_source)
            self.qualify(reversed_source, second)
            self.equal(buffers(first), buffers(second))

    def test_real_vertex_and_width_updates_equal_full_and_restore_remote(self):
        for mode in (0, 3):
            with self.subTest(mode=mode):
                source = angular_island_chain_network(model, mode, float32=True)
                backend = PythonMeshBackend(source.settings)
                initial = backend.prepare(source)
                saved = buffers(initial)
                remote = initial.products[model.Key('road', 99)]
                backend.commit(initial)
                for action in ('vertex', 'width'):
                    sections = []
                    for section in source.sections:
                        if section.id == (29 if action == 'vertex' else 17):
                            if action == 'vertex':
                                points = list(section.points)
                                x, y, z = points[24]
                                points[24] = (float(np.float32(x+.1)), y, z)
                                section = replace(section, points=tuple(points))
                            else:
                                section = replace(section, width=float(np.float32(section.width+.1)))
                        sections.append(section)
                    changed = replace(source, sections=tuple(sections))
                    local = backend.prepare(changed)
                    self.qualify(changed, local)
                    full = PythonMeshBackend(changed.settings).prepare(changed)
                    self.equal(buffers(local), buffers(full))
                    self.assertIs(local.products[model.Key('road', 99)], remote)
                    self.assertNotIn(model.Key('placement', 999), local.rebuilt)
                    for placement in source.furniture:
                        key = model.Key('placement', placement.id)
                        self.assertEqual(local.products[key], full.products[key])
                    backend.commit(local)
                    restored = backend.prepare(source)
                    self.equal(saved, buffers(restored))
                    backend.commit(restored)

    def test_two_and_four_islands_are_independent_components(self):
        for islands in (2, 4):
            for mode in (0, 3):
                with self.subTest(islands=islands, mode=mode):
                    source = angular_island_chain_network(model, mode, float32=True, islands=islands)
                    self.qualify(source, PythonMeshBackend(source.settings).prepare(source), islands)


if __name__ == '__main__':
    unittest.main()
