"""Full normal backend, whole interfaces and explicit Update for merged runs."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_merged_entrances import network
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.mesh import Mesh
from tested_python_streets.ownership import surface_key
from tested_python_streets.bundle_merged_plan import MergedBundlePlan


def source_for(mode, exaggeration=10.):
    source = network(mode)
    return replace(source, settings=replace(source.settings, exaggeration=exaggeration),
                   furniture=(*source.furniture, model.Placement(904, 301, 2., offset=-2.5)))


def snapshot(result):
    return {key: (part.mesh.xyz.copy(), part.mesh.quads.copy()) for key, part in result.products.items()
            if getattr(part, 'mesh', None) is not None}


def edges(mesh):
    e = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    unique, counts = np.unique(e, axis=0, return_counts=True)
    return unique, unique[counts == 1], counts


class MergedBundleTests(unittest.TestCase):
    def qualify(self, source, result):
        plan = result.products[model.Key('bundle_plan', 13)]
        self.assertIsInstance(plan, MergedBundlePlan)
        self.assertEqual(set(plan.source_ids), {17, 29, 117, 129})
        self.assertEqual({s.id for s in source.sections if s.street_id == 800}, {300, 301, 302})
        self.assertEqual(source.bundles[0].interior_section_ids, (301,))
        for sid in plan.source_ids:
            road = result.products[model.Key('road', sid)]
            self.assertIsNone(getattr(road, 'mesh', None))
            self.assertEqual(road.metadata['composite_owner'], 'bundle:13')
        consumed = result.products[model.Key('road', 301)]
        self.assertIsNone(getattr(consumed, 'mesh', None))
        self.assertEqual(consumed.metadata['consumed_by'], ['intersection:500'])
        self.assertIn(model.Key('placement', 904), result.products)
        product = result.products[surface_key(model.Key('intersection', 500))]
        item = result.products[model.Key('placement', 904)]
        self.assertAlmostEqual(item.position[2], product.height(np.array([item.position[:2]]))[0], places=10)
        saved, pavement, grass_edges = [], set(), []
        datum = source.sections[0].datum
        for key, part in result.products.items():
            if getattr(part, 'mesh', None) is None:
                continue
            mesh = Mesh(part.mesh.xyz.astype(np.float32), part.mesh.quads)
            unique, boundary, counts = edges(mesh)
            self.assertLessEqual(counts.max(), 2)
            self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
            self.assertGreaterEqual(mesh.xyz[:, 2].min(), datum-1e-6)
            self.assertLessEqual(mesh.xyz[:, 2].max(), datum+source.settings.exaggeration*.080001+1e-6)
            if source.settings.exaggeration:
                physical = Mesh(mesh.xyz.copy(), mesh.quads)
                physical.xyz[:, 2] = datum+(physical.xyz[:, 2]-datum)/source.settings.exaggeration
                self.assertLessEqual(physical.slope_statistics()['maximum'], .1+1e-6)
            interface = {tuple(sorted((tuple(mesh.xyz[a]), tuple(mesh.xyz[b])))) for a, b in boundary}
            if key.kind == 'separator':
                grass_edges.extend(interface)
            else:
                pavement.update(interface)
            if key != model.Key('road', 99):
                saved.append(mesh)
        self.assertTrue(all(edge in pavement for edge in grass_edges))
        points = np.vstack([mesh.xyz for mesh in saved])
        offset = np.cumsum([0, *[len(mesh.xyz) for mesh in saved]])
        points, inverse = np.unique(points, axis=0, return_inverse=True)
        faces = np.vstack([inverse[offset[i]:offset[i+1]][mesh.quads] for i, mesh in enumerate(saved)])
        whole = Mesh(points, faces)
        unique, boundary, counts = edges(whole)
        self.assertLessEqual(counts.max(), 2)
        self.assertEqual(len(points)-len(unique)+len(faces), 1)
        degree = np.bincount(boundary.ravel(), minlength=len(points))
        self.assertTrue(np.all(degree[np.unique(boundary)] == 2))
        parents = np.arange(len(points))
        def root(at):
            while parents[at] != at:
                parents[at] = parents[parents[at]]
                at = parents[at]
            return at
        for a, b in unique:
            parents[root(a)] = root(b)
        self.assertEqual(len({root(at) for at in range(len(points))}), 1)
        total = sum(len(result.products[model.Key(kind, 13)].mesh.quads)
                    for kind in ('bundle_pavement', 'separator'))
        self.assertLess(total, 8000)

    def test_normal_backend_generates_connected_bundles_in_both_modes(self):
        for mode in (0, 3):
            source = source_for(mode)
            original = source.to_dict()
            with self.subTest(mode=mode):
                result = PythonMeshBackend(source.settings).prepare(source)
                self.qualify(source, result)
                self.assertEqual(source.to_dict(), original)

    def test_update_full_restore_and_unrelated_identity(self):
        source = source_for(3)
        engine = PythonMeshBackend(source.settings)
        initial = engine.prepare(source)
        old = snapshot(initial)
        engine.commit(initial)
        sections = []
        for section in source.sections:
            points = np.asarray(section.points).copy()
            if section.id == 29:
                points[len(points)//2, 0] += .025
            sections.append(replace(section, points=tuple(map(tuple, points)),
                width=section.width+.1 if section.id == 17 else section.width))
        edited = replace(source, sections=tuple(sections))
        local = engine.prepare(edited)
        self.qualify(edited, local)
        full = PythonMeshBackend(edited.settings).prepare(edited)
        for key, (xyz, quads) in snapshot(local).items():
            np.testing.assert_array_equal(xyz, full.products[key].mesh.xyz)
            np.testing.assert_array_equal(quads, full.products[key].mesh.quads)
        for item in edited.furniture:
            key = model.Key('placement', item.id)
            self.assertEqual(local.products[key], full.products[key])
        self.assertIs(local.products[model.Key('road', 99)], initial.products[model.Key('road', 99)])
        self.assertIs(local.products[model.Key('placement', 999)], initial.products[model.Key('placement', 999)])
        engine.commit(local)
        restored = engine.prepare(source)
        for key, (xyz, quads) in old.items():
            np.testing.assert_array_equal(xyz, restored.products[key].mesh.xyz)
            np.testing.assert_array_equal(quads, restored.products[key].mesh.quads)

    def test_zero_exaggeration_uses_preserved_physical_caps(self):
        source = source_for(3, 0.)
        result = PythonMeshBackend(source.settings).prepare(source)
        self.qualify(source, result)
        for xyz, quads in snapshot(result).values():
            np.testing.assert_array_equal(xyz[:, 2], 0.)


if __name__ == '__main__':
    unittest.main()
