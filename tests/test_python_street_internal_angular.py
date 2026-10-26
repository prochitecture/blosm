"""Curved Bundle separators interrupted by a transverse same-identity Street."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_angular_fixtures import internal_angular_state_network
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.mesh import Mesh
from test_python_street_angular_bundles import buffers
from test_python_street_angular_island_chain import edges, components


def network(mode=3, *, exaggeration=10., datum=0., resolution=1):
    source = internal_angular_state_network(model, mode)
    return replace(source, sections=tuple(replace(section,
        points=tuple(map(tuple, np.asarray(section.points, dtype=np.float32).astype(float))),
        rise=float(np.float32(section.rise)), datum=datum) for section in source.sections),
        furniture=(*source.furniture, model.Placement(904, 301, 5.)),
        settings=replace(source.settings, exaggeration=exaggeration, resolution=resolution))


class InternalAngularTests(unittest.TestCase):
    def qualify(self, source, result):
        self.assertEqual({section.id for section in source.sections if section.street_id == 800}, {300, 301, 302})
        self.assertEqual(source.bundles[0].interior_section_ids, (301,))
        grass = result.products[model.Key('separator', 13)]
        self.assertEqual(len(grass.metadata['runs']), 2)
        self.assertEqual(components(edges(grass.mesh)[0]), 2)
        paving = set()
        for key, part in result.products.items():
            if getattr(part, 'mesh', None) is None:
                continue
            self.assertLessEqual(part.mesh.slope_statistics(source.settings.exaggeration)['maximum'], .1+1e-8)
            published = Mesh(part.mesh.xyz.astype(np.float32), part.mesh.quads)
            self.assertLessEqual(published.validate()['maximum_edge_ratio'], 8.)
            self.assertGreaterEqual(float(published.xyz[:, 2].min()), source.sections[0].datum-1e-6)
            self.assertLessEqual(float(published.xyz[:, 2].max()),
                source.sections[0].datum+source.settings.exaggeration*source.sections[0].rise+1e-6)
            if key.kind == 'separator':
                continue
            paving.update(tuple(sorted((tuple(published.xyz[a]), tuple(published.xyz[b]))))
                          for a, b in edges(published)[1])
        xyz = grass.mesh.xyz.astype(np.float32)
        for a, b in edges(grass.mesh)[1]:
            self.assertIn(tuple(sorted((tuple(xyz[a]), tuple(xyz[b])))), paving)
        self.assertLess(len(grass.mesh.quads), 6000*source.settings.resolution**2)
        self.assertIn(model.Key('placement', 904), result.products)
        return grass

    def test_exact_complete_published_interfaces_and_source_street_identity(self):
        source = network(0)
        saved = source.to_dict()
        self.qualify(source, PythonMeshBackend(source.settings).prepare(source))
        self.assertEqual(source.to_dict(), saved)

    def test_smooth_complete_published_interfaces_and_source_street_identity(self):
        source = network(3)
        saved = source.to_dict()
        self.qualify(source, PythonMeshBackend(source.settings).prepare(source))
        self.assertEqual(source.to_dict(), saved)

    def test_combined_control_width_edit_equals_full_restore_and_keeps_remote(self):
        for mode in (0, 3):
            with self.subTest(mode=mode):
                source = network(mode)
                backend = PythonMeshBackend(source.settings)
                initial = backend.prepare(source)
                self.qualify(source, initial)
                original = buffers(initial)
                backend.commit(initial)
                sections = []
                for section in source.sections:
                    if section.id == 29:
                        points = np.asarray(section.points).copy()
                        points[len(points)//2, 0] += .1
                        section = replace(section, points=tuple(map(tuple, points)))
                    if section.id == 17:
                        section = replace(section, width=section.width+.1)
                    sections.append(section)
                changed = replace(source, sections=tuple(sections))
                actual = backend.prepare(changed)
                self.qualify(changed, actual)
                expected = PythonMeshBackend(changed.settings).prepare(changed)
                reference = buffers(expected)
                for key, (xyz, quads) in buffers(actual).items():
                    np.testing.assert_array_equal(xyz, reference[key][0])
                    np.testing.assert_array_equal(quads, reference[key][1])
                self.assertIs(actual.products[model.Key('road', 99)], initial.products[model.Key('road', 99)])
                self.assertNotIn(model.Key('placement', 999), actual.rebuilt)
                for key in actual.products:
                    if key.kind == 'placement':
                        self.assertEqual(actual.products[key], expected.products[key])
                backend.commit(actual)
                restored = backend.prepare(source)
                for key, (xyz, quads) in original.items():
                    np.testing.assert_array_equal(xyz, buffers(restored)[key][0])
                    np.testing.assert_array_equal(quads, buffers(restored)[key][1])

    def test_zero_exaggeration_retains_exact_shared_edges_on_translated_datum(self):
        source = network(3, exaggeration=0., datum=1.7)
        result = PythonMeshBackend(source.settings).prepare(source)
        self.qualify(source, result)
        self.assertTrue(all(np.all(part.mesh.xyz[:, 2] == 1.7) for part in result.products.values()
                            if getattr(part, 'mesh', None) is not None))


if __name__ == '__main__':
    unittest.main()
