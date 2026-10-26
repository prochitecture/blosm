"""Complete exterior profiles and stable source ownership across compounds."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import compound_seam_fixture
from test_python_street_dimensional_bundles import meshes, reverse_sections
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.ownership import PLAN, surface_key
from tested_python_streets.providers import FurnitureProvider
from tested_python_streets.compound_ports import CompoundPortPlanner
from tested_python_streets.paths import Path


class PythonCompoundPortTests(unittest.TestCase):
    def qualify(self, network):
        original = network.to_dict()
        result = PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(original, network.to_dict())
        root = model.Key('intersection', 1)
        group = result.products[PLAN].by_site[root]
        self.assertIsNotNone(group.physical_layout)
        self.assertEqual(group.members, (root, model.Key('intersection', 2)))
        self.assertEqual({claim.section_id for claim in group.claims}, {1, 6})
        self.assertEqual(result.products[model.Key('seam_layout', 6)].absorbed_by, root)
        self.assertIsNone(result.products[model.Key('seam_mesh', 6)])
        self.assertIsNone(result.products[model.Key('junction_mesh', 2)])
        for sid in (1, 6):
            self.assertIsNone(getattr(result.products[model.Key('road', sid)], 'mesh', None))
        part = result.products[model.Key('junction_mesh', 1)]
        self.assertEqual(part.owner, model.Key('compound', 1))
        self.assertTrue(part.metadata['physical_exterior_ports'])
        self.assertEqual(part.metadata['member_junctions'], ['intersection:1', 'intersection:2'])
        self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
        xyz, quads = part.mesh.xyz, part.mesh.quads
        datum = network.section_index[2].datum
        ceiling = datum + network.settings.exaggeration * max(network.section_index[s].rise for s in (2, 3, 4, 5))
        self.assertGreaterEqual(float(xyz[:, 2].min()), datum - 1e-8)
        self.assertLessEqual(float(xyz[:, 2].max()), ceiling + 1e-8)
        edges, counts = np.unique(np.sort(np.stack((quads, np.roll(quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1),
                                  axis=0, return_counts=True)
        boundary = xyz[np.unique(edges[counts == 1])]
        mouths = []
        for sid in (2, 3, 4, 5):
            road = result.products[model.Key('road', sid)].mesh
            end = next(end for section, end in group.physical_layout.mouths if section == sid)
            row = road.xyz[road.seams[end]]
            distance = np.linalg.norm(row[:, None] - xyz[None], axis=-1)
            self.assertLessEqual(float(distance.min(axis=1).max()), 1e-8)
            ids = distance.argmin(axis=1)
            for first, last in zip(ids, ids[1:]):
                hits = np.flatnonzero(np.all(edges == sorted((int(first), int(last))), axis=1))
                self.assertEqual(len(hits), 1)
                self.assertEqual(counts[hits[0]], 1)
            profile = result.products[model.Key('profile', sid)]
            np.testing.assert_allclose(row[:, 2], profile.datum + network.settings.exaggeration
                * profile.value(np.linspace(-1., 1., len(row))), atol=1e-8, rtol=0.)
            mouths.extend(row)
        distance = np.linalg.norm(boundary[:, None] - np.asarray(mouths)[None], axis=-1).min(axis=1)
        curb = boundary[distance > 1e-8]
        self.assertGreater(len(curb), 4)
        np.testing.assert_allclose(curb[:, 2], datum, atol=1e-8, rtol=0.)
        context = type('Context', (), {'__getitem__': lambda self, key: result.products[key]})()
        for anchor in network.furniture:
            self.assertEqual(FurnitureProvider.compound_owner(context, anchor), root)
            position = np.asarray(result.products[model.Key('placement', anchor.id)].position)
            height = result.products[surface_key(root)].height(position[None, :2])[0]
            self.assertAlmostEqual(position[2], height + anchor.height, places=8)
        self.assertIn(PLAN, result.dependencies[model.Key('seam_layout', 6)])
        return result

    def test_modes_dimensions_rises_and_source_authoring_directions(self):
        for mode in (0, 3):
            for width in (6.05, 6.4):
                for rise in (.08, .12):
                    for ids in (set(), {1, 6}, {1, 2, 3, 4, 5, 6}):
                        with self.subTest(mode=mode, width=width, rise=rise, reversed=ids):
                            result = self.qualify(reverse_sections(compound_seam_fixture(model,
                                mode=mode, width=width, rise=rise), ids))
                            self.assertLess(len(result.products[model.Key('junction_mesh', 1)].mesh.quads), 1000)

    def test_radius_ownership_common_datum_and_resolution(self):
        n = compound_seam_fixture(model)
        first, second = n.intersections
        first = first.with_corner_radius(model.End(2, 0), 1.5)
        second = second.with_corner_radius(model.End(4, 0), 2.25)
        n = replace(n, intersections=(first, second))
        angle = .713
        matrix = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        for resolution in (1, 2, 4):
            sections = []
            for section in n.sections:
                p = np.asarray(section.points).copy()
                p[:, :2] = p[:, :2] @ matrix.T + (1234., -5678.)
                sections.append(replace(section, points=tuple(map(tuple, p)), datum=37.))
            transformed = replace(n, sections=tuple(sections), settings=replace(n.settings, resolution=resolution))
            with self.subTest(resolution=resolution):
                result = self.qualify(transformed)
                layout = result.products[PLAN].by_site[model.Key('intersection', 1)].physical_layout
                radii = {model.End(c.first.source.id, c.first.end): c.radius for c in layout.corners}
                self.assertEqual(radii[model.End(2, 0)], 1.5)
                self.assertEqual(radii[model.End(4, 0)], 2.25)

    def test_released_cache_update_full_equality_and_exact_restoration(self):
        n = compound_seam_fixture(model)
        engine = PythonMeshBackend(n.settings)
        first = engine.prepare(n)
        initial = meshes(first)
        engine.commit(first)
        remote = first.products[model.Key('road', 99)]
        engine.release_published_meshes()
        changed = replace(n, sections=tuple(replace(s, width=6.2, rise=.09) if s.id == 4 else s for s in n.sections))
        local = engine.prepare(changed)
        full = self.qualify(changed)
        self.assertIs(local.products[model.Key('road', 99)], remote)
        for key, (xyz, faces) in meshes(local).items():
            np.testing.assert_array_equal(xyz, full.products[key].mesh.xyz)
            np.testing.assert_array_equal(faces, full.products[key].mesh.quads)
        for anchor in n.furniture:
            key = model.Key('placement', anchor.id)
            self.assertEqual(local.products[key], full.products[key])
        engine.commit(local)
        engine.release_published_meshes()
        restored = engine.prepare(n)
        for key, (xyz, faces) in meshes(restored).items():
            np.testing.assert_array_equal(xyz, initial[key][0])
            np.testing.assert_array_equal(faces, initial[key][1])

    def test_unrelated_compound_does_not_invalidate_local_seam_or_surface(self):
        n = compound_seam_fixture(model)
        remote = compound_seam_fixture(model)
        sections = tuple(replace(s, id=s.id + 100, street_id=s.street_id + 100,
            points=tuple((x + 150., y + 150., z) for x, y, z in s.points)) for s in remote.sections if s.id != 99)
        junctions = tuple(replace(j, id=j.id + 100, ends=tuple(model.End(e.section_id + 100, e.end) for e in j.ends))
                          for j in remote.intersections)
        connections = tuple(replace(c, id=c.id + 100, ends=tuple(model.End(e.section_id + 100, e.end) for e in c.ends))
                            for c in remote.connections)
        n = replace(n, sections=(*n.sections, *sections), intersections=(*n.intersections, *junctions),
                    connections=(*n.connections, *connections))
        engine = PythonMeshBackend(n.settings)
        first = engine.prepare(n)
        engine.commit(first)
        changed = replace(n, sections=tuple(replace(s, width=6.2) if s.id == 104 else s for s in n.sections))
        update = engine.prepare(changed)
        self.assertIn(PLAN, update.rebuilt)
        for key in (model.Key('seam_layout', 6), model.Key('junction_mesh', 1), model.Key('mouths', 1), model.Key('mouths', 2)):
            self.assertNotIn(key, update.rebuilt)
            self.assertIs(first.products[key], update.products[key])

    def test_consumed_source_cannot_jump_a_concave_gap_or_exceed_support_budget(self):
        # The 10 cm notch lies strictly between 50 cm source samples. Both
        # segment endpoints are inside: whole-segment containment is required.
        section = model.Section(1, 1, ((1., 5., 0.), (5., 5., 0.)), width=.1, interpolation=0)
        network = model.Network('support_guard', (section,))
        path = Path(section)
        context = type('Context', (), {'network': network, '__getitem__': lambda self, key: path})()
        loop = type('Loop', (), {'xy': np.array(((0., 0.), (6., 0.), (6., 6.), (2.3, 6.),
                                               (2.3, 2.), (2.2, 2.), (2.2, 6.), (0., 6.)))})()
        with self.assertRaisesRegex(ValueError, 'exits'):
            CompoundPortPlanner(context, (), {}).check_support({1}, loop)
        path = Path(replace(section, points=((1., 5., 0.), (1000000., 5., 0.))))
        with self.assertRaisesRegex(ValueError, 'sample budget'):
            CompoundPortPlanner(context, (), {}).check_support({1}, loop)


if __name__ == '__main__':
    unittest.main()
