"""Full external profiles when a junction owns an internal crown-rise seam."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import absorbed_transition_fixture_networks
from test_python_street_dimensional_bundles import meshes, reverse_sections
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.ownership import footprint_key, surface_key
from tested_python_streets.providers import FurnitureProvider


def profile_network(*, rise=.12, width=6., mode=3, datum=0., inner_rise=.08):
    n = absorbed_transition_fixture_networks(model)[0]
    return replace(n, sections=tuple(replace(s, interpolation=mode, datum=datum,
        rise=rise if s.id == 2 else inner_rise if s.id == 1 else s.rise,
        width=width if s.id == 2 else s.width) for s in n.sections))


class PythonAbsorbedProfileTests(unittest.TestCase):
    def qualify(self, network, *, curb_tolerance=1e-10):
        before = network.to_dict()
        result = PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(before, network.to_dict())
        layout = result.products[footprint_key(model.Key('intersection', 1))]
        self.assertEqual(layout.absorbed_transitions, (model.Key('connection', 6),))
        self.assertIsNone(result.products[model.Key('seam_mesh', 6)])
        self.assertIsNone(getattr(result.products[model.Key('road', 1)], 'mesh', None))
        pavement = result.products[model.Key('junction_mesh', 1)].mesh
        self.assertLessEqual(pavement.validate()['maximum_edge_ratio'], 8.)
        all_mouths = []
        for sid in (2, 3, 4):
            end = next(end for section, end in layout.mouths if section == sid)
            road = result.products[model.Key('road', sid)].mesh
            row = road.xyz[road.seams[end]]
            profile = result.products[model.Key('profile', sid)]
            q = np.linspace(-1., 1., len(row))
            np.testing.assert_allclose(row[:, 2], profile.datum
                + network.settings.exaggeration * profile.value(q), atol=1e-12, rtol=0.)
            distance = np.linalg.norm(row[:, None] - pavement.xyz[None], axis=-1)
            self.assertLess(float(distance.min(axis=1).max()), 1e-8)
            all_mouths.extend(row[:, :2])
        datum = network.section_index[2].datum
        ceiling = datum + network.settings.exaggeration * max(network.section_index[s].rise for s in (2, 3, 4))
        self.assertGreaterEqual(float(pavement.xyz[:, 2].min()), datum - 1e-8)
        self.assertLessEqual(float(pavement.xyz[:, 2].max()), ceiling + 1e-8)
        faces = pavement.quads
        edges = np.sort(np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        rim = pavement.xyz[np.unique(unique[counts == 1])]
        mouth_distance = np.linalg.norm(rim[:, None, :2] - np.asarray(all_mouths)[None], axis=-1).min(axis=1)
        curb = rim[mouth_distance > 1e-8]
        self.assertGreater(len(curb), 4)
        np.testing.assert_allclose(curb[:, 2], datum, atol=curb_tolerance, rtol=0.)
        for corner in layout.corners:
            if corner.center is not None:
                self.assertAlmostEqual(corner.radius, 2., places=10)
                distance = np.linalg.norm(rim[:, :2] - corner.center, axis=1)
                self.assertGreaterEqual(np.count_nonzero(abs(distance - 2.) < 1e-8), 3)
        context = type('Context', (), {'__getitem__': lambda self, k: result.products[k]})()
        for anchor in network.furniture[:2]:
            self.assertEqual(FurnitureProvider.compound_owner(context, anchor), model.Key('intersection', 1))
            position = np.asarray(result.products[model.Key('placement', anchor.id)].position)
            height = result.products[surface_key(model.Key('intersection', 1))].height(position[None, :2])[0]
            self.assertAlmostEqual(position[2], height + anchor.height, places=10)
        return result

    def test_lower_higher_flat_profiles_both_modes_widths_and_authoring_directions(self):
        for mode in (0, 3):
            for width in (6., 6.05):
                for rise in (0., .04, .12, .18):
                    for reversed_ids in (set(), {1}, {2}, {1, 2}):
                        with self.subTest(mode=mode, width=width, rise=rise, reverse=reversed_ids):
                            self.qualify(reverse_sections(profile_network(rise=rise, width=width, mode=mode), reversed_ids))

    def test_datum_rotation_resolution_and_consumed_high_inner_crown(self):
        for resolution in (1, 2, 4):
            n = profile_network(rise=.04, inner_rise=.18, datum=37.)
            angle = .713
            matrix = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
            sections = []
            for s in n.sections:
                p = np.asarray(s.points).copy()
                p[:, :2] = p[:, :2] @ matrix.T + (1234., -5678.)
                sections.append(replace(s, points=tuple(map(tuple, p))))
            # Use the production height guard for kilometre-translated inputs;
            # the origin-centred cases above retain the tighter zero-curb check.
            self.qualify(replace(n, sections=tuple(sections), settings=replace(n.settings, resolution=resolution)),
                         curb_tolerance=1e-8)

    def test_rise_update_release_full_equivalence_ownership_reversal_and_restoration(self):
        n = profile_network()
        backend = PythonMeshBackend(n.settings)
        first = backend.prepare(n)
        saved = meshes(first)
        original_anchors = {p.id: first.products[model.Key('placement', p.id)] for p in n.furniture}
        backend.commit(first)
        remote = first.products[model.Key('road', 99)]
        backend.release_published_meshes()
        changed = replace(n, sections=tuple(replace(s, rise=.18) if s.id == 2 else s for s in n.sections))
        update = backend.prepare(changed)
        self.assertIs(update.products[model.Key('road', 99)], remote)
        full = PythonMeshBackend(changed.settings).prepare(changed)
        for key in update.rebuilt:
            a, b = update.products[key], full.products[key]
            if getattr(a, 'mesh', None) is not None:
                np.testing.assert_array_equal(a.mesh.xyz, b.mesh.xyz)
                np.testing.assert_array_equal(a.mesh.quads, b.mesh.quads)
        backend.commit(update)
        backend.release_published_meshes()
        apart = replace(n, sections=tuple(replace(s, points=((0., 0., 0.), (6., 0., 0.), (12., 0., 0.))) if s.id == 1
            else replace(s, points=((12., 0., 0.), *s.points[1:])) if s.id == 2 else s for s in n.sections))
        separate = backend.prepare(apart)
        self.assertIsNotNone(separate.products[model.Key('seam_mesh', 6)])
        backend.commit(separate)
        restored = backend.prepare(n)
        for key, (xyz, quads) in saved.items():
            if getattr(restored.products[key], 'mesh', None) is not None:
                np.testing.assert_array_equal(restored.products[key].mesh.xyz, xyz)
                np.testing.assert_array_equal(restored.products[key].mesh.quads, quads)
        for identity, pose in original_anchors.items():
            self.assertEqual(restored.products[model.Key('placement', identity)], pose)


if __name__ == '__main__':
    unittest.main()
