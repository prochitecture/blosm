"""Attached overlaps preserve actual mouth profiles, caps and source ownership."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import attached_overlap_fixture_networks
from test_python_street_bundles import boundary_edges
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_ends import PathRange
from tested_python_streets.paths import Path
from tested_python_streets.providers import Part
from tested_python_streets.profiles import plateau_extension


def meshes(update):
    return {key: (p.mesh.xyz.copy(), p.mesh.quads.copy()) for key, p in update.products.items()
            if isinstance(p, Part) and p.mesh is not None}


class PythonAttachedUnionTests(unittest.TestCase):
    def source(self):
        return attached_overlap_fixture_networks(model)[0]

    def equal(self, first, second):
        self.assertEqual(first.keys(), second.keys())
        for key in first:
            np.testing.assert_array_equal(first[key][0], second[key][0])
            np.testing.assert_array_equal(first[key][1], second[key][1])

    def test_original_range_polynomials_and_controls_are_not_reinterpolated(self):
        for mode in (0, 3):
            source = replace(self.source().sections[0], interpolation=mode)
            path = Path(source)
            view = PathRange(path, 4.75, path.length - 4.75)
            stations = np.linspace(0., view.length, 23)
            for a, b in zip(view.frame(stations), path.frame(stations + 4.75)):
                np.testing.assert_array_equal(a, b)
            self.assertIs(view.section, source)

    def test_translated_plateau_rounding_does_not_hide_a_real_height_conflict(self):
        datum, height = 37., 1.2
        peak = datum + height
        rounded = peak
        for _ in range(4):
            rounded = np.nextafter(rounded, -np.inf)
        boundary = np.array(((0., 0., datum), (1., 0., peak), (1., 0., rounded),
                             (1., 1., peak), (0., 1., datum)))
        xy = np.array(((.5, .5), (1., 0.), (1., 1.)))
        values, _ = plateau_extension(xy, boundary, datum, height)
        self.assertTrue(np.all(np.isfinite(values)))
        self.assertLessEqual(float(values.max()), peak)
        self.assertGreaterEqual(float(values.min()), datum)
        conflict = boundary.copy()
        conflict[2, 2] -= 1e-5
        with self.assertRaisesRegex(ValueError, 'conflicting heights'):
            plateau_extension(xy, conflict, datum, height)

    def test_body_and_rounded_caps_have_whole_interfaces_and_bounded_profiles(self):
        network = self.source()
        original = network.to_dict()
        update = PythonMeshBackend(network.settings).prepare(network)
        plan = update.products[model.Key('bundle_plan', 7)]
        self.assertEqual(plan.mode, 'union')
        self.assertEqual(set(plan.ranges), {1, 2})
        self.assertNotIn(model.Key('road', 1), meshes(update))
        pavement = update.products[model.Key('bundle_pavement', 7)].mesh
        grass = update.products[model.Key('separator', 7)]
        self.assertLessEqual(pavement.validate()['maximum_edge_ratio'], 8.)
        self.assertLessEqual(grass.mesh.validate()['maximum_edge_ratio'], 8.)
        self.assertEqual(grass.metadata['welded_cap_body_vertices'], 34)
        self.assertEqual(grass.metadata['islands'], 2)
        self.assertLess(grass.mesh.quads.shape[0], 3000)
        self.assertLess(pavement.slope_statistics(10.)['maximum'], .1)
        junctions = np.concatenate([boundary_edges(p.mesh) for key, p in update.products.items()
            if key.kind == 'junction_mesh' and getattr(p, 'mesh', None) is not None])
        matched = 0
        for edge in boundary_edges(pavement):
            axial = (edge[:, :2] - plan.origin) @ plan.axis
            if not (np.all(abs(axial-plan.first)<1e-8) or np.all(abs(axial-plan.last)<1e-8)):
                continue
            error = np.minimum(np.max(abs(junctions-edge), axis=(1, 2)), np.max(abs(junctions[:, ::-1]-edge), axis=(1, 2)))
            self.assertLess(float(error.min()), 1e-8)
            matched += 1
        self.assertEqual(matched, 64)
        self.assertTrue(update.products[model.Key('road', 1)].metadata['junction_spans'])
        for identity in (81, 82, 83):
            self.assertTrue(np.all(np.isfinite(update.products[model.Key('placement', identity)].position)))
        self.assertEqual(network.to_dict(), original)

    def test_position_width_and_coordinated_radius_edits_match_full_build_and_restore(self):
        network = self.source()
        engine = PythonMeshBackend(network.settings)
        initial = engine.prepare(network)
        engine.commit(initial)
        saved = meshes(initial)
        remote = initial.products[model.Key('road', 99)]
        points = list(network.sections[0].points)
        points[3] = (0., 2.2, 0.)
        edits = (replace(network, sections=(replace(network.sections[0], points=tuple(points)), *network.sections[1:])),
                 replace(network, sections=(replace(network.sections[0], width=6.2), *network.sections[1:])),
                 replace(network, intersections=tuple(replace(j, radius=3.) if j.id in (1, 2) else j for j in network.intersections)))
        for changed in edits:
            update = engine.prepare(changed)
            fresh = PythonMeshBackend(changed.settings).prepare(changed)
            self.equal(meshes(update), meshes(fresh))
            self.assertIs(update.products[model.Key('road', 99)], remote)
            self.assertNotIn(model.Key('placement', 99), update.rebuilt)
            engine.commit(update)
            restored = engine.prepare(network)
            self.equal(meshes(restored), saved)
            engine.commit(restored)

    def test_union_separated_switch_and_released_cache_rebuild_preserve_source_ids(self):
        network = self.source()
        engine = PythonMeshBackend(network.settings)
        initial = engine.prepare(network)
        engine.commit(initial)
        saved = meshes(initial)
        separated = replace(network, sections=tuple(replace(s, points=tuple(
            (x, (5. if s.id == 1 else -5.) if i == 3 else y, z) for i, (x, y, z) in enumerate(s.points)))
            if s.id in (1, 2) else s for s in network.sections))
        update = engine.prepare(separated)
        self.assertEqual(update.products[model.Key('bundle_plan', 7)].mode, 'union')
        self.assertEqual(update.products[model.Key('separator', 7)].metadata['islands'], 1)
        pavement = update.products[model.Key('bundle_pavement', 7)].mesh
        edges = np.sort(np.stack((pavement.quads, np.roll(pavement.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        self.assertEqual(len(pavement.xyz)-len(np.unique(edges, axis=0))+len(pavement.quads), 2)
        self.assertNotIn(model.Key('road', 1), meshes(update))
        self.equal(meshes(update), meshes(PythonMeshBackend(separated.settings).prepare(separated)))
        engine.commit(update)
        restored = engine.prepare(network)
        self.equal(meshes(restored), saved)
        engine.commit(restored)
        engine.release_published_meshes()
        self.assertIsNone(engine.updates.products[model.Key('bundle_plan', 7)]._pavement_part)
        rebuilt = engine.prepare(network, force=True)
        self.equal(meshes(rebuilt), saved)

    def test_modes_resolutions_reverse_member_and_common_datum_keep_valid_surfaces(self):
        for mode in (0, 3):
            for resolution in (1, 2, 3, 4):
                for reverse in (False, True):
                    with self.subTest(mode=mode, resolution=resolution, reverse=reverse):
                        network = self.source()
                        sections = tuple(replace(s, interpolation=mode, datum=37.) for s in network.sections)
                        network = replace(network, sections=sections, settings=replace(network.settings,
                            resolution=resolution, exaggeration=0. if resolution == 2 else 10.))
                        if reverse:
                            lower = sections[1]
                            length = Path(lower).length
                            network = replace(network,
                                sections=(sections[0], replace(lower, points=lower.points[::-1]), *sections[2:]),
                                intersections=tuple(replace(j, ends=tuple(model.End(e.section_id, 1-e.end)
                                    if e.section_id == 2 else e for e in j.ends)) for j in network.intersections),
                                bundles=(replace(network.bundles[0], members=(network.bundles[0].members[0],
                                    model.BundleMember((model.Traversal(2, False),)))),),
                                furniture=tuple(replace(p, station=length-p.station, side='RIGHT')
                                    if p.section_id == 2 else p for p in network.furniture))
                        result = PythonMeshBackend(network.settings).prepare(network)
                        self.assertEqual(result.products[model.Key('bundle_plan', 7)].mode, 'union')
                        for xyz, faces in meshes(result).values():
                            self.assertEqual(faces.shape[1], 4)
                            self.assertGreaterEqual(float(xyz[:, 2].min()), 37.-1e-8)
                            self.assertLessEqual(float(xyz[:, 2].max()), 37.+.12*network.settings.exaggeration+1e-8)


if __name__ == '__main__':
    unittest.main()
