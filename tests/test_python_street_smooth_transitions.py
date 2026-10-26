"""Sharp Smooth dimensional seams retain profiles and all update contracts."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from _python_street_fixtures import sharp_smooth_seam_fixture
from test_python_street_dimensional_bundles import meshes, reverse_sections
from test_python_street_fitted_transitions import network as sharp_network


class PythonSmoothTransitionTests(unittest.TestCase):
    def qualify(self, source):
        authored = source.to_dict()
        result = PythonMeshBackend(source.settings).prepare(source)
        self.assertEqual(source.to_dict(), authored)
        plan = result.products[model.Key('seam_layout', 8)].transition
        self.assertNotIn('exact_miter', plan.certification)
        self.assertTrue(any(c.curvature_intervals > 1 for c in plan.curbs))
        self.assertTrue(all(c.curvature_bound <= .5+1e-8 for c in plan.curbs))
        self.assertLessEqual(plan.certification['quality']['maximum_edge_ratio'], 8.)
        self.assertLessEqual(plan.certification['physical_slopes']['maximum'], plan.certification['maximum_slope']+1e-8)
        self.assertLessEqual(plan.certification['maximum_entrance_longitudinal_slope'], .03+1e-8)
        mesh = result.products[model.Key('seam_mesh', 8)].mesh
        for section in source.sections:
            self.assertLessEqual(result.products[model.Key('road', section.id)].mesh.validate()['maximum_edge_ratio'], 8.)
        q = np.linspace(-1., 1., plan.across+1)
        for i, approach in enumerate(plan.approaches):
            side = -q if i == 0 else q
            expected = approach.point(plan.distances[i], approach.offset+.5*approach.width*side)
            row = mesh.xyz[mesh.seams[i]]
            np.testing.assert_allclose(row[:, :2], expected, atol=1e-10, rtol=0.)
            np.testing.assert_allclose(row[:, 2], approach.profile.datum
                +source.settings.exaggeration*approach.profile.value(side), atol=1e-12, rtol=0.)
        dense = np.linspace(0., 1., 12001)
        for curb in plan.curbs:
            samples = curb(plan.parameters)
            chord = np.column_stack([np.interp(dense, plan.parameters, samples[:, i]) for i in (0, 1)])
            self.assertLessEqual(float(np.linalg.norm(curb(dense)-chord, axis=1).max()), source.settings.curve_tolerance)
        return result

    def equal(self, a, b):
        self.assertEqual(a.keys(), b.keys())
        for key in a:
            np.testing.assert_array_equal(a[key][0], b[key][0])
            np.testing.assert_array_equal(a[key][1], b[key][1])

    def test_both_turn_directions_keep_the_whole_curve_radius(self):
        source = sharp_smooth_seam_fixture(model)
        for sign in (-1., 1.):
            with self.subTest(sign=sign):
                sections = tuple(replace(s, points=tuple((x, sign*y, z) for x, y, z in s.points)) for s in source.sections)
                result = self.qualify(replace(source, sections=sections))
                self.assertEqual(len(result.products[model.Key('seam_mesh', 8)].mesh.quads), 960)

    def test_reverse_orientation_datum_resolution_and_zero_exaggeration(self):
        original = sharp_smooth_seam_fixture(model)
        for resolution in (1, 2, 4):
            with self.subTest(resolution=resolution):
                source = reverse_sections(original, {1, 2})
                source = replace(source, sections=tuple(replace(s, datum=37.) for s in source.sections),
                    settings=replace(source.settings, resolution=resolution, exaggeration=0. if resolution == 2 else 10.))
                self.qualify(source)

    def test_local_full_furniture_release_and_restoration(self):
        source = sharp_smooth_seam_fixture(model)
        backend = PythonMeshBackend(source.settings)
        first = backend.prepare(source)
        saved = meshes(first)
        pose = first.products[model.Key('placement', 85)]
        remote = first.products[model.Key('road', 99)]
        backend.commit(first)
        backend.release_published_meshes()
        sections = []
        for section in source.sections:
            if section.id == 1:
                section = replace(section, width=section.width+.02)
            elif section.id == 2:
                points = list(section.points)
                points[-1] = (points[-1][0]+.02, *points[-1][1:])
                section = replace(section, points=tuple(points))
            sections.append(section)
        changed = replace(source, sections=tuple(sections))
        local = backend.prepare(changed)
        self.assertIs(local.products[model.Key('road', 99)], remote)
        self.assertNotEqual(local.products[model.Key('placement', 85)], pose)
        # The retained remote mesh belongs to Blender after buffer release.
        # Combine its previously verified buffer with rebuilt numerical parts.
        self.equal({**saved, **meshes(local)}, meshes(PythonMeshBackend(source.settings).prepare(changed)))
        backend.commit(local)
        restored = backend.prepare(source)
        self.equal({**saved, **meshes(restored)}, saved)
        self.assertEqual(restored.products[model.Key('placement', 85)], pose)
        backend.commit(restored)
        backend.release_published_meshes()
        self.equal(meshes(backend.prepare(source, force=True)), saved)

    def test_unsafe_source_offsets_still_fail_before_surface_construction(self):
        source = sharp_network(length=30.)
        source = replace(source, sections=tuple(replace(s, interpolation=3) for s in source.sections))
        with self.assertRaisesRegex(ValueError, 'Smooth offset curvature could fold'):
            PythonMeshBackend(source.settings).prepare(source)

    def test_visible_road_sampling_uses_certified_trimmed_curvature(self):
        source = sharp_smooth_seam_fixture(model)
        result = self.qualify(source)
        for sid in (1, 2):
            with self.subTest(section=sid):
                path = result.products[model.Key('path', sid)]
                road = result.products[model.Key('road', sid)]
                first, last = road.metadata['start'], road.metadata['end']
                bound = path.visible_curvature_bound(first, last)
                self.assertLess(bound, path.curvature_bound)
                stations = path.samples(first, last, source.settings.spacing, local_curvature=True)
                self.assertLess(len(stations), len(path.samples(first, last, source.settings.spacing)))
                dense = np.linspace(first, last, 12001)
                index = np.clip(np.searchsorted(stations, dense)-1, 0, len(stations)-2)
                for sign in (-1., 1.):
                    offset = path.section.offset+sign*.5*path.section.width
                    points, samples = path.point(dense, offset), path.point(stations, offset)
                    a, edge = samples[index], samples[index+1]-samples[index]
                    f = np.clip(np.sum((points-a)*edge, axis=1)/np.sum(edge*edge, axis=1), 0., 1.)
                    error = np.linalg.norm(points-a-f[:, None]*edge, axis=1)
                    self.assertLessEqual(float(error.max()), source.settings.curve_tolerance)

    def test_shared_boundary_samples_match_pavement_after_buffer_release(self):
        from tested_python_streets.bundles import _road_rim
        for reversed_ids in (set(), {1, 2}):
            with self.subTest(reversed_ids=reversed_ids):
                source = reverse_sections(sharp_smooth_seam_fixture(model), reversed_ids)
                backend = PythonMeshBackend(source.settings)
                result = backend.prepare(source)
                context = type('Context', (), {'network': source, '__getitem__': lambda self,k: result.products[k]})()
                expected = {}
                for sid in (1, 2):
                    road = result.products[model.Key('road', sid)]
                    self.assertTrue(road.metadata['local_curvature_sampling'])
                    rows = road.mesh.xyz.reshape(-1, road.metadata['transverse_intervals']+1, 3)
                    for side in (-1, 1):
                        expected[sid, side] = rows[:, 0 if side == -1 else -1].copy()
                for released in (False, True):
                    if released:
                        backend.commit(result)
                        backend.release_published_meshes()
                    for (sid, side), points in expected.items():
                        rim = _road_rim(context, sid, side)
                        lo = int(np.argmin(np.linalg.norm(rim-points[0], axis=1)))
                        hi = int(np.argmin(np.linalg.norm(rim-points[-1], axis=1)))
                        np.testing.assert_array_equal(rim[lo:hi+1], points)


if __name__ == '__main__':
    unittest.main(verbosity=2)
