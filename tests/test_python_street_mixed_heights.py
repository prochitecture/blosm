"""Sharp mixed joins reconcile heights without moving their profiles or curbs."""
from dataclasses import replace
import math
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.drainage import BoundaryDrainage
from tested_python_streets.mixed_heights import MixedHeightCorrection
from test_python_street_uniform_mixed import uniform_source
from test_python_street_dimensional_bundles import meshes, reverse_sections


class PythonMixedHeightTests(unittest.TestCase):
    def build(self, source):
        authored = source.to_dict()
        backend = PythonMeshBackend(source.settings)
        result = backend.prepare(source)
        self.assertEqual(source.to_dict(), authored)
        plan = result.products[model.Key('seam_layout', 8)].transition
        self.assertIsNotNone(plan.height_corrections)
        self.assertFalse(plan.height_corrections.flags.writeable)
        self.assertTrue(plan.certification['transverse_drainage'])
        mesh = plan.mesh(backend.patterns, 1.)
        original = replace(plan, height_corrections=None).mesh(backend.patterns, 1.)
        np.testing.assert_array_equal(mesh.xyz[:, :2], original.xyz[:, :2])
        np.testing.assert_array_equal(mesh.quads, original.quads)
        shape = (plan.rows + 1, plan.across + 1)
        grid = np.arange(len(mesh.xyz)).reshape(shape)
        fixed = np.unique(np.r_[grid[:2].ravel(), grid[-2:].ravel(),
            grid[:, (0, plan.across // 2, plan.across)].ravel()])
        np.testing.assert_array_equal(mesh.xyz[fixed], original.xyz[fixed])
        heights = mesh.xyz[:, 2].reshape(shape)
        self.assertGreaterEqual(float(np.diff(heights[:, :plan.across // 2 + 1], axis=1).min()), -1e-10)
        self.assertLessEqual(float(np.diff(heights[:, plan.across // 2:], axis=1).max()), 1e-10)
        outlets = np.zeros(len(mesh.xyz), bool)
        outlets[grid[:, (0, plan.across)].ravel()] = True
        self.assertLessEqual(BoundaryDrainage().fill(mesh, outlets)[1]['maximum_fill'], 1e-10)
        self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
        self.assertLessEqual(float(mesh.quad_slopes().max()), plan.certification['maximum_slope'] + 1e-8)
        self.assertLessEqual(plan.certification['maximum_source_entrance_longitudinal_bound'], .03 + 1e-8)
        for t in (0., 1.):
            q = np.linspace(-1., 1., 601)
            np.testing.assert_allclose(plan.surface(t, q)[1], plan.datum + .08 * (1. - q*q), atol=1e-14, rtol=0.)
        return backend, result, plan

    def test_both_signed_turns_and_mode_orders_at_all_resolutions(self):
        for resolution in (1, 2, 4):
            for angle in (-110., 110.):
                for modes in ((0, 3), (3, 0)):
                    with self.subTest(resolution=resolution, angle=angle, modes=modes):
                        self.build(uniform_source(length=30., angle=angle, modes=modes, resolution=resolution))

    def test_reversed_rotated_translated_source_and_common_datum(self):
        source = reverse_sections(uniform_source(length=30., angle=110., datum=37.), {1, 2})
        theta = .713
        rotation = np.array(((math.cos(theta), -math.sin(theta)), (math.sin(theta), math.cos(theta))))
        sections = []
        for section in source.sections:
            xyz = np.asarray(section.points).copy()
            xyz[:, :2] = xyz[:, :2] @ rotation.T + (1234., -5678.)
            sections.append(replace(section, points=tuple(map(tuple, xyz))))
        self.build(replace(source, sections=tuple(sections)))

    def test_local_full_edits_furniture_queries_release_and_exact_restoration(self):
        source = uniform_source(length=30., angle=110.)
        source = replace(source, furniture=(model.Placement(85, 1, 24., 'RIGHT', -2.),
                                           model.Placement(99, 99, 10., 'RIGHT', -1.)))
        backend, initial, plan = self.build(source)
        saved = meshes(initial)
        remote = initial.products[model.Key('road', 99)]
        furniture = initial.products[model.Key('placement', 85)]
        backend.commit(initial)
        second = source.sections[1]
        points = list(second.points)
        points[-1] = (points[-1][0] + .025, points[-1][1], points[-1][2])
        edited = replace(source, sections=(replace(source.sections[0], width=6.02),
            replace(second, points=tuple(points)), source.sections[2]))
        local = backend.prepare(edited)
        full = PythonMeshBackend(edited.settings).prepare(edited)
        for key, (xyz, quads) in meshes(local).items():
            np.testing.assert_array_equal(xyz, meshes(full)[key][0])
            np.testing.assert_array_equal(quads, meshes(full)[key][1])
        self.assertEqual(local.products[model.Key('placement', 85)], full.products[model.Key('placement', 85)])
        self.assertIs(local.products[model.Key('road', 99)], remote)
        backend.commit(local)
        restored = backend.prepare(source)
        for key, (xyz, quads) in meshes(restored).items():
            np.testing.assert_array_equal(xyz, saved[key][0])
            np.testing.assert_array_equal(quads, saved[key][1])
        self.assertEqual(restored.products[model.Key('placement', 85)], furniture)
        backend.commit(restored)
        product = restored.products[model.Key('seam_surface', 8)]
        xy = product.pavement.mesh.xyz[:, :2].copy()
        expected = product.pavement.mesh.xyz[:, 2].copy()
        np.testing.assert_allclose(product.height(xy), expected, atol=1e-9, rtol=0.)
        backend.release_published_meshes()
        np.testing.assert_allclose(product.height(xy), expected, atol=1e-9, rtol=0.)
        moved = replace(source, furniture=(replace(source.furniture[0], station=24.01), source.furniture[1]))
        query = backend.prepare(moved)
        fresh = PythonMeshBackend(moved.settings).prepare(moved)
        self.assertEqual(query.rebuilt, (model.Key('placement', 85),))
        self.assertEqual(query.products[model.Key('placement', 85)], fresh.products[model.Key('placement', 85)])

    def test_projection_exhaustion_preserves_last_committed_build(self):
        source = uniform_source()
        backend = PythonMeshBackend(source.settings)
        initial = backend.prepare(source)
        backend.commit(initial)
        manifest = backend.manifest()
        with patch.object(MixedHeightCorrection, 'max_iterations', 0):
            with self.assertRaisesRegex(ValueError, 'iteration budget'):
                backend.prepare(uniform_source(length=30., angle=110.))
        self.assertEqual(backend.manifest(), manifest)
        self.assertIs(backend.prepare(source).products[model.Key('road', 99)], initial.products[model.Key('road', 99)])


if __name__ == '__main__':
    unittest.main(verbosity=2)
