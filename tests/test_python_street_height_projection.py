import gc
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import weakref
from dataclasses import replace
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.slope_projection import SlopeProjection, corner_matrices
from tested_python_streets.surface_query import QuadHeightField
from tested_python_streets.ownership import CompoundProduct


def peak_grid(datum=0.):
    pattern = PatternLibrary().grid(2, 2)
    xyz = np.column_stack((2. * pattern.uv, np.full(len(pattern.uv), datum)))
    xyz[4, 2] += .08
    fixed = np.ones(len(xyz), bool)
    fixed[4] = False
    return Mesh(xyz, pattern.quads), fixed


class PythonHeightProjectionTests(unittest.TestCase):
    def test_optional_feasibility_stop_keeps_physical_constraints_and_default_projection(self):
        mesh, fixed = peak_grid()
        result, report = SlopeProjection().solve(mesh, fixed, datum=0., ceiling=.08,
            maximum_slope=.049, feasibility_slope=.05)
        self.assertEqual(report['convergence'], 'feasible_surface')
        self.assertLessEqual(result.quad_slopes().max(), .05)
        np.testing.assert_array_equal(result.xyz[fixed], mesh.xyz[fixed])
        self.assertGreaterEqual(float(result.xyz[:, 2].min()), 0.)
        self.assertLessEqual(float(result.xyz[:, 2].max()), .08)
        for bound in (.048, float('nan'), float('inf')):
            with self.assertRaisesRegex(ValueError, 'Feasibility slope'):
                SlopeProjection().solve(mesh, fixed, datum=0., ceiling=.08,
                    maximum_slope=.049, feasibility_slope=bound)

    def test_disabled_projection_cannot_be_reenabled_by_the_publication_guard(self):
        from tested_python_streets.publication_slopes import PublicationSlopeGuard
        mesh, fixed = peak_grid()
        mesh.xyz[:, :2] = mesh.xyz[:, :2] * .1 + (12345., -12345.)
        mesh.xyz[4, 2] = .01 / np.sqrt(2.)
        field = SimpleNamespace(datum=0., ceiling=.08, metadata={})
        self.assertLessEqual(float(mesh.quad_slopes().max()), .1 + 1e-10)
        rounded = PublicationSlopeGuard.physical_buffer(mesh, 0., 1.)
        self.assertGreater(float(rounded.quad_slopes().max()), .1 + 1e-7)
        with patch('tested_python_streets.publication_slopes.SlopeProjection.solve',
                   side_effect=AssertionError('Disabled projection was called')):
            with self.assertRaisesRegex(ValueError, 'projection is disabled'):
                PublicationSlopeGuard().correct(mesh, field, fixed, exaggeration=1., limit=.1,
                    max_cells=100, check_joins=lambda candidate: None, allow_projection=False)

    def test_one_free_peak_has_the_analytical_minimum_displacement_solution(self):
        original, fixed = peak_grid()
        expected = .05 / np.sqrt(2.)
        result, report = SlopeProjection().solve(original, fixed, datum=0., ceiling=.08, maximum_slope=.05)
        self.assertAlmostEqual(result.xyz[4, 2], expected, places=10)
        np.testing.assert_array_equal(result.xyz[fixed], original.xyz[fixed])
        np.testing.assert_array_equal(result.xyz[:, :2], original.xyz[:, :2])
        np.testing.assert_array_equal(result.quads, original.quads)
        self.assertEqual(original.xyz[4, 2], .08)
        self.assertLessEqual(result.quad_slopes().max(), .05 + 1e-9)
        self.assertLess(report['iterations'], 1000)

    def test_corner_constraints_match_the_independent_nine_point_slope_audit(self):
        mesh, fixed = peak_grid()
        mesh.xyz[4, :2] = .8, .9
        matrices = corner_matrices(mesh)
        values = np.linalg.norm(np.einsum('fcij,fj->fci', matrices, mesh.xyz[mesh.quads, 2]), axis=-1).max(axis=1)
        np.testing.assert_allclose(values, mesh.quad_slopes(), atol=1e-14, rtol=0.)

    def test_repeated_runs_and_common_datum_translation_are_deterministic(self):
        mesh, fixed = peak_grid()
        solver = SlopeProjection()
        a, _ = solver.solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.05)
        b, _ = solver.solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.05)
        np.testing.assert_array_equal(a.xyz, b.xyz)
        raised, fixed = peak_grid(12.5)
        c, _ = solver.solve(raised, fixed, datum=12.5, ceiling=12.58, maximum_slope=.05)
        np.testing.assert_allclose(a.xyz[:, 2], c.xyz[:, 2] - 12.5, atol=1e-10, rtol=0.)

    def test_feasible_data_retains_its_buffers_and_bad_fixed_data_is_rejected(self):
        mesh, fixed = peak_grid()
        unchanged, report = SlopeProjection(max_iterations=0).solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.2)
        self.assertIs(unchanged, mesh)
        self.assertEqual(report['iterations'], 0)
        with self.assertRaisesRegex(ValueError, 'infeasible'):
            SlopeProjection().solve(mesh, np.ones(len(fixed), bool), datum=0., ceiling=.08, maximum_slope=.05)
        with self.assertRaisesRegex(ValueError, 'iteration budget'):
            SlopeProjection(max_iterations=0).solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.05)

    def test_projection_budgets_and_boundary_envelopes_are_checked(self):
        for kwargs in ({'max_cells': 0}, {'max_iterations': -1}, {'memory_bytes': 0}, {'max_colors': 0}):
            with self.assertRaises(ValueError):
                SlopeProjection(**kwargs)
        mesh, fixed = peak_grid()
        with self.assertRaisesRegex(ValueError, 'memory budget'):
            SlopeProjection(memory_bytes=4096).solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.05)
        with self.assertRaisesRegex(ValueError, 'Boolean'):
            SlopeProjection().solve(mesh, np.flatnonzero(fixed), datum=0., ceiling=.08, maximum_slope=.05)
        mesh.xyz[0, 2] = -.01
        with self.assertRaisesRegex(ValueError, 'crown envelope'):
            SlopeProjection().solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.05)

    def test_a_boundary_gradient_equal_to_the_limit_is_not_mistaken_for_infeasibility(self):
        mesh, fixed = peak_grid()
        mesh.xyz[:, 2] = .05 * mesh.xyz[:, 0]
        mesh.xyz[4, 2] += .08
        result, _ = SlopeProjection().solve(mesh, fixed, datum=0., ceiling=.2, maximum_slope=.05)
        np.testing.assert_allclose(result.xyz[:, 2], .05 * result.xyz[:, 0], atol=1e-10, rtol=0.)


class PythonAdjustedSurfaceQueryTests(unittest.TestCase):
    def test_standalone_junction_anchor_queries_survive_publication_buffer_release(self):
        from _python_street_fixtures import fixture_networks
        from tested_python_streets import model
        from tested_python_streets.backend import PythonMeshBackend
        network = fixture_networks(model)[0]
        network = replace(network, furniture=(model.Placement(81, 0, 1., offset=-3.),))
        backend = PythonMeshBackend(network.settings)
        prepared = backend.prepare(network)
        before = prepared.products[model.Key('placement', 81)]
        self.assertGreater(before.position[2], .1)
        backend.commit(prepared)
        backend.release_published_meshes()
        changed = replace(network, furniture=(replace(network.furniture[0], station=1.2),))
        edited = backend.prepare(changed)
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        self.assertEqual(edited.rebuilt, (model.Key('placement', 81),))
        self.assertEqual(edited.products[model.Key('placement', 81)], fresh.products[model.Key('placement', 81)])

    def test_inverse_queries_match_bilinear_heights_on_varied_convex_quads(self):
        random = np.random.default_rng(731)
        for _ in range(12):
            xy = np.array(((0., 0.), (3., 0.), (3., 2.), (0., 2.))) + random.uniform(-.2, .2, (4, 2))
            height = random.uniform(0., .08, 4)
            mesh = Mesh(np.column_stack((xy, height)), np.array(((0, 1, 2, 3),)))
            field = QuadHeightField(mesh, datum=0., ceiling=.08)
            uv = np.vstack((random.uniform(0., 1., (40, 2)), ((0., 0.), (1., 0.), (1., 1.), (0., 1.))))
            u, v = uv.T
            weights = np.column_stack(((1. - u) * (1. - v), u * (1. - v), u * v, (1. - u) * v))
            np.testing.assert_allclose(field(weights @ xy), weights @ height, atol=1e-11, rtol=0.)
            offset = np.array((1_000_000., -2_000_000.))
            moved = QuadHeightField(Mesh(np.column_stack((xy + offset, height)), mesh.quads), datum=0., ceiling=.08)
            np.testing.assert_allclose(moved(weights @ xy + offset), weights @ height, atol=1e-9, rtol=0.)

    def test_queries_use_adjusted_heights_and_survive_release_and_exaggeration(self):
        original, fixed = peak_grid()
        adjusted, _ = SlopeProjection().solve(original, fixed, datum=0., ceiling=.08, maximum_slope=.05)
        field = QuadHeightField(adjusted, datum=0., ceiling=.08)
        expected = .05 / np.sqrt(2.)
        np.testing.assert_allclose(field(np.array(((1., 1.), (.75, .25)))), (expected, .75 * .25 * expected), atol=1e-10)
        reference = weakref.ref(adjusted)
        adjusted.xyz[:, 2] *= 10.
        del adjusted
        gc.collect()
        self.assertIsNone(reference())
        self.assertAlmostEqual(field(np.array(((1., 1.),)))[0], expected)

    def test_owner_excludes_separator_holes_before_strict_pavement_queries(self):
        pattern = PatternLibrary().grid(3, 3)
        mesh = Mesh(np.column_stack((pattern.uv, np.full(len(pattern.uv), .08))), np.delete(pattern.quads, 4, axis=0))
        field = QuadHeightField(mesh, datum=0., ceiling=.08)
        loops = (SimpleNamespace(hole=False, xy=np.array(((0., 0.), (1., 0.), (1., 1.), (0., 1.)))),
                 SimpleNamespace(hole=True, xy=np.array(((1/3, 1/3), (1/3, 2/3), (2/3, 2/3), (2/3, 1/3)))))
        owner = CompoundProduct(None, SimpleNamespace(loops=loops), field, {}, {}, {}, {}, None, None, 10.)
        np.testing.assert_allclose(owner.height(np.array(((.1, .1), (.5, .5), (2., 2.)))), (.8, .8, 0.), atol=1e-14)
        with self.assertRaisesRegex(ValueError, 'outside pavement'):
            field(np.array(((.5, .5),)))

    def test_invalid_queries_and_spatial_budgets_are_explicit(self):
        mesh, fixed = peak_grid()
        field = QuadHeightField(mesh, datum=0., ceiling=.08)
        self.assertEqual(field(np.empty((0, 2))).shape, (0,))
        for query in (np.ones((2, 3)), np.array(((float('nan'), 0.),))):
            with self.assertRaises(ValueError):
                field(query)
        with self.assertRaisesRegex(ValueError, 'outside pavement'):
            field(np.array(((4., 4.),)))
        with self.assertRaisesRegex(ValueError, 'spatial/memory budget'):
            QuadHeightField(mesh, datum=0., ceiling=.08, max_entries=1)


if __name__ == '__main__':
    unittest.main()
