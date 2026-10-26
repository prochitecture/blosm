"""Fixed-outline caps, coordinated seams, source fidelity and bounded reuse."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_bundle_fixtures import bundle_networks
from test_python_street_bundles import boundary_edges
from test_python_street_surface_quality import context_for
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_axis import AxialPlan
from tested_python_streets.fibre_strips import CappedFibreStrips
from tested_python_streets.mesh import cross


def source(mode=3):
    return next(n for n in bundle_networks(model, interpolation=mode) if n.uuid.endswith('isolated_contact_probe'))


def distance(first, second):
    points = np.concatenate((first[:, 0], first[:, 1], first.mean(axis=1)))
    a, delta = second[:, 0], second[:, 1] - second[:, 0]
    result = 0.
    for start in range(0, len(points), 128):
        offset = points[start:start+128, None] - a
        t = np.clip(np.sum(offset * delta, axis=-1) / np.sum(delta * delta, axis=1), 0., 1.)
        result = max(result, float(np.linalg.norm(offset-t[..., None]*delta, axis=-1).min(axis=1).max()))
    return result


class PythonCappedFibreTests(unittest.TestCase):
    def test_bad_grass_chart_is_replaced_without_changing_its_outline_or_plateau(self):
        network = source()
        original = network.to_dict()
        context = context_for(network)
        plan = AxialPlan(context, network.bundles[0])
        previous = plan.separator_chart(context)
        result = plan.separator(context).mesh
        self.assertGreater(previous.validate()['maximum_edge_ratio'], 8.)
        self.assertLessEqual(result.validate()['maximum_edge_ratio'], 8.)
        self.assertLess(len(result.quads), .3 * len(previous.quads))
        old, new = boundary_edges(previous), boundary_edges(result)
        self.assertLess(distance(old[..., :2], new[..., :2]), 1e-9)
        self.assertLess(distance(new[..., :2], old[..., :2]), 1e-9)
        for mesh in (previous, result):
            self.assertTrue(np.all(mesh.xyz[:, 2] == .12 * network.settings.exaggeration))
        area = lambda mesh:.5 * np.sum(cross(mesh.xyz[mesh.quads, :2], np.roll(mesh.xyz[mesh.quads, :2], -1, axis=1)))
        self.assertAlmostEqual(area(previous), area(result), places=9)
        self.assertGreaterEqual(plan.separator_correction['minimum_ordinary_sine'], .03)
        self.assertLess(plan.separator_correction['physical_cusp_minimum_sine'], .03)
        self.assertEqual(network.to_dict(), original)

    def test_added_curb_rows_are_whole_shared_edges_of_the_pavement(self):
        network = source()
        update = PythonMeshBackend(network.settings).prepare(network)
        plan = update.products[model.Key('bundle_plan', 7)]
        grass = update.products[model.Key('separator', 7)].mesh
        pavement = update.products[model.Key('bundle_pavement', 7)].mesh
        self.assertLessEqual(pavement.validate()['maximum_edge_ratio'], 8.)
        edges = boundary_edges(grass)
        axial = (edges[..., :2] - plan.origin) @ plan.axis
        free = np.all(abs(axial - plan.first) < 1e-8, axis=1) | np.all(abs(axial - plan.last) < 1e-8, axis=1)
        others = boundary_edges(pavement)
        for edge in edges[~free]:
            delta = np.max(abs(others-edge), axis=(1, 2))
            reverse = np.max(abs(others[:, ::-1]-edge), axis=(1, 2))
            self.assertLess(float(np.minimum(delta, reverse).min()), 1e-8)

    def test_cache_release_reconstructs_the_same_chart_without_reinserting_rows(self):
        network = source()
        context = context_for(network)
        plan = AxialPlan(context, network.bundles[0])
        before = plan.separator(context).mesh
        peaks = plan.peak_a.points.copy(), list(plan.peak_a.keys)
        plan.release_geometry()
        self.assertIsNone(plan._separator_cache)
        after = plan.separator(context).mesh
        np.testing.assert_array_equal(after.xyz, before.xyz)
        np.testing.assert_array_equal(after.quads, before.quads)
        np.testing.assert_array_equal(plan.peak_a.points, peaks[0])
        self.assertEqual(plan.peak_a.keys, peaks[1])
        with self.assertRaisesRegex(ValueError, 'budget'):
            CappedFibreStrips(max_cells=10).build(context, plan)
        np.testing.assert_array_equal(plan.peak_a.points, peaks[0])

    def test_resolution_interpolation_and_datum_keep_valid_bounded_geometry(self):
        for mode in (0, 3):
            for resolution in (1, 2, 3, 4):
                with self.subTest(mode=mode, resolution=resolution):
                    network = source(mode)
                    network = replace(network, sections=tuple(replace(s, datum=37.) for s in network.sections),
                        settings=replace(network.settings, resolution=resolution, exaggeration=0. if resolution == 2 else 10.))
                    context = context_for(network)
                    plan = AxialPlan(context, network.bundles[0])
                    mesh = plan.separator(context).mesh
                    self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
                    self.assertLess(len(mesh.quads), 10000)
                    self.assertTrue(np.all(mesh.xyz[:, 2] == 37. + .12 * network.settings.exaggeration))


if __name__ == '__main__':
    unittest.main()
