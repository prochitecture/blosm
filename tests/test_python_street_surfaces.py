from dataclasses import dataclass, replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import compound_cross_network, short_junction_pair
from test_python_street_partition import inputs
from tested_python_streets import model
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.conditions import BoundaryData, PlateauEdge
from tested_python_streets.mesh import Mesh
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.refinement import QuadRefiner
from tested_python_streets.surface_fields import BoundaryJetField
from tested_python_streets.validation import validate_overlaps


@dataclass
class LinearRule:
    first: float
    last: float

    def sample(self, t):
        t = np.asarray(t)
        return self.first + t * (self.last - self.first), np.tile((.04, 0.), t.shape + (1,))


def boundary_data():
    xy = np.array(((0., 0.), (1., 0.), (1., 1.), (0., 1.)))
    height = .02 + .04 * xy[:, 0]
    rules = tuple(LinearRule(a, b) for a, b in zip(height, np.roll(height, -1)))
    return BoundaryData(xy, height, np.tile((.04, 0.), (4, 1)), rules, (), 0., .1, False)


class PythonCompoundSurfaceTests(unittest.TestCase):
    def test_field_boundary_values_tangent_planes_and_batched_queries(self):
        boundary = boundary_data()
        field = BoundaryJetField((boundary,))
        points = np.column_stack((np.linspace(.1, .9, 21), np.zeros(21)))
        np.testing.assert_allclose(field(points), .02 + .04 * points[:, 0], atol=1e-12)
        shifted = points + (0., 1e-5)
        np.testing.assert_allclose((field(shifted) - field(points)) / 1e-5, 0., atol=1e-6)
        shifted = points + (1e-5, 0.)
        np.testing.assert_allclose((field(shifted) - field(points)) / 1e-5, .04, atol=1e-6)
        random = np.random.default_rng(16).uniform(0., 1., (1001, 2))
        small = BoundaryJetField((boundary,), scratch_bytes=4096)
        np.testing.assert_allclose(field(random), small(random), atol=1e-15, rtol=0.)
        self.assertTrue(np.all((field(random) >= 0.) & (field(random) <= .1)))
        self.assertLess(small.batch_rows, field.batch_rows)

    def test_constant_field_and_invalid_queries_are_explicit(self):
        b = boundary_data()
        constant = replace(b, rules=(PlateauEdge(.04),) * 4)
        field = BoundaryJetField((constant,))
        np.testing.assert_allclose(field(np.array(((0., 0.), (.5, .5), (1., 1.)))), .04, atol=1e-14)
        self.assertEqual(field(np.empty((0, 2))).shape, (0,))
        for xy in (np.ones((2, 3)), np.array(((float('nan'), 0.),))):
            with self.assertRaises(ValueError):
                field(xy)
        with self.assertRaisesRegex(ValueError, 'one elevation envelope'):
            BoundaryJetField((b, replace(b, ceiling=.2)))
        with self.assertRaises(ValueError):
            BoundaryJetField(())

    def test_local_refinement_keeps_boundary_and_has_no_t_junctions(self):
        patterns = PatternLibrary()
        grid = patterns.grid(3, 3)
        mesh = Mesh(np.column_stack((grid.uv, np.full(len(grid.uv), .02))), grid.quads)
        before_edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(before_edges, axis=0, return_counts=True)
        rim = np.unique(unique[counts == 1])
        mesh.seams = {0: rim}
        refined, remap = QuadRefiner(patterns).refine(mesh, np.array((4,)), lambda xy: np.full(len(xy), .02))
        np.testing.assert_array_equal(refined.xyz[remap], mesh.xyz)
        np.testing.assert_array_equal(refined.seams[0], remap[rim])
        edges = np.sort(np.stack((refined.quads, np.roll(refined.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        self.assertEqual(len(refined.xyz) - len(unique) + len(refined.quads), 1)
        self.assertEqual(int(np.sum(counts == 1)), 12)
        self.assertTrue(set(counts).issubset({1, 2}))
        self.assertLess(refined.validate()['maximum_edge_ratio'], 8.)
        validate_overlaps({model.Key('refined', 1): refined})
        with self.assertRaisesRegex(ValueError, 'cell budget'):
            QuadRefiner(patterns, max_cells=4).refine(mesh, np.array((4,)), lambda xy: np.zeros(len(xy)))

    def test_compound_surfaces_preserve_full_profiles_holes_and_elevation(self):
        for network in (compound_cross_network(model), short_junction_pair(model)):
            with self.subTest(network=network.uuid):
                group, context = inputs(network)
                result = CompoundSurfaceBuilder(context.patterns).build(group, context)
                self.assertLess(result.metadata['total_quads'], 5000)
                self.assertLessEqual(result.metadata['quality']['maximum_edge_ratio'], 8.)
                self.assertLessEqual(result.metadata['physical_slopes']['maximum'], result.metadata['measured_slope_limit'])
                self.assertGreaterEqual(float(result.mesh.xyz[:, 2].min()), 0.)
                self.assertLessEqual(float(result.mesh.xyz[:, 2].max()), .08 * network.settings.exaggeration)
                actual_group = result.chart.partition.group
                layouts = {p.region.owner: p.layout for p in actual_group.footprints}
                from tested_python_streets.conditions import BoundaryConditions
                data = [BoundaryConditions().build(actual_group, loop, layouts) for loop in actual_group.loops]
                for vertex, ref in result.boundary_refs.items():
                    xy, height, gradient = data[ref.loop].edge_samples(ref.edge, ref.fraction)
                    np.testing.assert_allclose(result.mesh.xyz[vertex, :2], xy, atol=1e-8, rtol=0.)
                    self.assertAlmostEqual(result.mesh.xyz[vertex, 2], height * network.settings.exaggeration)
                meshes = {model.Key('compound', 1): result.mesh}
                meshes.update({model.Key('hole', key): mesh for key, mesh in result.holes.items()})
                validate_overlaps(meshes)

    def test_exaggeration_changes_only_output_elevation(self):
        network = short_junction_pair(model)
        first_group, first_context = inputs(network)
        first = CompoundSurfaceBuilder(first_context.patterns).build(first_group, first_context)
        other = replace(network, settings=replace(network.settings, exaggeration=1.))
        group, context = inputs(other)
        second = CompoundSurfaceBuilder(context.patterns).build(group, context)
        np.testing.assert_array_equal(first.mesh.quads, second.mesh.quads)
        np.testing.assert_array_equal(first.mesh.xyz[:, :2], second.mesh.xyz[:, :2])
        np.testing.assert_allclose(first.mesh.xyz[:, 2], 10. * second.mesh.xyz[:, 2], atol=1e-14)
        self.assertAlmostEqual(first.metadata['physical_slopes']['maximum'], second.metadata['physical_slopes']['maximum'])
        raised = replace(network, sections=tuple(replace(s, datum=12.5) for s in network.sections))
        group, context = inputs(raised)
        shifted = CompoundSurfaceBuilder(context.patterns).build(group, context)
        np.testing.assert_array_equal(first.mesh.quads, shifted.mesh.quads)
        np.testing.assert_array_equal(first.mesh.xyz[:, :2], shifted.mesh.xyz[:, :2])
        np.testing.assert_allclose(shifted.mesh.xyz[:, 2], first.mesh.xyz[:, 2] + 12.5, atol=1e-10, rtol=0.)


if __name__ == '__main__':
    unittest.main()
