"""Certified physical rims with holes retain their boundaries and quad guards."""
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch
import math
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import compound_cross_network
from test_python_street_partition import inputs
from tested_python_streets import model
from tested_python_streets.compounds import CompoundGroup, RimSegment, RimSpan, rim_label
from tested_python_streets.mesh import publication_coordinates
from tested_python_streets.partition import CompoundPartitioner, CompoundQuadChart
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.physical_hole_cuts import PhysicalHoleCuts, PhysicalHoleCutError
from tested_python_streets.regions import BoundaryLoop
from tested_python_streets.sampling import AnalyticalBoundarySampler
from tested_python_streets.validation import validate_overlaps
from test_python_street_source_cycles import triangle, prepare


def rectangle(owner, bounds, clockwise=False):
    x0, y0, x1, y1 = bounds
    corners = np.array(((x0, y0), (x1, y0), (x1, y1), (x0, y1)))
    if clockwise:
        corners = corners[::-1]
    points, labels = [], []
    for i, (a, b) in enumerate(zip(corners, np.roll(corners, -1, axis=0))):
        count = max(2, int(math.ceil(np.linalg.norm(b-a)/.5)))
        points.extend(np.linspace(a, b, count+1)[:-1])
        labels.extend(RimSpan(owner, RimSegment('curb', i, j/count, (j+1)/count)) for j in range(count))
    return BoundaryLoop(np.array(points), tuple(labels))


def holed_region(holes, angle=0., offset=(0., 0.)):
    owner = model.Key('intersection', 999)
    loops = (rectangle(owner, (-6., -5., 6., 5.)),
             *tuple(rectangle(owner, bounds, True) for bounds in holes))
    c, s = math.cos(angle), math.sin(angle)
    rotation = np.array(((c, -s), (s, c)))
    loops = tuple(BoundaryLoop(loop.xy @ rotation.T + offset, loop.spans) for loop in loops)
    return CompoundGroup((owner,), (), (), loops, object())


class PhysicalHolePartitionTests(unittest.TestCase):
    def test_one_and_two_holes_keep_exact_rims_area_and_conforming_positive_quads(self):
        context = SimpleNamespace(network=SimpleNamespace(settings=SimpleNamespace(curve_tolerance=.005)))
        cases = (((-2., -1., 2., 1.),), ((-2., -.078, 2., .078),),
                 ((-4., -1., -2., 1.), (2., -1., 4., 1.)),
                 ((-4., -1., -3., 1.), (-.5, -.078, .5, .078), (3., -1., 4., 1.)))
        for holes in cases:
            for angle, offset in ((0., (0., 0.)), (.713, (1200., -1800.))):
                with self.subTest(holes=len(holes), angle=angle):
                    group = holed_region(holes, angle, offset)
                    with patch.object(AnalyticalBoundarySampler, 'build', lambda self, value: value):
                        result = CompoundQuadChart(PatternLibrary()).build(group, context)
                    self.assertEqual(len(result.holes), len(holes))
                    self.assertTrue(result.metadata['physical_hole_partition'])
                    self.assertLess(result.partition.area_error, 1e-7)
                    published = result.mesh.with_positions(publication_coordinates(result.mesh.xyz, 1., physical=True))
                    self.assertLessEqual(published.validate()['maximum_edge_ratio'], 8.)
                    meshes = {model.Key('compound', 999): result.mesh}
                    for li, hole in result.holes.items():
                        self.assertLessEqual(hole.mesh.validate()['maximum_edge_ratio'], 8.)
                        meshes[model.Key('grass', li)] = hole.mesh
                        np.testing.assert_array_equal(result.mesh.xyz[result.mesh.seams[li]][::-1],
                                                      hole.mesh.xyz[hole.mesh.seams[0]])
                    validate_overlaps(meshes)
                    for li, loop in enumerate(group.loops):
                        boundary = result.mesh.xyz[result.mesh.seams[li], :2]
                        for point in loop.xy:
                            self.assertTrue(np.any(np.all(boundary == point, axis=1)))
                    for node, ref in result.boundary_refs.items():
                        loop = group.loops[ref.loop]
                        expected = (1.-ref.fraction)*loop.xy[ref.edge] + ref.fraction*loop.xy[(ref.edge+1) % len(loop.xy)]
                        np.testing.assert_allclose(result.mesh.xyz[node, :2], expected, atol=5e-10, rtol=0.)

    def test_source_curves_and_prescribed_profile_rows_survive_physical_hole_dispatch(self):
        group, context = inputs(compound_cross_network(model))
        group = replace(group, physical_layout=object())
        result = CompoundQuadChart(PatternLibrary()).build(group, context)
        self.assertEqual(len(result.holes), 1)
        self.assertLessEqual(result.metadata['quality']['maximum_edge_ratio'], 8.)
        sources = {footprint.region.owner: footprint for footprint in group.footprints}
        for ref in result.boundary_refs.values():
            if 1e-10 < ref.fraction < 1.-1e-10:
                loop = result.partition.group.loops[ref.loop]
                self.assertEqual(rim_label(loop.spans[ref.edge], sources).kind, 'curb')
        # Existing source-ownership cuts remain authoritative without a
        # certified replacement layout, including their original aliases.
        ordinary, ordinary_context = inputs(compound_cross_network(model))
        partition = CompoundPartitioner().build(ordinary, ordinary_context)
        self.assertEqual(len(partition.blocks), 4)
        self.assertTrue(all(cut.source_sections for cut in partition.cuts))

    def test_hole_cut_and_node_budgets_remain_explicit(self):
        group = holed_region(((-4., -1., -2., 1.), (2., -1., 4., 1.)))
        with self.assertRaisesRegex(ValueError, 'hole budget'):
            PhysicalHoleCuts(maximum_holes=1).cuts(group, None)
        with self.assertRaisesRegex(ValueError, 'cut budget'):
            PhysicalHoleCuts(maximum_cuts=2).cuts(group, None)
        with self.assertRaisesRegex(ValueError, 'graph-node budget'):
            CompoundPartitioner(max_nodes=4).build(group, None)
        with self.assertRaisesRegex(ValueError, 'certified exterior'):
            PhysicalHoleCuts().cuts(replace(group, physical_layout=None), None)
        with self.assertRaisesRegex(ValueError, 'require original source corridors'):
            PhysicalHoleCuts(hole_bridges=True)

    def test_finite_cut_exhaustion_is_distinct_from_invalid_physical_inputs(self):
        group = holed_region(((-2., -1., 2., 1.),))
        policy = PhysicalHoleCuts()
        with patch.object(policy, '_cuts', side_effect=ValueError('uncut hole')) as cuts:
            with self.assertRaises(PhysicalHoleCutError):
                policy.cuts(group, None)
        self.assertEqual(cuts.call_count, 4)
        with self.assertRaises(ValueError) as caught:
            policy.cuts(replace(group, physical_layout=None), None)
        self.assertNotIsInstance(caught.exception, PhysicalHoleCutError)

    def test_source_corridors_with_two_hole_connections_keep_exact_physical_rims(self):
        context, original = prepare(triangle(2.))
        group = AnalyticalBoundarySampler(spacing=context.network.settings.spacing,
            tolerance=context.network.settings.curve_tolerance).build(original)
        partition = CompoundPartitioner(policy=PhysicalHoleCuts(source_corridors=True,
            hole_bridges=True), require_visibility=False).build(group, context)
        self.assertLess(partition.area_error, 1e-7)
        self.assertGreaterEqual(len(partition.cuts), 4)
        self.assertEqual(sum(loop.hole for loop in partition.group.loops), 1)
        sources = {source.region.owner: source for source in group.footprints}
        for li, loop in enumerate(group.loops):
            boundary = partition.xy[list(partition.boundary_loops[li])]
            for point in loop.xy:
                self.assertTrue(np.any(np.all(boundary == point, axis=1)))
        for ref in partition.boundary_refs.values():
            if 1e-10 < ref.fraction < 1.-1e-10:
                self.assertEqual(rim_label(group.loops[ref.loop].spans[ref.edge], sources).kind, 'curb')

    def test_computational_cut_endpoints_can_reuse_existing_physical_vertices(self):
        group = holed_region(((-2., -1., 2., 1.),))
        cuts = PhysicalHoleCuts(existing_vertices=True).cuts(group, None)
        self.assertTrue(cuts)
        for cut in cuts:
            for ref in (cut.first, cut.last):
                self.assertEqual(ref.fraction, 0.)
                self.assertTrue(0 <= ref.edge < len(group.loops[ref.loop].xy))


if __name__ == '__main__':
    unittest.main()
