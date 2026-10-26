from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import compound_cross_network, short_junction_pair
from test_python_street_compounds import footprint_context
from tested_python_streets import model
from tested_python_streets.compounds import CompoundPlanner
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.partition import BoundaryRef, CompoundPartitioner, CompoundQuadChart, CutPolicy, JunctionCut
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.validation import validate_overlaps


def inputs(network):
    context = footprint_context(network)
    layouts = {site.key: IntersectionProvider.footprint(context, site) for site in network.intersections}
    return CompoundPlanner().build(layouts)[0], context


class PythonCompoundPartitionTests(unittest.TestCase):
    def test_divided_cross_has_four_source_owned_blocks_and_keeps_its_hole(self):
        network = compound_cross_network(model)
        group, context = inputs(network)
        result = CompoundPartitioner().build(group, context)
        self.assertEqual(len(result.blocks), 4)
        self.assertEqual(len(result.cuts), 4)
        self.assertEqual(len(result.boundary_loops), 2)
        self.assertLess(result.area_error, 1e-8)
        self.assertEqual({owner for block in result.blocks for owner in block.owners}, set(group.members))
        self.assertTrue(all(len(block.owners) == 1 for block in result.blocks))
        self.assertTrue(all(np.all(np.isfinite(block.center)) for block in result.blocks))
        for loop in group.loops:
            for point in loop.xy:
                self.assertTrue(np.any(np.all(result.xy == point, axis=1)))

    def test_shared_cut_samples_match_in_reverse_and_resolve_every_block_parity(self):
        group, context = inputs(compound_cross_network(model))
        result = CompoundPartitioner().build(group, context)
        sampled = result.sample(.7)
        shared = {}
        for block, rim in zip(result.blocks, sampled.rims):
            self.assertEqual(len(rim) % 2, 0)
            at = 0
            for edge_index, forward in block.edges:
                edge = result.edges[edge_index]
                count = sampled.cut_intervals.get(edge_index, 1)
                sequence = tuple(rim[(at + j) % len(rim)] for j in range(count + 1))
                self.assertEqual(sequence[0], edge.first if forward else edge.last)
                self.assertEqual(sequence[-1], edge.last if forward else edge.first)
                if edge.cut is not None:
                    shared.setdefault(edge_index, []).append(sequence)
                at += count
            self.assertEqual(at, len(rim))
        for sequences in shared.values():
            self.assertEqual(len(sequences), 2)
            self.assertEqual(sequences[0], sequences[1][::-1])
        self.assertTrue(all(len(loop) % 2 == 0 for loop in sampled.boundary_loops))
        sources = {p.region.owner: p for p in group.footprints}
        for node, ref in sampled.boundary_refs.items():
            if 1e-10 < ref.fraction < 1. - 1e-10:
                span = group.loops[ref.loop].spans[ref.edge]
                self.assertEqual(sources[span.owner].segments[span.segment].kind, 'curb')

    def test_short_two_site_link_produces_two_compatible_blocks(self):
        group, context = inputs(short_junction_pair(model))
        partition = CompoundPartitioner().build(group, context)
        self.assertEqual(len(partition.blocks), 2)
        self.assertEqual(len(partition.cuts), 1)
        self.assertEqual({b.owners[0] for b in partition.blocks}, set(group.members))
        sampled = partition.sample(.7)
        self.assertTrue(all(len(rim) % 2 == 0 for rim in sampled.rims))

    def test_small_source_bend_retains_semantic_cut_ownership(self):
        network = compound_cross_network(model)
        section = network.sections[1]
        a, b = section.points
        moved = replace(section, points=(a, (.5 * (a[0] + b[0]), .12 + .5 * (a[1] + b[1]), 0.), b))
        changed = replace(network, sections=(network.sections[0], moved, *network.sections[2:]))
        before, context_before = inputs(network)
        after, context_after = inputs(changed)
        first = CompoundPartitioner().build(before, context_before)
        second = CompoundPartitioner().build(after, context_after)
        self.assertEqual([cut.section_id for cut in first.cuts], [cut.section_id for cut in second.cuts])
        self.assertEqual({b.owners for b in first.blocks}, {b.owners for b in second.blocks})
        self.assertEqual([s.id for s in network.sections], [s.id for s in changed.sections])

    def test_cuts_cannot_fill_the_separator_hole_and_budgets_are_explicit(self):
        group, context = inputs(compound_cross_network(model))
        class CrossHole(CutPolicy):
            def cuts(self, group, context):
                points = group.loops[0].xy
                first = min(np.flatnonzero(points[:, 0] < 0.), key=lambda i: abs(points[i, 1]))
                last = min(np.flatnonzero(points[:, 0] > 0.), key=lambda i: abs(points[i, 1]))
                claim = group.claims[0]
                return (JunctionCut(claim.section_id, claim.first_owner, claim.last_owner,
                                    BoundaryRef(0, int(first), 0.), BoundaryRef(0, int(last), 0.)),)
        with self.assertRaisesRegex(ValueError, 'pavement'):
            CompoundPartitioner(policy=CrossHole()).build(group, context)
        with self.assertRaisesRegex(ValueError, 'cut budget'):
            CompoundPartitioner(max_cuts=1).build(group, context)
        with self.assertRaisesRegex(ValueError, 'graph-node budget'):
            CompoundPartitioner(max_nodes=4).build(group, context)

    def test_compound_pavement_and_separator_are_conforming_positive_quad_regions(self):
        network = compound_cross_network(model)
        group, context = inputs(network)
        result = CompoundQuadChart(PatternLibrary()).build(group, context)
        self.assertEqual(len(result.holes), 1)
        self.assertLess(result.metadata['total_quads'], 10000)
        self.assertLessEqual(result.metadata['quality']['maximum_edge_ratio'], 8.)
        meshes = {model.Key('compound', 1): result.mesh}
        for li, hole in result.holes.items():
            self.assertEqual(hole.mesh.material, 'grass')
            meshes[model.Key('grass', li)] = hole.mesh
            np.testing.assert_array_equal(result.mesh.xyz[result.mesh.seams[li], :2][::-1],
                                          hole.mesh.xyz[hole.mesh.seams[0], :2])
        validate_overlaps(meshes)
        for node, ref in result.boundary_refs.items():
            loop = result.partition.group.loops[ref.loop].xy
            expected = (1. - ref.fraction) * loop[ref.edge] + ref.fraction * loop[(ref.edge + 1) % len(loop)]
            np.testing.assert_allclose(result.mesh.xyz[node, :2], expected, atol=1e-12, rtol=0.)

    def test_two_site_compound_needs_no_separator_object(self):
        group, context = inputs(short_junction_pair(model))
        result = CompoundQuadChart(PatternLibrary()).build(group, context)
        self.assertEqual(result.holes, {})
        self.assertEqual(len(result.metadata['blocks']), 2)
        validate_overlaps({model.Key('compound', 1): result.mesh})


if __name__ == '__main__':
    unittest.main()
