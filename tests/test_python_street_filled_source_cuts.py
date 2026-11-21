"""A filled source cycle can use original corridor cuts without inventing a hole."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_source_cycles import triangle, prepare
from test_python_street_physical_visibility_cuts import contained_owner
from tested_python_streets.compounds import rim_label
from tested_python_streets.four_side_charts import quality
from tested_python_streets.mesh import publication_coordinates
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.physical_hole_cuts import PhysicalHoleCuts
from tested_python_streets.sampling import AnalyticalBoundarySampler


class FilledSourceCutsTests(unittest.TestCase):
    def test_certified_filled_cycle_keeps_exact_source_rim_and_prescribed_profiles(self):
        context, original = prepare(triangle(8.))
        sampled = AnalyticalBoundarySampler(spacing=context.network.settings.spacing,
            tolerance=context.network.settings.curve_tolerance).build(original)
        self.assertFalse(PhysicalHoleCuts.eligible(original))
        self.assertTrue(PhysicalHoleCuts.filled_eligible(original))
        result = CompoundQuadChart(context.patterns)._filled_source(
            original, context, sampled, 'unqualified ordinary seed')
        self.assertEqual(len(result.partition.group.loops), 1)
        self.assertFalse(result.holes)
        self.assertEqual(result.partition.group.claims, original.claims)
        self.assertTrue(result.metadata['filled_source_partition'])
        self.assertGreater(result.metadata['source_width_visibility']['total_cells'], 0)
        self.assertLess(result.partition.area_error, 1e-7)
        boundary = result.mesh.xyz[result.mesh.seams[0], :2]
        for point in sampled.loops[0].xy:
            self.assertTrue(np.any(np.all(boundary == point, axis=1)))
        sources = {source.region.owner: source for source in result.partition.group.footprints}
        for node, ref in result.boundary_refs.items():
            loop = result.partition.group.loops[ref.loop]
            expected = (1.-ref.fraction)*loop.xy[ref.edge]+ref.fraction*loop.xy[(ref.edge+1) % len(loop.xy)]
            np.testing.assert_allclose(result.mesh.xyz[node, :2], expected, atol=5e-10, rtol=0.)
            if 1e-10 < ref.fraction < 1.-1e-10:
                self.assertEqual(rim_label(loop.spans[ref.edge], sources).kind, 'curb')
        stored = quality(result.mesh.with_positions(publication_coordinates(result.mesh.xyz, 1., physical=True)))
        self.assertLessEqual(stored['maximum_edge_ratio'], 8.)
        self.assertGreaterEqual(stored['minimum_corner_sine'], .03)

    def test_filled_cut_policy_requires_explicit_source_union_qualification(self):
        group, context = contained_owner()
        ordinary = replace(group, physical_layout=SimpleNamespace())
        self.assertFalse(PhysicalHoleCuts.filled_eligible(ordinary))
        with self.assertRaisesRegex(ValueError, 'certified exterior'):
            PhysicalHoleCuts(allow_filled=True).cuts(ordinary, context)
        incomplete = replace(group, physical_layout=SimpleNamespace(source_region_union=True,
            internal_profile_fragments={'complete_source_width': {}}))
        self.assertFalse(PhysicalHoleCuts.filled_eligible(incomplete))
        certified = replace(group, physical_layout=SimpleNamespace(source_cycle_proofs={}))
        with self.assertRaisesRegex(ValueError, 'certified exterior'):
            PhysicalHoleCuts().cuts(certified, context)

    def test_filled_source_retry_still_rejects_an_incomplete_original_strip(self):
        group, context = contained_owner(outside=True)
        group = replace(group, physical_layout=SimpleNamespace(source_cycle_proofs={}))
        with self.assertRaisesRegex(ValueError, 'complete source strip'):
            CompoundQuadChart(context.patterns)._filled_source(group, context, group, 'unqualified seed')


if __name__ == '__main__':
    unittest.main()
