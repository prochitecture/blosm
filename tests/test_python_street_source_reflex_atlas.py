"""Source normals and true outer-curb bisectors are computational only."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.source_reflex_atlas import SourceReflexCuts,eligible,build
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.compounds import rim_label
from tested_python_streets.four_side_charts import quality
from tested_python_streets.mesh import publication_coordinates
from test_python_street_source_cycles import triangle,prepare


class SourceReflexAtlasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        network=triangle(8.)
        network=replace(network,sections=tuple(replace(s,points=(s.points[0],(0.,-30.,0.)))
            if s.id==4 else s for s in network.sections))
        cls.context,cls.group=prepare(network)

    def test_combined_cuts_use_original_source_normals_and_actual_outer_curbs(self):
        context,group=self.context,self.group;before=context.network.to_dict()
        sources={f.region.owner:f for f in group.footprints}
        labels=[[rim_label(span,sources) for span in loop.spans] for loop in group.loops]
        policy=SourceReflexCuts(source_corridors=True,allow_filled=True)
        cuts=policy._source_cuts(group,context,labels)
        self.assertEqual(len(cuts),policy.source_normal_count+policy.outer_reflex_count)
        self.assertGreater(policy.source_normal_count,0);self.assertGreater(policy.outer_reflex_count,0)
        self.assertLessEqual(len(cuts),128)
        self.assertTrue(all(c.first.loop==0 and c.last.loop==0 for c in cuts))
        self.assertEqual(context.network.to_dict(),before)

    def test_missing_source_proof_and_real_holes_do_not_enter_the_filled_family(self):
        context,group=self.context,self.group
        self.assertTrue(eligible(group))
        self.assertFalse(eligible(replace(group,physical_layout=None)))
        self.assertFalse(eligible(replace(group,claims=())))
        _,hole=prepare(triangle(2.))
        self.assertFalse(eligible(hole))
        with self.assertRaisesRegex(ValueError,'source-certified'):
            build(CompoundQuadChart(context.patterns),hole,context,[])

    def test_small_complete_source_pool_keeps_rim_and_stored_quad_contract(self):
        context,group=self.context,self.group
        result=build(CompoundQuadChart(context.patterns),group,context,[])
        self.assertEqual(result.partition.group.members,group.members)
        self.assertEqual(result.partition.group.claims,group.claims)
        report=result.metadata['source_normal_outer_reflex_recovery']
        self.assertFalse(report['physical_boundary_changed'])
        self.assertEqual(report['final_width_proof']['maximum_error'],.002)
        stored=quality(result.mesh.with_positions(publication_coordinates(result.mesh.xyz,1.,physical=True)))
        self.assertLessEqual(stored['maximum_edge_ratio'],8.)
        self.assertGreaterEqual(stored['minimum_corner_sine'],.03)


if __name__=='__main__':unittest.main()
