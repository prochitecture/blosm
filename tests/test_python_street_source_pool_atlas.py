"""Only a fresh, explicitly measured new input class changes atlas priority."""
from copy import copy,deepcopy
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.source_pool_atlas import admission,build
from tested_python_streets.source_pool_budget import certify
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.model import Settings
from tested_python_streets.patterns import PatternLibrary
from test_python_street_source_pool_budget import ring
from test_python_street_source_cycles import prepare
from test_python_street_source_rim_normalization import fixture


class SourcePoolAtlasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context,cls.group=prepare(ring())
        certify(cls.context,cls.group)

    def test_matching_fresh_budget_and_width_proof_preserve_the_accepted_plan(self):
        context,group=self.context,self.group
        budget=group.physical_layout.extended_source_pool
        proof=group.physical_layout.source_cycle_full_width_proof
        snapshot=context.network.to_dict()
        result=admission(group,context)
        self.assertEqual(result,budget)
        self.assertIs(group.physical_layout.extended_source_pool,budget)
        self.assertIs(group.physical_layout.source_cycle_full_width_proof,proof)
        self.assertEqual(context.network.to_dict(),snapshot)

    def test_stale_work_counts_or_downgraded_width_proof_cannot_change_priority(self):
        context,group=self.context,self.group
        layout=copy(group.physical_layout);layout.extended_source_pool=dict(layout.extended_source_pool,source_segments=1)
        with self.assertRaisesRegex(ValueError,'matching fresh measured-work'):
            admission(replace(group,physical_layout=layout),context)
        layout=copy(group.physical_layout);layout.source_cycle_full_width_proof=deepcopy(layout.source_cycle_full_width_proof)
        layout.source_cycle_full_width_proof['total_cells']=1
        with self.assertRaisesRegex(ValueError,'matching fresh complete-width'):
            admission(replace(group,physical_layout=layout),context)

    def test_missing_new_admission_raw_ordinary_and_eight_sites_keep_old_routes(self):
        context,group=self.context,self.group
        layout=copy(group.physical_layout);del layout.extended_source_pool
        self.assertIsNone(admission(replace(group,physical_layout=layout),context))
        self.assertIsNone(admission(replace(group,physical_layout=None),context))
        self.assertIsNone(admission(replace(group,members=group.members[:8]),context))
        layout=copy(group.physical_layout);layout.source_cycle_full_width_proof=None
        self.assertIsNone(admission(replace(group,physical_layout=layout),context))
        # A real already-positive ordinary exterior follows its unchanged chart.
        old=fixture();ctx=SimpleNamespace(network=SimpleNamespace(settings=Settings()),patterns=PatternLibrary())
        chart=CompoundQuadChart(ctx.patterns)
        expected=chart._build(old,ctx,require_quality=True)
        with patch('tested_python_streets.source_pool_atlas.build',side_effect=AssertionError('old route changed')):
            actual=chart.build(old,ctx)
        np.testing.assert_array_equal(actual.mesh.xyz,expected.mesh.xyz)
        np.testing.assert_array_equal(actual.mesh.quads,expected.mesh.quads)

    def test_priority_uses_same_original_corridors_and_retains_the_full_admission_metadata(self):
        context,group=self.context,self.group
        certificate=admission(group,context)
        sentinel=SimpleNamespace(metadata={'original_source_cycle_atlas':{'maximum_corridor_attempts':3}})
        chart=CompoundQuadChart(context.patterns)
        with patch('tested_python_streets.source_cycle_atlas.build',return_value=sentinel) as atlas:
            result=chart.build(group,context)
        self.assertIs(result,sentinel)
        self.assertIs(atlas.call_args.args[1],group)
        self.assertEqual(result.metadata['measured_source_pool_priority']['admission'],certificate)
        self.assertTrue(result.metadata['measured_source_pool_priority']['original_corridor_order_unchanged'])
        self.assertEqual(result.metadata['measured_source_pool_priority']['cell_budget'],50000)


if __name__=='__main__':unittest.main()
