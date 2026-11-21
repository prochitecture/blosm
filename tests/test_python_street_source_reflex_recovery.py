"""A late compound atlas must prove complete original consumed widths."""
from dataclasses import replace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from test_python_street_deferred_compounds import fixture
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.compound_chart_recovery import recover_source_reflex
from tested_python_streets.four_side_charts import quality
from tested_python_streets.mesh import publication_coordinates
from tested_python_streets.ownership import PLAN
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.updates import BuildContext


class SourceReflexRecoveryTests(unittest.TestCase):
    def test_complete_compound_keeps_source_widths_and_existing_successes_first(self):
        network=fixture();engine=PythonMeshBackend(network.settings)
        with patch('tested_python_streets.compound_chart_recovery.recover_source_reflex') as late:
            prepared=engine.prepare(network)
        late.assert_not_called()
        group=prepared.products[PLAN].groups[0]
        context=BuildContext(network,prepared.products,engine.patterns,
            prepared.task_fingerprints,prepared.source_fingerprints)
        original=network.to_dict()
        chart=recover_source_reflex(CompoundQuadChart(engine.patterns,spacing=network.settings.spacing),
            group,context,[dict(reason='earlier atlas families exhausted')])
        report=chart.metadata['fixed_rim_compound_recovery']
        self.assertTrue(report['original_width_proof'])
        self.assertTrue(report['final_width_proof'])
        self.assertEqual(network.to_dict(),original)
        for mesh in (chart.mesh,chart.mesh.with_positions(publication_coordinates(chart.mesh.xyz,1.,physical=True))):
            measured=quality(mesh)
            self.assertLessEqual(measured['maximum_edge_ratio'],8.)
            self.assertGreaterEqual(measured['minimum_corner_sine'],.03)
        # A geometrically clipped consumed source cannot borrow the chart's
        # good grading as permission to omit part of its carriageway width.
        loop=group.loops[0];center=np.mean(loop.xy,axis=0)
        clipped=replace(group,loops=(replace(loop,xy=center+.1*(loop.xy-center)),))
        with self.assertRaises(ValueError):
            recover_source_reflex(CompoundQuadChart(engine.patterns),clipped,context,[])


if __name__=='__main__':unittest.main()
