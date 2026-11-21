"""The last source-normal seed family keeps fixed source identity and limits."""
from unittest.mock import patch
import unittest
import numpy as np
import _python_street_package
from test_python_street_interior_quad_refinement import distorted_grid
from tested_python_streets.mesh import Mesh
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.source_normal_energy import EnergySourceNormalChart,_energy_attempt
from tested_python_streets.validation import validate_overlaps
from tested_python_streets.model import Key


class SourceNormalEnergyTests(unittest.TestCase):
    def test_positive_bad_candidate_repairs_without_moving_source_rim(self):
        mesh=distorted_grid();before=mesh.xyz.copy()
        chart=EnergySourceNormalChart(PatternLibrary(),max_cells=1000)
        result,correction,energy,physical,stored=_energy_attempt(chart,mesh,mesh,1,True)
        self.assertLessEqual(stored['maximum_edge_ratio'],8.)
        self.assertGreaterEqual(stored['minimum_corner_sine'],.03)
        self.assertGreater(energy['steps'],0)
        np.testing.assert_array_equal(result.seams[0],mesh.seams[0])
        np.testing.assert_array_equal(result.xyz[mesh.seams[0]],before[mesh.seams[0]])
        np.testing.assert_array_equal(mesh.xyz,before)
        validate_overlaps({Key('pavement',0):result})

    def test_an_old_qualified_chart_never_uses_the_extra_energy_family(self):
        mesh=Mesh(np.array(((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.))),np.array(((0,1,2,3),)))
        chart=EnergySourceNormalChart(PatternLibrary())
        with patch('tested_python_streets.source_normal_energy._energy_attempt_steps') as recovery:
            result,details=chart._correct(mesh)
        self.assertIs(result,mesh)
        self.assertEqual(chart.energy_recoveries,[])
        recovery.assert_not_called()

    def test_wrong_source_rim_nonplanar_or_over_budget_candidate_cannot_qualify(self):
        mesh=distorted_grid();chart=EnergySourceNormalChart(PatternLibrary(),max_cells=1000)
        xyz=mesh.xyz.copy();xyz[mesh.seams[0][0],0]-=.1
        with self.assertRaisesRegex(ValueError,'physical rim'):
            _energy_attempt(chart,mesh,mesh.with_positions(xyz),1,True)
        xyz=mesh.xyz.copy();xyz[0,2]=.1
        with self.assertRaisesRegex(ValueError,'bounded planar'):
            _energy_attempt(chart,mesh,mesh.with_positions(xyz),1,True)
        with self.assertRaisesRegex(ValueError,'bounded planar'):
            _energy_attempt(EnergySourceNormalChart(PatternLibrary(),max_cells=4),mesh,mesh,1,True)


if __name__=='__main__':unittest.main()
