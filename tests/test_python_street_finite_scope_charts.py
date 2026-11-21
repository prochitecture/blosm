"""The finite scope collar atlas retains its exact boundary-first contract."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.finite_scope_charts import fixed_rim_chart
from tested_python_streets.mesh import Mesh
from tested_python_streets.patterns import PatternLibrary


def boundary():
    corners=np.array(((0.,0.),(7.25,0.),(5.8,7.),(0.,7.)))
    counts=(38,12,26,8)
    return np.vstack([np.linspace(a,b,n+1)[:-1] for a,b,n in
        zip(corners,np.roll(corners,-1,axis=0),counts)])


class FiniteScopeChartTests(unittest.TestCase):
    def test_graded_atlas_keeps_exact_boundary_first_and_stored_quality(self):
        for shift in ((0.,0.),(-670.,-635.)):
            rim=boundary()+shift
            result=fixed_rim_chart(rim,PatternLibrary())
            np.testing.assert_array_equal(result.mesh.xyz[:len(rim),:2],rim)
            np.testing.assert_array_equal(result.mesh.seams[0],np.arange(len(rim)))
            self.assertTrue(result.metadata['local_collars'])
            self.assertTrue(result.metadata['all_direct_quads'])
            self.assertLessEqual(result.mesh.validate()['maximum_edge_ratio'],8.)
            stored=Mesh(result.mesh.xyz.astype(np.float32).astype(float),result.mesh.quads)
            self.assertLessEqual(stored.validate()['maximum_edge_ratio'],8.)

    def test_invalid_rim_and_allocation_are_rejected(self):
        with self.assertRaises(ValueError):fixed_rim_chart(boundary()[:-1],PatternLibrary())
        with self.assertRaises(ValueError):fixed_rim_chart(boundary(),PatternLibrary(),max_candidates=7)
        with self.assertRaises(ValueError):fixed_rim_chart(boundary(),PatternLibrary(),max_cells=4)


if __name__=='__main__':unittest.main()
