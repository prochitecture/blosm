"""Conforming uniform quad refinement and physical boundary constraints."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.quad_subdivision import UniformQuadSubdivision


class PythonUniformQuadSubdivisionTests(unittest.TestCase):
    def mesh(self):
        return Mesh(np.array(((0.,0.,5.), (2.,0.,5.), (2.,1.,5.), (0.,1.,5.),
                              (4.,0.,5.), (4.,1.,5.))), np.array(((0,1,2,3), (1,4,5,2))))

    def field(self, xy):
        return 5.+.2*np.sin(np.pi*xy[:,0]/4.)*np.sin(np.pi*xy[:,1])

    def test_reversed_shared_edges_are_conforming_at_each_factor(self):
        patterns = PatternLibrary()
        for factor in (2, 3, 4):
            with self.subTest(factor=factor):
                result = UniformQuadSubdivision(patterns).build(self.mesh(), factor, self.field, datum=5., ceiling=5.2)
                self.assertEqual(len(result.quads), 2*factor*factor)
                edges = np.sort(np.stack((result.quads,np.roll(result.quads,-1,axis=1)),axis=-1).reshape(-1,2),axis=1)
                unique, counts = np.unique(edges, axis=0, return_counts=True)
                self.assertEqual(len(result.xyz)-len(unique)+len(result.quads),1)
                self.assertLessEqual(int(counts.max()),2)
                self.assertEqual(int(np.sum(counts==1)),6*factor)
                self.assertEqual(len(result.xyz), (2*factor+1)*(factor+1))
                self.assertTrue(np.all(result.xyz[unique[counts==1],2] == 5.))
                self.assertGreater(float(result.xyz[:,2].max()),5.1)
                self.assertLessEqual(float(result.xyz[:,2].max()),5.2)

    def test_original_vertices_stay_exact_and_pattern_is_reused(self):
        source = self.mesh()
        patterns = PatternLibrary()
        builder = UniformQuadSubdivision(patterns)
        first = builder.build(source, 2, self.field, datum=5., ceiling=5.2)
        for point in source.xyz:
            self.assertTrue(np.any(np.all(first.xyz==point,axis=1)))
        second = builder.build(source, 2, self.field, datum=5., ceiling=5.2)
        np.testing.assert_array_equal(first.xyz, second.xyz)
        np.testing.assert_array_equal(first.quads, second.quads)
        self.assertGreater(patterns.stats()['hits'],0)

    def test_cell_budget_is_checked_before_field_evaluation(self):
        called = []
        with self.assertRaisesRegex(ValueError,'cell budget'):
            UniformQuadSubdivision(PatternLibrary(),max_cells=31).build(self.mesh(),4,
                lambda xy:called.append(True),datum=5.,ceiling=5.2)
        self.assertFalse(called)

    def test_identity_and_invalid_factors(self):
        source = self.mesh()
        builder = UniformQuadSubdivision(PatternLibrary())
        self.assertIs(builder.build(source,1,self.field,datum=5.,ceiling=5.2),source)
        for factor in (0,5,1.5):
            with self.assertRaises(ValueError):
                builder.build(source,factor,self.field,datum=5.,ceiling=5.2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
