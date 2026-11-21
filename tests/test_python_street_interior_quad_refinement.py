"""Physical-rim and direct-quad contracts of the bounded atlas fallback."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.interior_quad_refinement import InteriorQuadRefinement,InteriorQuadRefinementError
from tested_python_streets.mesh import Mesh,publication_coordinates
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.quality_correction import QuadShapeCorrection


def distorted_grid():
    xy=np.array([(x,y) for y in range(5) for x in range(5)],float)
    xy[[11,12,13],1]=1.05
    faces=np.array([(5*y+x,5*y+x+1,5*(y+1)+x+1,5*(y+1)+x) for y in range(4) for x in range(4)])
    rim=np.array([0,1,2,3,4,9,14,19,24,23,22,21,20,15,10,5])
    return Mesh(np.c_[xy,np.zeros(len(xy))],faces,seams={0:rim})


class InteriorQuadRefinementTests(unittest.TestCase):
    def test_distorted_interior_becomes_qualified_direct_quads_with_exact_rim_and_vertex_lineage(self):
        mesh=distorted_grid();before=mesh.xyz.copy();faces=mesh.quads.copy()
        _,initial,_=QuadShapeCorrection().evaluate(mesh.xyz[:,:2],mesh.quads)
        self.assertGreater(initial['maximum_edge_ratio'],8.)
        result,mapping,report=InteriorQuadRefinement(PatternLibrary()).improve(mesh)
        self.assertTrue(report['qualified'])
        self.assertGreater(len(result.quads),len(mesh.quads))
        self.assertEqual(result.quads.shape[1],4)
        self.assertEqual(len(mapping),len(mesh.xyz))
        self.assertEqual(len(set(mapping)),len(mapping))
        np.testing.assert_array_equal(result.xyz[result.seams[0]],before[mesh.seams[0]])
        np.testing.assert_array_equal(result.seams[0],mapping[mesh.seams[0]])
        np.testing.assert_array_equal(mesh.xyz,before);np.testing.assert_array_equal(mesh.quads,faces)
        self.assertEqual(report['physical_boundary_edges_split'],0)
        stored=publication_coordinates(result.xyz,1.,physical=True)
        _,quality,_=QuadShapeCorrection().evaluate(stored[:,:2],result.quads)
        self.assertTrue(QuadShapeCorrection().qualified(quality))

    def test_qualified_input_is_unchanged_and_an_immutable_skinny_single_quad_stays_unsupported(self):
        mesh=Mesh(np.array(((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.))),np.array(((0,1,2,3),)))
        result,mapping,report=InteriorQuadRefinement(PatternLibrary()).improve(mesh)
        self.assertIs(result,mesh);np.testing.assert_array_equal(mapping,np.arange(4));self.assertFalse(report['passes'])
        thin=Mesh(mesh.xyz*np.array((20.,1.,1.)),mesh.quads)
        with self.assertRaises(InteriorQuadRefinementError) as caught:
            InteriorQuadRefinement(PatternLibrary()).improve(thin)
        self.assertIn('immutable',caught.exception.metadata['refinement_failure'])
        self.assertFalse(caught.exception.metadata['qualified'])

    def test_allocation_and_nonplanar_inputs_fail_without_changing_the_rim(self):
        mesh=distorted_grid();before=mesh.xyz.copy()
        with self.assertRaises(InteriorQuadRefinementError) as caught:
            InteriorQuadRefinement(PatternLibrary(),max_cells=len(mesh.quads)).improve(mesh)
        self.assertIn('cell budget',caught.exception.metadata['refinement_failure'])
        np.testing.assert_array_equal(mesh.xyz,before)
        moved=mesh.xyz.copy();moved[0,2]=.1
        with self.assertRaisesRegex(ValueError,'planar'):
            InteriorQuadRefinement(PatternLibrary()).improve(Mesh(moved,mesh.quads))
        with self.assertRaisesRegex(ValueError,'allocation'):
            InteriorQuadRefinement(PatternLibrary(),max_passes=4)


if __name__=='__main__':unittest.main()
