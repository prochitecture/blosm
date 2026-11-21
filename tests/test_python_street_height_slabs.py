"""Directional constraints are optional and preserve the physical height box."""
import unittest
from types import SimpleNamespace
import numpy as np
import _python_street_package
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.entrance_height_recovery import entrance_slabs,recover
from tested_python_streets.height_slabs import HeightSlabs
from tested_python_streets.mesh import Mesh
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.slope_projection import SlopeProjection,ProjectionIterationBudgetError
from tested_python_streets.surface_query import QuadHeightField
from test_python_street_height_projection import peak_grid


class HeightSlabTests(unittest.TestCase):
    def test_one_free_height_has_the_known_slab_projection_and_none_is_unchanged(self):
        mesh,fixed=peak_grid()
        slab=HeightSlabs(mesh.quads[:1],np.array(((0.,0.,1.,-1.),)),np.array((.02,)))
        kwargs=dict(datum=0.,ceiling=.08,maximum_slope=.1)
        original,old_report=SlopeProjection().solve(mesh,fixed,**kwargs)
        unchanged,none_report=SlopeProjection().solve(mesh,fixed,height_slabs=None,**kwargs)
        np.testing.assert_array_equal(original.xyz,unchanged.xyz)
        self.assertEqual(old_report,none_report)
        result,report=SlopeProjection().solve(mesh,fixed,height_slabs=slab,**kwargs)
        self.assertAlmostEqual(result.xyz[4,2],.02,places=10)
        self.assertLessEqual(slab.violation(result.xyz[:,2]),1e-9)
        self.assertLessEqual(float(result.quad_slopes().max()),.1+1e-9)
        np.testing.assert_array_equal(result.xyz[fixed],mesh.xyz[fixed])
        np.testing.assert_array_equal(result.xyz[:,:2],mesh.xyz[:,:2])
        self.assertEqual(report['directional_height_rows'],1)
        horizontal=HeightSlabs(mesh.quads[:1],np.array(((0.,0.,1.,-1.),)),np.array((0.,)))
        level,_=SlopeProjection().solve(mesh,fixed,height_slabs=horizontal,**kwargs)
        self.assertAlmostEqual(level.xyz[4,2],0.,places=12)
        with self.assertRaises(ProjectionIterationBudgetError):
            SlopeProjection(max_iterations=0).solve(mesh,fixed,height_slabs=slab,**kwargs)

    def test_complete_entrance_edge_and_spill_checks_follow_the_coupled_solve(self):
        mesh,fixed=peak_grid();mesh.xyz[3,2]=.04
        conditions=(SimpleNamespace(xy=mesh.xyz[:,:2].copy(),ports=(SimpleNamespace(vertices=(0,3,6)),)),)
        field=QuadHeightField(mesh,datum=0.,ceiling=.08)
        builder=CompoundSurfaceBuilder(PatternLibrary())
        result,report=recover(builder,mesh,field,fixed,conditions,.1,100)
        slabs=entrance_slabs(result,fixed,conditions,.03)
        self.assertLessEqual(slabs.violation(result.xyz[:,2]),1e-9)
        self.assertLessEqual(report['joins']['maximum_entrance_longitudinal_slope'],.03)
        self.assertEqual(report['rounds'][-1]['drainage']['maximum_fill'],0.)
        np.testing.assert_array_equal(result.xyz[fixed],mesh.xyz[fixed])
        np.testing.assert_array_equal(result.quads,mesh.quads)

    def test_fixed_profile_infeasibility_and_allocation_are_rejected(self):
        mesh=Mesh(np.array(((0.,0.,0.),(1.,0.,.04),(1.,1.,.04),(0.,1.,0.))),np.array(((0,1,2,3),)))
        slab=HeightSlabs(mesh.quads,np.array(((-1.,1.,0.,0.),)),np.array((.01,)))
        with self.assertRaisesRegex(ValueError,'Fixed profile heights'):
            SlopeProjection().solve(mesh,np.ones(4,bool),datum=0.,ceiling=.08,maximum_slope=.1,height_slabs=slab)
        with self.assertRaisesRegex(ValueError,'mesh/memory budget'):
            SlopeProjection(memory_bytes=4096).solve(mesh,np.zeros(4,bool),datum=0.,ceiling=.08,
                maximum_slope=.1,height_slabs=HeightSlabs(np.repeat(mesh.quads,8,axis=0),
                    np.tile((-1.,1.,0.,0.),(8,1)),np.full(8,.01)))
        with self.assertRaises(ValueError):HeightSlabs(mesh.quads,np.zeros((1,4)),np.array((float('nan'),)))


if __name__=='__main__':unittest.main()
