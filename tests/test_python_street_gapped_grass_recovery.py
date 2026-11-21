"""Stored-only topology and height recovery retain all actual grass edges."""
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_positive_quad_patches import fixture, boundary_edges
import test_python_street_finite_free_closure as closure_fixture
from tested_python_streets.gapped_grass_recovery import recover
from tested_python_streets.mesh import Mesh
from tested_python_streets.fixed_grass_heights import boundary_mask
from tested_python_streets.publication_slopes import PublicationSlopeGuard


class GappedGrassRecoveryTests(unittest.TestCase):
    field = SimpleNamespace(datum=0., ceiling=.08)

    @staticmethod
    def constrained_ear(angle=0.,shift=(0.,0.)):
        xy=np.array(((-1.,0.),(0.,.05),(1.,0.),(1.05,-1.),
            (1.,-2.),(0.,-2.05),(-1.,-2.),(-1.05,-1.),(0.,-1.)))
        rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
        xyz=np.c_[xy@rotation.T+shift,np.full(len(xy),.08)]
        xyz[1,2]=.074
        return Mesh(xyz,np.array(((8,2,1,0),(8,4,3,2),(8,6,5,4),(8,0,7,6))), 'grass')

    def test_fixed_height_ear_uses_direct_quads_without_changing_the_physical_rim(self):
        for angle,shift in ((0.,(0.,0.)),(.713,(1200.,-1800.))):
            with self.subTest(angle=angle):
                mesh=self.constrained_ear(angle,shift);before=mesh.xyz.copy()
                fixed=boundary_mask(mesh);keys=[('source',i) for i in range(len(mesh.xyz))]
                self.assertGreater(PublicationSlopeGuard.fixed_gradient_floor(mesh,fixed,0.),.1)
                output,names,report=recover(mesh,keys,self.field,namespace='fixed ear')
                np.testing.assert_array_equal(mesh.xyz,before)
                np.testing.assert_array_equal(output.xyz[:len(before),:2],before[:,:2])
                np.testing.assert_array_equal(output.xyz[:len(before)][fixed],before[fixed])
                np.testing.assert_array_equal(boundary_edges(output)[0],boundary_edges(mesh)[0])
                self.assertEqual(names[:len(keys)],keys)
                self.assertFalse(report['topology_patch']['triangles_created'])
                self.assertEqual(output.quads.shape[1],4)
                self.assertLessEqual(len(output.quads),len(mesh.quads)+1)
                self.assertLessEqual(output.with_positions(output.xyz.astype(np.float32)).validate()['maximum_edge_ratio'],8.)
                self.assertLessEqual(report['height_reconciliation']['maximum_stored_physical_slope'],.1+1e-7)
                self.assertLessEqual(report['height_reconciliation']['drainage']['maximum_fill'],1e-6)

    def test_height_ear_cannot_increase_the_cell_budget_or_change_an_impossible_boundary(self):
        mesh=self.constrained_ear();keys=[('source',i) for i in range(len(mesh.xyz))]
        with self.assertRaises(ValueError):recover(mesh,keys,self.field,namespace='budget',max_cells=4)
        impossible=mesh.with_positions(mesh.xyz.copy());impossible.xyz[0,2]=0.
        impossible.xyz[:,:2]*=.2  # A physical boundary edge alone now exceeds ten percent.
        with self.assertRaises(ValueError):recover(impossible,keys,self.field,namespace='impossible')

    def test_existing_success_returns_the_same_mesh_and_keys(self):
        mesh = Mesh(np.array(((0., 0., .04), (1., 0., .04), (1., 1., .04), (0., 1., .04))),
            np.array(((0, 1, 2, 3),)), 'grass')
        keys = [('source', i) for i in range(len(mesh.xyz))]
        output, names, details = recover(mesh, keys, self.field, namespace='old')
        self.assertIs(output, mesh)
        self.assertIs(names, keys)
        self.assertIsNone(details)

    def test_stored_singular_patch_keeps_exact_boundary_and_no_source_geometry_changes(self):
        mesh = fixture(shift=(1200., 1800.))
        mesh.xyz[1, 1] -= 1e-7
        mesh.xyz[:, 2] = .04
        mesh.material = 'grass'
        keys = [('source', i) for i in range(len(mesh.xyz))]
        source = mesh.xyz.copy()
        output, names, report = recover(mesh, keys, self.field, namespace='stored')
        self.assertTrue(report['topology_patch']['admitted_stored_failures'])
        self.assertFalse(report['physical_boundary_changed'])
        np.testing.assert_array_equal(mesh.xyz, source)
        np.testing.assert_array_equal(output.xyz[:len(source)], source)
        np.testing.assert_array_equal(boundary_edges(output)[0], boundary_edges(mesh)[0])
        self.assertEqual(names[:len(keys)], keys)
        self.assertEqual(output.quads.shape[1], 4)
        stored = output.with_positions(output.xyz.astype(np.float32))
        self.assertLessEqual(stored.validate()['maximum_edge_ratio'], 8.)
        self.assertLessEqual(stored.quad_slopes().max(), .1+1e-7)

    def test_height_only_recovery_preserves_every_xy_quad_and_source_boundary(self):
        mesh, fixed = closure_fixture.FiniteFreeClosureTests.height_seed()
        keys = [('source', i) for i in range(len(mesh.xyz))]
        output, names, report = recover(mesh, keys, self.field, namespace='height')
        self.assertIsNone(report['topology_patch'])
        np.testing.assert_array_equal(output.xyz[:, :2], mesh.xyz[:, :2])
        np.testing.assert_array_equal(output.xyz[fixed], mesh.xyz[fixed])
        np.testing.assert_array_equal(output.quads, mesh.quads)
        self.assertEqual(names, keys)
        self.assertLessEqual(report['height_reconciliation']['maximum_stored_physical_slope'], .1+1e-7)


if __name__ == '__main__':unittest.main()
