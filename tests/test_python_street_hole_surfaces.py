"""Closed road blocks and median holes retain their actual curb heights."""
from types import SimpleNamespace
from dataclasses import replace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.hole_surfaces import build, ConstantHoleHeight
from tested_python_streets.mesh import Mesh
from tested_python_streets.surface_query import QuadHeightField
from tested_python_streets.ownership import CompoundProduct


def grid(nx=8, ny=6, *, datum=0., boundary_height=lambda x,y:0.):
    xy = np.array([(float(i),float(j)) for j in range(ny+1) for i in range(nx+1)])
    quads = np.array([(j*(nx+1)+i,j*(nx+1)+i+1,(j+1)*(nx+1)+i+1,(j+1)*(nx+1)+i)
                      for j in range(ny) for i in range(nx)])
    rim = np.array([*range(nx+1),*(j*(nx+1)+nx for j in range(1,ny+1)),
                    *(ny*(nx+1)+i for i in range(nx-1,-1,-1)),
                    *(j*(nx+1) for j in range(ny-1,0,-1))])
    template = Mesh(np.column_stack((xy,np.zeros(len(xy)))),quads,'grass',{0:rim})
    uv = np.arange(len(quads)*8,dtype=float).reshape(-1,4,2)/100.
    template = template.with_uv_layer('Authored UV',uv)
    xyz = template.xyz.copy()
    xyz[:,2] = datum+boundary_height(xy[:,0],xy[:,1])
    pavement = Mesh(xyz,quads,seams={1:rim[::-1]})
    return template,pavement


class SourceHoleSurfaceTests(unittest.TestCase):
    def test_grass_corner_can_fit_a_gentle_rim_with_independent_adjacent_quads(self):
        xyz = np.array([[-1., 0., 0.], [0., -.02, 0.], [1., 0., 0.], [0., 1., 0.]])
        face = np.array([[0, 1, 2, 3]])
        template = Mesh(xyz, face, 'grass', {0:face[0]}, {'Map':xyz[face,:2]})
        raised = xyz.copy()
        raised[:,2] = (0., .005, .025, .01)
        pavement = Mesh(raised, face, seams={1:face[0,::-1]})
        with self.assertRaisesRegex(ValueError, 'slope limit'):
            build(template, pavement, 1, datum=0., ceiling=.08, max_cells=4)
        result = build(template, pavement, 1, datum=0., ceiling=.08, max_cells=5)
        self.assertEqual(len(result.mesh.quads), 5)
        np.testing.assert_array_equal(result.mesh.xyz[face[0]], raised)
        np.testing.assert_allclose(result.mesh.uv_layers['Map'], result.mesh.xyz[result.mesh.quads,:2])
        self.assertLessEqual(result.mesh.quad_slopes().max(), .100001)
        self.assertLessEqual(result.mesh.validate()['maximum_edge_ratio'], 8.)

    def test_zero_curb_block_stays_zero_instead_of_growing_a_crown_mound(self):
        template,pavement=grid()
        before=template.xyz.copy()
        result=build(template,pavement,1,datum=0.,ceiling=.08)
        np.testing.assert_array_equal(result.mesh.xyz[:,2],0.)
        np.testing.assert_array_equal(result.field(np.array(((1.,1.),(4.,3.)))),0.)
        self.assertIsInstance(result.field,ConstantHoleHeight)
        np.testing.assert_array_equal(template.xyz,before)
        np.testing.assert_array_equal(result.mesh.quads,template.quads)
        np.testing.assert_array_equal(result.mesh.uv_layers['Authored UV'],template.uv_layers['Authored UV'])

    def test_real_constant_plateau_keeps_old_buffers_and_physical_queries(self):
        for datum in (0.,25.):
            template,pavement=grid(datum=datum,boundary_height=lambda x,y:.08+np.zeros_like(x))
            value=float(pavement.xyz[0,2])
            for exaggeration in (0.,1.,10.):
                result=build(template,pavement,1,datum=datum,ceiling=datum+.08,exaggeration=exaggeration)
                expected=template.xyz.copy();expected[:,2]=value
                np.testing.assert_array_equal(result.mesh.xyz,expected)
                np.testing.assert_array_equal(result.field(np.array(((4.,3.),))),value)

    def test_mixed_boundary_retains_every_curb_and_fixed_slope_limits(self):
        template,pavement=grid(boundary_height=lambda x,y:.01*x)
        for exaggeration in (0.,1.,10.):
            result=build(template,pavement,1,datum=0.,ceiling=.08,exaggeration=exaggeration)
            rim=template.seams[0]
            np.testing.assert_array_equal(result.mesh.xyz[rim],pavement.xyz[pavement.seams[1]][::-1])
            np.testing.assert_array_equal(result.mesh.xyz[:,:2],template.xyz[:,:2])
            np.testing.assert_array_equal(result.mesh.quads,template.quads)
            self.assertTrue(np.all((result.mesh.xyz[:,2]>=0.)&(result.mesh.xyz[:,2]<=.08)))
            self.assertLessEqual(float(result.mesh.quad_slopes().max()),.1+1e-7)
            self.assertLessEqual(result.metadata['reconciliation']['maximum_stored_physical_slope'],.1+1e-7)
            np.testing.assert_allclose(result.field(result.mesh.xyz[rim,:2]),result.mesh.xyz[rim,2],atol=1e-10,rtol=0.)
            np.testing.assert_array_equal(result.mesh.uv_layers['Authored UV'],template.uv_layers['Authored UV'])

    def test_real_source_cycle_product_keeps_zero_hole_and_queries_after_release(self):
        from test_python_street_source_cycles import triangle,prepare
        from tested_python_streets.conditions import BoundaryConditions
        original=BoundaryConditions.build
        for exaggeration in (0.,10.):
            network=triangle(2.)
            network=replace(network,settings=replace(network.settings,exaggeration=exaggeration))
            context,group=prepare(network)
            def source_conditions(self,value,loop,layouts):
                return original(self,value,loop,layouts,context=context)
            with patch.object(BoundaryConditions,'build',source_conditions):
                product=CompoundProduct.build(context,group.members[0],group)
            self.assertIsNotNone(product.physical_grass)
            np.testing.assert_array_equal(product.physical_grass.xyz[:,2],0.)
            np.testing.assert_array_equal(product.grass.mesh.xyz[:,2],0.)
            li=next(i for i,loop in enumerate(product.group.loops) if loop.hole)
            np.testing.assert_array_equal(product.physical_rims[li][:,2],0.)
            point=product.physical_grass.xyz[product.physical_grass.quads[0],:2].mean(0)[None,:]
            np.testing.assert_array_equal(product.height(point),0.)
            product.release_geometry()
            np.testing.assert_array_equal(product.height(point),0.)

    def test_mismatched_or_missing_boundary_and_cell_budget_are_fatal(self):
        template,pavement=grid()
        with self.assertRaisesRegex(ValueError,'cell budget'):
            build(template,pavement,1,datum=0.,ceiling=.08,max_cells=1)
        pavement.seams[1]=pavement.seams[1][1:]
        with self.assertRaisesRegex(ValueError,'match exactly'):
            build(template,pavement,1,datum=0.,ceiling=.08)
        template,pavement=grid();template.seams[0]=template.seams[0][1:]
        pavement.seams[1]=pavement.seams[1][:-1]
        with self.assertRaisesRegex(ValueError,'omits a physical boundary'):
            build(template,pavement,1,datum=0.,ceiling=.08)

    def test_owner_queries_actual_hole_surface_after_release_and_exaggeration(self):
        template,pavement=grid(3,3,boundary_height=lambda x,y:.02*x)
        field=QuadHeightField(Mesh(pavement.xyz,np.delete(pavement.quads,4,axis=0)),datum=0.,ceiling=.08)
        hole_xy=np.array(((1.,1.),(2.,1.),(2.,2.),(1.,2.)))
        hole=Mesh(np.column_stack((hole_xy,.02*hole_xy[:,0])),np.array(((0,1,2,3),)))
        hole_field=QuadHeightField(hole,datum=0.,ceiling=.08)
        loops=(SimpleNamespace(hole=False,xy=np.array(((0.,0.),(3.,0.),(3.,3.),(0.,3.)))),
               SimpleNamespace(hole=True,xy=hole_xy[::-1]))
        for exaggeration in (0.,1.,10.):
            owner=CompoundProduct(None,SimpleNamespace(loops=loops),field,{}, {}, {}, {},None,None,
                                  exaggeration,hole_fields={1:hole_field})
            owner.release_geometry()
            xy=np.array(((.5,.5),(1.5,1.5),(1.,1.5),(5.,5.)))
            np.testing.assert_allclose(owner.height(xy),exaggeration*np.array((.01,.03,.02,0.)),atol=1e-12)
        owner.hole_fields=None
        with self.assertRaisesRegex(ValueError,'missing its actual'):
            owner.height(np.array(((1.5,1.5),)))

    def test_height_policy_epoch_rebuilds_surface_without_rebuilding_source_geometry(self):
        from _python_street_fixtures import fixture_networks
        from tested_python_streets import model
        from tested_python_streets.backend import PythonMeshBackend
        network=fixture_networks(model)[0]
        backend=PythonMeshBackend(network.settings)
        initial=backend.prepare(network);backend.commit(initial)
        sid=network.sections[0].id
        path=initial.products[model.Key('path',sid)]
        surfaces={key:product for key,product in initial.products.items() if key.kind=='junction_surface'}
        self.assertTrue(surfaces)
        with patch('tested_python_streets.hole_surfaces.SOURCE_HOLE_HEIGHT_POLICY','next-height-policy'):
            changed=backend.prepare(network)
        self.assertIs(changed.products[model.Key('path',sid)],path)
        for key,product in surfaces.items():
            self.assertIn(key,changed.rebuilt)
            if product is not None:
                np.testing.assert_array_equal(changed.products[key].pavement.mesh.xyz,product.pavement.mesh.xyz)
                np.testing.assert_array_equal(changed.products[key].pavement.mesh.quads,product.pavement.mesh.quads)


if __name__=='__main__':
    unittest.main()
