"""Explicit normal Product construction retains source heights and body rows."""
from dataclasses import replace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from test_python_street_source_cycles import triangle,prepare
from test_python_street_publication_entrances import fixture as publication_fixture
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.ownership import CompoundProduct
from tested_python_streets.height_slabs import HeightSlabs
from tested_python_streets.source_body_heights import prepare as body_guides
from tested_python_streets.source_body_references import SourceBodyUnavailable
from tested_python_streets.publication_slopes import PublicationSlopeGuard
from tested_python_streets.surface_query import QuadHeightField
from tested_python_streets.patterns import PatternLibrary


def actual_boundary(mesh):
    edges=np.sort(np.stack((mesh.quads,np.roll(mesh.quads,-1,axis=1)),axis=-1).reshape(-1,2),axis=1)
    edges,counts=np.unique(edges,axis=0,return_counts=True)
    return np.unique(edges[counts==1])


class SourceHeightProductTests(unittest.TestCase):
    def test_retained_crossfall_uses_only_the_owned_sections_and_physical_width(self):
        from tested_python_streets.model import Key
        from tested_python_streets.profiles import Crown
        context, group = prepare(triangle(8.))
        retained = group.claims[0].section_id
        context.products[Key('profile', retained)] = Crown(.32, .5)
        # A separate approach is outside the retained-profile calculation.
        context.products[Key('profile', 4)] = Crown(1., .5)
        builder = CompoundSurfaceBuilder(context.patterns, source_heights=True)
        self.assertAlmostEqual(builder.retained_crossfall(group, context), .16)
        for half in (-1, 1):
            context.products[Key('profile', retained)] = Crown(.32, .5, half)
            self.assertAlmostEqual(builder.retained_crossfall(group, context), .08)
        self.assertEqual(builder.retained_crossfall(replace(group, claims=()), context), 0.)

    def test_normal_product_preserves_actual_body_targets_and_rebuilds_its_physical_query(self):
        context,group=prepare(triangle(8.))
        before=context.network.to_dict()
        product=CompoundProduct.build(context,group.members[0],group,source_heights=True)
        self.assertEqual(context.network.to_dict(),before)
        self.assertTrue(product.pavement.metadata['source_height_conditions'])
        guide=body_guides(context,product.group,product.physical_pavement,actual_boundary(product.physical_pavement))
        self.assertGreater(len(guide.core_nodes),0)
        self.assertFalse(guide.physical_fixed[guide.core_nodes].any())
        np.testing.assert_allclose(product.physical_pavement.xyz[guide.core_nodes,2],guide.core_targets,atol=1e-8,rtol=0.)
        points=product.physical_pavement.xyz[guide.core_nodes,:2]
        np.testing.assert_allclose(product.field(points),guide.core_targets,atol=1e-8,rtol=0.)
        expected=product.height(points)
        product.release_geometry()
        np.testing.assert_array_equal(product.height(points),expected)

    def test_normal_product_zero_hole_and_body_physical_heights_survive_zero_and_tenfold_exaggeration(self):
        physical=[];roundoff=[];fixed=None
        for exaggeration in (0.,10.):
            network=triangle(2.)
            network=replace(network,settings=replace(network.settings,exaggeration=exaggeration))
            context,group=prepare(network)
            product=CompoundProduct.build(context,group.members[0],group,source_heights=True)
            physical.append(product.physical_pavement.xyz.copy())
            roundoff.append(product.pavement.metadata['publication_roundoff'])
            guide=body_guides(context,product.group,product.physical_pavement,actual_boundary(product.physical_pavement))
            fixed=guide.physical_fixed
            np.testing.assert_allclose(product.physical_pavement.xyz[guide.core_nodes,2],
                                       guide.core_targets,atol=1e-8,rtol=0.)
            np.testing.assert_allclose(product.field(product.physical_pavement.xyz[guide.core_nodes,:2]),
                                       guide.core_targets,atol=1e-8,rtol=0.)
            self.assertIsNotNone(product.physical_grass)
            np.testing.assert_array_equal(product.physical_grass.xyz[:,2],0.)
            np.testing.assert_array_equal(product.grass.mesh.xyz[:,2],0.)
            li=next(i for i,loop in enumerate(product.group.loops) if loop.hole)
            np.testing.assert_array_equal(product.physical_rims[li][:,2],0.)
            rim=product.physical_rims[li]
            np.testing.assert_allclose(product.hole_fields[li](rim[:,:2]),rim[:,2],atol=1e-12,rtol=0.)
            point=product.physical_grass.xyz[product.physical_grass.quads[0],:2].mean(0)[None,:]
            np.testing.assert_array_equal(product.height(point),0.)
            mesh=product.physical_pavement
            np.testing.assert_array_equal(product.pavement.mesh.xyz[:,2],exaggeration*mesh.xyz[:,2])
            product.release_geometry()
            np.testing.assert_array_equal(product.height(point),0.)
        np.testing.assert_array_equal(physical[0][:,:2],physical[1][:,:2])
        np.testing.assert_array_equal(physical[0][fixed],physical[1][fixed])
        self.assertIsNone(roundoff[0])
        self.assertIsNotNone(roundoff[1])
        # Stored-coordinate certification may adjust free physical heights.
        # It must retain the original-scale source body targets and boundary.
        self.assertLessEqual(float(abs(physical[0][:,2]-physical[1][:,2]).max()),
                             roundoff[1]['maximum_height_change']+1e-8)

    def test_default_policy_is_unchanged_and_invalid_or_unsupported_opt_in_does_not_fall_back(self):
        context,group=prepare(triangle(8.))
        implicit=CompoundProduct.build(context,group.members[0],group)
        context,group=prepare(triangle(8.))
        explicit=CompoundProduct.build(context,group.members[0],group,source_heights=False)
        np.testing.assert_array_equal(implicit.pavement.mesh.xyz,explicit.pavement.mesh.xyz)
        np.testing.assert_array_equal(implicit.pavement.mesh.quads,explicit.pavement.mesh.quads)
        self.assertNotIn('source_height_conditions',implicit.pavement.metadata)
        self.assertNotIn('source_body_heights',implicit.pavement.metadata)
        for value in (1,None,'source'):
            with self.assertRaises(ValueError):CompoundSurfaceBuilder(context.patterns,source_heights=value)
            with self.assertRaises(ValueError):CompoundProduct.build(context,group.members[0],group,source_heights=value)
        with patch('tested_python_streets.source_body_heights.prepare_steps',
                   side_effect=SourceBodyUnavailable('Original paired field has no certified finite coordinates')):
            with self.assertRaisesRegex(SourceBodyUnavailable,'certified finite coordinates'):
                CompoundProduct.build(context,group.members[0],group,source_heights=True)

    def test_publication_fallback_retains_exact_body_rows_and_physical_only_outlets(self):
        mesh,_,fixed,conditions,rejected=publication_fixture()
        mesh.xyz[6,2]=mesh.xyz[4,2]
        rejected.xyz[6,2]=mesh.xyz[6,2]
        uv=np.arange(len(mesh.quads)*8,dtype=float).reshape(-1,4,2)/10.
        mesh=mesh.with_uv_layer('Authored UV',uv)
        field=QuadHeightField(mesh,datum=0.,ceiling=.08)
        body=HeightSlabs(np.array(((4,6,7,8),)),np.array(((1.,-1.,0.,0.),)),np.array((0.,)))
        builder=CompoundSurfaceBuilder(PatternLibrary())
        guard=PublicationSlopeGuard();ordinary=guard.correct_steps;calls=[];masks=[]
        def publication(candidate,query,mask,**kwargs):
            calls.append(kwargs);masks.append(mask.copy())
            if len(calls)==1:kwargs['check_joins'](rejected)
            return ordinary(candidate,query,mask,**kwargs)
        with patch.object(guard,'correct_steps',side_effect=publication),patch(
                'tested_python_streets.compound_surfaces.PublicationSlopeGuard',return_value=guard):
            result,query,report=builder.correct_publication(mesh,field,fixed,conditions,
                exaggeration=1.,limit=.1,max_cells=100,additional_slabs=body)
        self.assertEqual(len(calls),2)
        for call,mask in zip(calls,masks):
            rows=call['height_slabs']
            found=[i for i,row in enumerate(rows.vertices) if np.array_equal(row,body.vertices[0])]
            self.assertEqual(len(found),1)
            np.testing.assert_array_equal(rows.coefficients[found],body.coefficients)
            np.testing.assert_array_equal(mask,fixed)
            self.assertFalse(mask[4])
        self.assertLessEqual(body.violation(result.xyz[:,2]),1e-8)
        self.assertEqual(report['directional_entrance_recovery']['additional_height_rows'],1)
        self.assertLessEqual(report['stored_joins']['maximum_entrance_longitudinal_slope'],.03)
        np.testing.assert_array_equal(result.xyz[fixed],mesh.xyz[fixed])
        np.testing.assert_array_equal(result.xyz[:,:2],mesh.xyz[:,:2])
        np.testing.assert_array_equal(result.quads,mesh.quads)
        np.testing.assert_array_equal(result.uv_layers['Authored UV'],mesh.uv_layers['Authored UV'])
        np.testing.assert_allclose(query(result.xyz[:,:2]),result.xyz[:,2],atol=1e-10,rtol=0.)


if __name__=='__main__':unittest.main()
