"""Only an actual publication-join failure enters directional reconciliation."""
from dataclasses import dataclass
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder,PublicationEntranceSlopeError
from tested_python_streets.publication_slopes import PublicationSlopeGuard
from tested_python_streets.surface_query import QuadHeightField
from tested_python_streets.patterns import PatternLibrary
from test_python_street_height_projection import peak_grid


@dataclass
class Boundary:
    xy:object
    ports:tuple


def fixture():
    mesh,fixed=peak_grid()
    mesh.xyz[:,:2]=mesh.xyz[:,:2]*.1+(12345.,-12345.)
    mesh.xyz[4,2]=.0002
    field=QuadHeightField(mesh,datum=0.,ceiling=.08)
    boundary=Boundary(mesh.xyz[:,:2].copy(),(SimpleNamespace(vertices=(0,1,2)),))
    rejected=mesh.with_positions(mesh.xyz.copy());rejected.xyz[4,2]=.0069
    return mesh,field,fixed,(boundary,),rejected


class PublicationEntranceTests(unittest.TestCase):
    def test_midpoint_rows_equal_the_authoritative_join_check_on_skew_quads(self):
        from tested_python_streets.entrance_height_recovery import entrance_slabs
        mesh,field,fixed,conditions,rejected=fixture()
        rejected.xyz[4,:2]+=(.012,-.018)
        rows=entrance_slabs(rejected,fixed,conditions,0.,at_corners=False)
        actual=rows.violation(rejected.xyz[:,2])
        expected=CompoundSurfaceBuilder.projection_checks(rejected,fixed,conditions)
        self.assertAlmostEqual(actual,expected['maximum_entrance_longitudinal_slope'],places=13)
        self.assertEqual(len(rows.vertices),2)
        corners=entrance_slabs(rejected,fixed,conditions,0.)
        self.assertEqual(len(corners.vertices),4)

    def test_one_failed_join_reconciles_real_slabs_without_moving_source_geometry(self):
        mesh,field,fixed,conditions,rejected=fixture()
        builder=CompoundSurfaceBuilder(PatternLibrary());guard=PublicationSlopeGuard()
        ordinary=guard.correct_steps;calls=[]
        def publication(candidate,query,mask,**kwargs):
            calls.append(kwargs)
            if len(calls)==1:kwargs['check_joins'](rejected)
            return ordinary(candidate,query,mask,**kwargs)
        before=mesh.xyz.copy()
        with patch.object(guard,'correct_steps',side_effect=publication),patch(
                'tested_python_streets.compound_surfaces.PublicationSlopeGuard',return_value=guard):
            result,query,report=builder.correct_publication(mesh,field,fixed,conditions,
                exaggeration=1.,limit=.1,max_cells=100)
        self.assertEqual(len(calls),2)
        self.assertNotIn('height_slabs',calls[0]);self.assertIn('height_slabs',calls[1])
        self.assertGreater(report['failed_publication_joins']['maximum_entrance_longitudinal_slope'],.03)
        self.assertLessEqual(report['stored_joins']['maximum_entrance_longitudinal_slope'],.03)
        self.assertEqual(report['stored_joins']['isolated_interior_minima'],0)
        self.assertEqual(report['maximum_reconciliation_attempts'],1)
        self.assertTrue(report['exact_boundary'] and report['exact_xy'] and report['exact_quads'])
        self.assertGreater(report['maximum_height_change'],0.)
        self.assertEqual(query.metadata['publication_roundoff'],report)
        np.testing.assert_array_equal(result.xyz[fixed],before[fixed])
        np.testing.assert_array_equal(result.xyz[:,:2],before[:,:2])
        np.testing.assert_array_equal(result.quads,mesh.quads)
        np.testing.assert_array_equal(mesh.xyz,before)
        np.testing.assert_allclose(query(result.xyz[:,:2]),result.xyz[:,2],atol=1e-10,rtol=0.)

    def test_unrelated_publication_failure_and_changed_fixed_geometry_still_reject(self):
        mesh,field,fixed,conditions,rejected=fixture();builder=CompoundSurfaceBuilder(PatternLibrary())
        with patch('tested_python_streets.compound_surfaces.PublicationSlopeGuard.correct_steps',
                side_effect=ValueError('unrelated publication failure')):
            with self.assertRaisesRegex(ValueError,'unrelated publication failure'):
                builder.correct_publication(mesh,field,fixed,conditions,exaggeration=1.,limit=.1,max_cells=100)
        rejected.xyz[0,2]=.0001
        def changed(candidate,query,mask,**kwargs):kwargs['check_joins'](rejected)
        with patch('tested_python_streets.compound_surfaces.PublicationSlopeGuard.correct_steps',side_effect=changed):
            with self.assertRaisesRegex(ValueError,'unchanged feasible source geometry'):
                builder.correct_publication(mesh,field,fixed,conditions,exaggeration=1.,limit=.1,max_cells=100)

    def test_disabled_height_projection_never_enables_directional_retry(self):
        mesh,field,fixed,conditions,rejected=fixture()
        builder=CompoundSurfaceBuilder(PatternLibrary(),project_heights=False)
        def failed(candidate,query,mask,**kwargs):kwargs['check_joins'](rejected)
        with patch('tested_python_streets.compound_surfaces.PublicationSlopeGuard.correct_steps',side_effect=failed):
            with self.assertRaises(PublicationEntranceSlopeError):
                builder.correct_publication(mesh,field,fixed,conditions,exaggeration=1.,limit=.1,max_cells=100)


if __name__=='__main__':unittest.main()
