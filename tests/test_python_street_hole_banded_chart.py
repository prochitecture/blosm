"""Computational grass bands preserve exact rims and storage contracts."""
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.hole_banded_chart import build_banded_hole, audit
from tested_python_streets.mesh import Mesh
from tested_python_streets.model import Key
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.validation import validate_overlaps
from tested_python_streets.partition import CompoundQuadChart, GrassChartGradingError


def boundary():
    first=np.linspace(0.,80.,45)
    last=np.linspace(80.,0.,91)
    right=np.c_[2.+4.*np.sin(np.pi*first/80.),first]
    left=np.c_[-2.+4.*np.sin(np.pi*last/80.),last]
    top=np.linspace(0.,np.pi,17)[1:]
    bottom=np.linspace(np.pi,2*np.pi,17)[1:-1]
    return np.vstack((right,np.c_[2.*np.cos(top),80.+2.*np.sin(top)],left[1:],
                      np.c_[2.*np.cos(bottom),2.*np.sin(bottom)]))


class HoleBandedChartTests(unittest.TestCase):
    def test_curved_unequal_banks_keep_all_vertices_and_form_one_stored_quad_disk(self):
        original=boundary()
        for angle,offset in ((0.,(0.,0.)),(.713,(1200.,-1800.))):
            rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
            rim=original@rotation.T+offset
            output=build_banded_hole(PatternLibrary(),rim,.7,10000)
            report=audit(output.mesh,rim)
            self.assertLessEqual(report['stored_quality']['maximum_edge_ratio'],8.)
            self.assertTrue(report['exact_physical_polygon'])
            self.assertEqual(report['protected_vertices'],len(rim))
            self.assertEqual(output.mesh.material,'grass')
            validate_overlaps({Key('grass',0):output.mesh})
        np.testing.assert_array_equal(original,boundary())

    def test_grass_acceptance_keeps_small_fixed_tip_angles_and_rejects_lost_rim(self):
        rim=np.array(((0.,0.),(1.,-.005),(2.,0.),(1.,.005)))+(1200.,-1800.)
        mesh=Mesh(np.c_[rim,np.zeros(4)],np.array(((0,1,2,3),)),material='grass',seams={0:np.arange(4)})
        report=audit(mesh,rim)
        self.assertLess(report['stored_quality']['minimum_corner_sine'],.03)
        self.assertEqual(report['grass_quality_contract'],'positive corners and edge ratio <=8')
        changed=rim.copy();changed[1,0]+=.001
        with self.assertRaisesRegex(ValueError,'original physical rim vertex'):
            audit(mesh,changed)

    def test_invalid_topology_and_finite_allocations_reject(self):
        rim=boundary()
        for spacing,cells in ((float('nan'),50000),(.7,3),(.7,50001)):
            with self.assertRaises(ValueError):build_banded_hole(PatternLibrary(),rim,spacing,cells)
        with self.assertRaises(ValueError):build_banded_hole(PatternLibrary(),rim[::-1],.7,10000)
        crossed=np.array(((0.,0.),(2.,2.),(0.,2.),(2.,0.)))
        with self.assertRaises(ValueError):build_banded_hole(PatternLibrary(),crossed,.7,10000)
        with self.assertRaises(ValueError):build_banded_hole(PatternLibrary(),rim,.7,4)

    def test_later_old_physical_family_keeps_precedence_over_new_grass_recovery(self):
        from test_python_street_source_cycles import triangle,prepare
        context,group=prepare(triangle(2.))
        mesh=Mesh(np.array(((0.,0.,0.),(2.,0.,0.),(2.,1.,0.),(0.,1.,0.))),
                  np.array(((0,1,2,3),)))
        accepted=SimpleNamespace(mesh=mesh,metadata={},holes={})
        calls=[]
        def old_chart(chart,*args,**kwargs):
            yield
            calls.append(chart.recover_grass)
            if len(calls)==1:
                policy=chart.partitioner.policy
                raise GrassChartGradingError('old first grass family failed',family=(
                    policy.source_corridors,policy.hole_bridges,policy.existing_vertices))
            return accepted
        with patch.object(CompoundQuadChart,'_build_steps',old_chart), \
                patch('tested_python_streets.physical_hole_recovery.build_steps') as source_recovery, \
                patch('tested_python_streets.hole_banded_chart.retry_failed_family') as grass_recovery, \
                patch('tested_python_streets.paired_interval_atlas.build_steps') as paired_recovery:
            result=CompoundQuadChart(context.patterns).build(group,context)
        self.assertIs(result,accepted)
        self.assertEqual(calls,[False,False])
        source_recovery.assert_not_called()
        grass_recovery.assert_not_called()
        paired_recovery.assert_not_called()


if __name__=='__main__':unittest.main()
