"""Necessary incoming feasibility keeps immutable data and the original box."""
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.height_slabs import HeightSlabs
from tested_python_streets.source_height_chart_admission import entrance_interval_floor


class SourceHeightChartAdmissionTests(unittest.TestCase):
    def test_necessary_admission_uses_no_generic_interior_or_polynomial_field(self):
        from test_python_street_ordinary_source_prefix_atlas import straight_scope_fixture
        from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
        from tested_python_streets.partition import CompoundQuadChart
        from tested_python_streets.source_height_chart_admission import SourceHeightChartAdmission
        context,group=straight_scope_fixture()
        chart=CompoundQuadChart(context.patterns).build(group,context)
        builder=CompoundSurfaceBuilder(context.patterns,source_heights=True)
        with patch('tested_python_streets.surface_fields.BoundaryJetField',
                   side_effect=AssertionError('Necessary admission must not construct a fitted field')):
            report=SourceHeightChartAdmission(builder,context)(group,chart)
        self.assertFalse(report['solver_ran'])
        self.assertTrue(report['mandatory_final_surface_checks'])

    def test_immutable_entrance_violation_is_not_mistaken_for_body_freedom(self):
        rows=HeightSlabs(np.array(((0,1,2,3),)),np.array(((1.,-1.,0.,0.),)),np.array((.03,)))
        heights=np.array((.04,0.,.08,.08));immutable=np.ones(4,bool)
        self.assertGreater(entrance_interval_floor(rows,heights,immutable,0.,.08)[0],.03)
        immutable[1]=False
        self.assertEqual(entrance_interval_floor(rows,heights,immutable,0.,.08)[0],0.)

    def test_free_height_box_can_prove_impossibility_but_never_fix_a_free_value(self):
        rows=HeightSlabs(np.array(((0,1,2,3),)),np.array(((1.,-1.,-1.,0.),)),np.array((.03,)))
        heights=np.array((.2,0.,.04,0.));immutable=np.array((True,True,False,False))
        before=heights.copy()
        self.assertGreater(entrance_interval_floor(rows,heights,immutable,0.,.08)[0],.11)
        self.assertEqual(entrance_interval_floor(rows,heights,immutable,0.,.2)[0],0.)
        np.testing.assert_array_equal(heights,before)

    def test_outward_rounding_does_not_reject_the_configured_limit_itself(self):
        rows=HeightSlabs(np.array(((0,1,2,3),)),np.array(((1.,-1.,0.,0.),)),np.array((.03,)))
        heights=np.array((.03,0.,.08,.08));immutable=np.ones(4,bool)
        self.assertLessEqual(entrance_interval_floor(rows,heights,immutable,0.,.08)[0],.03)
        for mask in ((True,False,False,True),(False,False,True,True)):
            coefficients=np.array(((3.,-4.,2.,-1.),))
            rows=HeightSlabs(np.array(((0,1,2,3),)),coefficients,np.array((.03,)))
            lower=entrance_interval_floor(rows,heights,np.array(mask),0.,.08)[0]
            free=np.flatnonzero(~np.array(mask))
            for a in np.linspace(0.,.08,31):
                for b in np.linspace(0.,.08,31):
                    value=heights.copy();value[free]=(a,b)
                    self.assertLessEqual(lower,abs(float(coefficients[0]@value))+1e-15)


if __name__=='__main__':unittest.main()
