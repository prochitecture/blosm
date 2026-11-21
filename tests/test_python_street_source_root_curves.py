"""Portable geometric proofs for the opt-in native-curvature root family."""
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.continuations import WidthContinuation,_evaluate,ContinuationReachError
from tested_python_streets.curb_envelope import CurbEnvelope
from tested_python_streets.source_root_curves import NativeCurvatureCore,ShiftedRootBell,bounded_bell
from tested_python_streets.mesh import cross


def fixture(reverse=False,angle=0.,translation=(0.,0.)):
    first,last=np.array((0.,0.)),np.array((16.,1.))
    va,vb=np.array((1.,.12)),np.array((1.,-.04))
    curvatures=(-.08,-.005)
    if reverse:first,last,va,vb=last,first,-vb,-va;curvatures=tuple(-k for k in curvatures[::-1])
    rot=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    base=WidthContinuation(first@rot.T+translation,last@rot.T+translation,va@rot.T,vb@rot.T,2.,refine_curvature=True)
    return base,curvatures,np.array((0.,-1.))@rot.T


class SourceRootCurveTests(unittest.TestCase):
    def test_native_core_retains_positions_first_jets_and_geometric_curvature(self):
        base,targets,_=fixture();before=base.controls.tobytes();core=NativeCurvatureCore(base,targets,2.)
        np.testing.assert_array_equal(core.controls[[0,1,4,5]],base.controls[[0,1,4,5]])
        np.testing.assert_array_equal(core.first_derivative[[0,-1]],base.first_derivative[[0,-1]])
        for i,target in enumerate(targets):
            velocity=core.first_derivative[0 if i==0 else -1]
            acceleration=core.second_derivative[0 if i==0 else -1]
            self.assertAlmostEqual(float(cross(velocity,acceleration)/np.linalg.norm(velocity)**3),target,places=12)
        self.assertEqual(before,base.controls.tobytes())
        sample=np.linspace(0.,1.,4097)
        self.assertLessEqual(float(np.linalg.norm(core(sample)-base(sample),axis=1).max()),core.base_displacement_bound+1e-12)

    def test_bell_is_bounded_and_has_zero_endpoint_jets(self):
        base,targets,normal=fixture();core=NativeCurvatureCore(base,targets,2.)
        for index in (3,6):
            curve=ShiftedRootBell(core,normal,.6,2.,index=index,maximum_distance=1.75)
            t=np.linspace(0.,1.,4097)
            self.assertLessEqual(float(np.linalg.norm(curve(t)-core(t),axis=1).max()),.6+1e-11)
            np.testing.assert_allclose(curve(np.array((0.,1.))),core(np.array((0.,1.))),atol=1e-12,rtol=0.)
            np.testing.assert_allclose(curve.derivative(np.array((0.,1.))),core.derivative(np.array((0.,1.))),atol=1e-11,rtol=0.)
            np.testing.assert_allclose(_evaluate(curve.second_derivative,np.array((0.,1.))),
                _evaluate(core.second_derivative,np.array((0.,1.))),atol=2e-10,rtol=0.)
            self.assertAlmostEqual(float(ShiftedRootBell.bell(index/9.,index)),1.)
            self.assertLessEqual(curve.curvature_intervals,128)
            velocity=curve.derivative(t);acceleration=_evaluate(curve.second_derivative,t)
            curvature=np.abs(cross(velocity,acceleration))/np.linalg.norm(velocity,axis=1)**3
            self.assertLessEqual(float(curvature.max()),curve.curvature_bound+1e-12)
            self.assertLessEqual(curve.curvature_bound,.5*(1.+1e-8))

    def test_reverse_and_rigid_transforms_preserve_the_same_geometric_family(self):
        base,k,n=fixture();curve=ShiftedRootBell(NativeCurvatureCore(base,k,2.),n,.6,2.,index=6,maximum_distance=1.75)
        other,other_k,other_n=fixture(reverse=True)
        reversed_curve=ShiftedRootBell(NativeCurvatureCore(other,other_k,2.),other_n,.6,2.,index=3,maximum_distance=1.75)
        t=np.linspace(0.,1.,129)
        np.testing.assert_allclose(reversed_curve(t),curve(1.-t),atol=2e-11,rtol=0.)
        angle=.713;shift=np.array((1200.,-1800.));rot=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
        other,other_k,other_n=fixture(angle=angle,translation=shift)
        moved=ShiftedRootBell(NativeCurvatureCore(other,other_k,2.),other_n,.6,2.,index=6,maximum_distance=1.75)
        np.testing.assert_allclose(moved(t),curve(t)@rot.T+shift,atol=3e-11,rtol=0.)

    def test_degree_six_default_constructor_is_unchanged(self):
        base,_,n=fixture();normal=CurbEnvelope(base,n,.4,2.)
        explicit=bounded_bell(base,n,.4,2.,index=None,maximum_distance=1.75)
        self.assertIs(type(explicit),CurbEnvelope)
        np.testing.assert_array_equal(explicit.controls,normal.controls)
        np.testing.assert_array_equal(explicit.first_derivative,normal.first_derivative)
        np.testing.assert_array_equal(explicit.second_derivative,normal.second_derivative)
        self.assertEqual(explicit.curvature_bound,normal.curvature_bound)

    def test_invalid_caps_radii_shapes_and_nonfinite_jets_are_refused(self):
        base,k,n=fixture()
        for radius in (0.,-1.,float('nan'),True):
            with self.assertRaises(ValueError):NativeCurvatureCore(base,k,radius)
        with self.assertRaises(ValueError):NativeCurvatureCore(base,(float('inf'),0.),2.)
        for index in (0,4,9,True):
            with self.assertRaises(ValueError):ShiftedRootBell(base,n,.6,2.,index=index,maximum_distance=1.75)
        for distance in (-1.,1.750001,float('nan'),True):
            with self.assertRaises(ValueError):bounded_bell(base,n,distance,2.,index=6,maximum_distance=1.75)
        with self.assertRaises(ContinuationReachError):NativeCurvatureCore(base,(-2.,0.),2.)


if __name__=='__main__':unittest.main()
