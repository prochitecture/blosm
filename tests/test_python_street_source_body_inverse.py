"""The entire original normal-root bracket remains rigorously monotone."""
from types import SimpleNamespace
import math
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.profiles import Crown
from tested_python_streets.source_body_references import (
    OriginalBodyReference,normal_root_distance_bound)


class AnalyticalArc:
    """Independent unit-speed reference, not a replacement native Path solver."""
    def __init__(self,*,angle=0.,shift=(0.,0.),reverse=False,width=5.):
        self.rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
        self.shift=np.asarray(shift);self.reverse=reverse
        self.length=4.;self.radius=4.;self.curvature_bound=.25;self.exact=False
        self.section=SimpleNamespace(id=1,width=width,offset=0.)
    def frame(self,station):
        at=np.asarray(station);value=self.length-at if self.reverse else at
        theta=value/self.radius
        center=np.stack((self.radius*np.sin(theta),self.radius*(1.-np.cos(theta))),axis=-1)
        tangent=np.stack((np.cos(theta),np.sin(theta)),axis=-1)*(-1. if self.reverse else 1.)
        normal=np.stack((tangent[...,1],-tangent[...,0]),axis=-1)
        return center@self.rotation.T+self.shift,tangent@self.rotation.T,normal@self.rotation.T
    def point(self,station,offset):
        center,_,normal=self.frame(station)
        return center+np.asarray(offset)[...,None]*normal


class OriginalBodyInverseTests(unittest.TestCase):
    def test_previously_positive_certificate_never_uses_the_retry_bound(self):
        path=AnalyticalArc()
        reference=OriginalBodyReference(path,Crown(.08,0.),0.,4.,np.linspace(0.,4.,21))
        xy=path.point(np.array((1.5,2.5)),np.array((-.1,.1)))
        with patch('tested_python_streets.source_body_references.normal_root_distance_bound',
                   side_effect=AssertionError('An old positive certificate must remain unchanged')):
            _,_,valid,proof=reference.inverse(xy)
        self.assertTrue(valid.all())
        self.assertNotIn('normal_root_distance_bound',proof)

    def test_chord_error_bound_contains_whole_bracket_and_proves_previously_coarse_case(self):
        path=AnalyticalArc();low=np.array((0.,1.));high=np.array((3.,4.))
        points=path.point(np.array((1.5,2.5)),np.full(2,2.))
        bound=normal_root_distance_bound(points,path,low,high)
        old=np.linalg.norm(points-path.frame(low)[0],axis=1)+(high-low)
        self.assertTrue(np.all(1.-path.curvature_bound*old<0.))
        self.assertTrue(np.all(1.-path.curvature_bound*bound>0.))
        for t in np.linspace(0.,1.,1001):
            actual=np.linalg.norm(points-path.frame(low+t*(high-low))[0],axis=1)
            self.assertTrue(np.all(actual<=bound))

    def test_original_bracket_root_is_accurate_under_rigid_transform_and_reversal(self):
        for angle,shift,reverse in ((0.,(0.,0.),False),(.713,(300.,-700.),False),(.713,(300.,-700.),True)):
            path=AnalyticalArc(angle=angle,shift=shift,reverse=reverse)
            at=np.array((1.5,2.5));q=np.full(2,-.8 if reverse else .8)
            reference=OriginalBodyReference(path,Crown(.08,0.),0.,4.,np.arange(5,dtype=float))
            xy=path.point(at,.5*path.section.width*q)
            station,lateral,valid,proof=reference.inverse(xy)
            self.assertTrue(valid.all());self.assertGreater(proof['minimum_normal_root_jacobian'],0.)
            np.testing.assert_allclose(station,at,atol=1e-11,rtol=0.)
            np.testing.assert_allclose(lateral,q,atol=1e-11,rtol=0.)
            self.assertLess(proof['maximum_inverse_residual'],1e-10)

    def test_finite_extrapolation_and_genuinely_folded_reference_are_still_rejected(self):
        path=AnalyticalArc();reference=OriginalBodyReference(path,Crown(.08,0.),0.,4.,np.arange(5,dtype=float))
        _,_,valid,_=reference.inverse(path.point(np.array((-.5,4.5)),np.zeros(2)))
        self.assertFalse(valid.any())
        folded=AnalyticalArc(width=9.)
        reference=OriginalBodyReference(folded,Crown(.08,0.),0.,4.,np.arange(5,dtype=float))
        with self.assertRaises(ValueError):reference.inverse(folded.point(np.array((2.,)),np.array((0.,))))


if __name__=='__main__':unittest.main()
