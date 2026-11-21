"""Finite curvature-matched root curves; existing continuations stay first.

These are geometric G2 endpoint jets, with the original endpoint speeds. They
do not prescribe the parameterization of the native bank outside the contact.
The added bells have zero position/first/second jets at their endpoints.
"""
import math
import numpy as np

from .continuations import WidthContinuation,ContinuationReachError,_evaluate
from .curb_envelope import CurbEnvelope
from .mesh import cross


def _number(value):
    return (not isinstance(value,(bool,np.bool_)) and isinstance(value,(int,float,np.integer,np.floating))
            and math.isfinite(value))


def _finish(curve,controls,minimum_radius):
    """The existing Bernstein certificate, with the actual polynomial degree."""
    n=len(controls)-1
    curve.controls=controls
    curve.first_derivative=n*np.diff(controls,axis=0)
    curve.second_derivative=n*(n-1)*np.diff(controls,n=2,axis=0)
    axis=controls[-1]-controls[0];length=float(np.linalg.norm(axis))
    if length<=1e-8:raise ValueError('Original root curve has no finite chord')
    axis/=length
    if np.min(curve.first_derivative@axis)<=1e-8:
        raise ContinuationReachError('Original root curve lacks positive complete chord progress')
    numerator=np.zeros(2*n-2)
    for i,velocity in enumerate(curve.first_derivative):
        for j,acceleration in enumerate(curve.second_derivative):
            numerator[i+j]+=(math.comb(n-1,i)*math.comb(n-2,j)
                /math.comb(2*n-3,i+j))*cross(velocity,acceleration)
    curve.curvature_bound,curve.curvature_intervals=curve._refined_curvature(numerator,minimum_radius)
    curve.maximum_speed=float(np.linalg.norm(curve.first_derivative,axis=1).max())
    curve.maximum_second=float(np.linalg.norm(curve.second_derivative,axis=1).max())
    nodes,weights=np.polynomial.legendre.leggauss(8)
    curve.length=float(.5*np.sum(weights*np.linalg.norm(
        _evaluate(curve.first_derivative,.5*(nodes+1.)),axis=1)))
    for array in (curve.controls,curve.first_derivative,curve.second_derivative):array.flags.writeable=False


class NativeCurvatureCore(WidthContinuation):
    """Quintic Hermite core with the unchanged endpoint positions and speeds."""
    def __init__(self,base,curvatures,minimum_radius):
        old=np.asarray(base.controls)
        if (old.shape!=(6,2) or not np.isfinite(old).all() or not isinstance(curvatures,(tuple,list,np.ndarray)) or len(curvatures)!=2
                or not all(_number(k) for k in curvatures)
                or not _number(minimum_radius) or minimum_radius<=0.):
            raise ValueError('Native root curvature needs two finite original endpoint jets')
        controls=old.copy()
        for endpoint,control,velocity in ((0,2,base.first_derivative[0]),(1,3,base.first_derivative[-1])):
            speed=float(np.linalg.norm(velocity))
            if speed<=1e-8:raise ValueError('Native root endpoint has no finite tangent speed')
            tangent=velocity/speed;left=np.array((-tangent[1],tangent[0]))
            acceleration=float(curvatures[endpoint])*speed**2*left
            controls[control]=(2.*old[1]-old[0] if endpoint==0 else 2.*old[-2]-old[-1])+acceleration/20.
        self.native_curvatures=tuple(float(k) for k in curvatures)
        # A rigorous convex-hull displacement bound relative to the old core;
        # this is separate from the added bell's physical displacement cap.
        self.base_displacement_bound=float(np.linalg.norm(controls-old,axis=1).max())
        _finish(self,controls,minimum_radius)


class ShiftedRootBell(WidthContinuation):
    """One normalized degree-nine Bernstein bell, peaked at 1/3 or 2/3."""
    def __init__(self,base,outward,distance,minimum_radius,*,index,maximum_distance):
        old=np.asarray(base.controls);normal=np.asarray(outward,float)
        if (type(index) is not int or index not in (3,6) or old.shape!=(6,2)
                or not np.isfinite(old).all() or normal.shape!=(2,) or not np.isfinite(normal).all()
                or abs(np.linalg.norm(normal)-1.)>1e-8 or not _number(distance) or distance<0.
                or not _number(maximum_distance) or maximum_distance<0. or distance>maximum_distance
                or not _number(minimum_radius) or minimum_radius<=0.):
            raise ValueError('Invalid bounded original root bell')
        controls=old.copy()
        for degree in range(6,10):
            elevated=np.empty((degree+1,2));elevated[0],elevated[-1]=controls[0],controls[-1]
            for i in range(1,degree):
                elevated[i]=(i/degree)*controls[i-1]+(1.-i/degree)*controls[i]
            controls=elevated
        peak=index/9.
        maximum=math.comb(9,index)*peak**index*(1.-peak)**(9-index)
        controls[index]+=(distance/maximum)*normal
        self.envelope_distance=float(distance);self.bell_index=index;self.bell_degree=9
        _finish(self,controls,minimum_radius)

    @staticmethod
    def bell(t,index):
        if type(index) is not int or index not in (3,6):raise ValueError('Original root bell needs an authored finite shape')
        peak=index/9.;t=np.asarray(t)
        return t**index*(1.-t)**(9-index)/(peak**index*(1.-peak)**(9-index))


def bounded_bell(base,outward,distance,minimum_radius,*,index,maximum_distance):
    """Leave the degree-six constructor and its successful arithmetic intact."""
    if (not _number(distance) or not _number(maximum_distance) or distance<0.
            or maximum_distance<0. or distance>maximum_distance):
        raise ValueError('Original root bell exceeds its finite displacement cap')
    if index is None:return CurbEnvelope(base,outward,distance,minimum_radius)
    return ShiftedRootBell(base,outward,distance,minimum_radius,index=index,maximum_distance=maximum_distance)
