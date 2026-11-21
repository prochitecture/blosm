"""A bounded C2 curb envelope for a completely owned finite dead-end pad."""
import math
import numpy as np
from .continuations import WidthContinuation, _evaluate
from .mesh import cross


class CurbEnvelope(WidthContinuation):
    """Add one zero-jet Bernstein bell to an existing quintic curb.

    Elevating its degree preserves the original curve. The added degree-six
    basis has zero position, first and second derivatives at both ends, so
    complete incoming profiles and their tangent planes remain authoritative.
    A positive speed projection and interval Bernstein bounds certify radius.
    """
    def __init__(self, base, outward, distance, minimum_radius):
        normal = np.asarray(outward, float)
        old = np.asarray(base.controls, float)
        if (normal.shape != (2,) or not np.isfinite(normal).all()
                or abs(np.linalg.norm(normal)-1.) > 1e-8
                or not math.isfinite(distance) or distance < 0.
                or not math.isfinite(minimum_radius) or minimum_radius <= 0.
                or old.shape != (6, 2) or not np.isfinite(old).all()
                or np.linalg.norm(old[-1]-old[0]) <= 1e-8):
            raise ValueError('Invalid finite curb envelope')
        controls = np.empty((7, 2))
        controls[0], controls[-1] = old[0], old[-1]
        for i in range(1, 6):
            controls[i] = (i/6.)*old[i-1]+(1.-i/6.)*old[i]
        controls[3] += 3.2*distance*normal  # 64 t³(1-t)³ has unit peak.
        self.controls = controls
        self.first_derivative = 6.*np.diff(controls, axis=0)
        self.second_derivative = 30.*np.diff(controls, n=2, axis=0)
        axis = controls[-1]-controls[0]
        axis /= np.linalg.norm(axis)
        if np.min(self.first_derivative @ axis) <= 1e-8:
            raise ValueError('Finite pad envelope turns back')
        coefficients = np.zeros(10)
        for i, velocity in enumerate(self.first_derivative):
            for j, acceleration in enumerate(self.second_derivative):
                coefficients[i+j] += (math.comb(5, i)*math.comb(4, j)
                    /math.comb(9, i+j))*cross(velocity, acceleration)
        self.curvature_bound, self.curvature_intervals = self._refined_curvature(coefficients, minimum_radius)
        self.maximum_speed = float(np.linalg.norm(self.first_derivative, axis=1).max())
        self.maximum_second = float(np.linalg.norm(self.second_derivative, axis=1).max())
        nodes, weights = np.polynomial.legendre.leggauss(8)
        self.length = float(.5*np.sum(weights*np.linalg.norm(_evaluate(self.first_derivative, .5*(nodes+1.)), axis=1)))
        self.envelope_distance = float(distance)
        for values in (self.controls, self.first_derivative, self.second_derivative):
            values.flags.writeable = False
