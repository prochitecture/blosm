"""Certified smooth curb continuation across a nearly straight width change."""
import math
import numpy as np
from .mesh import cross


def _evaluate(controls, t):
    t = np.asarray(t)
    degree = len(controls) - 1
    weights = np.stack([math.comb(degree, i) * t ** i * (1. - t) ** (degree - i)
                        for i in range(degree + 1)], axis=-1)
    return weights @ controls


class WidthContinuation:
    """Quintic Hermite curb, with endpoint tangents and a radius lower bound.

    A through-width transition is not a circular turning fillet. The curve has
    zero endpoint second derivatives. Bernstein numerator bounds and a positive
    speed projection certify its curvature over the complete parameter interval.
    """
    def __init__(self, first, last, start_tangent, end_tangent, minimum_radius, *, refine_curvature=False):
        first, last = np.asarray(first, dtype=float), np.asarray(last, dtype=float)
        a, b = np.asarray(start_tangent, dtype=float), np.asarray(end_tangent, dtype=float)
        if (not all(v.shape == (2,) and np.all(np.isfinite(v)) for v in (first, last, a, b))
                or not math.isfinite(minimum_radius) or minimum_radius <= 0.
                or min(np.linalg.norm(a), np.linalg.norm(b)) <= 1e-12):
            raise ValueError('Invalid width-continuation geometry')
        a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
        axis = a + b
        if np.linalg.norm(axis) <= 1e-8:
            raise ValueError('Width-continuation endpoint tangents turn back')
        axis /= np.linalg.norm(axis)
        length = float((last - first) @ axis)
        if length <= 1e-8:
            raise ValueError('Width continuation has no positive longitudinal reach')
        va, vb = length * a, length * b
        self.controls = np.array((first, first + va / 5., first + 2. * va / 5.,
                                  last - 2. * vb / 5., last - vb / 5., last))
        self.first_derivative = 5. * np.diff(self.controls, axis=0)
        self.second_derivative = 20. * np.diff(self.controls, n=2, axis=0)
        lower_speed = float(np.min(self.first_derivative @ axis))
        if lower_speed <= 1e-8:
            raise ValueError('Width continuation needs more longitudinal reach to avoid turning back')
        coefficients = np.zeros(8)
        for i, velocity in enumerate(self.first_derivative):
            for j, acceleration in enumerate(self.second_derivative):
                coefficients[i + j] += (math.comb(4, i) * math.comb(3, j)
                    / math.comb(7, i + j)) * cross(velocity, acceleration)
        self.curvature_bound = float(np.max(abs(coefficients)) / lower_speed ** 3)
        self.curvature_intervals = 1
        if self.curvature_bound > (1. + 1e-8) / minimum_radius:
            if not refine_curvature:
                raise ValueError('Width continuation needs more longitudinal reach for the required minimum radius')
            self.curvature_bound, self.curvature_intervals = self._refined_curvature(coefficients, minimum_radius)
        self.maximum_speed = float(np.max(np.linalg.norm(self.first_derivative, axis=1)))
        self.maximum_second = float(np.max(np.linalg.norm(self.second_derivative, axis=1)))
        nodes, weights = np.polynomial.legendre.leggauss(8)
        self.length = float(.5 * np.sum(weights * np.linalg.norm(_evaluate(self.first_derivative, .5 * (nodes + 1.)), axis=1)))
        for values in (self.controls, self.first_derivative, self.second_derivative):
            values.flags.writeable = False

    @staticmethod
    def _halves(values):
        """de Casteljau restriction at one half, retaining derivative units."""
        work = values
        left, right = [work[0]], [work[-1]]
        while len(work) > 1:
            work = .5 * (work[:-1] + work[1:])
            left.append(work[0])
            right.append(work[-1])
        return np.asarray(left), np.asarray(right[::-1])

    def _refined_curvature(self, numerator, minimum_radius):
        """Certify the same polynomial on at most 128 dyadic intervals.

        A global speed projection can underestimate speed near a sharp bend.
        On each restricted interval its midpoint velocity supplies a local
        axis. Bernstein hulls bound both that projection and v cross a over
        the entire interval, not just samples. The original global positive
        projection still proves forward progress; no tangent/control changes
        or smaller radius are permitted by this refinement.
        """
        eps = np.finfo(float).eps
        speed_error = 128. * eps * max(1., float(abs(self.first_derivative).max()))
        numerator_error = 128. * eps * max(1., float(abs(numerator).max()))
        pending = [(self.first_derivative, numerator, 0)]
        maximum, intervals = 0., 0
        while pending:
            velocity, coefficients, depth = pending.pop()
            axis = _evaluate(velocity, np.asarray(.5))
            axis /= np.linalg.norm(axis)
            lower = float(np.min(velocity @ axis)) - speed_error
            upper = float(np.max(abs(coefficients))) + numerator_error
            bound = upper / lower ** 3 if lower > 0. else math.inf
            if bound <= (1. + 1e-8) / minimum_radius:
                maximum = max(maximum, bound)
                intervals += 1
                continue
            if depth == 7:
                raise ValueError('Width continuation needs more longitudinal reach for the required minimum radius')
            a, b = self._halves(velocity)
            c, d = self._halves(coefficients)
            pending.extend(((a, c, depth + 1), (b, d, depth + 1)))
        return maximum, intervals

    def __call__(self, t):
        return _evaluate(self.controls, t)

    def derivative(self, t):
        return _evaluate(self.first_derivative, t)

    def intervals(self, spacing, tolerance):
        # Linear interpolation error on a parameter interval is <= M dt² / 8.
        # The speed bound also controls the physical length of each mesh edge.
        if not all(math.isfinite(v) and v > 0. for v in (spacing, tolerance)):
            raise ValueError('Width-continuation sampling needs positive finite limits')
        return max(1, int(math.ceil(self.maximum_speed / spacing)),
                   int(math.ceil(math.sqrt(self.maximum_second / (8. * tolerance)))))
