"""Circular contacts on physical curbs, independently of profile miter frames.

A miter frame preserves an Exact road's offset boundary, but its transverse
vector need not have unit length or be perpendicular to that boundary. Radius
offsets therefore use the curb tangent's unit normal, not the profile vector.
"""
import numpy as np
from .corner_radii import CornerContactError
from .paths import right, unit


class CircularCurbContact:
    iterations = 32
    line_search = (1., .5, .25, .125, .0625)
    tolerance = 1e-8

    def __init__(self, first, second, radius, side):
        self.approaches = first, second
        self.distances = (first.offset - .5 * first.width,
                          second.offset + .5 * second.width)
        self.radii = -side * radius, side * radius

    def center(self, index, station):
        approach = self.approaches[index]
        point, tangent = approach.jet(float(station), self.distances[index])
        return point + self.radii[index] * right(unit(tangent))

    def jet(self, index, station):
        approach = self.approaches[index]
        step = 1e-3
        a, b = station - step, min(station + step, approach.length)
        point = self.center(index, station)
        return point, (self.center(index, b) - self.center(index, a)) / (b - a)

    def solve(self, initial, lower, upper):
        current = np.clip(np.asarray(initial, dtype=float), lower, upper)
        for _ in range(self.iterations):
            jets = [self.jet(i, station) for i, station in enumerate(current)]
            delta = jets[0][0] - jets[1][0]
            if np.linalg.norm(delta) < self.tolerance:
                contacts = tuple(a.point(float(s), d) for a, s, d in
                                 zip(self.approaches, current, self.distances))
                return current, jets[0][0], contacts
            matrix = np.column_stack((jets[0][1], -jets[1][1]))
            if abs(np.linalg.det(matrix)) < 1e-9:
                break
            step = np.linalg.solve(matrix, delta)
            for factor in self.line_search:
                candidate = np.clip(current - factor * step, lower, upper)
                error = self.center(0, candidate[0]) - self.center(1, candidate[1])
                if np.linalg.norm(error) < np.linalg.norm(delta):
                    current = candidate
                    break
            else:
                break
        raise CornerContactError('Requested-radius physical-curb contact did not converge')
