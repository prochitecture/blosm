"""Call-local workspaces for the fixed-rim quad shape descent.

These buffers change allocation and scatter overhead, not the shape energy or
line search. Keep corner/contribution/face accumulation order: regrouping the
floating-point sums by vertex can change subsequent accepted steps.
"""
import numpy as np
from .mesh import cross


class QuadEnergy:
    def __init__(self, faces, power):
        self.faces, self.power = faces, power
        self.corners = np.empty((len(faces), 4, 2))
        self.edge = np.empty_like(self.corners)
        self.next_edge = np.empty_like(self.corners)
        self.jac = np.empty((len(faces), 4))
        self.numerator = np.empty_like(self.jac)
        self.value = float('inf')
        self.rejected_corner = None

    def evaluate(self, points):
        # A corner rejected on a previous trial is a cheap exact witness. Read
        # its current coordinates/connectivity: it is never a cached geometric
        # value or an approximate acceptance test. Only finite evaluations have
        # usable scratch buffers; an accepted proposal still evaluates all faces.
        if self.rejected_corner is not None:
            face, corner = self.rejected_corner
            ids = self.faces[face]
            a = points[ids[(corner + 1) % 4]] - points[ids[corner]]
            b = points[ids[(corner + 2) % 4]] - points[ids[(corner + 1) % 4]]
            if a[0] * b[1] - a[1] * b[0] <= 1e-10:
                self.value = float('inf')
                return self.value
        np.take(points, self.faces, axis=0, out=self.corners)
        np.subtract(self.corners[:, 1:], self.corners[:, :-1], out=self.edge[:, :-1])
        np.subtract(self.corners[:, 0], self.corners[:, -1], out=self.edge[:, -1])
        self.next_edge[:, :-1] = self.edge[:, 1:]
        self.next_edge[:, -1] = self.edge[:, 0]
        self.jac[:] = cross(self.edge, self.next_edge)
        if self.jac.min() <= 1e-10:
            self.rejected_corner = divmod(int(np.argmin(self.jac)), 4)
            self.value = float('inf')
            return self.value
        # Preserve the original squared-edge addition and reduction order.
        # next_edge is scratch after evaluation; QuadGradient restores it.
        np.square(self.edge, out=self.corners)
        np.square(self.next_edge, out=self.next_edge)
        np.add(self.corners, self.next_edge, out=self.corners)
        # Exactly two positive squared-coordinate terms: the same one addition
        # can use a ufunc stream instead of a separate two-item reduction at
        # every corner. The shape sum below keeps its original reduction order.
        np.add(self.corners[..., 0], self.corners[..., 1], out=self.numerator)
        shape = self.numerator / self.jac
        self.value = float(np.sum(shape) if self.power == 1 else np.sum(shape ** 2) / 4.)
        return self.value


class QuadGradient:
    def __init__(self, faces, vertices):
        self.targets = np.stack([
            np.stack((faces[:, i], faces[:, (i + 1) % 4], faces[:, (i + 2) % 4]))
            for i in range(4)]).reshape(-1).astype(np.intp, copy=False)
        self.contributions = np.empty((4, 3, len(faces), 2))
        # Long contiguous coordinate streams let NumPy process whole corners
        # per ufunc loop, rather than broadcasting a denominator every two
        # values. Only storage changes; the scatter stream remains unchanged.
        self.right_a = np.empty((2, len(faces), 4))
        self.right_b = np.empty_like(self.right_a)
        self.a = np.empty_like(self.right_a)
        self.b = np.empty_like(self.right_a)
        self.da = np.empty_like(self.right_a)
        self.db = np.empty_like(self.right_a)
        self.term = np.empty_like(self.right_a)
        self.square = np.empty((1, len(faces), 4))
        self.factor = np.empty_like(self.square)
        self.gradient = np.zeros((vertices, 2))

    def evaluate(self, energy):
        energy.next_edge[:, :-1] = energy.edge[:, 1:]
        energy.next_edge[:, -1] = energy.edge[:, 0]
        a, b = self.a, self.b
        a[:] = energy.edge.transpose(2, 0, 1)
        b[:] = energy.next_edge.transpose(2, 0, 1)
        self.right_a[0], self.right_a[1] = a[1], -a[0]
        self.right_b[0], self.right_b[1] = b[1], -b[0]
        denominator = energy.jac[None]
        numerator = energy.numerator[None]
        # Preserve the original multiply/divide/add order, including both
        # divisions. Algebraically factoring these expressions changes the
        # floating-point gradient and can change later line-search decisions.
        da, db, term = self.da, self.db, self.term
        np.square(denominator, out=self.square)
        np.multiply(2., a, out=da)
        np.divide(da, denominator, out=da)
        np.multiply(numerator, self.right_b, out=term)
        np.divide(term, self.square, out=term)
        np.subtract(da, term, out=da)
        np.multiply(2., b, out=db)
        np.divide(db, denominator, out=db)
        np.multiply(numerator, self.right_a, out=term)
        np.divide(term, self.square, out=term)
        np.add(db, term, out=db)
        if energy.power == 2:
            np.multiply(2., denominator, out=self.factor)
            np.divide(numerator, self.factor, out=self.factor)
            np.multiply(da, self.factor, out=da)
            np.multiply(db, self.factor, out=db)
        np.negative(da.transpose(2, 1, 0), out=self.contributions[:, 0])
        np.subtract(da, db, out=term)
        self.contributions[:, 1] = term.transpose(2, 1, 0)
        self.contributions[:, 2] = db.transpose(2, 1, 0)
        # Both coordinate streams retain the original corner/contribution/face
        # accumulation order. minlength also clears unreferenced vertices.
        values = self.contributions.reshape(-1, 2)
        self.gradient[:, 0] = np.bincount(self.targets, weights=values[:, 0], minlength=len(self.gradient))
        self.gradient[:, 1] = np.bincount(self.targets, weights=values[:, 1], minlength=len(self.gradient))
        return self.gradient
