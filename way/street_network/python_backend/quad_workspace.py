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
        np.sum(self.corners, axis=-1, out=self.numerator)
        shape = self.numerator / self.jac
        self.value = float(np.sum(shape) if self.power == 1 else np.sum(shape ** 2) / 4.)
        return self.value


class QuadGradient:
    def __init__(self, faces, vertices):
        self.targets = np.stack([
            np.stack((faces[:, i], faces[:, (i + 1) % 4], faces[:, (i + 2) % 4]))
            for i in range(4)]).reshape(-1).astype(np.intp, copy=False)
        self.contributions = np.empty((4, 3, len(faces), 2))
        self.right_a = np.empty((len(faces), 4, 2))
        self.right_b = np.empty_like(self.right_a)
        self.gradient = np.zeros((vertices, 2))

    def evaluate(self, energy):
        a, b = energy.edge, energy.next_edge
        b[:, :-1] = a[:, 1:]
        b[:, -1] = a[:, 0]
        self.right_a[..., 0], self.right_a[..., 1] = a[..., 1], -a[..., 0]
        self.right_b[..., 0], self.right_b[..., 1] = b[..., 1], -b[..., 0]
        denominator = energy.jac[..., None]
        numerator = energy.numerator[..., None]
        da = 2. * a / denominator - numerator * self.right_b / denominator ** 2
        db = 2. * b / denominator + numerator * self.right_a / denominator ** 2
        if energy.power == 2:
            factor = numerator / (2. * denominator)
            da, db = da * factor, db * factor
        self.contributions[:, 0] = -da.transpose(1, 0, 2)
        self.contributions[:, 1] = (da - db).transpose(1, 0, 2)
        self.contributions[:, 2] = db.transpose(1, 0, 2)
        # Both coordinate streams retain the original corner/contribution/face
        # accumulation order. minlength also clears unreferenced vertices.
        values = self.contributions.reshape(-1, 2)
        self.gradient[:, 0] = np.bincount(self.targets, weights=values[:, 0], minlength=len(self.gradient))
        self.gradient[:, 1] = np.bincount(self.targets, weights=values[:, 1], minlength=len(self.gradient))
        return self.gradient
