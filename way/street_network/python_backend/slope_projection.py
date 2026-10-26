"""Bounded NumPy projections of fixed-chart heights onto slope constraints.

XY geometry, quad connectivity and all boundary heights stay fixed. The four
corner gradients bound the whole bilinear quad. Dykstra corrections combine
small local Euclidean projections with the common-datum/crown height box.
"""
import math
import numpy as np
from .mesh import Mesh, cross


def corner_matrices(mesh):
    xy = mesh.xyz[mesh.quads, :2]
    xy = xy - xy[:, :1]
    matrices = []
    for u, v in ((0., 0.), (1., 0.), (1., 1.), (0., 1.)):
        du = np.array((v - 1., 1. - v, v, -v))
        dv = np.array((u - 1., -u, u, 1. - u))
        a, b = np.einsum('fij,i->fj', xy, du), np.einsum('fij,i->fj', xy, dv)
        determinant = cross(a, b)
        if np.min(determinant) <= 1e-12:
            raise ValueError('Slope projection needs positive quad Jacobians')
        gx = (b[:, 1, None] * du - a[:, 1, None] * dv) / determinant[:, None]
        gy = (-b[:, 0, None] * du + a[:, 0, None] * dv) / determinant[:, None]
        matrices.append(np.stack((gx, gy), axis=1))
    return np.stack(matrices, axis=1)


def face_colors(faces, fixed, maximum):
    used = [0] * len(fixed)
    groups = []
    for index, face in enumerate(faces):
        free = [int(v) for v in face if not fixed[v]]
        mask = 0
        for vertex in free:
            mask |= used[vertex]
        bit = (mask + 1) & ~mask
        color = bit.bit_length() - 1
        if color >= maximum:
            raise ValueError('Slope projection exceeds its face-color budget')
        while len(groups) <= color:
            groups.append([])
        groups[color].append(index)
        for vertex in free:
            used[vertex] |= bit
    return tuple(np.asarray(group, dtype=np.int32) for group in groups)


class SlopeProjection:
    def __init__(self, *, max_iterations=1000, max_cells=50000,
                 memory_bytes=96 * 1024 * 1024, max_colors=64,
                 height_tolerance=1e-10, slope_tolerance=1e-9):
        if (type(max_iterations) is not int or max_iterations < 0
                or type(max_cells) is not int or max_cells < 1
                or type(memory_bytes) is not int or memory_bytes < 4096
                or type(max_colors) is not int or max_colors < 1
                or not all(math.isfinite(v) and v > 0. for v in (height_tolerance, slope_tolerance))):
            raise ValueError('Invalid slope-projection budgets')
        self.max_iterations, self.max_cells = max_iterations, max_cells
        self.memory_bytes, self.max_colors = memory_bytes, max_colors
        self.height_tolerance, self.slope_tolerance = height_tolerance, slope_tolerance

    @staticmethod
    def project(x, coefficient, free_coefficient, eigenvalues, eigenvectors, limit):
        """Project independent four-height rows onto one gradient-norm ball."""
        gradient = np.einsum('nij,nj->ni', coefficient, x)
        selected = np.flatnonzero(np.linalg.norm(gradient, axis=1) > limit + max(1e-14, limit * 1e-12))
        result = x.copy()
        if not len(selected):
            return result
        a = free_coefficient[selected]
        values, vectors = eigenvalues[selected], eigenvectors[selected]
        g = np.einsum('nji,nj->ni', vectors, gradient[selected])
        largest = values[:, 1]
        rank_zero = largest <= 1e-20
        if np.any(rank_zero):
            raise ValueError('Fixed boundary heights make the quad slope limit infeasible')
        rank_one = values[:, 0] <= largest * 1e-12
        delta = np.zeros_like(g)
        if np.any(rank_one):
            rows = np.flatnonzero(rank_one)
            normal2 = g[rows, 0] ** 2
            if np.any(normal2 > limit ** 2 + 1e-12):
                raise ValueError('Fixed boundary gradient makes the slope limit infeasible')
            radius = np.sqrt(np.maximum(0., limit ** 2 - normal2))
            change = np.clip(g[rows, 1], -radius, radius) - g[rows, 1]
            delta[rows, 1] = change / largest[rows]
        if np.any(~rank_one):
            rows = np.flatnonzero(~rank_one)
            w, h = values[rows], g[rows]
            norm = np.linalg.norm(h, axis=1)
            lower = np.zeros(len(rows))
            upper = (norm / limit - 1.) / w[:, 0]
            mu = (norm / limit - 1.) / w[:, 1]
            for _ in range(40):
                denominator = 1. + mu[:, None] * w
                component = h / denominator
                residual = np.sum(component ** 2, axis=1) - limit ** 2
                done = abs(residual) <= max(1e-16, limit ** 2 * 1e-12)
                if np.all(done):
                    break
                lower = np.where((~done) & (residual > 0.), mu, lower)
                upper = np.where((~done) & (residual <= 0.), mu, upper)
                derivative = -2. * np.sum(w * h ** 2 / denominator ** 3, axis=1)
                trial = mu - residual / derivative
                inside = (trial > lower) & (trial < upper) & np.isfinite(trial)
                mu = np.where(done, mu, np.where(inside, trial, .5 * (lower + upper)))
            else:
                raise ValueError('Local slope projection exceeded its scalar-root budget')
            delta[rows] = -mu[:, None] * h / (1. + mu[:, None] * w)
        world = np.einsum('nij,nj->ni', vectors, delta)
        result[selected] += np.einsum('nij,ni->nj', a, world)
        return result

    def solve(self, mesh, fixed, *, datum, ceiling, maximum_slope, minimum_heights=None,
              feasibility_slope=None):
        """Project to a slope ball, optionally stopping at a certified feasible mesh.

        By default the displacement must converge as well as the constraints.
        A caller needing feasibility rather than the nearest height vector may
        supply a larger acceptance slope and project towards a stricter inner
        target. Fixed vertices and the height box still hold at every return;
        the caller must independently certify its physical surface constraints.
        """
        mesh.validate()
        fixed = np.asarray(fixed)
        if fixed.dtype != np.bool_ or fixed.shape != (len(mesh.xyz),):
            raise ValueError('Slope projection needs a Boolean fixed-vertex mask')
        if (not all(map(math.isfinite, (datum, ceiling, maximum_slope)))
                or ceiling < datum or maximum_slope <= 0.):
            raise ValueError('Invalid slope-projection envelope')
        if feasibility_slope is not None and (not math.isfinite(feasibility_slope)
                or feasibility_slope < maximum_slope):
            raise ValueError('Feasibility slope must be finite and at least the projection target')
        lower = datum
        if minimum_heights is not None:
            lower = np.asarray(minimum_heights, dtype=float)
            if (lower.shape != (len(mesh.xyz),) or not np.all(np.isfinite(lower))
                    or np.any(lower < datum) or np.any(lower > ceiling)):
                raise ValueError('Invalid per-vertex minimum heights')
            if np.any(mesh.xyz[fixed, 2] < lower[fixed] - 1e-10):
                raise ValueError('Fixed boundary heights violate minimum heights')
            lower = lower.copy()
            lower[fixed] = mesh.xyz[fixed, 2]
        estimated = len(mesh.quads) * 2048 + len(mesh.xyz) * 96
        if len(mesh.quads) > self.max_cells or estimated > self.memory_bytes:
            raise ValueError('Slope projection exceeds its mesh/memory budget')
        initial = mesh.xyz[:, 2]
        if np.any(initial[fixed] < datum - 1e-10) or np.any(initial[fixed] > ceiling + 1e-10):
            raise ValueError('Fixed boundary heights exceed the crown envelope')
        coefficient = corner_matrices(mesh)
        faces = mesh.quads
        def gradients(z):
            return np.linalg.norm(np.einsum('fcij,fj->fci', coefficient, z[faces]), axis=-1)
        before = float(gradients(initial).max())
        if before <= maximum_slope and np.all((initial >= lower) & (initial <= ceiling)):
            return mesh, dict(iterations=0, maximum=before, maximum_height_change=0., colors=0,
                              estimated_bytes=estimated, construction='fixed-chart height projection')
        free = ~fixed[faces]
        movable = coefficient * free[:, None, None, :]
        gram = np.einsum('fcij,fckj->fcik', movable, movable)
        values, vectors = np.linalg.eigh(gram)
        values = np.maximum(values, 0.)
        groups = face_colors(faces, fixed, self.max_colors)
        corrections = np.zeros((len(faces), 4, 4))
        box_correction = np.zeros(len(initial))
        heights = initial.copy()
        aim = maximum_slope
        maximum, movement = before, float('inf')
        for iteration in range(self.max_iterations):
            previous = heights.copy()
            for group in groups:
                ids = faces[group]
                for corner in range(4):
                    x = heights[ids] + corrections[group, corner]
                    projected = self.project(x, coefficient[group, corner], movable[group, corner],
                        values[group, corner], vectors[group, corner], aim)
                    projected[~free[group]] = initial[ids][~free[group]]
                    corrections[group, corner] = (x - projected) * free[group]
                    heights[ids] = projected
            x = heights + box_correction
            heights = np.clip(x, lower, ceiling)
            heights[fixed] = initial[fixed]
            box_correction = x - heights
            box_correction[fixed] = 0.
            maximum = float(gradients(heights).max())
            movement = float(np.max(abs(heights - previous)))
            converged = (maximum <= maximum_slope + self.slope_tolerance and movement <= self.height_tolerance
                         if feasibility_slope is None else maximum <= feasibility_slope)
            if converged:
                xyz = mesh.xyz.copy()
                xyz[:, 2] = heights
                return Mesh(xyz, faces, mesh.material, dict(mesh.seams)), dict(iterations=iteration + 1,
                    initial_maximum=before, maximum=maximum,
                    maximum_height_change=float(np.max(abs(heights - initial))),
                    squared_displacement=float(np.sum((heights - initial) ** 2)),
                    colors=len(groups), estimated_bytes=estimated,
                    **({} if feasibility_slope is None else dict(convergence='feasible_surface')),
                    construction='fixed-chart height projection')
        raise ValueError(f'Height projection did not converge within its iteration budget: slope {maximum:.6f}, step {movement:.3g}')
