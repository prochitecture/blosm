"""Bounded planar quad correction aimed at the actual publication limits.

Physical rim vertices stay exact. Positive corner Jacobians are mandatory on
every accepted step. Feasibility is measured, never inferred from proxy energy.
"""
import math
import numpy as np
from .mesh import Mesh, cross


class QuadShapeCorrection:
    def __init__(self, *, max_iterations=2500, max_cells=50000, maximum_ratio=8., target_ratio=7.75,
                 minimum_sine=.03, target_sine=.05):
        if (type(max_iterations) is not int or max_iterations < 1 or type(max_cells) is not int or max_cells < 1
                or not all(math.isfinite(v) for v in (maximum_ratio, target_ratio, minimum_sine, target_sine))
                or not 1. < target_ratio < maximum_ratio or not 0. < minimum_sine < target_sine < 1.):
            raise ValueError('Invalid quad-shape correction bounds')
        self.max_iterations, self.max_cells = max_iterations, max_cells
        self.maximum_ratio, self.target_ratio = maximum_ratio, target_ratio
        self.minimum_sine, self.target_sine = minimum_sine, target_sine

    def evaluate(self, xy, faces, *, gradient=False):
        edge = np.roll(xy[faces], -1, axis=1) - xy[faces]
        lengths = np.linalg.norm(edge, axis=-1)
        next_edge = np.roll(edge, -1, axis=1)
        jac = cross(edge, next_edge)
        if jac.min() <= 1e-10 or lengths.min() <= 1e-12:
            return math.inf, None, None
        low, high = lengths.argmin(axis=1), lengths.argmax(axis=1)
        rows = np.arange(len(faces))
        a, b = lengths[rows, high], lengths[rows, low]
        ratios = a / b
        error = np.maximum(0., ratios - self.target_ratio)
        sine = jac / (lengths * np.roll(lengths, -1, axis=1))
        angle_error = np.maximum(0., self.target_sine - sine)
        value = float(np.sum(error ** 2) + 64. * np.sum(angle_error ** 2))
        report = dict(maximum_edge_ratio=float(ratios.max()), minimum_corner_sine=float(sine.min()),
                      cells_above_ratio_limit=int(np.sum(ratios > self.maximum_ratio)), loss=value)
        if not gradient:
            return value, report, None
        result = np.zeros_like(xy)
        ge = np.zeros_like(edge)
        ge[rows, high] += (2. * error / (a * b))[:, None] * edge[rows, high]
        ge[rows, low] -= (2. * error * a / b ** 3)[:, None] * edge[rows, low]
        right = lambda v: np.stack((v[..., 1], -v[..., 0]), axis=-1)
        factor = -128. * angle_error
        denominator = lengths * np.roll(lengths, -1, axis=1)
        ga = (right(next_edge) / denominator[..., None]
              - sine[..., None] * edge / lengths[..., None] ** 2) * factor[..., None]
        gb = (-right(edge) / denominator[..., None]
              - sine[..., None] * next_edge / np.roll(lengths, -1, axis=1)[..., None] ** 2) * factor[..., None]
        ge += ga + np.roll(gb, 1, axis=1)
        for i in range(4):
            np.add.at(result, faces[:, i], -ge[:, i])
            np.add.at(result, faces[:, (i + 1) % 4], ge[:, i])
        return value, report, result

    def qualified(self, report):
        return (report is not None and report['maximum_edge_ratio'] <= self.maximum_ratio
                and report['minimum_corner_sine'] >= self.minimum_sine)

    def improve(self, mesh, *, fixed_vertices=None):
        mesh.validate()
        if len(mesh.quads) > self.max_cells:
            raise ValueError('Quad-shape correction exceeds its cell budget')
        if np.ptp(mesh.xyz[:, 2]) > 1e-10:
            raise ValueError('Quad-shape correction needs a planar input before height construction')
        extra = np.asarray(() if fixed_vertices is None else fixed_vertices)
        if (extra.ndim != 1 or (extra.size and (extra.dtype.kind not in 'iu'
                or np.any(extra < 0) or np.any(extra >= len(mesh.xyz))))):
            raise ValueError('Quad-shape correction needs valid fixed guide vertex indices')
        extra = np.unique(extra.astype(np.int32))
        xy, faces = mesh.xyz[:, :2].copy(), mesh.quads
        initial_loss, initial, _ = self.evaluate(xy, faces)
        if self.qualified(initial):
            return mesh, dict(qualified=True, iterations=0, accepted_steps=0, final=initial,
                              initial=initial, exact_fixed_rim=True, fixed_guides=len(extra),
                              stop='already_qualified')
        edges = np.sort(np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        if counts.max() > 2:
            raise ValueError('Quad-shape correction requires manifold edge incidence')
        fixed = np.zeros(len(xy), bool)
        fixed[unique[counts == 1].ravel()] = True
        fixed[extra] = True
        src, dst = np.r_[unique[:, 0], unique[:, 1]], np.r_[unique[:, 1], unique[:, 0]]
        degree = np.bincount(dst, minlength=len(xy))
        if np.any(degree == 0):
            raise ValueError('Quad-shape correction needs incident faces at every vertex')
        best, best_score, accepted = xy.copy(), math.inf, 0
        history, stop, evaluations = [], 'iteration_budget', 0
        def score(report):
            return max(report['maximum_edge_ratio'] / self.maximum_ratio,
                       self.minimum_sine / report['minimum_corner_sine'])
        for iteration in range(self.max_iterations):
            value, report, gradient = self.evaluate(xy, faces, gradient=True)
            if score(report) < best_score:
                best, best_score = xy.copy(), score(report)
            if iteration % 100 == 0:
                history.append(dict(iteration=iteration, **report))
            # A margin avoids publishing a success at the floating-point limit.
            if (report['maximum_edge_ratio'] <= .5 * (self.maximum_ratio + self.target_ratio)
                    and report['minimum_corner_sine'] >= .5 * (self.minimum_sine + self.target_sine)):
                stop = 'qualified'
                break
            length2 = np.sum((xy[src] - xy[dst]) ** 2, axis=1)
            mean2 = np.bincount(dst, weights=length2, minlength=len(xy)) / degree
            move = -gradient * mean2[:, None] / degree[:, None]
            move[fixed] = 0.
            norm = np.linalg.norm(move, axis=1)
            cap = .1 * np.sqrt(mean2)
            move *= np.minimum(1., np.divide(cap, norm, out=np.ones_like(norm), where=norm > 0.))[:, None]
            for line in range(24):
                candidate = xy + 2. ** -line * move
                trial, _, _ = self.evaluate(candidate, faces)
                evaluations += 1
                if trial < value * (1. - 1e-12):
                    xy = candidate
                    accepted += 1
                    break
            else:
                stop = 'line_search_stalled'
                break
        _, last, _ = self.evaluate(xy, faces)
        if score(last) < best_score:
            best = xy.copy()
        _, final, _ = self.evaluate(best, faces)
        result = Mesh(np.column_stack((best, mesh.xyz[:, 2])), faces, mesh.material, dict(mesh.seams))
        result.validate()
        return result, dict(qualified=self.qualified(final), iterations=iteration + 1,
            accepted_steps=accepted, line_search_evaluations=evaluations, initial=initial, initial_loss=initial_loss,
            final=final, history=history, exact_fixed_rim=bool(np.array_equal(best[fixed], mesh.xyz[fixed, :2])),
            fixed_guides=len(extra),
            stop=stop, construction='fixed-rim edge-ratio and angle violation correction')
