"""Experimental direct edge-ratio/angle feasibility, with a fixed physical rim."""
import numpy as np
from tested_python_streets.mesh import Mesh, cross


class TargetedQuality:
    def __init__(self, *, iterations=2500, ratio=7.75, sine=.05):
        self.iterations, self.ratio, self.sine = iterations, ratio, sine

    def improve(self, mesh, progress=None):
        xy, f = mesh.xyz[:, :2].copy(), mesh.quads
        edges = np.sort(np.stack((f, np.roll(f, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        fixed = np.zeros(len(xy), bool);fixed[unique[counts == 1].ravel()] = True
        src, dst = np.r_[unique[:, 0], unique[:, 1]], np.r_[unique[:, 1], unique[:, 0]]
        degree = np.bincount(dst, minlength=len(xy))

        def evaluate(p, gradient=False):
            edge = np.roll(p[f], -1, axis=1) - p[f]
            lengths = np.linalg.norm(edge, axis=-1)
            next_edge = np.roll(edge, -1, axis=1)
            jac = cross(edge, next_edge)
            if jac.min() <= 1e-10 or lengths.min() <= 1e-12:
                return np.inf, None, None
            low, high = lengths.argmin(axis=1), lengths.argmax(axis=1)
            rows = np.arange(len(f))
            a, b = lengths[rows, high], lengths[rows, low]
            ratios = a / b
            error = np.maximum(0., ratios - self.ratio)
            sine = jac / (lengths * np.roll(lengths, -1, axis=1))
            angle_error = np.maximum(0., self.sine - sine)
            value = float(np.sum(error ** 2) + 64. * np.sum(angle_error ** 2))
            report = dict(maximum_ratio=float(ratios.max()), minimum_sine=float(sine.min()),
                          above_ratio_8=int(np.sum(ratios > 8.)), value=value)
            if not gradient:
                return value, report, None
            grad = np.zeros_like(p)
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
                np.add.at(grad, f[:, i], -ge[:, i])
                np.add.at(grad, f[:, (i + 1) % 4], ge[:, i])
            return value, report, grad

        initial, report, grad = evaluate(xy, True)
        history = []
        best, best_score = xy.copy(), np.inf
        for it in range(self.iterations):
            value, report, grad = evaluate(xy, True)
            if it % 100 == 0 or it == self.iterations - 1:
                record = dict(iteration=it, **report)
                history.append(record)
                if progress:progress(record)
            score = max(report['maximum_ratio'] / 8., .03 / report['minimum_sine'])
            if score < best_score:
                best, best_score = xy.copy(), score
            if report['maximum_ratio'] <= 7.9 and report['minimum_sine'] >= .04:
                break
            length2 = np.sum((xy[src] - xy[dst]) ** 2, axis=1)
            mean2 = np.bincount(dst, weights=length2, minlength=len(xy)) / degree
            move = -grad * mean2[:, None] / degree[:, None]
            move[fixed] = 0.
            norm = np.linalg.norm(move, axis=1)
            cap = .1 * np.sqrt(mean2)
            move *= np.minimum(1., np.divide(cap, norm, out=np.ones_like(norm), where=norm > 0.))[:, None]
            for line in range(24):
                candidate = xy + 2. ** -line * move
                trial, _, _ = evaluate(candidate)
                if trial < value * (1. - 1e-12):
                    xy = candidate
                    break
            else:
                break
        result = Mesh(np.column_stack((best, mesh.xyz[:, 2])), f, mesh.material, dict(mesh.seams))
        return result, dict(iterations=it + 1, initial_loss=initial, final=evaluate(best)[1],
                            history=history, exact_fixed_rim=bool(np.array_equal(best[fixed], mesh.xyz[fixed, :2])))
