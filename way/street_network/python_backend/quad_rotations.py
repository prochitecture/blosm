"""Bounded rotations of shared edges between two existing quad cells.

Each alternative is two quads over the same six-vertex boundary. No triangular
faces, new vertices or boundary edits are involved. This is a local grading
repair; a badly parameterized region still needs another chart family.
"""
import math
import numpy as np
from .mesh import Mesh, cross


_PAIRS = np.array((((1, 2, 3, 4), (1, 4, 5, 0)),
                   ((2, 3, 4, 5), (2, 5, 0, 1))), dtype=np.int32)
_PAIRS.flags.writeable = False


def _metrics(xy, faces):
    points = xy[faces]
    edge = np.roll(points, -1, axis=-2) - points
    jacobian = cross(edge, np.roll(edge, -1, axis=-2))
    lengths = np.linalg.norm(edge, axis=-1)
    ratio = np.divide(lengths.max(axis=-1), lengths.min(axis=-1),
                      out=np.full(lengths.shape[:-1], np.inf), where=lengths.min(axis=-1) > 1e-12)
    numerator = lengths ** 2 + np.roll(lengths ** 2, -1, axis=-1)
    energy = np.divide(numerator, jacobian, out=np.full_like(jacobian, np.inf), where=jacobian > 1e-10)
    return jacobian.min(axis=-1), ratio, energy.sum(axis=-1)


class QuadEdgeRotations:
    def __init__(self, *, maximum_ratio=8., max_passes=4, max_candidates=1024, max_cells=50000):
        if (not math.isfinite(maximum_ratio) or maximum_ratio < 1.
                or type(max_passes) is not int or max_passes < 0
                or type(max_candidates) is not int or max_candidates < 1
                or type(max_cells) is not int or max_cells < 1):
            raise ValueError('Invalid quad-rotation budgets')
        self.maximum_ratio, self.max_passes = maximum_ratio, max_passes
        self.max_candidates, self.max_cells = max_candidates, max_cells

    def proposals(self, mesh):
        faces, xy = mesh.quads, mesh.xyz[:, :2]
        _, ratios, energy = _metrics(xy, faces)
        edge = np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2)
        unique, inverse, counts = np.unique(np.sort(edge, axis=1), axis=0,
                                            return_inverse=True, return_counts=True)
        if counts.max() > 2:
            raise ValueError('Quad rotations require manifold shared edges')
        internal = np.flatnonzero(counts == 2)
        order = np.argsort(inverse, kind='stable')
        offsets = np.r_[0, np.cumsum(counts)]
        pairs = order[np.column_stack((offsets[internal], offsets[internal] + 1))] // 4
        worst = ratios[pairs].max(axis=1)
        active = np.flatnonzero(worst > min(7.5, self.maximum_ratio))
        active = sorted(active, key=lambda i: (-round(float(worst[i]), 8), int(i)))[:self.max_candidates]
        rings, retained = [], []
        for at in active:
            first, last = map(int, pairs[at])
            a, b = map(int, unique[internal[at]])
            f, g = tuple(map(int, faces[first])), tuple(map(int, faces[last]))
            start = next(i for i in range(4) if {f[i], f[(i + 1) % 4]} == {a, b})
            f = f[start:] + f[:start]
            a, b = f[:2]
            start = g.index(b)
            g = g[start:] + g[:start]
            if g[1] != a:
                raise ValueError('Shared quad edges must have opposite orientations')
            ring = (a, g[2], g[3], b, f[2], f[3])
            if len(set(ring)) == 6:
                rings.append(ring)
                retained.append((first, last))
        if not rings:
            return ()
        pairs = np.asarray(retained)
        replacements = np.asarray(rings)[:, _PAIRS]
        jacobian, ratios_new, energy_new = _metrics(xy, replacements)
        before = ratios[pairs].max(axis=1)[:, None]
        after = ratios_new.max(axis=2)
        valid = ((jacobian.min(axis=2) > 1e-10) & (after < before * (1. - 1e-5))
                 & (energy_new.sum(axis=2) <= energy[pairs].sum(axis=1)[:, None] * (1. + 1e-12)))
        rows, variants = np.nonzero(valid)
        options = [(float(before[row, 0] - after[row, variant]), int(pairs[row, 0]),
                    int(pairs[row, 1]), int(variant), replacements[row, variant])
                   for row, variant in zip(rows, variants)]
        return sorted(options, key=lambda item: (-round(item[0], 8), item[1], item[2], item[3]))

    def improve(self, mesh, relaxation):
        """Return the best measured ratio, never a worse final iterate.

        Relaxation must retain physical rim vertices. The existing fixed-rim
        quality pass is supplied by the caller, keeping this topology utility
        independent of the chart constructor and its height field.
        """
        initial = mesh.validate()
        if len(mesh.quads) > self.max_cells:
            raise ValueError('Quad rotations exceed their cell budget')
        best, best_ratio, history = mesh, initial['maximum_edge_ratio'], []
        current = Mesh(mesh.xyz.copy(), mesh.quads.copy(), mesh.material, dict(mesh.seams))
        for _ in range(self.max_passes):
            if best_ratio <= self.maximum_ratio:
                break
            faces, used, accepted = current.quads.copy(), set(), 0
            for score, first, last, variant, replacement in self.proposals(current):
                if first in used or last in used:
                    continue
                faces[first], faces[last] = replacement
                used.update((first, last))
                accepted += 1
            if not accepted:
                break
            rotated = Mesh(current.xyz.copy(), faces, mesh.material, dict(mesh.seams))
            ratio = rotated.validate()['maximum_edge_ratio']
            if ratio < best_ratio:
                best, best_ratio = rotated, ratio
            current, report = relaxation.improve(rotated)
            ratio = current.validate()['maximum_edge_ratio']
            if ratio < best_ratio:
                best, best_ratio = current, ratio
            history.append(dict(rotations=accepted, ratio=ratio, relaxation=report))
        return best, dict(passes=len(history), history=history,
                          initial_ratio=initial['maximum_edge_ratio'], best_ratio=best_ratio,
                          construction='two-quad edge rotations with fixed rim')
