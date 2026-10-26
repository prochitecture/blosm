"""Compact numerical mesh payloads and independent local geometry guards."""
from dataclasses import dataclass, field
import numpy as np


def cross(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def publication_coordinates(values, scale, *, physical=False):
    """Use the same float32 scene coordinates for writing and validation."""
    if not np.isfinite(scale) or scale <= 0.:
        raise ValueError('Publication needs a finite positive scene scale')
    with np.errstate(over='ignore', invalid='ignore'):
        result = np.asarray(np.asarray(values, dtype=float) / scale, dtype=np.float32)
    if not np.all(np.isfinite(result)):
        raise ValueError('Coordinates exceed Blender float32 publication range')
    return result.astype(float) * scale if physical else result


@dataclass
class Mesh:
    xyz: np.ndarray
    quads: np.ndarray
    material: str = 'pavement'
    seams: dict = field(default_factory=dict)

    def __post_init__(self):
        self.xyz = np.ascontiguousarray(self.xyz, dtype=np.float64)
        self.quads = np.ascontiguousarray(self.quads, dtype=np.int32)

    @property
    def nbytes(self):
        return self.xyz.nbytes + self.quads.nbytes

    @property
    def bounds(self):
        return tuple(map(float, self.xyz[:, :2].min(axis=0))), tuple(map(float, self.xyz[:, :2].max(axis=0)))

    def validate(self):
        if (self.xyz.ndim != 2 or self.xyz.shape[1] != 3 or not len(self.xyz)
                or not np.all(np.isfinite(self.xyz)) or self.quads.ndim != 2 or self.quads.shape[1] != 4
                or not len(self.quads) or self.quads.min() < 0 or self.quads.max() >= len(self.xyz)):
            raise ValueError('Invalid numerical mesh buffers')
        points = self.xyz[self.quads, :2]
        edges = np.roll(points, -1, axis=1) - points
        jacobians = cross(edges, np.roll(edges, -1, axis=1))
        if np.min(jacobians) <= 1e-10:
            raise ValueError('A generated quad is collapsed, inverted or concave')
        edge_lengths = np.linalg.norm(edges, axis=-1)
        ratio = np.max(edge_lengths, axis=1) / np.min(edge_lengths, axis=1)
        return dict(vertices=len(self.xyz), quads=len(self.quads),
                    minimum_corner=float(jacobians.min()), maximum_edge_ratio=float(ratio.max()))

    def quad_slopes(self, exaggeration=1.):
        """Independent physical gradients on a 3x3 grid in each bilinear quad.

        These are measured mesh slopes, not a proof about unsampled analytical
        surfaces. Batches bound scratch memory for large imported networks.
        Inspection exaggeration is removed before solving each 2D tangent plane.
        No triangle mesh is produced or used for this diagnostic.
        """
        if exaggeration <= 0.:
            return np.zeros(len(self.quads))
        maxima = np.zeros(len(self.quads))
        for start in range(0, len(self.quads), 4096):
            p = self.xyz[self.quads[start:start + 4096]].copy()
            p[:, :, 2] /= exaggeration
            measured = np.zeros(len(p))
            for u in (0., .5, 1.):
                for v in (0., .5, 1.):
                    a = (1. - v) * (p[:, 1] - p[:, 0]) + v * (p[:, 2] - p[:, 3])
                    b = (1. - u) * (p[:, 3] - p[:, 0]) + u * (p[:, 2] - p[:, 1])
                    determinant = cross(a[:, :2], b[:, :2])
                    if np.any(abs(determinant) <= 1e-12):
                        raise ValueError('A quad tangent plane is singular during slope measurement')
                    gx = (a[:, 2] * b[:, 1] - b[:, 2] * a[:, 1]) / determinant
                    gy = (a[:, 0] * b[:, 2] - b[:, 0] * a[:, 2]) / determinant
                    measured = np.maximum(measured, np.hypot(gx, gy))
            maxima[start:start + len(p)] = measured
        return maxima

    def slope_statistics(self, exaggeration=1.):
        maxima = self.quad_slopes(exaggeration)
        return dict(maximum=float(maxima.max()), percentile_99=float(np.quantile(maxima, .99)), samples_per_quad=9)


class Assembler:
    """Weld only declared parametric keys, never coincident unrelated inputs."""
    def __init__(self):
        self.xyz = []
        self.quads = []
        self.lookup = {}

    def add(self, points, faces, keys=None):
        ids = []
        for i, point in enumerate(points):
            key = keys[i] if keys is not None else None
            at = self.lookup.get(key) if key is not None else None
            if at is None:
                at = len(self.xyz)
                self.xyz.append(point)
                if key is not None:
                    self.lookup[key] = at
            elif np.linalg.norm(np.asarray(self.xyz[at]) - point) > 1e-8:
                raise ValueError('Two patches disagree on an explicitly shared seam')
            ids.append(at)
        self.quads.extend(np.asarray(ids, dtype=np.int32)[faces])

    def mesh(self, material='pavement'):
        return Mesh(np.asarray(self.xyz), np.asarray(self.quads), material)
