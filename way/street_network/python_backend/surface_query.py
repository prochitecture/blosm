"""Compact, bounded height queries on an adjusted direct-quad surface."""
import math
import numpy as np
from .mesh import cross


def inverse_quad(points, query, tolerance):
    """Invert the positive bilinear XY map without constructing triangles."""
    origin = points[0]
    a, b = points[1] - origin, points[3] - origin
    c = points[0] - points[1] + points[2] - points[3]
    delta = query - origin
    aa, bb, cc = -float(cross(a, c)), float(cross(delta, c) - cross(a, b)), float(cross(delta, b))
    scale = max(abs(float(cross(a, b))), np.linalg.norm(a) * np.linalg.norm(c), 1e-30)
    if abs(aa) <= 1e-14 * scale:
        roots = [-cc / bb] if abs(bb) > 1e-20 else []
    else:
        discriminant = bb * bb - 4. * aa * cc
        if discriminant < -1e-12 * max(bb * bb + abs(4. * aa * cc), 1e-30):
            roots = []
        else:
            q = -.5 * (bb + math.copysign(math.sqrt(max(0., discriminant)), bb))
            roots = [q / aa]
            if abs(q) > 1e-30:
                roots.append(cc / q)
    shortest = float(np.linalg.norm(np.roll(points, -1, axis=0) - points, axis=1).min())
    parameter_tolerance = max(1e-10, 4. * tolerance / shortest)
    for u in roots:
        direction = b + c * u
        length2 = float(direction @ direction)
        if length2 <= 1e-24:
            continue
        v = float((delta - a * u) @ direction / length2)
        if min(u, v) < -parameter_tolerance or max(u, v) > 1. + parameter_tolerance:
            continue
        u, v = float(np.clip(u, 0., 1.)), float(np.clip(v, 0., 1.))
        if np.linalg.norm(a * u + b * v + c * u * v - delta) <= 8. * tolerance:
            return u, v
    raise ValueError('Projected surface query could not invert a containing quad')


class QuadHeightField:
    """Keep physical heights for attachment queries after mesh publication.

    Copies exclude materials, Blender objects and temporary projection arrays.
    Queries must lie on the supplied pavement; the owner handles outside/holes.
    """
    def __init__(self, mesh, *, datum, ceiling, metadata=None, max_cells=50000,
                 max_entries=1_000_000, memory_bytes=96 * 1024 * 1024, tolerance=1e-7):
        mesh.validate()
        if (not all(map(math.isfinite, (datum, ceiling, tolerance))) or ceiling < datum or tolerance <= 0.
                or type(max_cells) is not int or max_cells < 1
                or type(max_entries) is not int or max_entries < 1
                or type(memory_bytes) is not int or memory_bytes < 4096):
            raise ValueError('Invalid projected-query budgets')
        if len(mesh.quads) > max_cells:
            raise ValueError('Projected surface query exceeds its cell budget')
        if np.any(mesh.xyz[:, 2] < datum - 1e-10) or np.any(mesh.xyz[:, 2] > ceiling + 1e-10):
            raise ValueError('Projected query heights exceed the crown envelope')
        self.datum, self.ceiling, self.tolerance = datum, ceiling, tolerance
        self.origin = mesh.xyz[:, :2].min(axis=0).copy()
        self.xy = np.ascontiguousarray(mesh.xyz[:, :2] - self.origin)
        self.heights = mesh.xyz[:, 2].copy()
        self.quads = mesh.quads.copy()
        points = self.xy[self.quads]
        self.minimum = points.min(axis=1) - tolerance
        self.maximum = points.max(axis=1) + tolerance
        span = (self.maximum - self.minimum).max(axis=1)
        self.cell_size = max(float(np.median(span)), float(span.max()) / 32., 1e-6)
        lo = np.floor(self.minimum / self.cell_size).astype(np.int64)
        hi = np.floor(self.maximum / self.cell_size).astype(np.int64)
        sizes = np.prod(hi - lo + 1, axis=1)
        entries = int(sizes.sum())
        estimated = entries * 176 + self.xy.nbytes + self.heights.nbytes + self.quads.nbytes + len(self.quads) * 64
        if sizes.max() > 4096 or entries > max_entries or estimated > memory_bytes:
            raise ValueError('Projected surface query exceeds its spatial/memory budget')
        cells = {}
        for face, (first, last) in enumerate(zip(lo, hi)):
            for x in range(int(first[0]), int(last[0]) + 1):
                for y in range(int(first[1]), int(last[1]) + 1):
                    cells.setdefault((x, y), []).append(face)
        self.cells = {key: np.asarray(value, dtype=np.int32) for key, value in cells.items()}
        self.metadata = dict(metadata or {}, height_query='adjusted bilinear quad surface',
                             query_cells=len(cells), query_entries=entries, query_estimated_bytes=estimated)
        for array in (self.origin, self.xy, self.heights, self.quads, self.minimum, self.maximum, *self.cells.values()):
            array.flags.writeable = False

    def __call__(self, xy):
        return self.evaluate(xy)

    def evaluate(self, xy, *, default=None):
        xy = np.asarray(xy, dtype=float)
        if xy.ndim != 2 or xy.shape[1] != 2 or not np.all(np.isfinite(xy)):
            raise ValueError('Projected surface queries need finite XY rows')
        if default is not None and not math.isfinite(default):
            raise ValueError('Surface query default height must be finite')
        result = np.empty(len(xy))
        tested = 0
        lower, upper = self.minimum.min(axis=0), self.maximum.max(axis=0)
        for index, point in enumerate(xy - self.origin):
            if np.any(point < lower) or np.any(point > upper):
                if default is None:
                    raise ValueError('Projected surface query is outside pavement')
                result[index] = default
                continue
            slot = tuple(np.floor(point / self.cell_size).astype(np.int64))
            candidates = self.cells.get(slot, ())
            found = False
            for start in range(0, len(candidates), 512):
                ids = np.asarray(candidates[start:start + 512], dtype=np.int32)
                tested += len(ids)
                if tested > 5_000_000:
                    raise ValueError('Projected surface queries exceed their work budget')
                valid = np.all((self.minimum[ids] <= point) & (self.maximum[ids] >= point), axis=1)
                ids = ids[valid]
                if not len(ids):
                    continue
                points = self.xy[self.quads[ids]]
                edge = np.roll(points, -1, axis=1) - points
                inside = np.all(cross(edge, point - points) >= -self.tolerance * np.linalg.norm(edge, axis=-1), axis=1)
                for local in np.flatnonzero(inside):
                    face = int(ids[local])
                    u, v = inverse_quad(points[local], point, self.tolerance)
                    height = self.heights[self.quads[face]]
                    weights = np.array(((1. - u) * (1. - v), u * (1. - v), u * v, (1. - u) * v))
                    result[index] = np.clip(weights @ height, self.datum, self.ceiling)
                    found = True
                    break
                if found:
                    break
            if not found:
                if default is None:
                    raise ValueError('Projected surface query is outside pavement')
                result[index] = default
        return result
