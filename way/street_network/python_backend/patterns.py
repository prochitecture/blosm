"""Bounded LRU of immutable quad connectivity and unit parameter coordinates.

Patterns contain no road dimensions, fixture IDs, world coordinates or heights.
Their reuse is independent of the current geometric construction algorithm.
"""
from collections import OrderedDict
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Grid:
    uv: np.ndarray
    quads: np.ndarray
    rows: int
    columns: int

    @property
    def nbytes(self):
        return self.uv.nbytes + self.quads.nbytes


class PatternLibrary:
    def __init__(self, max_bytes=8 * 1024 * 1024):
        self.max_bytes = max_bytes
        self.bytes = self.hits = self.misses = 0
        self._cache = OrderedDict()

    def grid(self, rows, columns):
        if type(rows) is not int or type(columns) is not int or min(rows, columns) < 1:
            raise ValueError('Quad subdivision counts must be positive integers')
        if (rows + 1) * (columns + 1) > 2_000_000:
            raise ValueError('Requested pattern exceeds the per-part vertex budget')
        key = ('grid', rows, columns)
        if key in self._cache:
            self.hits += 1
            self._cache.move_to_end(key)
            return self._cache[key]
        self.misses += 1
        u, v = np.meshgrid(np.linspace(0., 1., rows + 1), np.linspace(0., 1., columns + 1), indexing='ij')
        ids = np.arange((rows + 1) * (columns + 1), dtype=np.int32).reshape(rows + 1, columns + 1)
        # For XY maps, the caller's first and second parameter derivatives
        # determine winding. Mesh validation checks the resulting Jacobians.
        quads = np.stack((ids[:-1, :-1], ids[1:, :-1], ids[1:, 1:], ids[:-1, 1:]), axis=-1).reshape(-1, 4)
        uv = np.column_stack((u.ravel(), v.ravel()))
        uv.flags.writeable = quads.flags.writeable = False
        result = Grid(uv, quads, rows, columns)
        if result.nbytes <= self.max_bytes:
            while self._cache and self.bytes + result.nbytes > self.max_bytes:
                self.bytes -= self._cache.popitem(last=False)[1].nbytes
            self._cache[key] = result
            self.bytes += result.nbytes
        return result

    def stats(self):
        return dict(hits=self.hits, misses=self.misses, bytes=self.bytes,
                    limit=self.max_bytes, patterns=len(self._cache))

    def _remember(self, key, result):
        result.uv.flags.writeable = result.quads.flags.writeable = False
        if result.nbytes <= self.max_bytes:
            while self._cache and self.bytes + result.nbytes > self.max_bytes:
                self.bytes -= self._cache.popitem(last=False)[1].nbytes
            self._cache[key] = result
            self.bytes += result.nbytes
        return result

    def _lookup(self, key):
        if key in self._cache:
            self.hits += 1
            self._cache.move_to_end(key)
            return self._cache[key]
        self.misses += 1
        return None

    def strip(self, rows, first, last):
        """Direct quad strip with independently sampled end profiles.

        Local 1:3 transition blocks add two boundary intervals at a time.
        Subdivisions terminate in the boundary-adjacent cells instead of
        propagating through the entire mesh. A monotone parameter map makes
        both final profile partitions uniform. No triangles are constructed.
        """
        if first == last:
            return self.grid(rows, first)
        if (any(type(v) is not int or v < 1 for v in (rows, first, last))
                or (first - last) % 2):
            raise ValueError('A transition strip needs compatible boundary parity')
        if (rows + 1) * (min(first, last) + 1) + 2 * abs(first - last) > 2_000_000:
            raise ValueError('Requested transition exceeds the per-part vertex budget')
        key = ('strip', rows, first, last)
        saved = self._lookup(key)
        if saved is not None:
            return saved
        base = self.grid(rows, min(first, last))
        points = list(base.uv.copy())
        faces = list(map(tuple, base.quads))
        target_u = 0. if first > last else 1.
        owners = {}
        for fi, face in enumerate(faces):
            for a, b in zip(face, (*face[1:], face[0])):
                if points[a][0] == target_u and points[b][0] == target_u:
                    owners[tuple(sorted((a, b)))] = fi
        for _ in range(abs(first - last) // 2):
            edge = max(owners, key=lambda e: abs(points[e[1]][1] - points[e[0]][1]))
            owner = owners.pop(edge)
            face = faces[owner]
            at = next(i for i in range(4) if tuple(sorted((face[i], face[(i + 1) % 4]))) == edge)
            f = face[at:] + face[:at]
            a0, a3, b1, b0 = f
            a1, a2 = len(points), len(points) + 1
            points.extend(((2. * points[a0] + points[a3]) / 3.,
                           (points[a0] + 2. * points[a3]) / 3.))
            p, q = len(points), len(points) + 1
            points.extend((.5 * (points[a1] + (2. * points[b0] + points[b1]) / 3.),
                           .5 * (points[a2] + (points[b0] + 2. * points[b1]) / 3.)))
            added = ((a0, a1, p, b0), (a1, a2, q, p),
                     (a2, a3, b1, q), (p, q, b1, b0))
            indices = (owner, len(faces), len(faces) + 1, len(faces) + 2)
            faces[owner] = added[0]
            faces.extend(added[1:])
            for fi, a, b in zip(indices[:3], (a0, a1, a2), (a1, a2, a3)):
                owners[tuple(sorted((a, b)))] = fi
        uv = np.asarray(points)
        remapped = []
        for side in (0., 1.):
            boundary = np.unique(uv[uv[:, 0] == side, 1])
            remapped.append(np.interp(uv[:, 1], boundary, np.linspace(0., 1., len(boundary))))
        uv[:, 1] = (1. - uv[:, 0]) * remapped[0] + uv[:, 0] * remapped[1]
        result = Grid(uv, np.asarray(faces, dtype=np.int32), rows, max(first, last))
        self._positive(result)
        return self._remember(key, result)

    def staged_strip(self, rows, first, last, *, middle, end_rows):
        """Retain both end partitions while supplying extra interior quad rows.

        Three cached direct-quad strips share complete internal profiles.
        Curb and guide coordinates still have exactly rows + 1 uniform samples.
        """
        values = (rows, first, last, middle, end_rows)
        if (any(type(value) is not int or value < 1 for value in values)
                or 2 * end_rows >= rows or (first - middle) % 2 or (last - middle) % 2):
            raise ValueError('Staged strips need positive counts, compatible parity and a nonempty middle')
        if (rows + 1) * (max(first, last, middle) + 1) + 2 * (abs(first - middle) + abs(last - middle)) > 2_000_000:
            raise ValueError('Requested staged strip exceeds the per-part vertex budget')
        key = ('staged_strip', *values)
        saved = self._lookup(key)
        if saved is not None:
            return saved
        counts = (end_rows, rows - 2 * end_rows, end_rows)
        profiles = (first, middle, middle, last)
        points, faces, lookup, offset = [], [], {}, 0
        for block, count in enumerate(counts):
            pattern = self.strip(count, profiles[block], profiles[block + 1])
            uv = pattern.uv.copy()
            uv[:, 0] = (offset + count * uv[:, 0]) / rows
            ids = []
            for point in uv:
                address = tuple(np.round(point, 12))
                if address not in lookup:
                    lookup[address] = len(points)
                    points.append(point)
                ids.append(lookup[address])
            faces.extend(np.asarray(ids, dtype=np.int32)[pattern.quads])
            offset += count
        result = Grid(np.asarray(points), np.asarray(faces, dtype=np.int32), rows, max(profiles))
        self._positive(result)
        return self._remember(key, result)

    @staticmethod
    def _positive(pattern):
        points = pattern.uv[pattern.quads]
        edges = np.roll(points, -1, axis=1) - points
        other = np.roll(edges, -1, axis=1)
        determinants = edges[..., 0] * other[..., 1] - edges[..., 1] * other[..., 0]
        if np.min(determinants) <= 1e-13:
            raise ValueError('Direct quad transition pattern is not convex')

    def triangle(self, rows, end):
        """Three directly constructed quad blocks in a canonical triangle.

        Two sides have `rows` intervals; the opposite profile has `end`.
        Multiples of four leave compatible even partitions on all half-sides.
        Coordinates are barycentric XY in [(0,0),(1,0),(1,1)], not a triangle
        face: the final connectivity contains quads only.
        """
        if any(type(v) is not int or v < 4 or v % 4 for v in (rows, end)):
            raise ValueError('Collapsed profile patterns require counts divisible by four')
        key = ('triangle', rows, end)
        saved = self._lookup(key)
        if saved is not None:
            return saved
        a, b = rows // 2, end // 2
        corners = np.array(((0., 0.), (1., 0.), (1., 1.)))
        mids = .5 * (corners + np.roll(corners, -1, axis=0))
        center = np.mean(corners, axis=0)
        points, quads, lookup = [], [], {}
        for side in range(3):
            pattern = self.grid(a, a) if side == 0 else self.strip(a, b, a)
            uv = pattern.uv if side != 1 else pattern.uv[:, ::-1]
            faces = pattern.quads if side != 1 else pattern.quads[:, ::-1]
            u, v = uv.T
            p0, p1, p2, p3 = corners[side], mids[side], center, mids[side - 1]
            xy = ((1. - u) * (1. - v))[:, None] * p0 + (u * (1. - v))[:, None] * p1
            xy += (u * v)[:, None] * p2 + ((1. - u) * v)[:, None] * p3
            ids = []
            for point in xy:
                address = tuple(np.round(point, 12))
                if address not in lookup:
                    lookup[address] = len(points)
                    points.append(point)
                ids.append(lookup[address])
            quads.extend(np.asarray(ids)[faces])
        result = Grid(np.asarray(points), np.asarray(quads, dtype=np.int32), rows, end)
        self._positive(result)
        return self._remember(key, result)

    def quad(self, bottom, right, top, left):
        """A four-sided direct-quad atlas retaining every requested rim sample.

        Boundary partitions remain explicit (possibly nonuniform) parameters.
        Callers map them to their existing source rim vertices, so refinements
        do not create geometric cracks or force neighboring meshes to resample.
        """
        counts = (bottom, right, top, left)
        if (any(type(v) is not int or v < 1 for v in counts)
                or (bottom - top) % 2 or (right - left) % 2):
            raise ValueError('Quad rim counts require compatible opposite-edge parity')
        if ((min(bottom, top) + 1) * (min(left, right) + 1)
                + 2 * (abs(bottom - top) + abs(right - left)) > 2_000_000):
            raise ValueError('Requested quad atlas exceeds the per-part vertex budget')
        if bottom == top and left == right:
            return self.grid(bottom, left)
        key = ('quad', *counts)
        saved = self._lookup(key)
        if saved is not None:
            return saved
        base = self.grid(min(bottom, top), min(left, right))
        points, faces = list(base.uv.copy()), list(map(tuple, base.quads))
        owners = {}

        def side_of(a, b):
            if a[1] == b[1] == 0.:
                return 0
            if a[0] == b[0] == 1.:
                return 1
            if a[1] == b[1] == 1.:
                return 2
            if a[0] == b[0] == 0.:
                return 3
            return None

        def record(fi):
            face = faces[fi]
            for a, b in zip(face, (*face[1:], face[0])):
                side = side_of(points[a], points[b])
                if side is not None:
                    owners[tuple(sorted((a, b)))] = (side, fi)

        for fi in range(len(faces)):
            record(fi)
        for side, target in enumerate(counts):
            current = sum(s == side for s, fi in owners.values())
            for _ in range((target - current) // 2):
                edges = [edge for edge, (s, fi) in owners.items() if s == side]
                edge = max(edges, key=lambda e: np.linalg.norm(points[e[1]] - points[e[0]]))
                owner = owners[edge][1]
                face = faces[owner]
                for x, y in zip(face, (*face[1:], face[0])):
                    owners.pop(tuple(sorted((x, y))), None)
                at = next(i for i in range(4) if tuple(sorted((face[i], face[(i + 1) % 4]))) == edge)
                a0, a3, b1, b0 = face[at:] + face[:at]
                a1, a2 = len(points), len(points) + 1
                points.extend(((2. * points[a0] + points[a3]) / 3., (points[a0] + 2. * points[a3]) / 3.))
                p, q = len(points), len(points) + 1
                points.extend((.5 * (points[a1] + (2. * points[b0] + points[b1]) / 3.),
                               .5 * (points[a2] + (points[b0] + 2. * points[b1]) / 3.)))
                added = ((a0, a1, p, b0), (a1, a2, q, p),
                         (a2, a3, b1, q), (p, q, b1, b0))
                indices = (owner, len(faces), len(faces) + 1, len(faces) + 2)
                faces[owner] = added[0]
                faces.extend(added[1:])
                for fi in indices:
                    record(fi)
        result = Grid(np.asarray(points), np.asarray(faces, dtype=np.int32), max(bottom, top), max(left, right))
        self._positive(result)
        return self._remember(key, result)

    @staticmethod
    def boundary(pattern, side):
        """Return a rim's vertex IDs in increasing u or v, not winding order."""
        axis, value = ((1, 0.), (0, 1.), (1, 1.), (0, 0.))[side]
        ids = np.flatnonzero(pattern.uv[:, axis] == value)
        return ids[np.argsort(pattern.uv[ids, 1 - axis])]

    def band(self, rows, first, last, *, middle=4):
        """Corner coordinates with a full, half, or collapsed end profile."""
        if first and last:
            return self.strip(rows, first, last)
        key = ('band', rows, first, last, middle)
        saved = self._lookup(key)
        if saved is not None:
            return saved
        if first == last:
            if rows % 8:
                raise ValueError('A doubly capped band needs rows divisible by eight')
            left = self.band(rows // 2, 0, middle)
            right = self.band(rows // 2, middle, 0)
            points, faces, lookup = [], [], {}
            for offset, pattern in ((0., left), (.5, right)):
                uv = pattern.uv.copy()
                uv[:, 0] = offset + .5 * uv[:, 0]
                ids = []
                for point in uv:
                    address = tuple(np.round(point, 12))
                    if address not in lookup:
                        lookup[address] = len(points)
                        points.append(point)
                    ids.append(lookup[address])
                faces.extend(np.asarray(ids)[pattern.quads])
            return self._remember(key, Grid(np.asarray(points), np.asarray(faces, np.int32), rows, middle))
        triangle = self.triangle(rows, max(first, last))
        u = triangle.uv[:, 0]
        v = np.divide(triangle.uv[:, 1], u, out=np.full_like(u, .5), where=u > 1e-12)
        if first:
            u = 1. - u
        quads = triangle.quads if not first else triangle.quads[:, ::-1]
        return self._remember(key, Grid(np.column_stack((u, v)), quads, rows, max(first, last)))

    def set_budget(self, max_bytes):
        if type(max_bytes) is not int or max_bytes < 0:
            raise ValueError('Pattern memory budget must be a nonnegative integer')
        self.max_bytes = max_bytes
        while self._cache and self.bytes > max_bytes:
            self.bytes -= self._cache.popitem(last=False)[1].nbytes
