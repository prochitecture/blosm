"""Bounded point classification for one unchanged source-width proof.

Y bins only select candidates. Final ray and distance arithmetic, epsilon and
nearest-edge tie order remain those of regions._classify. The caller owns the
polygon's immutable lifetime; construct a new classifier for every new proof.
"""
import math
import numpy as np
from .regions import _classify, _classify_edges


class _PreparedPolygon:
    maximum_edges = 20000
    maximum_bins = 1024
    maximum_entries = 131072
    maximum_bytes = 4*1024*1024
    arithmetic_limit = math.sqrt(np.finfo(float).max)/16.

    def __init__(self, polygon, remaining_bytes):
        self.polygon, self.bytes, self.ready = polygon, 0, False
        a = np.asarray(polygon)
        if (a.ndim != 2 or a.shape[1] != 2 or not 3 <= len(a) <= self.maximum_edges
                or a.dtype != np.float64):
            return
        # Check the bounded original edge arrays before allocating them.
        budget = min(remaining_bytes, self.maximum_bytes)
        base_bytes = len(a)*88
        if base_bytes > budget or not np.isfinite(a).all():
            return
        self.a, self.b = a, np.roll(a, -1, axis=0)
        self.edge = self.b-a
        self.length2 = np.sum(self.edge*self.edge, axis=1)
        self.low, self.high = np.minimum(a, self.b), np.maximum(a, self.b)
        self.magnitude = max(1., float(np.max(abs(a))), float(np.max(abs(self.edge))))
        if (np.any(self.length2 <= 0.) or not np.isfinite(self.length2).all()
                or self.magnitude > self.arithmetic_limit):
            return
        self.y0, self.y1 = float(a[:, 1].min()), float(a[:, 1].max())
        height = self.y1-self.y0
        rounding = 512.*np.finfo(float).eps*self.magnitude
        if not height > 8.*rounding or not math.isfinite(height):
            return
        self.count = min(self.maximum_bins, max(32, 8*math.isqrt(len(a))))
        self.scale = self.count/height
        # One outward bin protects floor/index arithmetic; it changes no test.
        first = np.clip(np.floor((self.low[:, 1]-self.y0)*self.scale).astype(np.int64)-1, 0, self.count-1)
        last = np.clip(np.floor((self.high[:, 1]-self.y0)*self.scale).astype(np.int64)+1, 0, self.count-1)
        self.entries = int(np.sum(last-first+1))
        planned_bytes = base_bytes+4*self.entries+128*self.count
        if self.entries > self.maximum_entries or planned_bytes > budget:
            return
        self.bins = tuple(np.flatnonzero((first <= i) & (last >= i)).astype(np.int32)
                          for i in range(self.count))
        self.bytes, self.ready = planned_bytes, True

    def classify(self, points, epsilon):
        points = np.asarray(points)
        if (points.ndim != 2 or points.shape[1] != 2 or points.dtype != np.float64
                or not np.isfinite(points).all() or not math.isfinite(epsilon) or epsilon <= 0.
                or epsilon > self.arithmetic_limit
                or (len(points) and float(np.max(abs(points))) > self.arithmetic_limit)):
            return _classify(points, self.polygon, epsilon)
        inside = np.zeros(len(points), bool)
        boundary = np.full(len(points), -1, np.int32)
        if len(points) >= 8:
            # All crossings for these rays, and every possible epsilon hit,
            # lie in this Y interval. Preserve original edge order for ties.
            batch = max(1, min(64, 32768 // len(self.a)))
            for at in range(0, len(points), batch):
                p = points[at:at+batch]
                reach = epsilon + 512.*np.finfo(float).eps*max(self.magnitude, float(np.max(abs(p))))
                ids = np.flatnonzero((self.high[:, 1] >= p[:, 1].min()-reach)
                                     & (self.low[:, 1] <= p[:, 1].max()+reach))
                if len(ids):
                    hit, near = _classify_edges(p, self.a[ids], self.b[ids], epsilon)
                    inside[at:at+len(p)] = hit
                    boundary[at:at+len(p)] = np.where(near < 0, -1, ids[near])
            return inside, boundary
        for i, p in enumerate(points):
            rounding = 512.*np.finfo(float).eps*max(self.magnitude, float(np.max(abs(p))))
            reach = epsilon+rounding
            first_y, last_y = float(p[1]-reach), float(p[1]+reach)
            if last_y < self.y0 or first_y > self.y1:
                continue
            first = max(0, min(self.count-1, int(math.floor((first_y-self.y0)*self.scale))))
            last = max(0, min(self.count-1, int(math.floor((last_y-self.y0)*self.scale))))
            ids = (self.bins[first] if first == last else
                   np.unique(np.concatenate(self.bins[first:last+1])))
            if not len(ids):
                continue
            a, b, edge = self.a[ids], self.b[ids], self.edge[ids]
            crosses = (a[:, 1] > p[1]) != (b[:, 1] > p[1])
            offset = np.divide((p[1]-a[:, 1])*edge[:, 0], edge[:, 1],
                               out=np.zeros(len(ids)), where=crosses)
            intersections = a[:, 0]+offset
            inside[i] = np.count_nonzero(crosses & (p[0] < intersections)) % 2 == 1
            possible = np.all(p >= self.low[ids]-reach, axis=1) & np.all(p <= self.high[ids]+reach, axis=1)
            near = ids[possible]
            if not len(near):
                continue
            delta = p-self.a[near]
            fraction = np.clip(np.sum(delta*self.edge[near], axis=-1)/self.length2[near], 0., 1.)
            distance2 = np.sum((delta-fraction[:, None]*self.edge[near])**2, axis=-1)
            nearest = int(np.argmin(distance2))
            # Sorted ORIGINAL edge IDs preserve the baseline argmin tie.
            if distance2[nearest] <= epsilon**2:
                boundary[i] = near[nearest]
        return inside, boundary


class SourceProofClassifier:
    """No global state, eviction or reuse after a polygon edit.

    Buffer addresses distinguish equal-shape source polygons. Strong references
    pin those allocations while fresh views of the same buffer share an index.
    Declined preparations are remembered within the same bounded inventory.
    """
    maximum_polygons = 16
    maximum_cache_bytes = 8*1024*1024

    def __init__(self):
        self._cache = {}
        self.allocated_bytes = 0
        self.preparations = 0

    def __call__(self, points, polygon, epsilon):
        if not isinstance(polygon, np.ndarray):
            return _classify(points, polygon, epsilon)
        key = (polygon.__array_interface__['data'][0], polygon.shape, polygon.strides, polygon.dtype.str)
        entry = self._cache.get(key)
        if entry is None:
            if len(self._cache) >= self.maximum_polygons:
                return _classify(points, polygon, epsilon)
            prepared = _PreparedPolygon(polygon, self.maximum_cache_bytes-self.allocated_bytes)
            self.preparations += 1
            if not prepared.ready:
                prepared = None
            else:
                self.allocated_bytes += prepared.bytes
            # Pin declined buffers too, so their address cannot be recycled.
            entry = self._cache[key] = (polygon, prepared)
        owner, prepared = entry
        return (_classify(points, polygon, epsilon) if prepared is None else
                prepared.classify(points, epsilon))
