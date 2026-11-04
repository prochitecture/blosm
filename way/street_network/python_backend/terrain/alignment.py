"""Conservative, bounded-grade, no-cut placement of an entire street corridor.

No optimizer or SciPy dependency: a graph majorant coordinates shared heights,
then monotone Hermite curves provide C1 references with horizontal endpoints.
Clearance is checked on a configurable sampling of both triangulations of each
quad. This is a sampled terrain contract, not an engineering earthworks solver.
"""
from dataclasses import dataclass
import heapq
import math
import numpy as np


@dataclass(frozen=True)
class PlacementSettings:
    clearance: float = .02
    maximum_grade: float = .10
    control_spacing: float = 4.
    sample_spacing: float = 1.
    maximum_lift: float = 25.

    def validate(self):
        values = (self.clearance, self.maximum_grade, self.control_spacing,
                  self.sample_spacing, self.maximum_lift)
        if not all(math.isfinite(v) for v in values) or self.clearance < 0:
            raise ValueError('Invalid terrain placement settings')
        if min(values[1:]) <= 0 or self.maximum_grade > 1.:
            raise ValueError('Terrain grade, spacing and maximum lift must be positive; grade <= 100%')


@dataclass
class ElevationReference:
    identity: str
    stations: np.ndarray

    def __post_init__(self):
        self.stations = np.asarray(self.stations, dtype=float)
        if (not len(self.stations) or not np.isfinite(self.stations).all()
                or np.any(np.diff(self.stations) <= 0)):
            raise ValueError('Elevation reference stations must be finite and strictly increasing')


@dataclass
class CorridorComponent:
    identity: str
    reference: str
    xyz: np.ndarray
    quads: np.ndarray
    stations: np.ndarray
    datum: float = 0.
    clearance: float = 0.
    # Connections are explicit shared baseline IDs; arbitrary incidental XY
    # contact must not connect bridges or unrelated networks.
    contacts: tuple = ()  # (reference station, shared contact ID)

    def __post_init__(self):
        self.xyz = np.asarray(self.xyz, dtype=float)
        self.quads = np.asarray(self.quads, dtype=np.int32)
        self.stations = np.asarray(self.stations, dtype=float)
        if (self.xyz.shape != (len(self.stations), 3) or self.quads.ndim != 2
                or self.quads.shape[1] != 4 or not np.isfinite(self.xyz).all()
                or not np.isfinite(self.stations).all() or not math.isfinite(self.datum)
                or not math.isfinite(self.clearance) or self.clearance < 0):
            raise ValueError('Invalid street corridor component')
        if not len(self.xyz) or (self.quads.size and (self.quads.min() < 0 or self.quads.max() >= len(self.xyz))):
            raise ValueError('Corridor component has no vertices or invalid face indices')


class HeightCurve:
    def __init__(self, stations, heights):
        self.stations, self.heights = np.asarray(stations), np.asarray(heights)
        self.derivatives = np.zeros(len(stations))
        if len(stations) > 2:
            h = np.diff(self.stations)
            d = np.diff(self.heights) / h
            same = d[:-1] * d[1:] > 0
            w1, w2 = 2 * h[1:] + h[:-1], h[1:] + 2 * h[:-1]
            inner = np.zeros(len(d) - 1)
            inner[same] = (w1[same] + w2[same]) / (w1[same] / d[:-1][same] + w2[same] / d[1:][same])
            limit = 2 * np.minimum(np.abs(d[:-1]), np.abs(d[1:]))
            self.derivatives[1:-1] = np.clip(inner, -limit, limit)

    def interval(self, s):
        return np.clip(np.searchsorted(self.stations, s, side='right') - 1, 0, max(0, len(self.stations) - 2))

    def __call__(self, s):
        s = np.asarray(s)
        if len(self.stations) == 1:
            return np.full(s.shape, self.heights[0])
        i = self.interval(s)
        h = self.stations[i + 1] - self.stations[i]
        t = np.clip((s - self.stations[i]) / h, 0., 1.)
        return ((2*t**3 - 3*t**2 + 1)*self.heights[i] + (t**3 - 2*t**2 + t)*h*self.derivatives[i]
                + (-2*t**3 + 3*t**2)*self.heights[i+1] + (t**3 - t**2)*h*self.derivatives[i+1])

    def maximum_grade(self):
        if len(self.stations) == 1:
            return 0.
        d = np.diff(self.heights) / np.diff(self.stations)
        a = -6*d + 3*(self.derivatives[:-1] + self.derivatives[1:])
        b = 6*d - 4*self.derivatives[:-1] - 2*self.derivatives[1:]
        t = np.divide(-b, 2*a, out=np.zeros_like(a), where=np.abs(a) > 1e-15)
        t = np.clip(t, 0., 1.)
        return float(max(np.max(np.abs(self.derivatives)), np.max(np.abs(a*t*t + b*t + self.derivatives[:-1]))))


def quad_samples(xyz, quads, spacing):
    """Bounded batches of samples; use lower of the two possible quad diagonals.

    Yield face IDs and XYZ with the lower interpolated relief. This includes
    all vertices, edges and interiors and is conservative for either Blender
    tessellation. Refinement changes sampling, never the published quad mesh.
    """
    for start in range(0, len(quads), 512):
        q = xyz[quads[start:start+512]]
        size = np.max(np.linalg.norm(q[..., :2] - np.roll(q[..., :2], 1, axis=1), axis=2), axis=1)
        counts = np.maximum(2, np.ceil(size / spacing).astype(int))
        if np.max(counts, initial=0) > 64:
            raise ValueError('Terrain sampling requires more than 64 subdivisions of a quad; increase sampling spacing or refine the input')
        for n in np.unique(counts):
            ids = np.flatnonzero(counts == n)
            v = np.linspace(0., 1., n+1)
            u, v = np.meshgrid(v, v)
            u, v = u.ravel(), v.ravel()
            # Bilinear XY; calculate barycentric Z against each actual planar
            # triangle in XY, so nonrectangular quads are handled correctly.
            w = np.array([(1-u)*(1-v), u*(1-v), u*v, (1-u)*v]).T
            points = np.einsum('sj,fjk->fsk', w, q[ids])
            zz = []
            for triangles in (((0,1,2), (0,2,3)), ((0,1,3), (1,2,3))):
                height = np.full(points.shape[:2], np.inf)
                for tri in triangles:
                    a, b, c = (q[ids, j] for j in tri)
                    ab, ac = b[:,:2]-a[:,:2], c[:,:2]-a[:,:2]
                    p = points[:,:,:2]-a[:,None,:2]
                    det = ab[:,0]*ac[:,1]-ab[:,1]*ac[:,0]
                    if np.any(det <= 0):
                        raise ValueError('Terrain placement requires convex, counterclockwise quads')
                    wb = (p[:,:,0]*ac[:,None,1]-p[:,:,1]*ac[:,None,0])/det[:,None]
                    wc = (ab[:,None,0]*p[:,:,1]-ab[:,None,1]*p[:,:,0])/det[:,None]
                    z = a[:,None,2]+wb*(b-a)[:,None,2]+wc*(c-a)[:,None,2]
                    valid = (wb >= -1e-7) & (wc >= -1e-7) & (wb+wc <= 1+1e-7)
                    height = np.minimum(height, np.where(valid, z, np.inf))
                zz.append(height)
            points[:,:,2] = np.minimum(*zz)
            if not np.isfinite(points).all():
                raise ValueError('Could not sample a terrain footprint quad')
            yield start+ids, points


@dataclass
class PlacementResult:
    curves: dict
    coordinates: dict
    report: dict


def fit(references, components, terrain, settings=PlacementSettings()):
    """terrain(xy[N,2]) returns original-terrain Z[N], or raises on any miss."""
    settings.validate()
    refs = {r.identity: r for r in references}
    if len(refs) != len(references):
        raise ValueError('Duplicate elevation reference')
    if len({c.identity for c in components}) != len(components):
        raise ValueError('Duplicate corridor component')
    if any(c.reference not in refs for c in components):
        raise ValueError('Corridor component refers to an unknown elevation reference')
    node_ids, parent, floors = {}, [], []
    for key, ref in refs.items():
        node_ids[key] = np.arange(len(parent), len(parent)+len(ref.stations))
        parent.extend(range(len(parent), len(parent)+len(ref.stations)))
        floors.extend([-np.inf]*len(ref.stations))
    parent, floors = np.asarray(parent), np.asarray(floors)
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    def join(a, b):
        a, b = root(a), root(b)
        parent[max(a,b)] = min(a,b)
    contacts, mappings, sample_count = {}, {}, 0
    for c in components:
        ref = refs[c.reference]
        if np.min(c.stations) < ref.stations[0]-1e-6 or np.max(c.stations) > ref.stations[-1]+1e-6:
            raise ValueError('Component stations lie outside its elevation reference: '+c.identity)
        for s, contact in c.contacts:
            j = int(np.argmin(np.abs(ref.stations-s)))
            if abs(ref.stations[j]-s) > 1e-6:
                raise ValueError('A contact must lie on an elevation control')
            node = node_ids[c.reference][j]
            if contact in contacts:
                join(node, contacts[contact])
            contacts[contact] = node
        lo = np.clip(np.searchsorted(ref.stations, c.stations, side='right')-1, 0, max(0,len(ref.stations)-2))
        hi = np.minimum(lo+1, len(ref.stations)-1)
        mappings[c.identity] = (node_ids[c.reference][lo], node_ids[c.reference][hi])
        for ids, samples in quad_samples(c.xyz, c.quads, settings.sample_spacing):
            height = np.asarray(terrain(samples[:,:,:2].reshape(-1,2))).reshape(samples.shape[:2])
            if not np.isfinite(height).all():
                raise ValueError('Terrain does not cover the complete corridor: '+c.identity)
            required = np.max(height - (samples[:,:,2]-c.datum), axis=1) + settings.clearance+c.clearance
            # Raise all controls influencing every vertex of this face. The
            # final face, not only the analytical reference, meets these bounds.
            for nodes in mappings[c.identity]:
                np.maximum.at(floors, nodes[c.quads[ids]].ravel(), np.repeat(required, 4))
            sample_count += height.size
    roots = np.array([root(i) for i in range(len(parent))])
    merged = np.full(len(parent), -np.inf)
    np.maximum.at(merged, roots, floors)
    edges = [[] for _ in parent]
    for key, ref in refs.items():
        ns = roots[node_ids[key]]
        for a,b,distance in zip(ns[:-1], ns[1:], np.diff(ref.stations)):
            # Capped Hermite derivatives are bounded by twice the largest
            # adjacent secant. This reserve guarantees the requested grade.
            cost = .5*settings.maximum_grade*distance
            if a != b:
                edges[a].append((b,cost)); edges[b].append((a,cost))
    queue = [(-h,i) for i,h in enumerate(merged) if np.isfinite(h)]
    heapq.heapify(queue)
    while queue:
        negative, a = heapq.heappop(queue)
        if -negative < merged[a]-1e-12:
            continue
        for b,cost in edges[a]:
            value = merged[a]-cost
            if value > merged[b]+1e-12:
                merged[b] = value
                heapq.heappush(queue,(-value,b))
    curves = {key: HeightCurve(ref.stations, merged[roots[node_ids[key]]]) for key,ref in refs.items()}
    if any(not np.isfinite(c.heights).all() for c in curves.values()):
        raise ValueError('An elevation reference has no terrain constraints')
    grade = max((c.maximum_grade() for c in curves.values()), default=0.)
    if grade > settings.maximum_grade+1e-9:
        raise ValueError('Terrain reference exceeds the requested longitudinal grade')
    coordinates, maximum_lift = {}, 0.
    for c in components:
        xyz = c.xyz.copy()
        xyz[:,2] += curves[c.reference](c.stations)-c.datum
        ground = np.asarray(terrain(xyz[:,:2]))
        lift = float(np.max(xyz[:,2]-ground))
        maximum_lift = max(maximum_lift, lift)
        if lift > settings.maximum_lift + settings.clearance+c.clearance:
            raise ValueError(f'{c.identity}: terrain fit lifts the surface {lift:.2f} m above ground; limit is {settings.maximum_lift:.2f} m. Reduce the protected region or revise the alignment/terrain settings')
        coordinates[c.identity] = xyz
    return PlacementResult(curves, coordinates, dict(references=len(refs), samples=sample_count,
        maximum_reference_grade=grade, maximum_surface_lift=maximum_lift,
        sample_spacing=settings.sample_spacing, clearance=settings.clearance,
        policy='sampled no-cut; unchanged original terrain'))
