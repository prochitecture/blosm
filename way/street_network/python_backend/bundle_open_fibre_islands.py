"""Protected-bank direct quad strips for one-ended separator openings.

Only computational rows and interior vertices change. Every original curb
vertex survives, and additional boundary samples subdivide original chords.
The exported banks let the pavement builder use the identical whole edges.
"""
import math
import numpy as np
from .bundle_angular_open_islands import clip_ray
from .four_side_charts import quality
from .hole_fibre_chart import _pair_rows
from .mesh import Assembler, Mesh, cross
from .quality_correction import QuadShapeCorrection
from .regions import signed_area
from .validation import validate_overlaps


class _Bank:
    def __init__(self, stations, points):
        self.stations, self.points = stations, points

    def __call__(self, stations):
        return np.stack([np.interp(stations, self.stations, self.points[:, k]) for k in (0, 1)], axis=-1)


class _Frame:
    def __init__(self, banks, closed, index, spacing, denominator, *, width_columns=False):
        origin = .5*(banks[0][0]+banks[1][0])
        axis = .5*(banks[0][-1]+banks[1][-1])-origin
        length = float(np.linalg.norm(axis))
        if length <= 1e-8:
            raise ValueError('Open separator has no finite longitudinal extent')
        axis /= length
        left = np.array((-axis[1], axis[0]))
        first, second, self.pairing_cost = _pair_rows(*banks, origin, axis, left)
        if np.any((first[1:-1]-second[1:-1]) @ left <= 0.):
            raise ValueError('Open separator banks change their physical lateral order')
        original = .5*((first-origin) @ axis+(second-origin) @ axis)
        if np.any(np.diff(original) <= 1e-8):
            raise ValueError('Open separator needs strictly ordered matched rows')
        self.rays = tuple(_Bank(original, points) for points in (first, second))
        self.closed, self.index, self.spacing, self.denominator = closed, index, spacing, denominator
        self.width_columns = width_columns
        self.first, self.last = map(float, original[[0, -1]])
        tip_front = min(original[1]-original[0] if closed[0] else np.inf,
                        original[-1]-original[-2] if closed[1] else np.inf)
        self.transition_distance = max(.5*spacing, 1.5*tip_front)
        rows = []
        for a, b in zip(original, original[1:]):
            rows.append(a)
            cursor = a
            while cursor < b:
                if (a == self.first and closed[0]) or (b == self.last and closed[1]):
                    break
                distance = self.distance(cursor)
                width = float(np.linalg.norm(self.rays[0](cursor)-self.rays[1](cursor)))
                shear = min(distance/denominator, .28*spacing)
                across = 1 if distance < self.transition_distance else 3
                step = min(spacing, max(4.*width/across, 4.*shear/across, 1e-4))
                if b-cursor <= 1.5*step:
                    break
                cursor += step
                rows.append(cursor)
                if len(rows) > 8192:
                    raise ValueError('Open separator exceeds its finite row budget')
        rows.append(self.last)
        self.rows = np.asarray(rows)

    def distance(self, x):
        return np.minimum(x-self.first if self.closed[0] else np.inf,
                          self.last-x if self.closed[1] else np.inf)

    def columns(self, row):
        if self.width_columns:
            station = self.rows[row]
            width = float(np.linalg.norm(self.rays[0](station)-self.rays[1](station)))
            if width < .9*self.spacing:
                return 2
        return 2 if self.distance(self.rows[row]) < self.transition_distance else 6

    def address(self, row, v):
        return 'open fibre island', self.index, 'row', row, round(float(v), 12)

    def point(self, x, v):
        x, v = np.broadcast_arrays(np.asarray(x), np.asarray(v))
        source = x.copy()
        width = np.linalg.norm(self.rays[0](x)-self.rays[1](x), axis=-1)
        fade = np.maximum(0., 1.-(width/(.6*self.spacing))**2)
        cap = ((x == self.rows[1]) & self.closed[0]) | ((x == self.rows[-2]) & self.closed[1])
        fade = np.where(cap, np.maximum(fade, .5), fade)
        delta = np.zeros_like(x)
        if self.closed[0]:
            delta += np.minimum((x-self.first)/self.denominator, .28*self.spacing)*np.clip(
                (self.last-x)/self.spacing, 0., 1.)
        if self.closed[1]:
            delta -= np.minimum((self.last-x)/self.denominator, .28*self.spacing)*np.clip(
                (x-self.first)/self.spacing, 0., 1.)
        shifted = x+delta*fade*(1.-abs(2.*v-1.))
        # Interpolate both actual Cartesian banks continuously, including when
        # their protected matched points have different longitudinal positions.
        result = (1.-v[..., None])*self.rays[0](shifted)+v[..., None]*self.rays[1](shifted)
        for side, ray in ((0., self.rays[0]), (1., self.rays[1])):
            selected = v == side
            if np.any(selected):
                result[selected] = ray(source[selected])
        return result


def _assemble(frame, context, height, max_cells):
    assembly, tip_keys = Assembler(), []

    def add(xy, faces, keys):
        if len(assembly.quads)+len(faces) > max_cells:
            raise ValueError('Open separator exceeds its finite quad budget')
        assembly.add(np.c_[xy, np.full(len(xy), height)], faces, keys)

    for start, closed in ((True, frame.closed[0]), (False, frame.closed[1])):
        if not closed:
            continue
        row, root = (1, 0) if start else (len(frame.rows)-2, len(frame.rows)-1)
        profile = frame.point(frame.rows[row], np.array((0., .5, 1.)))
        tip = frame.point(np.array((frame.rows[root],)), np.array((.5,)))[0]
        order = (2, 1, 0) if start else (0, 1, 2)
        keys = [frame.address(root, .5), *(frame.address(row, (0., .5, 1.)[j]) for j in order)]
        add(np.vstack((tip, profile[list(order)])), np.array(((0, 1, 2, 3),)), keys)
        tip_keys.append(keys[0])
    first = 1 if frame.closed[0] else 0
    last = len(frame.rows)-2 if frame.closed[1] else len(frame.rows)-1
    for row in range(first, last):
        # Each half has a continuous sheared chart. Its existing 1:3 template
        # never crosses the midpoint cusp of the full computational fibre.
        for half in (0, 1):
            pattern = context.patterns.strip(1, frame.columns(row)//2, frame.columns(row+1)//2)
            u, v = pattern.uv.T
            across = .5*(half+v)
            x = (1.-u)*frame.rows[row]+u*frame.rows[row+1]
            keys = [frame.address(row if s == 0. else row+1, t) if s in (0., 1.) else
                    ('open fibre island', frame.index, 'inside', half, row, vi)
                    for vi, (s, t) in enumerate(zip(u, across))]
            add(frame.point(x, across), pattern.quads[:, ::-1], keys)
    return assembly.mesh('grass'), assembly, tip_keys


def _bank_proof(original, points, keys):
    positions = []
    for point in original:
        hits = np.flatnonzero(np.all(points == point, axis=1))
        if len(hits) != 1:
            raise ValueError('Open separator changed an original protected curb vertex')
        positions.append(int(hits[0]))
    if positions != sorted(positions) or positions[0] != 0 or positions[-1] != len(points)-1:
        raise ValueError('Open separator changed original curb order or endpoints')
    extras, error = [], 0.
    for edge, (first, last) in enumerate(zip(positions, positions[1:])):
        p, q = original[edge:edge+2]
        delta = q-p
        samples = points[first:last+1]
        t = (samples-p) @ delta/(delta @ delta)
        error = max(error, float(np.linalg.norm(samples-(p+t[:, None]*delta), axis=1).max()))
        if np.any(np.diff(t) <= 0.) or t[0] < -1e-9 or t[-1] > 1.+1e-9:
            raise ValueError('Open separator samples leave an original finite curb chord')
        extras.extend(dict(edge=edge, fraction=float(t[i]), key=keys[first+i]) for i in range(1, last-first))
    if error > 1e-8:
        raise ValueError('Open separator changed an original physical curb chord')
    return dict(original_vertex_positions=positions, added_chord_refs=extras, maximum_chord_error=error)


class OpenFibreIsland:
    def __init__(self, domain, *, max_cells=5000):
        if type(max_cells) is not int or not 1 <= max_cells <= 5000:
            raise ValueError('Open separator requires a finite existing island cell budget')
        self.domain, self.max_cells = domain, max_cells

    def build(self, first, last, roots, index):
        domain, context = self.domain, self.domain.context
        spacing = context.network.settings.spacing
        if (not math.isfinite(first+last+spacing) or first >= last or spacing <= 0.
                or (first in roots)+(last in roots) != 1):
            raise ValueError('Open fibres require one finite physical tip and one separated cut')
        banks = tuple(clip_ray(ray, first, last, roots) for ray in domain.rays[1:3])
        closed = first in roots, last in roots
        tip_vectors = tuple(bank[1]-bank[0] if closed[0] else bank[-2]-bank[-1] for bank in banks)
        intrinsic = float(abs(float(cross(*tip_vectors)))/
                          (np.linalg.norm(tip_vectors[0])*np.linalg.norm(tip_vectors[1])))
        if not math.isfinite(intrinsic) or intrinsic <= 0.:
            raise ValueError('Open separator needs a resolved positive source-tip sector')
        failures = []
        candidates = [(denominator, False) for denominator in (7., 7.5, 7.8)]
        candidates += [(denominator, True) for denominator in (7., 7.5, 7.8)]
        for denominator, width_columns in candidates:
            try:
                frame = _Frame(banks, closed, index, spacing, denominator, width_columns=width_columns)
                mesh, assembly, tips = _assemble(frame, context, domain.field.ceiling, self.max_cells)
                correction = QuadShapeCorrection(max_cells=self.max_cells, max_iterations=1000,
                    minimum_sine=min(.03, .5*intrinsic), target_sine=min(.05, .75*intrinsic))
                mesh, grading = correction.improve(mesh)
                double = quality(mesh)
                stored = quality(Mesh(mesh.xyz.astype(np.float32), mesh.quads, 'grass'))
                if max(double['maximum_edge_ratio'], stored['maximum_edge_ratio']) > 8.:
                    raise ValueError('Open fibres exceed the unchanged stored edge-ratio limit')
                lookup = assembly.lookup
                exported = []
                for side in (0, 1):
                    keys = [frame.address(row, float(side)) for row in range(len(frame.rows))]
                    if closed[0]: keys[0] = tips[0]
                    if closed[1]: keys[-1] = tips[0]
                    points = mesh.xyz[[lookup[key] for key in keys], :2].copy()
                    exported.append((points, keys))
                proofs = [_bank_proof(bank, *out) for bank, out in zip(banks, exported)]
                edges = np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2)
                unique, inverse, counts = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True)
                if len(mesh.xyz)-len(unique)+len(mesh.quads) != 1 or not set(counts).issubset({1, 2}):
                    raise ValueError('Open fibres do not form a conforming quad disk')
                boundary = edges[counts[inverse] == 1]
                following = {int(a): int(b) for a, b in boundary}
                if len(following) != len(boundary) or set(following) != set(following.values()):
                    raise ValueError('Open fibres need one manifold physical boundary at each vertex')
                at, visited = int(boundary[0, 0]), []
                while at not in visited:
                    visited.append(at); at = following[at]
                if at != visited[0] or len(visited) != len(boundary):
                    raise ValueError('Open fibres lost their one complete physical boundary')
                source_polygon = np.vstack((banks[1], banks[0][::-1]))
                origin = source_polygon[0]
                expected = signed_area(source_polygon-origin)
                actual = signed_area(mesh.xyz[visited, :2]-origin)
                if expected <= 0. or abs(actual-expected) > max(1e-8, expected*1e-10):
                    raise ValueError('Open fibres changed the original physical area')
                cuts = []
                for end, row in ((0, 0), (1, len(frame.rows)-1)):
                    if closed[end]: row = 1 if end == 0 else len(frame.rows)-2
                    values = np.linspace(1., 0., frame.columns(row)+1)
                    keys = [frame.address(row, v) for v in values]
                    cuts.append((mesh.xyz[[lookup[key] for key in keys], :2].copy(), keys))
                overlap = validate_overlaps({domain.bundle.key: mesh}, cell_size=2.)
                self.keys = [None]*len(mesh.xyz)
                for key, at in lookup.items(): self.keys[at] = key
                self.cuts = tuple(cuts)
                report = dict(construction='protected open separator with split fibre quad strips',
                    quality=double, stored_quality=stored, closed_tips=closed, original_tip_sine=intrinsic,
                    rows=len(frame.rows), maximum_columns=6, pairing_cost=frame.pairing_cost,
                    cap_shear_denominator=denominator, correction=grading, bank_sampling=proofs,
                    physical_area=expected, exact_physical_area=True, physical_curbs_unchanged=True,
                    conforming_disk=True, overlap=overlap, fibre_candidate_failures=failures,
                    tip_error=0., all_direct_quads=True,
                    flat_cut_intervals=[len(cuts[end][0])-1 if not closed[end] else None for end in (0, 1)])
                if width_columns:
                    report['width_adaptive_columns'] = dict(minimum_columns=2, maximum_columns=6,
                        threshold=.9*spacing, original_banks_unchanged=True,
                        selection='after all established column families failed')
                return mesh, exported, report
            except ValueError as error:
                failures.append(str(error))
        raise ValueError('Open separator fibre charts remain unqualified: '+str(failures))
