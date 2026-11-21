"""Use one certified polyline for an original bank shared by a corner and a port."""
from dataclasses import replace
from hashlib import sha256
import numpy as np

from .compound_curbs import OffsetCurb
from .compound_dimensional_corners import DimensionalCompositeCurb, DimensionalCorner, _patch
from .source_root_ports import NativeRolePort


def _snap_samples(values, knots, points, result):
    """Use the canonical coordinate at an exact sample, including the last knot."""
    at = np.asarray(values)
    indices = np.minimum(np.searchsorted(knots, at), len(knots)-1)
    return np.where((at == knots[indices])[..., None], points[indices], result)


class SharedOffsetCurb(OffsetCurb):
    def __init__(self, original, stations, points, error):
        self.approach, self.offset = original.approach, original.offset
        if original.parameters[-1] < original.parameters[0]:
            stations, points = stations[::-1], points[::-1]
        self.parameters = stations.copy()
        self.points = points.copy()
        self.cumulative = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
        self.length = float(self.cumulative[-1])
        self.knots = self.cumulative/self.length
        self.error = error
        for data in (self.parameters, self.points, self.cumulative, self.knots):
            data.flags.writeable = False

    def __call__(self, t):
        xy = self.approach.point(np.interp(t, self.knots, self.parameters), self.offset)
        return _snap_samples(t, self.knots, self.points, xy)

    def sample(self, sampler, first, last):
        if first == last:
            return np.array([first])
        if first > last:
            return self.sample(sampler, last, first)[::-1]
        if (sampler.tolerance + 1e-14 < self.error
                or np.max(np.diff(self.cumulative)) > sampler.spacing + 1e-10):
            raise ValueError('The shared bank needs finer sampling on both source polygons')
        at = np.r_[first, self.knots[(self.knots > first) & (self.knots < last)], last]
        if len(at) > sampler.max_samples:
            raise ValueError('Shared bank exceeds the requested sample budget')
        return at


def _native(traversal, at):
    return at if traversal.forward else traversal.path.length-at


def _bounds(traversal, offset, at, scale):
    """Bounds for the sagitta of each chord on the actual offset curve."""
    native = _native(traversal, at)
    path = traversal.path
    reach = abs(offset)
    bounds = []
    for a, b, length in zip(native, native[1:], np.diff(at)):
        k = 0. if path.exact else path.visible_curvature_bound(min(a,b), max(a,b))
        if reach*k >= 1.:
            raise ValueError('Shared bank has an irregular offset')
        bounds.append(k/(1.-reach*k) * ((1.+reach*k)*length)**2 / 8.)
    return np.asarray(bounds) + 256.*np.finfo(float).eps*max(1., scale)


def _grid(piece, curve_stations, port, sources, settings, tolerance):
    arm = piece.approach
    traversal, = arm.path.pieces
    low, high = sorted((float(piece.parameters[0]), float(piece.parameters[-1])))
    if not 0. <= low < high <= port.layout.outward_stations[0]:
        return None
    coordinates = port.layout.boundaries[0].source_coordinates
    old = coordinates.stations
    port_stations = np.r_[low, old[(low < old) & (old < high)], high]
    curve_stations = np.unique(np.r_[low, curve_stations, high])
    if np.array_equal(port_stations, curve_stations):
        return None
    scale = float(np.max(np.abs(port.region.xy)))
    target = min(tolerance, coordinates.maximum_chord_error,
                 _bounds(traversal, piece.offset, port_stations, scale).max(),
                 _bounds(traversal, piece.offset, curve_stations, scale).max())
    hard = {low, high}
    for value in _native(traversal, traversal.path.station(traversal.path.knots)):
        if low < value < high:
            hard.add(float(value))
    for source in sources:
        for corner in source.layout.corners:
            for a, station in zip((corner.first, corner.second), corner.stations):
                if a is arm and low < station < high:
                    hard.add(float(station))
    hard = sorted(hard)
    chunks = []
    count = 1
    for a, b in zip(hard, hard[1:]):
        first, last = sorted((_native(traversal, a), _native(traversal, b)))
        knots = traversal.path.samples(first, last, min(.5, settings.spacing), tolerance=.99*target,
                                       local_curvature=True, max_samples=20000)
        outward = knots.copy() if traversal.forward else traversal.path.length-knots[::-1]
        outward[0], outward[-1] = a, b
        count += len(outward)-1
        if count > 20000:
            raise ValueError('Shared bank exceeds its source sampling budget')
        chunks.append(outward[:-1])
    at = np.concatenate([*chunks, [high]])
    for _ in range(18):
        bounds = _bounds(traversal, piece.offset, at, scale)
        split = bounds > target
        if not np.any(split):
            if np.any(np.linalg.norm(np.diff(arm.point(at,piece.offset),axis=0),axis=1) <= 1e-8):
                raise ValueError('Distinct shared bank stations collapse to one position')
            return at, float(bounds.max())
        if len(at) + int(split.sum()) > 20000:
            break
        at = np.sort(np.r_[at, ((at[:-1]+at[1:])*.5)[split]])
    raise ValueError('Shared bank exceeded its analytical refinement budget')


def normalize(context, sources, partitions, settings):
    """Rebuild only the two representations of the same finite native interval."""
    from .compound_cycles import _port_at_stations
    from .intersections import Approach
    from .routes import Route
    replacements, grids, bindings, signatures = {}, dict(partitions), [], []
    occupied = {}
    for source in sources:
        for index, curve in enumerate(source.layout.boundaries):
            if type(curve) is not DimensionalCompositeCurb:
                continue
            corner, = source.layout.corners
            record = getattr(corner, 'dimensional_corner', None)
            if not isinstance(record, DimensionalCorner) or record.source is not source:
                continue
            pieces = list(curve.pieces)
            labels = [s for s in source.segments if s.kind=='curb' and s.index==index]
            parameters = np.asarray(sorted({u for label in labels for u in (label.first,label.last)}))
            for i, piece in enumerate(curve.pieces):
                if type(piece) is not OffsetCurb:
                    continue
                arm = piece.approach
                if (type(arm) is not Approach or type(arm.path) is not Route or len(arm.path.pieces)!=1
                        or arm.offset!=0. or abs(piece.offset)!=.5*arm.width):
                    continue
                matches = [p for p in sources if len(p.layout.approaches)==1 and p.layout.approaches[0] is arm]
                if len(matches)!=1:
                    continue
                owner = matches[0].region.owner
                port = replacements.get(owner, matches[0])
                NativeRolePort.build(context, port, arm)
                low, high = curve.lengths[i:i+2]/curve.length
                local = parameters[(parameters>=low) & (parameters<=high)]*curve.length-curve.lengths[i]
                shared = _grid(piece,np.interp(local,piece.cumulative,piece.parameters),port,sources,
                               settings,curve.prepared_tolerance)
                if shared is None:
                    continue
                stations, bound = shared
                a,b = stations[[0,-1]]
                if any(max(a,x)<min(b,y) for x,y in occupied.get(owner,())):
                    raise ValueError('Overlapping shared bank grids require a common subdivision')
                occupied.setdefault(owner,[]).append((a,b))
                old = port.layout.boundaries[0].source_coordinates.stations
                full = np.r_[old[old<a], stations, old[old>b]]
                port = _port_at_stations(owner,arm,port.layout.outward_stations[0],port.layout.radial,
                                         port.layout.boundaries[0].source_coordinates.maximum_chord_error,full)
                bank = port.layout.boundaries[0 if piece.offset>0. else 1]
                xy = bank.points if piece.offset>0. else bank.points[::-1]
                pieces[i] = SharedOffsetCurb(piece,stations,xy[np.searchsorted(full,stations)],bound)
                replacements[owner] = port
                traversal, = arm.path.pieces
                sid = traversal.section.id
                native = np.sort(_native(traversal,stations))
                prior = grids[sid]
                grids[sid] = np.unique(np.r_[prior[(prior<native[0]) | (prior>native[-1])],native])
                signatures.append((source.region.owner.text(),owner.text(),index,i,sid,arm.end,piece.offset,
                                   float(a),float(b),len(stations),bound,sha256(pieces[i].points.tobytes()).hexdigest(),
                                   sha256(stations.tobytes()).hexdigest()))
            if all(a is b for a,b in zip(pieces,curve.pieces)):
                continue
            updated = DimensionalCompositeCurb(curve.corner,pieces,curve.prepared_tolerance)
            updated.shared_bank_points = tuple(((updated.lengths[i]+p.knots*p.length)/updated.length,p.points)
                                               for i,p in enumerate(pieces) if isinstance(p,SharedOffsetCurb))
            rebuilt = _patch(source.region.owner,corner,updated,tuple(corner.stations),
                             replace(settings,curve_tolerance=curve.prepared_tolerance))
            rebuilt.layout.original_junction = source.layout.original_junction
            replacements[source.region.owner] = rebuilt
            bindings.append((corner,replace(record,source=rebuilt,curve=updated)))
    result = tuple(replacements.get(s.region.owner,s) for s in sources)
    if sum(len(s.region.xy) for s in result)>20000:
        raise ValueError('Shared source polygons exceed the existing sample budget')
    return result,grids,bindings,tuple(signatures)
