"""Piecewise polynomial width transitions retaining Exact source centerlines.

Only the transition interior receives common miter frames. Source paths and
their full outer cut profiles stay authoritative and are never rewritten.
"""
import math
import numpy as np
from .continuations import _evaluate
from .paths import right, unit


def exact_path(path):
    value = getattr(path, 'exact', None)
    return bool(value) if value is not None else bool(getattr(path, 'pieces', ())) and all(
        exact_path(piece.path) for piece in path.pieces)


def split(controls, t):
    work = np.asarray(controls).copy()
    first, last = [work[0]], [work[-1]]
    while len(work) > 1:
        work = (1.-t)*work[:-1]+t*work[1:]
        first.append(work[0]); last.append(work[-1])
    return np.asarray(first), np.asarray(last[::-1])


def restricted(controls, first, last):
    left, _ = split(controls, last)
    return split(left, first/last)[1] if first else left


class PolynomialCurb:
    def __init__(self, controls):
        self.controls = np.array(controls, dtype=float, copy=True)
        degree = len(self.controls)-1
        self.first_derivative = degree*np.diff(self.controls, axis=0)
        self.second_derivative = degree*(degree-1)*np.diff(self.controls, n=2, axis=0)
        self.maximum_speed = float(np.linalg.norm(self.first_derivative, axis=1).max())
        self.maximum_second = float(np.linalg.norm(self.second_derivative, axis=1).max())
        for values in (self.controls, self.first_derivative, self.second_derivative):
            values.flags.writeable = False

    def __call__(self, t):
        return _evaluate(self.controls, t)

    def intervals(self, spacing, tolerance):
        count = max(2, int(math.ceil(self.maximum_speed/spacing)),
                    int(math.ceil(math.sqrt(self.maximum_second/(8.*tolerance)))))
        return count+count % 2


class PiecewiseCurb:
    """Shared mandatory knots with independent Bernstein polynomial pieces."""
    def __init__(self, knots, segments):
        self.knots = np.array(knots, dtype=float, copy=True)
        self.segments = tuple(segments)
        if (len(self.knots) != len(self.segments)+1 or self.knots[0] != 0.
                or self.knots[-1] != 1. or np.any(np.diff(self.knots) <= 0.)):
            raise ValueError('Piecewise curb needs positive intervals covering zero to one')
        self.knots.flags.writeable = False

    def evaluate(self, t, derivative=False):
        t = np.asarray(t)
        if not np.all(np.isfinite(t)) or np.any(t < -1e-12) or np.any(t > 1.+1e-12):
            raise ValueError('Piecewise curb parameter is outside zero to one')
        flat = np.clip(t.ravel(), 0., 1.)
        indices = np.clip(np.searchsorted(self.knots, flat, side='right')-1, 0, len(self.segments)-1)
        result = np.empty((len(flat), 2))
        for i in np.unique(indices):
            mask = indices == i
            span = self.knots[i+1]-self.knots[i]
            local = (flat[mask]-self.knots[i])/span
            piece = self.segments[i]
            result[mask] = (_evaluate(piece.first_derivative, local)/span if derivative else piece(local))
        return result.reshape(t.shape+(2,))

    def __call__(self, t):
        return self.evaluate(t)

    def derivative(self, t):
        return self.evaluate(t, derivative=True)

    def sampling_parameters(self, spacing, tolerance, *, max_intervals=200000):
        return paired_parameters((self,), spacing, tolerance, max_intervals=max_intervals)


def paired_parameters(curbs, spacing, tolerance, *, max_intervals):
    if not all(math.isfinite(v) and v > 0. for v in (spacing, tolerance)):
        raise ValueError('Piecewise curb sampling needs positive finite limits')
    knots = curbs[0].knots
    if any(not np.array_equal(c.knots, knots) for c in curbs):
        raise ValueError('Paired curbs need the same mandatory parameter knots')
    counts = [max(c.segments[i].intervals(spacing, tolerance) for c in curbs) for i in range(len(knots)-1)]
    if sum(counts) > max_intervals:
        raise ValueError('Piecewise curb samples exceed their allocation budget')
    return np.r_[np.concatenate([np.linspace(a, b, count+1)[:-1]
                                for a,b,count in zip(knots,knots[1:],counts)]), 1.]


def outward_stations(approach, distance):
    path = approach.path
    if getattr(path, 'outward', False):
        controls = np.concatenate([path.stations[i]+(piece.path.knots if piece.forward else
            piece.path.length-piece.path.knots) for i,piece in enumerate(path.pieces)])
    else:
        controls = path.knots if approach.end == 0 else path.length-path.knots
    return np.unique(np.r_[0., controls[(controls > 1e-10) & (controls < distance-1e-10)], distance])


class ExactTransitionGeometry:
    """An immutable-centerline ruled chart across one dimensional seam."""
    @staticmethod
    def source_points(approaches, distances):
        first, last = approaches
        stations = tuple(outward_stations(a, d) for a,d in zip(approaches, distances))
        p0, _, n0 = first.frame(stations[0])
        p1, _, n1 = last.frame(stations[1])
        points = np.vstack((p0[::-1], p1[1:]))
        normals = np.vstack((-n0[::-1], n1[1:]))
        metric = np.r_[distances[0]-stations[0][::-1], distances[0]+stations[1][1:]]
        return points, normals, metric, len(stations[0])-1

    @classmethod
    def required(cls, approaches, distances):
        points, _, _, _ = cls.source_points(approaches, distances)
        directions = unit(np.diff(points, axis=0))
        return bool(np.any(np.linalg.norm(np.diff(directions, axis=0), axis=1) > 1e-8))

    @classmethod
    def build(cls, approaches, distances, *, coherent_collinear=False):
        if not all(exact_path(a.path) for a in approaches):
            raise ValueError('Exact miter transition needs two Exact source paths')
        points, normals, metric, seam = cls.source_points(approaches, distances)
        if len(points) > 4097:
            raise ValueError('Exact transition exceeds its source-piece budget')
        weld_error = float(np.linalg.norm(approaches[0].origin-approaches[1].origin))
        if weld_error > 1e-5:
            raise ValueError('Exact seam endpoint copies disagree')
        points[seam] = .5*(approaches[0].origin+approaches[1].origin)
        directions = unit(np.diff(points, axis=0))
        dot = np.sum(directions[:-1]*directions[1:], axis=1)
        if np.any(dot < -.95):
            raise ValueError('Exact transition has a near-reversal miter')
        normals[1:-1] = (right(directions[:-1])+right(directions[1:]))/(1.+dot)[:,None]
        if coherent_collinear:
            # Collinear authoring controls do not define additional turns.
            # Interpolate their derived frames between real bends and the
            # authoritative outer cuts, so resampling a straight source does
            # not force the shared miter into a microscopic final interval.
            fixed = np.r_[True,np.linalg.norm(np.diff(directions,axis=0),axis=1)>1e-8,True]
            for axis in (0,1):
                normals[:,axis] = np.interp(metric,metric[fixed],normals[fixed,axis])
        reach = float(np.linalg.norm(normals, axis=1).max())
        if reach > 5.:
            raise ValueError('Exact transition exceeds the supported miter reach')
        # Remove only redundant affine P/N controls. The ownership seam stays
        # compulsory, and every actual centerline/frame bend remains present.
        keep = []
        for index in range(len(points)):
            keep.append(index)
            while len(keep) >= 3 and keep[-2] != seam:
                a,b,c = keep[-3:]
                t = (metric[b]-metric[a])/(metric[c]-metric[a])
                if (np.linalg.norm(points[b]-(1.-t)*points[a]-t*points[c]) > 1e-10
                        or np.linalg.norm(normals[b]-(1.-t)*normals[a]-t*normals[c]) > 1e-12):
                    break
                keep.pop(-2)
        points, normals, metric = points[keep], normals[keep], metric[keep]
        knots = metric/sum(distances)
        knots[0], knots[-1] = 0., 1.
        curbs = []
        first, last = approaches
        for q in (-1., 1.):
            a, b = -first.offset+q*.5*first.width, last.offset+q*.5*last.width
            distance = np.array((a,a,a,b,b,b))
            offsets = _evaluate(distance, knots)
            shared = points+normals*offsets[:,None]
            shared[0] = first.point(distances[0],first.offset-q*.5*first.width)
            shared[-1] = last.point(distances[1],last.offset+q*.5*last.width)
            segments = []
            for i,(lo,hi) in enumerate(zip(knots,knots[1:])):
                values = restricted(distance, float(lo), float(hi))
                controls = np.linspace(points[i],points[i+1],7)
                for k in range(7):
                    if k < 6:controls[k] += (6-k)/6.*normals[i]*values[k]
                    if k:controls[k] += k/6.*normals[i+1]*values[k-1]
                controls[0], controls[-1] = shared[i], shared[i+1]
                segments.append(PolynomialCurb(controls))
            curbs.append(PiecewiseCurb(knots,segments))
        metadata = dict(exact_miter=True, polynomial_pieces=len(knots)-1,
                        maximum_miter_reach=reach, endpoint_reconciliation=weld_error,
                        mandatory_knots=knots.tolist())
        if coherent_collinear:
            metadata['coherent_collinear_frames'] = True
        return tuple(curbs), metadata
