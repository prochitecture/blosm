"""Simple source fillet patches with internal contact cross-sections.

Positive contacts close through the original centreline inside its road strip.
Negative auxiliary contacts retain the exact offset lead back to the original
cross-section. Radial node-to-arc spokes are never physical patch seams.
"""
from types import SimpleNamespace
import math
import numpy as np
from .compounds import RimSegment, SampledFootprint
from .regions import Region, RegionUnion, signed_area


class _PreparedCurve:
    continuation = None

    def __init__(self, points):
        self.points = np.asarray(points, dtype=float)
        self.contact_reference_points = self.points
        self.stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(self.points, axis=0), axis=1))]
        self.length = float(self.stations[-1])
        self.corner = SimpleNamespace(center=None)
        if self.length <= 1e-8:
            raise ValueError('Corner patch has a collapsed prepared curve')

    def __call__(self, t):
        station = np.asarray(t)*self.length
        return np.stack([np.interp(station, self.stations, self.points[:, k]) for k in (0, 1)], axis=-1)

    def certified_parameters(self, sampler, first, last):
        knots = sorted({first, last, *(float(s/self.length) for s in self.stations
                                     if first*self.length < s < last*self.length)})
        values = []
        for a, b in zip(knots, knots[1:]):
            values.extend(sampler.regular(a, b, sampler.spacing/self.length)[:-1])
        return np.r_[values, last]


def _stations(approach, station, spacing, tolerance, maximum):
    if station <= 1e-7:
        return np.array((0., station))
    path = approach.path
    if station > approach.length+1e-7:
        raise ValueError('Corner patch contact exceeds its original source approach')
    if not getattr(path, 'outward', False):
        first, last = ((0., station) if approach.end == 0 else (path.length-station, path.length))
        values = path.samples(first, last, spacing, tolerance=tolerance, max_samples=maximum)
        return values if approach.end == 0 else path.length-values[::-1]
    values = [0., station]
    for i, piece in enumerate(path.pieces):
        begin = float(path.stations[i])
        extent = min(station, float(path.stations[i+1]))-begin
        if extent <= 1e-7:
            break
        first, last = ((0., extent) if piece.forward else (piece.path.length-extent, piece.path.length))
        sampled = piece.path.samples(first, last, spacing, tolerance=tolerance, max_samples=maximum)
        if not piece.forward:
            sampled = piece.path.length-sampled[::-1]
        values.extend(begin+sampled)
        if len(values) > maximum:
            raise ValueError('Corner patch exceeds its sampling budget')
    return np.array(sorted(set(values)))


def corner_patch(owner, corner, spacing, tolerance, *, max_samples=8192):
    if (not all(math.isfinite(value) and value > 0. for value in (spacing, tolerance))
            or type(max_samples) is not int or max_samples < 4 or corner.center is None):
        raise ValueError('Invalid circular corner patch sampling request')
    if np.linalg.norm(corner.first.origin-corner.second.origin) > 1e-7:
        raise ValueError('Local source corner patch requires one original source node')
    points, kinds = [np.asarray(corner.first.origin)], []
    arc_phases,lead_tags=[np.nan],[None]

    def append(chain, kind, *, phases=None, lead_role=None, lead_stations=None):
        for index,point in enumerate(chain):
            phase=np.nan if phases is None else float(phases[index])
            lead=None if lead_role is None else (lead_role,float(lead_stations[index]))
            if np.linalg.norm(point-points[-1]) > 1e-8:
                points.append(np.asarray(point))
                kinds.append(kind)
                arc_phases.append(phase);lead_tags.append(lead)
                if len(points) > max_samples+1:
                    raise ValueError('Corner patch exceeds its sampling budget')
            else:
                # Retain the producer's semantic identity on the vertex
                # already kept by the unchanged duplicate handling.
                if phases is not None:arc_phases[-1]=phase
                if lead is not None:lead_tags[-1]=lead

    first, last = corner.stations
    at = _stations(corner.first, first, spacing, tolerance, max_samples)
    if first >= 0.:
        append(corner.first.frame(at)[0], 'internal')
    else:
        lead = corner.first.point(at, corner.first.offset-.5*corner.first.width)
        append(lead[:1], 'internal',lead_role=0,lead_stations=at[:1])
        append(lead[1:], 'curb',lead_role=0,lead_stations=at[1:])
    count = max(2, int(math.ceil(abs(corner.sweep)*corner.radius/min(
        spacing, math.sqrt(8.*tolerance*corner.radius)))))
    if count+1 > max_samples:
        raise ValueError('Corner patch exceeds its sampling budget')
    angle = np.linspace(corner.start_angle, corner.start_angle+corner.sweep, count+1)
    arc = corner.center+corner.radius*np.c_[np.cos(angle), np.sin(angle)]
    phases=np.linspace(0.,1.,count+1)
    append(arc[:1], 'internal' if first >= 0. else 'curb',phases=phases[:1])
    append(arc[1:], 'curb',phases=phases[1:])
    at = _stations(corner.second, last, spacing, tolerance, max_samples)[::-1]
    if last >= 0.:
        append(corner.second.frame(at)[0], 'internal')
    else:
        lead = corner.second.point(at, corner.second.offset+.5*corner.second.width)
        append(lead, 'curb',lead_role=1,lead_stations=at)
        append((corner.second.origin,), 'internal')
    if np.linalg.norm(points[-1]-points[0]) > 1e-8:
        append((points[0],), 'internal')
    points = np.asarray(points[:-1])
    arc_phases=np.asarray(arc_phases[:-1]);lead_tags=lead_tags[:-1]
    if len(kinds) != len(points):
        raise ValueError('Corner patch has unresolved closing incidence')
    if signed_area(points) < 0.:
        points = points[::-1]
        kinds = list(np.roll(np.asarray(kinds)[::-1], -1))
        arc_phases=arc_phases[::-1];lead_tags=lead_tags[::-1]
    # Start at a semantic change so no prepared curve wraps its parameters.
    start = next((i for i in range(len(points)) if kinds[i] != kinds[i-1]), 0)
    points = np.roll(points, -start, axis=0)
    kinds = kinds[start:]+kinds[:start]
    arc_phases=np.roll(arc_phases,-start);lead_tags=lead_tags[start:]+lead_tags[:start]
    region = Region(owner, points)
    RegionUnion().build((region,))
    curves, labels = [], []
    first = 0
    while first < len(points):
        last = first+1
        while last < len(points) and kinds[last] == kinds[first]:
            last += 1
        curve = _PreparedCurve(np.vstack((points[first:last], points[last % len(points)])))
        if kinds[first]=='curb':
            from .source_circular_intervals import declare
            indices=np.r_[np.arange(first,last),last%len(points)]
            declaration=declare(corner,curve,arc_phases[indices],
                [lead_tags[index] for index in indices],tolerance)
            if declaration is not None:curve.original_circle_interval=declaration
        index = len(curves)
        curves.append(curve)
        labels.extend(RimSegment(kinds[first], index, float(a/curve.length), float(b/curve.length))
                      for a, b in zip(curve.stations, curve.stations[1:]))
        first = last
    layout = SimpleNamespace(approaches=(), boundaries=tuple(curves), rows=len(points), mouths={}, corners=(corner,))
    return SampledFootprint(region, tuple(labels), layout)
