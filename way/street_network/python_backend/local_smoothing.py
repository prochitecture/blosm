"""Finite, tangent cubic bends with straight source runs between them.

All controls here are disposable evaluation data. Source vertices and IDs are
never edited. A bend fits a physical radius, a source-deviation envelope and
the finite neighboring segments; shrinking a radius below its safe minimum
is not an available fallback.
"""
import math
from dataclasses import dataclass
import numpy as np


def segment_distance(points, first, last):
    delta = last - first
    t = np.clip((points - first) @ delta / (delta @ delta), 0., 1.)
    return np.linalg.norm(points - first - t[..., None] * delta, axis=-1)


def bezier_segment_envelope(controls, first, last, limit, depth=3):
    """A bounded de Casteljau proof, including finite segment endpoints."""
    if float(segment_distance(controls, first, last).max()) <= limit + 1e-10:
        return True
    if depth == 0:
        return False
    levels = [controls]
    for unused in range(3):
        levels.append(.5 * (levels[-1][:-1] + levels[-1][1:]))
    left = np.array([v[0] for v in levels])
    right = np.array([v[-1] for v in levels[::-1]])
    return (bezier_segment_envelope(left, first, last, limit, depth - 1)
            and bezier_segment_envelope(right, first, last, limit, depth - 1))


def within_source_envelope(coefficients, points, limit):
    """Convex distance to a finite segment bounds its complete Bezier hull."""
    p, a, b, c = np.moveaxis(coefficients, 1, 0)
    controls = np.stack((p, p + a / 3., p + (2. * a + b) / 3., p + a + b + c), axis=1)
    return all(bezier_segment_envelope(q, first, last, limit)
               for q, first, last in zip(controls, points, points[1:]))


@dataclass(frozen=True)
class SharedEndpointBend:
    """One directed half of a single bend across an explicit source seam."""
    controls: object
    span: float
    setback: float
    radius: float


def shared_endpoint_bends(first, first_end, second, second_end, limit):
    """Round a sparse compatible seam without interpolating its sharp vertex.

    Callers use a canonical incidence order, so both paths receive the exact
    same midpoint and tangent. The source records themselves remain untouched.
    Gentle span-wide curves retain their established construction.
    """
    a = np.asarray(first.points, dtype=float)[:, :2]
    b = np.asarray(second.points, dtype=float)[:, :2]
    corner = a[0 if first_end == 0 else -1]
    other_corner = b[0 if second_end == 0 else -1]
    if np.linalg.norm(corner - other_corner) > 1e-6:
        return None  # The existing incidence validator reports disconnected copies.
    before = a[1 if first_end == 0 else -2]
    after = b[1 if second_end == 0 else -2]
    lengths = np.array((np.linalg.norm(corner - before), np.linalg.norm(after - corner)))
    u, v = (corner - before) / lengths[0], (after - corner) / lengths[1]
    turn = math.atan2(float(u[0] * v[1] - u[1] * v[0]), float(u @ v))
    # The old normalized endpoint handle's transverse bow is 4 L sin(a)/27.
    # Localize both sides together if either adjacent span exceeds its envelope.
    if 4. / 27. * lengths.max() * math.sin(abs(turn) * .5) <= limit:
        return None
    if abs(turn) >= math.radians(170.):
        raise ValueError('Shared Smooth bend has a near-reversal')
    reach = max(abs(s.offset) + .5 * s.width for s in (first, second))
    factor = math.tan(abs(turn) * .5)
    radius = min(2. * reach, .99 * limit / (1. - math.cos(abs(turn) * .5)),
                 .4 * lengths.min() / factor)
    if radius < 1.5 * reach - 1e-10:
        raise ValueError('Shared Smooth bend needs more source length for the carriageway width')
    setback = radius * factor
    start, finish = corner - setback * u, corner + setback * v
    center = start + math.copysign(radius, turn) * np.array((-u[1], u[0]))
    def piece(p, q, t, w, angle):
        handle = 4. / 3. * radius * math.tan(abs(angle) * .25)
        return np.array((p, p + handle * t, q - handle * w, q)), 3. * handle
    if abs(turn) <= math.pi * .5 + 1e-12:
        controls, span = piece(start, finish, u, v, turn)
        levels = [controls]
        for unused in range(3):
            levels.append(.5 * (levels[-1][:-1] + levels[-1][1:]))
        left = np.array([level[0] for level in levels])
        right = np.array([level[-1] for level in levels[::-1]])
        span *= .5
    else:
        angle = turn * .5
        rotation = np.array(((math.cos(angle), -math.sin(angle)),
                             (math.sin(angle), math.cos(angle))))
        midpoint, tangent = center + rotation @ (start - center), rotation @ u
        left, span = piece(start, midpoint, u, tangent, angle)
        right, unused = piece(midpoint, finish, tangent, v, angle)
    for controls, p, q in ((left, before, corner), (right, corner, after)):
        if not bezier_segment_envelope(controls, p, q, limit):
            raise ValueError('Shared Smooth bend exceeds the source-deviation limit')
    # Canonical traversal enters first's end, then exits second's start.
    return (SharedEndpointBend(left if first_end == 1 else left[::-1], span, setback, radius),
            SharedEndpointBend(right if second_end == 0 else right[::-1], span, setback, radius))


class LocalCubicSmoothing:
    """Width-relative default policy, independent of long source chord lengths."""
    def __init__(self, points, reach, deviation_limit, endpoint_directions=None, endpoint_bends=None):
        self.source = np.asarray(points, dtype=float)
        self.reach = float(reach)
        self.deviation_limit = float(deviation_limit)
        delta = np.diff(self.source, axis=0)
        self.lengths = np.linalg.norm(delta, axis=1)
        self.directions = delta / self.lengths[:, None]
        self.endpoints = endpoint_directions or {}
        self.endpoint_bends = endpoint_bends or {}

    def build(self):
        directions, lengths = self.directions, self.lengths
        turns = np.arctan2(directions[:-1, 0] * directions[1:, 1]
                          - directions[:-1, 1] * directions[1:, 0],
                          np.sum(directions[:-1] * directions[1:], axis=1))
        if np.any(abs(turns) >= math.radians(170.)):
            raise ValueError('Local Smooth bend has a near-reversal')
        active = abs(turns) > 1e-7
        factor = np.tan(abs(turns) * .5)
        # Reserve finite inner-curb reach before the independent certificate.
        minimum = np.where(active, 1.5 * self.reach * factor, 0.)
        preferred = np.where(active, 2. * self.reach * factor, 0.)
        deviation_radius = .99 * self.deviation_limit / np.maximum(1. - np.cos(abs(turns) * .5), 1e-30)
        preferred = np.minimum(preferred, deviation_radius * factor)
        if np.any(active & (preferred < minimum - 1e-10)):
            raise ValueError('Local Smooth bend cannot fit the source-deviation limit at its safe radius')
        # The whole bend lies inside its trim triangle. Bounding its distance
        # from the corner vertex also bounds distance from the source polyline.
        # The tighter half-Bezier hull check below certifies actual deviation.
        setback = np.minimum(preferred, .8 * np.minimum(lengths[:-1], lengths[1:]))
        setback = np.maximum(setback, minimum)
        endpoint_spans = np.zeros(2)
        endpoint_tangents = [directions[0], directions[-1]]
        for end, value in self.endpoints.items():
            if end in self.endpoint_bends:
                continue
            tangent = np.asarray(value, dtype=float)
            tangent = tangent / np.linalg.norm(tangent)
            endpoint_tangents[end] = tangent
            direction = directions[0 if end == 0 else -1]
            sine = abs(float(tangent[0] * direction[1] - tangent[1] * direction[0]))
            if float(tangent @ direction) <= 0.:
                raise ValueError('Local Smooth endpoint tangent turns back')
            if sine > 1e-8:
                # For this endpoint Hermite cubic the transverse displacement
                # is D sin(angle) u(1-u)^2, whose maximum is 4D sin(angle)/27.
                # Fit the preferred extent to that bound before constructing
                # and independently certifying the actual Bezier hull.
                endpoint_spans[end] = max(.02 * self.reach,
                    min(12. * self.reach * sine, 6.6 * self.deviation_limit / sine))
        # Neighboring bends share the finite segment budget. Only the excess
        # over independently required minimum radii may be reduced.
        for i, length in enumerate(lengths):
            ids = [j for j in (i - 1, i) if 0 <= j < len(setback)]
            required = float(minimum[ids].sum())
            endpoint_ids = ([0] if i == 0 else []) + ([1] if i == len(lengths) - 1 else [])
            endpoint = float(endpoint_spans[endpoint_ids].sum())
            fixed = sum(self.endpoint_bends[end].setback for end in endpoint_ids if end in self.endpoint_bends)
            available = .9 * length - required - fixed
            if endpoint > max(0., available):
                endpoint_spans[endpoint_ids] *= max(0., available) / endpoint
                endpoint = float(endpoint_spans[endpoint_ids].sum())
            budget = .9 * length - endpoint - fixed
            if required > budget + 1e-10:
                raise ValueError('Local Smooth bends need more source length for the carriageway width')
            excess = float((setback[ids] - minimum[ids]).sum())
            if excess > 0. and float(setback[ids].sum()) > budget:
                setback[ids] = minimum[ids] + (setback[ids] - minimum[ids]) * max(0., (budget - required) / excess)
        pieces = []
        def add(controls, span=None):
            controls = np.asarray(controls)
            if np.linalg.norm(controls[-1] - controls[0]) < 1e-9:
                return
            pieces.append((controls, float(span if span is not None else np.linalg.norm(controls[-1] - controls[0]))))
        def line(first, last):
            add(np.linspace(first, last, 4))
        cursor = self.source[0]
        if 0 in self.endpoint_bends:
            bend = self.endpoint_bends[0]
            add(bend.controls, bend.span)
            cursor = bend.controls[-1]
        elif endpoint_spans[0]:
            target = cursor + endpoint_spans[0] * directions[0]
            span = endpoint_spans[0]
            controls = np.array((cursor, cursor + span / 3. * endpoint_tangents[0],
                                 target - span / 3. * directions[0], target))
            if not bezier_segment_envelope(controls, cursor, target, self.deviation_limit):
                raise ValueError('Local Smooth endpoint exceeds the source-deviation limit')
            add(controls, span)
            cursor = target
        bends = []
        for i, (turn, distance) in enumerate(zip(turns, setback)):
            corner = self.source[i + 1]
            if not active[i]:
                line(cursor, corner)
                cursor = corner
                continue
            first = corner - distance * directions[i]
            last = corner + distance * directions[i + 1]
            line(cursor, first)
            radius = distance / factor[i]
            sign = math.copysign(1., turn)
            center = first + sign * radius * np.array((-directions[i, 1], directions[i, 0]))
            count = max(1, int(math.ceil((abs(turn) - 1e-12) / (math.pi * .5))))
            step = turn / count
            handle = 4. / 3. * radius * math.tan(abs(step) * .25)
            radial = first - center
            tangent = directions[i]
            for part in range(count):
                angle = step * (part + 1)
                rotation = np.array(((math.cos(angle), -math.sin(angle)),
                                     (math.sin(angle), math.cos(angle))))
                end_radial = rotation @ radial
                end_tangent = rotation @ directions[i]
                a = first if part == 0 else pieces[-1][0][-1]
                b = last if part == count - 1 else center + end_radial
                controls = np.array((a, a + handle * tangent, b - handle * end_tangent, b))
                # Restrict at 1/2 before testing the two finite source legs.
                levels = [controls]
                for unused in range(3):
                    levels.append(.5 * (levels[-1][:-1] + levels[-1][1:]))
                halves = (np.array([v[0] for v in levels]), np.array([v[-1] for v in levels[::-1]]))
                for half in halves:
                    bounds = [float(segment_distance(half, a0, b0).max()) for a0, b0 in
                              ((self.source[i], corner), (corner, self.source[i + 2]))]
                    if min(bounds) > self.deviation_limit + 1e-9:
                        raise ValueError('Local Smooth bend exceeds the source-deviation limit')
                add(controls, 3. * handle)
                tangent = end_tangent
            bends.append(dict(vertex=i + 1, radius=float(radius), setback=float(distance)))
            cursor = last
        end = (self.endpoint_bends[1].controls[0] if 1 in self.endpoint_bends else
               self.source[-1] - endpoint_spans[1] * directions[-1])
        line(cursor, end)
        if 1 in self.endpoint_bends:
            bend = self.endpoint_bends[1]
            add(bend.controls, bend.span)
        elif endpoint_spans[1]:
            span = endpoint_spans[1]
            controls = np.array((end, end + span / 3. * directions[-1],
                                 self.source[-1] - span / 3. * endpoint_tangents[1], self.source[-1]))
            if not bezier_segment_envelope(controls, end, self.source[-1], self.deviation_limit):
                raise ValueError('Local Smooth endpoint exceeds the source-deviation limit')
            add(controls, span)
        controls = np.asarray([p[0] for p in pieces])
        points = np.concatenate((controls[:, 0], controls[-1:, 3]))
        coefficients = np.stack((controls[:, 0], 3. * (controls[:, 1] - controls[:, 0]),
            3. * (controls[:, 2] - 2. * controls[:, 1] + controls[:, 0]),
            controls[:, 3] - 3. * controls[:, 2] + 3. * controls[:, 1] - controls[:, 0]), axis=1)
        return points, np.array([p[1] for p in pieces]), coefficients, bends
