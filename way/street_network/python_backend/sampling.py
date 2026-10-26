"""Redistribute compound curb samples on their original analytical curves.

Polygon clipping splits chords at incidental intersections and coincident lead
samples. Those tiny pieces are not physical road features. Retain semantic
ownership changes and full incoming profiles, refine curve intersections, then
sample each retained analytical curb interval without inheriting these pieces.
"""
from dataclasses import replace
import math
import numpy as np
from .compounds import RimSegment, RimSpan, rim_label
from .regions import BoundaryLoop, Region, RegionUnion
from .mesh import cross


class CertifiedCurbSampler:
    """Sample each real curve using its own curvature and physical spacing."""
    def __init__(self, tolerance, spacing, max_samples):
        if (not all(math.isfinite(v) and v > 0. for v in (tolerance, spacing))
                or type(max_samples) is not int or max_samples < 4):
            raise ValueError('Invalid certified curb sampling limits')
        self.tolerance, self.spacing, self.max_samples = tolerance, spacing, max_samples

    def regular(self, first, last, step):
        count = max(1, int(math.ceil((last - first) / step)))
        if count + 1 > self.max_samples:
            raise ValueError('Certified curb sampling exceeds its sample budget')
        return np.linspace(first, last, count + 1)

    def lead_metrics(self, approach, parameters, cumulative, distance, base, length):
        route = approach.path
        if getattr(route, 'outward', False):
            pieces = [(piece.path, float(route.stations[i]), piece.forward) for i, piece in enumerate(route.pieces)]
        else:
            pieces = [(route, 0., approach.end == 0)]
        fixed, curvature, speed = [], 0., 1.
        delta = np.diff(cumulative)
        ratio = float(np.max(np.divide(abs(np.diff(parameters)), delta,
                                      out=np.zeros_like(delta), where=delta > 1e-12)))
        first, last = min(parameters[0], parameters[-1]), max(parameters[0], parameters[-1])
        for path, offset, forward in pieces:
            if offset > last or offset + path.length < max(0., first):
                continue
            if path.exact:
                knots = path.stations if forward else path.length - path.stations
                stations = np.array([float(offset + s) for s in knots if first < offset + s < last])
                xp, fp = (parameters, cumulative) if parameters[-1] >= parameters[0] else (parameters[::-1], cumulative[::-1])
                fixed.extend((base + np.interp(stations, xp, fp)) / length)
                derivative = np.diff(path.normals, axis=0) / path.chords[:, None]
                maximum_speed = 1. + abs(distance) * float(np.linalg.norm(derivative, axis=1).max())
                speed = max(speed, ratio * maximum_speed)
            else:
                k = path.curvature_bound
                curvature = max(curvature, k / (1. - abs(distance) * k))
                speed = max(speed, ratio * (1. + abs(distance) * k))
        return curvature, speed, fixed

    def sample(self, curve, first, last):
        if curve.continuation is not None:
            count = curve.continuation.intervals(self.spacing, self.tolerance)
            parameters = self.regular(first, last, 1. / count)
        elif curve.corner.center is None:
            delta = curve.last - curve.first
            axis = -curve.corner.first.direction
            longitudinal = float(delta @ axis) * axis
            lateral = np.linalg.norm(delta - longitudinal)
            speed = np.linalg.norm(longitudinal) + 1.875 * lateral
            second = 6. * lateral  # >= max |quintic smoothstep second derivative|.
            step = min(self.spacing / max(speed, 1e-12),
                       math.sqrt(8. * self.tolerance / max(second, 1e-12)))
            parameters = self.regular(first, last, step)
        else:
            first_distance, last_distance = first * curve.length, last * curve.length
            curvature, speed, fixed = 0., 1., [first, last]
            for part in range(3):
                lo = max(first_distance, curve.lengths[part])
                hi = min(last_distance, curve.lengths[part + 1])
                if hi <= lo + 1e-10:
                    continue
                if part == 1:
                    curvature = max(curvature, 1. / curve.corner.radius)
                else:
                    approach, parameters, cumulative, offset = curve.leads[0 if part == 0 else 1]
                    k, v, knots = self.lead_metrics(approach, parameters, cumulative, offset,
                                                   curve.lengths[part], curve.length)
                    curvature, speed = max(curvature, k), max(speed, v)
                    fixed.extend(t for t in knots if first < t < last)
            step = self.spacing / speed
            if curvature > 1e-12:
                step = min(step, math.sqrt(8. * self.tolerance / curvature) / speed,
                           math.pi / (2. * curvature * speed))
            fixed = sorted(set(float(t) for t in fixed))
            values = []
            for a, b in zip(fixed, fixed[1:]):
                values.extend(self.regular(a, b, step / curve.length)[:-1])
            parameters = np.r_[values, last]
        # A refined ownership contact can differ from a lead/arc contact by
        # picometres. Retain the authoritative run endpoints and remove only
        # internal numerical duplicates, in physical rather than parameter units.
        for _ in range(4):
            points = curve(parameters)
            duplicates = np.flatnonzero(np.linalg.norm(np.diff(points, axis=0), axis=1) <= 1e-7)
            if not len(duplicates):
                break
            remove = {int(i) if i + 1 == len(parameters) - 1 else int(i + 1) for i in duplicates}
            if 0 in remove or len(parameters) - len(remove) < 2:
                raise ValueError('Retained analytical curb has no positive physical span')
            parameters = np.delete(parameters, sorted(remove))
        else:
            raise ValueError('Analytical curb contacts have unresolved numerical duplicates')
        if (len(parameters) - 1) % 2:
            points = curve(parameters)
            at = int(np.argmax(np.linalg.norm(np.diff(points, axis=0), axis=1)))
            parameters = np.insert(parameters, at + 1, .5 * (parameters[at] + parameters[at + 1]))
        if len(parameters) > self.max_samples:
            raise ValueError('Certified curb sampling exceeds its sample budget')
        return parameters


class AnalyticalBoundarySampler:
    def __init__(self, *, tolerance=.005, spacing=.7, max_samples=16384):
        if (not math.isfinite(tolerance) or tolerance <= 0.
                or not math.isfinite(spacing) or spacing <= 0.
                or type(max_samples) is not int or max_samples < 4):
            raise ValueError('Invalid analytical boundary sampling limits')
        self.tolerance, self.max_samples = tolerance, max_samples
        self.curves = CertifiedCurbSampler(tolerance, spacing, max_samples)

    @staticmethod
    def _circle_contact(curves, initial, lower, upper):
        if any(curve.corner.center is None for curve in curves):
            return None
        for curve, t in zip(curves, initial):
            if not curve.lengths[1] - 1e-7 <= t * curve.length <= curve.lengths[2] + 1e-7:
                return None
        a, b = (curve.corner for curve in curves)
        delta = b.center - a.center
        distance = float(np.linalg.norm(delta))
        candidates = []
        if distance <= 1e-8 and abs(a.radius - b.radius) <= 1e-8:
            # Coincident arcs have no unique Newton root. Keep a nearby common
            # point in both finite arc intervals, including their endpoints.
            midpoint = .5 * (curves[0](initial[0]) + curves[1](initial[1])) - a.center
            if np.linalg.norm(midpoint) > 1e-12:
                candidates.append(a.center + a.radius * midpoint / np.linalg.norm(midpoint))
            for corner in (a, b):
                for angle in (corner.start_angle, corner.start_angle + corner.sweep):
                    candidates.append(corner.center + corner.radius * np.array((math.cos(angle), math.sin(angle))))
        elif distance > 1e-8:
            along = (a.radius ** 2 - b.radius ** 2 + distance ** 2) / (2. * distance)
            height2 = a.radius ** 2 - along ** 2
            if height2 < -1e-8 * max(1., a.radius):
                return None
            axis = delta / distance
            base = a.center + along * axis
            side = math.sqrt(max(0., height2)) * np.array((-axis[1], axis[0]))
            candidates.extend((base + side, base - side))
        solutions = []
        for point in candidates:
            parameters = []
            for curve in curves:
                corner = curve.corner
                if abs(corner.sweep) <= 1e-12:
                    break
                start = np.array((math.cos(corner.start_angle), math.sin(corner.start_angle)))
                vector = point - corner.center
                angle = math.atan2(float(cross(start, vector)), float(start @ vector))
                fractions = [(angle + k * math.tau) / corner.sweep for k in (-1, 0, 1)]
                valid = [u for u in fractions if -1e-8 <= u <= 1. + 1e-8]
                if not valid:
                    break
                u = float(np.clip(valid[0], 0., 1.))
                parameters.append((curve.lengths[1] + u * (curve.lengths[2] - curve.lengths[1])) / curve.length)
            if len(parameters) == 2:
                parameters = np.asarray(parameters)
                if (np.all(parameters >= lower - 1e-9) and np.all(parameters <= upper + 1e-9)
                        and np.linalg.norm(curves[0](parameters[0]) - curves[1](parameters[1])) <= 1e-7):
                    solutions.append((float(np.sum((parameters - initial) ** 2)), parameters))
        return min(solutions, key=lambda item: item[0])[1] if solutions else None

    def _check_contact_budget(self, curves, current, spans, layouts):
        for span, curve, t in zip(spans, curves, current):
            points = curve(np.linspace(0., 1., layouts[span.owner].rows + 1))
            vectors = np.diff(points, axis=0)
            length2 = np.sum(vectors ** 2, axis=1)
            if np.any(length2 <= 1e-16):
                raise ValueError('Analytical contact has a collapsed source curb chord')
            offset = curve(t) - points[:-1]
            fraction = np.clip(np.sum(offset * vectors, axis=1) / length2, 0., 1.)
            distance = np.min(np.linalg.norm(offset - fraction[:, None] * vectors, axis=1))
            if distance > self.tolerance + 1e-7:
                raise ValueError('Analytical curb contact exceeds its source approximation budget')

    def _contact(self, first, last, labels, layouts):
        curves = [layouts[span.owner].boundaries[label.index]
                  for span, label in zip((first, last), labels)]
        initial = np.array((labels[0].last, labels[1].first))
        # Solve locally around the clipping contact. Another intersection of the
        # same finite curves must never be substituted for this one.
        reach = np.array([2. / layouts[span.owner].rows for span in (first, last)])
        lower, upper = np.maximum(0., initial - reach), np.minimum(1., initial + reach)
        circular = self._circle_contact(curves, initial, lower, upper)
        if circular is not None:
            self._check_contact_budget(curves, circular, (first, last), layouts)
            return circular
        current = initial.copy()
        for _ in range(16):
            delta = curves[0](current[0]) - curves[1](current[1])
            error = float(np.linalg.norm(delta))
            if error <= 1e-9:
                # A shallow intersection can move along the curb by more than
                # the chord's normal error. Bound distance to each original
                # sampled curve, rather than freeze its old polygonal contact.
                self._check_contact_budget(curves, current, (first, last), layouts)
                return current
            derivatives = []
            for curve, t in zip(curves, current):
                a, b = max(0., t - 1e-5), min(1., t + 1e-5)
                derivatives.append((curve(b) - curve(a)) / (b - a))
            matrix = np.column_stack((derivatives[0], -derivatives[1]))
            if abs(np.linalg.det(matrix)) <= 1e-10:
                break
            step = np.linalg.solve(matrix, delta)
            for power in range(10):
                candidate = np.clip(current - step * 2. ** -power, lower, upper)
                if np.linalg.norm(curves[0](candidate[0]) - curves[1](candidate[1])) < error:
                    current = candidate
                    break
            else:
                break
        raise ValueError('Analytical curb contact needs a different local intersection branch')

    def build(self, group):
        sources = {p.region.owner: p for p in group.footprints}
        layouts = {owner: source.layout for owner, source in sources.items()}
        if any(layout is None for layout in layouts.values()):
            raise ValueError('Analytical compound sampling needs the source footprint layouts')
        loops, total = [], 0
        for loop in group.loops:
            labels = [rim_label(span, sources) for span in loop.spans]
            contacts = {}
            for at, last in enumerate(loop.spans):
                first = loop.spans[at - 1]
                if first.owner == last.owner:
                    continue
                pair = (labels[at - 1], labels[at])
                if any(label.kind != 'curb' for label in pair):
                    raise ValueError('Compound clipping retained a partial incoming profile')
                contacts[at] = self._contact(first, last, pair, layouts)

            def token(index):
                index %= len(labels)
                return loop.spans[index].owner, labels[index].kind, labels[index].index

            starts = [i for i in range(len(labels)) if token(i) != token(i - 1)]
            if not starts:
                raise ValueError('Closed compound boundary has no finite analytical curb intervals')
            points, spans = [], []
            for index, start in enumerate(starts):
                end = starts[(index + 1) % len(starts)]
                owner, label = loop.spans[start].owner, labels[start]
                if label.kind == 'mouth':
                    indices = range(start, end if end > start else end + len(labels))
                    for at in indices:
                        at %= len(labels)
                        points.append(loop.xy[at])
                        spans.append(RimSpan(owner, labels[at]))
                    total += len(range(start, end if end > start else end + len(labels)))
                else:
                    first, last = label.first, labels[(end - 1) % len(labels)].last
                    if start in contacts:
                        first = float(contacts[start][1])
                    if end in contacts:
                        last = float(contacts[end][0])
                    if last <= first:
                        raise ValueError('Retained analytical curb interval reverses or collapses')
                    layout = layouts[owner]
                    parameters = self.curves.sample(layout.boundaries[label.index], first, last)
                    count = len(parameters) - 1
                    if total + count > self.max_samples:
                        raise ValueError('Analytical compound boundary exceeds its sample budget')
                    sampled = layout.boundaries[label.index](parameters)
                    points.extend(sampled[:-1])
                    spans.extend(RimSpan(owner, RimSegment('curb', label.index, float(a), float(b)))
                                 for a, b in zip(parameters, parameters[1:]))
                    total += count
            if total > self.max_samples:
                raise ValueError('Analytical compound boundary exceeds its sample budget')
            xy = np.asarray(points)
            # Check independently from the radial chart's visibility test.
            checked = xy[::-1] if loop.hole else xy
            RegionUnion().build((Region(group.members[0], checked),))
            xy.flags.writeable = False
            loops.append(BoundaryLoop(xy, tuple(spans)))
        if sum(not loop.hole for loop in loops) != 1 or len(loops) != len(group.loops):
            raise ValueError('Analytical sampling changed compound boundary topology')
        return replace(group, loops=tuple(loops))
