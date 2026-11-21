"""Finite original-profile and declared Bundle references for opt-in heights.

References choose no source geometry or incidence. They retain the canonical
Path/Profile/Bundle objects and reject evaluations outside their owned metric
domains. The caller still owns mesh, drainage and publication qualification.
"""
from dataclasses import dataclass
import math
import numpy as np
from .model import Key
from .mesh import Mesh
from .surface_query import QuadHeightField
from .profiles import memberships, ProfileProvider
from .bundle_axis import AxialChartUnavailable


class SourceBodyUnavailable(ValueError):
    def __init__(self, message, *, scope=None):
        self.scope = dict(scope or {}, reason=message)
        super().__init__(message)


class NoPairedBodyDomain(SourceBodyUnavailable):
    """Finite declared ranges are disjoint; no paired field is requested."""


def normal_root_distance_bound(points, path, low, high):
    """Bound distance to the entire original arc-length bracket.

    A unit-speed curve with ||C''|| <= k differs from its endpoint chord by
    at most k*(high-low)^2/8. Norm convexity bounds distance to that chord by
    the larger endpoint distance. The bracket and global curvature bound are
    unchanged; this is tighter than the old endpoint-distance-plus-length
    estimate, not permission to accept a nonmonotone normal root.
    """
    first = path.frame(low)[0]
    last = path.frame(high)[0]
    distance = np.maximum(np.linalg.norm(points-first, axis=1),
                          np.linalg.norm(points-last, axis=1))
    return distance+path.curvature_bound*(high-low)**2/8.


@dataclass(frozen=True)
class BodyGuideLimits:
    maximum_samples: int = 20000
    maximum_cells: int = 8192
    maximum_pairs: int = 2000000
    maximum_field_pairs: int = 256
    maximum_classification_work: int = 100000000
    maximum_polygon_vertices: int = 20000
    maximum_mesh_cells: int = 50000
    source_tolerance: float = .00025
    maximum_source_error: float = .002

    def __post_init__(self):
        for name, ceiling in (('maximum_samples', 20000), ('maximum_cells', 8192),
                ('maximum_pairs', 2000000), ('maximum_field_pairs', 256),
                ('maximum_classification_work', 100000000),
                ('maximum_polygon_vertices', 20000), ('maximum_mesh_cells', 50000)):
            value = getattr(self, name)
            if type(value) is not int or not 1 <= value <= ceiling:
                raise ValueError('Source-body work limits must remain positive and bounded')
        if (not math.isfinite(self.source_tolerance) or not 0. < self.source_tolerance <= .00025
                or not math.isfinite(self.maximum_source_error)
                or not 0. < self.maximum_source_error <= .002):
            raise ValueError('Source-body references require unchanged finite source tolerances')


class BodyWork:
    """Shared pre-allocation/source work allowance for one helper invocation."""
    def __init__(self, limits):
        self.limits = limits
        self.used = dict(samples=0, cells=0, pairs=0, field_pairs=0, classification_work=0)

    def charge(self, name, count):
        if type(count) is not int or count < 0:
            raise ValueError('Source-body work must use nonnegative integer counts')
        bound = getattr(self.limits, 'maximum_'+name)
        if self.used[name]+count > bound:
            raise SourceBodyUnavailable('Source-body '+name+' budget exhausted', scope=self.used)
        self.used[name] += count

    def samples(self, path, first, last, spacing):
        remaining = self.limits.maximum_samples-self.used['samples']
        if remaining < 3:
            raise SourceBodyUnavailable('Source-body samples budget exhausted', scope=self.used)
        values = path.samples(first, last, spacing, tolerance=self.limits.source_tolerance,
                              max_samples=remaining)
        self.charge('samples', len(values))
        return values


@dataclass(frozen=True, eq=False)
class OriginalBodyReference:
    path: object
    profile: object
    first: float
    last: float
    stations: object

    @classmethod
    def build(cls, context, section_id, first, last, work):
        section = context.network.section_index[section_id]
        path, profile = context[Key('path', section_id)], context[Key('profile', section_id)]
        expected = ProfileProvider.build(section, memberships(context.network).get(section_id),
                                         context.network.section_index)
        if path.section is not section or profile != expected:
            raise SourceBodyUnavailable('Source-body reference changed its original Path or Profile')
        if (not all(map(math.isfinite, (first, last))) or not 0. <= first < last <= path.length):
            raise ValueError('Source-body reference must use a positive finite owned original range')
        at = work.samples(path, first, last, min(.5, context.network.settings.spacing)).copy()
        at.flags.writeable = False
        return cls(path, profile, float(first), float(last), at)

    @property
    def section_id(self):
        return self.path.section.id

    @property
    def dependencies(self):
        return (Key('path', self.section_id), Key('profile', self.section_id))

    def inverse(self, points, *, max_cells=8192):
        """Invert only this finite original strip; no nearest-source choice."""
        points = np.asarray(points, float)
        if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
            raise ValueError('Original body coordinates require finite XY rows')
        at, path = self.stations, self.path
        count = len(at)
        left = path.point(at, path.section.offset-.5*path.section.width)
        right = path.point(at, path.section.offset+.5*path.section.width)
        xyz = np.c_[np.vstack((left, right)), np.r_[at, at]]
        faces = np.array([(count+i, count+i+1, i+1, i) for i in range(count-1)], np.int32)
        query = QuadHeightField(Mesh(xyz, faces), datum=self.first, ceiling=self.last,
                                max_cells=max_cells, tolerance=1e-10)
        station = query.evaluate(points, default=-1.)
        accepted = station >= self.first
        ids = np.flatnonzero(accepted)
        proof = dict(original_section=self.section_id, interval=(self.first, self.last),
                     cells=len(faces), maximum_inverse_residual=0.)
        if not path.exact and len(ids):
            cells = np.clip(np.searchsorted(at, station[ids], side='right')-1, 0, len(at)-2)
            low, high = at[np.maximum(0, cells-1)], at[np.minimum(len(at)-1, cells+2)]
            distance = np.linalg.norm(points[ids]-path.frame(low)[0], axis=1)+(high-low)
            jacobian = 1.-path.curvature_bound*distance
            rejected = jacobian <= 0.
            if rejected.any():
                # Keep every already-positive old certificate and its output
                # report. Only a failed coarse bound uses the tighter bound
                # over the same complete bracket and same global curvature.
                tight = normal_root_distance_bound(points[ids[rejected]], path, low[rejected], high[rejected])
                jacobian[rejected] = 1.-path.curvature_bound*tight
                proof['normal_root_distance_bound'] = 'same_bracket_arc_length_chord_curvature'
            if np.min(jacobian) <= 0.:
                raise SourceBodyUnavailable('Original body inverse lacks a monotone normal-root certificate')
            f0 = np.sum((points[ids]-path.frame(low)[0])*path.frame(low)[1], axis=1)
            f1 = np.sum((points[ids]-path.frame(high)[0])*path.frame(high)[1], axis=1)
            valid = (f0 >= -1e-8) & (f1 <= 1e-8)
            accepted[ids[~valid]] = False
            ids, low, high = ids[valid], low[valid], high[valid]
            for _ in range(56):
                middle = .5*(low+high)
                center, tangent, _ = path.frame(middle)
                positive = np.sum((points[ids]-center)*tangent, axis=1) > 0.
                low, high = np.where(positive, middle, low), np.where(positive, high, middle)
            station[ids] = .5*(low+high)
            proof['minimum_normal_root_jacobian'] = float(jacobian.min())
        q = np.full(len(points), np.nan)
        if len(ids):
            center, _, normal = path.frame(station[ids])
            offset = np.sum((points[ids]-center)*normal, axis=1)/np.sum(normal*normal, axis=1)
            q[ids] = 2.*(offset-path.section.offset)/path.section.width
            error = np.linalg.norm(path.point(station[ids], path.section.offset+.5*path.section.width*q[ids])
                                   -points[ids], axis=1)
            proof['maximum_inverse_residual'] = float(error.max())
            if error.max() > 1e-7:
                raise SourceBodyUnavailable('Original body inverse failed its physical source residual')
            accepted[ids[abs(q[ids]) > 1.+1e-10]] = False
        return station, q, accepted, proof

    def values(self, points):
        station, q, accepted, _ = self.inverse(points)
        if not accepted.all() or np.any(station < self.first) or np.any(station > self.last):
            raise SourceBodyUnavailable('Original body field would extrapolate its finite source range')
        return self.profile.datum+self.profile.value(q)


def _linear_progress(frame, path, forward):
    from .bundles import CenterlineRange
    from .bundle_boundaries import _certify_increasing
    direction = 1. if forward else -1.
    if path.exact:
        lower = min(float(np.min(np.diff(path.point(path.knots,
            path.section.offset+side*.5*path.section.width), axis=0)@frame.axis*direction))
            for side in (-1., 1.))
        if lower <= 1e-9:
            raise AxialChartUnavailable('Original offset curbs have no positive linear progress')
        return lower
    reach = abs(path.section.offset)+.5*path.section.width
    if reach*path.curvature_bound >= 1.:
        raise AxialChartUnavailable('Original source offsets lack regular linear coordinates')
    return min(_certify_increasing(3.*np.diff(c, axis=0)@frame.axis*direction,
               coordinate='source-body linear') for c in CenterlineRange.build(path, 0., path.length).curves)


def declared_frame(context, bundle):
    """Established angular coordinates first, one source-directed linear retry."""
    from .bundle_angular_internal import shared_frame
    from .finite_terminal_entrances import LinearFrame
    paths = [(context[Key('path', v.section_id)], v.forward)
             for member in bundle.members for v in member.path]
    try:
        frame = shared_frame(context, bundle)
        lower = min(frame.certify(path, forward=forward) for path, forward in paths)
        return frame, dict(family='existing_shared_angular', progress_lower=lower)
    except AxialChartUnavailable as error:
        angular_failure = str(error)
    path, forward = paths[0]
    center, tangent, _ = path.frame(.5*path.length)
    frame = LinearFrame(center, tangent*(1. if forward else -1.))
    lower = min(_linear_progress(frame, path, forward) for path, forward in paths)
    return frame, dict(family='source_directed_linear', progress_lower=lower,
                       established_angular_failure=angular_failure)


@dataclass(frozen=True, eq=False)
class PairedBodyReference:
    bundle: object
    originals: tuple
    forward: tuple
    frame: object
    field: object
    first: float
    last: float
    frame_proof: object

    @classmethod
    def build(cls, context, bundle, originals, frame, frame_proof):
        if not any(bundle is value for value in context.network.bundles) or len(originals) != 2:
            raise SourceBodyUnavailable('Paired body reference needs its actual declared Bundle identity')
        forward = []
        for role, source in enumerate(originals):
            visits = [v for v in bundle.members[role].path if v.section_id == source.section_id]
            if len(visits) != 1 or source.path is not context[Key('path', source.section_id)]:
                raise SourceBodyUnavailable('Paired body roles must retain unique original Bundle visits')
            forward.append(visits[0].forward)
            if source.profile.oriented(int(not visits[0].forward)).half != (1 if role == 0 else -1):
                raise SourceBodyUnavailable('Paired body reference lost its original facing half-profile')
        a, b = originals
        if abs(a.profile.datum-b.profile.datum) > 1e-10 or abs(a.profile.rise-b.profile.rise) > 1e-10:
            raise SourceBodyUnavailable('Paired body reference requires its original common datum and rise')
        curves = []
        for source, direct in zip(originals, forward):
            at = source.stations if direct else source.stations[::-1]
            sign = 1. if direct else -1.
            curves.extend(source.path.point(at, source.path.section.offset+side*sign*.5*source.path.section.width)
                          for side in (-1., 1.))
        from .bundle_angular_curves import AngularCrown
        from .finite_terminal_entrances import LinearCrown
        field = (LinearCrown if hasattr(frame, 'axis') else AngularCrown)(
            frame, curves, datum=a.profile.datum, rise=a.profile.rise)
        first = max(float(curve.stations[0]) for curve in field.curves)
        last = min(float(curve.stations[-1]) for curve in field.curves)
        if first >= last:
            raise NoPairedBodyDomain('Paired body fields have no complete finite common fibre range')
        return cls(bundle, tuple(originals), tuple(forward), frame, field, first, last, dict(frame_proof))

    @property
    def dependencies(self):
        return tuple(sorted({self.bundle.key, *(key for member in self.bundle.members
            for visit in member.path for key in (Key('path', visit.section_id), Key('profile', visit.section_id)))}))

    def values(self, points):
        station = self.frame.local(points)[..., 0]
        if np.any(station < self.first) or np.any(station > self.last):
            raise SourceBodyUnavailable('Paired body field would extrapolate an original curb domain')
        return self.field.values(points)

    def whole_body_intervals(self, core_ranges):
        """Intersect the complete, influence-free fibre intervals of both banks."""
        common = [(self.first, self.last)]
        for source, forward in zip(self.originals, self.forward):
            ranges = core_ranges[source]
            offsets = source.path.section.offset+.5*source.path.section.width*np.array((-1., 1.))
            intersections = []
            index = 0
            for interval in ranges if forward else reversed(ranges):
                first, last = interval if forward else interval[::-1]
                low = float(self.frame.local(source.path.point(first, offsets))[:, 0].max())
                high = float(self.frame.local(source.path.point(last, offsets))[:, 0].min())
                if high <= low:
                    continue
                while index < len(common) and common[index][1] <= low:
                    index += 1
                for i in range(index, len(common)):
                    a, b = common[i]
                    if a >= high:
                        break
                    intersections.append((max(a, low), min(b, high)))
            common = intersections
        return tuple(common)
