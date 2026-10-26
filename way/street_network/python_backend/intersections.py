"""Rounded corner bands and a directly constructed cyclic Coons quad core.

The topology follows the canonical cyclic patterns used during GN development:
each neighboring pair of half-profiles owns one corner band; every arm owns one
four-sided central block. Numerical curves replace GN field evaluation. Actual
geometry must pass the same local Jacobian guards before it can be published.
"""
from dataclasses import dataclass, replace
from collections import OrderedDict
import math
import numpy as np
from .mesh import Assembler, cross
from .model import End, Key
from .paths import right
from .providers import Part
from .updates import FeatureProvider, Task
from .routes import Route, RouteTopology, route_key
from .continuations import WidthContinuation
from .corner_radii import CornerContactError, CornerRadiusPolicy, fitted_radius_report


class CornerConstructionError(ValueError):
    """A failed circular contact with the actual finite source approaches."""
    def __init__(self, message, junction, approaches):
        super().__init__(message)
        self.junction = junction
        self.approaches = tuple(approaches)


class ApproachExtentError(ValueError):
    def __init__(self, junction, extents):
        super().__init__('Requested intersection footprint consumes its complete compatible approach route')
        self.junction = junction
        self.approaches = tuple(a for a, station in extents)
        self.required_stations = {End(a.source.id, a.end): float(station) for a, station in extents}


def smooth(t):
    return t * t * t * (10. + t * (-15. + 6. * t))


def half_intervals(settings):
    # A shared transverse parameter partition keeps both ends of every road
    # compatible without propagating a width edit through the whole network.
    # Physical edge-ratio guards still reject unsuitable dimensional extremes.
    return max(4, 4 * int(math.ceil(3. / (4. * settings.spacing))))


class Approach:
    frame_cache_limit = 16

    def __init__(self, path, end, profile):
        self.path, self.end = path, end
        self.source = path.section
        self.length = path.length
        self.offset = path.section.offset * (1. if end == 0 else -1.)
        self.width = path.section.width
        self.profile = profile.oriented(end)
        mouth = getattr(path, 'mouth_piece', None)
        if mouth is not None:
            # The incidence keeps its original identity. Its external port is
            # on a later, explicitly connected Section wholly outside the site.
            self.width = mouth.section.width
            self.offset = mouth.section.offset * (1. if mouth.forward else -1.)
            self.profile = path.mouth_profile
        self._frames = OrderedDict()
        self.origin, self.direction, _ = self.frame(0.)

    def frame(self, station):
        station = np.asarray(station)
        # Mouth conditions repeatedly query the same station and its two
        # finite-difference neighbors. Cache scalar frames only; vectorized
        # sampling and Exact miter arithmetic keep their established path.
        key = ((station.dtype.str, station.tobytes()) if station.ndim == 0
               and station.dtype.kind in 'fiu' and np.isfinite(station) else None)
        if key is not None and key in self._frames:
            self._frames.move_to_end(key)
            return tuple(v.copy() for v in self._frames[key])
        clipped = np.clip(station, 0., self.path.length)
        if getattr(self.path, 'outward', False):
            p, tangent, normal = self.path.frame(clipped)
        else:
            global_station = clipped if self.end == 0 else self.path.length - clipped
            p, tangent, normal = self.path.frame(global_station)
            if self.end == 1:
                tangent, normal = -tangent, -normal
        # A bounded tangent extension is auxiliary corner construction inside
        # the intersection. It changes none of the authored centerline controls.
        p = p + np.minimum(station, 0.)[..., None] * tangent
        if key is not None:
            self._frames[key] = tuple(v.copy() for v in (p, tangent, normal))
            if len(self._frames) > self.frame_cache_limit:
                self._frames.popitem(last=False)
        return p, tangent, normal

    def point(self, station, distance):
        point, _, normal = self.frame(station)
        return point + np.asarray(distance)[..., None] * normal

    def jet(self, station, distance):
        point = self.point(station, distance)
        step = 1e-4
        # Offset derivatives include Exact miter interpolation; assuming a
        # unit normal derivative here would change Exact boundary geometry.
        a, b = station - step, min(station + step, self.length)
        derivative = (self.point(b, distance) - self.point(a, distance)) / (b - a)
        return point, derivative

    def physical_gradient(self, station, q):
        """Incoming world-XY height gradient, before inspection exaggeration.

        An ordinary Section keeps its profile constant along the path. A
        derived entrance whose width/crown changes along its course overrides
        this method with the gradient of that actual incoming height field.
        The profile value and derivative alone cannot describe such a field.
        """
        q = np.asarray(q)
        normal = self.frame(station)[2]
        velocity = self.jet(station, self.offset + .5 * self.width * q)[1]
        covector = right(velocity)
        denominator = np.sum(covector * (.5 * self.width * normal), axis=-1)
        if np.any(abs(denominator) <= 1e-10):
            raise ValueError('Incoming boundary profile has a singular tangent frame')
        return self.profile.derivative(q)[..., None] * covector / denominator[..., None]


@dataclass
class Corner:
    first: Approach
    second: Approach
    stations: tuple
    center: object
    radius: float
    start_angle: float
    sweep: float
    gap: float
    continuation: bool = False
    requested_radius: float | None = None
    radius_mode: str = 'STRICT'
    radius_attempts: tuple = ()


def _corner(first, second, radius, continuation_angle=0.):
    gap = math.atan2(float(cross(first.direction, second.direction)),
                     float(first.direction @ second.direction)) % math.tau
    if gap < math.radians(3.):
        raise ValueError('Two intersection approaches have indistinguishable outward directions')
    mismatch = abs(first.offset + second.offset + .5 * (second.width - first.width)) > 1e-6
    if abs(gap - math.pi) < 1e-6 or (mismatch and abs(gap - math.pi) <= continuation_angle):
        return Corner(first, second, (0., 0.), None, radius, 0., 0., gap, continuation=mismatch)
    side = 1. if gap < math.pi else -1.
    distances = (first.offset - first.width * .5 - side * radius,
                 second.offset + second.width * .5 + side * radius)
    p0 = first.point(0., distances[0])
    p1 = second.point(0., distances[1])
    matrix = np.column_stack((first.direction, -second.direction))
    if abs(np.linalg.det(matrix)) < 1e-8:
        raise ValueError('Corner offset curves have a singular initial contact')
    initial = np.linalg.solve(matrix, p1 - p0)
    lower = -4. * max(first.width + radius, second.width + radius)
    upper = np.array((first.length, second.length))
    current = np.clip(initial, lower, upper)
    solved = False
    for _ in range(32):
        jets = [a.jet(float(s), d) for a, s, d in zip((first, second), current, distances)]
        delta = jets[0][0] - jets[1][0]
        if np.linalg.norm(delta) < 1e-8:
            solved = True
            break
        jacobian = np.column_stack((jets[0][1], -jets[1][1]))
        if abs(np.linalg.det(jacobian)) < 1e-9:
            break
        step = np.linalg.solve(jacobian, delta)
        accepted = False
        for factor in (1., .5, .25, .125, .0625):
            candidate = np.clip(current - factor * step, lower, upper)
            error = first.point(candidate[0], distances[0]) - second.point(candidate[1], distances[1])
            if np.linalg.norm(error) < np.linalg.norm(delta):
                current = candidate
                accepted = True
                break
        if not accepted:
            break
    miter_frame = first.source.interpolation == 0 or second.source.interpolation == 0
    if not solved and not miter_frame:
        raise CornerContactError('Requested-radius curved corner has no converged contact in the supported branch')
    center = first.point(current[0], distances[0])
    contacts = (first.point(current[0], first.offset - first.width * .5),
                second.point(current[1], second.offset + second.width * .5))
    vectors = np.array(contacts) - center
    radius_error = np.max(abs(np.linalg.norm(vectors, axis=1) - radius))
    # Exact miter vectors are profile coordinates, not radius normals. Keep
    # existing qualified contacts unchanged; retry on the physical curbs when
    # a miter-scaled radius or a failed initial solve needs the true normals.
    miter_tangency_error = 0.
    if miter_frame:
        for approach, station, distance, vector in zip((first, second), current,
                (first.offset - first.width * .5, second.offset + second.width * .5), vectors):
            tangent = approach.jet(float(station), distance)[1]
            miter_tangency_error = max(miter_tangency_error, abs(float(vector @ tangent)) / np.linalg.norm(tangent))
    if not solved or radius_error > 1e-5 or miter_tangency_error > 1e-7:
        from .curb_contact import CircularCurbContact
        current, center, contacts = CircularCurbContact(first, second, radius, side).solve(current, lower, upper)
        vectors = np.array(contacts) - center
    if np.max(abs(np.linalg.norm(vectors, axis=1) - radius)) > 1e-5:
        raise ValueError('Corner contacts do not preserve the prescribed radius')
    start = math.atan2(vectors[0, 1], vectors[0, 0])
    sweep = math.atan2(float(cross(vectors[0], vectors[1])), float(vectors[0] @ vectors[1]))
    if sweep * side > 1e-7 or abs(sweep) > math.pi + 1e-7:
        raise CornerContactError('Corner contact selects the wrong circular sector')
    return Corner(first, second, tuple(map(float, current)), center, radius, start, sweep, gap)


class BoundaryCurve:
    def __init__(self, corner, stations, tolerance):
        self.corner = corner
        self.mouths = stations
        self.first = corner.first.point(stations[0], corner.first.offset - corner.first.width * .5)
        self.last = corner.second.point(stations[1], corner.second.offset + corner.second.width * .5)
        self.continuation = None
        if corner.center is None:
            if corner.continuation:
                first_velocity = corner.first.jet(stations[0], corner.first.offset - corner.first.width * .5)[1]
                last_velocity = corner.second.jet(stations[1], corner.second.offset + corner.second.width * .5)[1]
                self.continuation = WidthContinuation(self.first, self.last, -first_velocity, last_velocity, corner.radius)
                self.lengths = np.array((0., self.continuation.length))
                return
            self.lengths = np.array((0., float(np.linalg.norm(self.last - self.first))))
            return
        self.contacts = (corner.first.point(corner.stations[0], corner.first.offset - corner.first.width * .5),
                         corner.second.point(corner.stations[1], corner.second.offset + corner.second.width * .5))
        # These tables describe the actual offset curves along their finite
        # leads. Circular coordinates are evaluated analytically below.
        self.leads = []
        lengths = []
        for index, approach in enumerate((corner.first, corner.second)):
            a, b = ((stations[0], corner.stations[0]) if index == 0 else
                    (corner.stations[1], stations[1]))
            parameters = np.linspace(a, b, max(2, int(abs(b - a) / .2) + 2))
            distance = approach.offset + (-.5 if index == 0 else .5) * approach.width
            xy = approach.point(parameters, distance)
            cumulative = np.r_[0., np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
            self.leads.append((approach, parameters, cumulative, distance))
            lengths.append(cumulative[-1])
        self.lengths = np.r_[0., np.cumsum((lengths[0], abs(corner.sweep) * corner.radius, lengths[1]))]

    @property
    def length(self):
        return float(self.lengths[-1])

    def __call__(self, t):
        t = np.asarray(t)
        corner = self.corner
        if corner.center is None:
            if self.continuation is not None:
                return self.continuation(t)
            # For opposite approaches, preserve endpoint tangents even if the
            # widths differ; this is a smooth back curb, with zero profile rise.
            h = smooth(t)
            direction = -corner.first.direction
            delta = self.last - self.first
            longitudinal = np.dot(delta, direction) * direction
            lateral = delta - longitudinal
            return self.first + t[..., None] * longitudinal + h[..., None] * lateral
        station = t * self.length
        result = np.empty(t.shape + (2,))
        first = station <= self.lengths[1]
        last = station >= self.lengths[2]
        arc = ~(first | last)
        for mask, index, offset in ((first, 0, 0.), (last, 1, self.lengths[2])):
            if np.any(mask):
                approach, parameters, cumulative, distance = self.leads[index]
                if cumulative[-1] < 1e-10:
                    result[mask] = self.contacts[index]
                else:
                    s = np.interp(station[mask] - offset, cumulative, parameters)
                    result[mask] = approach.point(s, distance)
        if np.any(arc):
            q = (station[arc] - self.lengths[1]) / (self.lengths[2] - self.lengths[1])
            angle = corner.start_angle + q * corner.sweep
            result[arc] = corner.center + corner.radius * np.stack((np.cos(angle), np.sin(angle)), axis=-1)
        return result


@dataclass(frozen=True)
class Mouth:
    station: float
    transverse_intervals: int


@dataclass
class JunctionLayout:
    approaches: tuple
    corners: tuple
    boundaries: tuple
    outward_stations: tuple
    mouths: dict
    rows: int
    radial: int
    center: np.ndarray
    guide_weights: list
    mesh_strategy: str = 'bands'
    transition: object = None
    absorbed_transitions: tuple = ()
    absorbed_by: object = None


class IntersectionProvider(FeatureProvider):
    source_kinds = frozenset(('intersection', 'connection'))
    def tasks(self, network):
        from .ownership import PLAN, JunctionOwnership, footprint_key, surface_key
        from .transitions import required, product_key, output_key, TransitionProduct, neighboring_junctions
        topology = RouteTopology(network)
        attachments = topology.attachments(network)
        from .bundle_merged_entrances import MergedEntrancePreparation, preparation_key
        merged_candidates = MergedEntrancePreparation.candidates(network)
        merged_sites = {site: preparation_key(bundle) for bundle, crossing, runs in merged_candidates
                        for site in crossing.sites}
        if len(merged_sites) != sum(len(crossing.sites) for bundle, crossing, runs in merged_candidates):
            raise ValueError('Two Bundles claim the same merged physical entrance')
        yield from MergedEntrancePreparation.tasks(network, merged_candidates)
        for junction in (*network.intersections, *network.connections):
            transition_needed = required(network, junction)
            neighbors = neighboring_junctions(attachments, junction) if transition_needed else ()
            if len(junction.ends) > 2:
                for end in junction.ends:
                    candidate = topology.trace(end)
                    deps = tuple(sorted({*(Key(kind, visit.section_id) for visit in candidate.traversals
                                           for kind in ('path', 'profile')), *candidate.seams}))
                    signature = None
                    if candidate.terminal_site is not None and candidate.terminal_site not in deps:
                        geometry_deps = deps
                        deps = tuple(sorted((*deps, candidate.terminal_site)))
                        # The far incidence establishes an occupied endpoint;
                        # its radii do not change this route's path geometry.
                        signature = lambda context, keys=geometry_deps: [
                            (key.text(), context.digest(key)) for key in keys]
                    yield Task(route_key(end), deps, lambda context, source=candidate: Route.build(context, source),
                               signature=signature)
                deps = (junction.key, *(route_key(end) for end in junction.ends),
                        *(Key('profile', end.section_id) for end in junction.ends))
            else:
                deps = (junction.key, *(Key('path', end.section_id) for end in junction.ends),
                        *(Key('profile', end.section_id) for end in junction.ends),
                        *(footprint_key(j.key) for j in neighbors))
            if len(junction.ends) > 2:
                if junction.key in merged_sites:
                    deps = (*deps, merged_sites[junction.key])
                yield Task(footprint_key(junction.key), deps,
                    lambda context, source=junction, key=merged_sites.get(junction.key):
                        context[key].layout if key is not None and context[key] is not None
                        else self.footprint(context, source))
                yield Task(junction.layout_key, (footprint_key(junction.key), PLAN),
                    lambda context, source=junction: JunctionOwnership.layout(context, source),
                    signature=lambda context, source=junction: JunctionOwnership.layout_signature(context, source))
            else:
                signature = None
                if len(neighbors) > 1:
                    from .junction_transitions import JunctionTransitions
                    geometry_deps = deps
                    deps = (*deps, PLAN)
                    signature = lambda context, keys=geometry_deps, nearby=neighbors: JunctionTransitions.seam_signature(
                        context, keys, nearby)
                yield Task(junction.layout_key, deps,
                    lambda context, source=junction, nearby=neighbors: self.layout(context, source, transition_neighbors=nearby),
                    signature=signature)
            if transition_needed:
                from .bundle_boundaries import candidate_owner
                from .transitions import published_pavement, pavement_signature
                bundle = candidate_owner(network, junction)
                transition = product_key(junction.key)
                yield Task(transition, (junction.layout_key,),
                    lambda context, source=junction: TransitionProduct.build(context, context[source.layout_key].transition)
                        if context[source.layout_key].transition is not None else None)
                yield Task(output_key(junction.key), (transition,
                    *((Key('bundle_plan', bundle.id),) if bundle is not None else ())),
                    lambda context, owner=junction.key, candidate=bundle: published_pavement(context, owner, candidate),
                    signature=lambda context, owner=junction.key, candidate=bundle: pavement_signature(context,owner,candidate))
            elif junction.mesh_key is not None:
                deps = (junction.layout_key, *((surface_key(junction.key),) if len(junction.ends) > 2 else ()))
                yield Task(junction.mesh_key, deps,
                    lambda context, source=junction: JunctionOwnership.pavement(context, source))
        yield from JunctionOwnership.tasks(network)

    @staticmethod
    def footprint(context, junction, transition_neighbors=()):
        """Solve only prescribed curbs and source trims, before choosing a mesh.

        Compound planning needs these physical extents even when the local
        corner-band guide family cannot parameterize the resulting surface.
        """
        approaches = [Approach(context[route_key(end) if len(junction.ends) > 2 else Key('path', end.section_id)], end.end,
                              context[Key('profile', end.section_id)]) for end in junction.ends]
        from .junction_transitions import JunctionTransitions
        if len(junction.ends) > 2:
            return JunctionTransitions.footprint(context, junction, approaches)
        absorbed = JunctionTransitions.absorbed(context, junction, approaches, transition_neighbors)
        if absorbed is not None:
            return absorbed
        return IntersectionProvider._footprint(context, junction, approaches, transition_neighbors)

    @staticmethod
    def _footprint(context, junction, approaches, transition_neighbors=()):
        if max(np.linalg.norm(a.origin - approaches[0].origin) for a in approaches) > 1e-5:
            raise ValueError(f'Intersection {junction.id} endpoint copies disagree; move connected copies together')
        return IntersectionProvider.curved_boundary(context, junction, approaches, transition_neighbors)

    @staticmethod
    def curved_boundary(context, junction, approaches, transition_neighbors=()):
        """Physical curb contacts and complete profiles, independent of origin.

        Local intersections validate their common endpoint in _footprint.
        A compound may use distinct original endpoint origins after its source
        incidence and consumed-span ownership have been established.
        """
        approaches.sort(key=lambda a: math.atan2(a.direction[1], a.direction[0]) % math.tau)
        if len(approaches) == 2:
            first, second = approaches
            from .transitions import mixed_bend_required
            if (abs(first.width - second.width) > 1e-6 or abs(first.offset + second.offset) > 1e-6
                    or abs(first.profile.rise - second.profile.rise) > 1e-7
                    or abs(first.profile.datum - second.profile.datum) > 1e-8
                    or first.profile.half + second.profile.half != 0
                    or mixed_bend_required(first.source,first.end,second.source,second.end)):
                from .transitions import SectionTransition, TransitionClearance
                clearance = TransitionClearance.build(context, approaches, transition_neighbors)
                transition = SectionTransition.solve(context, junction, approaches, clearance)
                radial = half_intervals(context.network.settings)
                mouths = {(a.source.id, a.end): Mouth(d if a.end == 0 else a.length - d, 2 * radial)
                          for a, d in zip(approaches, transition.distances)}
                return JunctionLayout(tuple(approaches), (), (), transition.distances, mouths,
                    transition.rows, radial, first.origin.copy(), [0., 0.], 'transition', transition)
            if np.linalg.norm(first.direction + second.direction) > 1e-5:
                raise ValueError('Exact bent Section seams need an explicit common miter')
            radial = half_intervals(context.network.settings)
            mouths = {(a.source.id, a.end): Mouth(0. if a.end == 0 else a.length, 2 * radial) for a in approaches}
            return JunctionLayout(tuple(approaches), (), (), (0., 0.), mouths, 0, radial,
                                  first.origin.copy(), [0., 0.])
        corners = []
        def entrance_fits(corner):
            return all(IntersectionProvider.mouth_fits(a, IntersectionProvider.mouth_station(
                a, max(0., s) + min(.25, a.width * .05), context.network.settings))
                for a, s in zip((corner.first, corner.second), corner.stations))
        for i, a in enumerate(approaches):
            b = approaches[(i + 1) % len(approaches)]
            radius, mode = junction.corner_radius(End(a.source.id, a.end))
            try:
                corners.append(CornerRadiusPolicy.solve(a, b, radius, mode,
                    math.radians(context.network.settings.continuation_angle_degrees),
                    accept=entrance_fits if mode == 'STYLE' else None))
            except ValueError as exc:
                message = f'Corner between Section {a.source.id} end {a.end} and '
                message += (
                    f'Section {b.source.id} end {b.end}, radius {radius:g} m, '
                    f'available approaches {a.length:.3f}/{b.length:.3f} m: {exc}')
                if isinstance(exc, CornerContactError):
                    raise CornerConstructionError(message, junction, (a, b)) from exc
                raise ValueError(message) from exc
        stations = [max(0., corners[i].stations[0], corners[i - 1].stations[1])
                    for i in range(len(approaches))]
        # A small lead permits independent side contacts and keeps complete
        # incoming profiles at the boundary of the generated intersection.
        settings = context.network.settings
        stations = [IntersectionProvider.mouth_station(a, s + min(.25, a.width * .05), settings)
                    for s, a in zip(stations, approaches)]
        if any(not IntersectionProvider.mouth_fits(a, s)
               for s, a in zip(stations, approaches)):
            raise ApproachExtentError(junction, [(a, s) for s, a in zip(stations, approaches)
                                                if s >= a.length - 1e-6])
        boundaries = [BoundaryCurve(corner, (stations[i], stations[(i + 1) % len(stations)]),
                                    context.network.settings.curve_tolerance)
                      for i, corner in enumerate(corners)]
        settings = context.network.settings
        radial = half_intervals(settings)
        rows = 2 * radial
        for boundary in boundaries:
            corner = boundary.corner
            chord = settings.spacing
            if corner.center is not None:
                chord = min(chord, 2. * math.sqrt(max(1e-12,
                    2. * corner.radius * settings.curve_tolerance - settings.curve_tolerance ** 2)))
            rows = max(rows, int(math.ceil(boundary.length / chord)))
            if boundary.continuation is not None:
                rows = max(rows, boundary.continuation.intervals(settings.spacing, settings.curve_tolerance))
        rows = 4 * int(math.ceil(rows / 4.))
        if any(a.profile.intervals(-1, radial) == 0 and
               approaches[(i + 1) % len(approaches)].profile.intervals(1, radial) == 0
               for i, a in enumerate(approaches)):
            rows = 8 * int(math.ceil(rows / 8.))
        mouths = {}
        for approach, distance in zip(approaches, stations):
            for identity, station in approach.path.trims(distance).items():
                mouths[identity] = Mouth(float(station), 2 * radial)
        layout = JunctionLayout(tuple(approaches), tuple(corners), tuple(boundaries), tuple(stations),
                                mouths, rows, radial, np.zeros(2), [0.] * len(approaches))
        return layout

    @staticmethod
    def mouth_fits(approach, station):
        return (station < approach.length - 1e-6 or
                (getattr(approach.path, 'terminal_free', False) and abs(station - approach.length) <= 1e-9))

    @staticmethod
    def mouth_station(approach, station, settings):
        """Absorb a sub-cell source remnant while retaining its full profile.

        Circle contacts stay fixed; only the curb lead extends to the
        next compatible Section seam or genuine free end. An incompatible
        continuation or another physical junction may not be consumed. The
        bound is the length below which the normal two-row road grid would
        exceed edge ratio eight.
        """
        if not isinstance(approach.path, Route):
            return station
        station = max(station, getattr(approach.path, 'minimum_mouth', 0.))
        boundaries = approach.path.stations[1:]
        if not approach.path.terminal_free:
            boundaries = boundaries[:-1]
        at = int(np.searchsorted(boundaries, station))
        if at == len(boundaries):
            return station
        target = float(boundaries[at])
        gap = target - station
        limit = min(.5 * settings.spacing, approach.width / (8. * half_intervals(settings)))
        if 0. <= gap <= limit:
            return target
        return station

    @staticmethod
    def layout(context, junction, footprint=None, transition_neighbors=()):
        layout = footprint if footprint is not None else IntersectionProvider.footprint(context, junction, transition_neighbors)
        # Mesh guide selection must not mutate the cached physical footprint,
        # which also belongs to compound planning and other feature providers.
        layout = replace(layout, center=layout.center.copy(), guide_weights=list(layout.guide_weights))
        if not layout.corners:
            return layout
        approaches, boundaries, radial = layout.approaches, layout.boundaries, layout.radial
        for side, boundary in enumerate(boundaries):
            pattern = IntersectionProvider.pattern(context, layout, side)
            t, v = pattern.uv.T
            faces = pattern.quads
            collapsed = (approaches[side].profile.intervals(-1, radial) == 0 or
                         approaches[(side + 1) % len(approaches)].profile.intervals(1, radial) == 0)
            # Chord-first selection can leave a very thin radial band at an
            # oblique corner: positive quads alone do not prevent steep heights.
            # Prefer a broad curved guide, retaining bounded fallback choices.
            # A collapsed cap needs the fully tangent variant at its tip.
            for weight in ((1., .75, .5, .25, 0.) if collapsed else (.75, .5, .25, 0.)):
                layout.guide_weights[side] = weight
                xy, z = IntersectionProvider.band_surface(layout, side, t, v)
                points = xy[faces]
                edges = np.roll(points, -1, axis=1) - points
                ceiling = max(a.profile.rise for a in approaches)
                if (np.min(cross(edges, np.roll(edges, -1, axis=1))) > 1e-10
                        and np.min(z) >= -1e-8 and np.max(z) <= ceiling + 1e-8):
                    break
            else:
                raise ValueError('Corner band needs a different compatible guide pattern')
        layout.center = np.mean([IntersectionProvider.guide(layout, i, np.asarray(.5))
                                 for i in range(len(approaches))], axis=0)
        return layout

    @staticmethod
    def guide_controls(layout, side):
        first, second = layout.approaches[side], layout.approaches[(side + 1) % len(layout.approaches)]
        a = first.point(layout.outward_stations[side], first.offset + first.profile.peak * .5 * first.width)
        b = second.point(layout.outward_stations[(side + 1) % len(layout.approaches)],
                         second.offset + second.profile.peak * .5 * second.width)
        h0, h1 = .7 * layout.outward_stations[side], .7 * layout.outward_stations[(side + 1) % len(layout.approaches)]
        collapse0, collapse1 = first.profile.peak == -1., second.profile.peak == 1.
        # A cap's two bounding curves must approach its tip at comparable
        # speeds. Otherwise a valid continuum map can still fold a finite quad
        # across its very thin last wedge. Keep the noncollapsed handle, which
        # places the guide inside the physical circular curb.
        if collapse0 and not collapse1:
            h0 = min(h0, layout.boundaries[side].length / 3.)
        if collapse1 and not collapse0:
            h1 = min(h1, layout.boundaries[side].length / 3.)
        p = a - h0 * first.frame(layout.outward_stations[side])[1]
        q = b - h1 * second.frame(layout.outward_stations[(side + 1) % len(layout.approaches)])[1]
        return a, p, q, b

    @staticmethod
    def guide(layout, side, t):
        a, p, q, b = IntersectionProvider.guide_controls(layout, side)
        t = np.asarray(t)[..., None]
        # The canonical pattern's chord-guide family has a nonsingular core.
        # A pair of longitudinally tangent crown guides creates a zero-angle
        # cusp, which can invert Coons cells after even a small source edit.
        # Guides are internal coordinates, not the prescribed physical curbs.
        chord = (1. - t) * a + t * b
        curved = (1. - t) ** 3 * a + 3. * (1. - t) ** 2 * t * p + 3. * (1. - t) * t ** 2 * q + t ** 3 * b
        weight = layout.guide_weights[side]
        return (1. - weight) * chord + weight * curved

    @staticmethod
    def height(layout, side, t):
        a = layout.approaches[side].profile.rise
        b = layout.approaches[(side + 1) % len(layout.approaches)].profile.rise
        return a + (b - a) * smooth(np.asarray(t))

    @staticmethod
    def pattern(context, layout, side):
        first = layout.approaches[side].profile.intervals(-1, layout.radial)
        last = layout.approaches[(side + 1) % len(layout.approaches)].profile.intervals(1, layout.radial)
        return context.patterns.band(layout.rows, first, last, middle=layout.radial)

    @classmethod
    def band_surface(cls, layout, side, t, v):
        """Return physical XY and profile rise, including source tangent planes."""
        a, b = layout.approaches[side], layout.approaches[(side + 1) % len(layout.approaches)]
        outer, guide = layout.boundaries[side](t), cls.guide(layout, side, t)
        xy = (1. - v[:, None]) * outer + v[:, None] * guide
        qa, qb = a.profile.coordinate(-1, v), b.profile.coordinate(1, v)
        mix = smooth(t)
        z = (1. - mix) * a.profile.value(qa) + mix * b.profile.value(qb)
        chord = cls.guide(layout, side, np.asarray(1.)) - cls.guide(layout, side, np.asarray(0.))
        weight = layout.guide_weights[side]
        anchor_a, handle_a, handle_b, anchor_b = cls.guide_controls(layout, side)
        derivative_a = (1. - weight) * chord + weight * 3. * (handle_a - anchor_a)
        derivative_b = (1. - weight) * chord + weight * 3. * (anchor_b - handle_b)
        boundary = layout.boundaries[side]
        epsilon = 1e-6
        outer_a = (boundary(np.asarray(epsilon)) - boundary(np.asarray(0.))) / epsilon
        outer_b = (boundary(np.asarray(1.)) - boundary(np.asarray(1. - epsilon))) / epsilon
        slopes = []
        for approach, q, parameter, d_outer, d_guide in (
            (a, qa, layout.outward_stations[side], outer_a, derivative_a),
            (b, qb, layout.outward_stations[(side + 1) % len(layout.approaches)], outer_b, derivative_b)):
            normal = approach.frame(parameter)[2]
            velocity = approach.jet(parameter, approach.offset + q * .5 * approach.width)[1]
            covector = right(velocity)
            denominator = np.sum(covector * (.5 * approach.width * normal), axis=-1)
            if np.min(abs(denominator)) < 1e-9:
                raise ValueError('Incoming profile frame is singular')
            gradient = approach.profile.derivative(q)[:, None] * covector / denominator[:, None]
            delta = (1. - v[:, None]) * d_outer + v[:, None] * d_guide
            slopes.append(np.sum(gradient * delta, axis=-1))
        z += t * (1. - t) ** 2 * slopes[0] + t * t * (t - 1.) * slopes[1]
        return xy, z

    @classmethod
    def build(cls, context, junction):
        return cls._build(context, junction, planar=False)

    @classmethod
    def planar_chart(cls, context, junction):
        """Preparation geometry only; its Z=0 mesh requires a height solver."""
        return cls._build(context, junction, planar=True)

    @classmethod
    def _build(cls, context, junction, *, planar):
        layout = context[junction.layout_key]
        if not layout.corners:
            # The two paths share a complete compatible profile. Its seam has
            # no area and therefore owns no duplicate pavement object.
            return None
        settings = context.network.settings
        n, m, h = layout.rows, layout.radial, layout.rows // 2
        sides = len(layout.approaches)
        datum = layout.approaches[0].profile.datum
        if any(abs(a.profile.datum - datum) > 1e-7 for a in layout.approaches):
            raise ValueError('Incoming minimum heights disagree')
        ceiling = max(a.profile.rise for a in layout.approaches)
        center_height = float(np.mean([a.profile.rise for a in layout.approaches]))
        assembly = Assembler()

        def key(owner, row):
            return ('guide', (owner + 1) % sides, 0) if row == n else ('guide', owner, row)

        def xyz(xy, z):
            if planar:
                return np.column_stack((xy, np.zeros(len(xy))))
            return np.column_stack((xy, datum + settings.exaggeration * z))

        for side in range(sides):
            previous = (side - 1) % sides
            pattern = cls.pattern(context, layout, side)
            t, v = pattern.uv.T
            if planar:
                xy = ((1. - v[:, None]) * layout.boundaries[side](t)
                      + v[:, None] * cls.guide(layout, side, t))
                z = np.zeros(len(xy))
            else:
                xy, z = cls.band_surface(layout, side, t, v)
                if np.min(z) < -1e-8 or np.max(z) > ceiling + 1e-8:
                    raise ValueError('Corner profile reconciliation exceeds its elevation envelope')
                z = np.clip(z, 0., ceiling)  # Remove only already-checked roundoff.
            first_count = layout.approaches[side].profile.intervals(-1, m)
            last_count = layout.approaches[(side + 1) % sides].profile.intervals(1, m)
            keys = [key(side, int(round(ti * n))) if vi > 1. - 1e-10 or
                    (ti < 1e-10 and first_count == 0) or (ti > 1. - 1e-10 and last_count == 0)
                    else ('band', side, i) for i, (ti, vi) in enumerate(pattern.uv)]
            assembly.add(xyz(xy, z), pattern.quads, keys)

            pattern = context.patterns.grid(h, h)
            u, v = pattern.uv.T
            a = cls.guide(layout, side, np.asarray(0.))
            b = cls.guide(layout, side, np.asarray(.5))
            d = cls.guide(layout, previous, np.asarray(.5))
            c = layout.center
            bottom = cls.guide(layout, side, u * .5)
            left = cls.guide(layout, previous, 1. - v * .5)
            top = (1. - u[:, None]) * d + u[:, None] * c
            right_side = (1. - v[:, None]) * b + v[:, None] * c
            bilinear = ((1. - u) * (1. - v))[:, None] * a + (u * (1. - v))[:, None] * b
            bilinear += ((1. - u) * v)[:, None] * d + (u * v)[:, None] * c
            xy = ((1. - v[:, None]) * bottom + v[:, None] * top
                  + (1. - u[:, None]) * left + u[:, None] * right_side - bilinear)
            # Smooth transfinite extension of guide heights. The explicit
            # envelope retains unequal entrance crowns and prevents overshoot.
            su, sv = smooth(u), smooth(v)
            ha = cls.height(layout, side, 0.)
            hb = cls.height(layout, side, .5)
            hd = cls.height(layout, previous, .5)
            hbnd = cls.height(layout, side, u * .5)
            hl = cls.height(layout, previous, 1. - v * .5)
            ht = (1. - su) * hd + su * center_height
            hr = (1. - sv) * hb + sv * center_height
            base = (1. - su) * (1. - sv) * ha + su * (1. - sv) * hb
            base += (1. - su) * sv * hd + su * sv * center_height
            z = np.clip((1. - sv) * hbnd + sv * ht + (1. - su) * hl + su * hr - base, 0., ceiling)
            keys = []
            for i in range(h + 1):
                for j in range(h + 1):
                    identity = ('core', side, i, j)
                    if i == h:
                        identity = ('spoke', side, j)
                    if j == h:
                        identity = ('spoke', previous, i)
                    if i == h and j == h:
                        identity = ('center',)
                    if i == 0:
                        identity = key(previous, n - j)
                    if j == 0:
                        identity = key(side, i)
                    keys.append(identity)
            assembly.add(xyz(xy, z), pattern.quads, keys)
        mesh = assembly.mesh()
        quality = mesh.validate()
        if planar:
            return Part(junction.key, mesh, dict(quality=quality, rows=n, radial=m,
                planar_only=True, construction='planar corner bands and direct cyclic Coons blocks'))
        if quality['maximum_edge_ratio'] > 8.:
            raise ValueError(f'Intersection {junction.id} quad edge ratio exceeds 8: '
                             f'{quality["maximum_edge_ratio"]:.3f}')
        return Part(junction.key, mesh, dict(quality=quality, rows=n, radial=m,
                    fitted_corner_radii=fitted_radius_report(junction.key, layout),
                    crown_ceiling=ceiling, guide_weights=list(layout.guide_weights),
                    approach_sections=[[piece.section.id for piece in a.path.pieces] for a in layout.approaches],
                    physical_slopes=mesh.slope_statistics(settings.exaggeration),
                    construction='corner bands and direct cyclic Coons blocks'))
