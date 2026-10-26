"""Early physical entrances for overlapping, curved Bundle members.

Preparation consumes original paths, profiles and incidence, never completed
surfaces. The two logical junctions and their consumed Street stay original.
"""
from dataclasses import dataclass
import math
import numpy as np
from .model import End, Key
from .updates import Task
from .bundle_axis import AxialChartUnavailable
from .bundle_angular_curves import AngularCrown
from .bundle_angular_internal import shared_frame
from .bundle_regions import interior_crossings
from .bundle_visible_runs import VisibleBundleRuns
from .corner_radii import CornerRadiusPolicy
from .intersections import Approach, BoundaryCurve, JunctionLayout, Mouth, half_intervals
from .paths import unit
from .routes import Piece, route_key
from .source_mouths import SourceMouth


def preparation_key(bundle):
    return Key('merged_entrances', bundle.id)


@dataclass(frozen=True)
class MergedCrossSection:
    rise: float
    datum: float
    peak: float
    half: int = 0

    def value(self, q):
        q = np.asarray(q)
        width = np.where(q <= self.peak, 1.+self.peak, 1.-self.peak)
        return self.rise*(1.-((q-self.peak)/width)**2)

    def derivative(self, q):
        q = np.asarray(q)
        width = np.where(q <= self.peak, 1.+self.peak, 1.-self.peak)
        return -2.*self.rise*(q-self.peak)/width**2

    def coordinate(self, side, v):
        return side+(self.peak-side)*np.asarray(v)

    def intervals(self, side, base):
        return max(2, 2*int(round(abs(self.peak-side)*base/2.)))


@dataclass(frozen=True)
class MemberInventory:
    """Metadata inventory only; deliberately exposes no Route.trims API."""
    pieces: tuple


class PairedAngularEntrance:
    maximum_samples = 20000

    def __init__(self, context, run, end, frame):
        self.frame_chart, self.end_index = frame, end
        self.original_ends = tuple(member.ends[end] for member in run.members)
        self.visits = tuple(member.visits[0] for member in run.members)
        self.paths = tuple(context[Key('path', visit.section_id)] for visit in self.visits)
        self.path = MemberInventory(tuple(Piece(path, True) for path in self.paths))
        curves, profiles, origins = [], [], []
        for role, (visit, path) in enumerate(zip(self.visits, self.paths)):
            profile = context[Key('profile', visit.section_id)].oriented(int(not visit.forward))
            if profile.half != (1 if role == 0 else -1):
                raise ValueError('Merged entrance lost its original half-profile orientation')
            profiles.append(profile)
            frame.certify(path, forward=visit.forward)
            samples = path.samples(0., path.length, context.network.settings.spacing,
                                   tolerance=.1*context.network.settings.curve_tolerance)
            if len(samples) > self.maximum_samples:
                raise ValueError('Merged entrance exceeds its source sampling budget')
            for side in (-1., 1.):
                q = side*(1 if visit.forward else -1)
                points = path.point(samples, path.section.offset+q*.5*path.section.width)
                curves.append(points if visit.forward else points[::-1])
            origins.append(path.frame(0. if self.original_ends[role].end == 0 else path.length)[0])
        if abs(profiles[0].rise-profiles[1].rise)>1e-9 or abs(profiles[0].datum-profiles[1].datum)>1e-9:
            raise ValueError('Merged Bundle entrance requires its original shared crown')
        self.field = AngularCrown(frame, curves, datum=profiles[0].datum, rise=profiles[0].rise)
        self.sign = 1. if end == 0 else -1.
        self.origin_station = float(np.mean(frame.local(origins)[:, 0]))
        self.length = (min(ray.stations[-1] for ray in self.field.curves)-self.origin_station
            if self.sign > 0. else self.origin_station-max(ray.stations[0] for ray in self.field.curves))
        high, low = self.bounds(0.)
        self.width, self.offset = float(high-low), 0.
        self.source, self.end = self.paths[0].section, self.original_ends[0].end
        self.origin = self._positive_point(np.asarray(0.), np.asarray(0.))
        self.direction = unit(self._positive_point(np.asarray(1e-4), np.asarray(0.))-self.origin)
        self.profile = self.cross_section(0.)

    def chart_station(self, station):
        return self.origin_station+self.sign*np.asarray(station)

    def bounds(self, station):
        ao, ai, bi, bo = self.field.edges(self.chart_station(station))
        return np.maximum(ao, bi), np.minimum(ai, bo)

    def gap(self, station):
        ao, ai, bi, bo = self.field.edges(self.chart_station(station))
        return np.maximum(ai, bo)-np.minimum(ao, bi)

    def _positive_point(self, station, distance):
        high, low = self.bounds(station)
        y = .5*(high+low)-self.sign*np.asarray(distance)*(high-low)/self.width
        return self.frame_chart.world(np.stack(np.broadcast_arrays(self.chart_station(station), y), axis=-1))

    def point(self, station, distance):
        station = np.asarray(station)
        return (self._positive_point(np.clip(station, 0., self.length), np.asarray(distance))
                +np.minimum(station, 0.)[..., None]*self.direction)

    def jet(self, station, distance):
        a, b = np.asarray(station)-1e-4, np.minimum(np.asarray(station)+1e-4, self.length)
        return self.point(station, distance), (self.point(b, distance)-self.point(a, distance))/(b-a)[..., None]

    def frame(self, station):
        point, tangent = self.jet(station, 0.)
        normal = (self.point(station, .5*self.width)-self.point(station, -.5*self.width))/self.width
        return point, unit(tangent), normal

    def physical_gradient(self, station, q):
        return self.field.gradient(self.point(station, .5*self.width*np.asarray(q)))

    def cross_section(self, station):
        if self.gap(station) >= -1e-7:
            raise AxialChartUnavailable('Merged entrance requires actual pavement overlap at its full mouth')
        first, last = self.field.peaks(self.chart_station(station))
        if abs(first-last)>1e-9:
            raise AxialChartUnavailable('A separator-bearing entrance needs its split profile policy')
        high, low = self.bounds(station)
        peak = -self.sign*(2.*first-high-low)/(high-low)
        if not -1.+1e-8 < peak < 1.-1e-8:
            raise ValueError('Merged entrance requires a positive width on both crown halves')
        return MergedCrossSection(self.field.rise, self.field.datum, float(peak))

    def corner_end(self, side):
        role = 0 if (side < 0.) == (self.sign > 0.) else 1
        return self.original_ends[role]

    def source_mouths(self, station):
        """Intersect each certified original centreline with the radial cut."""
        target = float(self.chart_station(station))
        peak_u = .5*(self.cross_section(station).peak+1.)
        aliases = []
        for role, (path, visit, endpoint) in enumerate(zip(self.paths, self.visits, self.original_ends)):
            low, high = 0., path.length
            first, last = self.frame_chart.local(path.frame(np.array((low, high)))[0])[:, 0]
            if not min(first, last)-1e-8 <= target <= max(first, last)+1e-8:
                raise ValueError('Merged radial cut lies outside an original member centreline')
            for _ in range(56):
                middle = .5*(low+high)
                value = float(self.frame_chart.local(path.frame(middle)[0])[0])
                if (value < target) == (last > first):
                    low = middle
                else:
                    high = middle
            at = .5*(low+high)
            if abs(float(self.frame_chart.local(path.frame(at)[0])[0])-target)>1e-7:
                raise ValueError('Merged source-mouth station did not converge')
            outer_u = float(role if self.sign > 0. else 1-role)
            factor = (1 if endpoint.end == 0 else -1)*(1 if role == 0 else -1)*(1 if visit.forward else -1)
            aliases.append(SourceMouth(endpoint, at, outer_u, peak_u, .5*(1-factor), .5*(1+factor)))
        return tuple(aliases)


@dataclass(frozen=True)
class MergedEntranceSpecification:
    bundle_id: int
    crossing: object
    frame: object
    layout: object
    group: object
    centerline_proof_checks: int

    @classmethod
    def build(cls, context, bundle, crossing, runs):
        try:
            frame = shared_frame(context, bundle)
            approaches = [PairedAngularEntrance(context, run, end, frame)
                for run in runs for end, token in enumerate(run.boundaries)
                if token == ('crossing', crossing.street_id)]
        except AxialChartUnavailable:
            return None  # Retain the established separated/internal provider.
        from .bundles import certify_disjoint
        proof_checks = certify_disjoint(tuple(context[Key('path', visit.section_id)] for visit in bundle.members[0].path),
                                        tuple(context[Key('path', visit.section_id)] for visit in bundle.members[1].path))
        sites = {site.key: site for site in context.network.intersections}
        exterior = {end for key in crossing.sites for end in sites[key].ends
                    if end.section_id not in crossing.sections}
        paired = {end for approach in approaches for end in approach.original_ends}
        for endpoint in sorted(exterior-paired, key=lambda end: (end.section_id, end.end)):
            approaches.append(Approach(context[route_key(endpoint)], endpoint.end,
                                       context[Key('profile', endpoint.section_id)]))
        approaches.sort(key=lambda approach: math.atan2(approach.direction[1], approach.direction[0]) % math.tau)
        by_end = {end: site for site in context.network.intersections for end in site.ends}
        corners = []
        for index, first in enumerate(approaches):
            second = approaches[(index+1)%len(approaches)]
            endpoint = (first.corner_end(-1.) if isinstance(first, PairedAngularEntrance)
                        else End(first.source.id, first.end))
            corners.append(CornerRadiusPolicy.solve(first, second, *by_end[endpoint].corner_radius(endpoint)))
        stations = tuple(max(0., corner.stations[0], corners[index-1].stations[1])+.25
                         for index, corner in enumerate(corners))
        mouths = {}
        radial = half_intervals(context.network.settings)
        from .source_mouths import source_mouths
        for approach, station in zip(approaches, stations):
            if station >= approach.length-1e-6:
                raise ValueError('Merged exterior consumes an entire original approach')
            if isinstance(approach, PairedAngularEntrance):
                approach.profile = approach.cross_section(station)
            for alias in source_mouths(approach, station):
                mouths[alias.endpoint.section_id, alias.endpoint.end] = Mouth(alias.station, 2*radial)
            if hasattr(approach.path, 'trims'):
                for identity, at in approach.path.trims(station).items():
                    mouths[identity] = Mouth(float(at), 2*radial)
        for site in crossing.sites:
            for endpoint in sites[site].ends:
                if endpoint.section_id in crossing.sections:
                    length = context[Key('path', endpoint.section_id)].length
                    mouths[endpoint.section_id, endpoint.end] = Mouth(length if endpoint.end == 0 else 0., 2*radial)
        boundaries = tuple(BoundaryCurve(corner, (stations[index], stations[(index+1)%len(stations)]),
            context.network.settings.curve_tolerance) for index, corner in enumerate(corners))
        layout = JunctionLayout(tuple(approaches), tuple(corners), boundaries, stations, mouths,
            4*int(math.ceil(max(2*radial, max(b.length/context.network.settings.spacing for b in boundaries))/4.)),
            radial, np.zeros(2), [0.]*len(approaches), mesh_strategy='boundary_chart')
        from .bundle_merged_rim import merged_footprint
        from .compounds import CompoundGroup, SharedClaim
        from .regions import RegionUnion
        owner = min(crossing.sites)
        sampled = merged_footprint(owner, layout, context.network.settings)
        loops = RegionUnion().build((sampled.region,))
        if len(loops) != 1 or loops[0].hole:
            raise ValueError('Merged entrance requires one simple actual exterior')
        from .compound_ports import CompoundPortPlanner
        CompoundPortPlanner(context, tuple(sites.values()), {}).check_support(crossing.sections, loops[0])
        claims = tuple(SharedClaim(sid, by_end[End(sid, 0)].key, by_end[End(sid, 1)].key,
                                  context[Key('path', sid)].length, 0.) for sid in crossing.sections)
        group = CompoundGroup(tuple(sorted(crossing.sites)), claims, (sampled,), loops,
                              physical_layout=sampled.layout)
        return cls(bundle.id, crossing, frame, layout, group, proof_checks)


class MergedEntrancePreparation:
    @staticmethod
    def candidates(network):
        result = []
        for bundle in network.bundles:
            if len(bundle.members) != 2 or not bundle.interior_section_ids or not all(m.complete for m in bundle.members):
                continue
            crossings = interior_crossings(network, bundle)
            if len(crossings) != 1 or len(crossings[0].path) != 1:
                continue
            try:
                runs = VisibleBundleRuns.build(network, bundle, crossings)
            except AxialChartUnavailable:
                continue
            if len(runs) == 2 and all(len(member.visits) == 1 for run in runs for member in run.members):
                result.append((bundle, crossings[0], runs))
        return tuple(result)

    @classmethod
    def tasks(cls, network, candidates):
        sites = {site.key: site for site in network.intersections}
        for bundle, crossing, runs in candidates:
            ends = {end for key in crossing.sites for end in sites[key].ends}
            ids = {*bundle.sections, *crossing.sections}
            deps = tuple(sorted({bundle.key, *crossing.sites,
                *(Key(kind, sid) for sid in ids for kind in ('path', 'profile')),
                *(route_key(end) for end in ends if end.section_id not in bundle.sections)}))
            yield Task(preparation_key(bundle), deps,
                lambda context, b=bundle, c=crossing, r=runs: MergedEntranceSpecification.build(context, b, c, r))
