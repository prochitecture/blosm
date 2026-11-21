"""Opt-in preparation of a junction on one member of overlapping pavement.

There is deliberately no production candidate/discovery hook. A caller must
supply the actual original site, both adjacent near-member visits and the
explicit opposite Bundle visit. Physical cuts restrict evaluation intervals;
they create no source Section, Street, persistent node or refitted Path.
"""
from dataclasses import dataclass, replace
from types import SimpleNamespace
import math
import numpy as np

from .bundle_axis import AxialChartUnavailable
from .bundle_merged_entrances import MemberInventory
from .finite_terminal_entrances import FinitePairedEntrance, LinearFrame, LinearCrown
from .model import End, Key
from .paths import unit
from .routes import Piece


@dataclass(frozen=True)
class SourceEvaluationInterval:
    section_id: int
    first: float
    last: float
    forward: bool

    def __post_init__(self):
        if not np.isfinite((self.first, self.last)).all() or not 0. <= self.first < self.last:
            raise ValueError('A source evaluation interval needs positive original metric support')


@dataclass(frozen=True)
class DerivedSourcePortKey:
    owner: Key
    section_id: int
    end: int

    def __post_init__(self):
        if self.end not in (0, 1):
            raise ValueError('Derived physical source orientation must be zero or one')


def _station_at(path, frame, target, interval):
    low, high = interval
    first, last = frame.local(path.frame(np.array((low, high)))[0])[:, 0]
    if not min(first, last)-1e-8 <= target <= max(first, last)+1e-8:
        raise AxialChartUnavailable('A paired interior fibre leaves its original source support')
    for _ in range(56):
        middle = .5*(low+high)
        value = float(frame.local(path.frame(middle)[0])[0])
        if (value < target) == (last > first):
            low = middle
        else:
            high = middle
    station = .5*(low+high)
    if abs(float(frame.local(path.frame(station)[0])[0])-target) > 1e-8:
        raise AxialChartUnavailable('An interior source fibre station did not converge')
    return float(station)


class InteriorPairedEntrance(FinitePairedEntrance):
    """Bounded common fibre with one real incidence and one interior source cut."""
    def __init__(self, context, near_visit, opposite_visit, near_end, near_role, *,
                 source_tolerance=.00025, source_parameters=None, protected_source_parameters=None):
        from .bundles import CenterlineRange
        from .bundle_boundaries import _certify_increasing
        if near_role not in (0, 1) or near_end not in (0, 1):
            raise ValueError('An interior paired entrance needs an original near-side orientation')
        if not math.isfinite(source_tolerance) or not 0. < source_tolerance <= .00025:
            raise ValueError('Interior source chord tolerance must be at most a quarter millimetre')
        visits = [None, None]
        visits[near_role], visits[1-near_role] = near_visit, opposite_visit
        self.visits, self.near_role = tuple(visits), near_role
        self.protected_source_parameters = ({} if protected_source_parameters is None else
            {sid: tuple(sorted(set(map(float, values)))) for sid, values in protected_source_parameters.items()})
        self.paths = tuple(context[Key('path', v.section_id)] for v in visits)
        near = self.paths[near_role]
        self.end_index = int(near_end == int(near_visit.forward))
        self.sign = 1. if self.end_index == 0 else -1.
        near_station = 0. if near_end == 0 else near.length
        point, tangent, _ = near.frame(near_station)
        axis = tangent*(1. if near_visit.forward else -1.)
        opposite = self.paths[1-near_role]
        # The opposite visit was selected from explicit Bundle membership.
        # Projection only locates a finite evaluation range on that Path.
        provisional = LinearFrame(point, axis)
        anchor = _station_at(opposite, provisional, 0., (0., opposite.length))
        other_point, other_tangent, _ = opposite.frame(anchor)
        other_axis = other_tangent*(1. if opposite_visit.forward else -1.)
        if float(unit(axis)@unit(other_axis)) < math.cos(math.radians(30.)):
            raise AxialChartUnavailable('Interior paired source directions differ by more than thirty degrees')
        self.frame_chart = LinearFrame(.5*(point+other_point), axis+other_axis)
        origins = [None, None]
        origins[near_role], origins[1-near_role] = point, other_point
        if float((origins[0]-origins[1])@self.frame_chart.left) <= 1e-6:
            raise AxialChartUnavailable('Interior paired roles disagree with original lateral order')
        endpoints = [End(v.section_id, int(not v.forward) if self.end_index == 0 else int(v.forward))
                     for v in visits]
        endpoints[near_role] = End(near_visit.section_id, near_end)
        self.original_ends = tuple(endpoints)  # Opposite token is explicitly non-incidental below.
        self.original_incidence = tuple(role == near_role for role in (0, 1))
        self.source_run_sections = tuple(v.section_id for v in visits)
        self.path = MemberInventory(tuple(Piece(path, True) for path in self.paths))
        reach = 8.*max(path.section.width for path in self.paths)
        profiles, curves, samples_by_role, intervals = [], [], [], []
        for role, (path, visit) in enumerate(zip(self.paths, visits)):
            profile = context[Key('profile', visit.section_id)].oriented(int(not visit.forward))
            if profile.half != (1 if role == 0 else -1):
                raise ValueError('Interior paired entrance lost its original half-profile orientation')
            profiles.append(profile)
            if role == near_role:
                first, last = ((0., min(reach, path.length)) if near_end == 0
                               else (max(0., path.length-reach), path.length))
            else:
                first, last = max(0., anchor-reach), min(path.length, anchor+reach)
            intervals.append((first, last))
            spacing = min(.5, context.network.settings.spacing)
            if not path.exact:
                offset = abs(path.section.offset)+.5*path.section.width
                if offset*path.curvature_bound >= 1.:
                    raise AxialChartUnavailable('Interior source offsets need regular original curvature')
                curvature = path.curvature_bound/(1.-offset*path.curvature_bound)
                speed = 1.+offset*path.curvature_bound
                if curvature > 1e-12:
                    spacing = min(spacing, math.sqrt(8.*source_tolerance/curvature)/speed,
                                  math.pi/(16.*curvature*speed))
                for control in CenterlineRange.build(path, first, last).curves:
                    _certify_increasing(3.*np.diff(control, axis=0)@self.frame_chart.axis
                        *(1. if visit.forward else -1.), coordinate='interior paired axial')
            if source_parameters is None:
                samples = path.samples(first, last, spacing, tolerance=source_tolerance)
            else:
                from .paired_source_sampling import sample_source_parameters
                samples = sample_source_parameters(path, first, last, spacing, source_tolerance,
                    source_parameters.get(path.section.id, ()), maximum_samples=self.maximum_samples,
                    protected_stations=self.protected_source_parameters.get(path.section.id, ()))
                first, last = float(samples[0]), float(samples[-1])
                intervals[-1] = first, last
            if len(samples) > self.maximum_samples:
                raise ValueError('Interior paired entrance exceeds its source sampling budget')
            samples_by_role.append(samples if visit.forward else samples[::-1])
            for side in (-1., 1.):
                xy = path.point(samples, path.section.offset
                    +side*(1. if visit.forward else -1.)*.5*path.section.width)
                curves.append(xy if visit.forward else xy[::-1])
        if abs(profiles[0].rise-profiles[1].rise)+abs(profiles[0].datum-profiles[1].datum) > 1e-9:
            raise ValueError('Interior paired entrance requires the original shared crown')
        self.field = LinearCrown(self.frame_chart, curves, datum=profiles[0].datum, rise=profiles[0].rise)
        self.source_ranges, self.source_curve_stations = tuple(intervals), tuple(samples_by_role)
        self.source_curve_tolerance = source_tolerance
        self.origin_station = 0.
        signed = [self.sign*ray.stations[[0, -1]] for ray in self.field.curves]
        self.minimum_station = max(0., max(float(value.min()) for value in signed))+1e-5
        self.length = min(float(value.max()) for value in signed)
        if self.length <= self.minimum_station+1e-6:
            raise AxialChartUnavailable('Interior paired entrance has no common finite source fibre')
        high, low = self.bounds(self.minimum_station)
        self.width, self.offset = float(high-low), 0.
        self.source, self.end = near.section, near_end
        self.origin = self.frame_chart.world(np.array((0., .5*(high+low))))
        self.direction = self.sign*self.frame_chart.axis
        # Cross sections are finalized at the chosen outward fibre. A nearby
        # separated root does not authorize filling a genuine separator gap.
        self.profile = context[Key('profile', near_visit.section_id)].oriented(near_end)

    def source_mouths(self, station):
        aliases = super().source_mouths(station)
        return tuple(replace(alias, original_incidence=self.original_incidence[role])
                     for role, alias in enumerate(aliases))


@dataclass(frozen=True)
class OneSidedPairedSpecification:
    owner: Key
    near_role: int
    group: object
    entrances: tuple
    source_intervals: tuple
    opposite_cut_keys: tuple
    proof: dict

    @classmethod
    def build(cls, context, site, local_layout, near_visits, opposite_visit, near_role, *, source_parameters=None,
              protected_source_parameters=None):
        """Prepare only explicitly supplied original Bundle/source incidence."""
        from .compounds import CompoundGroup, SampledFootprint, RimSegment, RimSpan, sample_footprint, rim_label
        from .compound_cycles import PolylineCurb
        from .terminal_curb_ancestry import original_strip, clip, validate as validate_bank
        from .intersections import half_intervals, Mouth
        from .regions import Region, RegionUnion, BoundaryLoop
        if len(near_visits) != 2 or near_role not in (0, 1):
            raise ValueError('One-sided preparation needs both adjacent original near visits')
        before, after = near_visits
        endpoints = (End(before.section_id, int(before.forward)),
                     End(after.section_id, int(not after.forward)))
        if any(endpoint not in site.ends for endpoint in endpoints):
            raise ValueError('Near source visits do not meet at the declared original junction')
        if opposite_visit.section_id in {end.section_id for end in site.ends}:
            raise ValueError('The opposite source must have no authored incidence at this junction')
        stems = [i for i, approach in enumerate(local_layout.approaches)
                 if End(approach.source.id, approach.end) not in endpoints]
        if len(stems) != 1 or len(local_layout.approaches) != 3:
            raise AxialChartUnavailable('One-sided paired preparation requires one actual three-arm source junction')
        from .source_mouths import source_mouths
        protected = {} if protected_source_parameters is None else {
            sid: set(values) for sid, values in protected_source_parameters.items()}
        parameters = {} if source_parameters is None else {sid: set(values) for sid, values in source_parameters.items()}
        # These old complete profile rows must become internal. Their exact
        # original metric stations are semantic constraints, not arbitrary
        # sampling points; otherwise a coarse strip chord can expose a tiny
        # old mouth fragment even though the intended source strip covers it.
        for approach, station in zip(local_layout.approaches, local_layout.outward_stations):
            if End(approach.source.id, approach.end) not in endpoints:
                continue
            aliases = source_mouths(approach, station)
            alias = next((a for a in aliases if a.endpoint.section_id == approach.source.id), None)
            if alias is None:
                raise AxialChartUnavailable('An old near mouth crosses another original source interval')
            protected.setdefault(alias.endpoint.section_id, set()).add(alias.station)
        for sid, values in protected.items():
            parameters.setdefault(sid, set()).update(values)
        tolerance = min(.00025, .0125*context.network.settings.curve_tolerance)
        entrances = tuple(InteriorPairedEntrance(context, visit, opposite_visit, endpoint.end,
            near_role, source_tolerance=tolerance, source_parameters=parameters,
            protected_source_parameters=protected)
            for visit, endpoint in zip(near_visits, endpoints))
        # Both evaluations of the opposite source must use identical original
        # chord partitions where their finite support overlaps.
        from .paired_source_sampling import shared_source_parameters
        partitions = shared_source_parameters(entrances)
        entrances = tuple(InteriorPairedEntrance(context, visit, opposite_visit, endpoint.end,
            near_role, source_tolerance=tolerance, source_parameters=partitions,
            protected_source_parameters=protected)
            for visit, endpoint in zip(near_visits, endpoints))
        local = sample_footprint(site.key, local_layout)
        stations = []
        for entrance in entrances:
            along = entrance.sign*entrance.frame_chart.local(local.region.xy)[:, 0]
            station = max(entrance.minimum_station, float(along.max())+.25)
            if station >= entrance.length-1e-6:
                raise AxialChartUnavailable('A one-sided source junction has no finite outward common fibre')
            if float(entrance.gap(station)) >= -1e-7:
                # The original strips may merge farther along their finite
                # body. Own only those actual strips up to a merged fibre;
                # any intervening separator remains an explicit union hole.
                previous = station
                for trial in np.linspace(station, entrance.length-1e-5, 65)[1:]:
                    if float(entrance.gap(trial)) < -1e-7:
                        low, high = previous, float(trial)
                        for _ in range(40):
                            middle = .5*(low+high)
                            if float(entrance.gap(middle)) < -1e-7:
                                high = middle
                            else:
                                low = middle
                        station = high+.25
                        break
                    previous = float(trial)
                else:
                    raise AxialChartUnavailable('One-sided original source has no merged finite outward fibre')
                if station >= entrance.length-1e-6 or float(entrance.gap(station)) >= -1e-7:
                    raise AxialChartUnavailable('One-sided merged outward fibre lacks finite original source clearance')
            entrance.profile = entrance.cross_section(station)
            stations.append(station)
        planes = tuple(a.cut_plane(s) for a, s in zip(entrances, stations))
        operands = [local.region]
        source_polygons, intervals, root_caps = [], [], []
        # Original near leads cover the old mouth rows. Their original node
        # cross-sections remain internal in the actual footprint union.
        strips = [(entrance, near_role) for entrance in entrances]
        strips.append((entrances[0], 1-near_role))
        strip_edges = {}
        for index, (entrance, role) in enumerate(strips):
            curves = entrance.field.curves[2*role:2*role+2]
            root_caps.extend((np.array((curves[0].points[at], curves[1].points[at])) for at in (0, -1)))
            points, edges = original_strip(entrance, role)
            for center, normal in planes:
                points, edges = clip(points, edges, center, normal)
            source_polygons.append(points)
            key = Key('paired_source_strip', 3*site.id+index)
            operands.append(Region(key, points))
            strip_edges[key] = entrance, edges
        loops = RegionUnion().build(operands)
        if sum(not boundary.hole for boundary in loops) != 1:
            raise AxialChartUnavailable('One-sided paired source union needs one physical exterior')
        loop = next(boundary for boundary in loops if not boundary.hole)
        for boundary in loops:
            for span in boundary.spans:
                if span.owner == site.key:
                    label = local.segments[span.segment]
                    if label.kind == 'mouth' and label.index != stems[0]:
                        raise ValueError('One-sided source union exposes an internal original near-side mouth')
        from .cycle_terminal_ports import _on_segment
        for boundary in loops:
            for first, last in zip(boundary.xy, np.roll(boundary.xy, -1, axis=0)):
                if any(_on_segment(first, cap) and _on_segment(last, cap) for cap in root_caps):
                    raise ValueError('One-sided physical union exposes an internal original source root cap')
        stem_index = stems[0]
        approaches = (*entrances, local_layout.approaches[stem_index])
        distances = (*stations, local_layout.outward_stations[stem_index])
        radial = half_intervals(context.network.settings)
        mouth_planes = (*planes, (approaches[2].point(distances[2], approaches[2].offset),
                                  approaches[2].frame(distances[2])[1]))
        tokens = []
        for a, b in zip(loop.xy, np.roll(loop.xy, -1, axis=0)):
            hits = [i for i, (center, normal) in enumerate(mouth_planes)
                    if max(abs((a-center)@normal), abs((b-center)@normal)) < 1e-7]
            if len(hits) > 1:
                raise ValueError('A one-sided physical mouth has ambiguous fibre ownership')
            tokens.append(hits[0] if hits else None)
        starts = [i for i, token in enumerate(tokens) if token is not None and tokens[i-1] != token]
        if len(starts) != 3 or {tokens[i] for i in starts} != {0, 1, 2}:
            raise ValueError('One-sided physical union did not retain exactly three complete profiles')
        def retained_curve(points, spans):
            curve = PolylineCurb(points)
            curve.source_spans = tuple(spans)
            curve.original_junction_curves = tuple(local_layout.boundaries)
            edges = []
            for i, span in enumerate(spans):
                if span.owner == site.key:
                    edge = rim_label(span, {site.key: local})
                    if edge.kind != 'curb':
                        raise ValueError('One-sided curb exposes an original junction mouth')
                else:
                    entrance, tags = strip_edges[span.owner]
                    tag = tags[span.segment]
                    if tag.kind != 'curb':
                        raise ValueError('One-sided curb exposes an internal or cut edge')
                    edge = tag.bank.subspan(span.first, span.last)
                    validate_bank(entrance, edge, points[i:i+2])
                edges.append(edge)
            curve.one_sided_edges = tuple(edges)
            return curve

        points, labels, curbs, provenance = [], [], [], []
        for start in starts:
            mouth = tokens[start]
            end = (start+1) % len(loop.xy)
            while tokens[end] == mouth:
                end = (end+1) % len(loop.xy)
            approach, station = approaches[mouth], distances[mouth]
            method = getattr(approach, 'mouth_parameters', None)
            u = method(radial, approach.width/(2*radial)) if method else np.linspace(1., 0., 2*radial+1)
            profile = approach.point(station, approach.offset+approach.width*(u-.5))
            error = max(np.linalg.norm(profile[0]-loop.xy[start]), np.linalg.norm(profile[-1]-loop.xy[end]))
            if error > 1e-7:
                raise ValueError(f'One-sided source union clips an authoritative paired profile endpoint: mouth {mouth}, error {error:g} m')
            points.extend(profile[:-1])
            labels.extend(RimSegment('mouth', mouth, float(a), float(b)) for a, b in zip(u, u[1:]))
            chain, source_spans = [profile[-1]], []
            at = end
            while tokens[at] is None:
                source_spans.append(loop.spans[at])
                at = (at+1) % len(loop.xy)
                chain.append(loop.xy[at])
            if not source_spans:
                raise ValueError('One-sided physical profiles need an explicit retained curb')
            curve = retained_curve(np.asarray(chain), source_spans)
            index = len(curbs)
            curbs.append(curve)
            points.extend(curve.points[:-1])
            labels.extend(RimSegment('curb', index, float(a/curve.length), float(b/curve.length))
                          for a, b in zip(curve.stations, curve.stations[1:]))
            provenance.append(tuple(source_spans))
        # The public layout must order its approaches with the actual CCW
        # mouth sequence, independent of original member traversal direction.
        order = [tokens[start] for start in starts]
        inverse = {old: new for new, old in enumerate(order)}
        labels = tuple(replace(label, index=inverse[label.index]) if label.kind == 'mouth' else label for label in labels)
        ordered = tuple(approaches[i] for i in order)
        outward = tuple(distances[i] for i in order)
        hole_loops, separator_pairs = [], []
        for boundary in loops:
            if not boundary.hole:
                continue
            split = len(boundary.xy)//2
            hole_spans = []
            separator_pairs.append((len(curbs), len(curbs)+1))
            for chain, spans in ((boundary.xy[:split+1], boundary.spans[:split]),
                                 (np.vstack((boundary.xy[split:], boundary.xy[0])), boundary.spans[split:])):
                curve = retained_curve(chain, spans)
                index = len(curbs)
                curbs.append(curve)
                hole_spans.extend(RimSpan(site.key, RimSegment('curb', index, float(a/curve.length), float(b/curve.length)))
                    for a, b in zip(curve.stations, curve.stations[1:]))
            hole_loops.append(BoundaryLoop(boundary.xy, tuple(hole_spans)))
        mouths = {}
        for approach, station in zip(ordered, outward):
            from .source_mouths import source_mouths
            for alias in source_mouths(approach, station):
                if alias.original_incidence:
                    mouths[alias.endpoint.section_id, alias.endpoint.end] = Mouth(alias.station, 2*radial)
        layout = SimpleNamespace(approaches=ordered, outward_stations=outward, boundaries=tuple(curbs),
            mouths=mouths, corners=local_layout.corners, rows=len(points), radial=radial,
            center=local_layout.center.copy(), guide_weights=[0.]*len(ordered),
            mesh_strategy='boundary_chart', transition=None, absorbed_transitions=(), absorbed_by=None,
            finite_scope_ports=(), finite_free_pads=(),
            original_junction_layout=local_layout, original_junction=site,
            one_sided_separator_pairs=tuple(separator_pairs),
            one_sided_source_preparation=True, source_cycle_proofs={}, source_curve_tolerance=tolerance,
            one_sided_source_polygons=tuple(source_polygons), retained_source_spans=tuple(provenance))
        layout.one_sided_original_local_polygon = local.region.xy
        layout.one_sided_original_regions = tuple(operands)
        layout.one_sided_strip_sections = tuple(entrance.paths[role].section.id for entrance, role in strips)
        layout.one_sided_original_seams_hidden = True
        footprint = SampledFootprint(Region(site.key, np.asarray(points)), labels, layout)
        final_loops = (*RegionUnion().build((footprint.region,)), *hole_loops)
        group = CompoundGroup((site.key,), (), (footprint,), final_loops, layout)
        for entrance, station in zip(entrances, stations):
            alias = entrance.source_mouths(station)[near_role]
            path = entrance.paths[near_role]
            interval = (0., alias.station) if alias.endpoint.end == 0 else (alias.station, path.length)
            intervals.append(SourceEvaluationInterval(path.section.id, *interval,
                entrance.visits[near_role].forward))
        cuts = tuple(a.source_mouths(s)[1-near_role] for a, s in zip(entrances, stations))
        first, last = sorted(alias.station for alias in cuts)
        intervals.append(SourceEvaluationInterval(opposite_visit.section_id, first, last, opposite_visit.forward))
        layout.one_sided_source_intervals = tuple(intervals)
        keys = tuple(DerivedSourcePortKey(site.key, alias.endpoint.section_id, alias.endpoint.end) for alias in cuts)
        proof = dict(original_source_topology=True, no_refitted_source_path=True,
            retained_original_local_footprint=True, opposite_original_incidence=False,
            original_near_mouths_internal=True, original_root_caps_internal=True,
            source_sampling_tolerance=tolerance, merged_profile_endpoint_error_limit=1e-7,
            near_original_ends=[(e.section_id, e.end) for e in endpoints],
            opposite_cuts=[dict(section=key.section_id, end=key.end, station=alias.station)
                           for key, alias in zip(keys, cuts)])
        proof['continuous_source_width'] = certify(group)
        return cls(site.key, near_role, group, entrances, tuple(intervals), keys, proof)

    @staticmethod
    def candidates(network):
        return ()  # Deliberately disabled until body integration qualifies.


def certify(group):
    from .updates import drain_steps
    return drain_steps(certify_steps(group))


def certify_steps(group):
    """Continuously retain the original local footprint and clipped full strips."""
    from .finite_containment import segments_on_polygon_side
    from .prepared_source_classifier import SourceProofClassifier
    from .source_strip_containment import _distance
    from .regions import _classify
    exteriors = [loop.xy for loop in group.loops if not loop.hole]
    holes = [loop.xy for loop in group.loops if loop.hole]
    if len(exteriors) != 1 or not getattr(group.physical_layout, 'one_sided_source_preparation', False):
        raise ValueError('One-sided source proof requires one certified physical exterior')
    exterior = exteriors[0]
    layout = group.physical_layout
    if not getattr(layout, 'one_sided_original_seams_hidden', False):
        # Older opt-in transports predate this certificate. Recheck the actual
        # source labels before accepting their consolidated curb descriptors.
        from .compounds import sample_footprint
        from .regions import RegionUnion
        from .cycle_terminal_ports import _on_segment
        local = sample_footprint(layout.original_junction.key, layout.original_junction_layout)
        paired = [a for a in layout.approaches if hasattr(a, 'near_role')]
        near = {a.paths[a.near_role].section.id for a in paired}
        caps = []
        for a, role in [(a, a.near_role) for a in paired] + [(paired[0], 1-paired[0].near_role)]:
            curves = a.field.curves[2*role:2*role+2]
            caps.extend(np.array((curves[0].points[at], curves[1].points[at])) for at in (0, -1))
        for boundary in RegionUnion().build(layout.one_sided_original_regions):
            for span in boundary.spans:
                if span.owner != local.region.owner:
                    continue
                label = local.segments[span.segment]
                if (label.kind == 'mouth' and
                        local.layout.approaches[label.index].source.id in near):
                    raise ValueError('One-sided source union exposes an internal original near-side mouth')
            for first, last in zip(boundary.xy, np.roll(boundary.xy, -1, axis=0)):
                if any(_on_segment(first, cap) and _on_segment(last, cap) for cap in caps):
                    raise ValueError('One-sided physical union exposes an internal original source root cap')
        layout.one_sided_original_seams_hidden = True
    projection = max((row['maximum_projection_error'] for row in
                      getattr(layout, 'source_rim_normalization', ())), default=0.)
    hole_projection = max((row['maximum_projection_error'] for row in
                           getattr(layout, 'source_hole_normalization', ())), default=0.)
    budgeted = getattr(layout, 'provenance_rim_sampling', {})
    projection += budgeted.get('exterior_projection', 0.)
    hole_projection += budgeted.get('hole_projection', 0.)
    local_band = 1e-7+projection+hole_projection
    band = .002-layout.source_curve_tolerance-projection-hole_projection
    if band <= 0.:
        raise ValueError('One-sided source errors exhaust the independent complete-width budget')
    classifier = SourceProofClassifier()
    for polygon, tolerance in ((layout.one_sided_original_local_polygon, local_band),
                               *((polygon, band) for polygon in layout.one_sided_source_polygons)):
        if len(polygon) < 3:
            raise ValueError('Finite source footprint needs a planar polygon')
        sides = np.stack((polygon, np.roll(polygon, -1, axis=0)), axis=1)
        for at in range(0, len(sides), 32):
            yield
            first, last = sides[at:at+32].transpose(1, 0, 2)
            if not np.all(segments_on_polygon_side(first, last, exterior, inside=True,
                                                  tolerance=tolerance, classifier=classifier)):
                if polygon is layout.one_sided_original_local_polygon:
                    raise ValueError('One-sided physical exterior clips the actual original local junction')
                raise ValueError('One-sided physical exterior clips an original full-width source strip')
            for hole in holes:
                if not np.all(segments_on_polygon_side(first, last, hole, inside=False,
                                                      tolerance=tolerance, classifier=classifier)):
                    raise ValueError('One-sided physical hole cuts an original source region')
        for hole in holes:
            occupied, near = _classify(hole, polygon, 1e-9)
            if np.any((occupied | (near >= 0)) & (_distance(hole, sides) > tolerance)):
                raise ValueError('One-sided physical hole lies inside an original source region')
    return dict(original_local_footprint=True, clipped_original_full_width_strips=3,
                maximum_combined_error=.002, source_chord_error=layout.source_curve_tolerance,
                boundary_band=band, holes=len(holes), original_local_projection=projection,
                physical_hole_projection=hole_projection)
