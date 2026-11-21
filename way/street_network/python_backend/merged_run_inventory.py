"""Original visible spans and actual inner curb bridges of a merged Bundle.

Three-arm source junctions split pavement spans. They are never concatenated
as an area-free member seam. Actual complete incoming profiles and published
inner curb chains delimit those spans; only an explicitly unmeshed outer field
connector may cross a junction's owned area.
"""
from dataclasses import dataclass, replace
from types import SimpleNamespace
import numpy as np

from .bundle_angular_boundaries import AngularMemberBoundary
from .bundle_angular_curves import AngularCoordinates, AngularCrown, RayCurb
from .bundle_angular_internal import SampledAngularPort
from .bundle_angular_ports import AngularPort
from .bundle_angular_routes import AngularMemberRoute
from .bundle_axis import AxialChartUnavailable
from .bundle_visible_ports import VisibleMouth
from .model import BundleMember, End, Key
from .ownership import PLAN, CompoundProduct, StandaloneProduct, surface_key
from .separator_ports import PORTS, bind_physical_port
from .source_mouths import source_mouths
from .updates import BuildContext


class MergedRunNeedsPairedJunction(AxialChartUnavailable):
    """A one-sided original junction still overlaps the opposite carriageway."""
    def __init__(self, site, endpoint, gap):
        self.site, self.endpoint, self.gap = site, endpoint, float(gap)
        super().__init__('Merged run needs paired physical preparation at original '+site.text())


def inventory_frame(context, bundle, approach):
    """Preserve polar frames; certify a separate whole-run chart for prefixes."""
    if not hasattr(approach, 'source_ranges'):
        return approach.frame_chart
    paths = [(context[Key('path', v.section_id)], v.forward)
             for member in approach.run.members for v in member.visits]
    points = np.concatenate([path.points for path, _ in paths])
    vectors = []
    for member in approach.run.members:
        first, last = member.visits[0], member.visits[-1]
        a, b = context[Key('path', first.section_id)], context[Key('path', last.section_id)]
        vectors.append(b.frame(b.length if last.forward else 0.)[0]-a.frame(0. if first.forward else a.length)[0])
    axis = np.sum(vectors, axis=0)
    if np.linalg.norm(axis) > 1e-8:
        axis /= np.linalg.norm(axis)
        left = np.array((-axis[1], axis[0]))
        origin = .5*(points.min(axis=0)+points.max(axis=0))
        extent = max(float(np.ptp((points-origin)@axis)), float(np.ptp((points-origin)@left)),
                     max(path.section.width for path, _ in paths), 1.)
        # This chooses coordinates only. The original complete paths/offsets
        # must all prove positive progress before this frame can be returned.
        # A far transverse centre avoids a nearly longitudinal terminal fibre
        # produced by a least-squares circle on a nearly axial source run.
        for factor in (16., 64.):
            frame = object.__new__(AngularCoordinates)
            frame.center = origin+factor*extent*left
            frame.scale = factor*extent
            radial = origin-frame.center
            frame.angle = float(np.arctan2(radial[1], radial[0]))
            frame.direction = 1.
            try:
                for path, forward in paths: frame.certify(path, forward=forward)
            except AxialChartUnavailable:
                continue
            return frame
    from .bundle_angular_internal import shared_frame
    frame = shared_frame(context, bundle)
    for path, forward in paths: frame.certify(path, forward=forward)
    return frame


@dataclass
class SpanEnd:
    kind: str
    endpoint: End | None
    site: object
    port: object = None
    physical_xyz: object = None
    owner: object = None
    source_station: float | None = None
    merged_entrance: object = None
    registry_key: object = None
    approach_index: int | None = None
    prepared_mouth: object = None
    original_terminal: End | None = None


@dataclass
class SourceSpan:
    role: int
    visits: tuple
    route: object
    head: SpanEnd
    tail: SpanEnd
    source_intervals: tuple = ()
    curves: tuple = ()
    transition_owners: tuple = ()
    route_interval: tuple | None = None
    curve_support: tuple = ()
    original_visits: tuple = ()
    owned_source_intervals: tuple = ()


@dataclass(frozen=True)
class InnerBridge:
    role: int
    site: object
    owner: object
    previous_end: End
    next_end: End
    physical_xyz: np.ndarray
    angular_stations: np.ndarray


def source_spans(context, bundle, approach, *, frame=None):
    """Inventory complete original visits without requiring finished surfaces."""
    run = approach.run
    frame = inventory_frame(context, bundle, approach) if frame is None else frame
    return _partition_source_spans(context,bundle,run,frame,approach.original_ends,approach)


def physical_source_spans(context,bundle,run,*,frame):
    """Partition actual original visits with no assumed merged endpoint.

    Every three-or-more-arm authored incidence remains a physical boundary.
    Hydration subsequently resolves the real published profiles; this creates
    no paired owner, derived endpoint or synthetic source-mouth alias.
    """
    if len(run.members)!=2 or len(bundle.members)!=2 or not all(m.complete for m in bundle.members):
        raise ValueError('Physical source spans need two complete original member inventories')
    for role,member in enumerate(run.members):
        visits=tuple(member.visits);original=tuple(bundle.members[role].path)
        if not visits or not any(original[i:i+len(visits)]==visits
                                 for i in range(len(original)-len(visits)+1)):
            raise ValueError('Physical source spans need an ordered original member subsequence')
    return _partition_source_spans(context,bundle,run,frame,None,None)


def _partition_source_spans(context,bundle,run,frame,merged_ends,approach):
    original = {v.section_id for member in bundle.members for v in member.path}
    used = [v.section_id for member in run.members for v in member.visits]
    if len(used) != len(set(used)) or not set(used) <= original:
        raise ValueError('Merged run repeats or invents an original member Section')
    sites = {end: j for j in (*context.network.intersections, *context.network.connections) for end in j.ends}
    result = []
    for role, member in enumerate(run.members):
        groups, current = [], []
        for visit in member.visits:
            if current:
                previous = current[-1]
                end = End(previous.section_id, int(previous.forward))
                entry = End(visit.section_id, int(not visit.forward))
                site = sites.get(end)
                if site is None or entry not in site.ends:
                    raise ValueError('Merged span lacks an explicit original Section incidence')
                if len(site.ends) > 2:
                    groups.append(tuple(current)); current = []
            current.append(visit)
        groups.append(tuple(current))
        spans = []
        for visits in groups:
            route = AngularMemberRoute(context, BundleMember(visits), role)
            route.certify(frame)
            endpoints = (End(visits[0].section_id, int(not visits[0].forward)),
                         End(visits[-1].section_id, int(visits[-1].forward)))
            ends = []
            for endpoint in endpoints:
                site = sites.get(endpoint)
                merged = merged_ends is not None and endpoint == merged_ends[role]
                if not merged and (site is None or len(site.ends) < 3):
                    raise AxialChartUnavailable('Merged visible spans need actual physical outer junction ports')
                ends.append(SpanEnd('MERGED' if merged else 'PHYSICAL', endpoint,
                    site.key if site is not None else None, merged_entrance=approach if merged else None))
            spans.append(SourceSpan(role, visits, route, *ends))
        result.append(tuple(spans))
    if {v.section_id for spans in result for span in spans for v in span.visits} != set(used):
        raise ValueError('Merged span partition dropped an original source interval')
    return tuple(result)


def sample_sources(context, spans, approach, *, max_samples=200000, additional_entrances=()):
    """Extend complete original sampling outside a prepared finite prefix."""
    if type(max_samples) is not int or not 2 <= max_samples <= 200000:
        raise ValueError('Merged source samples need a bounded original Path allocation')
    from .paired_source_sampling import verify_shared_source_phase
    entrances = tuple(a for a in (approach, *additional_entrances)
                      if hasattr(a, 'source_curve_stations'))
    verify_shared_source_phase(entrances)
    supplied = {}
    for entrance in entrances:
        for path, at in zip(entrance.paths, entrance.source_curve_stations):
            if context[Key('path', path.section.id)] is not path:
                raise ValueError('Merged source sampling lost its original Path object')
            supplied.setdefault(path.section.id, []).append(np.sort(at))
    result = {}
    for role in spans:
        for span in role:
            for visit in span.visits:
                path = context[Key('path', visit.section_id)]
                full = path.samples(0., path.length, context.network.settings.spacing,
                    tolerance=.1*context.network.settings.curve_tolerance, max_samples=max_samples)
                local = supplied.get(visit.section_id, ())
                if local:
                    keep = np.ones(len(full), dtype=bool)
                    for at in local: keep &= (full < at[0]-1e-10) | (full > at[-1]+1e-10)
                    full = np.unique(np.concatenate((full[keep], *local)))
                full.flags.writeable = False
                result[visit.section_id] = full
                if 2*sum(len(values) for values in result.values()) > max_samples:
                    raise AxialChartUnavailable('Merged run exceeds its bounded source sample allocation')
    return result


def evaluation_spans(context, bundle, approach, specifications):
    """Split evaluation intervals around paired owners without new source Ends.

    These are preparation descriptors, not accepted pavement. A derived end
    has endpoint=None and its separate physical registry key. Its original
    Path/global station survives, and no original End resolver is consulted.
    """
    from .one_sided_paired_junctions import DerivedSourcePortKey
    spans = source_spans(context, bundle, approach)
    specifications = tuple(specifications)
    result = []
    for role, own in enumerate(spans):
        leaves = []
        for span in own:
            by_source = {v.section_id: (i, v) for i, v in enumerate(span.visits)}
            claimed = []
            for specification in specifications:
                if specification.near_role == role:
                    continue
                cut_sources = {key.section_id for key in specification.opposite_cut_keys}
                for interval in specification.source_intervals:
                    if interval.section_id not in by_source or interval.section_id not in cut_sources:
                        continue
                    index, visit = by_source[interval.section_id]
                    path = context[Key('path', interval.section_id)]
                    if (interval.forward != visit.forward or interval.first < 0.
                            or interval.last > path.length+1e-8):
                        raise ValueError('Derived span cut changes its original source traversal or finite support')
                    ends = []
                    for entrance in specification.entrances:
                        layout = specification.group.physical_layout
                        at = next(i for i, a in enumerate(layout.approaches) if a is entrance)
                        distance = layout.outward_stations[at]
                        aliases = [a for a in entrance.source_mouths(distance)
                                   if not a.original_incidence and a.endpoint.section_id == interval.section_id]
                        if len(aliases) != 1:
                            raise ValueError('Derived evaluation cut lacks one explicitly non-incidental source alias')
                        alias = aliases[0]
                        key = DerivedSourcePortKey(specification.owner, interval.section_id, alias.endpoint.end)
                        if key not in specification.opposite_cut_keys:
                            raise ValueError('Derived source cut lost its separate physical registry identity')
                        end = SpanEnd('DERIVED_MERGED', None, specification.owner,
                            owner=specification.owner, source_station=alias.station,
                            merged_entrance=entrance, registry_key=key)
                        station = (alias.station if visit.forward else path.length-alias.station)
                        ends.append((float(span.route.stations[index]+station), end))
                    ends.sort(key=lambda row: row[0])
                    lower = float(span.route.stations[index]+(interval.first if visit.forward else path.length-interval.last))
                    upper = float(span.route.stations[index]+(interval.last if visit.forward else path.length-interval.first))
                    if (len(ends) != 2 or abs(ends[0][0]-lower) > 1e-8 or abs(ends[1][0]-upper) > 1e-8):
                        raise ValueError('Paired source interval does not match its two authoritative derived cuts')
                    claimed.append((lower, upper, ends[0][1], ends[1][1]))
            claimed.sort(key=lambda row: row[0])
            if any(a[1] > b[0]+1e-8 for a, b in zip(claimed, claimed[1:])):
                raise AxialChartUnavailable('Overlapping paired source owners require combined physical preparation')
            first, head = 0., span.head
            segments = []
            for lower, upper, before, after in claimed:
                if lower > first+1e-8: segments.append((first, lower, head, before))
                first, head = upper, after
            if span.route.length > first+1e-8:
                segments.append((first, span.route.length, head, span.tail))
            for lower, upper, head, tail in segments:
                visits, intervals, begin = [], [], None
                for index, visit in enumerate(span.visits):
                    path = context[Key('path', visit.section_id)]
                    base = float(span.route.stations[index])
                    a, b = max(lower, base), min(upper, base+path.length)
                    if b <= a+1e-8: continue
                    if begin is None: begin = base
                    first_global, last_global = ((a-base, b-base) if visit.forward else
                                                (path.length-(b-base), path.length-(a-base)))
                    for end, active in ((head, abs(a-lower) <= 1e-9), (tail, abs(b-upper) <= 1e-9)):
                        if (active and end.registry_key is not None
                                and end.registry_key.section_id == visit.section_id):
                            if (end is head) == visit.forward: first_global = end.source_station
                            else: last_global = end.source_station
                    first_global, last_global = np.clip((first_global, last_global), 0., path.length)
                    if first_global <= 1e-9: first_global = 0.
                    if path.length-last_global <= 1e-9: last_global = path.length
                    visits.append(visit)
                    intervals.append(dict(section=visit.section_id, forward=visit.forward,
                                          first=float(first_global), last=float(last_global)))
                route = AngularMemberRoute(context, BundleMember(tuple(visits)), role)
                leaves.append(SourceSpan(role, tuple(visits), route, head, tail, tuple(intervals),
                                         route_interval=(lower-begin, upper-begin)))
            if abs(sum(b-a for a, b, _, _ in segments)+sum(b-a for a, b, _, _ in claimed)-span.route.length) > 1e-7:
                raise ValueError('Derived evaluation spans lost original metric source coverage')
        result.append(tuple(leaves))
    return tuple(result)


def _hydrate_end(context, span, end, index):
    """Attach one authoritative physical port without reinterpreting a token."""
    from .prepared_merged_profiles import PreparedMergedMouth
    registry = context[PORTS]
    if end.registry_key is not None:
        key = end.registry_key
        if end.endpoint is not None or end.kind != 'DERIVED_MERGED':
            raise ValueError('A derived span port cannot bind an original source End')
        visible = registry.resolve_derived(key.owner, key.section_id, key.end)
        if visible is None:
            raise AxialChartUnavailable('Derived merged span needs its accepted owner-scoped physical port')
        binding = bind_physical_port(context, visible, End(key.section_id, key.end))
    else:
        if end.endpoint is None:
            raise ValueError('An original span port needs its original source End')
        visible = registry.resolve(end.endpoint)
        if visible is None:
            raise AxialChartUnavailable('Original merged span needs its accepted physical port')
        binding = bind_physical_port(context, visible, end.endpoint)
    owner, product = binding.owner, binding.product
    layout, incoming_index = binding.layout, binding.approach_index
    incoming = layout.approaches[incoming_index]
    distance = layout.outward_stations[incoming_index]
    aliases = source_mouths(incoming, distance)
    if end.registry_key is not None:
        matches = [a for a in aliases if not a.original_incidence
                   and (a.endpoint.section_id, a.endpoint.end) == (key.section_id, key.end)]
    else:
        matches = [a for a in aliases if a.original_incidence and a.endpoint == end.endpoint]
    if len(matches) != 1:
        raise ValueError('A hydrated span lost its authoritative original source alias')
    alias = matches[0]
    if (end.owner is not None and end.owner != owner
            or end.source_station is not None and abs(end.source_station-alias.station) > 1e-8):
        raise ValueError('A hydrated span changed its prepared physical owner or original station')
    if len(aliases) == 2:
        if end.merged_entrance is not None and end.merged_entrance is not incoming:
            raise ValueError('Merged span hydration needs the same freshly prepared physical entrance')
        mouth = PreparedMergedMouth.from_product(product, incoming_index)
        if incoming.end_index != index:
            raise ValueError('Merged span mouth has the wrong original traversal direction')
        end.owner, end.source_station = owner, float(alias.station)
        end.kind = 'DERIVED_MERGED' if end.registry_key is not None else 'MERGED'
        end.merged_entrance, end.approach_index = incoming, incoming_index
        end.prepared_mouth, end.physical_xyz = mouth, mouth.xyz
    else:
        if end.registry_key is not None:
            raise ValueError('A derived paired span cannot become a singleton original mouth')
        _resolve(context, span, end, index)
        end.approach_index = incoming_index


def _owned_terminal(context, span, index):
    """Move only through the accepted owner's explicit original approach route."""
    from .compounds import CompoundGroup
    from .ownership import footprint_key
    from .regions import BoundaryLoop, signed_area
    from .source_strip_containment import certify
    endpoint = span.head if index == 0 else span.tail
    registry = context[PORTS]
    if (endpoint.kind != 'PHYSICAL' or endpoint.registry_key is not None
            or endpoint.endpoint is None or registry.resolve(endpoint.endpoint) is not None):
        return span
    original = endpoint.endpoint
    declared = next((j for j in context.network.intersections if j.key == endpoint.site), None)
    if declared is None or original not in declared.ends:
        raise ValueError('Moved terminal needs its declared original physical incidence')
    owner = _physical_owner(context, endpoint.site)
    product = context.products.get(surface_key(owner))
    if product is None or product.pavement.mesh is None:
        raise AxialChartUnavailable('Moved terminal needs its accepted original physical owner')
    raw = context[footprint_key(endpoint.site)]
    approaches = [a for a in raw.approaches if (a.source.id, a.end) == (original.section_id, original.end)]
    if len(approaches) != 1 or not hasattr(approaches[0].path, 'pieces'):
        raise AxialChartUnavailable('Moved terminal lacks its explicit original approach route')
    route = approaches[0].path
    indices = range(len(span.visits)) if index == 0 else range(len(span.visits)-1, -1, -1)
    for at in indices:
        visit = span.visits[at]
        candidate = End(visit.section_id, int(not visit.forward) if index == 0 else int(visit.forward))
        visible = registry.resolve(candidate)
        if visible is None: continue
        if _physical_owner(context, visible.site.key) != owner:
            raise AxialChartUnavailable('Moved terminal reaches a different physical owner before its complete profile')
        binding = bind_physical_port(context, visible, candidate)
        layout = binding.layout
        incoming = layout.approaches[binding.approach_index]
        distance = layout.outward_stations[binding.approach_index]
        aliases = source_mouths(incoming, distance)
        if len(aliases) != 1 or not aliases[0].original_incidence or aliases[0].endpoint != candidate:
            raise AxialChartUnavailable('Moved terminal needs one complete original physical profile')
        alias = aliases[0]
        positions = [i for i, piece in enumerate(route.pieces) if piece.section.id == visit.section_id]
        if len(positions) != 1:
            raise ValueError('Moved terminal is absent from its original finite approach route')
        position = positions[0]
        expected = tuple(span.visits[:at+1]) if index == 0 else tuple(reversed(span.visits[at:]))
        if (tuple(piece.section.id for piece in route.pieces[:position+1]) != tuple(v.section_id for v in expected)
                or any(piece.path is not context[Key('path', piece.section.id)]
                       or piece.forward != (v.forward if index == 0 else not v.forward)
                       for piece, v in zip(route.pieces[:position+1], expected))):
            raise ValueError('Moved terminal cannot skip an original source incidence or replace its controlled route')
        piece = route.pieces[position]
        outward = float(route.stations[position]+(alias.station if piece.forward else piece.path.length-alias.station))
        trims = route.trims(outward)
        if list(trims)[-1] != (candidate.section_id, candidate.end) or abs(list(trims.values())[-1]-alias.station) > 1e-9:
            raise ValueError('Moved terminal profile does not match the complete original route trims')
        intervals, claims = {}, []
        original_ids = {v.section_id: v for v in span.visits}
        for (sid, end), station in trims.items():
            if sid not in original_ids:
                raise ValueError('Moved terminal crosses an unowned original source approach')
            path = context[Key('path', sid)]
            first, last = (0., station) if end == 0 else (station, path.length)
            if last-first > 1e-9: intervals[sid] = first, last
        if isinstance(product, CompoundProduct):
            group = product.group
        else:
            adjacency, _ = _boundary(product.pavement.mesh)
            first = min(adjacency); nodes, previous, current = [first], first, adjacency[first][0]
            while current != first:
                nodes.append(current)
                following = next(node for node in adjacency[current] if node != previous)
                previous, current = current, following
                if len(nodes) > len(adjacency): raise ValueError('Moved terminal owner has no simple physical rim')
            if len(nodes) != len(adjacency):
                raise ValueError('Moved standalone terminal has an unsupported second physical rim')
            xy = product.pavement.mesh.xyz[nodes, :2]
            if signed_area(xy) < 0.: xy = xy[::-1]
            group = CompoundGroup((owner,), (), (), (BoundaryLoop(xy, ()),))
        proof = certify(context, group, intervals, intervals=intervals)
        for sid, (first, last) in intervals.items():
            claims.append(dict(section=sid, forward=original_ids[sid].forward, first=float(first), last=float(last),
                               owner=owner, proof='complete_original_approach_route_prefix', containment=proof,
                               source_path=context[Key('path', sid)]))
        endpoint.endpoint, endpoint.original_terminal = candidate, original
        endpoint.owner, endpoint.source_station = owner, float(alias.station)
        chosen = slice(at, None) if index == 0 else slice(None, at+1)
        visits = tuple(span.visits[chosen])
        rows = tuple(span.source_intervals[chosen]) if span.source_intervals else ()
        return replace(span, visits=visits, route=AngularMemberRoute(context, BundleMember(visits), span.role),
                       source_intervals=rows, original_visits=span.original_visits or span.visits,
                       owned_source_intervals=(*span.owned_source_intervals, *claims))
    raise AxialChartUnavailable('Original terminal has no owner-certified visible source interval')


def hydrate_evaluation_spans(context, bundle, approach, specifications, *, spans=None,
                             frame=None, max_samples=200000):
    """Hydrate visible leaves with real mouths and finite original curb reserve.

    The active source intervals stop at the actual physical mouth aliases.
    Curves retain the complete finite Paths of those visits so a full profile
    plane and its inward bridge can cut them without extrapolation. This is
    support geometry only: callers still clip pavement to the active mouths.
    No whole-run field or one-sided separation assumption is made here.
    """
    specifications = tuple(specifications)
    spans = evaluation_spans(context, bundle, approach, specifications) if spans is None else spans
    frame = inventory_frame(context, bundle, approach) if frame is None else frame
    candidates = [a for spec in specifications for a in spec.group.physical_layout.approaches]
    candidates.extend(end.merged_entrance for own in spans for span in own for end in (span.head, span.tail))
    additional = tuple({id(a): a for a in candidates if hasattr(a, 'source_curve_stations')}.values())
    samples = sample_sources(context, spans, approach, max_samples=max_samples,
                             additional_entrances=additional)
    result, count = [], 0
    for role, own in enumerate(spans):
        leaves = []
        for original in own:
            if original.role != role:
                raise ValueError('A hydrated source span changed its original Bundle role')
            members = tuple(bundle.members[role].path)
            starts = [i for i, visit in enumerate(members) if original.visits and visit == original.visits[0]]
            if len(starts) != 1 or members[starts[0]:starts[0]+len(original.visits)] != tuple(original.visits):
                raise ValueError('Hydrated leaves must retain consecutive original Bundle visits')
            span = replace(original, head=replace(original.head), tail=replace(original.tail))
            span = _owned_terminal(context, span, 0)
            span = _owned_terminal(context, span, 1)
            _hydrate_end(context, span, span.head, 0)
            _hydrate_end(context, span, span.tail, 1)
            rows, support = [], []
            for i, visit in enumerate(span.visits):
                path = context[Key('path', visit.section_id)]
                frame.certify(path, forward=visit.forward)
                row = (dict(span.source_intervals[i]) if span.source_intervals else
                       dict(section=visit.section_id, forward=visit.forward, first=0., last=path.length))
                if row['section'] != visit.section_id or row['forward'] != visit.forward:
                    raise ValueError('A hydrated source interval changed an original visit')
                for end, active in ((span.head, i == 0), (span.tail, i == len(span.visits)-1)):
                    if not active: continue
                    sid = end.registry_key.section_id if end.registry_key is not None else end.endpoint.section_id
                    if sid != visit.section_id:
                        raise AxialChartUnavailable('A moved physical mouth needs its explicit original visit inventory')
                    head = end is span.head
                    row['first' if head == visit.forward else 'last'] = end.source_station
                if not 0. <= row['first'] < row['last'] <= path.length:
                    raise AxialChartUnavailable('Hydrated merged span has no finite visible original interval')
                rows.append(row)
                support.append(dict(section=visit.section_id, forward=visit.forward, first=0., last=path.length))
            boundary = _SourceBoundary(context, BundleMember(span.visits), frame, samples,
                                       join_tolerance=span.route.join_tolerance)
            span.curves = tuple(boundary.points(side) for side in (-1., 1.))
            span.transition_owners = boundary.owners
            span.source_intervals, span.curve_support = tuple(rows), tuple(support)
            first, last = rows[0], rows[-1]
            span.route_interval = (first['first'] if first['forward'] else span.route.pieces[0].path.length-first['last'],
                span.route.stations[-2]+(last['last'] if last['forward'] else span.route.pieces[-1].path.length-last['first']))
            for curve in span.curves: curve.flags.writeable = False
            count += sum(len(curve) for curve in span.curves)
            if count > max_samples:
                raise AxialChartUnavailable('Hydrated source leaves exceed their finite curb sample allocation')
            leaves.append(span)
        result.append(tuple(leaves))
    return tuple(result)


def hydrate_physical_spans(context,bundle,spans,*,frame,max_samples=200000):
    """Resolve only real junction profiles before sampling original source bands.

    A published two-alias mouth may resolve to its actual prepared interface.
    It contributes its compulsory source phase before full-source sampling.
    Unpaired intermediate junctions keep their original physical ends.
    """
    prepared=[]
    for own in spans:
        rows=[]
        for original in own:
            if any(end.kind!='PHYSICAL' or end.endpoint is None or end.registry_key is not None
                   or end.prepared_mouth is not None or end.merged_entrance is not None
                   for end in (original.head,original.tail)):
                raise ValueError('Physical source hydration needs explicit original physical ends')
            span=replace(original,head=replace(original.head),tail=replace(original.tail))
            span=_owned_terminal(context,span,0);span=_owned_terminal(context,span,1)
            _hydrate_end(context,span,span.head,0);_hydrate_end(context,span,span.tail,1)
            rows.append(span)
        prepared.append(tuple(rows))
    return hydrate_evaluation_spans(context,bundle,None,(),spans=tuple(prepared),frame=frame,
                                    max_samples=max_samples)


def _boundary(mesh):
    edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    if len(counts) and counts.max() > 2:
        raise ValueError('Actual junction boundary is nonmanifold')
    rim = unique[counts == 1]
    adjacency = {}
    for a, b in rim:
        adjacency.setdefault(int(a), []).append(int(b)); adjacency.setdefault(int(b), []).append(int(a))
    if any(len(neighbors) != 2 for neighbors in adjacency.values()):
        raise ValueError('Actual junction rim has an unresolved boundary incidence')
    return adjacency, {tuple(edge) for edge in rim}


def _physical_owner(context, site):
    plan = context.products.get(PLAN)
    return (plan.owner(site) if plan is not None else None) or site


def _standalone_mouth(product, port, endpoint, owner):
    mesh = product.pavement.mesh
    if mesh is None:
        raise AxialChartUnavailable('Merged span needs hydrated actual junction mouth samples')
    adjacency, edges = _boundary(mesh)
    indices = np.array(sorted(adjacency), dtype=int)
    outer, inner = port.coordinates(np.array((0., 1.)))[0]
    delta = inner-outer; squared = float(delta@delta)
    offset = mesh.xyz[indices, :2]-outer
    u = offset@delta/squared
    error = np.linalg.norm(offset-u[:, None]*delta, axis=1)
    selected = (u >= -1e-10) & (u <= 1.+1e-10) & (error <= 1e-8)
    indices, u = indices[selected], np.clip(u[selected], 0., 1.)
    order = np.argsort(u); indices, u = indices[order], u[order]
    if len(u) < 2 or abs(u[0]) > 1e-10 or abs(u[-1]-1.) > 1e-10 or np.any(np.diff(u) <= 1e-10):
        raise AxialChartUnavailable('Actual standalone mesh lacks a complete original profile')
    if any(tuple(sorted((int(a), int(b)))) not in edges for a, b in zip(indices, indices[1:])):
        raise ValueError('Declared incoming profile is not one actual junction boundary chain')
    u[[0, -1]] = 0., 1.
    xyz = mesh.xyz[indices].copy()
    expected = port.coordinates(u)[2]
    if product.exaggeration:
        physical = port.profile.datum+(xyz[:, 2]-port.profile.datum)/product.exaggeration
        if np.max(abs(physical-expected)) > 1e-8:
            raise ValueError('Actual standalone profile disagrees with its original physical crown')
    xyz[:, 2] = expected
    u.flags.writeable = xyz.flags.writeable = False
    return VisibleMouth(endpoint, owner, port.station, u, xyz)


def _resolve(context, span, end, index):
    visit = span.visits[0 if index == 0 else -1]
    visible = context[PORTS].resolve(end.endpoint)
    if visible is None:
        raise AxialChartUnavailable('Original merged span end is consumed; its visible interval needs preparation')
    binding = bind_physical_port(context, visible, end.endpoint)
    layout = binding.layout
    incoming = layout.approaches[binding.approach_index]
    distance = layout.outward_stations[binding.approach_index]
    alias = next((a for a in source_mouths(incoming, distance) if a.endpoint == end.endpoint), None)
    if alias is None:
        raise ValueError('Merged span lost its original source-mouth alias')
    end.owner = binding.owner
    end.source_station = alias.station
    if end.kind == 'MERGED':
        return
    if len(source_mouths(incoming, distance)) != 1:
        raise AxialChartUnavailable('An intermediate paired mouth needs its merged span interface policy')
    path = context[Key('path', visit.section_id)]
    profile = context[Key('profile', visit.section_id)].oriented(int(not visit.forward))
    product = context[surface_key(end.owner)]
    mouth = VisibleMouth.from_compound(context, end.endpoint, span.role, visit.forward)
    if mouth is None:
        if not isinstance(product, StandaloneProduct):
            raise ValueError('Original merged span has no accepted physical junction product')
        provisional = AngularPort(span.role, index, path, visit.forward, alias.station,
            layout.mouths[end.endpoint.section_id, end.endpoint.end].transverse_intervals, profile)
        mouth = _standalone_mouth(product, provisional, end.endpoint, end.owner)
    end.port = SampledAngularPort(span.role, index, path, visit.forward, profile, mouth)
    end.physical_xyz = mouth.physical_xyz


def _standalone_cap(context, product, first, last):
    mesh = product.pavement.mesh
    if mesh is None:
        raise AxialChartUnavailable('Merged inner curb needs hydrated published junction boundary samples')
    adjacency, edges = _boundary(mesh)
    def node(point):
        hits = [i for i in adjacency if np.linalg.norm(mesh.xyz[i, :2]-point[:2]) <= 1e-8]
        if len(hits) != 1:
            raise ValueError('Declared inner peak is missing or ambiguous on the actual junction rim')
        return hits[0]
    start, finish = node(first.physical_xyz[-1]), node(last.physical_xyz[-1])
    mouths = set()
    for approach, station in zip(product.layout.approaches, product.layout.outward_stations):
        a, b = approach.point(station, approach.offset+np.array((-.5, .5))*approach.width)
        delta = b-a; squared = float(delta@delta)
        for x, y in edges:
            points = mesh.xyz[[x, y], :2]-a
            u = points@delta/squared
            if np.all((u >= -1e-10) & (u <= 1.+1e-10)) and np.max(np.linalg.norm(points-u[:, None]*delta, axis=1)) <= 1e-8:
                mouths.add((x, y))
    candidates = []
    for neighbor in adjacency[start]:
        chain, previous, current = [start], start, neighbor
        while True:
            if tuple(sorted((previous, current))) in mouths:
                break
            chain.append(current)
            if current == finish:
                candidates.append(chain); break
            following = next(i for i in adjacency[current] if i != previous)
            previous, current = current, following
            if current == start or len(chain) > len(adjacency):
                break
    if len(candidates) != 1:
        raise AxialChartUnavailable('Actual inner curb must be one published chain without another mouth')
    xyz = mesh.xyz[candidates[0]].copy()
    if product.exaggeration:
        xyz[:, 2] = first.port.profile.datum+(xyz[:, 2]-first.port.profile.datum)/product.exaggeration
    else:
        if abs(first.physical_xyz[-1, 2]-last.physical_xyz[-1, 2]) > 1e-10:
            raise AxialChartUnavailable('Zero-exaggeration standalone cap needs one original plateau boundary')
        xyz[:, 2] = first.physical_xyz[-1, 2]
    return xyz


def _inner_bridge(context, role, previous, following, frame):
    a, b = previous.tail, following.head
    if a.site != b.site or a.owner != b.owner or a.kind != 'PHYSICAL' or b.kind != 'PHYSICAL':
        raise ValueError('Merged inner bridge lacks one original physical junction owner')
    product = context[surface_key(a.owner)]
    if isinstance(product, CompoundProduct):
        if product.physical_rims is None:
            raise ValueError('Merged inner curb requires its accepted physical rim')
        physical = replace(product, rims=product.physical_rims, exaggeration=1.)
        xyz = physical.cap((a.endpoint.section_id, a.endpoint.end), (b.endpoint.section_id, b.endpoint.end))
    elif isinstance(product, StandaloneProduct):
        xyz = _standalone_cap(context, product, a, b)
    else:
        raise ValueError('Merged inner bridge has no accepted physical junction')
    if (np.linalg.norm(xyz[0]-a.physical_xyz[-1]) > 1e-8
            or np.linalg.norm(xyz[-1]-b.physical_xyz[-1]) > 1e-8):
        raise ValueError('Actual inner curb does not meet both authoritative profile peaks')
    stations = frame.local(xyz[:, :2])[:, 0]
    if np.any(np.diff(stations) <= 1e-9):
        raise AxialChartUnavailable('Actual inner curb is not monotone in this merged run chart')
    xyz, stations = xyz.copy(), stations.copy()
    xyz.flags.writeable = stations.flags.writeable = False
    return InnerBridge(role, a.site, a.owner, a.endpoint, b.endpoint, xyz, stations)


class _SourceBoundary(AngularMemberBoundary):
    def __init__(self, context, member, frame, source_samples, **kwargs):
        self.source_samples = source_samples
        super().__init__(context, member, frame, **kwargs)

    def _path_boundary(self, path, traversal, side):
        original = getattr(path, 'path', path)
        first = getattr(path, 'first', 0.); last = first+path.length
        self._certify_path(path, traversal)
        samples = self.source_samples[original.section.id]
        samples = np.r_[first, samples[(samples > first+1e-10) & (samples < last-1e-10)], last]
        q = side*(1 if traversal.forward else -1)
        points = original.point(samples, original.section.offset+q*.5*original.section.width)
        return points if traversal.forward else points[::-1]


@dataclass
class MergedRunInventory:
    bundle: object
    run: object
    frame: object
    field: object
    spans_by_role: tuple
    inner_bridges_by_role: tuple
    source_sections: tuple
    source_samples: dict
    proofs: dict
    computational_outer_bridges: tuple
    transition_owners: tuple = ()
    entrance: object = None

    @property
    def rays(self):
        return self.field.curves

    @classmethod
    def build(cls, context, bundle, approach, *, max_samples=200000):
        from .bundle_ends import PathRange
        from .bundle_angular_open_islands import gap_intervals
        frame = inventory_frame(context, bundle, approach)
        spans = source_spans(context, bundle, approach, frame=frame)
        run = approach.run
        source_samples = sample_sources(context, spans, approach, max_samples=max_samples)
        for role in spans:
            for span in role:
                _resolve(context, span, span.head, 0); _resolve(context, span, span.tail, 1)
        # A one-sided original port cannot be connected to a merged body by
        # pretending its full source profile is separated. Check against the
        # explicitly opposite member before constructing any field connector.
        for role, own in enumerate(spans):
            other = run.members[1-role]
            rays = []
            for visit in other.visits:
                path = context[Key('path', visit.section_id)]
                samples = source_samples[visit.section_id]
                q = (1 if 1-role == 0 else -1)*(1 if visit.forward else -1)
                points = path.point(samples, path.section.offset+q*.5*path.section.width)
                rays.append(RayCurb(frame, points if visit.forward else points[::-1]))
            for before, after in zip(own, own[1:]):
                for end in (before.tail, after.head):
                    station, transverse = frame.local(end.physical_xyz[-1, :2])
                    candidates = [ray for ray in rays if ray.stations[0] <= station <= ray.stations[-1]]
                    if len(candidates) != 1:
                        raise AxialChartUnavailable('Opposite original curb needs its physical junction interval')
                    opposite = float(frame.local(candidates[0](station))[1])
                    gap = float(transverse-opposite) if role == 0 else float(opposite-transverse)
                    if gap <= 1e-7:
                        raise MergedRunNeedsPairedJunction(end.site, end.endpoint, gap)
        bridges = tuple(tuple(_inner_bridge(context, role, a, b, frame) for a, b in zip(own, own[1:]))
                        for role, own in enumerate(spans))
        curves, outer_records, transition_owners = [], [], set()
        for role, own in enumerate(spans):
            parts = [[], []]
            for index, span in enumerate(own):
                products = dict(context.products)
                intervals = []
                for i, visit in enumerate(span.visits):
                    path = context[Key('path', visit.section_id)]; low, high = 0., path.length
                    for end, active in ((span.head, i == 0), (span.tail, i == len(span.visits)-1)):
                        if active and end.kind == 'PHYSICAL':
                            if end.endpoint.end == 0: low = max(low, end.source_station)
                            else: high = min(high, end.source_station)
                    if high <= low+1e-7:
                        raise AxialChartUnavailable('A merged pavement span is completely owned by its junctions')
                    intervals.append(dict(section=visit.section_id, forward=visit.forward, first=low, last=high))
                    if low > 0. or high < path.length: products[Key('path', visit.section_id)] = PathRange(path, low, high)
                span.source_intervals = tuple(intervals)
                view = BuildContext(context.network, products, context.patterns)
                boundary = _SourceBoundary(view, BundleMember(span.visits), frame, source_samples,
                    join_tolerance=span.route.join_tolerance)
                span.transition_owners = boundary.owners
                transition_owners.update(boundary.owners)
                values = tuple(boundary.points(side) for side in (-1., 1.))
                outer_index, inner_index = ((0, 1) if role == 0 else (1, 0))
                for end, at in ((span.head, 0), (span.tail, -1)):
                    if end.kind != 'PHYSICAL':continue
                    for side, node in ((outer_index, 0), (inner_index, -1)):
                        if np.linalg.norm(values[side][at]-end.physical_xyz[node, :2]) > 1e-8:
                            raise ValueError('Prepared source curb differs from its actual incoming profile endpoint')
                        values[side][at] = end.physical_xyz[node, :2]
                span.curves = values
                for side in (0, 1):
                    if parts[side]:
                        bridge = bridges[role][index-1]
                        if side == inner_index:
                            joined = bridge.physical_xyz[:, :2]
                        else:
                            joined = np.array((parts[side][-1][-1], values[side][0]))
                            outer_records.append(dict(role=role, site=bridge.site.text(),
                                policy='COMPUTATIONAL_ONLY_UNMESHED_JUNCTION_FIELD', xy=joined.tolist()))
                        if np.linalg.norm(parts[side][-1][-1]-joined[0]) > 1e-8 or np.linalg.norm(joined[-1]-values[side][0]) > 1e-8:
                            raise ValueError('Merged field lost an actual prepared curb endpoint')
                        parts[side].append(joined[1:])
                    parts[side].append(values[side] if not parts[side] else values[side][1:])
            curves.extend(np.concatenate(part) for part in parts)
        if sum(len(curve) for curve in curves) > max_samples:
            raise AxialChartUnavailable('Merged run field exceeds its finite source sample budget')
        profile = approach.field
        field = AngularCrown(frame, curves, datum=profile.datum, rise=profile.rise)
        domain = SimpleNamespace(frame=frame, rays=field.curves, field=field)
        for role in bridges:
            for bridge in role:
                intervals, roots = gap_intervals(domain, bridge.angular_stations[0], bridge.angular_stations[-1])
                if roots or not intervals or any(not positive for _, _, positive in intervals):
                    raise MergedRunNeedsPairedJunction(bridge.site, bridge.previous_end,
                        float(np.min(field.edges(bridge.angular_stations)[1]-field.edges(bridge.angular_stations)[2])))
        proofs = dict(source_sections_preserved=True, actual_physical_ports=True, actual_inner_curbs=True,
            same_original_frame=frame is approach.frame_chart, independent_run_frame=frame is not approach.frame_chart,
            actual_entrance_field_preserved=True, source_sampling_phase_preserved=True,
            one_sided_ports_separated=True, source_join_errors=[s.route.maximum_curb_discrepancy for own in spans for s in own])
        return cls(bundle, run, frame, field, spans, bridges, run.source_sections, source_samples,
                   proofs, tuple(outer_records), tuple(sorted(transition_owners)), approach)
