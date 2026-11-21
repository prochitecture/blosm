"""Source-only admission for a side owner coupled to a Bundle terminal.

This closes original physical ownership before any road or junction mesh is
built. A short original continuation is a candidate-selection rule only. Two
independent positive-overlap witnesses admit it; the ordinary compound builder
then owns all finite mouth, width, profile, incidence and boundary proofs.
"""
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np

from .model import End, Key, Traversal
from .routes import RouteTopology
from .updates import drain_steps


COUPLED_OWNER_POLICY = '2-visible-source-side-terminal-coupling-dependencies'


@dataclass(frozen=True)
class SourceSideTerminal:
    bundle: object
    side: object
    terminal: object
    continuation: tuple
    opposite: object
    end: int
    selection_reach: float


def _owners(groups):
    result = {}
    for group in groups:
        for site in group.members:
            if site in result:
                raise ValueError('Coupled source admission needs distinct existing original owners')
            result[site] = group
    return result


def _members(site, owners):
    group = owners.get(site)
    return tuple(group.members) if group is not None else (site,)


def _endpoint(visit, end):
    return End(visit.section_id, end if visit.forward else 1-end)


def candidates(context, bundle, owners):
    """Only a unique declared member suffix/prefix can reach a common owner."""
    if (len(bundle.members) != 2 or bundle.interior_section_ids
            or not all(member.complete for member in bundle.members)):
        return ()
    network = context.network
    topology = RouteTopology(network)
    physical = {site.key: site for site in network.intersections if len(site.ends) > 2}
    if any(network.section_index[sid].scope_ends or network.section_index[sid].selection_ends
           for sid in bundle.sections):
        return ()
    reach = 4.*max(network.section_index[sid].width for sid in bundle.sections)
    terminals = {}
    for end in (0, 1):
        sites = [topology.incidents.get(_endpoint(member.path[0 if end == 0 else -1], end))
                 for member in bundle.members]
        if (any(site not in physical for site in sites)
                or _members(sites[0], owners) != _members(sites[1], owners)):
            continue
        terminals[end] = sites[0]
    result = []
    for role, member in enumerate(bundle.members):
        other = bundle.members[1-role]
        for index, (before, after) in enumerate(zip(member.path, member.path[1:])):
            left, right = _endpoint(before, 1), _endpoint(after, 0)
            site = topology.incidents.get(left)
            if site not in physical or topology.incidents.get(right) != site:
                continue
            original = physical[site]
            # One real side road, not a crossing or an opposite member End.
            if (len(original.ends) != 3 or any(e.section_id in
                    {v.section_id for v in other.path} for e in original.ends)):
                continue
            for end, terminal in sorted(terminals.items()):
                if site in _members(terminal, owners):
                    continue
                visits = (tuple(member.path[index+1:]) if end else
                    tuple(Traversal(v.section_id, not v.forward) for v in member.path[:index+1][::-1]))
                if not visits or len(visits) > 4:
                    continue
                route = topology.trace(_endpoint(visits[0], 0))
                if (route.traversals != visits or route.terminal_scope or route.terminal_selection
                        or route.terminal_site not in _members(terminal, owners)):
                    continue
                length = sum(context[Key('path', visit.section_id)].length for visit in visits)
                if length > reach:
                    continue
                result.append(SourceSideTerminal(bundle, site, terminal, visits,
                    other.path[0 if end == 0 else -1], end, reach))
    # No choice between two potential terminal owners for the same side.
    counts = {}
    for item in result:counts[item.side] = counts.get(item.side, 0)+1
    return tuple(item for item in result if counts[item.side] == 1)


def visible_intervals(context, bundle, groups, layouts):
    """Original station ranges excluding all actual owner prefixes and bodies."""
    bounds = {sid: [0., context[Key('path', sid)].length] for sid in bundle.sections}
    consumed = {claim.section_id for group in groups for claim in group.claims}
    owners, seen = _owners(groups), set()
    for site, original in layouts.items():
        group = owners.get(site)
        layout = (group.physical_layout if group is not None and group.physical_layout is not None
                  else original)
        if id(layout) in seen:continue
        seen.add(id(layout))
        for (sid, end), mouth in getattr(layout, 'mouths', {}).items():
            if sid not in bounds:continue
            station = float(mouth.station)
            if not 0. <= station <= context[Key('path', sid)].length:
                raise ValueError('Coupled source admission lost an original finite mouth station')
            if end == 0:bounds[sid][0] = max(bounds[sid][0], station)
            else:bounds[sid][1] = min(bounds[sid][1], station)
    spans = []
    for role, member in enumerate(bundle.members):
        values = []
        for visit in member.path:
            first, last = bounds[visit.section_id]
            if visit.section_id in consumed or last-first <= 1e-8:
                continue
            row = dict(section=visit.section_id, forward=visit.forward, first=first, last=last)
            values.append(SimpleNamespace(role=role, visits=(visit,), source_intervals=(row,)))
        spans.append(tuple(values))
    return tuple(spans)


def analytic_body_overlap(context, bundle, spans, *, maximum_samples=20000, maximum_pairs=256):
    """Sampled chords propose; original analytic curb inversions prove overlap."""
    from .bundle_angular_internal import shared_frame
    from .bundle_angular_curves import AngularCrown
    from .bundle_axis import AxialChartUnavailable
    from .gapped_runtime_admission import visible_overlap
    from .gapped_source_validity import certify_visible_centerlines
    if (type(maximum_samples) is not int or not 2 <= maximum_samples <= 20000
            or type(maximum_pairs) is not int or not 1 <= maximum_pairs <= 256):
        raise ValueError('Coupled source overlap needs bounded aggregate samples and source pairs')
    if any(not values for values in spans):return None
    try:
        frame = shared_frame(context, bundle)
        for values in spans:
            for span in values:
                visit = span.visits[0]
                frame.certify(context[Key('path', visit.section_id)], forward=visit.forward)
    except AxialChartUnavailable:
        return None
    tolerance = min(.001, .025*context.network.settings.curve_tolerance)
    sampled, samples = {}, 0
    for role, values in enumerate(spans):
        for span in values:
            visit, row = span.visits[0], span.source_intervals[0]
            path = context[Key('path', visit.section_id)]
            remaining = maximum_samples-samples
            if remaining < 2:
                raise ValueError('Coupled source overlap exceeds its finite source sample budget')
            at = path.samples(row['first'], row['last'], .5, tolerance=tolerance, max_samples=remaining)
            samples += len(at)
            if samples > maximum_samples:
                raise ValueError('Coupled source overlap exceeds its finite source sample budget')
            if not visit.forward:at = at[::-1]
            sign = 1. if visit.forward else -1.
            sampled[role, visit.section_id] = tuple(path.point(at,
                path.section.offset+.5*path.section.width*side*sign) for side in (-1., 1.))
    pairs = 0
    for first in spans[0]:
        for last in spans[1]:
            if pairs >= maximum_pairs:
                raise ValueError('Coupled source overlap exceeds its finite original pair budget')
            pairs += 1
            curves = (*sampled[0, first.visits[0].section_id],
                      *sampled[1, last.visits[0].section_id])
            a, b = (frame.local(points)[:, 0] for points in (curves[1], curves[2]))
            if max(a[0], b[0]) >= min(a[-1], b[-1]):continue
            try:
                field = AngularCrown(frame, curves,
                    datum=context[Key('profile', first.visits[0].section_id)].datum,
                    rise=context[Key('profile', first.visits[0].section_id)].rise)
                witness = visible_overlap(context, SimpleNamespace(frame=frame, field=field,
                    spans_by_role=((first,), (last,))))
            except AxialChartUnavailable:
                # A proposal chart can be unavailable even though the
                # original finite sources themselves remain valid. This is
                # still discovery: no closure has been admitted or built.
                continue
            if witness is not None:
                witness['centerline_separation'] = certify_visible_centerlines(context, spans)
                witness['source_samples'] = samples
                witness['original_source_pairs_tested'] = pairs
                witness['maximum_original_source_pairs'] = maximum_pairs
                # Each existing analytic inverter considers at most64 proposed
                # fibres. The finite pair bound also bounds unsuccessful work.
                witness['maximum_analytic_witness_candidates'] = 64*maximum_pairs
                return witness
    return None


def side_footprint_overlap(context, request, spans, groups, layouts):
    """Positive area beyond both certified chord allowances, never a grass hole."""
    from .compounds import CompoundGroup, sample_footprint
    from .regions import RegionUnion
    from .sampling import AnalyticalBoundarySampler
    from .source_mouth_clearance import _positive_overlap
    from .mesh import cross
    owners = _owners(groups)
    group = owners.get(request.side)
    if group is None:
        footprint = sample_footprint(request.side, layouts[request.side])
        group = CompoundGroup((request.side,), (), (footprint,), RegionUnion().build((footprint.region,)))
    # The narrow witness belongs to pavement, not an unclassified hole interior.
    if len(group.loops) != 1 or group.loops[0].hole:return None
    tolerance = min(.001, .025*context.network.settings.curve_tolerance)
    refined = AnalyticalBoundarySampler(tolerance=tolerance, spacing=.5, max_samples=20000).build(group)
    if len(refined.loops) != 1 or refined.loops[0].hole:return None
    source = next((span for values in spans for span in values
                   if span.visits[0] == request.opposite), None)
    if source is None:return None
    row = source.source_intervals[0]
    path = context[Key('path', request.opposite.section_id)]
    remaining = 20000-sum(len(loop.xy) for loop in refined.loops)
    if remaining < 2:
        raise ValueError('Coupled side overlap exceeds its aggregate source sample budget')
    at = path.samples(row['first'], row['last'], .5, tolerance=tolerance, max_samples=remaining)
    high = path.point(at, path.section.offset+.5*path.section.width)
    low = path.point(at, path.section.offset-.5*path.section.width)
    polygon = refined.loops[0].xy
    band = 2.*tolerance+1e-7
    for index in range(len(at)-1):
        cell = np.array((high[index], high[index+1], low[index+1], low[index]))
        if np.any(cell.max(0) < polygon.min(0)) or np.any(cell.min(0) > polygon.max(0)):continue
        edge = np.roll(cell, -1, axis=0)-cell
        turns = cross(edge, np.roll(edge, -1, axis=0))
        if not (np.all(turns > 1e-12) or np.all(turns < -1e-12)):
            raise ValueError('Coupled source overlap needs unfolded original finite strip cells')
        if _positive_overlap(cell, polygon, tolerance=band):
            return dict(original_side=request.side.text(), side_owner_members=[s.text() for s in group.members],
                opposite_section=path.section.id, visible_original_interval=[row['first'], row['last']],
                original_cell_stations=[float(at[index]), float(at[index+1])],
                source_cell=cell.tolist(), positive_overlap_band=band,
                certified_curve_tolerance=tolerance, side_boundary_segments=len(polygon),
                completed_mesh_used=False, grass_holes_excluded=True)
    return None


def _close(planner, members, owners):
    return drain_steps(_close_steps(planner, members, owners))


def _close_steps(planner, members, owners):
    from .compound_preparation import CompoundExpansionRequired
    from .source_pool_budget import MAXIMUM_ORIGINAL_SITES, certify_steps
    history = []
    while True:
        members = set(member for site in members for member in _members(site, owners))
        if any(site not in planner.sites for site in members):
            raise ValueError('Coupled source closure needs an available actual original physical owner')
        if len(members) > MAXIMUM_ORIGINAL_SITES:
            raise ValueError('Coupled source closure exceeds its existing original source-pool budget')
        token = tuple(sorted(members))
        try:
            group = yield from planner.exterior_steps(token)
        except CompoundExpansionRequired as required:
            expanded = members | set(required.neighbors)
            if expanded == members:
                raise ValueError('Coupled source closure must add a new original physical neighbor') from required
            history.append(dict(members=[s.text() for s in token],
                original_neighbors=[s.text() for s in required.neighbors]))
            members = expanded
            continue
        if len(members) > planner.maximum_members:
            yield from certify_steps(planner.context, group)
        return group, history


def prepare(context, junctions, layouts, groups, *, excluded_bundle_ids=()):
    return drain_steps(prepare_steps(context, junctions, layouts, groups,
        excluded_bundle_ids=excluded_bundle_ids))


def prepare_steps(context, junctions, layouts, groups, *, excluded_bundle_ids=()):
    """Admit only positively coupled source owners; preserve all old successes."""
    from .compound_preparation import DeferredCompoundPlanner
    groups = tuple(groups)
    excluded = set(excluded_bundle_ids)
    planner = DeferredCompoundPlanner(context, junctions, layouts)
    joins = 0
    join_budget = max(0, len(layouts)-1)
    for bundle in sorted(context.network.bundles, key=lambda item:item.id):
        if bundle.id in excluded:continue
        while True:
            owners = _owners(groups)
            requests = candidates(context, bundle, owners)
            if len(requests) > 64:
                raise ValueError('Coupled source admission exceeds its finite original candidate budget')
            if not requests:break
            spans = visible_intervals(context, bundle, groups, layouts)
            body, admitted = None, False
            for request in requests:
                # Skip only discovery. Once both overlap witnesses admit a
                # source closure, normal geometry/source errors remain fatal.
                side = side_footprint_overlap(context, request, spans, groups, layouts)
                if side is None:continue
                if body is None:body = analytic_body_overlap(context, bundle, spans)
                if body is None:break
                seed = set(_members(request.side, owners)+_members(request.terminal, owners))
                group, history = yield from _close_steps(planner, seed, owners)
                covered = set(group.members)
                if any(covered.intersection(old.members) and not set(old.members) <= covered for old in groups):
                    raise ValueError('Coupled source closure must preserve every complete existing original owner')
                components = {tuple(_members(site, owners)) for site in covered}
                added_joins = len(components)-1
                if added_joins < 1 or joins+added_joins > join_budget:
                    raise ValueError('Coupled source closure must reduce its finite original owner forest')
                joins += added_joins
                proof = dict(policy=COUPLED_OWNER_POLICY, bundle=bundle.id,
                    actual_analytic_visible_body_overlap=body, actual_side_pavement_overlap=side,
                    original_continuation=[dict(section=v.section_id, forward=v.forward) for v in request.continuation],
                    original_continuation_length=sum(context[Key('path', v.section_id)].length for v in request.continuation),
                    declared_terminal_owner_members=[s.text() for s in _members(request.terminal, owners)],
                    candidate_selection_reach=request.selection_reach, selection_reach_is_not_geometry_proof=True,
                    source_expansions=history, final_members=[s.text() for s in group.members],
                    original_owner_forest_joins=added_joins, total_original_owner_forest_joins=joins,
                    maximum_original_owner_forest_joins=join_budget,
                    completed_products_used=False, original_source_geometry_changed=False)
                previous = tuple(record for old in groups if covered.intersection(old.members)
                    for record in getattr(old.physical_layout, 'coupled_source_owner_proofs', ()))
                group.physical_layout.coupled_source_owner_proofs = (*previous, proof)
                inherited = {sid for old in groups if covered.intersection(old.members)
                    for sid in getattr(old.physical_layout, 'coupled_source_dependencies', ())}
                group.physical_layout.coupled_source_dependencies = tuple(sorted(inherited | set(bundle.sections)))
                groups = tuple(old for old in groups if not covered.intersection(old.members))+(group,)
                # Rediscover even within this Bundle: another independent
                # terminal can remain coupled after the first finite closure.
                admitted = True
                break
            if not admitted:break
    return groups
