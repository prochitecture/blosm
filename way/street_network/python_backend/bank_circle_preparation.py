"""Round a clipped original corner against the visible bank of a consumed Section."""
from copy import copy
from dataclasses import dataclass, replace
import math
from types import SimpleNamespace
import numpy as np

from .bank_circle_contacts import BankCircleFillet, ContactSampling, OriginalCornerArc
from .compound_ports import CompoundPortPlanner
from .compounds import RimSegment, SampledFootprint, rim_label
from .cycle_dimensional_retry import certify_layout
from .displaced_corner_authority import original_parent_corner
from .finite_free_caps import PreparedCapCurb
from .intersections import Approach
from .model import Key
from .regions import Region, RegionUnion, signed_area
from .updates import drain_steps
from .retained_arc_fillets import CircularBoundaryGeometry, _mouth_points, _topology
from .sampling import AnalyticalBoundarySampler
from .source_bank_reconciliation import (NativeBankWitness, SourceBankDiscovery, SourceBankExposureError,
                                         SourceBankReconciliationError, fit_exposed_switch, visible_bank_stations)
from .source_circular_intervals import PreparedCircularInterval
from .source_owned_connector_fillets import NativeCurbContact
from .source_parent_chain_fillets import ParentChainUnavailable, offset_solution, radius_candidates
from .source_parent_chain_preparation import _future_guard, _refine_strip, ordered_curb_closure


@dataclass(frozen=True)
class VisibleCorner:
    parent: object
    index: int
    arc: OriginalCornerArc
    phases: tuple
    bank: NativeBankWitness
    stations: tuple
    outward: tuple


def visible_corner(context, group, proposal):
    if proposal.kind != 'BANK_TO_CLIPPED_ORIGINAL_CORNER_EXIT':
        return None
    discovery = SourceBankDiscovery(context, group)
    loop = group.loops[proposal.loop]
    at = proposal.vertex
    span = loop.spans[at]
    source = discovery.sources[span.owner]
    label = rim_label(span, discovery.sources)
    if source.region.owner.kind != 'cycle_corner' or label.kind != 'curb' or len(source.layout.corners) != 1:
        return None
    corner, = source.layout.corners
    prepared = source.layout.boundaries[label.index]
    declared = getattr(prepared, 'original_circle_interval', None)
    if not isinstance(declared, PreparedCircularInterval) or declared.prepared_bounds != (0., 1.):
        return None
    declared.validate(prepared, corner)
    site = getattr(source.layout, 'original_junction', None)
    node = next((v for v in context.network.intersections if v.key == site), None)
    if (node is None or site not in group.members or corner.center is None or corner.sweep >= 0.
            or not original_parent_corner(SimpleNamespace(context=context, group=group), node, corner)
            or corner.radius_mode != 'STYLE' or proposal.mode != 'STYLE'
            or any(mode != 'STYLE' for _, _, mode in proposal.sectors)):
        return None
    if any(hasattr(arm.path, 'mouth_piece') or getattr(arm.path, 'absorbed_sites', ())
           for arm in (corner.first, corner.second)):
        return None
    bank = proposal.first
    if bank.owner.kind != 'cycle_strip' or any(o.endpoint.section_id != bank.path.section.id for o in bank.origins):
        return None

    count = len(loop.spans)
    domain = visible_bank_stations(discovery,proposal,incoming=True)
    if domain is None:
        return None
    start, stop = domain
    native = tuple(sorted((start, stop)))
    outward = (start, stop) if bank.side == 1 else (bank.path.length-start, bank.path.length-stop)
    last = at
    for _ in range(count):
        following = (last+1) % count
        next_span = loop.spans[following]
        next_label = rim_label(next_span, discovery.sources)
        if next_span.owner != span.owner or (next_label.kind,next_label.index) != ('curb',label.index):
            break
        last = following
    else:
        raise SourceBankReconciliationError('A single original circle occupies a complete junction rim')
    tail = rim_label(loop.spans[last], discovery.sources)
    if abs(tail.last-1.) > 1e-9 or not 0. < label.first < tail.last <= 1.:
        return None
    view = OriginalCornerArc(span.owner, label.index, prepared, corner)
    return VisibleCorner(source, label.index, view, (label.first, tail.last), bank, native, outward)


def fit(context, visible, proposal):
    bank = visible.bank
    path = bank.path
    arm = Approach(path, int(bank.side == -1), context[Key('profile', path.section.id)])
    offset = path.section.offset + .5*path.section.width*bank.side
    low, high = visible.stations
    if path.exact and np.any((low < path.stations) & (path.stations < high)):
        raise SourceBankReconciliationError('A clipped-corner contact crosses an exact-polyline miter')
    curvature = 0. if path.exact else path.visible_curvature_bound(low, high)
    for radius in radius_candidates(proposal.radius, proposal.mode):
        if (radius/visible.arc.corner.radius >= 1.-1e-10 or (abs(offset)+radius)*curvature >= 1.-1e-10):
            continue
        try:
            return radius, offset_solution(visible.arc, visible.phases, arm, visible.outward, radius, bank_first=True)
        except ParentChainUnavailable:
            continue
    raise SourceBankReconciliationError('No allowed radius fits the visible Section bank and original circle')


def refine_parent(visible, contact, sweep, arc_count, *, crossing=None):
    """Insert a tangent contact and refine the interval between it and the crossing."""
    arc = visible.arc
    start = visible.phases[0] if crossing is None else crossing
    step = .25 * abs(sweep / arc.corner.sweep) / arc_count
    intervals = max(2, math.ceil(abs(contact-start)/step))
    source = visible.parent
    if len(source.segments) + intervals > 20000:
        raise SourceBankReconciliationError('Circular contact exceeds the original vertex budget')
    cuts = np.sort(np.linspace(start, contact, intervals+1)[1:])
    points, labels, references = [], [], []
    for point, label in zip(source.region.xy, source.segments):
        points.append(point)
        if (label.kind,label.index) != ('curb',visible.index):
            labels.append(label)
            continue
        references.append(len(points)-1)
        previous = label.first
        for value in cuts[np.searchsorted(cuts, label.first, side='right'):np.searchsorted(cuts, label.last)]:
            labels.append(replace(label, first=previous, last=value))
            previous = value
            references.append(len(points))
            points.append(arc(value))
        labels.append(replace(label, first=previous))
    points = np.asarray(points)
    view = copy(arc)
    view.contact_reference_points = points[[*references,(references[-1]+1) % len(points)]]
    view.contact_reference_points.flags.writeable = False
    angle = max(np.diff(arc.original_reference_parameters))*abs(arc.corner.sweep)
    sagitta = 2.*arc.corner.radius*math.sin(.25*angle)**2
    curves = list(source.layout.boundaries)
    curves[visible.index] = ContactSampling(view, min(.00025,sagitta), contact)
    layout = copy(source.layout)
    layout.boundaries, layout.rows = tuple(curves), len(points)
    region = Region(source.region.owner, points)
    RegionUnion().build((region,))
    return SampledFootprint(region, tuple(labels), layout)


def _bank_knot(strip, point, side):
    matches = []
    for index, curve in enumerate(strip.layout.boundaries):
        for i in np.flatnonzero(np.all(curve.points[:-1] == point, axis=1)):
            if curve.source_coordinates.sides[i] == side:
                matches.append((index, float(curve.stations[i]/curve.length)))
    if len(matches) != 1:
        raise SourceBankReconciliationError('Circular join needs one identifiable original bank contact')
    return matches[0]


def _union(areas):
    return drain_steps(_union_steps(areas))


def _union_steps(areas):
    return RegionUnion(max_segments=20000).build_steps(
        a.region for a in areas if not getattr(a.region,'boundary_only',False))


def closure(areas, first, last, intermediate=()):
    return drain_steps(closure_steps(areas, first, last, intermediate))


def closure_steps(areas, first, last, intermediate=()):
    sources = {area.region.owner: area for area in areas}
    try:
        loops = yield from _union_steps(areas)
        points, ancestry = ordered_curb_closure(loops, sources, first, last, intermediate)
    except ParentChainUnavailable as error:
        raise SourceBankReconciliationError(str(error)) from error
    return PreparedCapCurb(points), ancestry


def validate(group):
    records = getattr(group.physical_layout,'native_original_circle_fillets',())
    if not records:
        return group
    sources = {a.region.owner:a for a in group.footprints}
    for record in records:
        labels = [rim_label(span,sources) for loop in group.loops for span in loop.spans if span.owner == record.owner]
        covered = 0.
        for label in sorted(labels, key=lambda label:label.first):
            if (label.kind,label.index) != ('curb',0) or abs(label.first-covered) > 1e-9 or label.last <= label.first:
                raise SourceBankReconciliationError('Circular join has a gap or hidden closure on the physical rim')
            covered = label.last
        if abs(covered-1.) > 1e-9:
            raise SourceBankReconciliationError('Physical rim clips the circular join before its tangent contact')
    return group


def insert(planner, original, proposal, visible):
    context = planner.context
    settings = context.network.settings
    tolerance = min(.001, .05*settings.curve_tolerance)
    radius, solved = fit(context, visible, proposal)
    corner = CircularBoundaryGeometry(solved['center'], radius, solved['start_angle'], solved['sweep'])
    bank = visible.bank
    distance = solved['bank_outward_contact']
    station = distance if bank.side == 1 else bank.path.length-distance
    contact = solved['poly_contact']
    length = radius*abs(corner.sweep)
    preview_step = min(.2, settings.spacing, math.sqrt(8.*.00025*radius))
    parent = refine_parent(visible, contact, corner.sweep, max(2, math.ceil(length/preview_step)))
    sources = {s.region.owner:s for s in original.footprints}
    strip, first, _ = _refine_strip(sources[bank.owner], bank.path, station, bank.side)
    index, parameter = _bank_knot(strip, first, bank.side)
    last = visible.arc(contact)
    areas = tuple(strip if s.region.owner == bank.owner else parent if s is visible.parent else s
                  for s in original.footprints)
    internal, ancestry = closure(areas, (bank.owner,index,parameter,first),
                                 (parent.region.owner,visible.index,contact,last))
    nominal = np.asarray(corner.center) + radius*np.column_stack((
        np.cos([corner.start_angle, corner.start_angle+corner.sweep]),
        np.sin([corner.start_angle, corner.start_angle+corner.sweep])))
    precision = 32.*np.spacing(max(1., radius, float(np.max(np.abs(corner.center)))))
    if np.max(np.linalg.norm(nominal-np.asarray((first,last)), axis=1)) > precision:
        raise SourceBankReconciliationError('Circular join fails to meet its two analytical tangent contacts')
    records = getattr(original.physical_layout,'native_original_circle_fillets',())
    owner = Key('native_original_circle_fillet', len(records))
    if owner in sources:
        raise SourceBankReconciliationError('Circular join collides with an existing physical owner')
    step = min(.2, settings.spacing, math.sqrt(8.*tolerance*radius))
    count = max(2, math.ceil(length/step))
    if count + len(internal.points) + sum(len(a.region.xy) for a in areas) > 20000:
        raise SourceBankReconciliationError('Circular join exceeds the original vertex budget')
    knots = tuple(np.linspace(0., 1., count+1))
    record = BankCircleFillet(owner, corner, tuple(first), tuple(last),
        NativeCurbContact(bank.path,bank.side,station,bank.intervals), bank.owner,index,parameter,
        parent.region.owner,visible.index,visible.arc.original_prepared_curve,visible.arc.corner,
        contact,visible.phases,visible.stations,proposal.sectors,ancestry,knots,tolerance)
    polygon = np.vstack((record(knots)[:-1],internal.points[:-1]))
    if signed_area(polygon) <= 1e-10:
        raise SourceBankReconciliationError('Bank/circle join has no positive pavement area')
    labels = [RimSegment('curb',0,a,b) for a,b in zip(knots,knots[1:])]
    u = internal.stations/internal.length
    labels.extend(RimSegment('internal',1,a,b) for a,b in zip(u,u[1:]))
    layout = SimpleNamespace(approaches=(), outward_stations=(), mouths={}, corners=(),
                             boundaries=(record,internal), rows=len(polygon), generated_junction_area=record)
    wedge = SampledFootprint(Region(owner,polygon),tuple(labels),layout)
    RegionUnion().build((wedge.region,))
    areas = (*areas,wedge)
    physical = copy(original.physical_layout)
    physical.native_original_circle_fillets = (*records,record)
    group = replace(original, footprints=areas, physical_layout=physical, loops=_union(areas))
    consumed = {claim.section_id for claim in group.claims}
    validate(group)
    _topology(original, group, tolerance)
    certify_layout(context, group, consumed)
    sampler = AnalyticalBoundarySampler(tolerance=settings.curve_tolerance, spacing=settings.spacing, max_samples=20000)
    group = sampler.build(group)
    validate(group)
    _topology(original, group, tolerance)
    physical.source_cycle_full_width_proof = certify_layout(context, group, consumed)
    expected, actual = _mouth_points(original), _mouth_points(group)
    if (actual.keys() != expected.keys() or any(not np.array_equal(xy,actual[k]) for k,xy in expected.items())
            or CompoundPortPlanner.incomplete_ports(group) or group.members != original.members
            or group.claims != original.claims or group.physical_layout.mouths != original.physical_layout.mouths):
        raise SourceBankReconciliationError('Bank/circle join changed a complete mouth or original ownership')
    _future_guard(planner, group, record)
    return group


def _key(proposal):
    return proposal.site, proposal.kind, proposal.first.token, proposal.last.token


def recover_steps(reconciler, failure):
    """Keep successful original bank fits and retry only clipped-circle failures."""
    from .original_bank_handoffs import sample
    from .source_bank_preparation import _apply_steps, _record
    from .source_bank_reconciliation import SOURCE_BANK_POLICY
    planner = reconciler.planner
    context = planner.context
    base = sample(planner, reconciler.original, reconciler.partitions)
    selected = SourceBankDiscovery(context,base).discover()
    by_key = {_key(p):p for p in selected}
    if len(by_key) != len(selected) or len(selected) > reconciler.maximum_fillets:
        raise SourceBankReconciliationError('Bank/circle retry needs a bounded set of distinct source pairs')
    rejected = set()
    while True:
        key = _key(failure.proposal)
        if (key in rejected or key not in by_key or len(rejected) >= 16
                or visible_corner(context,base,by_key[key]) is None):
            raise failure
        rejected.add(key)
        fits = [(p,fit_exposed_switch(context,base,p)) for p in selected if _key(p) not in rejected]
        if not fits:
            current, proofs = base, []
            break
        try:
            current, _, proof = yield from _apply_steps(planner,base,reconciler.partitions,fits)
            proofs = [proof]
            break
        except SourceBankExposureError as next_failure:
            failure = next_failure
    used = set()
    while True:
        remaining = SourceBankDiscovery(context,current).discover()
        if not remaining:
            break
        proposal = next((p for p in remaining if _key(p) in rejected-used), None)
        if proposal is None:
            raise SourceBankReconciliationError('Original bank switches remain after the circular join')
        visible = visible_corner(context,current,proposal)
        if visible is None:
            raise SourceBankReconciliationError('The surviving bank fits changed the original circular parent')
        current = insert(planner,current,proposal,visible)
        used.add(_key(proposal))
    current.physical_layout.source_bank_reconciliation = dict(
        policy=SOURCE_BANK_POLICY, records=[_record(p,fit) for p,fit in fits], proofs=proofs,
        original_corners_unchanged=True, original_profiles_and_mouths_unchanged=True,
        prepared_owner_area_guard='PENDING_FINAL_PLAN_INVENTORY')
    return current


def corridor_chart(parent, group, context):
    from .updates import drain_steps
    return drain_steps(corridor_chart_steps(parent, group, context))


def corridor_chart_steps(parent, group, context):
    """Use one wholly owned undivided Section to lay out computational cuts."""
    from .paired_interval_atlas import _checks
    from .partition import CompoundPartitioner
    from .source_bank_reconciliation import owned_ranges
    from .source_cycle_atlas import OriginalCorridorCuts, SourceCycleEnergyChart
    ranges = owned_ranges(context, group)
    claimed = {claim.section_id for claim in group.claims}
    bundled = {sid for bundle in context.network.bundles for sid in bundle.sections}
    options = [(context[Key('path',sid)].length, sid) for sid in claimed-bundled
               if ranges.get(sid) == ((0.,context[Key('path',sid)].length),)]
    if not options:
        return None
    _, sid = max(options, key=lambda row:(row[0],-row[1]))
    policy = OriginalCorridorCuts(sid, spacing=parent.spacing, source_corridors=True, hole_bridges=True,
        computational_extensions=True, existing_vertices=True, source_ranges=ranges,
        cadence_multiplier=2, maximum_unpaired_span_widths=8.)
    chart = SourceCycleEnergyChart(parent.patterns, spacing=parent.spacing, max_cells=parent.max_cells,
        quality_iterations=min(100,parent.quality_iterations), four_sided=True, recover_grass=True,
        partitioner=CompoundPartitioner(policy=policy,require_visibility=False))
    result = yield from chart.build_steps(group,context)
    current = validate(result.partition.group)
    if current.members != group.members or current.claims != group.claims or owned_ranges(context,current) != ranges:
        raise SourceBankReconciliationError('Computational cuts for a circular join changed physical ownership')
    original, actual = _mouth_points(group), _mouth_points(current)
    if original.keys() != actual.keys() or any(not np.array_equal(point,actual[key]) for key,point in original.items()):
        raise SourceBankReconciliationError('Computational cuts for a circular join changed an incoming profile')
    width = certify_layout(context,current,claimed)
    checks = _checks(parent,result)
    height = parent.source_height_admission(current,result) if parent.source_height_admission else None
    result.metadata = dict(result.metadata, bank_circle_corridor=dict(section=sid, source_ranges=ranges,
        full_width=width, stored_checks=checks, height_admission=height))
    return result
