"""Fill tiny, generated-only gaps that cannot retain the corner-radius floor.

This is a physical pavement addition, explicitly authorized for residual road
corner overlaps. Native banks, incoming profiles, medians and authored corners
cannot authorize this operation. All surviving boundaries remain identical.
"""
from copy import copy
from dataclasses import dataclass, replace
from types import SimpleNamespace
import math
import numpy as np
from .model import Key
from .compounds import SampledFootprint, RimSegment, rim_label
from .regions import Region, signed_area
from .finite_free_caps import PreparedCapCurb
from .generated_junction_curbs import GeneratedJunctionArea

SOURCE_RESIDUAL_GAP_POLICY = '1-certified-small-generated-corner-gaps'


@dataclass(frozen=True)
class ResidualRoadGapFill(GeneratedJunctionArea):
    original_loop: int
    original_area: float
    radius_floor: float
    width_upper_bound: float
    perpendicular_width_upper_bound: float
    source_parents: tuple

    def report(self):
        return dict(policy=SOURCE_RESIDUAL_GAP_POLICY, original_loop=self.original_loop,
            original_area=self.original_area, auxiliary_radius_floor=self.radius_floor,
            certified_width_upper_bound=self.width_upper_bound,
            certified_perpendicular_width_upper_bound=self.perpendicular_width_upper_bound,
            common_parent_disks_box_area_bound=self.width_upper_bound*self.perpendicular_width_upper_bound,
            original_circle_parents=self.source_parents,
            native_banks_or_profiles_erased=False, changes_original_corner_radius=False)


def _candidate(context, group, li, sources):
    from .corner_radii import CornerRadiusPolicy
    from .displaced_corner_authority import original_parent_corner
    from .retained_arc_fillets import _original_circle, OriginalAnalyticalCircleCurve
    loop = group.loops[li]
    if not 3 <= len(loop.xy) <= 512 or signed_area(loop.xy) >= 0.:
        return None
    vectors = np.roll(loop.xy, -1, axis=0)-loop.xy
    turns = vectors[:, 0]*np.roll(vectors, -1, axis=0)[:, 1]-vectors[:, 1]*np.roll(vectors, -1, axis=0)[:, 0]
    lengths = np.linalg.norm(vectors, axis=1)
    if np.any(lengths <= 1e-10) or np.any(turns > 1e-10*lengths*np.roll(lengths, -1)):
        return None
    parents = {}; sections = set(); resolver = SimpleNamespace(context=context, group=group)
    sites = {node.key: node for node in context.network.intersections}
    floor = CornerRadiusPolicy.minimum_gore_radius
    for span in loop.spans:
        label = rim_label(span, sources); source = sources[span.owner]
        if (label.kind != 'curb' or span.owner.kind != 'cycle_corner'
                or isinstance(getattr(source.layout, 'generated_junction_area', None), GeneratedJunctionArea)):
            return None
        corner = _original_circle(source, label.index)
        node = sites.get(getattr(source.layout, 'original_junction', None))
        if (corner is None or node is None or node.key not in group.members
                or corner.radius_mode != 'STYLE' or corner.radius < floor
                or any(arm.profile.half for arm in (corner.first, corner.second))
                or not original_parent_corner(resolver, node, corner)):
            return None
        token = span.owner, label.index
        if token not in parents:
            prepared = source.layout.boundaries[label.index]
            if isinstance(prepared, OriginalAnalyticalCircleCurve):
                prepared = prepared.original_prepared_curve
            parents[token] = (corner, OriginalAnalyticalCircleCurve(span.owner, label.index, prepared, corner))
        curve = parents[token][1]
        phase = curve.circle_parameter(np.array((label.first, label.last)))
        if (not np.isfinite(phase).all() or not 0. <= min(phase) <= max(phase) <= 1.
                or corner.sweep*(phase[1]-phase[0]) >= 0.):
            return None
        for arm in (corner.first, corner.second):
            sections.add(arm.source.id)
            sections.update(piece.section.id for piece in getattr(arm.path, 'pieces', ()))
        for walk in getattr(corner, 'source_contact_walks', ()):
            sections.add(walk.original.source.id)
            sections.update(piece.section.id for piece in walk.prefix)
    if len(parents) < 2 or any(sections.intersection(bundle.sections) for bundle in context.network.bundles):
        return None
    circles = tuple(value[0] for value in parents.values())
    centers = np.asarray([corner.center for corner in circles]); radii = np.asarray([corner.radius for corner in circles])
    # This common-disk bound concerns the complete analytical arcs, not just
    # their sampled chords. The clockwise convex loop is on the inside of
    # every parent disk; native or differently oriented pieces cannot enter.
    if np.any(np.linalg.norm(loop.xy[:, None, :]-centers[None, :, :], axis=2) > radii+1e-8):
        return None
    directions = vectors/lengths[:, None]
    best = None
    for axis in directions:
        axes = np.array((axis, (-axis[1], axis[0])))
        projected = centers @ axes.T
        upper = np.min(projected+radii[:, None], axis=0)-np.max(projected-radii[:, None], axis=0)
        if (np.all(upper > 0.) and min(upper) < 2.*floor
                and float(np.prod(upper)) < math.pi*floor*floor):
            key = float(np.prod(upper)), float(min(upper))
            if best is None or key < best[0]:
                best = key, tuple(sorted(map(float, upper)))
    if best is None:
        return None
    original_area = abs(signed_area(loop.xy))
    if original_area >= math.pi*floor*floor:
        return None
    parent_records = tuple(dict(owner=owner.text(), index=index, original_radius=corner.radius,
        requested_radius=corner.requested_radius, mode=corner.radius_mode,
        original_junction=getattr(sources[owner].layout, 'original_junction').text())
        for (owner, index), (corner, _) in parents.items())
    return ResidualRoadGapFill(li, original_area, floor, *best[1], parent_records)


def reconcile(planner, group, *, maximum_fills=16):
    """Add only the exact old gap polygons; keep all other loops by identity."""
    if type(maximum_fills) is not int or not 1 <= maximum_fills <= 16:
        raise ValueError('Residual road-gap fill needs its finite work allowance')
    if len(group.loops) < 2:
        return group
    sources = {source.region.owner: source for source in group.footprints}
    selected = []
    for li in range(1, len(group.loops)):
        record = _candidate(planner.context, group, li, sources)
        if record is not None:
            if len(selected) >= maximum_fills:
                raise ValueError('Residual road-gap fill exceeds its finite hole allowance')
            selected.append(record)
    if not selected:
        return group
    if getattr(group.physical_layout, 'residual_road_gap_reconciliation', None):
        raise ValueError('Residual road-gap fill repeats a prepared local construction')
    additions = []
    for i, record in enumerate(selected):
        owner = Key('residual_road_gap_fill', i)
        if owner in sources:
            raise ValueError('Residual road-gap fill repeats an original physical owner')
        xy = group.loops[record.original_loop].xy[::-1].copy()
        curve = PreparedCapCurb(np.vstack((xy, xy[0])))
        labels = tuple(RimSegment('internal', 0, float(a/curve.length), float(b/curve.length))
                       for a, b in zip(curve.stations, curve.stations[1:]))
        layout = SimpleNamespace(approaches=(), outward_stations=(), corners=(), mouths={},
            boundaries=(curve,), rows=len(xy), generated_junction_area=record)
        additions.append(SampledFootprint(Region(owner, xy), labels, layout))
    if sum(len(source.region.xy) for source in (*group.area_footprints, *additions)) > 20000:
        raise ValueError('Residual road-gap fill exceeds the original source segment budget')
    removed = {record.original_loop for record in selected}
    loops = tuple(loop for li, loop in enumerate(group.loops) if li not in removed)
    # Adding an exact complete hole polygon to this physical area has exactly
    # these boundary components. No numerical re-union or boundary relocation
    # is needed: every surviving loop, span, contact and source point is kept.
    layout = copy(group.physical_layout)
    layout.residual_road_gap_reconciliation = dict(policy=SOURCE_RESIDUAL_GAP_POLICY,
        records=tuple(record.report() for record in selected),
        construction=tuple((record.original_loop, tuple(map(tuple, group.loops[record.original_loop].xy)))
                           for record in selected),
        removed_generated_holes=len(selected), added_pavement_area=sum(record.original_area for record in selected),
        original_outer_and_other_loops_identical=True, native_source_claims_and_mouths_unchanged=True,
        prepared_owner_area_guard='PENDING_FINAL_PLAN_INVENTORY')
    return replace(group, footprints=(*group.footprints, *additions), loops=loops, physical_layout=layout)


def signature(group):
    record = getattr(group.physical_layout, 'residual_road_gap_reconciliation', None)
    return (SOURCE_RESIDUAL_GAP_POLICY, record['construction']) if record else None
