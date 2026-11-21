"""Conservative evidence that visible original carriageways overlap.

The body field proposes stations only. Admission inverts the original analytic
offset curbs inside the actual visible source intervals, on complete fibres
strictly beyond both original cross-section cuts. Owned junction gaps and
prefixes cannot supply an overlap witness.
"""
from dataclasses import dataclass
import math
import numpy as np
from .model import Key


@dataclass(frozen=True)
class _VisibleRange:
    role: int
    path: object
    section: int
    forward: bool
    first: float
    last: float
    lower: float
    upper: float
    inner: float
    outer: float


def _ranges(context, inventory, margin):
    frame = inventory.frame
    result = []
    if len(inventory.spans_by_role) != 2:
        raise ValueError('Visible overlap requires both original Bundle source roles')
    for role, spans in enumerate(inventory.spans_by_role):
        records = []
        for span in spans:
            if span.role != role or len(span.visits) != len(span.source_intervals):
                raise ValueError('Visible overlap lost its actual ordered original source intervals')
            for visit, row in zip(span.visits, span.source_intervals):
                if row['section'] != visit.section_id or row['forward'] != visit.forward:
                    raise ValueError('Visible overlap changed its original source traversal')
                path = context[Key('path', visit.section_id)]
                first, last = float(row['first']), float(row['last'])
                if not 0. <= first < last <= path.length:
                    raise ValueError('Visible overlap cannot extrapolate its original finite source interval')
                frame.certify(path, forward=visit.forward)
                coordinate = (1. if role == 0 else -1.)*(1. if visit.forward else -1.)
                if context[Key('profile', visit.section_id)].peak != coordinate:
                    raise ValueError('Visible overlap requires the original inward half profile')
                inner = path.section.offset+.5*path.section.width*coordinate
                outer = path.section.offset-.5*path.section.width*coordinate
                start, finish = (first, last) if visit.forward else (last, first)
                # An original cross-section is affine in offset. Its angular
                # coordinate is monotone on that finite straight segment, so
                # both endpoints bound the ENTIRE transverse source cut.
                head = frame.local(path.point(start, np.array((outer, inner))))[:, 0]
                tail = frame.local(path.point(finish, np.array((outer, inner))))[:, 0]
                lower, upper = float(head.max())+margin, float(tail.min())-margin
                if lower < upper:
                    records.append(_VisibleRange(role, path, visit.section_id, visit.forward,
                        first, last, lower, upper, inner, outer))
        records.sort(key=lambda item: (item.lower, item.upper, item.section))
        if any(a.upper > b.lower+margin for a, b in zip(records, records[1:])):
            raise ValueError('Visible overlap has ambiguous original source interval ownership')
        result.append(tuple(records))
    return tuple(result)


def _intersection(frame, source, station, offset, margin):
    """Bounded monotone inversion on the original analytical offset curve."""
    first, last = source.first, source.last
    direction = 1. if source.forward else -1.
    endpoints = frame.local(source.path.point(np.array((first, last)), offset))[:, 0]
    values = direction*(endpoints-station)
    if values[0] >= 0. or values[1] <= 0.:
        raise ValueError('Visible overlap witness leaves its original analytical curb domain')
    for _ in range(64):
        middle = .5*(first+last)
        if middle == first or middle == last:
            break
        value = direction*(float(frame.local(source.path.point(middle, offset))[0])-station)
        if value < 0.:first = middle
        else:last = middle
    source_station = .5*(first+last)
    point = source.path.point(source_station, offset)
    local = frame.local(point)
    residual = float(abs(local[0]-station))
    if residual > margin:
        raise ValueError('Original analytical curb inversion did not resolve its common fibre')
    return dict(source_station=source_station, offset=offset, point=point.tolist(),
                transverse=float(local[1]), inversion_residual=residual)


def visible_overlap(context, inventory):
    """Return one certified original-PV overlap witness, otherwise None.

    Failure to find a witness preserves the old separated runtime policy; it
    does not certify separation. Full source/mesh checks still own publication.
    """
    tolerance = float(context.network.settings.curve_tolerance)
    if not math.isfinite(tolerance) or tolerance <= 0.:
        raise ValueError('Visible overlap needs the authored source curve allowance')
    frame = inventory.frame
    scale = max(1., float(frame.scale), float(abs(frame.center).max()))
    margin = max(1e-7, 512.*np.finfo(float).eps*scale)
    allowance = 2.*tolerance
    ranges = _ranges(context, inventory, margin)
    rays = inventory.field.curves[1:3]
    lower = max(float(ray.stations[0]) for ray in rays)
    upper = min(float(ray.stations[-1]) for ray in rays)
    first_index = last_index = candidates = 0
    while first_index < len(ranges[0]) and last_index < len(ranges[1]):
        a, b = ranges[0][first_index], ranges[1][last_index]
        first, last = max(a.lower, b.lower, lower), min(a.upper, b.upper, upper)
        if first < last:
            knots = np.unique(np.concatenate((np.asarray((first, last)),
                *(ray.stations[(ray.stations > first) & (ray.stations < last)] for ray in rays))))
            nodes = np.unique(np.r_[knots, .5*(knots[:-1]+knots[1:])])
            edges = inventory.field.edges(nodes)
            depth = edges[2]-edges[1]
            proposed = np.flatnonzero(depth > allowance+margin)
            proposed = proposed[np.argsort(-depth[proposed], kind='stable')]
            for index in proposed:
                if candidates >= 64:return None
                candidates += 1
                station = float(nodes[index])
                actual = []
                for source in (a, b):
                    actual.append(dict(section=source.section, forward=source.forward,
                        visible_interval=[source.first, source.last],
                        whole_transverse_fibre_range=[source.lower, source.upper],
                        outer=_intersection(frame, source, station, source.outer, margin),
                        inner=_intersection(frame, source, station, source.inner, margin)))
                ao, ai = actual[0]['outer']['transverse'], actual[0]['inner']['transverse']
                bo, bi = actual[1]['outer']['transverse'], actual[1]['inner']['transverse']
                if not ao > ai or not bi > bo:
                    raise ValueError('Original source fibres lost their positive carriageway widths')
                overlap = min(ao, bi)-max(ai, bo)
                if overlap > allowance+margin:
                    return dict(station=station, overlap_depth=overlap,
                        combined_source_allowance=allowance, source_curve_allowances=[tolerance, tolerance],
                        numerical_margin=margin, original_analytic_curbs=actual,
                        analytical_curbs_inverted=True, complete_transverse_fibres=True,
                        visible_original_intervals=True, owned_junction_gaps_excluded=True,
                        owned_prefixes_excluded=True, wholly_owned_query_spans_excluded=True,
                        proposal_candidates=candidates, sampled_field_is_proposal_only=True)
        if a.upper <= b.upper:first_index += 1
        else:last_index += 1
    return None
