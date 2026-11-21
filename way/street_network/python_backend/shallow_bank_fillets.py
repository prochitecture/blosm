"""Recover a fitted arc whose sagitta falls below the boolean union's welding distance."""
from dataclasses import replace
import math
import numpy as np

from .compounds import rim_label
from .regions import BoundaryIncidenceError, RegionUnion
from .updates import drain_steps
from .source_bank_reconciliation import NativeBankWitness, SourceBankDiscovery


def _edges(loop):
    return {(tuple(a), tuple(b), span) for a, b, span
            in zip(loop.xy, np.roll(loop.xy, -1, axis=0), loop.spans)}


def candidate_union(context, group, sources, fits):
    return drain_steps(candidate_union_steps(context, group, sources, fits))


def candidate_union_steps(context, group, sources, fits):
    areas = tuple(source.region for source in sources if not source.region.boundary_only)
    union = RegionUnion(max_segments=20000)
    try:
        loops = yield from union.build_steps(areas)
        return loops, None
    except BoundaryIncidenceError as error:
        failure = error
    if len(fits) != 1:
        raise failure
    proposal, fitted = fits[0]
    corner = fitted.corner
    if proposal.kind != 'DIRECT_ORIGINAL_BANK_SWITCH' or not 0. < abs(corner.sweep) < math.pi:
        raise failure
    generated = sources[-1]
    record = getattr(generated.layout, 'source_bank_fillet', None)
    if record is None or record.fitted is not fitted or record.proposal is not proposal:
        raise failure
    sagitta = 2. * corner.radius * math.sin(.25 * corner.sweep)**2
    precision = 128. * np.spacing(max(1., float(np.abs(corner.center).max()), corner.radius))
    if sagitta >= union.epsilon or sagitta/16. <= precision:
        raise failure
    # Keep the original union authoritative for every edge outside the fitted
    # arc. A smaller epsilon is allowed only to resolve its two true contacts.
    base = yield from union.build_steps(s.region for s in group.area_footprints)
    if len(base) != 1 or base[0].hole:
        raise failure
    loops = yield from RegionUnion(epsilon=sagitta/4., max_segments=union.max_segments).build_steps(areas)
    check = yield from RegionUnion(epsilon=sagitta/16., max_segments=union.max_segments).build_steps(areas)
    if (len(loops) != 1 or len(check) != 1 or loops[0].hole or check[0].hole
            or loops[0].spans != check[0].spans or not np.array_equal(loops[0].xy, check[0].xy)):
        raise failure
    before, after = _edges(base[0]), _edges(loops[0])
    removed, added = before-after, after-before
    if not removed or not added or any(span.owner != generated.region.owner for _, _, span in added):
        raise failure
    discovery = SourceBankDiscovery(context, replace(group, loops=base))
    banks = (proposal.first, proposal.last)
    for _, _, span in removed:
        ends = [discovery.describe(span, tail) for tail in (False, True)]
        if not any(all(isinstance(end, NativeBankWitness) and discovery.same_bank(end, bank)
                       for end in ends) for bank in banks):
            raise failure
    # Require an ordered copy of every emitted segment of the physical arc.
    # Interior construction seams and partially buried arcs are ineligible.
    inventory = {source.region.owner: source for source in sources}
    expected = [(i, label) for i, label in enumerate(generated.segments) if label.kind == 'curb']
    by_segment = {span.segment: (a, b, span) for a, b, span in added}
    if len(added) != len(expected) or len(by_segment) != len(added) or not expected:
        raise failure
    bank_points = [arm.point(at, arm.offset + side*.5*arm.width) for arm, at, side in
                   zip((corner.first, corner.second), corner.stations, (-1., 1.))]
    arc = []
    for i, label in expected:
        if i not in by_segment:
            raise failure
        a, b, span = by_segment[i]
        if rim_label(span, inventory) != label or arc and arc[-1] != a:
            raise failure
        if not arc:
            arc.append(a)
        arc.append(b)
    if not np.array_equal([arc[0], arc[-1]], bank_points):
        raise failure
    return loops, dict(owner=generated.region.owner.text(), sagitta=sagitta,
                       exact_contact_points=tuple(map(tuple, bank_points)),
                       replaced_bank_spans=len(removed), complete_arc_spans=len(added),
                       base_welding_distance=union.epsilon)
