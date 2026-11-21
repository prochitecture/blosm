"""Measured resource admission for larger, certified original-source pools."""
import numpy as np
from .model import Key
from .updates import drain_steps


MAXIMUM_ORIGINAL_SITES=16
MAXIMUM_SOURCE_SEGMENTS=20000
MAXIMUM_SOURCE_CELLS=20000
MAXIMUM_SOURCE_OPERANDS=256
MAXIMUM_EXTERNAL_PROFILES=16
MAXIMUM_HOLES=8


def certify(context,group):
    return drain_steps(certify_steps(context,group))


def certify_steps(context,group):
    """The ordinary eight-site path stays unchanged; larger pools prove work."""
    if not 8<len(group.members)<=MAXIMUM_ORIGINAL_SITES:
        raise ValueError('Extended source pool needs nine to sixteen original physical owners')
    physical={site.key for site in context.network.intersections if len(site.ends)>2}
    if len(set(group.members))!=len(group.members) or not set(group.members)<=physical:
        raise ValueError('Extended source pool lost an original physical owner')
    layout=group.physical_layout
    if (layout is None or not hasattr(layout,'source_cycle_unexposed_corners')
            or not group.claims or sum(not loop.hole for loop in group.loops)!=1):
        raise ValueError('Extended pool requires a certified complete original source-cycle exterior')
    if len(group.footprints)>MAXIMUM_SOURCE_OPERANDS:
        raise ValueError('Extended source pool exceeds its finite operand budget')
    segments=sum(len(source.region.xy) for source in group.footprints)
    boundary=sum(len(loop.xy) for loop in group.loops)
    if max(segments,boundary)>MAXIMUM_SOURCE_SEGMENTS:
        raise ValueError('Extended source pool exceeds its original finite segment budget')
    if len(layout.approaches)>MAXIMUM_EXTERNAL_PROFILES or sum(loop.hole for loop in group.loops)>MAXIMUM_HOLES:
        raise ValueError('Extended source pool exceeds its profile or physical-hole budget')
    consumed=set()
    for claim in group.claims:
        path=context[Key('path',claim.section_id)]
        if (claim.first_owner not in group.members or claim.last_owner not in group.members
                or abs(claim.first_trim-path.length)>1e-8 or abs(claim.last_trim)>1e-8
                or path.section.scope_ends or path.section.selection_ends):
            raise ValueError('Extended source pool requires complete original source ownership without scope cuts')
        consumed.add(claim.section_id)
    from .cycle_dimensional_retry import certify_layout_steps
    proof=yield from certify_layout_steps(context,group,consumed)
    if proof['total_cells']>MAXIMUM_SOURCE_CELLS or proof['maximum_error']>.002:
        raise ValueError('Extended source pool exceeds its complete-width proof budget')
    layout.source_cycle_full_width_proof=proof
    xy=np.vstack([source.region.xy for source in group.footprints])
    report=dict(policy='MEASURED_COMPLETE_ORIGINAL_SOURCE_POOL',ordinary_site_budget=8,
        original_sites=len(group.members),maximum_original_sites=MAXIMUM_ORIGINAL_SITES,
        consumed_sections=sorted(consumed),source_operands=len(group.footprints),
        maximum_source_operands=MAXIMUM_SOURCE_OPERANDS,source_segments=segments,
        boundary_segments=boundary,maximum_source_segments=MAXIMUM_SOURCE_SEGMENTS,
        complete_width_cells=proof['total_cells'],maximum_complete_width_cells=MAXIMUM_SOURCE_CELLS,
        maximum_source_error=proof['maximum_error'],complete_external_profiles=len(layout.approaches),
        maximum_external_profiles=MAXIMUM_EXTERNAL_PROFILES,physical_holes=sum(loop.hole for loop in group.loops),
        maximum_holes=MAXIMUM_HOLES,extent=np.ptp(xy,axis=0).tolist(),
        source_inputs_unchanged=True,downstream_quad_budget=50000,
        downstream_surface_and_stored_checks_required=True)
    layout.extended_source_pool=report
    return report
