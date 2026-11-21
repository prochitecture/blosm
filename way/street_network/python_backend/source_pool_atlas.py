"""Priority atlas only for the newly admitted measured large source pools."""
from copy import copy
from dataclasses import replace


def admission(group,context):
    """Recognize new measured pools without changing any old default route.

    An attached certificate alone is insufficient: original source ownership,
    full widths, finite work counts and the atlas's source contracts are checked
    again. A copied layout keeps this read-only test from overwriting the
    accepted PLAN's certificate before the comparison.
    """
    layout=group.physical_layout
    saved=getattr(layout,'extended_source_pool',None)
    if (context is None or len(group.members)<=8 or not isinstance(saved,dict)
            or saved.get('policy')!='MEASURED_COMPLETE_ORIGINAL_SOURCE_POOL'
            or not hasattr(layout,'source_cycle_unexposed_corners')):
        return None
    from .source_cycle_atlas import eligible
    if not eligible(group,context):return None
    from .source_pool_budget import certify
    candidate=replace(group,physical_layout=copy(layout))
    fresh=certify(context,candidate)
    if saved!=fresh:
        raise ValueError('Priority source-pool atlas requires its matching fresh measured-work certificate')
    if layout.source_cycle_full_width_proof!=candidate.physical_layout.source_cycle_full_width_proof:
        raise ValueError('Priority source-pool atlas requires its matching fresh complete-width proof')
    return fresh


def build(parent,group,context,certificate):
    from .source_cycle_atlas import build as original_corridors
    result=original_corridors(parent,group,context,[dict(
        policy='NEW_MEASURED_SOURCE_POOL_PRIORITY',ordinary_site_budget=8,
        original_sites=len(group.members),old_default_input_class=False)])
    result.metadata=dict(result.metadata,measured_source_pool_priority=dict(
        admission=certificate,old_eight_site_and_raw_ordinary_routes_unchanged=True,
        original_corridor_order_unchanged=True,physical_boundary_changed=False,
        complete_source_and_stored_checks_required=True,cell_budget=parent.max_cells))
    return result
