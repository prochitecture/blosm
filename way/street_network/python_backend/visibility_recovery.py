"""Source-certified visibility cuts after a radial partition has no kernel."""
import numpy as np
from .compound_chart_recovery import _finish_steps, recover_steps, recover_longitudinal_steps
from .compound_ports import CompoundPortPlanner
from .partition import BoundaryOwnershipCuts, CompoundPartitioner, CompoundQuadChart
from .sampling import AnalyticalBoundarySampler
from .source_strip_containment import certify_steps
from .updates import drain_steps


class _SourceVisibilityCuts(BoundaryOwnershipCuts):
    def __init__(self, proof):
        self.visibility_source_proof = proof

    def cuts(self, group, context):
        yield from self._physical_cuts(group)


def source_visibility(parent, group, context, previous):
    return drain_steps(source_visibility_steps(parent, group, context, previous))


def source_visibility_steps(parent, group, context, previous):
    if (context is None or group.physical_layout is not None or not group.claims
            or len(group.members) < 2 or len(group.loops) != 1 or group.loops[0].hole):
        raise ValueError('Source visibility recovery needs one original multi-junction compound exterior')
    if CompoundPortPlanner.incomplete_ports(group):
        raise ValueError('Source visibility recovery cannot replace an incomplete original external profile')
    sampled = AnalyticalBoundarySampler(tolerance=context.network.settings.curve_tolerance,
        spacing=parent.spacing, piecewise=parent.piecewise_sampling).build(group)
    sources = {claim.section_id for claim in group.claims}
    original_proof = yield from certify_steps(context, sampled, sources)
    chart = CompoundQuadChart(parent.patterns, spacing=parent.spacing, max_cells=parent.max_cells,
        quality_iterations=parent.quality_iterations, piecewise_sampling=parent.piecewise_sampling,
        partitioner=CompoundPartitioner(policy=_SourceVisibilityCuts(original_proof)))
    result = yield from chart._build_steps(group, context, require_quality=True)
    if (len(sampled.loops) != len(result.partition.group.loops) or
            any(not np.array_equal(first.xy, last.xy) or first.spans != last.spans
                for first, last in zip(sampled.loops, result.partition.group.loops))):
        raise ValueError('Source visibility recovery changed its sampled physical rim or ownership')
    result = yield from _finish_steps(chart, result, previous, policy='SOURCE_WIDTH_VISIBILITY')
    final_proof = yield from certify_steps(context, result.partition.group, sources)
    result.metadata = dict(result.metadata, source_width_visibility=final_proof,
        visibility_kernel_recovery=dict(previous_chart_failures=previous,
            original_width_proof=original_proof, final_width_proof=final_proof,
            physical_boundary_changed=False, original_profiles=True, cell_budget=parent.max_cells))
    return result


def recover_kernel_steps(parent, group, context, failure):
    previous = [dict(policy='ORIGINAL_PARTITION', reason=str(failure))]
    # These established constructions remain first, with their full assembled
    # quality checks. The source visibility policy is the last new alternative.
    for name, construct in (('LONGITUDINAL', parent._longitudinal_steps), ('REFLEX', parent._reflex_steps),
                            ('FIXED_RIM_COMPOUND', lambda g,c,p: recover_steps(parent,g,c,p)),
                            ('ALL_NOSES_LONGITUDINAL', lambda g,c,p: recover_longitudinal_steps(parent,g,c,p))):
        try:
            result = yield from construct(group, context, previous)
            result = yield from _finish_steps(parent, result, previous, policy=name)
        except ValueError as rejection:
            previous.append(dict(policy=name, reason=str(rejection)))
            continue
        result.metadata = dict(result.metadata, visibility_kernel_recovery=dict(
            previous_chart_failures=previous, policy=name, physical_boundary_changed=False,
            cell_budget=parent.max_cells))
        return result
    return (yield from source_visibility_steps(parent, group, context, previous))


def recover_splitter_steps(parent,group,context,failure):
    """Exhaust old charts before a source-certified outer-reflex combination."""
    previous=[dict(policy='ORIGINAL_VISIBILITY_CUTS',reason=str(failure))]
    try:
        return (yield from recover_kernel_steps(parent,group,context,failure))
    except ValueError as rejection:
        previous.append(dict(policy='EXISTING_KERNEL_RECOVERY',reason=str(rejection)))
    from .source_reflex_atlas import eligible,build
    if not eligible(group):raise failure
    try:
        return (yield from parent._filled_source_steps(group,context,group,str(failure)))
    except ValueError as rejection:
        previous.append(dict(policy='EXISTING_FILLED_SOURCE_NORMALS',reason=str(rejection)))
    from .source_exterior_sampling import eligible as sampling_eligible,recover_sampling
    if sampling_eligible(group,context):
        return recover_sampling(parent,group,context,failure)
    return build(parent,group,context,previous)
