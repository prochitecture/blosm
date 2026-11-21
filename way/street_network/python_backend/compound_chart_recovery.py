"""A bounded fixed-rim alternative after a radial compound chart fails.

Only computational cuts and interior vertices differ. The original compound
owns its source intervals before this module runs; this creates no new region.
"""
from .partition import CompoundQuadChart,CompoundPartitioner
from .reflex_cuts import ReflexCurbCuts
from .mesh import publication_coordinates
from .four_side_charts import quality
from .updates import drain_steps


class _FixedRimCompoundChart(CompoundQuadChart):
    def _repairable_leaf(self,group):
        # Eligibility is checked on the complete compound before sampling.
        # A poor positive computational leaf is permitted only as input to
        # full-region correction, never as a qualified output by itself.
        return True


def _finish_steps(parent,result,previous,**policy):
    double=quality(result.mesh)
    correction=None
    if double['maximum_edge_ratio']>8. or double['minimum_corner_sine']<.03:
        result.mesh,correction=yield from parent._correct_steps(result.mesh)
        double=quality(result.mesh)
    stored=quality(result.mesh.with_positions(publication_coordinates(result.mesh.xyz,1.,physical=True)))
    if (max(double['maximum_edge_ratio'],stored['maximum_edge_ratio'])>8.
            or min(double['minimum_corner_sine'],stored['minimum_corner_sine'])<.03):
        raise ValueError('Recovered compound chart fails original physical or stored quality')
    result.metadata=dict(result.metadata,quality=double,
        fixed_rim_compound_recovery=dict(previous_chart_failures=previous,
            physical_boundary_changed=False,stored_quality=stored,correction=correction,
            cell_budget=parent.max_cells,all_direct_quads=True,**policy))
    return result


def recover(parent,group,context,previous):
    return drain_steps(recover_steps(parent,group,context,previous))


def recover_steps(parent,group,context,previous):
    if (group.physical_layout is None or not group.claims or len(group.members)<2
            or len(group.loops)!=1 or group.loops[0].hole):
        raise ValueError('Fixed-rim recovery needs one prepared multi-junction compound exterior')
    chart=_FixedRimCompoundChart(parent.patterns,spacing=parent.spacing,
        partitioner=CompoundPartitioner(policy=ReflexCurbCuts(),require_visibility=False),
        max_cells=parent.max_cells,quality_iterations=min(100,parent.quality_iterations),
        four_sided=True,piecewise_sampling=parent.piecewise_sampling)
    result=yield from chart._build_steps(group,context,require_quality=True)
    result.metadata['reflex_curb_partition']=True
    return (yield from _finish_steps(chart,result,previous,source_claims=[claim.section_id for claim in group.claims],
        policy='REFLEX_COMPOUND'))


def recover_longitudinal(parent,group,context,previous):
    return drain_steps(recover_longitudinal_steps(parent,group,context,previous))


def recover_longitudinal_steps(parent,group,context,previous):
    """Keep every substantial fork nose when an ordinary atlas has failed."""
    from .longitudinal_cuts import LongitudinalCuts
    if len(group.loops)!=1 or group.loops[0].hole:
        raise ValueError('Multi-nose longitudinal recovery needs one physical exterior')
    chart=_FixedRimCompoundChart(parent.patterns,spacing=parent.spacing,
        partitioner=CompoundPartitioner(policy=LongitudinalCuts(parent.spacing,
            minimum_aspect=2.,all_noses=True),require_visibility=False),
        max_cells=parent.max_cells,quality_iterations=min(100,parent.quality_iterations),
        four_sided=True,piecewise_sampling=parent.piecewise_sampling)
    result=yield from chart._build_steps(group,context,require_quality=True)
    return (yield from _finish_steps(chart,result,previous,policy='ALL_NOSES_LONGITUDINAL',minimum_aspect=2.,
        original_profile_rows=True))


def recover_source_reflex(parent,group,context,previous):
    """Try the maintained compound atlas after other grading families fail."""
    if (group.physical_layout is None or not group.claims or len(group.members)<2
            or len(group.loops)!=1 or group.loops[0].hole):
        raise ValueError('Source-reflex recovery needs a prepared compound exterior')
    from .source_strip_containment import certify
    sources={claim.section_id for claim in group.claims}
    original=certify(context,group,sources)
    result=recover(parent,group,context,previous)
    proof=certify(context,result.partition.group,sources)
    result.metadata['source_width_visibility']=proof
    result.metadata['fixed_rim_compound_recovery']['original_width_proof']=original
    result.metadata['fixed_rim_compound_recovery']['final_width_proof']=proof
    return result
