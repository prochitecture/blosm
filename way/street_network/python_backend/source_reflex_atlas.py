"""Late exact-rim source normals combined with real outer curb bisectors."""
import numpy as np
from .physical_hole_cuts import PhysicalHoleCuts
from .reflex_cuts import ReflexCurbCuts
from .partition import CompoundPartitioner
from .source_normal_energy import EnergySourceNormalChart
from .cycle_dimensional_retry import certify_layout
from .four_side_charts import quality
from .mesh import publication_coordinates
from .sampling import AnalyticalBoundarySampler


class SourceReflexCuts(PhysicalHoleCuts):
    def _source_cuts(self,group,context,labels):
        original=super()._source_cuts(group,context,labels)
        outer=ReflexCurbCuts().cuts(group,context)
        cuts=[*original,*outer]
        if len(cuts)>self.maximum_cuts:
            raise ValueError('Combined source/reflex cuts exceed the original cut budget')
        self.source_normal_count,self.outer_reflex_count=len(original),len(outer)
        return cuts


def eligible(group):
    return (PhysicalHoleCuts.filled_eligible(group) and bool(group.claims)
        and bool(getattr(group.physical_layout,'source_cycle_full_width_proof',None)))


def build(parent,group,context,previous):
    if not eligible(group):
        raise ValueError('Source/reflex atlas needs one complete source-certified physical exterior')
    consumed={claim.section_id for claim in group.claims}
    before=certify_layout(context,group,consumed)
    sampled=AnalyticalBoundarySampler(tolerance=context.network.settings.curve_tolerance,
        spacing=parent.spacing,piecewise=parent.piecewise_sampling).build(group)
    policy=SourceReflexCuts(spacing=parent.spacing,source_corridors=True,allow_filled=True)
    chart=EnergySourceNormalChart(parent.patterns,spacing=parent.spacing,max_cells=parent.max_cells,
        quality_iterations=min(100,parent.quality_iterations),four_sided=True,
        partitioner=CompoundPartitioner(policy=policy,require_visibility=False))
    result=chart.build(group,context)
    current=result.partition.group
    if (current.members!=group.members or current.claims!=group.claims
            or len(sampled.loops)!=len(current.loops)
            or any(not np.array_equal(first.xy,last.xy) or first.spans!=last.spans
                   for first,last in zip(sampled.loops,current.loops))):
        raise ValueError('Source/reflex atlas changed its original physical rim or ownership')
    after=certify_layout(context,current,consumed)
    if result.holes or len(result.mesh.quads)>parent.max_cells:
        raise ValueError('Source/reflex atlas changed its hole topology or cell budget')
    result.mesh.validate()
    physical=quality(result.mesh)
    stored=quality(result.mesh.with_positions(publication_coordinates(result.mesh.xyz,1.,physical=True)))
    if (max(physical['maximum_edge_ratio'],stored['maximum_edge_ratio'])>8.
            or min(physical['minimum_corner_sine'],stored['minimum_corner_sine'])<.03):
        raise ValueError('Source/reflex atlas fails unchanged physical or stored quality')
    result.metadata=dict(result.metadata,filled_source_partition=True,source_width_visibility=after,
        source_normal_outer_reflex_recovery=dict(previous_chart_failures=previous,
            original_width_proof=before,final_width_proof=after,physical_boundary_changed=False,
            original_source_normal_cuts=policy.source_normal_count,
            actual_outer_reflex_cuts=policy.outer_reflex_count,maximum_cuts=policy.maximum_cuts,
            stored_quality=stored,physical_quality=physical,cell_budget=parent.max_cells,
            energy_recoveries=chart.energy_recoveries))
    return result
