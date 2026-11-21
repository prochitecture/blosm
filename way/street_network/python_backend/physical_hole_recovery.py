"""A final source-normal atlas, qualified by complete finite source strips."""
from .four_side_charts import quality
from .mesh import publication_coordinates
from .partition import CompoundQuadChart,CompoundPartitioner,QuadChartGradingError,GrassChartGradingError
from .physical_hole_cuts import PhysicalHoleCuts
from .source_strip_containment import certify
from .updates import drain_steps


def build(parent,group,context,previous):
    return drain_steps(build_steps(parent,group,context,previous))


def build_steps(parent,group,context,previous):
    if not PhysicalHoleCuts.eligible(group) or not group.claims:
        raise ValueError('Source-normal recovery needs a prepared compound with actual consumed source strips')
    sources={claim.section_id for claim in group.claims}
    original=certify(context,group,sources)
    chart=CompoundQuadChart(parent.patterns,spacing=parent.spacing,max_cells=parent.max_cells,
        quality_iterations=min(100,parent.quality_iterations),four_sided=True,
        partitioner=CompoundPartitioner(policy=PhysicalHoleCuts(spacing=parent.spacing,
            source_corridors=True,existing_vertices=True),require_visibility=False))
    energy_failure=None
    try:
        result=yield from chart.build_steps(group,context)
    except QuadChartGradingError as failure:
        if isinstance(failure,GrassChartGradingError) or failure.candidate is None or failure.candidate.material=='grass':
            raise
        from .source_normal_energy import EnergySourceNormalChart
        energy_failure=str(failure)
        chart=EnergySourceNormalChart(parent.patterns,spacing=parent.spacing,max_cells=parent.max_cells,
            quality_iterations=min(100,parent.quality_iterations),four_sided=True,recover_grass=True,
            partitioner=CompoundPartitioner(policy=PhysicalHoleCuts(spacing=parent.spacing,
                source_corridors=True,existing_vertices=True),require_visibility=False))
        result=yield from chart.build_steps(group,context)
    yield
    proof=certify(context,result.partition.group,sources)
    if len(result.mesh.quads)+sum(len(hole.mesh.quads) for hole in result.holes.values())>parent.max_cells:
        raise ValueError('Source-normal recovery exceeds the compound cell budget')
    checks=[]
    for mesh in (result.mesh,*(hole.mesh for hole in result.holes.values())):
        mesh.validate()
        double=quality(mesh)
        published=mesh.with_positions(publication_coordinates(mesh.xyz,1.,physical=True))
        published.validate()
        stored=quality(published)
        if (max(double['maximum_edge_ratio'],stored['maximum_edge_ratio'])>8.
                or (mesh.material!='grass'
                    and min(double['minimum_corner_sine'],stored['minimum_corner_sine'])<.03)):
            raise ValueError('Source-normal recovery fails unchanged physical or stored quad quality')
        checks.append(dict(material=mesh.material,physical=double,stored=stored))
    result.metadata=dict(result.metadata,physical_hole_partition=True,source_corridors=True,
        existing_cut_vertices=True,source_width_visibility=proof,
        source_normal_recovery=dict(previous_chart_failures=previous,original_width_proof=original,
            final_width_proof=proof,stored_checks=checks,physical_boundary_changed=False,
            original_source_paths=True,cell_budget=parent.max_cells,
            **(dict(previous_source_normal_failure=energy_failure,energy_recoveries=chart.energy_recoveries)
               if energy_failure is not None else {})))
    return result
