"""Late bounded charts aligned to actual unbundled original source corridors.

Choosing computational normals never groups roads laterally or changes a
physical rim. Complete original widths, profiles, holes and stored quad
quality independently qualify each candidate before height construction.
"""
import numpy as np
from .model import Key
from .routes import Route
from .physical_hole_cuts import PhysicalHoleCuts
from .partition import CompoundPartitioner,CompoundQuadChart,QuadChartGradingError
from .source_normal_energy import EnergySourceNormalChart,_energy_attempt_steps
from .four_side_charts import quality
from .mesh import publication_coordinates
from .cycle_dimensional_retry import certify_layout


def eligible(group,context):
    layout=group.physical_layout
    if (not PhysicalHoleCuts.eligible(group) or not group.claims or len(group.members)>16
            or not hasattr(layout,'source_cycle_unexposed_corners')
            or not getattr(layout,'source_cycle_full_width_proof',None)):
        return False
    ids={claim.section_id for claim in group.claims}
    if any(ids.intersection(bundle.sections) for bundle in context.network.bundles):return False
    for a in layout.approaches:
        if (not isinstance(a.path,Route) or a.profile.half
                or getattr(a.path,'absorbed_sites',()) or hasattr(a.path,'mouth_piece')
                or getattr(a.path,'deferred_profile_reserve',None) is not None
                or abs(a.width-a.source.width)>1e-8
                or abs(a.profile.rise-a.source.rise)>1e-10
                or abs(a.profile.datum-a.source.datum)>1e-10):return False
    return True


def _score(mesh):
    points=mesh.xyz[mesh.quads,:2]
    lengths=np.linalg.norm(np.roll(points,-1,axis=1)-points,axis=2)
    ratios=lengths.max(1)/lengths.min(1)
    return int((ratios>8.).sum()),float(ratios.max())


class OriginalCorridorCuts(PhysicalHoleCuts):
    def __init__(self,section_id,**kwargs):
        super().__init__(**kwargs);self.section_id=section_id

    def _source_inventory(self,group,context):
        if self.section_id not in {claim.section_id for claim in group.claims}:
            raise ValueError('Computational corridor must remain a consumed original source')
        return {self.section_id},set()


class SourceCycleEnergyChart(EnergySourceNormalChart):
    def _correct_steps(self,mesh,*,retry_rotations=False):
        try:
            return (yield from CompoundQuadChart._correct_steps(self,mesh,retry_rotations=retry_rotations))
        except QuadChartGradingError as original:
            if mesh.material=='grass' or original.candidate is None:raise
            best=original.candidate;failures=[]
            # Retain the existing two attempts, including their original seed.
            for power in (1,2):
                try:
                    result,correction,energy,physical,stored=yield from _energy_attempt_steps(
                        self,mesh,original.candidate,power,retry_rotations)
                except ValueError as failure:
                    failures.append(dict(power=power,reason=str(failure)))
                    candidate=getattr(failure,'candidate',None)
                    if candidate is not None and _score(candidate)<_score(best):best=candidate
                    continue
                return self._accepted(result,correction,energy,physical,stored,power,original,failures)
            if not 1<=_score(best)[0]<=128:raise original
            continuation=[]
            for power in (1,2,1,2):
                before=_score(best)
                try:
                    result,correction,energy,physical,stored=yield from _energy_attempt_steps(self,mesh,best,power,retry_rotations)
                except ValueError as failure:
                    row=dict(power=power,before=before,reason=str(failure))
                    candidate=getattr(failure,'candidate',None)
                    if candidate is not None:
                        row['after']=_score(candidate)
                        if _score(candidate)<_score(best):best=candidate
                    continuation.append(row);continue
                return self._accepted(result,correction,energy,physical,stored,power,original,failures,
                    dict(maximum_attempts=4,maximum_initial_bad_cells=128,
                         actual_attempts=len(continuation)+1,previous_attempts=continuation))
            original.details=dict(original.details or {},bounded_energy_continuation=continuation)
            raise original

    def _accepted(self,result,correction,energy,physical,stored,power,original,failures,continuation=None):
        details=dict(power=power,energy_iterations=500,correction_energy_iterations=300,
            energy=energy,previous_correction=original.details,previous_energy_failures=failures,
            physical_quality=physical,stored_quality=stored,fixed_source_rim=True,cell_budget=self.max_cells)
        if continuation is not None:details['bounded_continuation']=continuation
        self.energy_recoveries.append(details)
        return result,dict(correction,source_normal_energy_recovery=details)

    def recover_separator_attachment(self,mesh,refs,guides,previous,refinements):
        # Only obsolete interior guides are released. All actual physical rim
        # coordinates and IDs are mandatory fixed data of the general solver.
        parent=EnergySourceNormalChart(self.patterns,spacing=self.spacing,
            max_cells=self.max_cells,quality_iterations=300)
        fixed={key:(ids.copy(),mesh.xyz[ids].copy()) for key,ids in mesh.seams.items()}
        result,correction=parent._correct(mesh,retry_rotations=True)
        for key,(ids,xyz) in fixed.items():
            if not np.array_equal(result.seams[key],ids) or not np.array_equal(result.xyz[ids],xyz):
                raise ValueError('Source attachment correction changed its physical rim')
        physical=quality(result)
        stored=quality(result.with_positions(publication_coordinates(result.xyz,1.,physical=True)))
        if (max(physical['maximum_edge_ratio'],stored['maximum_edge_ratio'])>8.
                or min(physical['minimum_corner_sine'],stored['minimum_corner_sine'])<.03):
            raise ValueError('Source attachment correction fails unchanged stored quality')
        return result,refs,guides,dict(correction,qualified=True,
            source_attachment_global_correction=dict(previous_local=previous,previous_refinements=refinements,
                physical_rim_exact=True,interior_guides_released=True,
                physical=physical,stored=stored,energy=parent.energy_recoveries))


def build(parent,group,context,previous):
    if not eligible(group,context):
        raise ValueError('Original-corridor atlas needs a complete unbundled source-cycle width proof')
    sources={claim.section_id for claim in group.claims}
    before=certify_layout(context,group,sources)
    selected=sorted(sources,key=lambda sid:(-context[Key('path',sid)].length,sid))[:3]
    failures=[]
    # Keep every unchanged corridor attempt first. Only their complete failure
    # admits one bounded cut-phase alternative on those same three sources.
    for near_fraction in (0., .125):
        for sid in selected:
            policy=OriginalCorridorCuts(sid,spacing=parent.spacing,source_corridors=True,
                hole_bridges=True,computational_extensions=True,near_endpoint_fraction=near_fraction)
            chart=SourceCycleEnergyChart(parent.patterns,spacing=parent.spacing,max_cells=parent.max_cells,
                quality_iterations=min(100,parent.quality_iterations),four_sided=True,recover_grass=True,
                partitioner=CompoundPartitioner(policy=policy,require_visibility=False))
            try:
                result=chart.build(group,context)
                proof=certify_layout(context,result.partition.group,sources)
                meshes=(result.mesh,*(hole.mesh for hole in result.holes.values()))
                if sum(len(mesh.quads) for mesh in meshes)>parent.max_cells:
                    raise ValueError('Original-corridor atlas exceeds the unchanged compound cell budget')
                checks=[]
                for mesh in meshes:
                    mesh.validate();physical=quality(mesh)
                    published=mesh.with_positions(publication_coordinates(mesh.xyz,1.,physical=True))
                    published.validate();stored=quality(published)
                    if (max(physical['maximum_edge_ratio'],stored['maximum_edge_ratio'])>8.
                            or min(physical['minimum_corner_sine'],stored['minimum_corner_sine'])<.03):
                        raise ValueError('Original-corridor atlas fails unchanged physical or stored quad quality')
                    checks.append(dict(material=mesh.material,physical=physical,stored=stored))
            except ValueError as failure:
                row=dict(original_section=sid,reason=str(failure))
                if near_fraction:row['near_endpoint_fraction']=near_fraction
                failures.append(row);continue
            result.metadata=dict(result.metadata,physical_hole_partition=True,source_corridors=True,
                source_hole_bridges=True,source_width_visibility=proof,
                original_source_cycle_atlas=dict(original_corridor=sid,maximum_corridor_attempts=3,
                    source_length=context[Key('path',sid)].length,all_physical_source_claims_preserved=True,
                    lateral_grouping_inferred=False,physical_boundary_changed=False,
                    original_width_proof=before,final_width_proof=proof,stored_checks=checks,
                    previous_chart_failures=previous,previous_corridor_failures=failures,
                    computational_extensions=policy.extension_reports,
                    energy_recoveries=chart.energy_recoveries,cell_budget=parent.max_cells))
            if near_fraction:
                details=result.metadata['original_source_cycle_atlas']
                details.update(maximum_corridor_attempts=6,maximum_original_corridors=3,
                    maximum_endpoint_policies=2,actual_atlas_attempts=len(failures)+1,
                    near_endpoint_fraction=near_fraction,
                    endpoint_reuse_reports=policy.endpoint_reuse_reports,
                    physical_source_vertices_changed=False)
            return result
    raise ValueError('No qualified original source-cycle atlas: '+str(failures))
