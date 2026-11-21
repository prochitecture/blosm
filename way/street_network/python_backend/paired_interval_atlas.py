"""Late source-normal charts for both retained original Bundle intervals.

Opposed member normals can intersect into poor computational blocks. Select one
member's axes while retaining both complete physical intervals and all original
source claims. This changes no pavement, grass hole or incoming profile.
"""
from dataclasses import dataclass
import math
from .model import Key
from .cycle_dimensional_retry import certify_layout
from .four_side_charts import quality
from .mesh import publication_coordinates
from .partition import CompoundPartitioner
from .physical_hole_cuts import PhysicalHoleCuts
from .source_normal_energy import EnergySourceNormalChart
from .source_identity import matches_current_section

PAIRED_PREFIX_CHART_POLICY='1-finite-original-paired-prefix-normal-charts'


def eligible(group):
    layout=group.physical_layout
    if not PhysicalHoleCuts.eligible(group) or not getattr(layout,'partial_bundle_width_proof',None):
        return False
    roles={}
    for claim in getattr(layout,'partial_bundle_claims',()):
        if claim.bundle_end in ('HEAD','TAIL') and claim.sections:
            roles.setdefault((claim.bundle,claim.bundle_end),set()).add(claim.role)
    return any(value=={0,1} for value in roles.values())


def original_intervals(context,group):
    from .partial_bundle_claims import prepare
    proof=getattr(group.physical_layout,'partial_bundle_width_proof',{})
    if (not isinstance(proof,dict) or not proof.get('sections')
            or not 0.<proof.get('maximum_error',float('inf'))<=.002):
        raise ValueError('Paired-interval atlas needs a complete original source-width proof')
    stored=tuple(getattr(group.physical_layout,'partial_bundle_claims',()))
    source_ids={claim.section_id for claim in group.claims}
    bundles={bundle.id:bundle for bundle in context.network.bundles}
    accepted=[]
    for bid in sorted({claim.bundle for claim in stored}):
        if bid not in bundles:
            raise ValueError('Paired-interval atlas lost its original Bundle')
        own=tuple(claim for claim in stored if claim.bundle==bid)
        if (len(own)!=2 or {claim.role for claim in own}!={0,1}
                or len({claim.bundle_end for claim in own})!=1
                or own[0].bundle_end not in ('HEAD','TAIL') or not all(claim.sections for claim in own)):
            continue
        fresh=prepare(context,bundles[bid],source_ids,group.physical_layout.approaches,set(group.members))
        if set(fresh)!=set(own):
            raise ValueError('Paired-interval atlas changed an original source interval or half-profile cut')
        accepted.extend(own)
    if not accepted:
        raise ValueError('Paired-interval atlas needs both nonempty original terminal intervals')
    return tuple(accepted)


class PairedIntervalCuts(PhysicalHoleCuts):
    def __init__(self,claims,role,**kwargs):
        if role not in (0,1):raise ValueError('Paired-interval cuts need an original member role')
        super().__init__(**kwargs)
        self.claims,self.role=tuple(claims),role

    def _source_inventory(self,group,context):
        selected,paired=super()._source_inventory(group,context)
        excluded={sid for claim in self.claims if claim.role!=self.role for sid in claim.sections}
        return selected-excluded,paired


@dataclass(frozen=True)
class OwnedPairedPrefixScope:
    claims: tuple
    ranges: tuple
    dependencies: tuple
    width_proof: object

    @classmethod
    def build(cls,context,group,claims,proof):
        from .partial_bundle_claims import PartialBundleClaim
        if (len(claims)!=2 or any(type(claim) is not PartialBundleClaim for claim in claims)
                or len({claim.bundle for claim in claims})!=1 or {claim.role for claim in claims}!={0,1}
                or len({claim.bundle_end for claim in claims})!=1
                or claims[0].bundle_end not in ('HEAD','TAIL') or not all(claim.sections for claim in claims)):
            raise ValueError('Owned-prefix normals require one typed original paired terminal interval')
        if not isinstance(proof,dict) or not proof.get('sections'):
            raise ValueError('Owned-prefix normals need a fresh complete original width proof')
        rows={};count=0
        for row in proof['sections']:
            sid=row['section'];first,last=map(float,row['interval']);count+=1
            path=context[Key('path',sid)]
            if (not matches_current_section(context,path.section)
                    or not all(map(math.isfinite,(first,last))) or not 0.<=first<last<=path.length
                    or count>20000):
                raise ValueError('Owned-prefix normal ranges changed an actual finite original source')
            rows.setdefault(sid,[]).append((first,last))
        ranges=[]
        for sid,values in sorted(rows.items()):
            values=tuple(sorted(values))
            if any(a[1]>b[0] for a,b in zip(values,values[1:])):
                raise ValueError('Owned-prefix normal ranges must remain disjoint, never a hull')
            ranges.append((sid,values))
        bundle=next(bundle for bundle in context.network.bundles if bundle.id==claims[0].bundle)
        dependencies={bundle.key,*(Key(kind,sid) for sid in rows for kind in ('path','profile'))}
        dependencies.update(Key(kind,visit.section_id) for member in bundle.members for visit in member.path
                            for kind in ('path','profile'))
        return cls(tuple(claims),tuple(ranges),tuple(sorted(dependencies)),proof)


def signature(context,group):
    """Project both role candidates and the complete original frame inputs."""
    layout=group.physical_layout
    if not getattr(layout,'partial_bundle_width_proof',None) or not eligible(group):return None
    claims=tuple(getattr(layout,'partial_bundle_claims',()))
    accepted=[claim for claim in claims if claim.bundle_end in ('HEAD','TAIL') and claim.sections]
    if len(accepted)!=2 or len({claim.bundle for claim in accepted})!=1:return None
    proof=getattr(layout,'source_cycle_full_width_proof',None)
    if not proof:return None
    # Construction reruns the proof and original interval validation. The
    # signature records its exact current input range ledger without doing
    # source geometry work on every cached surface lookup.
    scope=OwnedPairedPrefixScope.build(context,group,accepted,proof)
    return (PAIRED_PREFIX_CHART_POLICY,scope.ranges,
        tuple((key.text(),context.digest(key)) for key in scope.dependencies),
        tuple((claim.bundle,claim.role,claim.sections,tuple((end.section_id,end.end) for end in claim.cuts),
               claim.bundle_end) for claim in scope.claims))


class OwnedPairedPrefixCuts(PairedIntervalCuts):
    def __init__(self,scope,role,**kwargs):
        if type(role) is not int or role not in (0,1) or not isinstance(scope,OwnedPairedPrefixScope):
            raise ValueError('Owned-prefix cuts need an explicit original role and finite proof scope')
        super().__init__(scope.claims,role,source_ranges=dict(scope.ranges),cadence_multiplier=2,
            maximum_unpaired_span_widths=8.,**kwargs)
        self.scope=scope

    def _source_inventory(self,group,context):
        selected,all_intervals=super()._source_inventory(group,context)
        bundle_id=self.scope.claims[0].bundle
        bundle=next(bundle for bundle in context.network.bundles if bundle.id==bundle_id)
        member={visit.section_id for visit in bundle.members[self.role].path}
        owned=set(self.source_ranges)
        selected=(selected|member)&owned
        # This flag requests every actual pavement intersection of a normal.
        # It is a computational interval mode, not lateral source ownership.
        return selected,all_intervals|(selected&member)


def _checks(parent,result):
    meshes=(result.mesh,*(hole.mesh for hole in result.holes.values()))
    if sum(len(mesh.quads) for mesh in meshes)>parent.max_cells:
        raise ValueError('Paired-interval atlas exceeds the unchanged compound cell budget')
    checks=[]
    for mesh in meshes:
        mesh.validate();physical=quality(mesh)
        stored_mesh=mesh.with_positions(publication_coordinates(mesh.xyz,1.,physical=True))
        stored_mesh.validate();stored=quality(stored_mesh)
        if (max(physical['maximum_edge_ratio'],stored['maximum_edge_ratio'])>8.
                or (mesh.material!='grass' and min(physical['minimum_corner_sine'],stored['minimum_corner_sine'])<.03)):
            raise ValueError('Paired-interval atlas fails unchanged physical or stored quad quality')
        checks.append(dict(material=mesh.material,physical=physical,stored=stored))
    return checks


def _owned_prefix_chart_steps(parent,group,context,claims,proof,previous,failures):
    from .source_cycle_atlas import SourceCycleEnergyChart
    scope=OwnedPairedPrefixScope.build(context,group,claims,proof)
    class OwnedPrefixChart(SourceCycleEnergyChart):
        def _fibre_hole(self,boundary,max_cells):
            from .hole_fibre_chart import build_fibre_hole
            return build_fibre_hole(self.patterns,boundary,self.spacing,max_cells,ruled_paired=True)
    rows=[]
    for role in (0,1):
        policy=OwnedPairedPrefixCuts(scope,role,spacing=parent.spacing,source_corridors=True,
            hole_bridges=True,computational_extensions=True,existing_vertices=True)
        chart=OwnedPrefixChart(parent.patterns,spacing=parent.spacing,max_cells=parent.max_cells,
            quality_iterations=min(100,parent.quality_iterations),four_sided=True,recover_grass=True,
            partitioner=CompoundPartitioner(policy=policy,require_visibility=False))
        try:
            result=yield from chart.build_steps(group,context)
            yield
            final=certify_layout(context,result.partition.group,{claim.section_id for claim in group.claims})
            if result.partition.group.members!=group.members or result.partition.group.claims!=group.claims:
                raise ValueError('Computational owned-prefix normals changed original physical ownership')
            checks=_checks(parent,result)
        except ValueError as failure:
            rows.append(dict(role=role,reason=str(failure)));continue
        result.metadata=dict(result.metadata,physical_hole_partition=True,source_corridors=True,
            source_hole_bridges=True,source_width_visibility=final,
            owned_paired_prefix_atlas=dict(policy=PAIRED_PREFIX_CHART_POLICY,normal_member_role=role,
                original_intervals=[claim.report() for claim in claims],source_ranges=dict(scope.ranges),
                dependencies=[key.text() for key in scope.dependencies],all_physical_source_claims_preserved=True,
                physical_boundary_changed=False,original_width_proof=proof,final_width_proof=final,
                stored_checks=checks,old_role_failures=failures,previous_chart_failures=previous,
                previous_owned_role_failures=rows,maximum_new_role_attempts=2,cadence_multiplier=2,
                maximum_bands_per_original_section=64,maximum_computational_cuts=policy.maximum_cuts,
                unpaired_span_reports=policy.source_cut_reports,computational_extensions=policy.extension_reports,
                energy_recoveries=chart.energy_recoveries,cell_budget=parent.max_cells))
        return result
    raise ValueError('No qualified finite owned-prefix paired atlas: '+str(rows))


def build(parent,group,context,previous):
    from .updates import drain_steps
    return drain_steps(build_steps(parent,group,context,previous))


def build_steps(parent,group,context,previous):
    if not eligible(group):raise ValueError('Paired-interval atlas needs certified original paired terminal ownership')
    claims=original_intervals(context,group)
    sources={claim.section_id for claim in group.claims}
    before=certify_layout(context,group,sources)
    failures=[]
    for role in (0,1):
        try:
            policy=PairedIntervalCuts(claims,role,spacing=parent.spacing,
                source_corridors=True,hole_bridges=True)
            chart=EnergySourceNormalChart(parent.patterns,spacing=parent.spacing,
                max_cells=parent.max_cells,quality_iterations=min(100,parent.quality_iterations),
                four_sided=True,recover_grass=True,
                partitioner=CompoundPartitioner(policy=policy,require_visibility=False))
            result=yield from chart.build_steps(group,context)
            yield
            proof=certify_layout(context,result.partition.group,sources)
            checks=_checks(parent,result)
        except ValueError as failure:
            failures.append(dict(role=role,reason=str(failure)));continue
        result.metadata=dict(result.metadata,physical_hole_partition=True,source_corridors=True,
            source_hole_bridges=True,source_width_visibility=proof,
            paired_interval_atlas=dict(normal_member_role=role,original_intervals=[claim.report() for claim in claims],
                all_physical_source_claims_preserved=True,physical_boundary_changed=False,
                original_width_proof=before,final_width_proof=proof,stored_checks=checks,
                previous_chart_failures=previous,previous_member_failures=failures,
                energy_recoveries=chart.energy_recoveries,cell_budget=parent.max_cells))
        return result
    if len(claims)==2 and len({claim.bundle for claim in claims})==1:
        return (yield from _owned_prefix_chart_steps(parent,group,context,claims,before,previous,failures))
    raise ValueError('No qualified original paired-interval atlas: '+str(failures))
