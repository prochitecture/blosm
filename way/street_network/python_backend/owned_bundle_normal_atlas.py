"""Late computational normal charts on proven original Bundle-owned ranges.

Unlike terminal PartialBundleClaims, this scope can include original incoming
prefixes with no completely consumed member Section. It grants no physical
ownership: the fresh continuous source-width proof and unchanged native ranges
bound all computational normals, while both original member roles remain real.
"""
from dataclasses import dataclass
import math
from .model import Key
from .physical_hole_cuts import PhysicalHoleCuts
from .source_bank_reconciliation import owned_ranges
from .source_identity import matches_current_section

SOURCE_OWNED_BUNDLE_NORMAL_POLICY='1-certified-owned-native-bundle-role-normals'


def _bundle(context,group,ranges):
    claimed={claim.section_id for claim in group.claims}
    candidates=[bundle for bundle in context.network.bundles if claimed.intersection(bundle.sections)]
    if len(candidates)!=1:return None
    bundle,=candidates
    if (len(bundle.members)!=2 or not all(member.complete for member in bundle.members)
            or not all(any(visit.section_id in ranges for visit in member.path) for member in bundle.members)):
        return None
    return bundle


def eligible(group,context):
    if not getattr(group,'claims',()) or not hasattr(group,'loops'):return False
    proof=getattr(getattr(group,'physical_layout',None),'source_cycle_full_width_proof',None)
    if (not 1<=len(getattr(group,'members',()))<=16 or not PhysicalHoleCuts.eligible(group)
            or not isinstance(proof,dict) or not proof.get('sections')
            or not 0.<proof.get('maximum_error',float('inf'))<=.002):return False
    if any(not isinstance(row,dict) or type(row.get('section')) is not int for row in proof['sections']):return False
    ranges={row['section'] for row in proof['sections']}
    return _bundle(context,group,ranges) is not None


@dataclass(frozen=True)
class OwnedBundleNormalScope:
    bundle: object
    ranges: tuple
    dependencies: tuple
    width_proof: object

    @classmethod
    def build(cls,context,group,proof):
        if (not PhysicalHoleCuts.eligible(group) or not group.claims or not 1<=len(group.members)<=16
                or not isinstance(proof,dict) or not proof.get('sections')
                or not 0.<proof.get('maximum_error',float('inf'))<=.002):
            raise ValueError('Owned Bundle normals need a continuous original source-width proof')
        rows={};count=0
        for row in proof['sections']:
            try:
                sid=row['section'];first,last=map(float,row['interval'])
            except (KeyError,TypeError,ValueError) as failure:
                raise ValueError('Owned Bundle normals need explicit original source interval rows') from failure
            if type(sid) is not int or sid not in context.network.section_index:
                raise ValueError('Owned Bundle normals lost an actual original source Section')
            path=context[Key('path',sid)]
            count+=1
            if (not matches_current_section(context,path.section) or count>20000
                    or not all(map(math.isfinite,(first,last))) or not 0.<=first<last<=path.length):
                raise ValueError('Owned Bundle normals changed an actual finite original source range')
            rows.setdefault(sid,[]).append((first,last))
        ranges={}
        for sid,values in sorted(rows.items()):
            values=tuple(sorted(values))
            if any(a[1]>b[0] for a,b in zip(values,values[1:])):
                raise ValueError('Owned Bundle normal ranges must stay disjoint, never a hull')
            ranges[sid]=values
        if ranges!=owned_ranges(context,group):
            raise ValueError('Owned Bundle normal ranges do not equal the actual prepared native ownership')
        bundle=_bundle(context,group,ranges)
        if bundle is None:raise ValueError('Owned Bundle normals need exactly one actual complete two-role Bundle')
        dependencies={bundle.key,*(Key(kind,sid) for sid in ranges for kind in ('path','profile'))}
        dependencies.update(Key(kind,visit.section_id) for member in bundle.members for visit in member.path
                            for kind in ('path','profile'))
        return cls(bundle,tuple(ranges.items()),tuple(sorted(dependencies)),proof)


class OwnedBundleNormalCuts(PhysicalHoleCuts):
    def __init__(self,scope,role,**kwargs):
        if type(scope) is not OwnedBundleNormalScope or type(role) is not int or role not in (0,1):
            raise ValueError('Owned Bundle normal cuts need a typed original scope and member role')
        super().__init__(source_ranges=dict(scope.ranges),cadence_multiplier=2,
            maximum_unpaired_span_widths=8.,**kwargs)
        self.scope,self.role=scope,role

    def _source_inventory(self,group,context):
        selected,paired=super()._source_inventory(group,context)
        bundle=self.scope.bundle
        if not any(original is bundle for original in context.network.bundles):
            raise ValueError('Owned Bundle normal cuts changed their original Bundle identity')
        member={visit.section_id for visit in bundle.members[self.role].path}
        other={visit.section_id for visit in bundle.members[1-self.role].path}
        selected=((selected-other)|member)&set(self.source_ranges)
        return selected,paired|(selected&member)


def signature(context,group):
    if not eligible(group,context):return None
    proof=group.physical_layout.source_cycle_full_width_proof
    try:scope=OwnedBundleNormalScope.build(context,group,proof)
    except ValueError:return None  # Optional late recovery cannot reject a valid earlier family.
    return (SOURCE_OWNED_BUNDLE_NORMAL_POLICY,scope.bundle.id,scope.ranges,
        tuple((key.text(),context.digest(key)) for key in scope.dependencies),
        tuple(tuple((visit.section_id,visit.forward) for visit in member.path) for member in scope.bundle.members))


def build(parent,group,context,previous):
    from .cycle_dimensional_retry import certify_layout
    from .source_cycle_atlas import SourceCycleEnergyChart
    from .partition import CompoundPartitioner
    from .paired_interval_atlas import _checks
    from .source_parent_chain_fillets import validate as validate_parent_chains
    if not eligible(group,context):raise ValueError('Owned Bundle normal atlas needs its original finite source scope')
    claimed={claim.section_id for claim in group.claims};proof=certify_layout(context,group,claimed)
    scope=OwnedBundleNormalScope.build(context,group,proof);failures=[]
    class ScopedNormalChart(SourceCycleEnergyChart):
        def _fibre_hole(self,boundary,max_cells):
            from .hole_fibre_chart import build_fibre_hole
            return build_fibre_hole(self.patterns,boundary,self.spacing,max_cells,ruled_paired=True)
    for role in (0,1):
        policy=OwnedBundleNormalCuts(scope,role,spacing=parent.spacing,source_corridors=True,
            hole_bridges=True,computational_extensions=True,existing_vertices=True)
        chart=ScopedNormalChart(parent.patterns,spacing=parent.spacing,max_cells=parent.max_cells,
            quality_iterations=min(100,parent.quality_iterations),four_sided=True,recover_grass=True,
            partitioner=CompoundPartitioner(policy=policy,require_visibility=False))
        try:
            result=chart.build(group,context);checks=_checks(parent,result)
            current=validate_parent_chains(result.partition.group)
            if current.members!=group.members or current.claims!=group.claims:
                raise ValueError('Computational Bundle normals changed original physical ownership')
            final=certify_layout(context,current,claimed)
            if owned_ranges(context,current)!=dict(scope.ranges):
                raise ValueError('Computational Bundle normals changed an original native profile plane')
            height_admission=(parent.source_height_admission(current,result)
                              if parent.source_height_admission is not None else None)
        except ValueError as failure:
            failures.append(dict(role=role,reason=str(failure)));continue
        result.metadata=dict(result.metadata,physical_hole_partition=True,source_corridors=True,
            source_hole_bridges=True,source_width_visibility=final,
            owned_bundle_normal_atlas=dict(policy=SOURCE_OWNED_BUNDLE_NORMAL_POLICY,bundle=scope.bundle.id,
                normal_member_role=role,source_ranges=dict(scope.ranges),
                dependencies=[key.text() for key in scope.dependencies],
                physical_boundary_changed=False,all_original_claims_and_profiles_preserved=True,
                original_width_proof=proof,final_width_proof=final,stored_checks=checks,
                source_height_admission=height_admission,maximum_role_attempts=2,
                cadence_multiplier=2,maximum_bands_per_original_section=64,
                maximum_computational_cuts=policy.maximum_cuts,
                maximum_unpaired_span_widths=8.,source_cut_reports=policy.source_cut_reports,
                computational_extensions=policy.extension_reports,energy_recoveries=chart.energy_recoveries,
                previous_chart_failures=previous,previous_role_failures=failures,cell_budget=parent.max_cells))
        return result
    raise ValueError('No qualified finite owned Bundle normal atlas: '+str(failures))
