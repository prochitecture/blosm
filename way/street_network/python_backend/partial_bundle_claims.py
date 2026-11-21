"""Original contiguous Bundle intervals wholly owned by a local junction pool.

This is a source-ownership permission, never a separator acceptance shortcut.
Visible Bundle runs still require their exact accepted junction curb caps,
centerline separation and final pavement/grass overlap qualification.
"""
from dataclasses import dataclass
from .model import End, Key


@dataclass(frozen=True)
class PartialBundleClaim:
    bundle: int
    role: int
    sections: tuple
    cuts: tuple
    bundle_end: str | None = None

    def report(self):
        return dict(bundle=self.bundle, member=self.role, sections=list(self.sections),
            policy=('RETAINED_'+self.bundle_end+'_HALF_PROFILE' if not self.sections else
                    'WHOLE_SOURCE_CYCLE_MEMBER' if self.bundle_end=='WHOLE' else
                    'PAIRED_'+self.bundle_end+'_MEMBER_INTERVAL' if self.bundle_end
                    else 'CONTIGUOUS_INTERNAL_MEMBER_INTERVAL'),
            complete_exterior_half_profiles=[dict(section=e.section_id,end=e.end) for e in self.cuts],
            separator_policy='accepted_junction_capped_visible_runs', endpoint_extrapolation=False)


def prepare(context, bundle, consumed, approaches, members):
    """Prove contiguous member intervals and surviving source incidences.

    Interior cuts can affect either member. A consumed Bundle end must include
    both corresponding lateral terminal intervals and their actual original
    terminal sites. No coordinate search creates a continuation, and every
    exterior cut retains its original half profile.
    """
    if len(bundle.members)!=2 or not all(member.complete for member in bundle.members):
        raise ValueError('Partial Bundle ownership needs two complete original members')
    active=[(role,[i for i,v in enumerate(member.path) if v.section_id in consumed])
            for role,member in enumerate(bundle.members)]
    active=[(role,indices) for role,indices in active if indices]
    incident={end:j.key for j in (*context.network.intersections,*context.network.connections) for end in j.ends}
    if set(bundle.sections)<=consumed:
        claims=[]
        for role,member in enumerate(bundle.members):
            first,last=member.path[0],member.path[-1]
            ends=(End(first.section_id,0 if first.forward else 1),
                  End(last.section_id,1 if last.forward else 0))
            if any(incident.get(end) not in members for end in ends):
                raise ValueError('Whole Bundle ownership requires all original terminal junctions')
            claims.append(PartialBundleClaim(bundle.id,role,tuple(v.section_id for v in member.path),(), 'WHOLE'))
        return tuple(claims)
    by_end={End(a.source.id,a.end):a for a in approaches}
    def validate_cut(end):
        approach=by_end.get(end)
        expected=context[Key('profile',end.section_id)].oriented(end.end)
        if (approach is None or incident.get(end) not in members or not expected.half
                or approach.profile!=expected or getattr(approach.path,'absorbed_sites',())
                or getattr(approach.path,'mouth_piece',None) is not None):
            raise ValueError('Partial Bundle interval lost an original exterior half-profile incidence')
    candidates=[]
    for role,indices in active:
        path=bundle.members[role].path
        first,last=indices[0],indices[-1]
        if indices!=list(range(first,last+1)) or len(indices)==len(path):
            raise ValueError('Partial Bundle ownership needs a contiguous interval with surviving sources')
        bundle_end='HEAD' if first==0 else 'TAIL' if last==len(path)-1 else None
        cuts=[]
        if first>0:
            before=path[first-1];cuts.append(End(before.section_id,1 if before.forward else 0))
        if last<len(path)-1:
            after=path[last+1];cuts.append(End(after.section_id,0 if after.forward else 1))
        for end in cuts:
            validate_cut(end)
        if bundle_end:
            visit=path[0 if bundle_end=='HEAD' else -1]
            terminal=End(visit.section_id,(0 if visit.forward else 1) if bundle_end=='HEAD'
                         else (1 if visit.forward else 0))
            if incident.get(terminal) not in members:
                raise ValueError('Partial Bundle end lacks its original terminal junction owner')
        candidates.append(PartialBundleClaim(bundle.id,role,tuple(path[i].section_id for i in indices),
            tuple(cuts),bundle_end))
    if len(candidates)==1 and candidates[0].bundle_end:
        # A long terminal Section on the opposite member need not be consumed
        # just because the short terminal Section of this member is internal.
        # Its original terminal must instead survive as a complete half port.
        claim=candidates[0];other=1-claim.role
        path=bundle.members[other].path
        visit=path[0 if claim.bundle_end=='HEAD' else -1]
        end=End(visit.section_id,(0 if visit.forward else 1) if claim.bundle_end=='HEAD'
                else (1 if visit.forward else 0))
        validate_cut(end)
        if incident.get(End(end.section_id,1-end.end)) in members:
            raise ValueError('Retained Bundle terminal actually belongs to the consumed source interval')
        candidates.append(PartialBundleClaim(bundle.id,other,(),(end,),claim.bundle_end))
    ends={claim.bundle_end for claim in candidates if claim.bundle_end}
    if ends and (len(candidates)!=2 or len(ends)!=1 or any(not c.bundle_end for c in candidates)):
        raise ValueError('Partial Bundle terminal ownership requires both corresponding lateral intervals')
    return tuple(candidates)


def certify_whole_separation(context, group, claims):
    """Entire consumed nonconverging pairs still need their source-curve proof."""
    from .bundles import CenterlineRange,certify_disjoint
    ids={claim.bundle for claim in claims if claim.bundle_end=='WHOLE'}
    existing=getattr(group.physical_layout,'source_cycle_proofs',{})
    for bundle in context.network.bundles:
        if bundle.id not in ids or bundle.id in existing:continue
        ranges=[tuple(CenterlineRange.build(context[Key('path',v.section_id)],0.,
            context[Key('path',v.section_id)].length) for v in member.path) for member in bundle.members]
        existing[bundle.id]=certify_disjoint(*ranges)
    if ids:group.physical_layout.source_cycle_proofs=existing


def certify(context, group, claims):
    """Require complete physical source coverage and actual surviving profiles."""
    if not claims:return group
    from .compound_coverage import CompoundSourceCoverage
    from .conditions import BoundaryConditions
    coverage=CompoundSourceCoverage(context,group,{})
    layouts={source.region.owner:source.layout for source in group.footprints}
    exposed={alias.endpoint for loop in group.loops
             for port in BoundaryConditions().build(group,loop,layouts).ports for alias in port.aliases}
    for claim in claims:
        if not set(claim.sections)<=coverage.fully_owned:
            raise ValueError('Partial Bundle interval lacks complete physical ownership')
        if any(end not in exposed or end.section_id in coverage.fully_owned for end in claim.cuts):
            raise ValueError('Partial Bundle interval hides an exterior half-profile cut')
    from .source_strip_containment import certify as certify_source_strips
    sections={sid for claim in claims for sid in claim.sections}
    whole={claim.bundle for claim in claims if claim.bundle_end=='WHOLE'}
    sections.update(sid for bundle in context.network.bundles if bundle.id in whole
                    for sid in bundle.interior_section_ids)
    width_proof=certify_source_strips(context,group,sections)
    group.physical_layout.partial_bundle_claims=tuple(claims)
    group.physical_layout.partial_bundle_width_proof=width_proof
    return group
