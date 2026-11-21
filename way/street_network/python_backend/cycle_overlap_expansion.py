"""Expand only real connected ground incidences after a finite gore fails."""
from .model import End
from .corner_radii import CornerRadiusPolicy
from .overlapping_approaches import overlap_certificate


def required(planner,site,first,second,radius,mode,members,internal,consumed):
    """Return a real external terminal and its bounded overlap proof, or None."""
    if mode!='STYLE':return None
    ends=(End(first.source.id,first.end),End(second.source.id,second.end))
    flags=tuple(end in internal for end in ends)
    if flags[0]==flags[1] or not set(ends)<=set(planner.sites[site].ends):return None
    inside,outside=(first,second) if flags[0] else (second,first)
    if any(a.profile.half or getattr(a.path,'absorbed_sites',()) or getattr(a.path,'mouth_piece',None) is not None
            or getattr(a.path,'minimum_mouth',0.)>1e-9 or getattr(a.path,'deferred_profile_reserve',None) is not None
            for a in (inside,outside)):return None
    sites=set(members)
    routes=[]
    for approach in (inside,outside):
        route=planner.topology.trace(End(approach.source.id,approach.end))
        if (route.terminal_scope or route.terminal_selection
                or tuple((p.section.id,p.forward) for p in approach.path.pieces)!=tuple((v.section_id,v.forward) for v in route.traversals)):
            return None
        routes.append(route)
    original,external=routes
    if (original.terminal_site not in sites or original.terminal_site==site
            or any(p.section.id not in consumed for p in inside.path.pieces)
            or external.terminal_site not in planner.sites or external.terminal_site in sites):return None
    opposite=original.traversals[-1]
    other=End(opposite.section_id,1 if opposite.forward else 0)
    if other not in internal or other not in planner.sites[original.terminal_site].ends:return None
    minimum=min(CornerRadiusPolicy.contact_candidates(first,second,radius,mode))
    proof=overlap_certificate(inside,outside,clearance=2.*minimum,refine=True)
    if proof is None:return None
    record=dict(policy='FINITE_INTERNAL_SOURCE_OVERLAP_EXPANSION',junction=site.text(),
        terminal_junction=external.terminal_site.text(),internal_section=inside.source.id,
        external_section=outside.source.id,minimum_style_radius=minimum,
        original_internal_sections=[p.section.id for p in inside.path.pieces],
        original_external_sections=[p.section.id for p in outside.path.pieces],**proof)
    return external.terminal_site,record
