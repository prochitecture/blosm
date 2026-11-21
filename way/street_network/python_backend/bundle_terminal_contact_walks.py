"""Finite inner-bank contacts across an owned full-profile terminal link.

A short full-profile road can join the two original ends of a paired Bundle.
The local circle may reach the opposite member beyond that consumed road.
Keep its original centreline as the construction closure; each contact keeps
its actual original half profile and width. No uniform-profile Route replaces
the consumed source or the opposite member.
"""
from dataclasses import replace
import math
import numpy as np
from .model import End,Key


def _original(context,approach):
    from .routes import RouteTopology
    path=approach.path;pieces=getattr(path,'pieces',())
    if (not pieces or len(pieces)>32 or approach.length<=0.
            or hasattr(approach,'original_approaches') or getattr(path,'absorbed_sites',())
            or getattr(path,'mouth_piece',None) is not None
            or getattr(path,'deferred_profile_reserve',None) is not None
            or abs(approach.offset)>1e-8):return False
    end=End(approach.source.id,approach.end)
    source=context.network.section_index[end.section_id]
    if (source is not pieces[0].section or abs(approach.width-source.width)>1e-8
            or end.end in (*source.scope_ends,*source.selection_ends)
            or approach.profile!=context[Key('profile',source.id)].oriented(end.end)):
        return False
    route=RouteTopology(context.network).trace(end)
    return tuple((p.section.id,p.forward) for p in pieces)==tuple(
        (v.section_id,v.forward) for v in route.traversals)


def _terminal(member,end):
    visit=member.path[0 if end==0 else -1]
    return End(visit.section_id,int(not visit.forward) if end==0 else int(visit.forward))


def contact_pair(context,first,second,members,consumed):
    """Return actual adjacent source walks, or decline this narrow policy."""
    from .profiles import memberships
    from .routes import RouteTopology
    from .compound_displaced_contacts import SourceContactWalk
    pair=(first,second);active=[i for i,a in enumerate(pair) if a.source.id in consumed]
    if len(active)!=1 or not all(_original(context,a) for a in pair):return None
    index=active[0];internal=pair[index];external=pair[1-index]
    if internal.profile.half or abs(external.profile.half)!=1:return None
    if (any(p.section.id not in consumed for p in internal.path.pieces)
            or any(p.section.id in consumed for p in external.path.pieces)):return None
    bundled=memberships(context.network);entry=bundled.get(external.source.id)
    if entry is None or internal.source.id in bundled:return None
    bundle,role,_=entry
    if len(bundle.members)!=2 or not all(member.complete for member in bundle.members):return None
    external_end=End(external.source.id,external.end)
    ends=[end for end in (0,1) if _terminal(bundle.members[role],end)==external_end]
    if len(ends)!=1:return None
    bundle_end=ends[0];target_end=_terminal(bundle.members[1-role],bundle_end)
    if target_end.section_id in consumed:return None
    incident={e:j for j in (*context.network.intersections,*context.network.connections) for e in j.ends}
    root=incident.get(external_end)
    internal_end=End(internal.source.id,internal.end)
    route=RouteTopology(context.network).trace(internal_end)
    far=next((j for j in context.network.intersections if j.key==route.terminal_site),None)
    if (root is None or internal_end not in root.ends or root.key not in members
            or far is None or far.key not in members or root.key==far.key or len(far.ends)<3
            or route.terminal_scope or route.terminal_selection or target_end not in far.ends):return None
    # Source incidence selects the far site and its facing bank's immediate
    # angular neighbour. Coordinates only certify the already selected seam.
    from .ownership import footprint_key
    original_layout=context.products.get(footprint_key(far.key))
    if original_layout is None:return None
    ordered=sorted(original_layout.approaches,key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
    source_ends=[End(a.source.id,a.end) for a in ordered]
    last=route.traversals[-1];reverse=End(last.section_id,int(last.forward))
    if set(source_ends)!=set(far.ends) or reverse not in source_ends:return None
    side=-1 if index==0 else 1
    following=ordered[(source_ends.index(reverse)+side)%len(ordered)]
    if End(following.source.id,following.end)!=target_end or not _original(context,following):return None
    if any(p.section.id in consumed for p in following.path.pieces):return None
    local_first=reverse if side>0 else target_end
    if far.corner_radius(local_first)[1]!='STYLE':return None
    contact=[first,second];contact[index]=following
    facing=[]
    for approach,curb_side in zip(contact,(-1.,1.)):
        if (abs(approach.profile.half)!=1 or approach.profile.peak!=curb_side
                or abs(float(approach.profile.derivative(curb_side)))>1e-10):return None
        facing.append(approach.profile.datum+float(approach.profile.value(curb_side)))
    if (abs(facing[0]-facing[1])>1e-10 or facing[0]-contact[0].profile.datum<=1e-10
            or abs(contact[0].profile.datum-contact[1].profile.datum)>1e-10
            or abs(internal.profile.datum-contact[0].profile.datum)>1e-10
            or abs(contact[0].profile.rise-contact[1].profile.rise)>1e-10
            or any(abs(a.width-internal.width)>1e-8 for a in contact)
            or np.linalg.norm(internal.frame(internal.length)[0]-following.origin)>1e-8
            or np.linalg.norm(first.origin-second.origin)>1e-8):return None
    walks=[SourceContactWalk(a,a) for a in pair]
    walks[index]=SourceContactWalk(internal,following,tuple(internal.path.pieces))
    report=dict(policy='ORIGINAL_OPPOSED_BUNDLE_TERMINAL_BANK_CONTACT',bundle=bundle.id,
        bundle_end=bundle_end,original_junction=root.key.text(),owned_terminal_junction=far.key.text(),
        consumed_original_prefix=[p.section.id for p in internal.path.pieces],
        original_half_profile_peak_height=facing[0],unchanged_original_profiles=True,
        declared_adjacent_source_incidence=True,endpoint_extrapolation=False)
    return tuple(walks),report


def retry(context,first,second,radius,mode,members,consumed):
    """Failed-only original STYLE circles, in their established radius order."""
    if mode!='STYLE':return None
    prepared=contact_pair(context,first,second,members,consumed)
    if prepared is None:return None
    walks,proof=prepared;a,b=(walk.approach for walk in walks)
    from .corner_radii import CornerRadiusPolicy,CornerContactError
    from .compound_contact_routes import regular_contact
    from .compound_displaced_contacts import patch
    from .intersections import _corner
    gap=math.atan2(first.direction[0]*second.direction[1]-first.direction[1]*second.direction[0],
        first.direction@second.direction)%math.tau
    side=1. if gap<math.pi else -1.;settings=context.network.settings
    tolerance=min(.001,.05*settings.curve_tolerance);attempted=[]
    for value in CornerRadiusPolicy.contact_candidates(first,second,radius,mode):
        attempted.append(float(value))
        try:
            corner=_corner(a,b,value,math.radians(settings.continuation_angle_degrees),side_override=side)
            if (corner.center is None or any(s<0. or s>p.length for p,s in zip((a,b),corner.stations))
                    or not all(regular_contact(p,s) for p,s in zip((a,b),corner.stations))):continue
            corner=replace(corner,requested_radius=radius,radius_mode=mode,
                radius_attempts=tuple(attempted),fit_policy=proof['policy'])
            corner.source_contact_walks=walks
            corner.source_contact_support=dict(**proof,radius=value,
                radius_attempts=list(attempted),walks=[walk.report(station)
                    for walk,station in zip(walks,corner.stations)])
            patch(Key('bundle_terminal_contact_proof',0),corner,settings.spacing,tolerance)
            return corner
        except (CornerContactError,ValueError):continue
    return None
