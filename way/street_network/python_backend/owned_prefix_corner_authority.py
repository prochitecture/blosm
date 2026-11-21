"""Logical corner policy with a typed, wholly owned prefix contact frame.

An OwnedPrefixRoute keeps its real junction End while its physical carriageway
frame belongs to a later original Section. Neither that later End nor its
profile is an incidence of the logical junction. This validator preserves both
identities and rechecks complete native support before exposing source ranges.
"""
from dataclasses import dataclass
import math
import numpy as np
from .model import End,Key
from .paths import Path
from .routes import Route,RouteTopology
from .intersections import Approach
from .owned_prefixes import OwnedPrefix,OwnedPrefixRoute,candidates
from .source_identity import matches_current_section

SOURCE_OWNED_PREFIX_CORNER_POLICY='1-real-prefix-logical-sector-physical-contact-authority'


class OwnedPrefixCornerAuthorityError(ValueError):
    """A recorded prefix lost its actual logical or finite physical proof."""


@dataclass(frozen=True)
class OwnedPrefixRangeProof:
    claim: OwnedPrefix
    ranges: tuple
    requested_station: float
    mouth_station: float
    final_mouth_ranges: tuple
    final_mouth_width_proof: object


@dataclass(frozen=True)
class OwnedPrefixCornerAuthority:
    junction: object
    logical_ends: tuple
    contact_ends: tuple
    source_walks: tuple
    prefix_range_proofs: tuple
    radius_attempts: tuple
    native_order: tuple
    proxy_order: tuple
    contact_width_proof: object


def _require(value,message):
    if not value:raise OwnedPrefixCornerAuthorityError(message)


def _route_pieces(context,route,traversals):
    _require(tuple((p.section.id,p.forward) for p in route.pieces)==tuple((v.section_id,v.forward) for v in traversals),
        'Owned prefix changed its original directed external source walk')
    _require(route.pieces and all(p.path is context[Key('path',p.section.id)]
        and matches_current_section(context,p.section) for p in route.pieces),
        'Owned prefix borrowed a different native Path or changed a source Section')
    stations=np.r_[0.,np.cumsum([p.path.length for p in route.pieces])]
    _require(np.array_equal(route.stations,stations) and route.length==float(stations[-1]),
        'Owned prefix changed its finite external Route metric')


def _prefix_arm(context,arm):
    _require(type(arm) is Approach and type(arm.path) is OwnedPrefixRoute,
        'Owned-prefix authority needs the explicit original route type')
    route=arm.path;claim=route.claim
    _require(type(claim) is OwnedPrefix,'Owned prefix lost its typed producer claim')
    prefix=context[Key('path',claim.endpoint.section_id)]
    _require(type(prefix) is Path and route.prefix is prefix
        and getattr(prefix,'owned_prefix',None) is claim and route.section is prefix.section
        and arm.source is prefix.section and arm.end==claim.endpoint.end
        and matches_current_section(context,prefix.section),
        'Owned prefix changed its canonical native Path/claim identity')
    _require(candidates(context.network).get(prefix.section.id)==claim,
        'Owned prefix no longer satisfies its original local producer contract')
    node=next((node for node in context.network.intersections if node.key==claim.junction),None)
    seams=[site for site in (*context.network.connections,*context.network.intersections)
           if site.key==claim.seam and len(site.ends)==2]
    _require(node is not None and claim.endpoint in node.ends and len(seams)==1
        and set(seams[0].ends)=={End(claim.endpoint.section_id,1-claim.endpoint.end),claim.external},
        'Owned prefix lost its real logical junction or two-ended seam')
    trace=RouteTopology(context.network).trace(claim.external)
    _route_pieces(context,route,trace.traversals)
    piece=route.pieces[0];external=End(piece.section.id,0 if piece.forward else 1)
    _require(external==claim.external and route.mouth_piece is piece
        and tuple(route.absorbed_sites)==(claim.seam,),
        'Owned prefix changed its actual external contact incidence')
    profile=context[Key('profile',piece.section.id)].oriented(external.end)
    offset=piece.section.offset*(1 if piece.forward else -1)
    _require(route.mouth_profile==profile and arm.profile==profile
        and arm.width==piece.section.width and arm.offset==offset and arm.length==route.length,
        'Owned prefix changed its actual external dimensions or profile')
    for p in route.pieces:
        end=0 if p.forward else 1
        _require(p.section.width==arm.width and p.section.offset*(1 if p.forward else -1)==arm.offset
            and context[Key('profile',p.section.id)].oriented(end)==profile,
            'Owned prefix crosses an incompatible external width or profile')
    first=prefix.frame(0. if claim.endpoint.end==0 else prefix.length)[0]
    last=prefix.frame(prefix.length if claim.endpoint.end==0 else 0.)[0]
    _require(np.linalg.norm(last-arm.origin)<=1e-8 and np.isfinite(first).all(),
        'Owned prefix does not reach its actual external frame')
    return claim,prefix,node,external


def _intervals(context,route,station):
    result={}
    for (sid,end),at in route.trims(station).items():
        path=context[Key('path',sid)];length=path.length
        _require(math.isfinite(at) and 0.<=at<=length,'Owned prefix extrapolated a native source interval')
        first,last=(0.,at) if end==0 else (at,length)
        if last>first:result.setdefault(sid,[]).append((float(first),float(last)))
    return {sid:tuple(sorted(rows)) for sid,rows in sorted(result.items())}


def _mouth_intervals(context,group):
    result={claim.section_id:[(0.,context[Key('path',claim.section_id)].length)] for claim in group.claims}
    for (sid,end),mouth in group.physical_layout.mouths.items():
        length=context[Key('path',sid)].length;at=mouth.station
        _require(math.isfinite(at) and 0.<=at<=length,'Prepared mouth escaped an actual native Section')
        low,high=(0.,at) if end==0 else (at,length)
        if high>low:result.setdefault(sid,[]).append((low,high))
    merged={}
    for sid,rows in result.items():
        intervals=[]
        for low,high in sorted(rows):
            if intervals and low<=intervals[-1][1]:intervals[-1]=(intervals[-1][0],max(high,intervals[-1][1]))
            else:intervals.append((low,high))
        merged[sid]=tuple(intervals)
    return merged


def _same_arm(a,b):
    if (type(a) is not Approach or type(b) is not Approach or a.source.id!=b.source.id
            or a.end!=b.end or a.width!=b.width or a.offset!=b.offset or a.profile!=b.profile
            or type(a.path) is not type(b.path)):return False
    if type(a.path) is OwnedPrefixRoute and a.path.claim is not b.path.claim:return False
    return (len(a.path.pieces)==len(b.path.pieces)
        and all(x.path is y.path and x.forward==y.forward for x,y in zip(a.path.pieces,b.path.pieces)))


def verified_owned_prefix_ranges(context,group,arm,station):
    """Fresh finite whole-prefix/final-mouth proof, without owned_ranges()."""
    from .compound_coverage import CompoundSourceCoverage
    from .compound_ports import CompoundPortPlanner
    from .source_strip_containment import certify_intervals
    claim,prefix,node,external=_prefix_arm(context,arm)
    _require(claim.junction in group.members and group.physical_layout is not None,
        'Owned prefix is outside its actual physical compound')
    matches=[(candidate,float(at)) for candidate,at in zip(group.physical_layout.approaches,
        group.physical_layout.outward_stations) if _same_arm(candidate,arm)]
    _require(len(matches)==1,'Owned prefix lost its exact final external Approach')
    _,mouth=matches[0]
    _require(math.isfinite(station) and 0.<=station<=mouth<=arm.length,
        'Owned-prefix contact exceeds its original final profile plane')
    final=_intervals(context,arm.path,mouth);requested=_intervals(context,arm.path,station)
    _require(final.get(prefix.section.id)==((0.,prefix.length),),
        'Owned-prefix range omitted its whole original Section')
    for (sid,end),at in arm.path.trims(mouth).items():
        stored=group.physical_layout.mouths.get((sid,end))
        _require(stored is not None and stored.station==at,
            'Owned prefix changed an original complete source/external mouth alias')
    _require(prefix.section.id in CompoundSourceCoverage(context,group,{}).fully_owned,
        'Owned prefix lacks complete original source-interval ownership')
    outer=tuple(loop for loop in group.loops if not loop.hole)
    _require(len(outer)==1,'Owned prefix needs its actual single exterior')
    try:
        CompoundPortPlanner(context,(node,),{}).check_support({prefix.section.id},outer[0])
        width=certify_intervals(context,group,final)
    except ValueError as failure:
        raise OwnedPrefixCornerAuthorityError('Owned prefix failed fresh complete native support: '+str(failure)) from failure
    _require(all(any(a<=low<=high<=b for a,b in final.get(sid,()))
        for sid,rows in requested.items() for low,high in rows),
        'Owned-prefix contact intervals escape its certified final mouth')
    return OwnedPrefixRangeProof(claim,tuple(requested.items()),float(station),mouth,tuple(final.items()),width)


def _ordinary_arm(context,arm):
    _require(type(arm) is Approach and type(arm.path) is Route and arm.path.pieces,
        'Prefix corner requires unchanged ordinary source arms beside its typed prefixes')
    _require(matches_current_section(context,arm.source) and arm.path.section is arm.source
        and arm.source is arm.path.pieces[0].section
        and arm.path.pieces[0].forward==(arm.end==0)
        and arm.width==arm.source.width and arm.offset==arm.source.offset*(1 if arm.end==0 else -1)
        and arm.profile==context[Key('profile',arm.source.id)].oriented(arm.end)
        and not hasattr(arm.path,'mouth_piece') and not getattr(arm.path,'absorbed_sites',()),
        'Prefix corner changed an ordinary logical source frame')
    _require(all(p.path is context[Key('path',p.section.id)] and matches_current_section(context,p.section)
        for p in arm.path.pieces),'Prefix corner borrowed an ordinary native Path')


def validate_owned_prefix_corner(context,group,node,corner):
    """Separate real logical sector/radius from proved actual contact frames."""
    from .compound_displaced_contacts import SourceContactWalk
    from .compound_contact_routes import regular_contact
    from .corner_radii import CornerRadiusPolicy
    from .ownership import footprint_key
    from .source_strip_containment import certify_intervals
    _require(any(node is actual for actual in context.network.intersections) and node.key in group.members,
        'Prefix corner lost its real original junction')
    support=getattr(corner,'source_contact_support',None);walks=getattr(corner,'source_contact_walks',())
    _require(isinstance(support,dict) and support.get('policy')=='EXPLICIT_WHOLE_OWNED_PREFIX_CONTACTS'
        and support.get('original_junction')==node.key.text() and support.get('endpoint_extrapolation') is False
        and type(walks) is tuple and len(walks)==2 and all(type(walk) is SourceContactWalk for walk in walks)
        and corner.center is not None,'Prefix corner lacks its explicit original producer authority')
    arms=(corner.first,corner.second);logical=[];contacts=[];range_proofs=[];prefix_count=0
    known=_mouth_intervals(context,group);contact_ranges={}
    for arm,walk,station in zip(arms,walks,corner.stations):
        _require(walk.approach is arm,'Prefix corner changed its recorded actual contact Approach')
        endpoint=End(arm.source.id,arm.end);logical.append(endpoint)
        _require(endpoint in node.ends and math.isfinite(station) and 0.<=station<=arm.length
            and regular_contact(arm,station),'Prefix corner left a real finite source incidence/contact')
        if type(arm.path) is OwnedPrefixRoute:
            claim,prefix,actual_node,external=_prefix_arm(context,arm);prefix_count+=1
            _require(actual_node is node and claim.endpoint==endpoint and len(walk.prefix)==1
                and walk.prefix[0].path is prefix and walk.prefix[0].forward==(endpoint.end==0)
                and walk.original.path is prefix and walk.original.source is prefix.section
                and walk.original.end==endpoint.end and walk.original.width==prefix.section.width
                and walk.original.offset==prefix.section.offset*(1 if endpoint.end==0 else -1)
                and walk.original.profile==context[Key('profile',prefix.section.id)].oriented(endpoint.end),
                'Prefix corner lost its real original centreline closure')
            proof=verified_owned_prefix_ranges(context,group,arm,station);range_proofs.append(proof);contacts.append(external)
            ranges=dict(proof.ranges)
        else:
            _ordinary_arm(context,arm)
            _require(walk.original is arm and not walk.prefix,'Ordinary prefix-corner arm acquired a fabricated source walk')
            contacts.append(endpoint);ranges=_intervals(context,arm.path,station)
            _require(all(any(a<=low<=high<=b for a,b in known.get(sid,()))
                for sid,rows in ranges.items() for low,high in rows),'Prefix corner ordinary contact exceeds actual prepared ownership')
        for sid,rows in ranges.items():contact_ranges.setdefault(sid,[]).extend(rows)
    _require(prefix_count>0 and logical[0]!=logical[1]
        and np.linalg.norm(walks[0].original.origin-walks[1].original.origin)<=1e-8,
        'Prefix corner lacks one real logical original junction')
    layout=context.products.get(footprint_key(node.key));ordered=getattr(layout,'approaches',())
    ends=tuple(End(arm.source.id,arm.end) for arm in ordered)
    _require(len(ends)==len(node.ends) and set(ends)==set(node.ends)
        and all(sum(_same_arm(arm,candidate) for candidate in ordered)==1 for arm in arms),
        'Prefix corner lost the one-to-one original producer Approach layout')
    proxy=tuple(End(arm.source.id,arm.end) for arm in sorted(ordered,
        key=lambda arm:math.atan2(arm.direction[1],arm.direction[0])%math.tau))
    _require(proxy[(proxy.index(logical[0])+1)%len(proxy)]==logical[1],
        'Prefix corner changed its actual proxy-layout cyclic sector')
    requested,mode=node.corner_radius(logical[0])
    _require(corner.requested_radius==requested and corner.radius_mode==mode,
        'Prefix corner changed its real logical radius policy')
    allowed=tuple(CornerRadiusPolicy.contact_candidates(*arms,requested,mode));attempts=tuple(corner.radius_attempts)
    _require(attempts and attempts==allowed[:len(attempts)] and corner.radius==attempts[-1]
        and (mode!='STRICT' or corner.radius==requested),'Prefix corner changed its bounded original radius candidates')
    a,b=(arm.direction for arm in arms)
    gap=math.atan2(float(a[0]*b[1]-a[1]*b[0]),float(a@b))%math.tau
    side=1. if gap<math.pi else -1.
    _require(math.isfinite(corner.sweep) and 0.<abs(corner.sweep)<=math.pi+1e-7,
        'Prefix corner selected a different finite circular sector')
    _require(abs(corner.gap-gap)<=1e-9 and corner.sweep*side<=1e-7,
        'Prefix corner changed its actual proxy-frame circular orientation')
    angles=corner.start_angle+corner.sweep*np.array((0.,1.));circle=corner.center+corner.radius*np.column_stack((np.cos(angles),np.sin(angles)))
    for i,(arm,station,side) in enumerate(zip(arms,corner.stations,(-1,1))):
        point,tangent=arm.jet(station,arm.offset+side*.5*arm.width);radial=point-corner.center
        _require(np.isfinite(point).all() and np.isfinite(tangent).all() and np.linalg.norm(tangent)>1e-10
            and np.linalg.norm(point-circle[i])<=1e-8 and abs(np.linalg.norm(radial)-corner.radius)<=1e-7
            and abs(float(radial@tangent))/np.linalg.norm(tangent)<=1e-7,
            'Prefix corner changed its actual finite circular tangencies')
    try:width=certify_intervals(context,group,contact_ranges)
    except ValueError as failure:
        raise OwnedPrefixCornerAuthorityError('Prefix corner contact support failed: '+str(failure)) from failure
    def native_angle(end):
        path=context[Key('path',end.section_id)];t=path.frame(0. if end.end==0 else path.length)[1]*(1 if end.end==0 else -1)
        return math.atan2(t[1],t[0])%math.tau
    native=tuple(sorted(node.ends,key=native_angle))
    return OwnedPrefixCornerAuthority(node.key,tuple(logical),tuple(contacts),walks,tuple(range_proofs),attempts,native,proxy,width)
