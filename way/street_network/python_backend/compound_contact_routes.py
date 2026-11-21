"""Finite contact searches across original incidences inside a compound.

These piecewise paths are construction supports, not smooth road seams. Their
original centerlines and offset frames stay unchanged. A positive fillet may
contact a later source Section; its complete external mouth must then lie past
that contact. Union/coverage proofs still reject exposed preparation seams.
"""
import math
import numpy as np
from .model import End,Key
from .routes import Piece,Route


def extend(context, approach, members, consumed, *, maximum_pieces=32):
    """Follow unique profile-compatible source continuations at owned nodes.

    Street IDs need not agree: a cross street can divide one physical through
    carriageway into separate Streets. Only an explicit original incidence may
    connect them. Ambiguous forks, sharp turns, scope/selection cuts, changed
    profiles and absorbed approaches cannot extend a contact search.
    """
    from .intersections import Approach
    if (approach.source.id not in consumed or not hasattr(approach.path,'pieces')
            or getattr(approach.path,'absorbed_sites',())
            or getattr(approach.path,'mouth_piece',None) is not None):return approach
    pieces=list(approach.path.pieces)
    if not pieces or any(p.section.id not in consumed for p in pieces):return approach
    incident={end:j for j in (*context.network.intersections,*context.network.connections) for end in j.ends}
    scoped={End(s.id,e) for s in context.network.sections for e in (*s.scope_ends,*s.selection_ends)}
    seen={p.section.id for p in pieces};through=[]
    cosine=math.cos(math.radians(20.))
    while len(pieces)<maximum_pieces:
        previous=pieces[-1]
        exit=End(previous.section.id,1 if previous.forward else 0)
        site=incident.get(exit)
        if site is None or site.key not in members or exit in scoped:break
        candidates=[]
        for end in site.ends:
            if end==exit or end.section_id in seen or end in scoped:continue
            source=context.network.section_index[end.section_id]
            if source.interpolation!=approach.source.interpolation:continue
            profile=context[Key('profile',source.id)].oriented(end.end)
            offset=source.offset*(1. if end.end==0 else -1.)
            if (profile!=approach.profile or abs(source.width-approach.width)>1e-8
                    or abs(offset-approach.offset)>1e-8):continue
            piece=Piece(context[Key('path',source.id)],end.end==0)
            a,ta,_=previous.frame(previous.path.length)
            b,tb,_=piece.frame(0.)
            if np.linalg.norm(a-b)>1e-8 or float(ta@tb)<cosine:continue
            candidates.append(piece)
        if len(candidates)!=1:break
        following=candidates[0];pieces.append(following);seen.add(following.section.id);through.append(site.key)
        if following.section.id not in consumed:break
    if not through:return approach
    path=Route(tuple(pieces))
    path.contact_extension_through=tuple(through)
    return Approach(path,approach.end,context[Key('profile',approach.source.id)])


def regular_contact(approach, station):
    """A root must use one original smooth frame, away from an internal seam."""
    return all(abs(station-float(s))>2e-4 for s in approach.path.stations[1:-1])


def outward_station(approach, contacts, *, direct_contacts=()):
    """Retain direct original coordinates before mapping other source aliases.

    A native cut converted back through a reversed Route can gain an ulp. A
    contact already expressed on this exact Approach supplies its authoritative
    original coordinate; only other native contacts require an inverse map.
    """
    maximum=0.
    direct={}
    for source,station in direct_contacts:
        if source is not approach:continue
        if not math.isfinite(station) or not 0.<=station<=approach.length:
            raise ValueError('Direct source contact exceeds its original finite Approach')
        maximum=max(maximum,float(station))
        for key,at in source.path.trims(station).items():
            direct.setdefault(key,set()).add(at)
    for i,piece in enumerate(approach.path.pieces):
        end=0 if piece.forward else 1
        for at in contacts.get((piece.section.id,end),()):
            if at in direct.get((piece.section.id,end),()):continue
            local=at if piece.forward else piece.path.length-at
            maximum=max(maximum,float(approach.path.stations[i]+local))
    return maximum
