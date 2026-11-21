"""Declared physical bank coordinates for an explicitly owned source prefix.

The logical Section remains the real junction incidence. Coordinate pieces and
profile jets belong to the unchanged external native Route, never to that
logical Section's different width/profile.
"""
import math
import numpy as np
from .model import End
from .paths import Path
from .intersections import Approach
from .owned_prefixes import OwnedPrefix,OwnedPrefixRoute


def validate_declaration(path,approach,context=None):
    """Exact producer type and source identities; no geometric bank inference."""
    from .source_curb_coordinates import CurbCoordinatesUnavailable
    def require(value,message):
        if not value:raise CurbCoordinatesUnavailable(message)
    require(type(path) is OwnedPrefixRoute and type(approach) is Approach
        and approach.path is path,'Prefix curb needs its exact original physical Approach')
    claim=path.claim;prefix=path.prefix
    require(type(claim) is OwnedPrefix and type(prefix) is Path
        and getattr(prefix,'owned_prefix',None) is claim
        and path.section is prefix.section and approach.source is prefix.section
        and approach.end==claim.endpoint.end and claim.endpoint.section_id==prefix.section.id,
        'Prefix curb changed its actual logical Path/claim identity')
    require(path.pieces and all(type(piece.path) is Path and piece.section is piece.path.section
        for piece in path.pieces),'Prefix curb needs declared original native Route pieces')
    first=path.pieces[0];external=End(first.section.id,0 if first.forward else 1)
    metric=np.r_[0.,np.cumsum([piece.path.length for piece in path.pieces])]
    require(external==claim.external and path.mouth_piece is first
        and tuple(path.absorbed_sites)==(claim.seam,)
        and np.array_equal(path.stations,metric) and path.length==float(metric[-1]),
        'Prefix curb changed its finite actual external Route')
    require(approach.width==first.section.width
        and approach.offset==first.section.offset*(1 if first.forward else -1)
        and approach.profile==path.mouth_profile and approach.length==path.length,
        'Prefix curb changed its actual external physical profile/frame')
    require(all(piece.section.width==approach.width
        and piece.section.offset*(1 if piece.forward else -1)==approach.offset
        for piece in path.pieces),'Prefix curb crosses incompatible native physical dimensions')
    if context is not None:
        from .owned_prefix_corner_authority import _prefix_arm
        _prefix_arm(context,approach)
    return True


def validate_curve(resolver,source,index,curve,first,last):
    """Fresh current-group proof before the existing native-coordinate rule.

    Ordinary Path/Route ledgers pass through unchanged. An explicit prefix
    ledger additionally proves whole-prefix and final-mouth support and that
    every declared physical bank piece lies in those finite native intervals.
    """
    from .source_curb_coordinates import CurbCoordinates
    from .source_curb_conditions import SourceCurbConditionUnavailable
    coordinates=getattr(curve,'source_coordinates',None)
    if type(coordinates) is not CurbCoordinates:
        raise SourceCurbConditionUnavailable('Prefix curb lost its original coordinate ledger')
    if type(coordinates.path) is not OwnedPrefixRoute:return None
    def require(value,message):
        if not value:raise SourceCurbConditionUnavailable(message)
    require(source.region.owner.kind=='cycle_port' and index in (0,1)
        and len(source.layout.approaches)==1
        and coordinates.approach is source.layout.approaches[0]
        and resolver.sources.get(source.region.owner) is source
        and source.layout.boundaries[index] is curve,
        'Prefix bank ledger lost its exact prepared source port')
    require(math.isfinite(first) and math.isfinite(last) and 0.<=first<last<=1.+1e-9,
        'Prefix bank needs a finite increasing source interval')
    approach=coordinates.approach
    validate_declaration(coordinates.path,approach,resolver.context)
    coordinates.validate(curve.points,resolver.context)
    station=float(source.layout.outward_stations[0]);values=coordinates.stations
    require(np.isfinite(station) and tuple(values[[0,-1]]) in ((0.,station),(station,0.))
        and (np.all(np.diff(values)>0.) or np.all(np.diff(values)<0.))
        and len(set(coordinates.sides))==1,
        'Prefix bank lost its complete declared monotone physical station map')
    from .owned_prefix_corner_authority import verified_owned_prefix_ranges
    # A resolver owns one immutable prepared group, like its existing rule and
    # native-knot caches. Prove this actual port once in that resolver, rather
    # than repeating its complete-width proof for every exposed chord. No proof
    # is reused across a new resolver, group, Approach or final profile plane.
    cache=getattr(resolver,'owned_prefix_source_support',None)
    if cache is None:cache=resolver.owned_prefix_source_support={}
    key=(source.region.owner,approach,station)
    cached=cache.get(key)
    if cached is None or cached[0] is not resolver.group:
        proof=verified_owned_prefix_ranges(resolver.context,resolver.group,approach,station)
        cache[key]=(resolver.group,proof)
    else:proof=cached[1]
    require(proof.requested_station==proof.mouth_station,
        'Prefix source port no longer reaches its actual complete profile plane')
    owned=dict(proof.ranges)
    for edge in range(len(values)-1):
        for piece in coordinates.edge_sources(edge):
            low,high=sorted((piece.first_station,piece.last_station))
            require(piece.authored_side in (-1,1) and any(a<=low<=high<=b
                for a,b in owned.get(piece.path.section.id,())),
                'Prefix bank left its certified actual native physical support')
    return proof
