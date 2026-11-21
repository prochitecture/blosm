"""Explicit short source prefixes wholly owned by one physical junction.

An authored dimensional seam immediately beside a junction can make a forced
shared endpoint tangent impossible. The original prefix still needs a certified
centerline. Its surface may be deferred only to the incident junction, whose
complete rounded exterior must contain that prefix and the absorbed seam.
"""
from dataclasses import dataclass
import math
from .model import End, Key
from .routes import Route, RouteTopology


@dataclass(frozen=True)
class OwnedPrefix:
    junction: object
    endpoint: End
    seam: object
    external: End


def candidates(network):
    """Bounded topology candidates; activation still requires a failed join."""
    from .profiles import memberships
    bundled=memberships(network)
    physical={end:j.key for j in network.intersections if len(j.ends)>2 for end in j.ends}
    result={}
    for seam in (*network.connections,*network.intersections):
        if len(seam.ends)!=2:continue
        for inner,outer in (seam.ends,seam.ends[::-1]):
            prefix=network.section_index[inner.section_id]
            external=network.section_index[outer.section_id]
            start=End(prefix.id,1-inner.end)
            if (start not in physical or prefix.id in bundled
                    or prefix.interpolation!=3 or external.interpolation!=3
                    or prefix.selection_ends or prefix.scope_ends
                    or abs(prefix.datum-external.datum)>1e-8):
                continue
            if external.id in bundled:
                bundle,role,_=bundled[external.id]
                member=bundle.members[role]
                first,last=member.path[0],member.path[-1]
                terminals={End(first.section_id,0 if first.forward else 1),
                           End(last.section_id,1 if last.forward else 0)}
                if not member.complete or outer not in terminals:
                    continue  # An interior Bundle seam needs extent ownership.
            length=sum(math.dist(a[:2],b[:2]) for a,b in zip(prefix.points,prefix.points[1:]))
            available=sum(math.dist(a[:2],b[:2]) for a,b in zip(external.points,external.points[1:]))
            if length>max(prefix.width,external.width) or available<=max(2.*length,external.width):
                continue
            if prefix.id in result:
                raise ValueError('Short source prefix has ambiguous physical junction ownership')
            result[prefix.id]=OwnedPrefix(physical[start],start,seam.key,outer)
    return result


class OwnedPrefixRoute(Route):
    """The real external path, with original consumed-prefix trim aliases."""
    def __init__(self, context, claim):
        prefix=context[Key('path',claim.endpoint.section_id)]
        external=Route.build(context,RouteTopology(context.network).trace(claim.external))
        super().__init__(external.pieces,terminal_free=external.terminal_free)
        self.section=prefix.section
        self.prefix=prefix
        self.claim=claim
        self.mouth_piece=self.pieces[0]
        self.mouth_profile=context[Key('profile',claim.external.section_id)].oriented(claim.external.end)
        self.absorbed_sites=(claim.seam,)

    def trims(self,distance):
        claim=self.claim
        result={(claim.endpoint.section_id,claim.endpoint.end):
            self.prefix.length if claim.endpoint.end==0 else 0.}
        result.update(super().trims(distance))
        return result


def approaches(context,junction,original):
    from .intersections import Approach
    result=[];claims=[]
    for approach in original:
        claim=getattr(context[Key('path',approach.source.id)],'owned_prefix',None)
        if claim is None or claim.junction!=junction.key or claim.endpoint.end!=approach.end:
            result.append(approach);continue
        path=OwnedPrefixRoute(context,claim)
        result.append(Approach(path,approach.end,context[Key('profile',approach.source.id)]))
        claims.append(claim)
    return result,tuple(claims)


def contact_walks(context,junction,corner):
    """Close displaced fillets through their actual wholly owned prefix.

    The external route retains its real frame, width and half profile. Its
    closure follows the original finite prefix centerline back to the actual
    junction; a chord between displaced origins is never invented.
    """
    import numpy as np
    from dataclasses import replace
    from .intersections import Approach
    from .routes import Piece
    from .compound_displaced_contacts import SourceContactWalk
    if np.linalg.norm(corner.first.origin-corner.second.origin)<=1e-7:
        return corner
    walks=[];changed=False
    seams={site.key:site for site in (*context.network.connections,*context.network.intersections)
           if len(site.ends)==2}
    for approach in (corner.first,corner.second):
        path=approach.path
        if not isinstance(path,OwnedPrefixRoute):
            walks.append(SourceContactWalk(approach,approach));continue
        claim=path.claim
        seam=seams.get(claim.seam)
        expected=End(approach.source.id,approach.end)
        if (claim.junction!=junction or claim.endpoint!=expected or seam is None
                or set(seam.ends)!={End(expected.section_id,1-expected.end),claim.external}
                or not path.pieces or End(path.pieces[0].section.id,int(not path.pieces[0].forward))!=claim.external):
            raise ValueError('Source corner lost its explicit owned-prefix incidence')
        original=Approach(path.prefix,approach.end,context[Key('profile',expected.section_id)])
        prefix=Piece(path.prefix,expected.end==0)
        if np.linalg.norm(prefix.frame(path.prefix.length)[0]-approach.origin)>1e-8:
            raise ValueError('Owned prefix does not reach its actual external source frame')
        walks.append(SourceContactWalk(original,approach,(prefix,)));changed=True
    if not changed:
        return corner
    if np.linalg.norm(walks[0].original.origin-walks[1].original.origin)>1e-8:
        raise ValueError('Owned-prefix contact closure lacks one original source junction')
    result=replace(corner)
    result.source_contact_walks=tuple(walks)
    result.source_contact_support=dict(policy='EXPLICIT_WHOLE_OWNED_PREFIX_CONTACTS',
        original_junction=junction.text(),endpoint_extrapolation=False,
        walks=[walk.report(station) for walk,station in zip(walks,corner.stations)])
    return result


def certify(context,junction,layout,claims,*,group=None):
    """Prove source containment and whole ownership before any surface exists."""
    if not claims:return
    from .compound_ports import CompoundPortPlanner
    from .compounds import sample_footprint
    from .regions import RegionUnion
    if group is None:
        footprint=sample_footprint(junction.key,layout)
        loops=RegionUnion().build((footprint.region,))
        if len(loops)!=1 or loops[0].hole:
            raise ValueError('Owned source prefix needs one simple junction exterior')
    else:
        if group.physical_layout is not layout or junction.key not in group.members:
            raise ValueError('Owned prefix certificate needs its actual prepared source group')
        loops=tuple(loop for loop in group.loops if not loop.hole)
        if len(loops)!=1:
            raise ValueError('Owned source prefix needs one simple junction exterior')
    consumed={claim.endpoint.section_id for claim in claims}
    for claim in claims:
        if group is not None:
            source=context[Key('path',claim.endpoint.section_id)]
            if (claim.junction not in group.members or getattr(source,'owned_prefix',None)!=claim):
                raise ValueError('Owned source prefix lost its real original junction claim')
        mouth=layout.mouths.get((claim.endpoint.section_id,claim.endpoint.end))
        length=context[Key('path',claim.endpoint.section_id)].length
        required=length if claim.endpoint.end==0 else 0.
        if mouth is None or abs(mouth.station-required)>1e-8:
            raise ValueError('Owned source prefix has no complete source-interval claim')
        if (claim.external.section_id,claim.external.end) not in layout.mouths:
            raise ValueError('Owned source prefix lost its complete external profile')
    CompoundPortPlanner(context,(junction,),{}).check_support(consumed,loops[0])
    if group is not None:
        from .compound_coverage import CompoundSourceCoverage
        from .source_strip_containment import certify as certify_widths
        if not consumed<=CompoundSourceCoverage(context,group,{}).fully_owned:
            raise ValueError('Prepared source group does not wholly own its original prefix')
        layout.owned_prefix_width_proof=certify_widths(context,group,consumed)
