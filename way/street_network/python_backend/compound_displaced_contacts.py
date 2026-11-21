"""Finite circles on actual source curbs across an owned dimensional corner.

Each contact uses its own original Approach and width. A source-incidence
curb walk supplies the internal centreline closure; no uniform-width Route,
invented node bridge or endpoint extension represents a width change.
"""
from dataclasses import dataclass,replace
from types import SimpleNamespace
import math
import numpy as np
from .model import End,Key
from .routes import RouteTopology
from .corner_radii import CornerRadiusPolicy,CornerContactError
from .compound_contact_routes import regular_contact
from .compounds import SampledFootprint,RimSegment
from .finite_free_caps import PreparedCapCurb
from .regions import Region,RegionUnion,signed_area


@dataclass(frozen=True)
class SourceContactWalk:
    original: object
    approach: object
    prefix: tuple=()

    @property
    def changed_width(self):
        return any(abs(piece.section.width-self.approach.width)>1e-8 for piece in self.prefix)

    @property
    def length(self):return sum(piece.path.length for piece in self.prefix)

    def report(self,station):
        return dict(original_incidence=dict(section=self.original.source.id,end=self.original.end),
            contact_incidence=dict(section=self.approach.source.id,end=self.approach.end),
            actual_contact_width=self.approach.width,contact_outward_station=float(station),
            original_centerline_closure=[dict(section=p.section.id,forward=p.forward,width=p.section.width)
                for p in self.prefix],actual_contact_aliases=[dict(section=sid,end=end,station=float(s))
                for (sid,end),s in self.approach.path.trims(station).items()],endpoint_extrapolation=False)


def _ordinary(approach):
    return (hasattr(approach.path,'pieces') and approach.length>0. and not approach.profile.half
        and abs(approach.offset)<1e-8 and not getattr(approach.path,'absorbed_sites',())
        and getattr(approach.path,'mouth_piece',None) is None
        and getattr(approach.path,'deferred_profile_reserve',None) is None
        and getattr(approach.path,'minimum_mouth',0.)<=1e-9)


def _walks(context,original,side,members,consumed,maximum_steps=8):
    from .ownership import footprint_key
    topology=RouteTopology(context.network)
    sites={j.key:j for j in context.network.intersections}
    result=[SourceContactWalk(original,original)];seen=set();prefix=[];current=original
    for _ in range(maximum_steps):
        if (not _ordinary(current) or any(p.section.id not in consumed or p.section.id in seen
                for p in current.path.pieces)):break
        route=topology.trace(End(current.source.id,current.end))
        if (route.terminal_site not in members or route.terminal_site not in sites
                or route.terminal_scope or route.terminal_selection
                or tuple((p.section.id,p.forward) for p in current.path.pieces)
                   !=tuple((v.section_id,v.forward) for v in route.traversals)):break
        last=route.traversals[-1];reverse=End(last.section_id,1 if last.forward else 0)
        layout=context.products.get(footprint_key(route.terminal_site))
        if layout is None:break
        approaches=tuple(a for a in layout.approaches if _ordinary(a))
        ordered=sorted(approaches,key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
        ends=[End(a.source.id,a.end) for a in ordered]
        if set(ends)!=set(sites[route.terminal_site].ends) or reverse not in ends:break
        at=ends.index(reverse);following=ordered[(at+(1 if side>0 else -1))%len(ordered)]
        local_first=reverse if side>0 else End(following.source.id,following.end)
        if sites[route.terminal_site].corner_radius(local_first)[1]!='STYLE':break
        if (following.source.id in seen or following.source.id in {p.section.id for p in current.path.pieces}
                or following.profile!=original.profile):break
        end=current.frame(current.length)[0]
        if np.linalg.norm(end-following.origin)>1e-8:break
        prefix.extend(current.path.pieces);seen.update(p.section.id for p in current.path.pieces)
        current=following
        result.append(SourceContactWalk(original,current,tuple(prefix)))
    return tuple(result)


def _centerline(walk,station,spacing,tolerance,maximum):
    from .corner_patch import _stations
    points=[];previous=None
    for piece in walk.prefix:
        at=piece.path.samples(0.,piece.path.length,spacing,tolerance=tolerance,max_samples=maximum)
        local=piece.frame(at)[0]
        if previous is not None and np.linalg.norm(previous-local[0])>1e-8:
            raise ValueError('Displaced source contact lacks an original centreline incidence')
        points.extend(local[:-1])
        previous=local[-1]
        if len(points)>maximum:raise ValueError('Displaced source patch exceeds its finite sample budget')
    at=_stations(walk.approach,station,spacing,tolerance,maximum)
    local=walk.approach.frame(at)[0]
    if points and np.linalg.norm(walk.prefix[-1].frame(walk.prefix[-1].path.length)[0]-local[0])>1e-8:
        raise ValueError('Displaced source contact leaves its actual source prefix')
    points.extend(local)
    if len(points)>maximum or np.linalg.norm(points[0]-walk.original.origin)>1e-8:
        raise ValueError('Displaced source contact lost its original source node')
    return np.asarray(points)


def patch(owner,corner,spacing,tolerance,*,max_samples=8192):
    """One simple circle patch closed only through explicit source centrelines."""
    if (not all(math.isfinite(value) and value>0. for value in (spacing,tolerance))
            or type(max_samples) is not int or max_samples<4):
        raise ValueError('Invalid bounded displaced source patch sampling request')
    walks=getattr(corner,'source_contact_walks',None)
    if walks is None or corner.center is None or len(walks)!=2:
        raise ValueError('Displaced patch needs two explicit actual-source contact walks')
    if np.linalg.norm(walks[0].original.origin-walks[1].original.origin)>1e-8:
        raise ValueError('Displaced patch lacks one original source junction')
    if any(not math.isfinite(s) or s<0. or s>a.length for a,s in zip((corner.first,corner.second),corner.stations)):
        raise ValueError('Displaced source patch cannot extrapolate an original contact')
    chains=[_centerline(w,s,min(.3,spacing),tolerance,max_samples) for w,s in zip(walks,corner.stations)]
    points=[chains[0][0]];kinds=[]
    def append(chain,kind):
        for point in chain:
            if np.linalg.norm(point-points[-1])<=1e-8:continue
            points.append(np.asarray(point));kinds.append(kind)
            if len(points)>max_samples:raise ValueError('Displaced source patch exceeds its finite sample budget')
    append(chains[0][1:],'internal')
    count=max(2,int(math.ceil(abs(corner.sweep)*corner.radius/min(spacing,math.sqrt(8.*tolerance*corner.radius)))))
    if count+len(points)>max_samples:raise ValueError('Displaced source patch exceeds its finite sample budget')
    angle=np.linspace(corner.start_angle,corner.start_angle+corner.sweep,count+1)
    arc=corner.center+corner.radius*np.column_stack((np.cos(angle),np.sin(angle)))
    append(arc[:1],'internal');append(arc[1:],'curb')
    append(chains[1][::-1],'internal')
    if np.linalg.norm(points[-1]-points[0])>1e-8:
        raise ValueError('Displaced source patch needs an invented closing bridge')
    points=np.asarray(points[:-1])
    if len(points)!=len(kinds):raise ValueError('Displaced source patch has unresolved boundary incidence')
    if signed_area(points)<0.:
        points=points[::-1];kinds=list(np.roll(np.asarray(kinds)[::-1],-1))
    start=next((i for i in range(len(points)) if kinds[i]!=kinds[i-1]),0)
    points=np.roll(points,-start,axis=0);kinds=kinds[start:]+kinds[:start]
    region=Region(owner,points)
    try:RegionUnion().build((region,))
    except ValueError as failure:
        support=getattr(corner,'source_contact_support',None)
        if (str(failure)!='A physical footprint boundary intersects or touches itself'
                or not isinstance(support,dict)
                or support.get('policy')!='EXPLICIT_WHOLE_OWNED_PREFIX_CONTACTS'):
            raise
        from .internal_contact_closure import normalize
        points,kinds,proof=normalize(points,kinds,corner,spacing,tolerance)
        corner.source_contact_support['internal_closure']=proof
        region=Region(owner,points)
    curves=[];labels=[];first=0
    while first<len(points):
        last=first+1
        while last<len(points) and kinds[last]==kinds[first]:last+=1
        curve=PreparedCapCurb(np.vstack((points[first:last],points[last%len(points)])))
        index=len(curves);curves.append(curve)
        labels.extend(RimSegment(kinds[first],index,float(a/curve.length),float(b/curve.length))
            for a,b in zip(curve.stations,curve.stations[1:]))
        first=last
    layout=SimpleNamespace(approaches=(),outward_stations=(),boundaries=tuple(curves),rows=len(points),mouths={},corners=(corner,))
    return SampledFootprint(region,tuple(labels),layout)


def retry(context,first,second,radius,mode,members,consumed):
    from .updates import drain_steps
    return drain_steps(retry_steps(context,first,second,radius,mode,members,consumed))


def retry_steps(context,first,second,radius,mode,members,consumed):
    if mode!='STYLE' or not all(_ordinary(a) for a in (first,second)):return None
    if np.linalg.norm(first.origin-second.origin)>1e-8:return None
    a=_walks(context,first,-1,members,consumed);b=_walks(context,second,1,members,consumed)
    pairs=[(x,y) for x in a for y in b if (x.prefix or y.prefix) and (x.changed_width or y.changed_width)]
    pairs=sorted(pairs,key=lambda pair:(pair[0].length+pair[1].length,
        pair[0].approach.source.id,pair[0].approach.end,pair[1].approach.source.id,pair[1].approach.end))[:16]
    if not pairs:return None
    from .intersections import _corner_steps
    gap=math.atan2(float(first.direction[0]*second.direction[1]-first.direction[1]*second.direction[0]),
        float(first.direction@second.direction))%math.tau
    side=1. if gap<math.pi else -1.;settings=context.network.settings
    tolerance=min(.001,.05*settings.curve_tolerance);attempted=[]
    for value in CornerRadiusPolicy.contact_candidates(first,second,radius,mode):
        attempted.append(float(value))
        for x,y in pairs:
            yield
            try:
                c=yield from _corner_steps(x.approach,y.approach,value,math.radians(settings.continuation_angle_degrees),side_override=side)
                if c.center is None or any(s<0. or s>a.length for a,s in zip((c.first,c.second),c.stations)):
                    continue
                if not all(regular_contact(a,s) for a,s in zip((c.first,c.second),c.stations)):continue
                c=replace(c,requested_radius=radius,radius_mode=mode,radius_attempts=tuple(attempted),
                    fit_policy='DISPLACED_ORIGINAL_DIMENSIONAL_SOURCE_CURBS')
                c.source_contact_walks=(x,y)
                c.source_contact_support=[walk.report(station) for walk,station in zip((x,y),c.stations)]
                patch(Key('displaced_contact_proof',0),c,settings.spacing,tolerance)
                return c
            except (CornerContactError,ValueError):continue
    return None
