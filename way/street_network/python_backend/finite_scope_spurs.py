"""A finite open scope cut beside a continuing ground carriageway.

The finite ground spur can end before even the smallest style gore fits.
Keep the continuing road's full incoming profile and its actual curb. A new
scope profile follows the spur's exact cut plane across the small remaining
gap to that curb. This is an explicit open grade handoff, not a tiny fillet,
a closed dead end, or an extension of the omitted elevated road.
"""
from dataclasses import replace
import math
import numpy as np
from .finite_containment import contains_polygon
from .finite_scope_ports import FiniteScopeAlias, FiniteScopeInventory
from .model import End, Key
from .paths import unit
from .source_mouths import SourceMouth


class ScopeSpurCrown:
    half=0

    def __init__(self,rise,datum,width,peak_u,outer_side):
        self.rise,self.datum,self.width=rise,datum,width
        self.peak_u,self.peak=peak_u,2.*peak_u-1.
        self.outer_side=outer_side

    def value(self,q):
        u=.5*(np.asarray(q)+1.)
        outer=(u<=self.peak_u) if self.outer_side<0 else (u>=self.peak_u)
        span=self.peak_u if self.outer_side<0 else 1.-self.peak_u
        t=(u-self.peak_u)/span
        old=self.rise*(1.-t*t)
        inner=(u-self.peak_u)/(1.-self.peak_u) if self.outer_side<0 else (self.peak_u-u)/self.peak_u
        h=inner**3*(10.+inner*(-15.+6.*inner))
        return np.where(outer,old,self.rise*(1.-h))

    def derivative(self,q):
        u=.5*(np.asarray(q)+1.)
        outer=(u<=self.peak_u) if self.outer_side<0 else (u>=self.peak_u)
        span=self.peak_u if self.outer_side<0 else 1.-self.peak_u
        old=-self.rise*(u-self.peak_u)/(span*span)
        length=(1.-self.peak_u) if self.outer_side<0 else self.peak_u
        t=(u-self.peak_u)/length if self.outer_side<0 else (self.peak_u-u)/length
        inner=(-1. if self.outer_side<0 else 1.)*.5*self.rise*30.*t*t*(1.-t)**2/length
        return np.where(outer,old,inner)

    def coordinate(self,side,v):return side+(self.peak-side)*np.asarray(v)
    def intervals(self,side,base):return max(2,2*int(round(abs(self.peak-side)*base/2.)))


class ScopeCurbLead:
    """The original wide curb, ending on the exact finite scope cut plane."""
    def __init__(self,corner,wide,start,finish,side,settings):
        self.corner,self.wide,self.start,self.finish,self.side=corner,wide,start,finish,side
        self.continuation=None
        samples=_source_stations(wide,max(start,finish),min(settings.spacing,.2),settings.curve_tolerance*.1)
        samples=samples[(samples>min(start,finish)+1e-10)&(samples<max(start,finish)-1e-10)]
        samples=np.r_[min(start,finish),samples,max(start,finish)]
        if finish<start:samples=samples[::-1]
        parameters=(samples-start)/(finish-start)
        # An even number of intervals is needed for a conforming quad rim.
        parameters=np.sort(np.r_[parameters,.5*(parameters[:-1]+parameters[1:])])
        self.initial_parameters=parameters
        points=self(parameters)
        self.length=float(np.linalg.norm(np.diff(points,axis=0),axis=1).sum())
        self.first,self.last=points[0],points[-1]

    def __call__(self,t):
        station=self.start+np.asarray(t)*(self.finish-self.start)
        return self.wide.point(station,self.wide.offset+self.side*.5*self.wide.width)


class FiniteScopeSpurPort:
    def __init__(self,spur,wide,cut_station,wide_station,cut_point,outer_side,closing_gap):
        self.source,self.end,self.length=spur.source,spur.end,0.
        self.original_approaches=(spur,)
        self.path=FiniteScopeInventory(spur.path.pieces)
        self.spur,self.wide=spur,wide
        self.cut_station,self.wide_station=cut_station,wide_station
        self.outer_side=outer_side
        self.closing_gap=closing_gap
        outer=spur.point(spur.length,outer_side*.5*spur.width)
        self.points=np.array((outer,cut_point) if outer_side<0 else (cut_point,outer))
        self.width=float(np.linalg.norm(self.points[1]-self.points[0]))
        self.offset=0.
        self.normal=unit(self.points[1]-self.points[0])
        self.origin=.5*(self.points[0]+self.points[1])
        self.direction=spur.direction.copy()
        original_center=spur.frame(spur.length)[0]
        peak=float((original_center-self.points[0])@self.normal/self.width)
        self.profile=ScopeSpurCrown(spur.profile.rise,spur.profile.datum,self.width,peak,outer_side)
        self.aliases=(FiniteScopeAlias(End(spur.source.id,spur.end),End(spur.source.id,1-spur.end),
            spur.length if spur.end==0 else 0.,spur.length),)

    def point(self,station,distance):
        if np.any(abs(np.asarray(station))>1e-10):
            raise ValueError('An open finite scope spur cannot be extrapolated')
        q=2.*np.asarray(distance)/self.width
        if np.any(abs(q)>1.+1e-10):raise ValueError('An open finite scope spur cannot leave its cut')
        return self.origin+np.asarray(distance)[...,None]*self.normal

    def physical_gradient(self,station,q):
        self.point(station,.5*self.width*np.asarray(q))
        return self.profile.derivative(q)[...,None]*self.normal/(.5*self.width)

    def mouth_parameters(self,radial,spacing):
        result=[]
        for a,b in ((0.,self.profile.peak_u),(self.profile.peak_u,1.)):
            n=max(2,int(math.ceil((b-a)*self.width/spacing)));n+=n%2
            result.extend(np.linspace(a,b,n+1)[:-1])
        return np.r_[result,1.][::-1]

    def source_mouths(self,station):
        self.point(station,0.)
        alias=self.aliases[0]
        return tuple(SourceMouth(alias.junction_end,alias.scope_station,float(side),self.profile.peak_u,
            float(side),.5) for side in (0,1))

    def report(self):
        alias=self.aliases[0]
        return dict(mode='OPEN_GROUND_SCOPE_SPUR',profile='original_outer_half_reconciled_inner_half',
            datum=self.profile.datum,rise=self.profile.rise,outline=self.points.tolist(),
            source_aliases=[dict(section=alias.junction_end.section_id,junction_end=alias.junction_end.end,
                scope_end=alias.scope_end.end,scope_station=alias.scope_station,consumed_interval=[0.,alias.length])],
            original_scope_point=self.spur.frame(self.spur.length)[0].tolist(),
            closing_gap=self.closing_gap,wide_section=self.wide.source.id,
            wide_curb_cut_station=self.cut_station,wide_complete_profile_station=self.wide_station,
            boundary_kind='OPEN_SCOPE_BOUNDARY',reconciled_inner_half=True,
            endpoint_extrapolation=False,source_containment='finite_ground_side_plus_original_wide_strip')

    def additional_claims(self,owner):
        from .compounds import SharedClaim
        claims=[]
        for index,piece in enumerate(self.wide.path.pieces):
            if self.wide_station>=self.wide.path.stations[index+1]-1e-9:
                claims.append(SharedClaim(piece.section.id,owner,owner,piece.path.length,0.))
        return claims

    def support_footprint(self,context,owner,layout):
        """Keep the original ground curb where an auxiliary lead falls inside.

        A nearly straight continuing source need not coincide with the other
        incoming road's negative tangent extension. Union its actual finite
        strip before choosing exposed curbs; never waive whole-strip support.
        """
        from .compounds import SampledFootprint,RimSegment
        from .regions import Region
        from .intersections import Corner
        settings=context.network.settings
        distance=self.wide_station
        stations=self.support_stations
        high=self.wide.point(stations,self.wide.offset+.5*self.wide.width)
        low=self.wide.point(stations,self.wide.offset-.5*self.wide.width)
        points=np.r_[high,low[::-1]]
        proxy=Corner(self.wide,self.wide,(0.,0.),None,1.,0.,0.,math.pi)
        curves=(ScopeCurbLead(proxy,self.wide,0.,distance,1.,settings),
                ScopeCurbLead(proxy,self.wide,distance,0.,-1.,settings))
        wide_index=layout.approaches.index(self.wide)
        u=stations/distance
        labels=[RimSegment('curb',0,float(a),float(b)) for a,b in zip(u,u[1:])]
        labels.append(RimSegment('mouth',wide_index,1.,0.))
        labels.extend(RimSegment('curb',1,float(1.-a),float(1.-b)) for a,b in zip(u[:0:-1],u[-2::-1]))
        labels.append(RimSegment('internal',0,0.,1.))
        source_layout=replace(layout,corners=(),boundaries=curves,finite_scope_ports=())
        return SampledFootprint(Region(Key('scope_ground_support',owner.id),points),tuple(labels),source_layout)

    def validate_group(self,context,loop):
        settings=context.network.settings
        stations=_source_stations(self.wide,self.wide_station,min(.2,settings.spacing),
            min(.00025,settings.curve_tolerance*.05))
        high=self.wide.point(stations,self.wide.offset+.5*self.wide.width)
        low=self.wide.point(stations,self.wide.offset-.5*self.wide.width)
        ribbon=np.r_[high,low[::-1]]
        # Actual original curb samples have a certified sub-millimetre chord
        # error. Surface curbs retain the ordinary published curve tolerance.
        tolerance=settings.curve_tolerance+1e-7
        if not contains_polygon(ribbon,loop.xy,tolerance=tolerance):
            raise ValueError('Open scope pool lost its original continuing ground strip')
        scope_origin=self.spur.frame(self.spur.length)[0]
        beyond=_clip_halfplane(loop.xy,scope_origin,self.spur.direction,tolerance)
        if len(beyond)>=3 and not contains_polygon(beyond,ribbon,tolerance=tolerance):
            raise ValueError('Open scope pool crosses the grade boundary outside its original wide ground strip')


def _source_stations(approach,distance,spacing,tolerance):
    stations=[]
    for index,piece in enumerate(approach.path.pieces):
        start=approach.path.stations[index]
        if distance<=start+1e-10:break
        length=min(piece.path.length,distance-start)
        path=piece.path
        if path.exact and len(path.section.points)!=2:
            raise ValueError('Open scope pool needs a separately certified Exact miter policy')
        k=0. if path.exact else path.curvature_bound
        reach=abs(approach.offset)+.5*approach.width
        if reach*k>=1.:raise ValueError('Open scope pool has a folded continuing curb')
        speed=1.+reach*k;curvature=k/(1.-reach*k)
        step=spacing
        if curvature>1e-12:step=min(step,math.sqrt(8.*tolerance/curvature)/speed)
        if length/step+len(stations)>8190:raise ValueError('Open scope pool exceeds its finite source sampling budget')
        a,b=(0.,length) if piece.forward else (path.length-length,path.length)
        local=path.samples(a,b,step,tolerance=tolerance)
        if not piece.forward:local=path.length-local[::-1]
        stations.extend((start+local)[:-1])
        stations.append(start+length)
    return np.asarray(sorted(set(stations)))


def _clip_halfplane(polygon,origin,normal,offset):
    values=(polygon-origin)@normal-offset
    output=[]
    for a,b,va,vb in zip(polygon,np.roll(polygon,-1,axis=0),values,np.roll(values,-1)):
        if va>=0.:output.append(a)
        if (va>=0.)!=(vb>=0.):output.append(a+(b-a)*va/(va-vb))
    points=np.asarray(output)
    if len(points)>1:points=points[np.r_[True,np.linalg.norm(np.diff(points,axis=0),axis=1)>1e-9]]
    if len(points)>1 and np.linalg.norm(points[-1]-points[0])<1e-9:points=points[:-1]
    return points


def _align_support_samples(layout,port,wide_first,settings):
    """Use identical chords on repeated pieces of the same original curb.

    Different independent samplings of one shallow curved source otherwise
    create artificial microscopic strips in a planar union. Keep all source
    knots and existing curb samples, then refine their common parameter set.
    """
    wide=port.wide;distance=port.wide_station
    outer=layout.boundaries[2 if wide_first else 1]
    if not hasattr(outer,'leads'):
        values=sorted({0.,distance,port.cut_station,*_source_stations(wide,distance,
            min(.2,settings.spacing),min(.00025,settings.curve_tolerance*.05))})
        values=np.asarray(values)
        values=np.sort(np.r_[values,.5*(values[:-1]+values[1:])])
        port.support_stations=values
        inner=layout.boundaries[0]
        values=values[(values>=port.cut_station-1e-9)&(values<=distance+1e-9)]
        inner.initial_parameters=np.sort((values-inner.start)/(inner.finish-inner.start))
        return
    lead_index=1 if wide_first else 0
    approach,parameters,cumulative,offset=outer.leads[lead_index]
    base=outer.lengths[2] if lead_index else 0.
    old=np.linspace(0.,1.,layout.rows+1)
    stations=list(_source_stations(wide,distance,min(.2,settings.spacing),
        min(.00025,settings.curve_tolerance*.05)))+[port.cut_station]
    nonwide=[]
    for t in old:
        at=t*outer.length-base
        if 0.<=at<=cumulative[-1]:
            s=float(np.interp(at,cumulative,parameters))
            if 0.<=s<=distance:stations.append(s);continue
        nonwide.append(float(t))
    stations=np.asarray(sorted(stations))
    stations=stations[np.r_[True,np.diff(stations)>1e-9]]
    stations=np.sort(np.r_[stations,.5*(stations[:-1]+stations[1:])])
    if len(stations)>2048:raise ValueError('Open scope curb alignment exceeds its source sample budget')
    port.support_stations=stations
    local=(np.interp(stations,parameters,cumulative) if parameters[-1]>parameters[0]
           else np.interp(stations,parameters[::-1],cumulative[::-1]))
    t=(base+local)/outer.length
    values=np.asarray(sorted([*nonwide,*t]));values=values[np.r_[True,np.diff(values)>1e-11]]
    if (len(values)-1)%2:
        at=0 if wide_first else len(values)-2
        values=np.sort(np.r_[values,.5*(values[at]+values[at+1])])
    outer.initial_parameters=values
    inner=layout.boundaries[0]
    values=stations[(stations>=port.cut_station-1e-9)&(stations<=distance+1e-9)]
    inner.initial_parameters=np.sort((values-inner.start)/(inner.finish-inner.start))


def from_failure(context,junction,approaches,failure):
    from .intersections import CornerConstructionError,Corner,BoundaryCurve,JunctionLayout,Mouth,half_intervals
    from .profiles import memberships
    from .corner_radii import CornerRadiusPolicy
    if not isinstance(failure,CornerConstructionError) or len(approaches)!=3:return None
    if any(junction.corner_radius(End(a.source.id,a.end))[1]!='STYLE' for a in approaches):return None
    bundled=memberships(context.network)
    if any(a.profile.half or a.source.id in bundled or getattr(a.path,'absorbed_sites',())
           or getattr(a.path,'deferred_profile_reserve',None) is not None
           or getattr(a.path,'minimum_mouth',0.)>1e-9 for a in approaches):return None
    scoped=[a for a in failure.approaches if len(a.path.pieces)==1
        and 1-a.end in a.source.scope_ends and not a.source.selection_ends
        and len(a.source.points)==2 and abs(a.offset)<1e-9]
    if len(scoped)!=1:return None
    spur=scoped[0];wide=next(a for a in failure.approaches if a is not spur)
    third=next(a for a in approaches if a not in failure.approaches)
    if (spur.width>wide.width or wide.length<3.*spur.length
            or wide.direction@spur.direction<math.cos(math.radians(30.))
            or wide.direction@third.direction>-math.cos(math.radians(15.))
            or any(np.linalg.norm(a.origin-spur.origin)>1e-8 for a in approaches)
            or max(a.profile.datum for a in approaches)-min(a.profile.datum for a in approaches)>1e-9):return None
    first,second=failure.approaches
    wide_first=first is wide
    wide_side=-1. if wide_first else 1.
    outer_side=-1. if wide_first else 1.
    center=spur.frame(spur.length)[0]
    def cut(s):return float((wide.point(s,wide.offset+wide_side*.5*wide.width)-center)@spur.direction)
    lo=0.;hi=min(wide.length,spur.length+2.*wide.width)
    if cut(lo)>=0. or cut(hi)<=0.:return None
    for _ in range(56):
        mid=.5*(lo+hi)
        if cut(mid)<0.:lo=mid
        else:hi=mid
    at=.5*(lo+hi)
    point=wide.point(at,wide.offset+wide_side*.5*wide.width)
    scope_normal=spur.frame(spur.length)[2]
    physical_half=float((point-center)@scope_normal)*(-outer_side)
    gap=physical_half-.5*spur.width
    requested,mode=junction.corner_radius(End(first.source.id,first.end))
    minimum=min(CornerRadiusPolicy.contact_candidates(first,second,requested,mode))
    if (gap<0. or gap>2.*minimum+1e-7 or physical_half<.9375*.5*spur.width
            or abs(float((point-center)@spur.direction))>1e-7):return None
    ordered=sorted(approaches,key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
    index=ordered.index(first)
    ordered=ordered[index:]+ordered[:index]
    if ordered[1] is not second:return None
    c1=CornerRadiusPolicy.solve(second,third,*junction.corner_radius(End(second.source.id,second.end)))
    c2=CornerRadiusPolicy.solve(third,first,*junction.corner_radius(End(third.source.id,third.end)))
    source_rim=np.array((spur.point(0.,-.5*spur.width),spur.point(spur.length,-.5*spur.width),
        spur.point(spur.length,.5*spur.width),spur.point(0.,.5*spur.width)))
    third_support=float(np.max((source_rim-third.origin)@third.direction))+.25
    distances=[0.,0.,max(0.,c1.stations[1],c2.stations[0],third_support)+min(.25,.05*third.width)]
    wi=0 if wide_first else 1
    distances[wi]=max(at+max(2.,wide.width),c2.stations[1] if wi==0 else c1.stations[0])+.25
    if distances[wi]>=.4*wide.length or distances[2]>=third.length-1e-6:return None
    settings=context.network.settings;radial=half_intervals(settings)
    port=FiniteScopeSpurPort(spur,wide,at,distances[wi],point,outer_side,gap)
    output=[port if a is spur else a for a in ordered]
    proxy=Corner(output[0],output[1],(0.,0.),None,requested,0.,0.,
        math.acos(float(np.clip(wide.direction@spur.direction,-1.,1.))))
    lead=ScopeCurbLead(proxy,wide,distances[wi] if wide_first else at,
        at if wide_first else distances[wi],wide_side,settings)
    # The other corners continue to use the original source curb contacts.
    d0=spur.length if first is spur else distances[0]
    d1=spur.length if second is spur else distances[1]
    boundaries=(lead,BoundaryCurve(c1,(d1,distances[2]),settings.curve_tolerance),
        BoundaryCurve(c2,(distances[2],d0),settings.curve_tolerance))
    rows=4*int(math.ceil(max(2*radial,*[b.length/min(settings.spacing,
        2.*math.sqrt(max(1e-12,2.*b.corner.radius*settings.curve_tolerance-settings.curve_tolerance**2)))
        for b in boundaries])/4.))
    if 3*(rows+2*radial)>8192:raise ValueError('Open scope pool exceeds its rim sample budget')
    mouths={}
    for a,d in ((wide,distances[wi]),(third,distances[2])):
        for identity,station in a.path.trims(d).items():mouths[identity]=Mouth(float(station),2*radial)
    alias=port.aliases[0]
    mouths[alias.junction_end.section_id,alias.junction_end.end]=Mouth(alias.scope_station,2*radial)
    distances[1-wi]=0.
    layout=JunctionLayout(tuple(output),(proxy,c1,c2),boundaries,tuple(distances),mouths,rows,radial,
        np.zeros(2),[0.]*3,mesh_strategy='boundary_chart',finite_scope_ports=(port,))
    _align_support_samples(layout,port,wide_first,settings)
    return layout
