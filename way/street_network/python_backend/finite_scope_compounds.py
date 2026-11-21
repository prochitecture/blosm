"""Compose an open finite scope pool with its adjacent source junction.

The adjacent narrow ground road may remain beside the scope spur until the
retained ground ends. Its facing curb then joins the *same original scope
plane*, rather than an impossible fillet beyond that plane. Only original
continuing road strips may cross the plane. Every consumed source incidence,
full strip and far incoming profile remains explicit and is certified here.
"""
from dataclasses import replace
from types import SimpleNamespace
import math
import numpy as np
from .compounds import CompoundGroup,SharedClaim,SampledFootprint,RimSegment,rim_label
from .corner_radii import CornerRadiusPolicy
from .finite_containment import contains_polygon
from .finite_free_caps import PreparedCapCurb
from .finite_scope_spurs import FiniteScopeSpurPort,_source_stations
from .model import End,Key
from .regions import Region,RegionUnion,signed_area


class ResampledScopeCurb(PreparedCapCurb):
    """Bounded projection onto a derived union polyline, with its provenance.

    The original analytical primitives remain in the physical layout and
    source curve references. This class records additional *polyline* error;
    it never claims exact equivalence to an analytical circle or offset curve.
    """
    def __init__(self,points,source_spans,source_curves,tolerance):
        original=np.asarray(points,dtype=float).copy()
        keep={0,len(original)-1};pending=[(0,len(original)-1)]
        while pending:
            first,last=pending.pop()
            if last-first<2:continue
            edge=original[last]-original[first]
            fraction=np.clip((original[first+1:last]-original[first])@edge/(edge@edge),0.,1.)
            error=np.linalg.norm(original[first+1:last]-original[first]-fraction[:,None]*edge,axis=1)
            at=int(np.argmax(error))
            if error[at]>tolerance:
                index=first+1+at;keep.add(index);pending.extend(((first,index),(index,last)))
        indices=sorted(keep);maximum=0.
        for first,last in zip(indices,indices[1:]):
            edge=original[last]-original[first]
            fraction=np.clip((original[first:last+1]-original[first])@edge/(edge@edge),0.,1.)
            maximum=max(maximum,float(np.linalg.norm(original[first:last+1]-original[first]-fraction[:,None]*edge,axis=1).max()))
        if maximum>tolerance+1e-12:raise ValueError('Composed scope curb exceeds its original polyline projection allowance')
        super().__init__(original[indices])
        self.original_polyline=original;self.original_polyline.flags.writeable=False
        self.source_spans,self.source_curves=tuple(source_spans),tuple(source_curves)
        self.maximum_projection_error,self.projection_allowance=maximum,tolerance
        parameters=self.stations/self.length
        if (len(parameters)-1)%2:
            at=int(np.argmax(np.diff(self.stations)))
            parameters=np.sort(np.r_[parameters,.5*(parameters[at]+parameters[at+1])])
        self.initial_parameters=parameters

    def report(self):
        return dict(reference='sampled_original_ground_union_polyline',original_points=len(self.original_polyline),
            retained_points=len(self.points),maximum_projection_error=self.maximum_projection_error,
            additional_error_allowance=self.projection_allowance,profile_and_scope_endpoints_exact=True,
            original_analytical_primitives_retained=True,
            source_curves=[dict(owner=owner.text(),index=index,curve_type=type(curve).__name__)
                for owner,index,curve in self.source_curves])


class ScopePlaneCuts:
    """Two finite computational cuts separate the continuing source tails."""
    def cuts(self,group,context):
        from .mesh import cross
        from .partition import BoundaryRef,JunctionCut
        port,=group.physical_layout.finite_scope_ports
        if not isinstance(port,ExpandedScopePort) or len(group.loops)!=1:
            raise ValueError('Scope-plane chart needs one composed finite exterior')
        loop=group.loops[0];points=loop.xy;delta=np.roll(points,-1,axis=0)-points
        inventory={p.region.owner:p for p in group.footprints}
        cuts=[]
        for start,direction in ((port.points[0],-port.normal),(port.points[1],port.normal)):
            vertex=int(np.argmin(np.linalg.norm(points-start,axis=1)))
            if np.linalg.norm(points[vertex]-start)>1e-8:
                raise ValueError('Scope-plane chart lost an original cut endpoint')
            denominator=cross(np.broadcast_to(direction,delta.shape),delta)
            offset=points-start;valid=abs(denominator)>1e-12
            t=np.divide(cross(offset,delta),denominator,out=np.full(len(points),np.inf),where=valid)
            u=np.divide(cross(offset,np.broadcast_to(direction,delta.shape)),denominator,
                out=np.full(len(points),np.inf),where=valid)
            hits=np.flatnonzero((t>1e-7)&(u>=-1e-10)&(u<=1.+1e-10))
            if not len(hits):raise ValueError('Scope-plane chart has no finite opposite curb')
            edge=int(hits[np.argmin(t[hits])]);fraction=float(np.clip(u[edge],0.,1.))
            if rim_label(loop.spans[edge],inventory).kind!='curb':
                raise ValueError('Scope-plane chart cannot split an original incoming profile')
            cuts.append(JunctionCut(None,group.members[0],group.members[0],
                BoundaryRef(0,vertex,0.),BoundaryRef(0,edge,fraction)))
        return tuple(cuts)


class ScopeCutCrown:
    """Two reconciled monotone halves with stationary crest and curb jets."""
    half=0

    def __init__(self,rise,datum,peak_u):
        if not 0.<peak_u<1.:raise ValueError('Expanded scope crest lies outside its finite cut')
        self.rise,self.datum,self.peak_u=rise,datum,peak_u
        self.peak=2.*peak_u-1.

    def _half(self,q):
        u=.5*(np.asarray(q)+1.)
        left=u<=self.peak_u
        span=np.where(left,self.peak_u,1.-self.peak_u)
        t=np.where(left,(self.peak_u-u)/span,(u-self.peak_u)/span)
        return left,span,np.clip(t,0.,1.)

    def value(self,q):
        _,_,t=self._half(q)
        return self.rise*(1.-t**3*(10.+t*(-15.+6.*t)))

    def derivative(self,q):
        left,span,t=self._half(q)
        return np.where(left,1.,-1.)*.5*self.rise*30.*t*t*(1.-t)**2/span

    def coordinate(self,side,v):return side+(self.peak-side)*np.asarray(v)
    def intervals(self,side,base):return max(2,2*int(round(abs(self.peak-side)*base/2.)))


class ExpandedScopePort(FiniteScopeSpurPort):
    def __init__(self,old,other,other_cut,other_station,points,members):
        self.__dict__.update(old.__dict__)
        self.other,self.other_cut,self.other_station=other,other_cut,other_station
        self.points=np.asarray(points)
        self.width=float(np.linalg.norm(self.points[1]-self.points[0]))
        self.normal=(self.points[1]-self.points[0])/self.width
        self.origin=.5*(self.points[0]+self.points[1])
        center=self.spur.frame(self.spur.length)[0]
        peak=float((center-self.points[0])@self.normal/self.width)
        self.profile=ScopeCutCrown(self.spur.profile.rise,self.spur.profile.datum,peak)
        self.members=tuple(members)
        self.closing_gaps=(peak*self.width-.5*self.spur.width,
            (1.-peak)*self.width-.5*self.spur.width)
        self.closing_gap=max(self.closing_gaps)

    def report(self):
        result=super().report()
        gaps=getattr(self,'closing_gaps',(self.profile.peak_u*self.width-.5*self.spur.width,
            (1.-self.profile.peak_u)*self.width-.5*self.spur.width))
        result.update(mode='COMPOSED_OPEN_GROUND_SCOPE',profile='reconciled_complete_crown',
            reconciled_inner_half=False,reconciled_both_halves=True,
            source_containment='finite_ground_side_plus_two_original_continuing_strips',
            member_junctions=[key.text() for key in self.members],
            second_continuing_section=self.other.source.id,
            second_curb_cut_station=self.other_cut,
            second_complete_profile_station=self.other_station,
            closing_gaps=list(gaps),original_scope_width=self.spur.width,
            original_scope_outline=[self.spur.point(self.spur.length,side*.5*self.spur.width).tolist()
                for side in (-1.,1.)])
        return result


def _cut(approach,side,center,direction):
    def value(s):return float((approach.point(s,approach.offset+side*.5*approach.width)-center)@direction)
    # The complete finite route must cross the plane only once. The ordinary
    # source chord tolerance also bounds possible offset-curve crossings.
    samples=_source_stations(approach,approach.length,.3,.00025)
    values=np.asarray([value(s) for s in samples])
    changes=np.flatnonzero((values[:-1]<0.)&(values[1:]>=0.))
    if values[0]>=0. or values[-1]<=0. or len(changes)!=1 or np.any(
            (values[:-1]>1e-7)&(values[1:]<-1e-7)):
        raise ValueError('Composed scope needs one original finite curb crossing')
    at=int(changes[0]);lo,hi=samples[at],samples[at+1]
    for _ in range(56):
        mid=.5*(lo+hi)
        if value(mid)<0.:lo=mid
        else:hi=mid
    station=.5*(lo+hi)
    return station,approach.point(station,approach.offset+side*.5*approach.width)


def _source_prefix(owner,approach,distance,radial,stations,*,internal=False):
    right=approach.point(stations,approach.offset+.5*approach.width)
    left=approach.point(stations,approach.offset-.5*approach.width)[::-1]
    u=np.linspace(1.,0.,2*radial+1)
    mouth=approach.point(distance,approach.offset+approach.width*(u-.5))
    curves=(PreparedCapCurb(right),PreparedCapCurb(left))
    xy=np.vstack((right[:-1],mouth[:-1],left))
    labels=[RimSegment('curb',0,float(a/curves[0].length),float(b/curves[0].length))
        for a,b in zip(curves[0].stations,curves[0].stations[1:])]
    labels.extend(RimSegment('internal' if internal else 'mouth',0,float(a),float(b))
        for a,b in zip(u,u[1:]))
    labels.extend(RimSegment('curb',1,float(a/curves[1].length),float(b/curves[1].length))
        for a,b in zip(curves[1].stations,curves[1].stations[1:]))
    labels.append(RimSegment('internal',0,0.,1.))
    layout=SimpleNamespace(approaches=(approach,),outward_stations=(distance,),boundaries=curves,
        radial=radial,rows=len(xy),mouths={},corners=())
    return SampledFootprint(Region(owner,xy),tuple(labels),layout)


def _physical_exterior(owner,layout,sources,loop,ports):
    """Retain every whole source profile before assigning one physical owner."""
    inventory={s.region.owner:s for s in sources}
    labels=[rim_label(span,inventory) for span in loop.spans]
    if any(label.kind not in ('curb','mouth') for label in labels):
        raise ValueError('Composed finite scope exposes an internal source cross-section')
    tokens=[ports.get((span.owner,label.index)) if label.kind=='mouth' else None
        for span,label in zip(loop.spans,labels)]
    if any(label.kind=='mouth' and token is None for label,token in zip(labels,tokens)):
        raise ValueError('Composed finite scope has an unmapped exposed mouth')
    starts={}
    for i,token in enumerate(tokens):
        if token is not None and token!=tokens[i-1]:
            if token in starts:raise ValueError('Composed finite scope splits an incoming profile')
            starts[token]=i
    if set(starts)!=set(range(len(layout.approaches))):
        raise ValueError('Composed finite scope hides an incoming profile')
    start=starts[0];xy=np.roll(loop.xy,-start,axis=0)
    spans=list(loop.spans[start:])+list(loop.spans[:start])
    labels=labels[start:]+labels[:start];tokens=tokens[start:]+tokens[:start]
    starts=sorted((i,token) for i,token in enumerate(tokens) if token is not None and token!=tokens[i-1])
    order=[token for _,token in starts];curves=[]
    for side,(begin,token) in enumerate(starts):
        finish=starts[side+1][0] if side+1<len(starts) else len(tokens)
        index=begin;value=1.
        while index<finish and tokens[index]==token:
            label=labels[index]
            if abs(label.first-value)>1e-8 or label.last>=label.first:
                raise ValueError('Composed finite scope breaks an original profile partition')
            value=label.last;index+=1
        if abs(value)>1e-8 or index>=finish:
            raise ValueError('Composed finite scope clips a profile or loses a physical curb')
        points=np.vstack((xy[index:finish],xy[finish%len(xy)]))
        provenance={}
        for span,label in zip(spans[index:finish],labels[index:finish]):
            provenance[span.owner,label.index]=inventory[span.owner].layout.boundaries[label.index]
        curve=ResampledScopeCurb(points,spans[index:finish],
            [(key[0],key[1],value) for key,value in sorted(provenance.items())],1e-6)
        curves.append(curve)
    physical=replace(layout,approaches=tuple(layout.approaches[i] for i in order),
        outward_stations=tuple(layout.outward_stations[i] for i in order),
        guide_weights=[0.]*len(order),boundaries=tuple(curves),
        rows=max(len(c.points)-1 for c in curves))
    physical.finite_scope_curb_sampling=tuple(curve.report() for curve in curves)
    from .compounds import sample_footprint
    return sample_footprint(owner,physical),physical


def _compose(context,sites,layouts,group,neighbor):
    from .compound_preparation import DeferredCompoundFootprint
    from .corner_patch import corner_patch
    from .intersections import JunctionLayout,Mouth
    from .profiles import memberships
    from .source_strip_containment import certify
    old_layout=group.physical_layout
    ports=old_layout.finite_scope_ports
    if len(ports)!=1 or not isinstance(ports[0],FiniteScopeSpurPort) or isinstance(ports[0],ExpandedScopePort):
        raise ValueError('Deferred finite composition needs one original open scope spur')
    old=ports[0];spur,wide=old.spur,old.wide
    source=group.members[0];deferred=layouts[neighbor]
    if not isinstance(deferred,DeferredCompoundFootprint) or len(deferred.approaches)!=3 or len(old_layout.approaches)!=3:
        raise ValueError('Deferred finite composition needs two original three-way source sites')
    internal=[a for a in old_layout.approaches if a is not old and a is not wide]
    if len(internal)!=1:raise ValueError('Deferred finite composition has ambiguous internal source')
    opposite=internal[0]
    matching=[a for a in deferred.approaches if a.source.id==opposite.source.id and a.end==1-opposite.end]
    if (len(matching)!=1 or len(opposite.path.pieces)!=1 or len(matching[0].path.pieces)!=1
            or source not in deferred.neighbors):
        raise ValueError('Deferred finite composition lacks both original connecting incidences')
    inner=matching[0]
    # Walk the original incidence cycle, rather than sorting displaced arms.
    from .compound_preparation import DeferredCompoundPlanner
    planner=DeferredCompoundPlanner(context,tuple(sites.values()),layouts)
    members=tuple(sorted((source,neighbor)))
    internal_ends={End(opposite.source.id,opposite.end),End(inner.source.id,inner.end)}
    order=planner.exterior_order(members,internal_ends)
    by_end={End(a.source.id,a.end):a for a in (*deferred.approaches,*old_layout.approaches)
        if End(a.source.id,a.end) not in internal_ends}
    scope_end=End(old.source.id,old.end);at=order.index(scope_end)
    before,after=by_end[order[at-1]],by_end[order[(at+1)%len(order)]]
    if before is wide:other=after;wide_side=-1.;other_side=1.
    elif after is wide:other=before;wide_side=1.;other_side=-1.
    else:raise ValueError('Composed finite scope lost its original continuing incidence')
    east=next(a for a in deferred.approaches if a is not inner and a is not other)
    bundled=memberships(context.network)
    for key in members:
        if any(mode!='STYLE' for _,mode in (sites[key].corner_radius(end) for end in sites[key].ends)):
            raise ValueError('Composed finite scope requires authored style radii')
    if any(a.profile.half or a.source.id in bundled or abs(a.offset)>1e-8
            or getattr(a.path,'absorbed_sites',()) or getattr(a.path,'minimum_mouth',0.)>1e-9
            or getattr(a.path,'deferred_profile_reserve',None) is not None
            for a in (spur,wide,other,east,inner,opposite)):
        raise ValueError('Composed finite scope needs ordinary original ground source strips')
    if (not spur.source.scope_ends or spur.source.selection_ends
            or 1-spur.end not in spur.source.scope_ends or len(spur.source.points)!=2
            or max(a.profile.datum for a in (spur,wide,other,east,inner))-
               min(a.profile.datum for a in (spur,wide,other,east,inner))>1e-8):
        raise ValueError('Composed finite scope needs the original straight ground scope incidence')
    center=spur.frame(spur.length)[0]
    cw,pw=_cut(wide,wide_side,center,spur.direction)
    co,po=_cut(other,other_side,center,spur.direction)
    # A normal physical gore must be tried first. This policy is only for a
    # finite original scope which cannot support any permitted contact.
    first,second=(spur,other) if before is wide else (other,spur)
    try:
        site=source if first is spur else neighbor
        CornerRadiusPolicy.solve(first,second,*sites[site].corner_radius(End(first.source.id,first.end)),displaced=True)
    except ValueError:pass
    else:raise ValueError('Composed open scope is unnecessary where a finite physical fillet fits')
    corners=[]
    local_by_site={source:{End(a.source.id,a.end):a for a in (wide,spur,opposite)},
                   neighbor:{End(a.source.id,a.end):a for a in deferred.approaches}}
    for key in members:
        ordered=sorted(local_by_site[key].values(),key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
        for a,b in zip(ordered,ordered[1:]+ordered[:1]):
            if a is spur or b is spur or (key==neighbor and ((a is inner and b is other) or (a is other and b is inner))):continue
            corners.append(CornerRadiusPolicy.solve(a,b,*sites[key].corner_radius(End(a.source.id,a.end))))
    dw=max(old.wide_station,cw+max(2.,wide.width)+.25)
    do=co+max(2.,other.width)+.25
    # The complete original root cross-sections must lie behind the external
    # east mouth. Their finite transverse projection matters at shallow forks.
    roots=np.vstack([a.point(0.,np.array((-.5,.5))*a.width) for a in (inner,other)])
    de=max(0.,float(np.max((roots-east.origin)@east.direction)),
        *[c.stations[i] for c in corners for i,a in enumerate((c.first,c.second)) if a is east])+.25
    distances={id(wide):dw,id(other):do,id(east):de,id(inner):inner.length,id(spur):spur.length}
    if any(d>=a.length-1e-6 for a,d in ((wide,dw),(other,do),(east,de))):
        raise ValueError('Composed finite scope consumes a complete external source route')
    points=(po,pw) if float((pw-po)@spur.frame(spur.length)[2])>0. else (pw,po)
    port=ExpandedScopePort(old,other,co,do,points,members)
    port.wide_station=dw
    radial=old_layout.radial;settings=context.network.settings
    samples={}
    for a in (wide,other,east,inner,spur):
        d=distances[id(a)]
        values=_source_stations(a,d,min(.2,settings.spacing),min(.00025,.05*settings.curve_tolerance))
        knots=[cw] if a is wide else [co] if a is other else []
        knots += [c.stations[i] for c in corners for i,x in enumerate((c.first,c.second)) if x is a and 0.<c.stations[i]<d]
        samples[id(a)]=np.asarray(sorted({*values,*knots}))
    sources=[];port_map={}
    for a in (east,wide,other,inner,spur):
        owner=Key('scope_compound_source',2*a.source.id+a.end)
        f=_source_prefix(owner,a,distances[id(a)],radial,samples[id(a)],internal=a in (inner,spur))
        if a not in (inner,spur):port_map[owner,0]=(east,wide,port,other).index(a)
        sources.append(f)
    for i,c in enumerate(corners):
        if c.center is not None:sources.append(corner_patch(Key('scope_compound_corner',i),c,
            min(.2,settings.spacing),min(.00025,.05*settings.curve_tolerance)))
    # The source centreline closing edges are internal, inside their complete
    # original strips. Only the facing wide curb and original cut plane expose.
    sw=samples[id(wide)];sw=sw[sw<=cw+1e-9]
    so=samples[id(other)];so=so[so<=co+1e-9][::-1]
    wc=wide.point(sw,wide.offset+wide_side*.5*wide.width)
    oc=other.frame(so)[0]
    ic=inner.frame(samples[id(inner)])[0]
    u=port.mouth_parameters(radial,port.width/(2*radial))
    cap=port.point(0.,port.width*(u-.5))
    # With the mirrored source rotation, traverse the prepared cut backwards
    # and reverse the whole region below, retaining decreasing mouth parameters.
    forward=np.linalg.norm(cap[0]-pw)<1e-7
    if not forward:cap=cap[::-1];u=u[::-1]
    xy=np.vstack((wc[:-1],cap,oc,ic[1:]))
    kinds=['curb']*(len(wc)-1)+['mouth']*(len(cap)-1)+['internal']*(len(oc)+len(ic))
    firsts=[0.]*(len(wc)-1)+list(u[:-1])+[0.]*(len(oc)+len(ic))
    lasts=[1.]*(len(wc)-1)+list(u[1:])+[1.]*(len(oc)+len(ic))
    from .regions import signed_area
    if signed_area(xy)<0.:
        xy=xy[::-1]
        kinds=list(np.roll(np.asarray(kinds)[::-1],-1))
        firsts,lasts=list(np.roll(np.asarray(lasts)[::-1],-1)),list(np.roll(np.asarray(firsts)[::-1],-1))
    # One prepared closed curve supplies exact parameters for retained curb
    # segments; internal parameter values have no physical boundary meaning.
    closed=PreparedCapCurb(np.vstack((xy,xy[0])))
    labels=tuple(RimSegment(kind,0,float(firsts[i]) if kind=='mouth' else float(closed.stations[i]/closed.length),
        float(lasts[i]) if kind=='mouth' else float(closed.stations[i+1]/closed.length)) for i,kind in enumerate(kinds))
    cut_owner=Key('a_scope_compound_cut',members[0].id)
    cut=SampledFootprint(Region(cut_owner,xy),labels,SimpleNamespace(approaches=(port,),outward_stations=(0.,),
        boundaries=(closed,),radial=radial,mouths={},corners=()))
    sources.append(cut);port_map[cut_owner,0]=2
    if sum(len(s.region.xy) for s in sources)>8192:
        raise ValueError('Composed finite scope exceeds its bounded source rim budget')
    # Half-plane membership of every polygon vertex proves every source chord
    # stays on the finite ground side. Only the two exact original continuing
    # strips are exempt; no free polygon may extend an elevated source.
    through={2*wide.source.id+wide.end,2*other.source.id+other.end}
    tolerance=min(.00025,.05*settings.curve_tolerance)+1e-7
    for s in sources:
        if s.region.owner.kind=='scope_compound_source' and s.region.owner.id in through:continue
        if np.max((s.region.xy-center)@spur.direction)>tolerance:
            raise ValueError('Composed scope crosses its grade cut outside original ground strips')
    loops=RegionUnion().build(s.region for s in sources)
    if len(loops)!=1 or loops[0].hole:raise ValueError('Composed finite scope needs one simple physical exterior')
    # Complete support from the internal through road can overtake the old
    # external fillet. A newly exposed intersection of its actual curb and
    # the narrow original curb still needs the permitted rounded contact.
    inventory={s.region.owner:s for s in sources}
    other_owner=Key('scope_compound_source',2*other.source.id+other.end)
    inner_owner=Key('scope_compound_source',2*inner.source.id+inner.end)
    switches=[]
    for i,span in enumerate(loops[0].spans):
        previous=loops[0].spans[i-1]
        if {previous.owner,span.owner}=={other_owner,inner_owner}:
            if all(rim_label(s,inventory).kind=='curb' for s in (previous,span)):switches.append(i)
    rounding_reports=[]
    if switches:
        if len(switches)!=1:raise ValueError('Composed scope needs an unambiguous supported curb contact')
        from .source_curb_rounding import intersection_stations,patch as rounding_patch
        hidden=next(c for c in corners if {id(c.first),id(c.second)}=={id(other),id(east)})
        previous=loops[0].spans[switches[0]-1]
        pair=(other,opposite) if previous.owner==other_owner else (opposite,other)
        rounded=CornerRadiusPolicy.solve(*pair,hidden.requested_radius or hidden.radius,hidden.radius_mode,displaced=True)
        if rounded.center is None or min(rounded.stations)<0.:
            raise ValueError('Composed scope needs a finite supported circular curb contact')
        rounded=replace(rounded,fit_policy='COMPOSED_SCOPE_ORIGINAL_CURB_CONTACT')
        contact=intersection_stations(rounded,loops[0].xy[switches[0]])
        oi=pair.index(other);ii=1-oi
        # Original curb chords and the circle share their exact contact knots.
        for approach,stations in ((other,(rounded.stations[oi],contact[oi])),
                (inner,(inner.length-rounded.stations[ii],inner.length-contact[ii]))):
            values=np.asarray(sorted({*samples[id(approach)],*stations}))
            samples[id(approach)]=values
            identity=Key('scope_compound_source',2*approach.source.id+approach.end)
            sources=[_source_prefix(identity,approach,distances[id(approach)],radial,values,
                internal=approach is inner) if f.region.owner==identity else f for f in sources]
        wedge=rounding_patch(Key('scope_compound_rounding',members[0].id),rounded,contact,
            min(.2,settings.spacing),min(.00025,.05*settings.curve_tolerance),
            partitions=tuple(samples[id(other)] if a is other else inner.length-samples[id(inner)][::-1] for a in pair))
        sources.append(wedge)
        rounding_reports.append(dict(policy='ROUND_EXPOSED_ORIGINAL_CURB_INTERSECTION',
            radius=rounded.radius,requested_radius=rounded.requested_radius,radius_mode=rounded.radius_mode,
            radius_attempts=list(rounded.radius_attempts),wedge_area=float(signed_area(wedge.region.xy)),
            actual_curbs=[dict(section=a.source.id,end=a.end,width=a.width,
                outward_contact=float(s),outward_intersection=float(q)) for a,s,q in zip(pair,rounded.stations,contact)],
            closure='ACTUAL_ORIGINAL_CURBS',endpoint_extrapolation=False))
        corners.append(rounded)
        loops=RegionUnion().build(s.region for s in sources)
        if len(loops)!=1 or loops[0].hole:raise ValueError('Rounded composed scope lost its simple physical exterior')
    mouths={}
    for a,d in ((east,de),(wide,dw),(other,do)):
        for identity,station in a.path.trims(d).items():mouths[identity]=Mouth(float(station),2*radial)
    alias=port.aliases[0];mouths[alias.junction_end.section_id,alias.junction_end.end]=Mouth(alias.scope_station,2*radial)
    mouths[inner.source.id,0]=Mouth(inner.length,2*radial);mouths[inner.source.id,1]=Mouth(0.,2*radial)
    for site,other_layout in layouts.items():
        if site in members or isinstance(other_layout,DeferredCompoundFootprint):continue
        for (sid,end),mouth in mouths.items():
            opposite_mouth=other_layout.mouths.get((sid,1-end))
            if opposite_mouth is not None and (mouth.station>opposite_mouth.station+1e-7 if end==0 else mouth.station<opposite_mouth.station-1e-7):
                raise ValueError('Composed finite scope reaches another original junction profile')
    physical=JunctionLayout((east,wide,port,other),tuple(corners),(),(de,dw,0.,do),mouths,8,radial,
        spur.origin.copy(),[0.]*4,mesh_strategy='boundary_chart',finite_scope_ports=(port,))
    footprint,physical=_physical_exterior(members[0],physical,sources,loops[0],port_map)
    physical.finite_scope_rounded_contacts=tuple(rounding_reports)
    loops=RegionUnion().build((footprint.region,))
    claims=[SharedClaim(inner.source.id,neighbor if inner.end==0 else source,
        source if inner.end==0 else neighbor,inner.length,0.)]
    claims.extend(replace(c,first_owner=members[0],last_owner=members[0]) for c in group.claims)
    result=CompoundGroup(members,tuple(claims),(footprint,),loops,physical)
    physical.finite_scope_source_proof=certify(context,result,{c.section_id for c in claims})
    for a in (east,wide,other):
        f=next(s for s in sources if s.region.owner==Key('scope_compound_source',2*a.source.id+a.end))
        if not contains_polygon(f.region.xy,loops[0].xy,tolerance=tolerance):
            raise ValueError('Composed finite scope lost an original continuing source prefix')
    from .compound_ports import CompoundPortPlanner
    if CompoundPortPlanner.incomplete_ports(result):raise ValueError('Composed finite scope clips an original full profile')
    return result


def compose(context,junctions,layouts,merged_groups):
    """Replace a finite group and its explicitly connected deferred neighbor."""
    from .compound_preparation import DeferredCompoundFootprint
    sites={j.key:j for j in junctions};result=[];consumed=set()
    for group in merged_groups:
        neighbors=[key for key,layout in layouts.items() if isinstance(layout,DeferredCompoundFootprint)
            and set(layout.neighbors).intersection(group.members)]
        if not neighbors or not getattr(group.physical_layout,'finite_scope_ports',()):result.append(group);continue
        if len(group.members)!=1 or len(neighbors)!=1 or consumed.intersection(neighbors):
            raise ValueError('Finite scope composition needs one unambiguous adjacent original source junction')
        result.append(_compose(context,sites,layouts,group,neighbors[0]));consumed.update(neighbors)
    return tuple(result)
