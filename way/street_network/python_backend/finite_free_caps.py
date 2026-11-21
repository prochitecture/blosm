"""One finite, rounded free cap with two original incoming profiles.

An authored short, wide free carriageway can exceed the reach of a two-curb
through pool. Give its *actual terminal plane* a zero-height curb side, keep
the two genuine incoming profiles, and union their original finite source
strips with certified rounded cap sides. No scope end is reclassified here.
"""
from dataclasses import dataclass,replace
from types import SimpleNamespace
import math
import numpy as np
from .compounds import SampledFootprint,RimSegment,CompoundGroup,SharedClaim,rim_label
from .continuations import WidthContinuation
from .corner_radii import CornerRadiusPolicy
from .finite_containment import contains_polygon
from .finite_scope_spurs import _source_stations
from .model import End,Key
from .regions import Region,RegionUnion


class PreparedCapCurb:
    continuation=None

    def __init__(self,points,corner=None):
        points=np.asarray(points,float)
        self.points=points[np.r_[True,np.linalg.norm(np.diff(points,axis=0),axis=1)>1e-9]]
        self.contact_reference_points=self.points
        self.stations=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(self.points,axis=0),axis=1))]
        self.length=float(self.stations[-1])
        self.corner=corner or SimpleNamespace(center=None)
        if self.length<=1e-8:raise ValueError('Finite free-cap curb has no positive length')

    def __call__(self,t):
        s=np.asarray(t)*self.length
        return np.stack([np.interp(s,self.stations,self.points[:,axis]) for axis in (0,1)],axis=-1)

    def certified_parameters(self,sampler,first,last):
        knots=sorted({first,last,*(float(s/self.length) for s in self.stations
            if first*self.length<s<last*self.length)})
        result=[]
        for a,b in zip(knots,knots[1:]):result.extend(sampler.regular(a,b,sampler.spacing/self.length)[:-1])
        return np.r_[result,last]


@dataclass
class OwnedFiniteFreeCap:
    pad: object
    approaches: tuple
    distances: tuple
    radius: float
    side_curves: tuple
    outer_corner: object
    cap_points: np.ndarray
    requested_radius: float=0.
    radius_attempts: tuple=()
    source_union_cap: bool=True

    def report(self):
        return dict(section=self.pad.source.id,end=self.pad.end,
            consumed_interval=[0.,self.pad.length],policy='FINITE_FREE_CAP_SOURCE_UNION',
            cap_radius=self.radius,cap_points=self.cap_points.tolist(),
            requested_radius=self.requested_radius,radius_attempts=list(self.radius_attempts),
            cap_side_curvature_bounds=[c.curvature_bound for c in self.side_curves],
            incoming_sections=[a.source.id for a in self.approaches],
            incoming_stations=list(self.distances),independent_terminal_profile=False,
            terminal_boundary='ZERO_HEIGHT_CURB_ON_ORIGINAL_FREE_END_PLANE',
            endpoint_extrapolation=False)


def _sample_curve(curve,settings):
    count=curve.intervals(min(settings.spacing,.3),min(.00025,.05*settings.curve_tolerance))
    count+=count%2
    if count>2048:raise ValueError('Finite free cap exceeds its bounded curve sample budget')
    return curve(np.linspace(0.,1.,count+1))


def _port(owner,approach,station,radial,stations):
    right=approach.point(stations,approach.offset+.5*approach.width)
    left=approach.point(stations,approach.offset-.5*approach.width)[::-1]
    u=np.linspace(1.,0.,2*radial+1)
    mouth=approach.point(station,approach.offset+approach.width*(u-.5))
    curves=(PreparedCapCurb(right),PreparedCapCurb(left))
    for curve, native_stations in zip(curves, (stations, stations[::-1])):
        if len(native_stations) != len(curve.points):
            raise ValueError('Finite free-cap bank lost an original station')
        curve.approach_stations = native_stations
    xy=np.vstack((right[:-1],mouth[:-1],left))
    labels=[]
    labels.extend(RimSegment('curb',0,float(a/curves[0].length),float(b/curves[0].length))
        for a,b in zip(curves[0].stations,curves[0].stations[1:]))
    labels.extend(RimSegment('mouth',0,float(a),float(b)) for a,b in zip(u,u[1:]))
    labels.extend(RimSegment('curb',1,float(a/curves[1].length),float(b/curves[1].length))
        for a,b in zip(curves[1].stations,curves[1].stations[1:]))
    labels.append(RimSegment('internal',0,0.,1.))
    layout=SimpleNamespace(approaches=(approach,),outward_stations=(station,),boundaries=curves,
        radial=radial,rows=len(xy),mouths={},corners=())
    return SampledFootprint(Region(owner,xy),tuple(labels),layout)


def _group(context,owner,layout):
    claim=next(c for c in layout.finite_free_pads if getattr(c,'source_union_cap',False))
    a,b=claim.approaches;da,db=claim.distances;ca,cb=claim.side_curves
    settings=context.network.settings;radial=layout.radial
    corner=claim.outer_corner
    stations=[]
    for approach,distance,contact in ((a,da,corner.stations[1]),(b,db,corner.stations[0])):
        values=_source_stations(approach,distance,min(.3,settings.spacing),min(.00025,.05*settings.curve_tolerance))
        values=np.array(sorted({*values,float(contact)}))
        if np.min(values)<-1e-9 or np.max(values)>distance+1e-9:
            raise ValueError('Finite free cap needs original positive external curb contacts')
        stations.append(values)
    pa=_port(Key('free_cap_port',2*a.source.id+a.end),a,da,radial,stations[0])
    pb=_port(Key('free_cap_port',2*b.source.id+b.end),b,db,radial,stations[1])
    cap=claim.cap_points
    first=_sample_curve(ca,settings);last=_sample_curve(cb,settings)
    capline=np.linspace(cap[0],cap[1],2*radial+1)
    capcurve=PreparedCapCurb(np.vstack((first[:-1],capline[:-1],last)))
    sb=stations[1][stations[1]>=corner.stations[0]-1e-9][::-1]
    sa=stations[0][stations[0]>=corner.stations[1]-1e-9]
    start=b.point(sb,b.offset-.5*b.width)
    finish=a.point(sa,a.offset+.5*a.width)
    tol=min(.00025,.05*settings.curve_tolerance)
    count=max(2,int(math.ceil(abs(corner.sweep)*corner.radius/min(settings.spacing,math.sqrt(8.*tol*corner.radius)))))
    count+=count%2
    theta=np.linspace(corner.start_angle,corner.start_angle+corner.sweep,count+1)
    arc=corner.center+corner.radius*np.stack((np.cos(theta),np.sin(theta)),axis=-1)
    outer=PreparedCapCurb(np.vstack((start[:-1],arc[:-1],finish)),corner)
    u=np.linspace(1.,0.,2*radial+1)
    ma=a.point(da,a.offset+a.width*(u-.5));mb=b.point(db,b.offset+b.width*(u-.5))
    xy=np.vstack((ma[:-1],capcurve.points[:-1],mb[:-1],outer.points[:-1]))
    labels=[RimSegment('internal',0,float(x),float(y)) for x,y in zip(u,u[1:])]
    labels.extend(RimSegment('curb',0,float(x/capcurve.length),float(y/capcurve.length))
        for x,y in zip(capcurve.stations,capcurve.stations[1:]))
    labels.extend(RimSegment('internal',1,float(x),float(y)) for x,y in zip(u,u[1:]))
    labels.extend(RimSegment('curb',1,float(x/outer.length),float(y/outer.length))
        for x,y in zip(outer.stations,outer.stations[1:]))
    derived=SimpleNamespace(approaches=(a,b),outward_stations=(da,db),boundaries=(capcurve,outer),
        radial=radial,rows=len(xy),mouths=layout.mouths,corners=(claim.outer_corner,),
        finite_free_pads=(claim,))
    # Original port owners win coincident far cross-sections. The prepared
    # pool's own closing cross-sections remain explicitly internal.
    pool=SampledFootprint(Region(Key('z_free_cap_pool',owner.id),xy),tuple(labels),derived)
    sources=(pa,pb,pool)
    if sum(len(s.region.xy) for s in sources)>8192:raise ValueError('Finite free-cap source union exceeds its rim budget')
    loops=RegionUnion().build(s.region for s in sources)
    if len(loops)!=1 or loops[0].hole:raise ValueError('Finite free cap needs one simple supported exterior')
    inventory={s.region.owner:s for s in sources}
    if any(rim_label(span,inventory).kind=='internal' for span in loops[0].spans):
        raise ValueError('Finite free cap exposes an internal source cross-section')
    pad=claim.pad
    rectangle=np.array((pad.point(0.,-.5*pad.width),pad.point(pad.length,-.5*pad.width),
        pad.point(pad.length,.5*pad.width),pad.point(0.,.5*pad.width)))
    if not contains_polygon(rectangle,loops[0].xy,tolerance=1e-7):
        raise ValueError('Finite free cap lost its complete original wide pad footprint')
    footprint,physical=_physical_exterior(owner,layout,sources,loops[0])
    loops=RegionUnion().build((footprint.region,))
    group=CompoundGroup((owner,),(SharedClaim(pad.source.id,owner,owner,pad.length,0.),),
        (footprint,),loops,physical)
    from .compound_ports import CompoundPortPlanner
    if CompoundPortPlanner.incomplete_ports(group):raise ValueError('Finite free cap clips an original incoming profile')
    return group


def _physical_exterior(owner,layout,sources,loop):
    """Consolidate the union into two complete physical ports and curb chains.

    Both original full profiles are checked before assigning the derived
    physical owner. Changing a source of an adjacent retained curb is not a
    clipped profile, but separate provisional owners cannot remain as the
    final physical-port boundary description.
    """
    inventory={s.region.owner:s for s in sources}
    labels=[rim_label(span,inventory) for span in loop.spans]
    port_owners={sources[0].region.owner:0,sources[1].region.owner:1}
    tokens=[port_owners[span.owner] if label.kind=='mouth' and span.owner in port_owners else None
        for span,label in zip(loop.spans,labels)]
    starts={}
    for i,token in enumerate(tokens):
        if token is not None and token!=tokens[i-1]:
            if token in starts:raise ValueError('Finite free cap splits an original incoming profile')
            starts[token]=i
    if set(starts)!={0,1}:raise ValueError('Finite free cap hides an original incoming profile')
    start=starts[0]
    xy=np.roll(loop.xy,-start,axis=0)
    labels=labels[start:]+labels[:start];tokens=tokens[start:]+tokens[:start]
    second=next(i for i,t in enumerate(tokens) if t==1)
    ends=[]
    for token,begin in ((0,0),(1,second)):
        index=begin;value=1.
        while index<len(tokens) and tokens[index]==token:
            label=labels[index]
            if abs(label.first-value)>1e-8 or label.last>=label.first:
                raise ValueError('Finite free cap breaks an incoming source parameter partition')
            value=label.last;index+=1
        if abs(value)>1e-8:raise ValueError('Finite free cap clips an original incoming profile')
        ends.append(index)
    if ends[0]>=second or ends[1]>=len(xy):raise ValueError('Finite free cap needs two positive physical curb chains')
    curves=(PreparedCapCurb(xy[ends[0]:second+1]),
            PreparedCapCurb(np.vstack((xy[ends[1]:],xy[0]))))
    # Keep the source of each exposed chord on the free-cap side. The other
    # boundary contains the through corner's circular fillet and its contacts.
    ordered_spans = loop.spans[start:] + loop.spans[:start]
    ancestry = []
    for span, label in zip(ordered_spans[ends[0]:second], labels[ends[0]:second]):
        if label.kind != 'curb':
            raise ValueError('The free-cap side contains an internal edge')
        if span.owner in port_owners and label.index in (0, 1):
            bank = inventory[span.owner].layout.boundaries[label.index]
            sa, sb = np.interp(np.array((label.first, label.last))*bank.length,
                               bank.stations, bank.approach_stations)
            ancestry.append((port_owners[span.owner], 1 if label.index == 0 else -1, sa, sb))
        elif span.owner == sources[2].region.owner and label.index == 0:
            ancestry.append(None)  # Generated side or the transverse free curb.
        else:
            break  # A different arrangement needs its own source description.
    if len(ancestry) == len(curves[0].points)-1:
        curves[0].free_cap_banks = tuple(ancestry)
    result=[]
    for i,label in enumerate(labels):
        if tokens[i] is not None:
            result.append(RimSegment('mouth',tokens[i],label.first,label.last))
        else:
            side=0 if i<second else 1
            at=i-ends[side]
            curve=curves[side]
            result.append(RimSegment('curb',side,float(curve.stations[at]/curve.length),
                float(curve.stations[at+1]/curve.length)))
    physical=replace(layout,boundaries=curves,rows=max(len(c.points)-1 for c in curves))
    return SampledFootprint(Region(owner,xy),tuple(result),physical),physical


def groups(context,layouts):
    return tuple(_group(context,key,layout) for key,layout in layouts.items()
        if any(getattr(c,'source_union_cap',False) for c in getattr(layout,'finite_free_pads',())))


def from_failure(context,junction,approaches,failure):
    from .intersections import CornerConstructionError,ApproachExtentError,JunctionLayout,Mouth,half_intervals
    from .profiles import memberships
    if not isinstance(failure,(CornerConstructionError,ApproachExtentError)) or len(approaches)!=3:return None
    if any(junction.corner_radius(End(a.source.id,a.end))[1]!='STYLE' for a in approaches):return None
    bundled=memberships(context.network)
    if any(a.source.id in bundled or a.profile.half or abs(a.offset)>1e-8
        or getattr(a.path,'absorbed_sites',()) or getattr(a.path,'minimum_mouth',0.)>1e-9
        or getattr(a.path,'deferred_profile_reserve',None) is not None for a in approaches):return None
    failed={a.source.id for a in failure.approaches}
    pads=[a for a in approaches if a.source.id in failed and getattr(a.path,'terminal_free',False)
        and len(a.path.pieces)==1 and len(a.source.points)==2
        and not a.source.scope_ends and not a.source.selection_ends
        and a.length<=.5*a.width]
    if len(pads)!=1:return None
    pad=pads[0]
    ordered=sorted(approaches,key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
    at=ordered.index(pad)
    a,b=ordered[at-1],ordered[(at+1)%3]
    if (max(x.profile.datum for x in approaches)-min(x.profile.datum for x in approaches)>1e-9
            or any(np.linalg.norm(x.origin-pad.origin)>1e-8 for x in approaches)):return None
    try:corner=CornerRadiusPolicy.solve(b,a,*junction.corner_radius(End(b.source.id,b.end)))
    except ValueError:return None
    if corner.center is None or min(corner.stations)<0.:return None
    settings=context.network.settings;radial=half_intervals(settings)
    cap,normal=pad.frame(pad.length)[::2]
    rectangle=np.array((pad.point(0.,-.5*pad.width),pad.point(pad.length,-.5*pad.width),
        pad.point(pad.length,.5*pad.width),pad.point(0.,.5*pad.width)))
    base_a=max(0.,corner.stations[1],float(np.max((rectangle-a.origin)@a.direction)))+.5
    base_b=max(0.,corner.stations[0],float(np.max((rectangle-b.origin)@b.direction)))+.5
    requested,mode=junction.corner_radius(End(pad.source.id,pad.end))
    radii=tuple(CornerRadiusPolicy.contact_candidates(pad,b,requested,mode))
    for radius_index,radius in enumerate(radii):
        for factor in (1.,1.5,2.,3.,4.,6.):
            da,db=base_a*factor,base_b*factor
            if da>=.4*a.length or db>=.4*b.length:continue
            ends=np.array((cap+(.5*pad.width+radius)*normal,cap-(.5*pad.width+radius)*normal))
            try:
                pa,va=a.jet(da,a.offset-.5*a.width);pb,vb=b.jet(db,b.offset+.5*b.width)
                curves=(WidthContinuation(pa,ends[0],-va,-normal,radius,refine_curvature=True),
                    WidthContinuation(ends[1],pb,-normal,vb,radius,refine_curvature=True))
                if np.max((curves[0].controls-pad.origin)@pad.direction)>pad.length+1e-8:
                    raise ValueError('Finite free-cap side crosses the original terminal plane')
                claim=OwnedFiniteFreeCap(pad,(a,b),(da,db),radius,curves,corner,ends,
                    requested_radius=requested,radius_attempts=tuple(radii[:radius_index+1]))
                mouths={}
                for x,d in ((a,da),(b,db)):
                    for identity,station in x.path.trims(d).items():mouths[identity]=Mouth(float(station),2*radial)
                mouths[pad.source.id,pad.end]=Mouth(pad.length if pad.end==0 else 0.,2*radial)
                layout=JunctionLayout((a,b),(corner,),(),(da,db),mouths,8,radial,
                    pad.origin.copy(),[0.,0.],mesh_strategy='boundary_chart',finite_free_pads=(claim,))
                _group(context,junction.key,layout)
                return layout
            except ValueError:continue
    return None
