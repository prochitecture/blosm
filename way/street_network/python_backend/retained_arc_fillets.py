"""Finite circular junction boundaries with exact original parent ancestry.

This module adds no network entity, source incidence, or height condition.
Every operand records its two actual retained analytical parent arcs.
"""
from dataclasses import dataclass,replace
from copy import copy
from types import SimpleNamespace
import math
import numpy as np
from .compounds import RimSpan,RimSegment,SampledFootprint,rim_label
from .regions import Region,RegionUnion,BoundaryLoop,SourceSpan,signed_area
from .finite_free_caps import PreparedCapCurb
from .paths import right
from .generated_junction_curbs import GeneratedJunctionArea
from .model import Key

SOURCE_ARC_POLICY = "5-original-terminal-chord-contact"

class RetainedArcReconciliationError(ValueError):
    """A finite retained-circle construction failed an original source gate."""

class RetainedArcFitUnavailable(ValueError):
    """A radius cannot fit these two finite visible circular intervals."""


@dataclass(frozen=True)
class CircularBoundaryGeometry:
    center: tuple
    radius: float
    start_angle: float
    sweep: float


def phase(corner,point,reference):
    angle=math.atan2(*(np.asarray(point)-corner.center)[::-1])
    values=[(angle-corner.start_angle+2.*math.pi*k)/corner.sweep for k in range(-2,3)]
    return min(values,key=lambda value:abs(value-reference))


def circular_point(corner,value):
    angle=corner.start_angle+corner.sweep*np.asarray(value)
    return corner.center+corner.radius*np.stack((np.cos(angle),np.sin(angle)),axis=-1)


@dataclass(frozen=True)
class RetainedOriginalArc:
    owner: object
    curve_index: int
    original_curve: object
    original_corner: object
    exposed_circle_parameters: tuple
    consumed_circle_parameters: tuple
    source_spans: tuple
    source_labels: tuple
    original_chord_contact_parameter: float
    original_chord_contact_xy: tuple
    analytical_contact_parameter: float
    analytical_contact_xy: tuple
    analytical_to_original_chord_distance: float

    def report(self):
        c=self.original_corner
        result=dict(owner=self.owner.text(),curve_index=self.curve_index,
            analytical_primitive='FINITE_RETAINED_ORIGINAL_CIRCLE',
            original_center=c.center.tolist(),original_radius=c.radius,
            original_requested_radius=c.requested_radius,original_radius_mode=c.radius_mode,
            original_contact_stations=c.stations,
            original_approaches=[dict(section=a.source.id,end=a.end,width=a.width,half=a.profile.half)
                for a in (c.first,c.second)],
            exposed_circle_parameters=self.exposed_circle_parameters,
            consumed_circle_parameters=self.consumed_circle_parameters,
            analytical_contact_parameter=self.analytical_contact_parameter,
            analytical_contact_xy=self.analytical_contact_xy,
            original_chord_contact_parameter=self.original_chord_contact_parameter,
            original_chord_contact_xy=self.original_chord_contact_xy,
            analytical_to_original_chord_distance=self.analytical_to_original_chord_distance,
            source_span_ancestry=[dict(owner=span.owner.text(),label=dict(kind=label.kind,index=label.index,
                first=label.first,last=label.last)) for span,label in zip(self.source_spans,self.source_labels)])
        declared=getattr(self.original_curve,'original_circle_interval',None)
        if declared is not None and declared.prepared_bounds!=(0.,1.):
            result.update(prepared_circle_bounds=declared.prepared_bounds,
                prepared_contact_parameter=float(declared.source_parameter(self.analytical_contact_parameter)),
                prepared_parameter_authority='EXPLICIT_ORIGINAL_PRODUCER_CIRCLE_KNOTS')
        return result


@dataclass(frozen=True)
class RetainedArcFillet(GeneratedJunctionArea):
    owner: object
    center: tuple
    radius: float
    start_angle: float
    sweep: float
    parent_arcs: tuple
    preferred_reference_radius: float
    attempted_radii: tuple
    explicit_candidate: bool
    original_crossing_xy: tuple
    loop_index: int
    join_index: int
    analytical_closing_intersection_xy: tuple
    minimum_derived_style_radius: float=.5
    source_reference_parameters: tuple=()
    source_chord_tolerance: float=.00025
    exact_retained_parent_contacts=True

    @property
    def length(self):return self.radius*abs(self.sweep)

    @property
    def curvature_bound(self):return 1./self.radius

    @property
    def corner(self):
        # This is a circular boundary geometry record, not an intersection
        # Corner or road Approach. It deliberately has no first/second/End.
        return CircularBoundaryGeometry(self.center,self.radius,self.start_angle,self.sweep)

    @property
    def lengths(self):return np.asarray((0.,0.,self.length,self.length))

    @property
    def leads(self):return ()

    @property
    def continuation(self):return None

    @property
    def contact_reference_points(self):
        if not self.source_reference_parameters:
            raise ValueError('Fillet sampling requires its explicit prepared source-reference knots')
        return self(np.asarray(self.source_reference_parameters))

    def certified_parameters(self,sampler,first,last):
        step=min(sampler.spacing,math.sqrt(8.*min(sampler.tolerance,self.source_chord_tolerance)*self.radius))
        low,high=sorted((first,last))
        regular=sampler.regular(low,high,step/self.length)
        knots=np.asarray(self.source_reference_parameters)
        values=np.unique(np.r_[regular,knots[(knots>low)&(knots<high)]])
        return values if last>=first else values[::-1]

    def __call__(self,t):
        angle=self.start_angle+self.sweep*np.asarray(t)
        return np.asarray(self.center)+self.radius*np.stack((np.cos(angle),np.sin(angle)),axis=-1)

    def tangent(self,t):
        angle=self.start_angle+self.sweep*np.asarray(t)
        return np.sign(self.sweep)*np.stack((-np.sin(angle),np.cos(angle)),axis=-1)

    def report(self):
        return dict(policy=SOURCE_ARC_POLICY,
            primitive='FINITE_ANALYTICAL_CIRCLE_BOUNDARY',owner=self.owner.text(),
            center=self.center,radius=self.radius,curvature_bound=self.curvature_bound,
            start_angle=self.start_angle,sweep=self.sweep,length=self.length,
            preferred_reference_radius=self.preferred_reference_radius,
            attempted_radii=self.attempted_radii,
            minimum_derived_style_radius=self.minimum_derived_style_radius,
            auxiliary_style_continuation=self.explicit_candidate,
            original_crossing_xy=self.original_crossing_xy,
            analytical_closing_intersection_xy=self.analytical_closing_intersection_xy,
            internal_closure='THE_TWO_FINITE_ORIGINAL_ANALYTICAL_ARCS',
            loop_index=self.loop_index,join_index=self.join_index,
            changes_original_radius=False,creates_approach_or_original_end=False,
            height_condition_defined=False,source_reference_parameters=self.source_reference_parameters,
            parent_arcs=[arc.report() for arc in self.parent_arcs])


class OriginalAnalyticalCircleCurve:
    """A sampler view of a validated original pure-circle prepared interval."""
    continuation=None
    leads=()

    def __init__(self,owner,index,prepared,original_corner):
        self.source_owner,self.source_curve_index=owner,index
        self.original_prepared_curve=prepared
        self.corner=original_corner
        self.circle_interval=None
        declared=getattr(prepared,'original_circle_interval',None)
        if declared is not None:
            from .source_circular_intervals import PreparedCircularInterval
            if type(declared) is not PreparedCircularInterval:
                raise RetainedArcReconciliationError('Original circle has an unknown prepared interval authority')
            declared.validate(prepared,original_corner)
            if declared.circle_parameters[-1]<declared.circle_parameters[0]:
                raise RetainedArcFitUnavailable('Retained circle needs its original clockwise prepared orientation')
            if declared.prepared_bounds!=(0.,1.):self.circle_interval=declared
        self.length=(prepared.length if self.circle_interval is not None
                     else original_corner.radius*abs(original_corner.sweep))
        bounds=self.circle_bounds
        self.lengths=np.asarray((0.,bounds[0]*self.length,bounds[1]*self.length,self.length))
        self.contact_reference_points=prepared.points
        self.original_reference_parameters=prepared.stations/prepared.length
        self.parameter_map_maximum_error=float(np.max(np.linalg.norm(
            self(self.original_reference_parameters)-prepared.points,axis=1)))
        self.points=prepared.points
        self.stations=np.asarray(self.original_reference_parameters)*self.length
        if self.parameter_map_maximum_error>1e-8:
            raise ValueError('Original prepared curb is not the declared uniform finite circular interval')

    @property
    def circle_bounds(self):
        return self.circle_interval.prepared_bounds if self.circle_interval is not None else (0.,1.)

    def circle_parameter(self,t):
        return self.circle_interval.circle_parameter(t) if self.circle_interval is not None else np.asarray(t)

    def source_parameter(self,t):
        return self.circle_interval.source_parameter(t) if self.circle_interval is not None else np.asarray(t)

    def __call__(self,t):
        if self.circle_interval is None:return circular_point(self.corner,t)
        at=np.asarray(t);flat=at.reshape(-1)
        low,high=self.circle_bounds;inside=(flat>=low)&(flat<=high)
        result=np.asarray(self.original_prepared_curve(flat)).copy()
        result[inside]=circular_point(self.corner,self.circle_parameter(flat[inside]))
        return result.reshape(at.shape+(2,))

    def certified_parameters(self,sampler,first,last):
        step=min(sampler.spacing,math.sqrt(8.*sampler.tolerance*self.corner.radius))
        if self.circle_interval is not None:
            bounds=(first,*(value for value in self.circle_bounds if first<value<last),last)
            values=[]
            for a,b in zip(bounds,bounds[1:]):
                if self.circle_bounds[0]<=.5*(a+b)<=self.circle_bounds[1]:
                    start,end=map(float,self.circle_parameter(np.asarray((a,b))))
                    phase=sampler.regular(start,end,step/(self.corner.radius*abs(self.corner.sweep)))
                    knots=np.asarray(self.source_parameter(phase));knots[0],knots[-1]=a,b
                else:
                    knots=self.original_prepared_curve.certified_parameters(sampler,a,b)
                values.extend(knots[:-1])
            return np.r_[values,last]
        return sampler.regular(first,last,step/self.length)


def _chord_contact(loop,first,last,corner,value,source_inventory,*,prepared_bounds=(0.,1.)):
    target=circular_point(corner,value);direction=(target-corner.center)/corner.radius
    for index in range(first,last+1):
        label=rim_label(loop.spans[index],source_inventory)
        if label.last<=prepared_bounds[0] or label.first>=prepared_bounds[1]:continue
        a,b=loop.xy[index],loop.xy[(index+1)%len(loop.xy)]
        edge=b-a;denominator=float(np.linalg.det(np.stack((direction,edge))))
        if abs(denominator)<1e-16:continue
        u=float(np.linalg.det(np.stack((a-corner.center,direction)))/denominator)
        if not -1e-10<=u<=1.+1e-10:continue
        u=float(np.clip(u,0.,1.));point=a+u*edge
        if float((point-corner.center)@direction)<=0.:continue
        parameter=label.first+u*(label.last-label.first)
        if not prepared_bounds[0]-1e-12<=parameter<=prepared_bounds[1]+1e-12:continue
        return index,point,float(parameter),target
    raise ValueError('Analytical contact has no retained source-chord ancestry')


def prepare_pair(group,loop_index,join_index,owner,radius,*,attempted_radii=(),explicit_candidate=False,
                 spacing,tolerance,max_samples=8192):
    if (not math.isfinite(radius) or radius<.5 or not all(math.isfinite(x) and x>0. for x in (spacing,tolerance))
            or type(max_samples) is not int or not 4<=max_samples<=8192):
        raise ValueError('Invalid bounded retained-arc STYLE fillet request')
    inventory={source.region.owner:source for source in group.footprints}
    if owner in inventory:raise ValueError('Derived fillet repeats an existing source operand key')
    original_join=join_index
    loop=group.loops[loop_index]
    before=(join_index-1)%len(loop.xy)
    previous_owner=loop.spans[before].owner
    for _ in range(len(loop.xy)):
        if loop.spans[(before-1)%len(loop.xy)].owner!=previous_owner:break
        before=(before-1)%len(loop.xy)
    shift=-before
    loop=replace(loop,xy=np.roll(loop.xy,shift,axis=0),spans=tuple(np.roll(np.asarray(loop.spans,dtype=object),shift)))
    join_index=(join_index+shift)%len(loop.xy)

    left,right_span=loop.spans[join_index-1],loop.spans[join_index]
    ll,rl=rim_label(left,inventory),rim_label(right_span,inventory)
    if left.owner==right_span.owner or ll.kind!='curb' or rl.kind!='curb':
        raise ValueError('Retained-arc pair is not an exposed physical curb switch')
    records=[]
    for span,label in ((left,ll),(right_span,rl)):
        source=inventory[span.owner];corners=tuple(source.layout.corners)
        if len(corners)!=1:raise ValueError('Retained arc lacks one explicit original circular record')
        c=corners[0]
        if c.center is None or c.radius<=0. or c.sweep>=0. or c.radius_mode not in ('STYLE','STRICT'):
            raise ValueError('Retained fillet requires two actual clockwise original circles')
        curve=source.layout.boundaries[label.index]
        if isinstance(curve,OriginalAnalyticalCircleCurve):curve=curve.original_prepared_curve
        records.append((source,c,curve,OriginalAnalyticalCircleCurve(span.owner,label.index,curve,c)))
    (_,a,_,aview),(_,b,_,bview)=records
    required=max((float(c.requested_radius or c.radius) for c in (a,b) if c.radius_mode=="STRICT"),default=0.)
    if radius<required:
        raise RetainedArcFitUnavailable("Derived boundary fillet violates an original STRICT minimum")
    d=float(np.linalg.norm(b.center-a.center));ra,rb=a.radius-radius,b.radius-radius
    if min(ra,rb)<=0. or not abs(ra-rb)<d<ra+rb:
        raise RetainedArcFitUnavailable('No internally tangent circle inside both parent disks')
    ex=(b.center-a.center)/d;along=(ra*ra-rb*rb+d*d)/(2.*d)
    h=math.sqrt(max(0.,ra*ra-along*along));middle=a.center+along*ex
    sampled_crossing=loop.xy[join_index]
    original_along=(a.radius*a.radius-b.radius*b.radius+d*d)/(2.*d)
    original_h=math.sqrt(max(0.,a.radius*a.radius-original_along*original_along))
    original_middle=a.center+original_along*ex
    crossing=min((original_middle+original_h*right(ex),original_middle-original_h*right(ex)),
        key=lambda xy:np.linalg.norm(xy-sampled_crossing))
    center=min((middle+h*right(ex),middle-h*right(ex)),key=lambda xy:np.linalg.norm(xy-crossing))
    contact_a=a.center+a.radius/ra*(center-a.center)
    contact_b=b.center+b.radius/rb*(center-b.center)
    aref=(float(aview.circle_parameter(np.clip(ll.last,*aview.circle_bounds)))
          if aview.circle_interval is not None else ll.last)
    bref=(float(bview.circle_parameter(np.clip(rl.first,*bview.circle_bounds)))
          if bview.circle_interval is not None else rl.first)
    pa=phase(a,crossing,aref);pb=phase(b,crossing,bref)
    fa,fb=phase(a,contact_a,pa),phase(b,contact_b,pb)
    before=join_index-1
    while before>0 and loop.spans[before-1].owner==left.owner:before-=1
    after=join_index
    while after+1<len(loop.spans) and loop.spans[after+1].owner==right_span.owner:after+=1
    low=phase(a,loop.xy[before],0.);high=phase(b,loop.xy[(after+1)%len(loop.xy)],1.)
    if aview.circle_interval is not None:
        low=float(aview.circle_parameter(max(aview.circle_bounds[0],rim_label(loop.spans[before],inventory).first)))
    if bview.circle_interval is not None:
        high=float(bview.circle_parameter(min(bview.circle_bounds[1],rim_label(loop.spans[after],inventory).last)))
    if not low+1e-9<fa<pa or not pb<fb<high-1e-9:
        raise RetainedArcFitUnavailable('A fillet contact leaves its retained finite analytical arc')
    start=math.atan2(*(contact_a-center)[::-1]);end=math.atan2(*(contact_b-center)[::-1])
    sweep=-((start-end)%(2.*math.pi))
    if not -math.pi<sweep<0.:raise RetainedArcFitUnavailable('Fillet selects the wrong circular sector')
    ia,qa,ua,ta=_chord_contact(loop,before,join_index-1,a,fa,inventory,prepared_bounds=aview.circle_bounds)
    ib,qb,ub,tb=_chord_contact(loop,join_index,after,b,fb,inventory,prepared_bounds=bview.circle_bounds)
    if max(np.linalg.norm(ta-qa),np.linalg.norm(tb-qb))>tolerance+1e-9:
        raise ValueError('Retained analytical contact exceeds its original chord error budget')
    source_a=tuple(loop.spans[ia:join_index]);source_b=tuple(loop.spans[join_index:ib+1])
    arcs=(RetainedOriginalArc(left.owner,ll.index,records[0][2],a,(low,pa),(fa,pa),source_a,
        tuple(rim_label(s,inventory) for s in source_a),ua,tuple(qa),fa,tuple(ta),float(np.linalg.norm(ta-qa))),
        RetainedOriginalArc(right_span.owner,rl.index,records[1][2],b,(pb,high),(pb,fb),source_b,
        tuple(rim_label(s,inventory) for s in source_b),ub,tuple(qb),fb,tuple(tb),float(np.linalg.norm(tb-qb))))
    preferred=min(a.requested_radius or a.radius,b.requested_radius or b.radius)
    descriptor=RetainedArcFillet(owner,tuple(center),radius,start,sweep,arcs,preferred,
        tuple(attempted_radii),bool(explicit_candidate),tuple(sampled_crossing),loop_index,original_join,tuple(crossing),
        source_chord_tolerance=float(tolerance))
    count=max(2,int(math.ceil(descriptor.length/min(spacing,math.sqrt(8.*tolerance*radius)))))
    # The closure follows the true parent circles. Their original sampled
    # chords remain untouched area operands. Matching a finer angular phase
    # on these internal arcs prevents the closure polygon crossing the new
    # tangent fillet merely because a parent chord cuts its own circle.
    angle_step=abs(sweep)/count
    chains=[]
    for parent,curve,view,first,last in ((b,records[1][2],bview,fb,pb),(a,records[0][2],aview,pa,fa)):
        angle=abs((last-first)*parent.sweep)
        step=min(.25*angle_step,math.sqrt(8.*tolerance/parent.radius))
        samples=max(2,int(math.ceil(angle/step)))
        values=[(float(value),circular_point(parent,value),False) for value in np.linspace(first,last,samples+1)]
        # Every old point must remain a shared point. Omitting these knots
        # leaves microscopic gaps between two differently phased polygons,
        # despite their common authoritative analytical circle.
        declared=view.circle_interval
        old=(zip(declared.circle_parameters,curve.points[np.asarray(declared.vertex_indices)]) if declared is not None
             else ((phase(parent,point,.5*(first+last)),point) for point in curve.points))
        for value,point in old:
            if min(first,last)<value<max(first,last) and abs(np.linalg.norm(point-parent.center)-parent.radius)<1e-8:
                values.append((float(value),point,True))
        values.sort(key=lambda row:row[0],reverse=last<first)
        merged=[]
        for row in values:
            if merged and abs(row[0]-merged[-1][0])<1e-12:
                if row[2]:merged[-1]=row
            else:merged.append(row)
        chains.append(np.asarray([row[1] for row in merged]))
    closure=np.vstack((chains[0][:-1],chains[1]))
    keep=np.r_[True,np.linalg.norm(np.diff(closure,axis=0),axis=1)>1e-10];closure=closure[keep]
    if count+len(closure)>max_samples:raise ValueError('Retained fillet exceeds its finite sample budget')
    parameters=np.linspace(0.,1.,count+1);arc=descriptor(parameters)
    closure[0]=arc[-1];closure[-1]=arc[0]
    internal=PreparedCapCurb(closure)
    points=np.vstack((arc[:-1],closure[:-1]))
    if signed_area(points)<=1e-10:raise ValueError('Fillet does not add a positive extensive pavement wedge')
    labels=tuple([RimSegment('curb',0,float(x),float(y)) for x,y in zip(parameters,parameters[1:])]
        +[RimSegment('internal',1,float(x/internal.length),float(y/internal.length))
          for x,y in zip(internal.stations,internal.stations[1:])])
    layout=SimpleNamespace(approaches=(),boundaries=(descriptor,internal),rows=len(points),mouths={},corners=(),
        retained_arc_fillets=(descriptor,),original_source_arc_ancestry=arcs,generated_junction_area=descriptor)
    footprint=SampledFootprint(Region(owner,points),labels,layout)
    RegionUnion().build((footprint.region,))
    return descriptor,footprint


def auxiliary_candidates(first, second):
    """Ordinary STYLE fitting first, then the existing rounded-gore tail."""
    from .corner_radii import CornerRadiusPolicy
    preferred=min(first.requested_radius or first.radius,second.requested_radius or second.radius)
    strict=max((c.requested_radius or c.radius for c in (first,second)
                if c.radius_mode=='STRICT'),default=0.)
    if strict:
        yield float(max(preferred,strict)),False
        return
    ordinary=tuple(CornerRadiusPolicy.candidates(preferred,'STYLE'))
    yield from ((float(radius),False) for radius in ordinary)
    previous=ordinary[-1]
    for _ in range(CornerRadiusPolicy.maximum_reductions):
        if previous<=CornerRadiusPolicy.minimum_gore_radius:return
        radius=max(CornerRadiusPolicy.minimum_gore_radius,previous*CornerRadiusPolicy.ratio)
        yield float(radius),True
        previous=radius
    if previous>CornerRadiusPolicy.minimum_gore_radius:
        yield float(CornerRadiusPolicy.minimum_gore_radius),True


def _original_circle(source,index):
    if source.region.owner.kind!='cycle_corner' or len(getattr(source.layout,'corners',()))!=1:
        return None
    corner=source.layout.corners[0]
    if corner.center is None or corner.radius<=0. or corner.sweep>=0.:
        return None
    return corner


def _circle(source,index,parameter=None):
    curve=source.layout.boundaries[index]
    if isinstance(curve,RetainedArcFillet):return curve.corner
    from .source_bank_reconciliation import PreparedBankFillet
    marker=getattr(source.layout,'generated_junction_area',None)
    if isinstance(marker,PreparedBankFillet) and hasattr(curve,'points'):
        corner=marker.fitted.corner
        if (corner.center is not None and corner.sweep<0.
                and np.max(abs(np.linalg.norm(np.asarray(curve.points)-corner.center,axis=1)-corner.radius))<=1e-8):
            return corner
    circle=_original_circle(source,index)
    declared=(curve.circle_interval if isinstance(curve,OriginalAnalyticalCircleCurve)
              else getattr(curve,'original_circle_interval',None))
    if circle is not None and declared is not None and parameter is not None:
        low,high=declared.prepared_bounds
        if not low-1e-12<=parameter<=high+1e-12:return None
    return circle


def _intersection(first,last,near):
    ca,cb=np.asarray(first.center),np.asarray(last.center)
    distance=float(np.linalg.norm(cb-ca))
    if distance<=1e-10:return None
    along=(first.radius**2-last.radius**2+distance**2)/(2.*distance)
    height2=first.radius**2-along**2
    if height2 < -1e-8:return None
    axis=(cb-ca)/distance;middle=ca+along*axis;height=math.sqrt(max(0.,height2))
    return min((middle+height*right(axis),middle-height*right(axis)),key=lambda point:np.linalg.norm(point-near))


def reflex_switches(group, *, tolerance):
    """Read actual analytical turns, excluding polygon phase/chord artifacts."""
    sources={source.region.owner:source for source in group.footprints}
    result=[]
    for li,loop in enumerate(group.loops):
        for at,last in enumerate(loop.spans):
            first=loop.spans[at-1]
            if first.owner==last.owner:continue
            labels=(rim_label(first,sources),rim_label(last,sources))
            if any(label.kind!='curb' for label in labels):continue
            circles=tuple(_circle(sources[span.owner],label.index,parameter)
                          for span,label,parameter in zip((first,last),labels,(labels[0].last,labels[1].first)))
            if any(circle is None for circle in circles):continue
            point=_intersection(*circles,loop.xy[at])
            if point is None or np.linalg.norm(point-loop.xy[at])>4.*tolerance+1e-7:continue
            tangents=tuple(right((point-np.asarray(circle.center))/circle.radius) for circle in circles)
            turn=math.atan2(float(np.linalg.det(np.stack(tangents))),float(tangents[0]@tangents[1]))
            if turn>=-math.radians(1e-4):continue
            if any(isinstance(sources[span.owner].layout.boundaries[label.index],RetainedArcFillet)
                   for span,label in zip((first,last),labels)):
                raise RetainedArcReconciliationError('A derived fillet creates another genuine analytical reflex switch')
            result.append((li,at,first.owner,last.owner,float(turn)))
    return tuple(result)


def _mouth_points(group):
    sources={source.region.owner:source for source in group.footprints};result={}
    for loop in group.loops:
        for at,span in enumerate(loop.spans):
            label=rim_label(span,sources)
            if label.kind!='mouth':continue
            for value,point in ((label.first,loop.xy[at]),(label.last,loop.xy[(at+1)%len(loop.xy)])):
                key=span.owner,label.index,round(value,12)
                if key in result and not np.array_equal(result[key],point):
                    raise RetainedArcReconciliationError('One original profile parameter names two physical points')
                result[key]=point
    return result


def _topology(original,current,tolerance):
    old=[loop for loop in original.loops if loop.hole]
    new=[loop for loop in current.loops if loop.hole]
    if (sum(not loop.hole for loop in current.loops)!=1 or len(old)!=len(new)
            or any(-loop.area<=1e-8 for loop in new)):
        raise RetainedArcReconciliationError('Retained-arc fill changes a genuine hole topology')
    from .regions import _classify
    matched=set()
    for loop in new:
        parents=[]
        for i,parent in enumerate(old):
            inside,boundary=_classify(loop.xy,parent.xy,2.*tolerance+1e-7)
            if np.all(inside | (boundary>=0)):parents.append(i)
        if len(parents)!=1 or parents[0] in matched:
            raise RetainedArcReconciliationError('Retained-arc fill loses a distinct original hole')
        matched.add(parents[0])


def _view(group):
    sources={source.region.owner:source for source in group.footprints}
    parents={};records={}
    for fillet in group.physical_layout.retained_arc_fillets:
        source=sources[fillet.owner]
        parameters=tuple([source.segments[0].first,*[s.last for s in source.segments if s.kind=='curb']])
        records[fillet.owner]=replace(fillet,source_reference_parameters=parameters)
        for parent in fillet.parent_arcs:
            token=parent.owner,parent.curve_index
            if token not in parents:
                parents[token]=OriginalAnalyticalCircleCurve(*token,parent.original_curve,parent.original_corner)
    footprints=[]
    for source in group.footprints:
        curves=list(source.layout.boundaries);changed=False
        for i,curve in enumerate(curves):
            if (source.region.owner,i) in parents:
                curves[i]=parents[source.region.owner,i];changed=True
            elif source.region.owner in records and i==0:
                curves[i]=records[source.region.owner];changed=True
        if changed:
            layout=copy(source.layout);layout.boundaries=tuple(curves)
            if source.region.owner in records:
                layout.retained_arc_fillets=(records[source.region.owner],)
                layout.generated_junction_area=records[source.region.owner]
            source=replace(source,layout=layout)
        footprints.append(source)
    physical=copy(group.physical_layout);physical.retained_arc_fillets=tuple(records.values())
    return replace(group,footprints=tuple(footprints),physical_layout=physical)


def shared_contact(sampler,first,last,labels,layouts):
    """Use explicit finite parent parameters, never nearest-curve ancestry."""
    spans=(first,last)
    curves=[layouts[span.owner].boundaries[label.index] for span,label in zip(spans,labels)]
    indices=[i for i,curve in enumerate(curves) if isinstance(curve,RetainedArcFillet)]
    if len(indices)!=1:raise RetainedArcReconciliationError('Retained contact needs exactly one derived fillet')
    fi=indices[0];pi=1-fi;fillet=curves[fi];parent_curve=curves[pi]
    if not isinstance(parent_curve,OriginalAnalyticalCircleCurve):
        raise RetainedArcReconciliationError('Retained contact lost its explicit analytical parent')
    matches=[(i,parent) for i,parent in enumerate(fillet.parent_arcs)
        if parent.owner==spans[pi].owner and parent.curve_index==labels[pi].index
        and parent.original_corner is parent_curve.corner
        and parent.original_curve is parent_curve.original_prepared_curve]
    if len(matches)!=1:raise RetainedArcReconciliationError('Retained contact changed its finite original ancestry')
    which,parent=matches[0];token=fillet.owner,which
    tokens=sampler.retained_arc_contact_tokens
    if token in tokens or len(tokens)>=32:
        raise RetainedArcReconciliationError('Retained contacts exceed their finite pair allocation')
    values=np.empty(2);values[fi]=float(which)
    values[pi]=float(parent_curve.source_parameter(parent.analytical_contact_parameter))
    low,high=sorted(parent.exposed_circle_parameters)
    if not low<=parent.analytical_contact_parameter<=high or not 0.<=values[fi]<=1.:
        raise RetainedArcReconciliationError('Retained contact leaves its original finite arc')
    initial=np.array((labels[0].last,labels[1].first))
    reach=np.array([2./layouts[span.owner].rows for span in spans])
    outside=(values<np.maximum(0.,initial-reach)-1e-9)|(values>np.minimum(1.,initial+reach)+1e-9)
    # The wedge's row count includes its internal closing polygon. The actual
    # terminal fillet chord supplies a local interval independent of that count.
    span=spans[fi]
    if outside[fi] and which==1-fi and isinstance(span,SourceSpan):
        knots=fillet.source_reference_parameters
        terminal=len(knots)-2 if which else 0
        if (len(knots)>=2 and span.segment==terminal and knots[0]==0. and knots[-1]==1.
                and 0.<=span.first<=span.last<=1.):
            a,b=knots[terminal:terminal+2]
            expected=a+(b-a)*np.array([span.first,span.last])
            clipped=np.array([labels[fi].first,labels[fi].last])
            if b>a and np.allclose(clipped,expected,atol=1e-12,rtol=0.):
                outside[fi]=False
    if np.any(outside):
        raise RetainedArcReconciliationError('Retained contact leaves the original ordinary contact neighborhood')
    if np.linalg.norm(curves[0](values[0])-curves[1](values[1]))>1e-7:
        raise RetainedArcReconciliationError('Explicit parent parameters do not name the same analytical contact')
    sampler._check_contact_budget(curves,values,spans,layouts)
    tolerance=fillet.source_chord_tolerance
    # The declared descriptor also records each original prepared normal error.
    for curve,value in zip(curves,values):
        points=np.asarray(curve.contact_reference_points);edges=np.diff(points,axis=0)
        length2=np.sum(edges*edges,axis=1)
        if np.any(length2<=1e-16):raise RetainedArcReconciliationError('Retained source chord collapsed')
        delta=curve(value)-points[:-1];u=np.clip(np.sum(delta*edges,axis=1)/length2,0.,1.)
        error=float(np.min(np.linalg.norm(delta-u[:,None]*edges,axis=1)))
        if error>tolerance+1e-7:
            raise RetainedArcReconciliationError('Retained contact exceeds an original prepared chord budget')
    tokens.add(token)
    sampler.retained_arc_contact_handoffs.append(dict(fillet=fillet.owner.text(),parent=parent.owner.text(),
        parameters=values.tolist(),proof='EXACT_PARENT_OBJECT_OWNER_INDEX_AND_FINITE_PARAMETERS'))
    return values


def canonicalize_tangent_fragments(group,sampler,*,context=None,maximum_proof_points=4096):
    """Resolve chord-phase re-entries at explicit analytical tangencies.

    An internally tangent fillet has exactly one contact with each original
    parent circle. A polygon chord can nevertheless clip a short disconnected
    fragment near that tangency. Only spans of the recorded consumed parents
    may be replaced; every other source switch remains a hard rejection.
    """
    fillets=tuple(getattr(group.physical_layout,'retained_arc_fillets',()))
    if len(fillets)>16:
        raise RetainedArcReconciliationError('Retained tangent resolution exceeds its finite pair allocation')
    if type(maximum_proof_points) is not int or not 3<=maximum_proof_points<=4096:
        raise ValueError('Hidden tangent proof needs its bounded original allocation')
    sources={source.region.owner:source for source in group.footprints}
    layouts={owner:source.layout for owner,source in sources.items()}
    loops=list(group.loops);journal=[];budget=dict(maximum=maximum_proof_points,used=0)
    checker=copy(sampler)
    checker.retained_arc_contact_tokens=set();checker.retained_arc_contact_handoffs=[]
    for fillet in fillets:
        uses=[li for li,loop in enumerate(loops) if any(span.owner==fillet.owner for span in loop.spans)]
        if not uses:
            from .retained_tangent_visibility import restore_hidden
            current=replace(group,loops=tuple(loops))
            current,record=restore_hidden(current,fillet,checker,context=context,budget=budget)
            loops=list(current.loops);journal.append(record)
            continue
        if len(uses)!=1:
            raise RetainedArcReconciliationError('Retained fillet fragments lost their single physical boundary')
        li=uses[0];loop=loops[li]
        ids=[i for i,span in enumerate(loop.spans) if span.owner==fillet.owner]
        runs=sum(loop.spans[i-1].owner!=fillet.owner for i in ids)
        if runs<=1:continue
        labels=[rim_label(span,sources) for span in loop.spans]
        begin=min(ids,key=lambda i:labels[i].first);end=max(ids,key=lambda i:labels[i].last)
        order=np.roll(np.arange(len(loop.spans)),-begin)
        stop=int(np.flatnonzero(order==end)[0])+1
        if stop>=len(order) or not set(ids)<=set(order[:stop]):
            raise RetainedArcReconciliationError('Retained tangent fragments reverse their finite arc order')
        previous,following=loop.spans[order[-1]],loop.spans[order[stop]]
        parents=fillet.parent_arcs
        if any(span.owner!=parent.owner or label.index!=parent.curve_index
               for span,label,parent in ((previous,labels[order[-1]],parents[0]),
                                          (following,labels[order[stop]],parents[1]))):
            raise RetainedArcReconciliationError('Retained tangent fragments meet a foreign source boundary')
        removed=[]
        for index in order[:stop]:
            span=loop.spans[index];label=labels[index]
            if span.owner==fillet.owner:continue
            matches=[parent for parent in parents
                     if span.owner==parent.owner and label.index==parent.curve_index]
            if len(matches)!=1 or label.kind!='curb':
                raise RetainedArcReconciliationError('Retained tangent resolution cannot replace a third source')
            parent=matches[0];low,high=sorted(parent.consumed_circle_parameters)
            curve=layouts[parent.owner].boundaries[parent.curve_index]
            actual=(tuple(map(float,curve.circle_parameter(np.asarray((label.first,label.last)))))
                    if isinstance(curve,OriginalAnalyticalCircleCurve) else (label.first,label.last))
            if (not isinstance(curve,OriginalAnalyticalCircleCurve)
                    or curve.corner is not parent.original_corner
                    or curve.original_prepared_curve is not parent.original_curve
                    or not low-1e-9<=actual[0]<actual[1]<=high+1e-9
                    or abs(np.linalg.norm(np.asarray(fillet.center)-curve.corner.center)
                           +fillet.radius-curve.corner.radius)>1e-7):
                raise RetainedArcReconciliationError('Retained re-entry is not its consumed internally tangent parent')
            removed.append((span.owner.text(),label.index,label.first,label.last))
        if not removed:
            raise RetainedArcReconciliationError('Retained tangent fragments have no consumed parent alias')
        # The full ordinary ancestry, local contact neighbourhood and BOTH
        # prepared chord budgets are checked independently at BOTH endpoints.
        shared_contact(checker,previous,loop.spans[order[0]],
                       (labels[order[-1]],labels[order[0]]),layouts)
        shared_contact(checker,loop.spans[order[stop-1]],following,
                       (labels[order[stop-1]],labels[order[stop]]),layouts)
        contacts=[float(layouts[parent.owner].boundaries[parent.curve_index].source_parameter(
                         parent.analytical_contact_parameter)) for parent in parents]
        if (labels[order[-1]].first>=contacts[0] or labels[order[stop]].last<=contacts[1]):
            raise RetainedArcReconciliationError('Retained tangent resolution consumes an original finite parent endpoint')
        kept=[loop.spans[i] for i in order[stop:]]
        kept[0]=RimSpan(kept[0].owner,replace(labels[order[stop]],first=contacts[1]))
        kept[-1]=RimSpan(kept[-1].owner,replace(labels[order[-1]],last=contacts[0]))
        spans=(RimSpan(fillet.owner,RimSegment('curb',0,0.,1.)),*kept)
        xy=np.vstack((fillet(0.),fillet(1.),loop.xy[order[stop+1:]]))
        if len(xy)!=len(spans):
            raise RetainedArcReconciliationError('Retained tangent resolution lost its boundary incidence')
        loops[li]=BoundaryLoop(xy,tuple(spans))
        journal.append(dict(fillet=fillet.owner.text(),consumed_parent_fragments=tuple(removed),
            original_vertices=len(loop.xy),canonical_vertices=len(xy),
            proof='EXACT_INTERNAL_CIRCLE_TANGENCY_CONSUMED_PARENT_AND_BOTH_SOURCE_CHORD_BUDGETS'))
    return (replace(group,loops=tuple(loops)),tuple(journal)) if journal else (group,())


def check_consumed_reservations(candidate,reserved):
    """Keep a positive original arc between exact-parent auxiliary joins."""
    for parent in candidate.parent_arcs:
        for previous in reserved:
            for other in previous.parent_arcs:
                if (parent.owner,parent.curve_index)!=(other.owner,other.curve_index):continue
                if (parent.original_curve is not other.original_curve
                        or parent.original_corner is not other.original_corner):
                    raise RetainedArcReconciliationError('Reserved parent interval changed its original analytical ancestry')
                a,b=sorted(parent.consumed_circle_parameters)
                c,d=sorted(other.consumed_circle_parameters)
                gap=max(a,c)-min(b,d)
                length=parent.original_corner.radius*abs(parent.original_corner.sweep)
                if gap*length<=1e-7:
                    raise RetainedArcFitUnavailable('Fillet consumes another reserved finite original parent interval')


def _original_parents(planner,group,sources,switch):
    from .displaced_corner_authority import original_parent_corner
    li,at,first,last,turn=switch;context=planner.context;parents=[]
    for owner,index in ((first,rim_label(group.loops[li].spans[at-1],sources).index),
                        (last,rim_label(group.loops[li].spans[at],sources).index)):
        source=sources[owner];corner=_original_circle(source,index)
        node=next((node for node in context.network.intersections
                   if node.key==getattr(source.layout,'original_junction',None)),None)
        if (corner is None or node is None or node.key not in group.members
                or not original_parent_corner(SimpleNamespace(context=context,group=group),node,corner)):
            raise RetainedArcReconciliationError('Retained parent lost its actual source radius/incidence authority')
        curve=source.layout.boundaries[index]
        if isinstance(curve,OriginalAnalyticalCircleCurve):curve=curve.original_prepared_curve
        OriginalAnalyticalCircleCurve(owner,index,curve,corner)
        parents.append(corner)
    return tuple(parents)


def reconcile(planner,group,*,maximum_pairs=16,maximum_batches=2):
    from .updates import drain_steps
    return drain_steps(reconcile_steps(planner, group, maximum_pairs=maximum_pairs, maximum_batches=maximum_batches))


def reconcile_steps(planner,group,*,maximum_pairs=16,maximum_batches=2):
    """Bounded extensive rounding of genuine exposed original-circle switches."""
    if (type(maximum_pairs) is not int or not 1<=maximum_pairs<=16
            or type(maximum_batches) is not int or not 1<=maximum_batches<=2):
        raise ValueError('Retained-arc closure needs its bounded original work allocation')
    context=planner.context;settings=context.network.settings
    tolerance=min(.001,.05*settings.curve_tolerance)
    original=group;mouths=_mouth_points(original);visited=set();journal=[]
    from .sampling import AnalyticalBoundarySampler
    from .compound_ports import CompoundPortPlanner
    from .cycle_dimensional_retry import certify_layout_steps
    for batch in range(maximum_batches):
        switches=reflex_switches(group,tolerance=tolerance)
        if not switches:return group
        if len(visited)+len(switches)>maximum_pairs:
            raise RetainedArcReconciliationError('Retained-arc closure exceeds its visited-pair budget')
        sources={source.region.owner:source for source in group.footprints}
        existing=tuple(getattr(group.physical_layout,'retained_arc_fillets',()))
        if len(existing)+len(switches)>maximum_pairs:
            raise RetainedArcReconciliationError('Retained-arc closure exceeds its finite total pair budget')
        parents_by_pair=[_original_parents(planner,group,sources,switch) for switch in switches]
        strict_prepared={};strict_reserved=[]
        # Reserve hard intervals before adjustable choices, but retain the
        # original publication order, owners and single STRICT radius attempt.
        for index,(switch,parents) in enumerate(zip(switches,parents_by_pair)):
            if not any(parent.radius_mode=='STRICT' for parent in parents):continue
            li,at,first,last,turn=switch
            radius,tail=next(auxiliary_candidates(*parents))
            owner=Key('retained_arc_fillet',len(existing)+index)
            if owner in sources:raise RetainedArcReconciliationError('Retained arc owner namespace is already occupied')
            try:
                prepared=prepare_pair(group,li,at,owner,radius,attempted_radii=(radius,),
                    explicit_candidate=tail,spacing=min(.2,settings.spacing),tolerance=tolerance)
                check_consumed_reservations(prepared[0],(*existing,*strict_reserved))
            except RetainedArcFitUnavailable as failure:
                raise RetainedArcReconciliationError('No allowed finite retained-circle fillet meets the original radius constraints (STRICT)') from failure
            strict_prepared[index]=prepared;strict_reserved.append(prepared[0])
        fillets=[];wedges=[]
        reservation_declines=[]
        for index,(switch,parents) in enumerate(zip(switches,parents_by_pair)):
            li,at,first,last,turn=switch
            token=first,last
            if token in visited:raise RetainedArcReconciliationError('Retained-arc closure repeats an original pair')
            visited.add(token)
            owner=Key('retained_arc_fillet',len(existing)+len(fillets))
            if owner in sources:raise RetainedArcReconciliationError('Retained arc owner namespace is already occupied')
            attempts=[];selected=None
            if index in strict_prepared:
                selected=strict_prepared[index];attempts=[selected[0].radius]
            else:
                for radius,tail in auxiliary_candidates(*parents):
                    if len(attempts)>=17:raise RetainedArcReconciliationError('Retained radius search exceeds its finite attempt budget')
                    attempts.append(radius)
                    try:
                        prepared=prepare_pair(group,li,at,owner,radius,attempted_radii=attempts,
                            explicit_candidate=tail,spacing=min(.2,settings.spacing),tolerance=tolerance)
                        try:check_consumed_reservations(prepared[0],(*existing,*strict_reserved,*fillets))
                        except RetainedArcFitUnavailable:
                            reservation_declines.append(dict(owner=owner.text(),radius=radius,
                                reason='OVERLAPPING_CONSUMED_ORIGINAL_PARENT_INTERVAL'))
                            raise
                        selected=prepared
                    except RetainedArcFitUnavailable:continue
                    break
            if selected is None:
                if len(switches)==1:
                    from .retained_circle_bank_join import recover_steps
                    repaired=yield from recover_steps(planner,group,switch)
                    if repaired is not None and not reflex_switches(repaired,tolerance=tolerance):
                        from .exposed_bank_joins import reconcile_steps as round_remaining_banks
                        return (yield from round_remaining_banks(planner,repaired,check_tangency=True))
                raise RetainedArcReconciliationError('No allowed finite retained-circle fillet meets the original radius constraints')
            fillet,wedge=selected
            fillets.append(fillet);wedges.append(wedge)
            journal.append(dict(previous=first.text(),following=last.text(),loop=li,
                analytical_turn_degrees=math.degrees(turn),attempted_radii=attempts,radius=fillet.radius))
        operands=(*group.area_footprints,*wedges)
        if sum(len(source.region.xy) for source in operands)>20000:
            raise RetainedArcReconciliationError('Retained-arc union exceeds its unchanged source segment budget')
        loops=yield from RegionUnion(max_segments=20000).build_steps(source.region for source in operands)
        physical=copy(group.physical_layout);physical.retained_arc_fillets=(*existing,*fillets)
        result=replace(group,footprints=(*group.footprints,*wedges),loops=loops,physical_layout=physical)
        _topology(original,result,tolerance)
        sampler=AnalyticalBoundarySampler(tolerance=settings.curve_tolerance,
            spacing=settings.spacing,max_samples=16384)
        sampler.source_curve_tolerance=settings.curve_tolerance
        view,tangent_fragments=canonicalize_tangent_fragments(_view(result),sampler,context=context)
        result=sampler.build(view)
        _topology(original,result,tolerance)
        inventory={source.region.owner:source for source in result.footprints}
        for fillet in result.physical_layout.retained_arc_fillets:
            labels=[rim_label(span,inventory) for loop in result.loops for span in loop.spans if span.owner==fillet.owner]
            if (not labels or any(label.kind!='curb' for label in labels)
                    or abs(min(label.first for label in labels))>1e-12
                    or abs(max(label.last for label in labels)-1.)>1e-12
                    or abs(sum(label.last-label.first for label in labels)-1.)>1e-12):
                raise RetainedArcReconciliationError('Whole source union clips a retained fillet before its true contacts')
        actual=_mouth_points(result)
        if mouths.keys()!=actual.keys() or any(not np.array_equal(point,actual[key]) for key,point in mouths.items()):
            raise RetainedArcReconciliationError('Retained-arc union changed an original complete profile')
        if CompoundPortPlanner.incomplete_ports(result):
            raise RetainedArcReconciliationError('Retained-arc union clips an original incoming profile')
        consumed={claim.section_id for claim in result.claims}
        CompoundPortPlanner(context,tuple(planner.sites.values()),{}).check_support(consumed,result.loops[0])
        result.physical_layout.source_cycle_full_width_proof=yield from certify_layout_steps(context,result,consumed)
        result.physical_layout.source_retained_arc_policy=dict(policy=SOURCE_ARC_POLICY,
            batches=batch+1,visited_pairs=len(visited),maximum_pairs=maximum_pairs,
            maximum_batches=maximum_batches,maximum_radius_attempts=17,journal=tuple(journal),
            exact_parent_contacts=tuple(sampler.retained_arc_contact_handoffs),
            original_logical_corners_unchanged=result.physical_layout.corners is original.physical_layout.corners)
        if tangent_fragments:
            result.physical_layout.source_retained_arc_policy['canonical_tangent_fragments']=tangent_fragments
        if reservation_declines:
            result.physical_layout.source_retained_arc_policy['consumed_interval_fit_declines']=tuple(reservation_declines)
        group=result
    if reflex_switches(group,tolerance=tolerance):
        raise RetainedArcReconciliationError('Retained-arc closure exhausted its finite batch budget')
    return group
