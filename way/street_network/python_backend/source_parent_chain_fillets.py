"""Finite rounding of a certified core→original circle→its native bank.

Only the observed consecutive source chain is eligible. The original circle
and certified polynomial stay unchanged; an additional circular boundary closes
along their existing physical source chain. No new road incidence is created.
"""
from dataclasses import dataclass
import math
import numpy as np
from .generated_junction_curbs import GeneratedJunctionArea,is_generated_junction_area
from .source_bank_reconciliation import SourceBankDiscovery,NativeBankWitness
from .source_owned_connector_fillets import NativeCurbContact,_policy
from .source_circular_intervals import PreparedCircularInterval
from .compound_dimensional_corners import DimensionalCompositeCurb,_CertifiedContinuation
from .continuations import WidthContinuation
from .compounds import rim_label
from .corner_radii import CornerRadiusPolicy
from .paths import Path,right
from .model import End

SOURCE_PARENT_CHAIN_POLICY='1-finite-certified-core-original-circle-native-bank'


class ParentChainUnavailable(ValueError):
    """A genuine observed chain has no admitted finite original-source fit."""


@dataclass(frozen=True)
class CertifiedCoreContact:
    owner: object
    index: int
    curve: object
    parameter: float
    visible_parameters: tuple
    original_corner: object
    original_site: object


@dataclass(frozen=True)
class OriginalCircleChain:
    owner: object
    index: int
    curve: object
    interval: PreparedCircularInterval
    corner: object
    original_site: object


@dataclass(frozen=True)
class PreparedParentChainFillet(GeneratedJunctionArea):
    owner: object
    center: tuple
    radius: float
    start_angle: float
    sweep: float
    first: tuple
    last: tuple
    core_contact: CertifiedCoreContact
    native_contact: NativeCurbContact
    parent_circle: OriginalCircleChain
    closure_ancestry: tuple
    preferred_radius: float
    radius_mode: str
    attempted_radii: tuple
    source_reference_parameters: tuple
    source_chord_tolerance: float
    regularity_proof: tuple

    @property
    def length(self):return self.radius*abs(self.sweep)
    @property
    def lengths(self):return np.asarray((0.,0.,self.length,self.length))
    @property
    def curvature_bound(self):return 1./self.radius
    @property
    def corner(self):
        from .retained_arc_fillets import CircularBoundaryGeometry
        return CircularBoundaryGeometry(self.center,self.radius,self.start_angle,self.sweep)
    @property
    def leads(self):return ()
    @property
    def continuation(self):return None
    @property
    def contact_reference_points(self):return self(np.asarray(self.source_reference_parameters))

    def __call__(self,t):
        values=np.asarray(t);angle=self.start_angle+self.sweep*values
        points=np.asarray(self.center)+self.radius*np.stack((np.cos(angle),np.sin(angle)),axis=-1)
        points=np.where((values==0.)[...,None],self.first,points)
        return np.where((values==1.)[...,None],self.last,points)

    def tangent(self,t):
        angle=self.start_angle+self.sweep*np.asarray(t)
        return np.sign(self.sweep)*np.stack((-np.sin(angle),np.cos(angle)),axis=-1)

    def certified_parameters(self,sampler,first,last):
        budget=min(sampler.tolerance,self.source_chord_tolerance)
        step=min(sampler.spacing,math.sqrt(8.*budget*self.radius))
        regular=sampler.regular(first,last,step/self.length)
        low,high=sorted((first,last));knots=np.asarray(self.source_reference_parameters)
        values=np.unique(np.r_[regular,knots[(knots>low)&(knots<high)]])
        return values if last>=first else values[::-1]

    def report(self):
        core=self.core_contact;native=self.native_contact;circle=self.parent_circle
        return dict(owner=self.owner.text(),policy=SOURCE_PARENT_CHAIN_POLICY,
            radius=self.radius,preferred_radius=self.preferred_radius,radius_mode=self.radius_mode,
            attempted_radii=self.attempted_radii,center=self.center,start_angle=self.start_angle,sweep=self.sweep,
            core_owner=core.owner.text(),core_index=core.index,core_parameter=core.parameter,
            core_visible_parameters=core.visible_parameters,core_site=core.original_site.text(),
            original_circle_owner=circle.owner.text(),original_circle_index=circle.index,
            original_circle_site=circle.original_site.text(),
            native_section=native.path.section.id,native_side=native.side,native_station=native.station,
            native_owned_intervals=native.owned_intervals,regularity_proof=self.regularity_proof,
            source_chord_tolerance=self.source_chord_tolerance,original_sources_unchanged=True,
            closure_ancestry=[(owner.text(),label.kind,label.index,label.first,label.last)
                              for owner,label in self.closure_ancestry])


@dataclass(frozen=True)
class _Candidate:
    group: object
    core_source: object
    core_index: int
    core_curve: object
    piece_index: int
    polynomial: object
    polynomial_interval: tuple
    core_visible_parameters: tuple
    core_site: object
    circle_source: object
    circle_index: int
    circle_curve: object
    circle: object
    circle_site: object
    bank: NativeBankWitness
    bank_arm: object
    bank_outward_interval: tuple
    radius: float
    mode: str


def _original_policy(planner,source,corner):
    from .junction_bank_transitions import original_corner
    site=getattr(source.layout,'original_junction',None)
    node=planner.sites.get(site)
    if node is None or site not in planner._chain_members:return None
    if not original_corner(planner,node,corner):return None
    endpoint=End(corner.first.source.id,corner.first.end)
    radius,mode=node.corner_radius(endpoint)
    return site,endpoint,float(radius),mode


def discover(planner,group,*,maximum_pairs=16):
    """Actual adjacent typed parents only; never infer a source from XY."""
    if type(maximum_pairs) is not int or not 1<=maximum_pairs<=16:raise ValueError('Invalid finite parent-chain pair budget')
    d=SourceBankDiscovery(planner.context,group);sources=d.sources;result=[]
    from types import SimpleNamespace
    authority=SimpleNamespace(context=planner.context,sites=planner.sites,_chain_members=group.members)
    claims={claim.section_id for claim in group.claims}
    for li,loop in enumerate(group.loops):
        for at,span in enumerate(loop.spans):
            prior=loop.spans[at-1]
            if prior.owner==span.owner:continue
            a,b=sources[prior.owner],sources[span.owner]
            if is_generated_junction_area(a) or is_generated_junction_area(b):continue
            labels=(rim_label(prior,sources),rim_label(span,sources))
            if any(label.kind!='curb' for label in labels):continue
            curve=a.layout.boundaries[labels[0].index];parent=b.layout.boundaries[labels[1].index]
            declared=getattr(parent,'original_circle_interval',None)
            if type(curve) is not DimensionalCompositeCurb or type(declared) is not PreparedCircularInterval:continue
            circle=declared.original_corner;declared.validate(parent,circle)
            policies=(_original_policy(authority,a,curve.corner),_original_policy(authority,b,circle))
            if any(value is None for value in policies):continue
            policy=_policy(policies)
            if policy is None:continue
            distance=labels[0].last*curve.length
            pi=int(np.clip(np.searchsorted(curve.lengths,distance,side='right')-1,0,len(curve.pieces)-1))
            piece=curve.pieces[pi]
            if type(piece) is not _CertifiedContinuation or type(piece.curve.continuation) is not WidthContinuation:continue
            poly=piece.curve.continuation
            if circle.center is None or abs(circle.sweep)<=1e-10:continue
            phase=float(declared.circle_parameter(labels[1].first))
            if not 0.<phase<1.:continue
            tangent=poly.derivative((distance-curve.lengths[pi])/piece.length)
            tangent/=np.linalg.norm(tangent)
            angle=circle.start_angle+phase*circle.sweep
            following=np.sign(circle.sweep)*np.array((-math.sin(angle),math.cos(angle)))
            turn=math.atan2(float(tangent[0]*following[1]-tangent[1]*following[0]),float(tangent@following))
            if turn>=-math.radians(1e-4):continue
            finish=at
            while ((loop.spans[(finish+1)%len(loop.spans)].owner,
                    rim_label(loop.spans[(finish+1)%len(loop.spans)],sources).index)==(span.owner,labels[1].index)):
                finish=(finish+1)%len(loop.spans)
                if finish==at:raise ParentChainUnavailable('Parent circle occupies a complete physical boundary')
            tail=rim_label(loop.spans[finish],sources)
            if abs(tail.last-declared.prepared_bounds[1])>1e-9:continue
            ni=(finish+1)%len(loop.spans);bankspan=loop.spans[ni];bank=d.describe(bankspan,False)
            original_bank=d.corner_bank(b,circle,False)
            if (not isinstance(bank,NativeBankWitness) or original_bank is None
                    or not d.same_bank(original_bank,bank) or bankspan.owner.kind!='cycle_strip'
                    or bank.path.section.id not in claims):continue
            end=ni
            while loop.spans[(end+1)%len(loop.spans)].owner==bankspan.owner:
                nextlabel=rim_label(loop.spans[(end+1)%len(loop.spans)],sources)
                if nextlabel.kind!='curb' or nextlabel.index!=rim_label(bankspan,sources).index:break
                end=(end+1)%len(loop.spans)
                if end==ni:raise ParentChainUnavailable('Native bank occupies a complete physical boundary')
            last=d.describe(loop.spans[end],True)
            if not isinstance(last,NativeBankWitness) or not d.same_bank(bank,last):continue
            arm=circle.second
            from .routes import Route
            if (type(arm.path) is not Route or len(arm.path.pieces)!=1
                    or arm.path.pieces[0].path is not bank.path
                    or arm.path.pieces[0].forward!=(arm.end==0)):continue
            first=float(circle.stations[1]);stop=float(last.station if arm.end==0 else arm.length-last.station)
            if not 0.<=first<stop<=arm.length:continue
            if not any(low<=min(original_bank.station,last.station)<=max(original_bank.station,last.station)<=high
                       for low,high in bank.intervals):continue
            before=(at-1)%len(loop.spans)
            while ((loop.spans[(before-1)%len(loop.spans)].owner,
                    rim_label(loop.spans[(before-1)%len(loop.spans)],sources).index)==(prior.owner,labels[0].index)):
                before=(before-1)%len(loop.spans)
                if before==(at-1)%len(loop.spans):raise ParentChainUnavailable('Certified core occupies a complete physical boundary')
            entry=rim_label(loop.spans[before],sources).first
            low=max(0.,(entry*curve.length-curve.lengths[pi])/piece.length)
            high=(distance-curve.lengths[pi])/piece.length
            if not 0.<=low<high<=1.:continue
            result.append(_Candidate(group,a,labels[0].index,curve,pi,poly,(low,high),
                (entry,labels[0].last),policies[0][0],b,labels[1].index,parent,circle,policies[1][0],
                bank,arm,(first,stop),*policy))
            if len(result)>maximum_pairs:raise ParentChainUnavailable('Parent-chain discovery exceeds its finite pair budget')
    return tuple(result)


def _regularity(candidate,radius):
    poly=candidate.polynomial;path=candidate.bank.path;arm=candidate.bank_arm
    first,last=candidate.bank_outward_interval
    low,high=sorted((first if arm.end==0 else path.length-first,last if arm.end==0 else path.length-last))
    offset=path.section.offset+.5*path.section.width*candidate.bank.side
    if path.exact:
        if np.any((path.stations>low)&(path.stations<high)):
            raise ParentChainUnavailable('A native parent-chain offset crosses an Exact miter')
        curvature=0.
    else:curvature=path.visible_curvature_bound(low,high)
    core=float(radius*poly.curvature_bound);native=float((abs(offset)+radius)*curvature)
    if core>=1.-1e-10 or native>=1.-1e-10:
        raise ParentChainUnavailable('A parent-chain normal-offset locus lacks a regularity certificate')
    return core,native


def bank_tangent(arm,station,interval):
    """Exact native offset direction, without differentiating translated XY."""
    from .routes import Route
    if type(arm.path) is Path:path=arm.path;direction=1 if arm.end==0 else -1
    elif type(arm.path) is Route and len(arm.path.pieces)==1:
        piece,=arm.path.pieces;path=piece.path;direction=1 if piece.forward else -1
    else:raise ParentChainUnavailable('Parent-chain bank needs one actual native Path')
    native=float(station if direction==1 else path.length-station)
    if path.exact:
        first,last=sorted(value if direction==1 else path.length-value for value in interval)
        if np.any((path.knots>first)&(path.knots<last)):
            raise ParentChainUnavailable('Parent-chain exact bank crosses an original miter')
        i=int(np.clip(np.searchsorted(path.knots,.5*(first+last),side='right')-1,0,len(path.chords)-1))
        offset=direction*(arm.offset+.5*arm.width)
        tangent=path.directions[i]+offset*(path.normals[i+1]-path.normals[i])/path.chords[i]
    else:tangent=path.frame(native)[1]
    length=float(np.linalg.norm(tangent))
    if length<=1e-10:raise ParentChainUnavailable('Parent-chain native offset has an undefined tangent')
    return direction*tangent/length


def offset_solution(polynomial,poly_interval,arm,bank_interval,radius,*,bank_first=False):
    """Bounded local two-offset solve; no extrapolation or unbounded search."""
    low,high=map(float,poly_interval);first,last=map(float,bank_interval)
    if (not 0.<=low<high<=1. or not 0.<=first<last<=arm.length
            or not math.isfinite(radius) or radius<=0.):raise ValueError('Invalid finite parent-chain contact domain')
    def pc(u):
        tangent=polynomial.derivative(u);tangent/=np.linalg.norm(tangent,axis=-1)[...,None]
        return polynomial(u)+radius*right(tangent)
    distance=arm.offset+.5*arm.width
    def bc(t):
        return arm.point(t,distance)+radius*right(bank_tangent(arm,float(t),bank_interval))
    us=np.linspace(low,high,33);ts=np.linspace(first,last,33)
    A=pc(us);B=np.array([bc(float(t)) for t in ts])
    i,j=np.unravel_index(np.argmin(np.linalg.norm(A[:,None]-B[None,:],axis=2)),(33,33))
    current=np.array((us[i],ts[j]))
    def derivative(function,t,a,b,step):
        x,y=max(a,t-step),min(b,t+step)
        return (function(y)-function(x))/(y-x)
    for iteration in range(32):
        delta=pc(current[0])-bc(current[1]);error=float(np.linalg.norm(delta))
        if error<1e-10:
            if not low<current[0]<high or not first<current[1]<last:break
            P=polynomial(current[0]);Q=arm.point(current[1],distance);C=pc(current[0])
            if bank_first:
                P, Q = Q, P
            start=math.atan2(*(P-C)[::-1]);stop=math.atan2(*(Q-C)[::-1]);sweep=-((start-stop)%math.tau)
            if not -math.pi<sweep<0.:break
            scale=max(1.,float(np.max(abs(C))),radius)
            trig=C+radius*np.array(((math.cos(start),math.sin(start)),(math.cos(start+sweep),math.sin(start+sweep))))
            if np.max(np.linalg.norm(trig-np.array((P,Q)),axis=1))>max(1e-10,32*np.spacing(scale)):break
            return dict(poly_contact=float(current[0]),bank_outward_contact=float(current[1]),
                center=tuple(C),first=tuple(P),last=tuple(Q),start_angle=start,sweep=sweep,
                iterations=iteration+1,residual=error)
        J=np.column_stack((derivative(pc,current[0],low,high,1e-5),
                           -derivative(bc,current[1],first,last,1e-3)))
        if abs(np.linalg.det(J))<1e-9:break
        step=np.linalg.solve(J,delta)
        for factor in (1.,.5,.25,.125,.0625):
            trial=np.clip(current-factor*step,(low,first),(high,last))
            if np.linalg.norm(pc(trial[0])-bc(trial[1]))<error:current=trial;break
        else:break
    raise ParentChainUnavailable('No finite contacts exist on the certified core and observed native bank')


def radius_candidates(radius,mode):
    """One shared, bounded policy for generation and source-rule validation."""
    if not math.isfinite(radius) or radius<=0. or mode not in ('STRICT','STYLE'):
        raise ValueError('Invalid original parent-chain radius policy')
    radii=(radius,) if mode=='STRICT' else tuple(CornerRadiusPolicy.candidates(radius,'STYLE'))
    if mode=='STYLE':
        radii=list(radii)
        while radii[-1]>.5:radii.append(max(.5,.8*radii[-1]))
    if len(radii)>17:raise ParentChainUnavailable('Parent-chain radius policy exceeds its finite allowance')
    return tuple(radii)


def fit(candidate):
    attempts=[]
    for radius in radius_candidates(candidate.radius,candidate.mode):
        attempts.append(radius)
        if len(attempts)>17:raise ParentChainUnavailable('Parent-chain radius search exceeds its finite allowance')
        try:
            regularity=_regularity(candidate,radius)
            solved=offset_solution(candidate.polynomial,candidate.polynomial_interval,candidate.bank_arm,
                                   candidate.bank_outward_interval,radius)
        except ParentChainUnavailable:continue
        return solved,radius,tuple(attempts),regularity
    raise ParentChainUnavailable('No allowed parent-chain radius meets its finite original STYLE/STRICT constraints')


def reconcile(planner,group):
    from .source_parent_chain_preparation import reconcile as prepare
    return prepare(planner,group)


def validate(group):
    """Recheck full arcs after any subsequent retained-circle preparation."""
    sources={source.region.owner:source for source in group.footprints}
    for area in group.footprints:
        record=getattr(area.layout,'source_parent_chain_fillet',None)
        if type(record) is not PreparedParentChainFillet:continue
        labels=[rim_label(span,sources) for loop in group.loops for span in loop.spans if span.owner==area.region.owner]
        if (not labels or any(label.kind!='curb' or label.index!=0 for label in labels)
                or abs(min(label.first for label in labels))>1e-12
                or abs(max(label.last for label in labels)-1.)>1e-12
                or abs(sum(label.last-label.first for label in labels)-1.)>1e-12):
            raise ParentChainUnavailable('Whole source union clips a parent-chain arc or exposes its internal closure')
    return group


def signature(group):
    record=getattr(group.physical_layout,'source_parent_chain_reconciliation',None)
    return (SOURCE_PARENT_CHAIN_POLICY,tuple(record['records'])) if record else None
