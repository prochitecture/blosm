"""Fit a single extensive, tangent join when two exposed original arcs cannot meet."""
from dataclasses import dataclass
import math
from types import SimpleNamespace

import numpy as np

from .arc_lead_joins import publish_join_steps
from .bank_circle_contacts import ContactSampling, OriginalCornerArc
from .bank_circle_preparation import closure_steps, refine_parent
from .updates import drain_steps
from .compounds import rim_label
from .corner_radii import CornerRadiusPolicy
from .displaced_corner_authority import original_parent_corner
from .exposed_bank_joins import BankJoin, _bank_arm, _sectors, refine_bank_contact
from .model import Key
from .paths import right
from .retained_arc_fillets import CircularBoundaryGeometry, _intersection, phase
from .sampling import AnalyticalContactError
from .source_bank_reconciliation import NativeBankWitness, SourceBankDiscovery, visible_bank_stations
from .source_parent_chain_fillets import ParentChainUnavailable, offset_solution


@dataclass(frozen=True)
class CircleContact:
    arc: OriginalCornerArc
    interval: tuple
    parameter: float

    @property
    def owner(self):
        return self.arc.source_owner

    @property
    def index(self):
        return self.arc.source_curve_index

    def signature(self):
        return (self.owner.text(),self.index,self.arc.corner.radius,self.interval,self.parameter)

    def check_curve(self,curve):
        if (not isinstance(curve,ContactSampling)
                or curve.curve.original_prepared_curve is not self.arc.original_prepared_curve
                or curve.curve.corner is not self.arc.corner or curve.contact!=self.parameter):
            raise AnalyticalContactError('Rounded join lost its retained original circular contact')

    def profile_value(self,resolver):
        from .source_curb_conditions import _corner_slopes, SourceCurbConditionUnavailable
        source=resolver.sources[self.owner]
        node=next((n for n in resolver.context.network.intersections if n.key==source.layout.original_junction),None)
        corner=self.arc.corner
        if (node is None or node.key not in resolver.group.members or not original_parent_corner(resolver,node,corner)
                or not self.interval[0]<self.parameter<self.interval[1]):
            raise SourceCurbConditionUnavailable('Rounded join lost its retained parent profile')
        a,b=(arm.profile.datum+float(arm.profile.value(side))
             for arm,side in ((corner.first,-1),(corner.second,1)))
        if abs(a-b)>1e-10:
            raise SourceCurbConditionUnavailable('A graded original circle needs its along-curb derivative')
        u=self.parameter
        first,last=_corner_slopes(corner)
        slope=first+(last-first)*u**3*(10.+u*(-15.+6.*u))
        tangent=self.arc.derivative(u)
        return a,-slope*right(tangent/np.linalg.norm(tangent))


def _original_run(discovery,loop,at,step):
    """Follow one original circle from the failed crossing to its adjacent Section."""
    spans=loop.spans
    start=spans[at]
    source=discovery.sources[start.owner]
    label=rim_label(start,discovery.sources)
    if label.kind!='curb' or len(getattr(source.layout,'corners',()))!=1:
        return None
    c,=source.layout.corners
    curve=source.layout.boundaries[label.index]
    declaration=getattr(curve,'original_circle_interval',None)
    site=discovery.sites.get(getattr(source.layout,'original_junction',None))
    if (declaration is None or declaration.prepared_bounds!=(0.,1.)
            or declaration.circle_parameters[0]!=0. or c.radius_mode!='STYLE'
            or c.sweep>=0. or site is None or site.key not in discovery.group.members
            or not original_parent_corner(SimpleNamespace(context=discovery.context,group=discovery.group),site,c)):
        return None
    arc=OriginalCornerArc(start.owner,label.index,curve,c)
    for _ in range(len(spans)):
        adjacent=(at+step)%len(spans)
        span=spans[adjacent]
        if span.owner!=start.owner:
            if span.owner.kind!='cycle_strip':
                return None
            bank=discovery.describe(span,step==-1)
            if not isinstance(bank,NativeBankWitness):
                return None
            boundary=label.first if step==-1 else label.last
            tangent=arc.derivative(boundary)
            tangent=tangent/np.linalg.norm(tangent)
            if np.linalg.norm(tangent-discovery.tangent(bank))>1e-6:
                return None
            return arc,bank,at if step==-1 else adjacent,boundary
        other=rim_label(span,discovery.sources)
        if ((other.kind,other.index)!=('curb',label.index)
                or abs((label.first-other.last) if step==-1 else (other.first-label.last))>1e-10):
            return None
        label,at=other,adjacent
    return None


def _insert_steps(planner,group,discovery,kept,buried,bank,domain,bounds,solution,radius,sectors,circle_first):
    u=solution['poly_contact']
    native=solution['bank_outward_contact']
    if bank.side==-1:
        native=bank.path.length-native
    source=discovery.sources[kept.source_owner]
    tolerance=min(.00025,.05*planner.context.network.settings.curve_tolerance)
    length=radius*abs(solution['sweep'])
    count=max(2,math.ceil(length/min(.2,planner.context.network.settings.spacing,math.sqrt(8*tolerance*radius))))
    parent=SimpleNamespace(arc=kept,parent=source,index=kept.source_curve_index,phases=bounds)
    refined=refine_parent(parent,u,solution['sweep'],count,crossing=bounds[1 if circle_first else 0])
    strip,point,contact=refine_bank_contact(discovery.sources[bank.owner],bank,domain,native)
    contacts=[CircleContact(kept,bounds,u),contact]
    points=[kept(u),point]
    if not circle_first:
        contacts.reverse()
        points.reverse()
    areas=tuple(refined if a is source else strip if a.region.owner==bank.owner else a for a in group.footprints)
    ends=[(c.owner,c.index,c.parameter,p) for c,p in zip(contacts,points)]
    hidden=(buried.source_owner,buried.source_curve_index)
    internal,_=yield from closure_steps(areas,*ends,intermediate=(hidden,))
    number=len(getattr(group.physical_layout,'original_bank_joins',()))
    corner=CircularBoundaryGeometry(solution['center'],radius,solution['start_angle'],solution['sweep'])
    record=BankJoin(Key('original_bank_join',number),corner,tuple(points[0]),tuple(points[1]),tuple(contacts),
                    buried.source_owner,buried.corner,sectors,tuple(np.linspace(0.,1.,count+1)),tolerance)
    candidate=yield from publish_join_steps(planner,group,areas,record,internal)
    for loop in candidate.loops:
        for span in loop.spans:
            if span.owner==hidden[0] and rim_label(span,discovery.sources).index==hidden[1]:
                raise ValueError('The proposed join leaves part of the intervening curb exposed')
    return candidate


def recover(planner,group,switch):
    return drain_steps(recover_steps(planner,group,switch))


def recover_steps(planner,group,switch):
    """Try both finite arc-to-bank joins only after the ordinary arc pair exhausted STYLE radii."""
    li,at,_,_,_=switch
    discovery=SourceBankDiscovery(planner.context,group)
    loop=group.loops[li]
    runs=(_original_run(discovery,loop,(at-1)%len(loop.spans),-1),_original_run(discovery,loop,at,1))
    if any(run is None for run in runs):
        return None
    arcs,banks,indices,outer_phases=zip(*runs)
    if discovery.same_bank(*banks):
        return None
    sectors=_sectors(discovery,banks)
    if not sectors:
        return None
    crossing=_intersection(arcs[0].corner,arcs[1].corner,loop.xy[at])
    if crossing is None:
        return None
    labels=[rim_label(loop.spans[i],discovery.sources) for i in (at-1,at)]
    bounds=((outer_phases[0],phase(arcs[0].corner,crossing,labels[0].last)),
            (phase(arcs[1].corner,crossing,labels[1].first),outer_phases[1]))
    if any(not 0.<=lo<hi<=1. for lo,hi in bounds):
        return None
    radii=CornerRadiusPolicy.candidates(min(row[2] for row in sectors),'STYLE')
    visible_domains = {}
    for radius in radii:
        for i,arc in enumerate(arcs):
            j=1-i
            bank=banks[j]
            if j not in visible_domains:
                query=SimpleNamespace(loop=li,vertex=indices[j],first=bank,last=bank)
                visible_domains[j]=visible_bank_stations(discovery,query,incoming=j==0)
            domain=visible_domains[j]
            if domain is None or radius>=arc.corner.radius:
                continue
            frame=_bank_arm(planner.context,bank,domain)
            if frame is None:
                continue
            arm,outward,k=frame
            if (abs(arm.offset+.5*arm.width)+radius)*k>=1.-1e-10:
                continue
            try:
                solution=offset_solution(arc,bounds[i],arm,outward,radius,bank_first=i==1)
            except ParentChainUnavailable:
                continue
            try:
                return (yield from _insert_steps(planner,group,discovery,arc,arcs[j],bank,domain,bounds[i],
                               solution,radius,sectors,i==0))
            except ValueError:
                continue
    return None
