"""Round two exposed Section banks separated by a clipped original style corner."""
from dataclasses import dataclass
from types import SimpleNamespace
import math

import numpy as np

from .arc_lead_joins import failed_profile_join, match_join, publish_join_steps
from .bank_circle_contacts import CircularJoin, OriginalCornerArc
from .bank_circle_preparation import _bank_knot, closure_steps
from .updates import drain_steps
from .compounds import rim_label
from .corner_radii import CornerRadiusPolicy
from .displaced_corner_authority import original_parent_corner
from .intersections import Approach
from .local_bank_contacts import DirectedBank
from .model import Key
from .retained_arc_fillets import CircularBoundaryGeometry
from .sampling import AnalyticalContactError
from .source_bank_reconciliation import SourceBankDiscovery, NativeBankWitness, visible_bank_stations, owned_ranges
from .source_curb_coordinates import CurbCoordinates
from .source_parent_chain_fillets import ParentChainUnavailable, offset_solution
from .source_parent_chain_preparation import _refine_strip


def _sectors(discovery, banks):
    """Read the style for the actual side of each original endpoint wherever it belongs."""
    records=[]
    for bank in banks:
        section=bank.path.section
        if (not bank.origins or section.scope_ends or section.selection_ends or
                discovery.context[Key('profile',section.id)].half):
            return ()
        for origin in bank.origins:
            end=origin.endpoint
            site=discovery.sites.get(origin.site)
            if (site is None or site.key not in discovery.group.members or end not in site.ends
                    or end.section_id!=section.id or bank.side*(1 if end.end==0 else -1)!=origin.outward_side):
                return ()
            def bearing(e):
                path=discovery.context[Key('path',e.section_id)]
                tangent=path.frame(0. if e.end==0 else path.length)[1]*(1 if e.end==0 else -1)
                return math.atan2(tangent[1],tangent[0])%math.tau
            ring=sorted(site.ends,key=bearing)
            sector=end if origin.outward_side==-1 else ring[ring.index(end)-1]
            radius,mode=site.corner_radius(sector)
            if mode!='STYLE' or not math.isfinite(radius) or radius<=0.:
                return ()
            records.append((site.key,sector,radius))
    return tuple(records)


@dataclass(frozen=True)
class BankContact:
    bank: NativeBankWitness
    station: float
    interval: tuple
    index: int
    parameter: float

    @property
    def owner(self):
        return self.bank.owner

    def signature(self):
        return (self.bank.path.section.id,self.bank.side,self.station,self.interval,self.parameter)

    def check_curve(self,curve):
        coordinates=getattr(curve,'source_coordinates',None)
        if (not isinstance(coordinates,CurbCoordinates) or coordinates.path is not self.bank.path
                or not np.any((coordinates.stations==self.station)&(coordinates.sides==self.bank.side))):
            raise AnalyticalContactError('A rounded bank join lost its original Section contact')

    def profile_value(self,resolver):
        from .source_curb_conditions import SourceCurbConditionUnavailable, _native_gradient
        bank=self.bank
        sid=bank.path.section.id
        if (resolver.context[Key('path',sid)] is not bank.path
                or owned_ranges(resolver.context,resolver.group).get(sid)!=bank.intervals
                or not self.interval[0]<self.station<self.interval[1]):
            raise SourceCurbConditionUnavailable('Rounded bank contact left its original owned interval')
        profile=resolver.context[Key('profile',sid)]
        return (profile.datum+float(profile.value(bank.side)),
                _native_gradient(bank.path,profile,self.station,bank.side,interval=self.interval))


@dataclass(frozen=True)
class BankJoin(CircularJoin):
    owner: object
    corner: CircularBoundaryGeometry
    first: tuple
    last: tuple
    contacts: tuple
    original_corner_owner: object
    original_corner: object
    sectors: tuple
    source_reference_parameters: tuple
    source_chord_tolerance: float

    def signature(self):
        return (self.owner.text(),self.original_corner_owner.text(),self.original_corner.radius,
                tuple(c.signature() for c in self.contacts),
                tuple((site.id,end.section_id,end.end,r) for site,end,r in self.sectors),
                self.corner.center,self.corner.radius,self.corner.sweep,self.source_reference_parameters)

    def contact(self,owner,index,curve):
        for i,c in enumerate(self.contacts):
            if (c.owner,c.index)==(owner,index):
                c.check_curve(curve)
                return float(i),c.parameter
        raise AnalyticalContactError('A rounded bank join meets an unrelated source')

    def height_rule(self,resolver,first,last):
        from .source_curb_conditions import _SourceRule, SourceCurbConditionUnavailable
        from .source_parent_chain_conditions import endpoint_slopes
        source=resolver.sources.get(self.original_corner_owner)
        sites={n.key:n for n in resolver.context.network.intersections if n.key in resolver.group.members}
        if (source is None or not any(c is self.original_corner for c in source.layout.corners)
                or not self.sectors or any(site not in sites or end not in sites[site].ends
                    or sites[site].corner_radius(end)!=(radius,'STYLE') for site,end,radius in self.sectors)
                or self.corner.radius not in CornerRadiusPolicy.candidates(min(s[2] for s in self.sectors),'STYLE')):
            raise SourceCurbConditionUnavailable('Rounded join lost its original corner or Section radius policy')
        levels,gradients=[],[]
        for c in self.contacts:
            self.contact(c.owner,c.index,resolver.sources[c.owner].layout.boundaries[c.index])
            z,gradient=c.profile_value(resolver)
            levels.append(z)
            gradients.append(gradient)
        slopes=endpoint_slopes(self,gradients)
        rule=_SourceRule(self,first,last,*levels,*slopes,corner_fragment=True)
        rule.source_anchor=False
        return rule


def _parent_runs(group):
    """Find single visible curb runs with an independent, finite source bank at either end."""
    for li,loop in enumerate(group.loops):
        for at,span in enumerate(loop.spans):
            if span.owner.kind!='cycle_corner' or loop.spans[at-1].owner.kind!='cycle_strip':
                continue
            after=(at+1)%len(loop.spans)
            while loop.spans[after].owner==span.owner:
                after=(after+1)%len(loop.spans)
            if loop.spans[after].owner.kind=='cycle_strip' and loop.spans[after].owner!=loop.spans[at-1].owner:
                yield li,at,after


def _case(discovery,li,at,after,*,check_tangency=False):
    group=discovery.group
    loop=group.loops[li];spans=loop.spans
    source=discovery.sources[spans[at].owner]
    labels=[rim_label(spans[i%len(spans)],discovery.sources)
            for i in range(at,at+(after-at)%len(spans))]
    first,last=labels[0],labels[-1]
    if (any((x.kind,x.index)!=('curb',first.index) for x in labels)
            or any(abs(a.last-b.first)>1e-10 for a,b in zip(labels,labels[1:]))
            or len(getattr(source.layout,'corners',()))!=1):
        return None
    c,=source.layout.corners
    curve=source.layout.boundaries[first.index]
    declaration=getattr(curve,'original_circle_interval',None)
    site=discovery.sites.get(getattr(source.layout,'original_junction',None))
    if (c.center is None or c.sweep>=0. or c.radius_mode!='STYLE' or site is None or site.key not in group.members
            or declaration is None
            or declaration.circle_parameters[0]!=0. or c.first.profile.half or c.second.profile.half):
        return None
    bounds=first.first,last.last
    low,high=declaration.prepared_bounds
    if bounds==(low,high) or not low<=bounds[0]<bounds[1]<=high:
        return None
    banks=discovery.describe(spans[at-1],True),discovery.describe(spans[after],False)
    if any(not isinstance(b,NativeBankWitness) for b in banks):
        return None
    angles=c.start_angle+c.sweep*declaration.circle_parameter(np.asarray(bounds))
    tangents=np.stack((np.sin(angles),-np.cos(angles)),axis=1)
    pairs=((discovery.tangent(banks[0]),tangents[0]),(tangents[1],discovery.tangent(banks[1])))
    turns=[math.atan2(a[0]*b[1]-a[1]*b[0],float(a@b)) for a,b in pairs]
    # The old chain's inward crossfall would point out of the pavement at
    # this switch. Tangencies elsewhere on the chain must survive.
    threshold=-math.radians(1e-4) if check_tangency else -.5*math.pi
    if min(turns)>=threshold or max(turns)>1e-6:
        return None
    if not original_parent_corner(SimpleNamespace(context=discovery.context,group=group),site,c):
        return None
    sections=_sectors(discovery,banks)
    if not sections:
        return None
    intervals=[]
    for i,bank in enumerate(banks):
        query=SimpleNamespace(loop=li,vertex=(at,after)[i],first=bank,last=bank)
        visible=visible_bank_stations(discovery,query,incoming=i==0)
        if visible is None:
            return None
        intervals.append(visible)
    arc=OriginalCornerArc(source.region.owner,first.index,curve,c)
    return banks,tuple(intervals),source,first.index,arc,bounds,turns,sections


def _bank_arm(context,bank,domain):
    path=bank.path
    low,high=sorted(domain)
    if path.exact and np.any((low<path.knots)&(path.knots<high)):
        return None
    curvature=0. if path.exact else path.visible_curvature_bound(low,high)
    arm=Approach(path,int(bank.side==-1),context[Key('profile',path.section.id)])
    outward=tuple(domain if bank.side==1 else (path.length-v for v in domain))
    return arm,outward,curvature


def refine_bank_contact(source,bank,domain,station):
    """Insert a native station without moving either original Section end cross-section."""
    strip,point,_=_refine_strip(source,bank.path,station,bank.side)
    index,parameter=_bank_knot(strip,point,bank.side)
    for i,end in getattr(source.layout,'original_root_caps',{}).items():
        matches=[j for j,e in strip.layout.original_root_caps.items() if e==end]
        if len(matches)!=1 or not np.array_equal(source.region.xy[[i,(i+1)%len(source.region.xy)]],
                strip.region.xy[[matches[0],(matches[0]+1)%len(strip.region.xy)]]):
            raise ValueError('Adding a bank contact moved a Section end cross-section')
    return strip,point,BankContact(bank,station,tuple(sorted(domain)),index,parameter)


def _publish_steps(planner,group,case,solution,radius,number):
    banks,domains,source,index,arc,bounds,turns,sectors=case
    context=planner.context
    native=(float(DirectedBank(banks[0],domains[0]).station(solution['poly_contact'])),
            float(solution['bank_outward_contact'] if banks[1].side==1 else
                  banks[1].path.length-solution['bank_outward_contact']))
    sources={s.region.owner:s for s in group.footprints}
    cuts=[];contacts=[]
    for bank,domain,station in zip(banks,domains,native):
        strip,p,contact=refine_bank_contact(sources[bank.owner],bank,domain,station)
        sources[bank.owner]=strip
        cuts.append((bank.owner,contact.index,contact.parameter,p))
        contacts.append(contact)
    areas=tuple(sources[s.region.owner] for s in group.footprints)
    internal,_=yield from closure_steps(areas,*cuts,intermediate=((source.region.owner,index),))
    corner=CircularBoundaryGeometry(solution['center'],radius,solution['start_angle'],solution['sweep'])
    tol=min(.00025,.05*context.network.settings.curve_tolerance)
    count=max(2,math.ceil(radius*abs(corner.sweep)/min(.2,context.network.settings.spacing,math.sqrt(8*tol*radius))))
    record=BankJoin(Key('original_bank_join',number),corner,tuple(cuts[0][-1]),tuple(cuts[1][-1]),
                    tuple(contacts),source.region.owner,arc.corner,sectors,tuple(np.linspace(0.,1.,count+1)),tol)
    candidate=yield from publish_join_steps(planner,group,areas,record,internal)
    if any(s.owner==source.region.owner and rim_label(s,sources).index==index for loop in candidate.loops for s in loop.spans):
        raise ValueError('The original corner is only partly covered by the new bank join')
    return candidate


def reconcile(planner,group,*,check_tangency=False):
    return drain_steps(reconcile_steps(planner,group,check_tangency=check_tangency))


def reconcile_steps(planner,group,*,check_tangency=False):
    """Repair failed profile contacts, or remaining clockwise turns after an arc-pair recovery."""
    for number in range(len(getattr(group.physical_layout,'original_bank_joins',())),17):
        discovery=None
        repaired=None
        failure=None
        unresolved=False
        for li,at,after in _parent_runs(group):
            if discovery is None:
                discovery=SourceBankDiscovery(planner.context,group)
            case=_case(discovery,li,at,after,check_tangency=check_tangency)
            if case is None:
                continue
            unresolved=True
            if not check_tangency:
                if failure is None:
                    from .source_curb_conditions import SourceJetDirectionError
                    failure=failed_profile_join(planner.context,group,SourceJetDirectionError)
                    if failure is None:
                        return group
                tolerance=planner.context.network.settings.curve_tolerance
                if not any(match_join(failure,group.loops[li],i,tolerance) for i in (at,after)):
                    continue
            banks,domains,source,index,arc,bounds,turns,sectors=case
            arms=[_bank_arm(planner.context,b,d) for b,d in zip(banks,domains)]
            if None in arms:
                continue
            radii=tuple(CornerRadiusPolicy.candidates(min(row[2] for row in sectors),'STYLE'))
            # Prefer a join retaining a visible piece of the original arc.
            # Leave such a case with the original-circle constructor.
            retained=False
            for i,(arm,outward,k) in enumerate(arms):
                if turns[i]>=-1e-6:
                    continue
                for r in radii:
                    if r>=arc.corner.radius or (abs(arm.offset+.5*arm.width)+r)*k>=1.-1e-10:
                        continue
                    try:
                        offset_solution(arc,bounds,arm,outward,r,bank_first=i==0)
                    except ParentChainUnavailable:
                        continue
                    retained=True
                    break
            if retained:
                continue
            if number==16:
                raise ValueError('Rounded bank joins exceed the finite junction contact budget')
            first=DirectedBank(banks[0],domains[0])
            arm,outward,_=arms[1]
            for radius in radii:
                if any((abs(a.offset+.5*a.width)+radius)*k>=1.-1e-10 for a,_,k in arms):
                    continue
                try:
                    solution=offset_solution(first,(0.,1.),arm,outward,radius)
                except ParentChainUnavailable:
                    continue
                try:
                    repaired=yield from _publish_steps(planner,group,case,solution,radius,number)
                except ValueError:
                    continue
                break
            if repaired is not None:
                break
        if repaired is None:
            if check_tangency and unresolved:
                raise ValueError('A clipped original circle still meets its Section bank at a clockwise turn')
            return group
        group=repaired
