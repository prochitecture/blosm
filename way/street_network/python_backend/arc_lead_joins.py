"""Round an exposed original circle against a finite auxiliary corner lead."""
from copy import copy
from dataclasses import dataclass, replace
from types import SimpleNamespace
import math

import numpy as np

from .bank_circle_contacts import CircularJoin, ContactSampling, OriginalCornerArc
from .bank_circle_preparation import _union_steps, closure_steps, refine_parent
from .compound_ports import CompoundPortPlanner
from .compounds import RimSegment, SampledFootprint, rim_label
from .cycle_dimensional_retry import certify_layout_steps
from .updates import drain_steps
from .displaced_corner_authority import original_parent_corner
from .model import Key
from .paths import right
from .regions import Region, signed_area
from .retained_arc_fillets import CircularBoundaryGeometry, _mouth_points, _topology, auxiliary_candidates, phase
from .sampling import AnalyticalBoundarySampler, AnalyticalContactError
from .source_circular_intervals import PreparedCircularInterval
from .source_parent_chain_preparation import _future_guard, _refine_core


@dataclass(frozen=True)
class Parent:
    source: object
    index: int
    curve: object
    declaration: PreparedCircularInterval
    exposed: tuple
    lead: object = None

    @property
    def owner(self):
        return self.source.region.owner


def _parent(group, loop, sources, index, incoming):
    span = loop.spans[index]
    if span.owner.kind != 'cycle_corner':
        return None
    label = rim_label(span, sources)
    source = sources[span.owner]
    if label.kind != 'curb' or getattr(source.layout, 'original_junction', None) not in group.members:
        return None
    curve = source.layout.boundaries[label.index]
    declaration = getattr(curve, 'original_circle_interval', None)
    if not isinstance(declaration, PreparedCircularInterval) or isinstance(curve, ContactSampling):
        return None
    corner = declaration.original_corner
    if corner.sweep >= 0. or corner.radius_mode != 'STYLE' or tuple(source.layout.corners) != (corner,):
        return None
    cursor = label.last if incoming else label.first
    if declaration.prepared_bounds == (0., 1.) and declaration.circle_parameters[0] == 0.:
        lead, domain = None, (0., 1.)
    else:
        matches = [lead for lead in declaration.auxiliary_leads if lead.prepared_bounds[0] < cursor < lead.prepared_bounds[1]]
        if len(matches) != 1:
            return None
        lead, = matches
        domain = lead.prepared_bounds
    low, high = label.first, label.last
    step = -1 if incoming else 1
    for _ in range(len(loop.spans)-1):
        index = (index+step) % len(loop.spans)
        following = loop.spans[index]
        if following.owner != span.owner:
            break
        other = rim_label(following, sources)
        if (other.kind,other.index) != ('curb',label.index):
            break
        if incoming and abs(other.last-low) < 1e-10:
            low = other.first
        elif not incoming and abs(other.first-high) < 1e-10:
            high = other.last
        else:
            break
    return Parent(source,label.index,curve,declaration,(max(low,domain[0]),min(high,domain[1])),lead)


def fit(circle, bank, circle_first):
    """Intersect the inner offset circle with the directed offset lead and bound both contact parameters."""
    c = circle.declaration.original_corner
    points = np.asarray(bank.lead.reference_points)
    delta = points[-1]-points[0]
    length = float(np.linalg.norm(delta))
    if length <= 1e-8:
        return None
    direction = delta/length
    stations = (points-points[0]) @ direction
    if np.any(np.diff(stations) <= 0.) or np.max(abs((points-points[0]) @ right(direction))) > 1e-9:
        return None
    for radius,_ in auxiliary_candidates(c,bank.declaration.original_corner):
        if radius >= c.radius:
            continue
        origin = points[0]+radius*right(direction)
        offset = origin-c.center
        middle = -float(offset@direction)
        squared = (c.radius-radius)**2-float(np.linalg.norm(offset+middle*direction))**2
        if squared <= 0.:
            continue
        for distance in (middle-math.sqrt(squared),middle+math.sqrt(squared)):
            if not 0. < distance < length:
                continue
            center = origin+distance*direction
            p = c.center+c.radius*(center-c.center)/(c.radius-radius)
            q = points[0]+distance*direction
            u = phase(c,p,np.mean(circle.exposed))
            v = float(np.interp(distance,stations,bank.lead.prepared_parameters))
            if not (circle.exposed[0] < u < circle.exposed[1] and bank.exposed[0] < v < bank.exposed[1]):
                continue
            first,last = (p,q) if circle_first else (q,p)
            start = math.atan2(*(first-center)[::-1])
            sweep = (math.atan2(*(last-center)[::-1])-start) % -math.tau
            if -math.pi < sweep < 0.:
                return CircularBoundaryGeometry(tuple(center),radius,start,sweep),u,v
    return None


@dataclass(frozen=True)
class ArcLeadJoin(CircularJoin):
    owner: object
    corner: CircularBoundaryGeometry
    first: tuple
    last: tuple
    circle: Parent
    bank: Parent
    circle_first: bool
    circle_parameter: float
    bank_parameter: float
    source_reference_parameters: tuple
    source_chord_tolerance: float

    def signature(self):
        return (self.owner.text(),self.circle.owner.text(),self.circle.index,self.bank.owner.text(),
                self.bank.index,self.bank.lead.original_approach.source.id,self.bank.lead.bank_side,
                self.circle_first,self.circle_parameter,self.bank_parameter,self.corner.center,
                self.corner.radius,self.corner.sweep,self.source_reference_parameters)

    def contact(self, owner, index, curve):
        if (owner,index) == (self.circle.owner,self.circle.index):
            if not isinstance(curve,ContactSampling) or curve.curve.original_prepared_curve is not self.circle.curve:
                raise AnalyticalContactError('Rounded join lost its original circle')
            return float(not self.circle_first),self.circle_parameter
        if ((owner,index) == (self.bank.owner,self.bank.index)
                and getattr(curve,'original_circle_interval',None) is self.bank.declaration):
            return float(self.circle_first),self.bank_parameter
        raise AnalyticalContactError('Rounded join lost its declared auxiliary lead')

    def height_rule(self,resolver,first,last):
        from .source_curb_conditions import SourceCurbConditionUnavailable, _SourceRule, _corner_slopes
        from .source_parent_chain_conditions import endpoint_slopes
        for parent in (self.circle,self.bank):
            parent_source = resolver.sources[parent.owner]
            self.contact(parent.owner,parent.index,parent_source.layout.boundaries[parent.index])
            parent.declaration.validate(parent.curve,parent.declaration.original_corner)
            node = next((n for n in resolver.context.network.intersections if n.key==parent_source.layout.original_junction),None)
            if node is None or node.key not in resolver.group.members or not original_parent_corner(resolver,node,parent.declaration.original_corner):
                raise SourceCurbConditionUnavailable('Arc/lead contact needs both original road profiles')
        c = self.circle.declaration.original_corner
        levels = [arm.profile.datum+float(arm.profile.value(side)) for arm,side in ((c.first,-1),(c.second,1))]
        if abs(levels[0]-levels[1]) > 1e-10:
            raise SourceCurbConditionUnavailable('A graded original arc needs its along-curb derivative')
        u = self.circle_parameter
        a,b = _corner_slopes(c)
        k = a+(b-a)*u**3*(10.+u*(-15.+6.*u))
        angle = c.start_angle+c.sweep*u
        circle_tangent = np.array((math.sin(angle),-math.cos(angle)))
        lead = self.bank.lead
        arm,side = lead.original_approach,lead.bank_side
        bank_gradient = arm.physical_gradient(float(lead.source_station(self.bank_parameter)),side)
        gradients = [-k*right(circle_tangent),bank_gradient]
        heights = [levels[0],arm.profile.datum+float(arm.profile.value(side))]
        if not self.circle_first:
            gradients.reverse()
            heights.reverse()
        slopes = endpoint_slopes(self,gradients)
        rule = _SourceRule(self,first,last,*heights,*slopes,corner_fragment=True)
        rule.source_anchor = False
        return rule


def _insert_steps(planner,group,circle,bank,circle_first,solution,number):
    corner,u,v = solution
    settings = planner.context.network.settings
    tolerance = min(.00025,.05*settings.curve_tolerance)
    length = corner.radius*abs(corner.sweep)
    count = max(2,math.ceil(length/min(.2,settings.spacing,math.sqrt(8*tolerance*corner.radius))))
    arc = OriginalCornerArc(circle.owner,circle.index,circle.curve,circle.declaration.original_corner)
    crossing = circle.exposed[1 if circle_first else 0]
    refined = refine_parent(SimpleNamespace(arc=arc,parent=circle.source,index=circle.index,phases=circle.exposed),
                            u,corner.sweep,count,crossing=crossing)
    knots_at_contact = [i for i,label in enumerate(bank.source.segments)
                        if label.kind=='curb' and label.index==bank.index and label.first==v]
    if len(knots_at_contact)==1:
        lead_source,q = bank.source,bank.curve(v)
        if np.linalg.norm(lead_source.region.xy[knots_at_contact[0]]-q)>1e-10:
            raise ValueError('Original auxiliary lead knot no longer meets the source curb')
    else:
        lead_source,_,q = _refine_core(bank.source,bank.index,v)
    contacts = [(circle.owner,circle.index,u,arc(u)),(bank.owner,bank.index,v,q)]
    if not circle_first:
        contacts.reverse()
    areas = tuple(refined if s is circle.source else lead_source if s is bank.source else s for s in group.footprints)
    internal,_ = yield from closure_steps(areas,*contacts)
    key = Key('original_lead_fillet',number)
    knots = tuple(np.linspace(0.,1.,count+1))
    record = ArcLeadJoin(key,corner,tuple(contacts[0][-1]),tuple(contacts[1][-1]),circle,bank,circle_first,u,v,knots,tolerance)
    return (yield from publish_join_steps(planner, group, areas, record, internal))


def publish_join(planner, group, areas, record, internal):
    return drain_steps(publish_join_steps(planner, group, areas, record, internal))


def publish_join_steps(planner, group, areas, record, internal):
    """Prove the complete generated curb, ordered closure, original profiles and source coverage."""
    settings = planner.context.network.settings
    knots = record.source_reference_parameters
    key,tolerance = record.owner,record.source_chord_tolerance
    if any(s.region.owner==key for s in areas):
        raise ValueError('Circular join duplicates an existing physical area')
    polygon = np.vstack((record(knots)[:-1],internal.points[:-1]))
    if signed_area(polygon)<=1e-10 or len(polygon)+sum(len(s.region.xy) for s in areas)>20000:
        raise ValueError('Arc/lead wedge needs positive pavement area within the source vertex budget')
    labels = [RimSegment('curb',0,a,b) for a,b in zip(knots,knots[1:])]
    f = internal.stations/internal.length
    labels.extend(RimSegment('internal',1,a,b) for a,b in zip(f,f[1:]))
    layout = SimpleNamespace(approaches=(),outward_stations=(),mouths={},corners=(),rows=len(polygon),
                             boundaries=(record,internal),generated_junction_area=record)
    wedge = SampledFootprint(Region(key,polygon),tuple(labels),layout)
    areas = (*areas,wedge)
    loops = yield from _union_steps(areas)
    result = replace(group,footprints=areas,loops=loops,physical_layout=copy(group.physical_layout))
    if isinstance(record,ArcLeadJoin):
        result.physical_layout.original_lead_joins = (*getattr(group.physical_layout,'original_lead_joins',()),record)
    else:
        result.physical_layout.original_bank_joins = (*getattr(group.physical_layout,'original_bank_joins',()),record)
    consumed = {c.section_id for c in group.claims}
    # The union operates on the original source polygons. The incoming rim
    # can already have been resampled at a different set of curve parameters.
    original_loops = yield from _union_steps(group.footprints)
    original_areas = replace(group, loops=original_loops)
    _topology(original_areas,result,tolerance)
    yield from certify_layout_steps(planner.context,result,consumed)
    result = AnalyticalBoundarySampler(tolerance=settings.curve_tolerance,spacing=settings.spacing,max_samples=20000).build(result)
    _topology(original_areas,result,tolerance)
    result.physical_layout.source_cycle_full_width_proof = yield from certify_layout_steps(planner.context,result,consumed)
    before,after = _mouth_points(group),_mouth_points(result)
    if before.keys()!=after.keys() or any(not np.array_equal(p,after[k]) for k,p in before.items()) or CompoundPortPlanner.incomplete_ports(result):
        raise ValueError('Arc/lead join changed an original entrance')
    inventory = {s.region.owner:s for s in result.footprints}
    for owner in (k for k in inventory if k.kind in ('original_lead_fillet','original_bank_join')):
        covered = 0.
        segments = [rim_label(s,inventory) for loop in result.loops for s in loop.spans if s.owner==owner]
        for label in sorted(segments,key=lambda s:s.first):
            if label.kind!='curb' or abs(label.first-covered)>1e-10 or label.last<=label.first:
                raise ValueError('An arc/lead join was clipped by another source operand')
            covered=label.last
        if abs(covered-1.)>1e-10:
            raise ValueError('The compound boundary lost an arc/lead tangent contact')
    _future_guard(planner,result)
    return result


def failed_profile_join(context,group,error_type):
    """Return the actual rejected join only after the ordinary source-profile rules run."""
    from .compound_surfaces import CompoundSurfaceBuilder
    try:
        CompoundSurfaceBuilder(context.patterns,source_heights=True).boundary_conditions(group,context)
    except error_type as exc:
        return getattr(exc,'join',None)
    except ValueError:
        pass
    return None


def match_join(failure,loop,at,tolerance):
    return (failure[:2]==(loop.spans[at-1].owner,loop.spans[at].owner)
            and np.linalg.norm(np.asarray(failure[2])-loop.xy[at])<=tolerance)


def reconcile(planner,group):
    return drain_steps(reconcile_steps(planner,group))


def reconcile_steps(planner,group):
    """Repair only physical switches between a declared straight lead and a different original circle."""
    count=len(getattr(group.physical_layout,'original_lead_joins',()))
    for number in range(count,17):
        failure=None
        sources = {s.region.owner:s for s in group.footprints}
        changed = False
        for loop in group.loops:
            for i,s in enumerate(loop.spans):
                previous = loop.spans[i-1]
                if previous.owner==s.owner or previous.owner.kind!='cycle_corner' or s.owner.kind!='cycle_corner':
                    continue
                p = loop.xy[i]-loop.xy[i-1]
                q = loop.xy[(i+1)%len(loop.xy)]-loop.xy[i]
                if p[0]*q[1]-p[1]*q[0] >= -1e-5*np.linalg.norm(p)*np.linalg.norm(q):
                    continue
                a = _parent(group,loop,sources,(i-1)%len(loop.xy),True)
                b = _parent(group,loop,sources,i,False)
                if a is None or b is None or (a.lead is None)==(b.lead is None):
                    continue
                if failure is None:
                    from .source_curb_conditions import SourceAnchorConflict
                    failure=failed_profile_join(planner.context,group,SourceAnchorConflict)
                    if failure is None:
                        return group
                if not match_join(failure,loop,i,planner.context.network.settings.curve_tolerance):
                    continue
                circle_first = a.lead is None
                circle,bank = (a,b) if circle_first else (b,a)
                for parent in (circle,bank):
                    parent.declaration.validate(parent.curve,parent.declaration.original_corner)
                solution = fit(circle,bank,circle_first)
                if solution is None:
                    continue
                if number==16:
                    raise ValueError('Too many distinct exposed circle/lead switches in one junction')
                try:
                    group = yield from _insert_steps(planner,group,circle,bank,circle_first,solution,number)
                except ValueError:
                    # A speculative patch must not replace the original
                    # source area when a hole, contact or full-width check rejects it.
                    continue
                changed = True
                break
            if changed:
                break
        if not changed:
            return group
