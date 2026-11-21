"""Failed-only native-curvature templates for an actual ordinary source root."""
from dataclasses import dataclass,replace
import hashlib
import math
import numpy as np

from .model import End,Key
from .mesh import cross
from .regions import RegionUnion
from .paths import right
from .continuations import WidthContinuation,ContinuationReachError,_evaluate
from .compound_root_envelopes import RootEnvelope
from .compound_dimensional_corners import DimensionalCompositeCurb,_CertifiedContinuation,_patch
from .free_pads import EnvelopeBoundary
from .generated_junction_curbs import GeneratedJunctionArea
from .source_root_curves import NativeCurvatureCore,ShiftedRootBell,bounded_bell
from .source_root_ports import OriginalRootPort,OriginalRootPortUnavailable,native_bank_curvature,refine_contact_port
from .source_root_witnesses import StagedRecipientWitness,StagedRecipientUnavailable,_group_token


SOURCE_ROOT_POLICY='2-native-cap-chains-and-endpoint-roots'


class OriginalRootTemplateUnavailable(ValueError):
    pass


@dataclass(frozen=True,eq=False)
class OriginalRootReference:
    network_uuid: str
    group: object
    site: object
    originals: tuple
    ports: tuple
    donor_owner: object
    source_objects: tuple
    group_token: tuple
    contact_bounds: tuple
    staged_witnesses: tuple

    @classmethod
    def build(cls,context,group,site,owner,donor,originals,maximum_stations,original_bounds):
        from .ordinary_singleton_sources import OriginalSingletonRootBounds
        from .routes import Route
        from .cycle_dimensional_retry import certify_layout
        originals=tuple(originals)
        if context.network.records().get(site.key)!=site:
            raise OriginalRootTemplateUnavailable('Native root template lost its actual original Junction')
        if (group is None or group.members!=(site.key,) or group.claims or any(loop.hole for loop in group.loops)
                or not 3<=len(originals)<=8 or set(site.ends)!={End(a.source.id,a.end) for a in originals}
                or len(set(site.ends))!=len(originals)
                or not any(c is donor for c in group.physical_layout.corners)
                or any(type(a.path) is not Route or a.profile.half or a.offset!=0.
                       or a.width!=a.source.width for a in originals)):
            raise OriginalRootTemplateUnavailable('Native root template needs one complete original full-profile junction exterior')
        actual={End(a.source.id,a.end):a for a in group.physical_layout.approaches}
        if set(actual)!=set(site.ends) or any(actual[End(a.source.id,a.end)] is not a for a in originals):
            raise OriginalRootTemplateUnavailable('Native root template changed the original approach inventory')
        if any(a.profile.datum!=originals[0].profile.datum for a in originals):
            raise OriginalRootTemplateUnavailable('Native root template needs its unchanged common original datum')
        if any(np.linalg.norm(a.origin-originals[0].origin)>1e-8
                or not np.array_equal(a.origin,a.frame(0.)[0])
                or not np.array_equal(a.direction,a.frame(0.)[1]) for a in originals):
            raise OriginalRootTemplateUnavailable('Native root template changed its actual common original origin or frame')
        own=tuple(End(a.source.id,a.end) for a in (donor.first,donor.second))
        requested,mode=site.corner_radius(own[0])
        preferred=donor.requested_radius if donor.requested_radius is not None else donor.radius
        ordered=sorted(originals,key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
        cyclic={(End(a.source.id,a.end),End(ordered[(i+1)%len(ordered)].source.id,ordered[(i+1)%len(ordered)].end))
            for i,a in enumerate(ordered)}
        if (own not in cyclic or donor.radius_mode!='STYLE' or mode!='STYLE'
                or preferred!=requested or donor.center is None or not math.pi*.5<donor.gap<math.pi*1.5):
            raise OriginalRootTemplateUnavailable('Native root template requires its actual half-plane STYLE donor')
        ports=[]
        for a in originals:
            values=[source for source in group.area_footprints if tuple(getattr(source.layout,'approaches',()))==(a,)]
            if len(values)!=1:raise OriginalRootTemplateUnavailable('Native root template lost an actual original source port')
            ports.append(OriginalRootPort.build(context,values[0],a))
        donors=[s for s in group.area_footprints if any(c is donor for c in getattr(s.layout,'corners',()))]
        if len(donors)!=1 or donors[0].region.owner!=owner:
            raise OriginalRootTemplateUnavailable('Native root template lost its actual replaced donor')
        port_owners={port.footprint.region.owner for port in ports}
        for source in group.area_footprints:
            if source.region.owner in port_owners:continue
            corners=getattr(source.layout,'corners',())
            if (getattr(source.layout,'original_junction',None)!=site.key or not corners
                    or any(any(actual.get(End(a.source.id,a.end)) is not a for a in (c.first,c.second)) for c in corners)):
                raise OriginalRootTemplateUnavailable('Native root template cannot borrow unowned physical support')
        staged=()
        if original_bounds is not None:
            if type(original_bounds) is not OriginalSingletonRootBounds or original_bounds.site!=site.key:
                raise OriginalRootTemplateUnavailable('Native root template has no typed original contact-plane request')
            expected=tuple(original_bounds.limits()[End(a.source.id,a.end)] for a in (donor.first,donor.second))
            if maximum_stations!=expected:
                raise OriginalRootTemplateUnavailable('Native root template changed an original donor bound')
            staged=StagedRecipientWitness.capture(context,group,original_bounds)
        elif maximum_stations is not None:
            raise OriginalRootTemplateUnavailable('Native root template cannot invent untyped contact bounds')
        by_end={port.endpoint:port for port in ports}
        limits=[]
        for i,a in enumerate((donor.first,donor.second)):
            end=End(a.source.id,a.end);port=by_end[end]
            value=(maximum_stations[i] if maximum_stations is not None else port.station-min(.25,.05*a.width))
            for witness in staged:
                if witness.endpoint==end:value=max(value,witness.contact)
            if not 0.<value<port.station:
                raise OriginalRootTemplateUnavailable('Native root template lost its finite profile lead')
            limits.append(float(value))
        certify_layout(context,group,())
        fixed=tuple(s.region for s in group.area_footprints if s.region.owner!=owner)
        # Invalid fixed operands must never be misclassified as bad candidates.
        RegionUnion(max_segments=20000).build(fixed)
        dependencies=tuple(sorted({Key(kind,p.section.id) for a in originals for p in a.path.pieces for kind in ('path','profile')}))
        return cls(context.network.uuid,group,site,originals,tuple(ports),owner,
            tuple((key,context[key]) for key in dependencies),_group_token(group),tuple(limits),staged)

    def validate(self,context):
        if (context.network.uuid!=self.network_uuid or _group_token(self.group)!=self.group_token
                or context.network.records().get(self.site.key)!=self.site
                or any(context[key] is not value for key,value in self.source_objects)
                or self.group.members!=(self.site.key,) or set(self.site.ends)!={port.endpoint for port in self.ports}):
            raise OriginalRootTemplateUnavailable('Native root reference is stale or belongs to another source graph')
        for port in self.ports:OriginalRootPort.build(context,port.footprint,port.approach)
        for witness in self.staged_witnesses:witness.validate_pre(context)

    def fixed_sources(self):
        return tuple(s for s in self.group.area_footprints if s.region.owner!=self.donor_owner)

    def contact_grid(self,context,port,stations,contact):
        """Use the actual circle station when a generated sample maps to that exact native cut."""
        from types import SimpleNamespace
        from .junction_bank_transitions import original_corner
        from .source_circular_intervals import PreparedCircularInterval
        arm=port.approach
        cut=arm.path.trims(contact)
        # Outward distance on a reversed Section can acquire a second floating
        # point representation after subtraction from its native length.
        i=int(np.searchsorted(stations,contact))
        duplicates=[]
        for direction in (-1,1):
            j=i+direction
            while 0<=j<len(stations) and arm.path.trims(float(stations[j]))==cut:
                duplicates.append(float(stations[j]))
                j+=direction
        if not duplicates:return stations,()
        self.validate(context)
        if not any(p.approach is arm and p.station==port.station and p.footprint.region.owner==port.footprint.region.owner
                   for p in self.ports):
            raise OriginalRootTemplateUnavailable('Coincident contact changed its original port')
        hard={0.,port.station}
        authenticated=False
        for source in self.group.area_footprints:
            for corner in source.layout.corners:
                contacts=[s for a,s in zip((corner.first,corner.second),corner.stations) if a is arm]
                hard.update(contacts)
                if (contact not in contacts or source.region.owner==self.donor_owner
                        or getattr(source.layout,'original_junction',None)!=self.site.key
                        or corner.center is None or corner.radius_mode!='STYLE'
                        or not any(c is corner for c in self.group.physical_layout.corners)
                        or not original_corner(SimpleNamespace(context=context),self.site,corner)):
                    continue
                for curve in source.layout.boundaries:
                    interval=getattr(curve,'original_circle_interval',None)
                    if type(interval) is PreparedCircularInterval and interval.original_corner is corner:
                        interval.validate(curve,corner)
                        authenticated=True
        if not authenticated:
            raise OriginalRootTemplateUnavailable('Coincident contact is missing its original retained circle')
        offsets=arm.offset+.5*arm.width*np.linspace(-1.,1.,5)
        for value in duplicates:
            route=arm.path
            if value in hard or not 0.<value<route.length:
                raise OriginalRootTemplateUnavailable('Coincident contact cannot replace a hard plane')
            piece_index=int(np.searchsorted(route.stations,value,side='right')-1)
            piece=route.pieces[piece_index]
            native=value-route.stations[piece_index]
            if not piece.forward:native=piece.path.length-native
            if (np.any(piece.path.knots==piece.path.parameter(native))
                    or not np.array_equal(arm.point(value,offsets),arm.point(contact,offsets))
                    or any(not np.array_equal(a,b) for a,b in zip(arm.frame(value),arm.frame(contact)))):
                raise OriginalRootTemplateUnavailable('Coincident contact cannot replace a source control, hard plane or distinct frame')
        return stations[~np.isin(stations,duplicates)],tuple(sorted(duplicates))

    def refine(self,context,contacts,source_vertices):
        fixed=[];reports=[]
        allocated=source_vertices+sum(len(s.region.xy) for s in self.fixed_sources())
        for source in self.fixed_sources():
            arms=getattr(source.layout,'approaches',())
            end=End(arms[0].source.id,arms[0].end) if len(arms)==1 else None
            if end in contacts:
                remaining=(20000-allocated)//2
                if remaining<3:raise OriginalRootTemplateUnavailable('Native root refinement exhausts its existing source allocation')
                previous=len(source.region.xy)
                source,details=refine_contact_port(context,source,contacts[end],maximum_samples=min(4096,remaining),root=self)
                allocated+=len(source.region.xy)-previous
                reports.append(details)
            fixed.append(source)
        if sum(len(s.region.xy) for s in fixed)>20000:
            raise OriginalRootTemplateUnavailable('Native root reference exceeds its existing source allocation')
        RegionUnion(max_segments=20000).build(s.region for s in fixed)
        return tuple(fixed),tuple(reports)

    def audit(self,candidate,source,fixed):
        by_end={port.endpoint:port for port in self.ports}
        for a,station,limit in zip((candidate.first,candidate.second),candidate.stations,self.contact_bounds):
            if station>limit or station+min(.25,.05*a.width)>by_end[End(a.source.id,a.end)].station:return False
        if len(source.region.xy)+sum(len(s.region.xy) for s in fixed)>20000:
            raise OriginalRootTemplateUnavailable('Native root audit exceeds its existing source allocation')
        try:loops=RegionUnion(max_segments=20000).build((source.region,*(s.region for s in fixed)))
        except ValueError as error:
            if str(error)!='Compound boundary has a point-only contact or unresolved numerical incidence':raise
            return False  # Fixed operands already passed their independent union.
        exposed=0.
        for loop in loops:
            for a in self.originals:
                cap=a.point(0.,np.array((-.5*a.width,.5*a.width)));axis=cap[1]-cap[0]
                length=float(np.linalg.norm(axis));delta=loop.xy-cap[0]
                close=abs(cross(delta,axis))/length<=1e-8;t=delta@axis/(axis@axis)
                low=np.maximum(0.,np.minimum(t,np.roll(t,-1)));high=np.minimum(1.,np.maximum(t,np.roll(t,-1)))
                if np.any(close&np.roll(close,-1)&((high-low)*length>1e-7)):return False
            for span in loop.spans:
                if loop.hole or span.owner!=source.region.owner:continue
                label=source.segments[span.segment]
                if label.kind=='curb':exposed+=abs((label.last-label.first)*(span.last-span.first))
        return abs(exposed-1.)<=1e-8 and sum(not loop.hole for loop in loops)==1 and not any(loop.hole for loop in loops)

    def signature(self,context):
        self.validate(context)
        return (self.site.key.text(),context.digest(self.site.key),tuple((a.source.id,a.end,a.width,a.offset,a.profile.rise,a.profile.datum)
            for a in self.originals),self.contact_bounds,
            tuple((key.text(),context.digest(key)) for key,_ in self.source_objects),
            tuple(w.signature(context) for w in self.staged_witnesses))


@dataclass(frozen=True)
class NativeRootEnvelope(RootEnvelope,GeneratedJunctionArea):
    reference: object
    contact_bounds: tuple
    staged_witnesses: tuple
    refinements: tuple


def prepare(context,group,owner,donor,settings,*,site,original_approaches,required_roots,
            maximum_stations=None,original_bounds=None):
    """One compulsory actual contact pair replaces an unusable coarse reach.

    The family spends one of the existing seventeen reach slots, with at most
    eight bell amplitudes and the established finite STYLE radius sequence.
    It is called only after the old source-root family has failed.
    """
    from .corner_radii import CornerRadiusPolicy
    reference=OriginalRootReference.build(context,group,site,owner,donor,original_approaches,
        maximum_stations,original_bounds)
    reference.validate(context)
    required_roots=tuple(required_roots)
    by_end={End(a.source.id,a.end):a for a in reference.originals}
    own={End(a.source.id,a.end) for a in (donor.first,donor.second)}
    if not 1<=len(required_roots)<=len(reference.originals)-2:
        raise OriginalRootTemplateUnavailable('Native root template needs a bounded actual exposed original root')
    for a,side in required_roots:
        end=End(a.source.id,a.end)
        if by_end.get(end) is not a or end in own or side not in (-1,1) or abs(float(a.profile.value(side)))>1e-10:
            raise OriginalRootTemplateUnavailable('Native root template changed its required original zero-height half')
        cap=a.point(0.,np.array((0.,side*.5*a.width)));axis=cap[1]-cap[0];length=float(np.linalg.norm(axis))
        exposed=False
        for loop in reference.group.loops:
            delta=loop.xy-cap[0];close=abs(cross(delta,axis))/length<=1e-8;t=delta@axis/(axis@axis)
            low=np.maximum(0.,np.minimum(t,np.roll(t,-1)));high=np.minimum(1.,np.maximum(t,np.roll(t,-1)))
            exposed|=bool(np.any(close&np.roll(close,-1)&((high-low)*length>1e-7)))
        if not exposed:raise OriginalRootTemplateUnavailable('Native root template cannot replace an already internal original half')
    stations=reference.contact_bounds
    contacts={End(a.source.id,a.end):s for a,s in zip((donor.first,donor.second),stations)}
    fixed=reference.fixed_sources();refined=None;refinement_reports=()
    preferred=donor.requested_radius if donor.requested_radius is not None else donor.radius
    radii=tuple(CornerRadiusPolicy.contact_candidates(donor.first,donor.second,preferred,'STYLE'))
    if len(radii)>17:raise OriginalRootTemplateUnavailable('Native root template radius allocation changed')
    width=min(donor.first.width,donor.second.width)
    support_width=max(a.width for a in reference.originals)
    required_width=max(a.width for a,side in required_roots)
    limit=max(.25*width,.5*max(0.,required_width-width))
    roots=tuple(a.point(0.,side*.5*a.width) for a,side in required_roots)+(
        donor.first.point(0.,-.5*donor.first.width),donor.second.point(0.,.5*donor.second.width))
    curvatures=(native_bank_curvature(donor.first,stations[0],-1,reverse=True),
                native_bank_curvature(donor.second,stations[1],1))
    attempts=0
    for radius in radii:
        first,va=donor.first.jet(stations[0],-.5*donor.first.width)
        last,vb=donor.second.jet(stations[1],.5*donor.second.width)
        try:
            old=WidthContinuation(first,last,-va,vb,radius,refine_curvature=True)
            base=NativeCurvatureCore(old,curvatures,radius)
        except ContinuationReachError:continue
        axis=(last-first)/np.linalg.norm(last-first);outward=right(axis);targets=[]
        for point in roots:
            projection=float((point-first)@axis)
            if not 0.<projection<float(np.linalg.norm(last-first)):continue
            low,high=0.,1.
            for _ in range(60):
                mid=.5*(low+high)
                if float((_evaluate(base.controls,mid)-first)@axis)<projection:low=mid
                else:high=mid
            t=.5*(low+high);targets.append((t,float((point-base(t))@outward)))
        for index,multiples in ((None,(0.,1.,2.,4.)),(6,(1.,4.)),(3,(1.,4.))):
            bell=lambda t:64.*t**3*(1.-t)**3 if index is None else float(ShiftedRootBell.bell(t,index))
            required=max((0.,*(distance/bell(t) for t,distance in targets)))
            amplitudes=sorted({0. if multiple==0. else required+multiple*settings.curve_tolerance
                for multiple in multiples if multiple==0. or required+multiple*settings.curve_tolerance<=limit})
            for amplitude in amplitudes:
                attempts+=1
                try:
                    envelope=bounded_bell(base,outward,amplitude,radius,index=index,maximum_distance=limit)
                    corner=replace(donor,stations=stations,center=None,radius=radius,continuation=True,
                        fit_policy='NATIVE_CURVATURE_ORIGINAL_ROOT_TEMPLATE')
                    curve=DimensionalCompositeCurb(corner,(_CertifiedContinuation(EnvelopeBoundary(corner,envelope)),),settings.curve_tolerance)
                    try:source=_patch(owner,corner,curve,stations,settings)
                    except ValueError as invalid:
                        if str(invalid) not in (
                                'Compound boundary has a point-only contact or unresolved numerical incidence',
                                'A physical footprint boundary intersects or touches itself',
                                'A physical footprint boundary overlaps or touches itself',
                                'A physical footprint boundary backtracks over itself'):
                            raise
                        continue  # Only this newly derived candidate is invalid.
                except ContinuationReachError:continue
                accepted=reference.audit(corner,source,fixed)
                if not accepted:
                    if refined is None:refined,refinement_reports=reference.refine(context,contacts,len(source.region.xy))
                    accepted=reference.audit(corner,source,refined)
                if not accepted:continue
                report=dict(policy='NATIVE_CURVATURE_ORIGINAL_ROOT_TEMPLATE',source_root_policy=SOURCE_ROOT_POLICY,
                    admission='EXPOSED_CYCLIC_ORIGINAL_ZERO_HEIGHT_ROOT',mode='STYLE',original_junction=site.key.text(),
                    first_section=donor.first.source.id,first_end=donor.first.end,
                    second_section=donor.second.source.id,second_end=donor.second.end,
                    preferred_minimum_radius=radius,requested_radius=preferred,curvature_bound=envelope.curvature_bound,
                    curvature_intervals=envelope.curvature_intervals,native_contact_curvatures=curvatures,
                    geometric_endpoint_curvature_matched=True,parametric_second_jets_prescribed=False,
                    base_displacement_convex_hull_bound=base.base_displacement_bound,
                    bell_maximum_displacement=amplitude,maximum_bell_displacement_m=limit,
                    total_displacement_from_zero_curvature_base_bound_m=base.base_displacement_bound+amplitude,
                    bell_degree=6 if index is None else 9,bell_index=index,
                    candidate_count=attempts,compulsory_actual_contact_candidate=True,
                    reach_candidates=1,maximum_reach_candidates=17,maximum_amplitudes_per_reach=8,
                    stations=stations,original_sources_unchanged=True,
                    complete_new_curb_exposed=True,all_original_root_caps_geometrically_internal=True,
                    native_port_contact_refinements=list(refinement_reports),
                    staged_recipient_baselines=[dict(w.proof) for w in reference.staged_witnesses])
                result=NativeRootEnvelope(source,curve,stations,report,corner,reference,stations,
                    reference.staged_witnesses,tuple(refinement_reports))
                source.layout.generated_junction_area=result
                return result
    raise ContinuationReachError('Native-curvature original root templates exhaust their unchanged finite allocation')


def refine_requested_port(context,source,records):
    """Normal source-cycle replay of the accepted template's proof sampling."""
    arms=getattr(source.layout,'approaches',())
    if len(arms)!=1:return source
    end=End(arms[0].source.id,arms[0].end);contacts={}
    for record in records:
        if not isinstance(record,NativeRootEnvelope):continue
        record.reference.validate(context)
        if record.refinements:
            for a,s in zip((record.corner.first,record.corner.second),record.stations):
                if End(a.source.id,a.end)==end:
                    # Shared knots can use the ordinary port grid. Replacing
                    # a generated station requires an unambiguous producer.
                    contacts[s]=record.reference if contacts.get(s,record.reference) is record.reference else None
    for contact,reference in sorted(contacts.items()):
        source,unused=refine_contact_port(context,source,contact,root=reference)
    from .bundle_crest_roots import CrestEnvelope
    for record in records:
        if isinstance(record,CrestEnvelope) and record.refined:
            record.reference.validate(context)
            for arm,station in zip((record.corner.first,record.corner.second),record.stations):
                if End(arm.source.id,arm.end)==end:
                    source=record.reference.refine(context,source,station)
    return source


def construction_signature(context,group):
    """Local source construction only; final neighbor inventory is separate."""
    records=[c.root_support_envelope for c in getattr(group.physical_layout,'corners',())
        if isinstance(getattr(c,'root_support_envelope',None),NativeRootEnvelope)]
    if not records:return None
    for record in records:
        record.reference.validate(context)
        if record.contact_bounds!=record.reference.contact_bounds or record.corner.stations!=record.stations:
            raise OriginalRootTemplateUnavailable('Constructed native root has forged contact bounds or planes')
        if getattr(record.source.layout,'generated_junction_area',None) is not record:
            raise OriginalRootTemplateUnavailable('Constructed native root lost its typed generated-area ownership')
        if (len(record.staged_witnesses)!=len(record.reference.staged_witnesses)
                or any(a is not b for a,b in zip(record.staged_witnesses,record.reference.staged_witnesses))):
            raise OriginalRootTemplateUnavailable('Constructed native root borrowed another staged recipient witness')
    return (SOURCE_ROOT_POLICY,tuple((r.reference.signature(context),r.stations,
        r.corner.radius,r.corner.requested_radius,r.corner.radius_mode,
        r.report['native_contact_curvatures'],r.report['bell_degree'],r.report['bell_index'],
        r.report['bell_maximum_displacement'],
        hashlib.sha256(np.asarray(r.curve.pieces[0].curve.continuation.controls).tobytes()).hexdigest())
        for r in records))
