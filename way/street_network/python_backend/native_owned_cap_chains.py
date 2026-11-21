"""Native C2 repair of an authenticated owned root-cap chain.

This is a separate authority from an adjacent original corner or a cyclic
root donor. It can cover two exposed full-profile caps, or the zero-height
outer banks of consecutive sections of a divided carriageway. The finite
bank chain between them and every original area operand are retained.
"""
from copy import copy
from dataclasses import dataclass,replace
import hashlib,math
import numpy as np
from .model import End,Key,fingerprint
from .compounds import rim_label
from .generated_junction_curbs import GeneratedJunctionArea
from .source_bank_reconciliation import CircularCornerWitness,SourceBankDiscovery,NativeBankWitness,fit_switch
from .source_root_ports import NativeRolePort,OriginalRootPort,native_bank_curvature
from .source_curb_rounding import FiniteSourceBank,fit_banks
from .compound_dimensional_corners import DimensionalCompositeCurb,_CertifiedContinuation,_patch
from .continuations import WidthContinuation,ContinuationReachError,_evaluate
from .source_root_curves import NativeCurvatureCore,ShiftedRootBell,bounded_bell
from .free_pads import EnvelopeBoundary
from .paths import right
from .regions import RegionUnion
from .corner_radii import CornerRadiusPolicy
from .corner_radii import CornerContactError
from .compound_ports import CompoundPortPlanner
from .cycle_dimensional_retry import certify_layout
from .native_c2_polygon_exposure import certify as certify_polynomial_exposure,NativeC2ExposureError

POLICY='1-authenticated-finite-owned-native-cap-chain-C2'


class OwnedCapChainAuthorityError(RuntimeError):
    """Original native authority changed; never treat this as a search miss."""


class OwnedCapChainCandidateError(ValueError):
    """One finite geometric candidate did not satisfy its unchanged gates."""


def _path_token(path):
    return (fingerprint(path.section),path.length,tuple((name,value.dtype.str,value.shape,
        hashlib.sha256(value.tobytes()).hexdigest()) for name,value in sorted(vars(path).items())
        if isinstance(value,np.ndarray)))


def _input_token(context,group):
    from .ordinary_mouth_clearance import _authenticate_native_half_profiles
    try:authenticated=_authenticate_native_half_profiles(context,group.physical_layout)
    except ValueError as error:
        raise OwnedCapChainAuthorityError('Owned cap chain lost canonical native source authority') from error
    if not authenticated:
        raise OwnedCapChainAuthorityError('Owned cap chain needs canonical native source roles')
    arms=group.physical_layout.approaches
    return (fingerprint(context.network.settings),tuple(group.members),tuple(group.claims),
        tuple((End(a.source.id,a.end),a.width,a.offset,fingerprint(a.profile),
            tuple((p.forward,_path_token(p.path)) for p in a.path.pieces)) for a in arms),
        tuple(group.physical_layout.outward_stations),tuple(sorted(group.physical_layout.mouths.items())))


@dataclass(frozen=True)
class OwnedNativeCapChainReference:
    original_group: object
    first: object
    last: object
    intermediate: tuple
    original_ports: tuple
    sectors: tuple
    radius: float
    mode: str
    contacts: tuple
    ledger: tuple
    input_token: tuple
    original_geometry: tuple

    def validate(self,context):
        if _input_token(context,self.original_group)!=self.input_token:
            raise OwnedCapChainAuthorityError('Owned cap chain changed its original Path/Profile/plane authority')
        port_type = NativeRolePort if self.first.profile.half else OriginalRootPort
        for reference in self.original_ports:
            if port_type.build(context,reference.footprint,reference.approach)!=reference:
                raise OwnedCapChainAuthorityError('Owned cap chain changed an original complete port')
        actual=tuple((s.region.owner,s.region.xy.tobytes(),s.segments)
            for s in self.original_group.area_footprints)
        if actual!=self.original_geometry:
            raise OwnedCapChainAuthorityError('Owned cap chain changed an original source operand')
        # Re-derive semantic chain, native side/ranges and finite ledger limits
        # from actual ports.  A typed record alone cannot authorize forged cuts.
        expected=_references(context,self.original_group)
        if len(expected)!=1:
            raise OwnedCapChainAuthorityError('Owned cap chain lost its unique actual boundary chain')
        other,=expected
        witness_token=lambda w:(id(w.path),w.side,w.station,w.intervals,w.origins,w.owner,w.index,w.proof)
        if (self.first is not other.first or self.last is not other.last
                or self.contacts!=other.contacts or self.ledger!=other.ledger
                or self.sectors!=other.sectors or self.radius!=other.radius or self.mode!=other.mode
                or tuple(map(witness_token,self.intermediate))!=tuple(map(witness_token,other.intermediate))):
            raise OwnedCapChainAuthorityError('Owned cap chain changed its actual native ranges/contact/sector authority')
        return self


class ContactPhasedNativeC2(DimensionalCompositeCurb):
    """Same analytical curve, certified cadence aligned to a real native grid."""
    def __init__(self,original,parameters,error):
        self.__dict__.update(original.__dict__)
        self.contact_phase_parameters=np.asarray(parameters,float)
        self.contact_phase_error_bound=float(error)
        if (self.contact_phase_parameters[0]!=0. or self.contact_phase_parameters[-1]!=1.
                or np.any(np.diff(self.contact_phase_parameters)<=0.)):
            raise OwnedCapChainCandidateError('Native C2 contact phase is not a complete partition')
        self.contact_phase_parameters.flags.writeable=False

    def certified_parameters(self,sampler,first,last):
        if self.contact_phase_error_bound>sampler.tolerance:
            raise OwnedCapChainCandidateError('Native C2 contact cadence needs a tighter chord certificate')
        values=np.r_[first,self.contact_phase_parameters[(self.contact_phase_parameters>first)
            &(self.contact_phase_parameters<last)],last]
        if len(values)>sampler.max_samples:
            raise OwnedCapChainCandidateError('Native C2 cadence exhausts the unchanged sample budget')
        return values


@dataclass(frozen=True)
class PreparedOwnedNativeCapChain(GeneratedJunctionArea):
    reference: OwnedNativeCapChainReference
    curve: object
    control_sha256: str
    report: dict


@dataclass(frozen=True)
class OwnedCapChainCarrier:
    """Derived carrier for one actual cap chain, never an original incidence."""
    first: object
    second: object
    stations: tuple
    radius: float
    requested_radius: float
    radius_mode: str
    radius_attempts: tuple
    owned_cap_chain_reference: object
    center: object = None
    continuation: bool = True
    median_cap: bool = False
    fit_policy: str = 'OWNED_NATIVE_CAP_CHAIN_C2'


def _port(group,arm):
    matches=[s for s in group.area_footprints if tuple(getattr(s.layout,'approaches',()))==(arm,)]
    if len(matches)!=1:raise OwnedCapChainAuthorityError('Owned cap chain lost its unique original source port')
    return matches[0]


def _references(context,group):
    if (len(group.members)!=1 or group.claims or len(group.loops)!=1 or group.loops[0].hole
            or not 3<=len(group.physical_layout.approaches)<=8):return ()
    token=_input_token(context,group);discovery=SourceBankDiscovery(context,group)
    loop=group.loops[0];inventory=discovery.sources;spans=loop.spans;n=len(spans)
    caps=[i for i,s in enumerate(spans) if rim_label(s,inventory).kind=='internal']
    # A root-cap chain has exactly two actual exposed transverse source caps.
    if len(caps)!=2:return ()
    site=discovery.sites[group.members[0]];arms={End(a.source.id,a.end):a for a in group.physical_layout.approaches}
    results=[]
    for start,end in (tuple(caps),tuple(reversed(caps))):
        walk=[];at=(start+1)%n
        while at!=end and len(walk)<=n:walk.append(at);at=(at+1)%n
        if not walk or len(walk)>n:continue
        first=discovery.describe(spans[start-1],True);last=discovery.describe(spans[(end+1)%n],False)
        if not isinstance(first,NativeBankWitness) or not isinstance(last,NativeBankWitness):continue
        if len(first.origins)!=1 or len(last.origins)!=1:continue
        aorigin,borigin=first.origins[0],last.origins[0]
        if (aorigin.site!=site.key or borigin.site!=site.key
                or aorigin.outward_side!=-1 or borigin.outward_side!=1):continue
        a,b=arms.get(aorigin.endpoint),arms.get(borigin.endpoint)
        if a is None or b is None or a is b:continue
        bundle_banks = bool(a.profile.half or b.profile.half)
        if bundle_banks:
            from .local_width_continuation import _continuous_separator_crest
            if not (a.profile.half and b.profile.half and _continuous_separator_crest(context,a,b)):continue
            # For regular Smooth offsets the original frame gives the analytic
            # tangent directly. Exact miters need a different derivative.
            if any(piece.path.exact for arm in (a,b) for piece in arm.path.pieces):continue
        if a.profile.datum!=b.profile.datum or a.profile.value(-1)!=0. or b.profile.value(1)!=0.:continue
        port_a,port_b=_port(group,a),_port(group,b)
        if (getattr(port_a.layout,'original_root_caps',{}).get(spans[start].segment)!=aorigin.endpoint
                or getattr(port_b.layout,'original_root_caps',{}).get(spans[end].segment)!=borigin.endpoint):continue
        port_type = NativeRolePort if bundle_banks else OriginalRootPort
        refs=(port_type.build(context,port_a,a),port_type.build(context,port_b,b))
        intermediate=[];valid=True
        for i in walk:
            for tail in (False,True):
                witness=discovery.describe(spans[i],tail)
                banks = (witness,)
                if bundle_banks and isinstance(witness,CircularCornerWitness):
                    if witness.site!=site.key or getattr(inventory[witness.owner].layout,'source_bank_fillet',None):
                        valid=False;break
                    banks=(witness.first,witness.last)
                if any(not isinstance(w,NativeBankWitness) or not w.origins or any(
                        o.site!=site.key or o.endpoint not in arms for o in w.origins) for w in banks):
                    valid=False;break
                intermediate.extend(banks)
            if not valid:break
        if not valid:continue
        middle_ends={o.endpoint for w in intermediate for o in w.origins}
        if not middle_ends or middle_ends.intersection((aorigin.endpoint,borigin.endpoint)):continue
        def angle(endpoint):
            arm=arms[endpoint];v=arm.frame(0.)[1];return math.atan2(v[1],v[0])%math.tau
        ring=tuple(sorted(site.ends,key=angle));i=ring.index(aorigin.endpoint);j=ring.index(borigin.endpoint)
        sectors=[];k=i
        while k!=j:
            sectors.append(ring[k]);k=(k+1)%len(ring)
        if set(sectors[1:])!=middle_ends:continue
        authored=tuple((endpoint,*site.corner_radius(endpoint)) for endpoint in sectors)
        strict={radius for _,radius,mode in authored if mode=='STRICT'}
        style=[radius for _,radius,mode in authored if mode=='STYLE']
        if len(strict)>1:continue
        radius=next(iter(strict)) if strict else min(style)
        mode='STRICT' if strict else 'STYLE'
        if style and radius>min(style):continue
        try:
            fit_banks(context,FiniteSourceBank(first.path,first.side,first.intervals),
                FiniteSourceBank(last.path,last.side,last.intervals),radius,mode,
                .5*(a.point(0.,-.5*a.width)+b.point(0.,.5*b.width)))
        except CornerContactError:pass
        else:continue  # Existing native circle remains the preferred construction.
        # Keep the real mouth reserve.  A genuine subsequent native switch
        # additionally bounds the last contact before its existing circle.
        limits=[refs[0].station-min(.25,.05*a.width),refs[1].station-min(.25,.05*b.width)]
        for proposal in discovery.discover():
            if (proposal.first.path is last.path and proposal.first.side==last.side
                    and any(o.endpoint==borigin.endpoint for o in proposal.first.origins)):
                fitted=fit_switch(context,proposal)
                arm=fitted.corner.first
                if arm.source.id!=b.source.id:continue
                contact=fitted.corner.stations[0]
                if arm.end!=b.end:contact=b.length-contact
                limits[1]=min(limits[1],contact-min(.25,.05*b.width))
        contacts=[];ledger=[]
        for arm,side,port,limit in ((a,-1,port_a,limits[0]),(b,1,port_b,limits[1])):
            bank=port.layout.boundaries[0 if side==1 else 1]
            values=bank.source_coordinates.stations
            allowed=values[(values>0.)&(values<=limit)]
            if not len(allowed):valid=False;break
            # A short cross street can intersect a Bundle bank close to the
            # root even if its mouth is far away. Use a local original bank
            # vertex at the scale of the declared radius for this join.
            station=float(min(allowed,key=lambda s:abs(s-radius)) if bundle_banks else max(allowed))
            point=arm.point(station,side*.5*arm.width)
            if not np.any(np.all(bank.points==point,axis=1)):
                raise OwnedCapChainAuthorityError('Owned cap contact is not an exact original bank vertex')
            contacts.append(station);ledger.append((arm.source.id,arm.end,side,limit,station,tuple(point)))
        if not valid:continue
        geometry=tuple((s.region.owner,s.region.xy.tobytes(),s.segments) for s in group.area_footprints)
        results.append(OwnedNativeCapChainReference(group,a,b,tuple(intermediate),refs,authored,
            radius,mode,tuple(contacts),tuple(ledger),token,geometry))
    if len(results)>1:raise OwnedCapChainAuthorityError('Owned cap chain has ambiguous original boundary authority')
    return tuple(results)


def _split(controls,t):
    work=np.asarray(controls).copy();left=[work[0].copy()];right=[work[-1].copy()]
    for _ in range(len(work)-1):
        work=(1.-t)*work[:-1]+t*work[1:];left.append(work[0].copy());right.append(work[-1].copy())
    return np.asarray(left),np.asarray(right[::-1])


def _restrict(controls,a,b):
    if b<1.:controls=_split(controls,b)[0]
    if a>0.:controls=_split(controls,a/b)[1]
    return controls


def _chord_bound(controls,parameters):
    maximum=0.
    for a,b in zip(parameters,parameters[1:]):
        values=_restrict(controls,float(a),float(b));n=len(values)-1
        line=values[0]+np.arange(n+1)[:,None]/n*(values[-1]-values[0])
        maximum=max(maximum,float(np.linalg.norm(values-line,axis=1).max()))
    return maximum


def _contact_phase(reference,curve,source,settings):
    core=curve.pieces[0].curve.continuation;last=reference.last;contact=reference.contacts[1]
    bank=_port(reference.original_group,last).layout.boundaries[0];grid=bank.source_coordinates.stations
    old=np.asarray(sorted({0.,1.,*(v for s in source.segments if s.kind=='curb' for v in (s.first,s.last))}))
    if contact < min(.5,.05*last.width):
        error=_chord_bound(core.controls,old)
        if error>settings.curve_tolerance:
            raise OwnedCapChainCandidateError('Native cap-chain original cadence exceeds the curve tolerance')
        return ContactPhasedNativeC2(curve,old,error),dict(original_parameters=len(old),
            phased_parameters=len(old),previous_chord_bound=error,phased_chord_bound=error,
            analytical_controls_unchanged=True)
    tangent=core.first_derivative[-1]/np.linalg.norm(core.first_derivative[-1]);endpoint=core.controls[-1]
    begin=max(i for i,s in enumerate(grid) if s<=contact-min(.5,.05*last.width))
    end=next(i for i,s in enumerate(grid) if s==contact);points=bank.points[begin:end+1]
    projection=(points-endpoint)@tangent
    if not np.all(np.diff(projection)>0.) or abs(projection[-1])>1e-12:
        raise OwnedCapChainCandidateError('Native cap-chain contact has no monotone common grid')
    mapped=[]
    for wanted in projection:
        if abs(wanted)<=1e-12:t=1.
        else:
            low,high=0.,1.
            for _ in range(60):
                middle=.5*(low+high)
                if float((curve(middle)-endpoint)@tangent)<wanted:low=middle
                else:high=middle
            t=.5*(low+high)
        mapped.append(t)
    mapped=np.asarray(mapped)
    if (_restrict(core.first_derivative,float(mapped[0]),1.)@tangent).min()<=1e-9:
        raise OwnedCapChainCandidateError('Native cap-chain generated phase loses monotonicity')
    parameters=np.unique(np.r_[old[old<mapped[0]],mapped])
    before,after=_chord_bound(core.controls,old),_chord_bound(core.controls,parameters)
    if after>before*(1.+1e-12) or after>settings.curve_tolerance:
        raise OwnedCapChainCandidateError('Native cap-chain generated phase worsens its chord certificate')
    return ContactPhasedNativeC2(curve,parameters,after),dict(original_parameters=len(old),
        phased_parameters=len(parameters),phase_start=float(mapped[0]),native_grid_indices=(begin,end),
        previous_chord_bound=before,phased_chord_bound=after,analytical_controls_unchanged=True)


def prepare(planner,group):
    context=planner.context;settings=context.network.settings;references=_references(context,group)
    if not references:return None
    reference, = references
    first, last = reference.first, reference.last
    curvatures = (native_bank_curvature(first, reference.contacts[0], -1, reverse=True),
                  native_bank_curvature(last, reference.contacts[1], 1))
    candidate, source, record = build_native_patch(
        context, group, reference, Key('owned_native_cap_chain', 0), settings,
        curvatures, OwnedCapChainCarrier, _contact_phase, POLICY)
    source.layout.owned_native_cap_chain = record
    physical = copy(candidate.physical_layout)
    physical.source_cycle_owned_native_cap_chains = (record,)
    physical.source_cycle_full_width_proof = record.report['complete_source_width']
    return replace(candidate, physical_layout=physical)


def build_native_patch(context, group, reference, owner, settings, curvatures,
                       make_carrier, contact_phase, policy):
    """Fit the two native bank jets and certify the whole curb and paved width."""
    reference.validate(context)
    first, last = reference.first, reference.last
    a, va = first.jet(reference.contacts[0], -.5 * first.width)
    b, vb = last.jet(reference.contacts[1], .5 * last.width)
    if isinstance(reference,OwnedNativeCapChainReference) and first.profile.half:
        # WidthContinuation consumes unit directions. The derivative of a
        # regular smooth parallel curve is a positive multiple of its tangent.
        va=first.frame(reference.contacts[0])[1]
        vb=last.frame(reference.contacts[1])[1]
    roots=(first.point(0.,-.5*first.width),last.point(0.,.5*last.width))
    limit=.25*min(first.width,last.width);attempts=[]
    radii=tuple(CornerRadiusPolicy.contact_candidates(first,last,reference.radius,reference.mode))
    if len(radii)>17:raise OwnedCapChainAuthorityError('Owned cap chain changed the existing finite radius policy')
    for radius in radii:
        try:base=NativeCurvatureCore(WidthContinuation(a,b,-va,vb,radius,refine_curvature=True),curvatures,radius)
        except ContinuationReachError as error:attempts.append((radius,str(error)));continue
        axis=(b-a)/np.linalg.norm(b-a);outward=right(axis);targets=[]
        for point in roots:
            projection=float((point-a)@axis)
            if not 0.<projection<float(np.linalg.norm(b-a)):continue
            low,high=0.,1.
            for _ in range(60):
                middle=.5*(low+high)
                if float((_evaluate(base.controls,middle)-a)@axis)<projection:low=middle
                else:high=middle
            t=.5*(low+high);targets.append((t,float((point-base(t))@outward)))
        for index,multiples in ((None,(0.,1.,2.,4.)),(6,(1.,4.)),(3,(1.,4.))):
            bell=(lambda t:64.*t**3*(1.-t)**3) if index is None else (lambda t:float(ShiftedRootBell.bell(t,index)))
            required=max((0.,*(distance/bell(t) for t,distance in targets)))
            amplitudes=sorted({0. if multiple==0. else required+multiple*settings.curve_tolerance
                for multiple in multiples if multiple==0. or required+multiple*settings.curve_tolerance<=limit})
            for amplitude in amplitudes:
                try:
                    core=bounded_bell(base,outward,amplitude,radius,index=index,maximum_distance=limit)
                    carrier=make_carrier(first,last,reference.contacts,radius,
                        reference.radius,reference.mode,tuple(r for r in radii if r>=radius),reference)
                    curve=DimensionalCompositeCurb(carrier,(_CertifiedContinuation(EnvelopeBoundary(carrier,core)),),settings.curve_tolerance)
                    if owner in {s.region.owner for s in group.footprints}:
                        raise OwnedCapChainAuthorityError('Owned cap chain generated owner namespace collided')
                    source=_patch(owner,carrier,curve,reference.contacts,settings)
                    if contact_phase is None:
                        phase = dict(normal_unchanged_analytical_sampling=True)
                    else:
                        curve, phase = contact_phase(reference, curve, source, settings)
                        source = _patch(owner, carrier, curve, reference.contacts, settings)
                    loops=RegionUnion(max_segments=20000).build(s.region for s in (*group.area_footprints,source))
                    candidate=replace(group,footprints=(*group.footprints,source),loops=loops)
                    inventory={s.region.owner:s for s in candidate.footprints}
                    spans=[rim_label(s,inventory) for loop in loops for s in loop.spans
                        if s.owner==owner and rim_label(s,inventory).kind=='curb' and not loop.hole]
                    exposure=sum(abs(s.last-s.first) for s in spans)
                    if (len(loops)!=1 or loops[0].hole or abs(exposure-1.)>1e-8
                            or any(rim_label(s,inventory).kind=='internal' for loop in loops for s in loop.spans)
                            or CompoundPortPlanner.incomplete_ports(candidate)):
                        raise OwnedCapChainCandidateError('Owned native cap-chain candidate is not continuously exposed with complete mouths')
                    continuous=certify_polynomial_exposure(core.controls,group.area_footprints)
                    width=certify_layout(context,candidate,())
                except (ContinuationReachError,OwnedCapChainCandidateError,NativeC2ExposureError) as error:
                    attempts.append((radius,index,amplitude,type(error).__name__,str(error)));continue
                reference.validate(context)
                report=dict(policy=policy,original_site=group.members[0].text(),
                    original_adjacent_donor_claim=False,non_donor_root_exposure_claim=False,
                    original_sectors=[(e.section_id,e.end,r,m) for e,r,m in reference.sectors],
                    intermediate_native_banks=sorted({w.token for w in reference.intermediate}),
                    original_ledger_contacts=reference.ledger,requested_radius=reference.radius,
                    fitted_radius=radius,curvature_bound=core.curvature_bound,
                    native_endpoint_curvatures=curvatures,base_displacement_bound=base.base_displacement_bound,
                    bell_amplitude=amplitude,bell_index=index,maximum_bell_displacement=limit,
                    full_curve_exposure=exposure,continuous_original_polygon_exposure=continuous,
                    phase=phase,complete_source_width=width,original_operands_retained=True,
                    original_paths_profiles_mouth_planes_unchanged=True,
                    prepared_neighbor_guard='PENDING_COMPLETE_EVENTUAL_OWNER_INVENTORY',attempts=attempts)
                record=PreparedOwnedNativeCapChain(reference,curve,hashlib.sha256(core.controls.tobytes()).hexdigest(),report)
                source.layout.generated_junction_area = record
                source.layout.original_junction = group.members[0]
                return candidate, source, record
    raise OwnedCapChainCandidateError('Authenticated owned native cap chain exhausted its existing radius/bell candidates: '+repr(attempts))


def construction_signature(context,group):
    records=getattr(group.physical_layout,'source_cycle_owned_native_cap_chains',())
    if not records:return None
    result=[]
    from .routes import route_key
    for record in records:
        if type(record) is not PreparedOwnedNativeCapChain:
            raise OwnedCapChainAuthorityError('Owned cap chain needs its exact typed construction record')
        ref=record.reference.validate(context);core=record.curve.pieces[0].curve.continuation
        sources=[s for s in group.area_footprints if getattr(s.layout,'owned_native_cap_chain',None) is record]
        if (len(sources)!=1 or getattr(sources[0].layout,'generated_junction_area',None) is not record
                or type(record.curve) is not ContactPhasedNativeC2
                or type(record.curve.corner) is not OwnedCapChainCarrier
                or record.curve.corner.owned_cap_chain_reference is not ref
                or record.curve.corner.first is not ref.first or record.curve.corner.second is not ref.last
                or record.curve.corner.stations!=ref.contacts):
            raise OwnedCapChainAuthorityError('Owned cap chain lost its real producer marker or derived carrier')
        if hashlib.sha256(core.controls.tobytes()).hexdigest()!=record.control_sha256:
            raise OwnedCapChainAuthorityError('Owned cap chain changed its certified polynomial geometry')
        parameters=record.curve.contact_phase_parameters
        if (parameters[0]!=0. or parameters[-1]!=1. or np.any(np.diff(parameters)<=0.)
                or _chord_bound(core.controls,parameters)>record.curve.contact_phase_error_bound+1e-12
                or record.curve.contact_phase_error_bound>context.network.settings.curve_tolerance):
            raise OwnedCapChainAuthorityError('Owned cap chain changed its certified generated cadence')
        dependencies={Key('settings',0),*group.members};section_ids=set()
        for arm in ref.original_group.physical_layout.approaches:
            dependencies.add(route_key(End(arm.source.id,arm.end)))
            section_ids.update(p.section.id for p in arm.path.pieces)
        for bundle in context.network.bundles:
            if section_ids.intersection(bundle.sections):
                dependencies.add(bundle.key);section_ids.update(bundle.sections)
        dependencies.update(Key(kind,sid) for sid in section_ids for kind in ('section','path','profile'))
        # Runtime identity checks are separate from deterministic cache inputs.
        public_inputs=(fingerprint(context.network.settings),tuple(k.text() for k in ref.original_group.members),
            tuple((a.source.id,a.end,a.width,a.offset,fingerprint(a.profile),
                tuple((p.forward,_path_token(p.path)) for p in a.path.pieces))
                for a in ref.original_group.physical_layout.approaches),
            tuple(ref.original_group.physical_layout.outward_stations),
            tuple((sid,end,fingerprint(mouth)) for (sid,end),mouth in sorted(ref.original_group.physical_layout.mouths.items())))
        result.append((POLICY,tuple((k.text(),context.digest(k)) for k in sorted(dependencies)),
            tuple((e.section_id,e.end,r,m) for e,r,m in ref.sectors),ref.contacts,
            tuple((s.region.owner.text(),hashlib.sha256(s.region.xy.tobytes()).hexdigest(),
                tuple((v.kind,v.index,v.first,v.last) for v in s.segments)) for s in ref.original_group.area_footprints),
            record.control_sha256,tuple(record.curve.contact_phase_parameters),
            record.curve.contact_phase_error_bound,public_inputs))
    return tuple(result)
