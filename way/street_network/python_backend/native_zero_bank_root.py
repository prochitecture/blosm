"""Source geometry for an original full zero-bank endpoint and native half bank.

The existing full/interior root providers are unchanged.  NativeRolePort
authenticates actual ProfileProvider roles; the endpoint provider below is
explicitly one-sided at a genuine original Path endpoint.
"""
from dataclasses import dataclass,replace
import hashlib,math
import numpy as np
from .model import End,Key,fingerprint
from .generated_junction_curbs import GeneratedJunctionArea
from .compound_root_envelopes import RootEnvelope
from .source_root_ports import native_bank_curvature
from .corner_radii import CornerRadiusPolicy,CornerContactError
from .corner_source_support import complete_source_roots
from .compounds import CompoundGroup,rim_label
from .compound_cycles import _canonical_samples,_port_footprint
from .corner_patch import corner_patch
from .regions import RegionUnion
from .cycle_dimensional_retry import certify_layout
from .intersections import IntersectionProvider
from .source_root_ports import NativeRolePort
from . import native_owned_cap_chains as common

POLICY='1-authenticated-original-full-endpoint-to-native-half-zero-bank-C2'


class NativeZeroBankAuthorityError(RuntimeError):pass


@dataclass(frozen=True)
class OriginalEndpointBankJet:
    port: object
    side: int

    def curvature(self,context,*,reverse=False):
        from .mesh import cross
        port=self.port;a=port.approach
        NativeRolePort.build(context,port.footprint,a)
        if (a.profile.half or self.side not in (-1,1) or a.profile.value(self.side)!=0.
                or len(a.path.pieces)!=1 or a.path.stations[0]!=0.):
            raise NativeZeroBankAuthorityError('Endpoint jet needs its actual complete zero-height original bank')
        piece=a.path.pieces[0];path=piece.path
        native=0. if piece.forward else path.length;parameter=path.parameter(native)
        if parameter not in (path.knots[0],path.knots[-1]):
            raise NativeZeroBankAuthorityError('Endpoint jet cannot borrow an interior authored knot or Route seam')
        if path.exact:return 0.
        _,velocity,acceleration=path.jet(parameter)
        k=float(cross(velocity,acceleration)/np.linalg.norm(velocity)**3)*(1. if piece.forward else -1.)
        denominator=1.+.5*self.side*a.width*k
        if not math.isfinite(denominator) or denominator<=0.:
            raise NativeZeroBankAuthorityError('Endpoint bank has no regular actual offset curvature')
        return k/denominator*(-1. if reverse else 1.)


def _base(context,site,raw,settings):
    paths={p.section.id:p.path for a in raw.approaches for p in a.path.pieces};cuts={sid:set() for sid in paths}
    for c in raw.corners:
        for a,s in zip((c.first,c.second),c.stations):
            for (sid,end),at in a.path.trims(max(0.,s)).items():cuts[sid].add(at)
    for (sid,end),mouth in raw.mouths.items():cuts[sid].add(mouth.station)
    grids=_canonical_samples(paths,cuts,settings.spacing,settings.curve_tolerance);areas=[]
    for i,c in enumerate(raw.corners):
        area=corner_patch(Key('zero_bank_original_corner',i),c,settings.spacing,settings.curve_tolerance)
        area.layout.original_junction=site;areas.append(area)
    for a,s in zip(raw.approaches,raw.outward_stations):
        areas.append(_port_footprint(Key('cycle_port',2*a.source.id+a.end),a,s,raw.radial,
            settings.spacing,settings.curve_tolerance,grids))
    loops=RegionUnion(max_segments=20000).build(s.region for s in areas)
    return CompoundGroup((site,),(),tuple(areas),loops,raw),grids


@dataclass(frozen=True)
class NativeZeroBankReference:
    original_group: object
    first: object
    last: object
    ports: tuple
    sectors: tuple
    radius: float
    mode: str
    contacts: tuple
    ledger: tuple
    bounds: object
    input_token: tuple
    geometry_token: tuple
    donor: object
    grids: object
    intermediate: tuple=()

    @classmethod
    def build(cls,planner,site,first,last,radius,mode,settings):
        from .ordinary_singleton_sources import OriginalSingletonRootBounds
        from .source_strip_containment import SourceStripContainmentError
        bounds=getattr(planner,'ordinary_singleton_root_bounds',None)
        if (type(bounds) is not OriginalSingletonRootBounds or bounds.site!=site or mode!='STYLE'
                or not isinstance(bounds.failure,SourceStripContainmentError)
                or bounds.layout is not planner.layouts.get(site)):return None
        raw=bounds.layout;context=planner.context
        if (first.profile.half or last.profile.half!=-1 or first.width!=last.width
                or first.profile.value(-1.)!=0. or last.profile.value(1.)!=0.
                or first.profile.datum!=last.profile.datum):return None
        if (len(first.path.pieces)!=1 or len(last.path.pieces)!=1
                or any(a.offset!=0. for a in raw.approaches)):return None
        group,grids=_base(context,site,raw,settings);token=common._input_token(context,group)
        actual_site=context.network.records()[site]
        ordered=sorted(raw.approaches,key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
        pairs=[(a,ordered[(i+1)%len(ordered)]) for i,a in enumerate(ordered)]
        if not any(a is first and b is last for a,b in pairs):return None
        if actual_site.corner_radius(End(first.source.id,first.end))!=(radius,mode):
            raise NativeZeroBankAuthorityError('Native zero bank changed its authored sector policy')
        donors=[c for c in raw.corners if c.first is first and c.second is last]
        if len(donors)!=1 or donors[0].center is None:return None
        donor,=donors
        if abs(donor.gap-math.pi)>math.radians(15.):return None
        refused=[]
        def accept(c):return complete_source_roots(c,settings.spacing,settings.curve_tolerance,
            support_approaches=ordered,audit_all_contacts=True,audit_facing_halves=True,on_refusal=refused.append)
        try:CornerRadiusPolicy.solve(first,last,radius,mode,math.radians(settings.continuation_angle_degrees),accept=accept)
        except CornerContactError:pass
        else:return None
        if not any(r.get('section')==first.source.id and r.get('end')==first.end for r in refused):return None
        ports=tuple(NativeRolePort.build(context,common._port(group,a),a) for a in raw.approaches)
        by_arm={id(p.approach):p for p in ports};half_port=by_arm[id(last)]
        # This is an EXISTING original corner's real native contact, including
        # its exact native parameter. No approximate equality of XY is used.
        choices=[(float(s),c) for c in raw.corners for a,s in zip((c.first,c.second),c.stations) if a is last and s>=0.]
        contact,producer=max(choices,key=lambda row:row[0])
        if IntersectionProvider.mouth_station(last,contact+min(.25,.05*last.width),context.network.settings)!=half_port.station:
            raise NativeZeroBankAuthorityError('Native zero bank lost its actual original complete profile lead')
        contacts=(0.,contact);ledger=[]
        for arm,side,station in ((first,-1,0.),(last,1,contact)):
            bank=by_arm[id(arm)].footprint.layout.boundaries[0 if side==1 else 1]
            matches=[i for i,s in enumerate(bank.source_coordinates.stations)
                if arm.path.trims(float(s))==arm.path.trims(station)]
            if len(matches)!=1:raise NativeZeroBankAuthorityError('Zero-bank contact lost its unique exact native source parameter')
            i,=matches
            ledger.append((arm.source.id,arm.end,side,station,float(bank.source_coordinates.stations[i]),
                tuple(sorted(arm.path.trims(station).items())),tuple(bank.points[i])))
        geometry=tuple((s.region.owner,s.region.xy.tobytes(),s.segments) for s in group.area_footprints)
        return cls(group,first,last,ports,((End(first.source.id,first.end),radius,mode),),radius,mode,
            contacts,tuple(ledger),bounds,token,geometry,donor,grids)

    def validate(self,context):
        from .ordinary_singleton_sources import OriginalSingletonRootBounds
        from .ownership import footprint_key
        if (type(self.bounds) is not OriginalSingletonRootBounds or self.bounds.layout is not self.original_group.physical_layout
                or self.bounds.site!=self.original_group.members[0]
                or context.products.get(footprint_key(self.bounds.site)) is not self.bounds.layout
                or common._input_token(context,self.original_group)!=self.input_token):
            raise NativeZeroBankAuthorityError('Zero-bank authority changed its original source/Profile/mouth bounds')
        if tuple((s.region.owner,s.region.xy.tobytes(),s.segments) for s in self.original_group.area_footprints)!=self.geometry_token:
            raise NativeZeroBankAuthorityError('Zero-bank authority changed its original area operands')
        for port in self.ports:NativeRolePort.build(context,port.footprint,port.approach)
        site=context.network.records()[self.bounds.site]
        if (site.corner_radius(End(self.first.source.id,self.first.end))!=(self.radius,self.mode)
                or self.sectors!=((End(self.first.source.id,self.first.end),self.radius,self.mode),)
                or self.mode!='STYLE' or not any(c is self.donor for c in self.bounds.layout.corners)):
            raise NativeZeroBankAuthorityError('Zero-bank authority changed its actual authored sector or donor')
        settings=context.network.settings
        rebuilt,_=_base(context,self.bounds.site,self.bounds.layout,
            replace(settings,curve_tolerance=min(.001,.05*settings.curve_tolerance)))
        rebuilt_token=tuple((s.region.owner,s.region.xy.tobytes(),s.segments) for s in rebuilt.area_footprints)
        if rebuilt_token!=self.geometry_token:
            raise NativeZeroBankAuthorityError('Zero-bank authority differs from rebuilt original source geometry')
        if self.contacts[0]!=0. or self.first.profile.half or self.last.profile.half!=-1:
            raise NativeZeroBankAuthorityError('Zero-bank authority lost its actual original full endpoint/native half role')
        expected=max(float(s) for c in self.bounds.layout.corners for a,s in zip((c.first,c.second),c.stations) if a is self.last)
        if self.contacts[1]!=expected or self.first.profile.value(-1.)!=0. or self.last.profile.value(1.)!=0.:
            raise NativeZeroBankAuthorityError('Zero-bank authority changed its actual native zero contact')
        for sid,end,side,station,ledger_station,native,point in self.ledger:
            arm=self.first if (sid,end)==(self.first.source.id,self.first.end) else self.last
            if tuple(sorted(arm.path.trims(station).items()))!=native or arm.path.trims(station)!=arm.path.trims(ledger_station):
                raise NativeZeroBankAuthorityError('Zero-bank authority forged an original native parameter alias')
            port=next(p for p in self.ports if p.approach is arm)
            bank=port.footprint.layout.boundaries[0 if side==1 else 1]
            indices=[i for i,s in enumerate(bank.source_coordinates.stations) if float(s)==ledger_station]
            if len(indices)!=1 or tuple(bank.points[indices[0]])!=point:
                raise NativeZeroBankAuthorityError('Zero-bank authority changed its exact native contact coordinate')
        return self

    def curvatures(self,context):
        self.validate(context)
        port=next(p for p in self.ports if p.approach is self.first)
        return OriginalEndpointBankJet(port,-1).curvature(context,reverse=True),native_bank_curvature(self.last,self.contacts[1],1)


def _carrier(first,last,stations,radius,requested,mode,attempts,reference):
    corner=replace(reference.donor,stations=stations,radius=radius,requested_radius=requested,
        radius_mode=mode,radius_attempts=attempts,center=None,continuation=True,
        fit_policy='ORIGINAL_ENDPOINT_NATIVE_ZERO_BANK_C2')
    corner.native_zero_bank_reference=reference
    return corner


@dataclass(frozen=True)
class NativeZeroBankEnvelope(RootEnvelope,GeneratedJunctionArea):
    reference: object
    control_sha256: str
    retained_sources: tuple


def prepare(planner,site,owner,first,last,radius,mode,settings):
    reference=NativeZeroBankReference.build(planner,site,first,last,radius,mode,settings)
    if reference is None:return None
    _, source, temporary = common.build_native_patch(
        planner.context, reference.original_group, reference, owner, settings,
        reference.curvatures(planner.context), _carrier, None, POLICY)
    old=next(s for s in reference.original_group.area_footprints if any(c is reference.donor for c in s.layout.corners))
    retained=replace(old,region=replace(old.region,owner=Key('retained_zero_bank_original_circle',2*first.source.id+first.end)))
    report=dict(temporary.report,policy=POLICY,original_endpoint_jet=True,
        original_circle_retained=True,all_native_profiles_and_original_mouths_preserved=True,
        original_native_contact_alias=reference.ledger,normal_analytical_sampling=True)
    envelope=NativeZeroBankEnvelope(source,temporary.curve,reference.contacts,report,temporary.curve.corner,
        reference,temporary.control_sha256,(retained,))
    source.layout.generated_junction_area=envelope;source.layout.native_zero_bank_root=envelope
    envelope.corner.root_envelope=envelope
    return envelope


def validate_group(context,group):
    records=tuple(getattr(c,'root_envelope',None) for c in group.physical_layout.corners
        if type(getattr(c,'root_envelope',None)) is NativeZeroBankEnvelope)
    if not records:return group
    inventory={s.region.owner:s for s in group.footprints}
    for record in records:
        record.reference.validate(context)
        core=record.curve.pieces[0].curve.continuation
        if hashlib.sha256(core.controls.tobytes()).hexdigest()!=record.control_sha256:
            raise NativeZeroBankAuthorityError('Native zero bank changed its certified polynomial controls')
        exposed=sum(abs(rim_label(s,inventory).last-rim_label(s,inventory).first)
            for l in group.loops for s in l.spans if not l.hole and s.owner==record.source.region.owner
            and rim_label(s,inventory).kind=='curb')
        if abs(exposed-1.)>1e-8:raise NativeZeroBankAuthorityError('Final native zero bank lost its complete generated curb')
    group.physical_layout.source_cycle_native_zero_bank_roots=records
    group.physical_layout.source_cycle_full_width_proof=certify_layout(context,group,())
    return group


def signature(context,group):
    records=getattr(group.physical_layout,'source_cycle_native_zero_bank_roots',())
    if not records:return None
    from .routes import route_key
    values=[]
    for record in records:
        if type(record) is not NativeZeroBankEnvelope:raise NativeZeroBankAuthorityError('Native zero bank lost its typed producer')
        ref=record.reference.validate(context)
        areas=[s for s in group.area_footprints if getattr(s.layout,'native_zero_bank_root',None) is record]
        if (len(areas)!=1 or getattr(areas[0].layout,'generated_junction_area',None) is not record
                or record.curve.corner is not record.corner or record.corner.root_envelope is not record
                or record.corner.first is not ref.first or record.corner.second is not ref.last
                or record.corner.stations!=ref.contacts):
            raise NativeZeroBankAuthorityError('Native zero bank changed its actual producer/marker/contact binding')
        if hashlib.sha256(record.curve.pieces[0].curve.continuation.controls.tobytes()).hexdigest()!=record.control_sha256:
            raise NativeZeroBankAuthorityError('Native zero bank changed its certified C2 geometry token')
        settings=context.network.settings
        rebuilt=common._patch(areas[0].region.owner,record.corner,record.curve,record.stations,
            replace(settings,curve_tolerance=min(.001,.05*settings.curve_tolerance)))
        if rebuilt.region.xy.tobytes()!=areas[0].region.xy.tobytes() or rebuilt.segments!=areas[0].segments:
            raise NativeZeroBankAuthorityError('Native zero bank source polygon differs from its exact declared curve sampling')
        inventory={s.region.owner:s for s in group.footprints}
        exposed=sum(abs(rim_label(s,inventory).last-rim_label(s,inventory).first)
            for loop in group.loops for s in loop.spans if not loop.hole and s.owner==areas[0].region.owner
            and rim_label(s,inventory).kind=='curb')
        if abs(exposed-1.)>1e-8 or any(s.region.owner not in inventory for s in record.retained_sources):
            raise NativeZeroBankAuthorityError('Native zero bank lost complete final exposure or a retained original circle')
        ids={p.section.id for a in ref.original_group.physical_layout.approaches for p in a.path.pieces}
        dependencies={*group.members,Key('settings',0)}
        for a in ref.original_group.physical_layout.approaches:dependencies.add(route_key(End(a.source.id,a.end)))
        for b in context.network.bundles:
            if ids.intersection(b.sections):dependencies.add(b.key);ids.update(b.sections)
        dependencies.update(Key(kind,sid) for sid in ids for kind in ('section','path','profile'))
        values.append((POLICY,tuple((k.text(),context.digest(k)) for k in sorted(dependencies)),ref.contacts,
            record.control_sha256,tuple((a.source.id,a.end,a.width,a.offset,fingerprint(a.profile),
                tuple(common._path_token(p.path) for p in a.path.pieces)) for a in ref.original_group.physical_layout.approaches),
            tuple(ref.original_group.physical_layout.outward_stations)))
    return tuple(values)
