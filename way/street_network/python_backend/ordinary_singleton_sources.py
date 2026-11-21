"""A failed ordinary junction can retain its finite original source union.

This is one original physical site and its owned prefixes. It creates neither
a graph cycle nor a consumed road, paired carriageway or new source incidence.
"""
from dataclasses import dataclass, replace
import numpy as np
from .model import Key
from .updates import drain_steps
from .source_strip_containment import SourceStripContainmentError


class OrdinarySourceSingletonUnavailable(ValueError):
    """The source union would change this singleton's original inventory."""


@dataclass(frozen=True)
class OriginalSingletonRootBounds:
    site: object
    layout: object
    failure: SourceStripContainmentError

    def limits(self):
        from .model import End
        return {End(a.source.id,a.end):s-min(.25,.05*a.width)
                for a,s in zip(self.layout.approaches,self.layout.outward_stations)}

    def preserve_mouths(self, approaches, stations):
        """Seed a singleton's search at or beyond each original profile plane."""
        from .intersections import IntersectionProvider
        original = _ends(self.layout)
        result = list(stations)
        for index, (arm, station) in enumerate(zip(approaches, stations)):
            entry = original.get((arm.source.id, arm.end))
            if entry is None or entry[0] is not arm:
                raise OrdinarySourceSingletonUnavailable('Singleton mouth lost its original Approach')
            previous = entry[1]
            if station < previous:
                if not 0. <= previous < arm.length or not IntersectionProvider.mouth_fits(arm, previous):
                    raise OrdinarySourceSingletonUnavailable('Original mouth falls outside its finite approach')
                result[index] = previous
        return result


def _ends(layout):
    if len(layout.approaches)!=len(layout.outward_stations):
        raise OrdinarySourceSingletonUnavailable('Ordinary source union lost a complete original End station')
    result={}
    for a,s in zip(layout.approaches,layout.outward_stations):
        key=a.source.id,a.end
        if key in result:
            raise OrdinarySourceSingletonUnavailable('Ordinary source union repeats an original End')
        result[key]=a,s
    return result


def _same_buffer(first,last):
    first,last=np.asarray(first),np.asarray(last)
    return first.dtype==last.dtype and first.shape==last.shape and first.tobytes()==last.tobytes()


def _same_mouth_representation(a,s,b,t):
    """Keep exact native source cuts despite a reversed Route round trip.

    A transverse section is affine in its original lateral coordinate. Exact
    center/tangent/normal frames establish every point, not a probe subset.
    The two frames and denominator used by Approach.jet establish its complete
    offset Jacobian too; no numerical proximity admits an actual moved mouth.
    """
    if (a.path is not b.path or a.profile is not b.profile or a.source is not b.source
            or a.width!=b.width or a.offset!=b.offset or a.end!=b.end):
        return False
    if s==t:return True  # Keep the original successful path unchanged.
    from .intersections import Approach
    from .routes import Route
    if (type(a) is not Approach or type(b) is not Approach or type(a.path) is not Route
            or not 0.<=s<=a.length or not 0.<=t<=b.length
            or a.path.trims(s)!=b.path.trims(t)):
        return False
    if any(getattr(getattr(arm,name),'__func__',None) is not getattr(Approach,name)
            for arm in (a,b) for name in ('frame','point','jet','physical_gradient')):
        return False
    # These are the exact finite-difference frames used by Approach.jet.
    first_s,last_s=s-1e-4,min(s+1e-4,a.length)
    first_t,last_t=t-1e-4,min(t+1e-4,b.length)
    if not _same_buffer(last_s-first_s,last_t-first_t):return False
    for u,v in ((s,t),(first_s,first_t),(last_s,last_t)):
        if not all(_same_buffer(x,y) for x,y in zip(a.frame(u),b.frame(v))):return False
    q=np.array((-1.,-.5,0.,.5,1.));distance=a.offset+.5*a.width*q
    if not all(_same_buffer(x,y) for x,y in zip(a.jet(s,distance),b.jet(t,distance))):return False
    return _same_buffer(a.physical_gradient(s,q),b.physical_gradient(t,q))


def dimensional_contact_candidate(planner,site,owner,corner,settings):
    """Try the original singleton's certified contact planes before a coarse search.

    A coarse continuation reach can skip an already feasible original plane.
    This one exact candidate has the same radius/curvature/source-width gates;
    failure retains the existing finite search and final mouth admission.
    """
    bounds=getattr(planner,'ordinary_singleton_root_bounds',None)
    if (type(bounds) is not OriginalSingletonRootBounds or bounds.site!=site
            or not isinstance(bounds.failure,SourceStripContainmentError)
            or bounds.layout is not planner.layouts.get(site)):
        return None
    old=_ends(bounds.layout);limits=bounds.limits();stations=[]
    from .model import End
    from .compound_dimensional_corners import prepare
    from .continuations import ContinuationReachError
    for approach in (corner.first,corner.second):
        key=approach.source.id,approach.end
        if key not in old:return None
        original,station=old[key]
        if (approach is not original or approach.profile.half
                or not _same_mouth_representation(original,station,approach,station)):
            return None
        stations.append(limits[End(*key)])
    try:
        record=prepare(owner,corner,settings,tuple(stations),fixed_stations=True)
    except ContinuationReachError:
        return None
    record.report['original_singleton_profile_planes']=[dict(section=approach.source.id,
        end=approach.end,profile_station=old[approach.source.id,approach.end][1],
        contact_station=contact) for approach,contact in zip((corner.first,corner.second),stations)]
    return record


def _dimensional_lead_reserve(context,group,a,s,b,t):
    """Admit one existing normal reserve beyond a proved wider-root contact.

    A root-anchored width transition can reach the old profile plane itself.
    Its normal source union still needs an exposed complete profile beyond
    that curb. This is real, recorded movement, never a station alias or a
    permission to advance a root-support donor.
    """
    from .intersections import Approach,IntersectionProvider
    from .compound_dimensional_corners import DimensionalCorner
    if (type(a) is not Approach or type(b) is not Approach
            or a.path is not b.path or a.profile is not b.profile or a.source is not b.source
            or a.end!=b.end or a.width!=b.width or a.offset!=b.offset
            or not 0.<=s<t<a.length or not IntersectionProvider.mouth_fits(b,t)):
        return None
    old_cuts,new_cuts=a.path.trims(s),b.path.trims(t)
    if set(old_cuts)!=set(new_cuts) or any(
            (new_cuts[key]<value if key[1]==0 else new_cuts[key]>value)
            for key,value in old_cuts.items()):
        return None
    for corner in group.physical_layout.corners:
        prepared=getattr(corner,'dimensional_corner',None)
        if (not isinstance(prepared,DimensionalCorner)
                or not any(source is prepared.source for source in group.area_footprints)
                or not any(owned is corner for owned in prepared.source.layout.corners)
                or prepared.report.get('policy')!='FINITE_ORIGINAL_WIDTH_CONTINUATION'
                or prepared.stations!=corner.stations):
            continue
        arms=corner.first,corner.second
        wide=int(arms[1].width>arms[0].width)
        if arms[wide].width<=arms[1-wide].width:
            continue
        contact=prepared.stations[wide]
        if not _same_mouth_representation(a,s,arms[wide],contact):
            continue
        reserve=min(.25,.05*a.width)
        expected=IntersectionProvider.mouth_station(a,contact+reserve,context.network.settings)
        if not _same_mouth_representation(a,expected,b,t):
            continue
        return dict(policy='EXISTING_FINITE_WIDTH_CONTINUATION_LEAD_RESERVE',
            section=a.source.id,end=a.end,original_route_station=s,
            prepared_wide_contact=contact,candidate_route_station=t,
            advance=t-s,normal_reserve=reserve,normal_mouth_station=expected,
            original_source_cuts=[dict(section=sid,end=end,old=value,new=new_cuts[sid,end])
                for (sid,end),value in old_cuts.items()],
            minimum_radius=corner.radius,curvature_bound=prepared.report['curvature_bound'],
            complete_original_profile_preserved=True,source_identity_preserved=True)
    return None


def _dimensional_recipient_lead(context,group,a,s,b,t):
    """Admit the existing narrow lead only after its old profile is clipped.

    This records an actual finite mouth advance. It neither treats it as a
    station alias nor changes the already certified continuation geometry.
    """
    from .intersections import Approach,IntersectionProvider
    from .routes import Route
    from .compound_dimensional_corners import DimensionalCorner,_CertifiedContinuation
    if (a is not b or type(a) is not Approach or type(a.path) is not Route
            or len(group.members)!=1 or group.claims
            or not 0.<=s<t<a.length or not IntersectionProvider.mouth_fits(a,t)
            or any(getattr(group.physical_layout,name,()) for name in
                   ('source_cycle_root_envelopes','source_cycle_root_support_envelopes'))):
        return None
    old_cuts,new_cuts=a.path.trims(s),a.path.trims(t)
    if set(old_cuts)!=set(new_cuts) or any(new_cuts[key]<value if key[1]==0 else new_cuts[key]>value
                                         for key,value in old_cuts.items()):
        return None
    for corner in group.physical_layout.corners:
        prepared=getattr(corner,'dimensional_corner',None)
        if (type(prepared) is not DimensionalCorner
                or not any(source is prepared.source for source in group.area_footprints)
                or not any(value is corner for value in prepared.source.layout.corners)
                or prepared.report.get('policy')!='FINITE_ORIGINAL_WIDTH_CONTINUATION'
                or prepared.stations!=corner.stations):
            continue
        arms=corner.first,corner.second
        narrow=int(arms[1].width<arms[0].width);wide=1-narrow
        if arms[narrow] is not a or arms[wide].width<=a.width:continue
        cores=[piece.curve for piece in prepared.curve.pieces if type(piece) is _CertifiedContinuation]
        if (len(cores)!=1 or not any(curve is prepared.curve for curve in prepared.source.layout.boundaries)
                or cores[0].corner.first is not arms[0] or cores[0].corner.second is not arms[1]
                or cores[0].corner.radius!=corner.radius):
            continue
        curvature=float(cores[0].continuation.curvature_bound)
        reserve=min(.25,.05*a.width);contact=prepared.stations[narrow]
        if (not s<contact<t
                or t!=IntersectionProvider.mouth_station(a,contact+reserve,context.network.settings)
                or prepared.report.get('minimum_radius')!=corner.radius
                or prepared.report.get('curvature_bound')!=curvature
                or not np.isfinite(curvature) or not 0.<=curvature<=1./corner.radius):
            continue
        if (a.profile!=context[Key('profile',a.source.id)].oriented(a.end)
                or any(piece.path is not context[Key('path',piece.section.id)] for piece in a.path.pieces)):
            continue
        ports=[source for source in group.area_footprints if len(getattr(source.layout,'approaches',()))==1
               and source.layout.approaches[0] is a and source.layout.outward_stations==(t,)]
        if len(ports)!=1:continue
        port=ports[0]
        from .source_curb_coordinates import CurbCoordinates
        grids={piece.section.id:[0.,piece.path.length] for piece in a.path.pieces}
        for curve in port.layout.boundaries:
            coordinates=getattr(curve,'source_coordinates',None)
            if type(coordinates) is not CurbCoordinates:return None
            coordinates.validate(curve.points,context)
            for edge in range(len(coordinates.stations)-1):
                for piece in coordinates.edge_sources(edge):
                    if piece.path.section.id not in grids:return None
                    grids[piece.path.section.id].extend((piece.first_station,piece.last_station))
        grids={sid:np.unique(values) for sid,values in grids.items()}
        from .compound_cycles import _port_footprint
        settings=context.network.settings
        previous=_port_footprint(port.region.owner,a,s,port.layout.radial,settings.spacing,
            min(.001,.05*settings.curve_tolerance),grids)
        sources=tuple(previous if source is port else source for source in group.area_footprints)
        if sum(len(source.region.xy) for source in sources)>20000:
            raise OrdinarySourceSingletonUnavailable('Dimensional recipient evidence exceeds its original source-sample budget')
        from .regions import RegionUnion
        from .incomplete_profile_leads import _exposed
        loops=RegionUnion(max_segments=20000).build(source.region for source in sources)
        old_measure=_exposed(previous,sources,loops)
        new_measure=_exposed(port,group.footprints,group.loops)
        if not 0.<=old_measure<1.-1e-8 or new_measure<1.-1e-8:return None
        return dict(policy='EXISTING_FINITE_WIDTH_CONTINUATION_RECIPIENT_LEAD',
            section=a.source.id,end=a.end,original_route_station=s,candidate_route_station=t,
            advance=t-s,prepared_narrow_contact=contact,normal_reserve=reserve,
            minimum_radius=corner.radius,curvature_bound=prepared.report['curvature_bound'],
            actual_old_exposed_profile_measure=old_measure,actual_new_exposed_profile_measure=new_measure,
            actual_old_profile_incomplete=True,complete_original_profile_preserved=True,
            original_source_cuts=[dict(section=sid,end=end,old=value,new=new_cuts[sid,end])
                for (sid,end),value in old_cuts.items()],
            original_path_profile_and_incidence_unchanged=True,no_root_donor_override=True)
    return None


def build(context,site,layout,layouts,junctions,failure):
    return drain_steps(build_steps(context, site, layout, layouts, junctions, failure))


def build_steps(context,site,layout,layouts,junctions,failure):
    """Use existing source patches after every ordinary width family failed."""
    from .ordinary_mouth_clearance import _eligible,ORDINARY_MOUTH_POLICY
    if not isinstance(failure,SourceStripContainmentError) or not _eligible(context,layout):
        raise OrdinarySourceSingletonUnavailable('Ordinary source union needs a real complete-prefix failure')
    if any(abs(a.offset)>1e-8 for a in layout.approaches):
        raise OrdinarySourceSingletonUnavailable('Ordinary source union needs centred original carriageways')
    old=_ends(layout)
    paths={piece.section.id:piece.path for a,s in old.values() for piece in a.path.pieces}
    if failure.section_id not in paths:
        raise OrdinarySourceSingletonUnavailable('Ordinary source union lost its failing original Section')
    from .compound_preparation import DeferredCompoundPlanner
    from .compound_cycles import ExposedPreparationSeam, source_cycle_group_steps
    planner=DeferredCompoundPlanner(context,junctions,layouts)
    planner.ordinary_singleton_root_bounds=OriginalSingletonRootBounds(site,layout,failure)
    from .local_width_continuation import fit_steps, wider_approach_steps
    from .ordinary_mouth_clearance import _canonical_candidate, _group, record_width_transition
    fitted = yield from fit_steps(context, site, layout, failure.section_id, _canonical_candidate,
                                  extrapolated=True)
    if fitted is None:
        try:
            group=yield from source_cycle_group_steps(planner,(site,),set(),set(),{},set())
        except ExposedPreparationSeam:
            yield
            fitted = yield from fit_steps(context, site, layout, failure.section_id, _canonical_candidate)
            if fitted is None:
                fitted = yield from wider_approach_steps(context, site, layout, failure.section_id, _canonical_candidate)
            if fitted is None:
                raise
    if fitted is not None:
        candidate, proof, record = fitted
        candidate = record_width_transition(candidate, site, proof, [], 0, record)
        return _group(site, candidate)
    physical=group.physical_layout
    if group.members!=(site,) or group.claims or physical is None:
        raise OrdinarySourceSingletonUnavailable('Ordinary source union changed its original single owner')
    new=_ends(physical)
    if set(new)!=set(old):
        raise OrdinarySourceSingletonUnavailable('Ordinary source union changed its original complete mouths')
    aliases=[];advances=[];expected_mouths={}
    for key,(a,s) in old.items():
        b,t=new[key]
        if not _same_mouth_representation(a,s,b,t):
            record=_dimensional_lead_reserve(context,group,a,s,b,t)
            if record is None:record=_dimensional_recipient_lead(context,group,a,s,b,t)
            if record is None:
                from .source_root_witnesses import admission as staged_recipient_admission
                record=staged_recipient_admission(context,group,a,s,b,t)
            if record is None:
                from .incomplete_profile_leads import admission
                record=admission(context,group,a,s,b,t)
            if record is None:
                raise OrdinarySourceSingletonUnavailable('Ordinary source union changed its original Path, Profile or station')
            advances.append(record)
        elif s!=t:
            aliases.append(dict(section=key[0],end=key[1],original_route_station=s,candidate_route_station=t,
                exact_original_source_stations=[dict(section=sid,end=end,station=value)
                    for (sid,end),value in a.path.trims(s).items()],
                complete_physical_frame_and_offset_jacobian_exact=True,physical_gradient_exact=True))
        for token,value in b.path.trims(t).items():
            if token not in layout.mouths:
                raise OrdinarySourceSingletonUnavailable('Ordinary source union added an original mouth incidence')
            mouth=replace(layout.mouths[token],station=float(value))
            if token in expected_mouths and expected_mouths[token]!=mouth:
                raise OrdinarySourceSingletonUnavailable('Ordinary source union has ambiguous original mouth cuts')
            expected_mouths[token]=mouth
    if physical.mouths!=expected_mouths:
        raise OrdinarySourceSingletonUnavailable('Ordinary source union changed its original complete mouths')
    for sid,path in paths.items():
        if path is not context[Key('path',sid)]:
            raise OrdinarySourceSingletonUnavailable('Ordinary source union changed its original Path ownership')
    for corner in layout.corners:
        if corner.radius_mode!='STRICT':continue
        matches=[c for c in physical.corners if (c.first.source.id,c.first.end,c.second.source.id,c.second.end)
            ==(corner.first.source.id,corner.first.end,corner.second.source.id,corner.second.end)]
        if (len(matches)!=1 or matches[0].radius!=corner.radius
                or matches[0].requested_radius!=corner.requested_radius or matches[0].radius_mode!='STRICT'):
            raise OrdinarySourceSingletonUnavailable('Ordinary source union changed an explicit STRICT radius')
    from .cycle_dimensional_retry import certify_layout
    from .sampling import AnalyticalBoundarySampler
    before=certify_layout(context,group,())
    settings=context.network.settings
    group=AnalyticalBoundarySampler(tolerance=settings.curve_tolerance,
        spacing=settings.spacing,max_samples=20000).build(group)
    if getattr(group.physical_layout, 'shallow_bank_fillet', None):
        # Keep the full native arc and original operands for profiles and
        # height evaluation. Its sub-millimetre samples need not be quad rows.
        from .source_rim_normalization import normalize_exterior
        group = normalize_exterior(group)
    after=certify_layout(context,group,())
    group.physical_layout.source_cycle_full_width_proof=after
    group.physical_layout.ordinary_complete_mouth_clearance=dict(policy=ORDINARY_MOUTH_POLICY,
        original_site=site.text(),correction='FINITE_ORIGINAL_SINGLE_JUNCTION_SOURCE_UNION',
        original_members=[site.text()],consumed_sections=[],
        original_corners_contacts_radii_preserved=False,strict_authored_radii_preserved=True,
        original_paths_profiles_incidence_preserved=True,original_mouth_planes_preserved=not advances,
        complete_source_grid=True,complete_source_width=after,
        route_station_representation_aliases=aliases,
        dimensional_lead_reserve_advances=[r for r in advances
            if r['policy']=='EXISTING_FINITE_WIDTH_CONTINUATION_LEAD_RESERVE'],
        dimensional_recipient_lead_advances=[r for r in advances
            if r['policy']=='EXISTING_FINITE_WIDTH_CONTINUATION_RECIPIENT_LEAD' and not r.get('admission_stage')],
        incomplete_profile_lead_advances=[r for r in advances
            if r['policy']=='TYPED_EXISTING_INCOMPLETE_PROFILE_LEAD_SEARCH'],
        before_canonical_complete_source_width=before,
        source_semantics='one_original_junction_and_its_owned_prefixes')
    staged=[r for r in advances if r.get('admission_stage')=='INDEPENDENT_PRE_ROOT_WIDTH_TRANSITION']
    if staged:group.physical_layout.ordinary_complete_mouth_clearance['staged_dimensional_recipient_lead_advances']=staged
    return group
