"""Physical holes of a compound's closed original source-incidence cycles.

An internal graph cycle does not itself prove a grass island: finite source
carriageways and their original local fillets can overlap and close it. Compute
the remaining holes from those physical regions, keeping their bounded curb
samples and provenance. This never connects Streets by their coordinates.
"""
from types import SimpleNamespace
import math
import numpy as np

from .compounds import RimSegment, SampledFootprint
from .model import End, Key
from .updates import drain_steps
from .regions import BoundaryIncidenceError, Region, RegionUnion, signed_area


class ExposedPreparationSeam(ValueError):
    """An internal preparation edge survived on the physical pavement rim."""


class PolylineCurb:
    """A retained bounded-error physical polyline, including every corner."""
    continuation = None

    def __init__(self, points):
        self.points = np.asarray(points)
        self.contact_reference_points = self.points
        self.stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(self.points, axis=0), axis=1))]
        self.length = float(self.stations[-1])
        self.corner = SimpleNamespace(center=None)
        if self.length <= 1e-8:
            raise ValueError('Source-cycle curb has no positive physical length')

    def __call__(self, t):
        station = np.asarray(t)*self.length
        return np.stack([np.interp(station, self.stations, self.points[:, k]) for k in (0, 1)], axis=-1)

    def certified_parameters(self, sampler, first, last):
        knots = sorted({first, last, *(float(s/self.length) for s in self.stations
                                     if first*self.length < s < last*self.length)})
        values = []
        for a, b in zip(knots, knots[1:]):
            values.extend(sampler.regular(a, b, sampler.spacing/self.length)[:-1])
        return np.r_[values, last]


def _footprint(owner, points, *, internal_edges=(),source_coordinates=None):
    points = np.asarray(points, dtype=float)
    internal={tuple(sorted((tuple(points[i]),tuple(points[(i+1)%len(points)])))) for i in internal_edges}
    keep = np.r_[True, np.linalg.norm(np.diff(points, axis=0), axis=1)>1e-8]
    indices=np.flatnonzero(keep)
    points = points[keep]
    if len(points)>2 and np.linalg.norm(points[0]-points[-1])<=1e-8:
        points=points[:-1]
        indices=indices[:-1]
    if len(points)<3 or abs(signed_area(points))<=1e-10:
        return None
    if signed_area(points)<0.:
        points=points[::-1]
        indices=indices[::-1]
    # Two finite analytical intervals avoid a wrapping parameter on a closed
    # curve. Their exact chord joins survive later contact redistribution.
    split = len(points)//2
    curves = (PolylineCurb(points[:split+1]), PolylineCurb(np.vstack((points[split:],points[0]))))
    if source_coordinates is not None:
        retained=source_coordinates.subset(indices)
        curves[0].source_coordinates=retained.subset(np.arange(split+1))
        curves[1].source_coordinates=retained.subset(np.r_[np.arange(split,len(points)),0])
    labels=[]
    for i, curve in enumerate(curves):
        labels.extend(RimSegment('curb',i,float(a/curve.length),float(b/curve.length))
                      for a,b in zip(curve.stations,curve.stations[1:]))
    if internal:
        from dataclasses import replace
        labels=[replace(label,kind='internal') if tuple(sorted((tuple(a),tuple(b)))) in internal else label
                for a,b,label in zip(points,np.roll(points,-1,axis=0),labels)]
    layout=SimpleNamespace(approaches=(),boundaries=curves,rows=len(points),mouths={},corners=())
    return SampledFootprint(Region(owner,points),tuple(labels),layout)


def shared_bundle_sites(network, bundle):
    """A complete two-sided source pair must share both original terminal sites."""
    if len(bundle.members)!=2 or bundle.interior_section_ids or not all(m.complete for m in bundle.members):
        return ()
    by_end={end:site for site in network.intersections for end in site.ends}
    terminals=[]; sites=set()
    for member in bundle.members:
        first,last=member.path[0],member.path[-1]
        ends=(End(first.section_id,0 if first.forward else 1),End(last.section_id,1 if last.forward else 0))
        owners=tuple(by_end.get(end) for end in ends)
        if any(owner is None for owner in owners):return ()
        terminals.append(tuple(owner.key for owner in owners))
        for visit in member.path:
            for end in (0,1):
                site=by_end.get(End(visit.section_id,end))
                if site is not None:sites.add(site.key)
    if terminals[0]!=terminals[1] or terminals[0][0]==terminals[0][1]:return ()
    return tuple(sorted(sites))


def _sampling_intervals(paths,consumed,approaches,stations):
    """Rendering-only owned intervals, unsuitable for future-source proofs."""
    result={sid:[(0.,paths[sid].length)] for sid in consumed}
    for approach,station in zip(approaches,stations):
        if hasattr(approach,'original_approaches'):
            # These finite original ranges also support the terminal's own
            # prepared footprint. They conservatively include its mouth.
            for path,interval in zip(approach.paths,approach.source_ranges):
                result.setdefault(path.section.id,[]).append(tuple(interval))
            continue
        for (sid,end),cut in approach.path.trims(station).items():
            result.setdefault(sid,[]).append((0.,cut) if end==0 else (cut,paths[sid].length))
    return result


class _CanonicalSamples(dict):
    """Distinguish complete source grids from rendering-only interval grids."""
    def __init__(self,*,complete_source_ranges):
        super().__init__()
        self.complete_source_ranges=complete_source_ranges


def _canonical_samples(paths,cuts,spacing,tolerance,*,intervals=None,max_samples=20000):
    """Keep the original compulsory-cut phase within a finite allocation.

    The default samples every original finite Path, including future unowned
    tails used by mouth search and continuous exclusion proofs. Explicit
    intervals are only for rebuilding current physical source operands.
    """
    from .paths import PathSampleBudgetExceeded
    try:
        return _sample_path_ranges(paths, cuts, spacing, tolerance, intervals, max_samples, False)
    except PathSampleBudgetExceeded:
        if all(path.exact for path in paths.values()):
            raise
        # An isolated tight bend can over-refine a whole long carriageway.
        # Repeat the complete ledger using certified per-piece curvature.
        return _sample_path_ranges(paths, cuts, spacing, tolerance, intervals, max_samples, True)


def _sample_path_ranges(paths, cuts, spacing, tolerance, intervals, max_samples, piecewise):
    from .paths import PathSampleBudgetExceeded
    if type(max_samples) is not int or max_samples<1:
        raise ValueError('Source-cycle sampling needs a positive finite allocation')
    result=_CanonicalSamples(complete_source_ranges=intervals is None);allocated=0
    for sid,path in paths.items():
        active=[]
        requested=((0.,path.length),) if intervals is None else intervals.get(sid,())
        for first,last in sorted(requested):
            if not (math.isfinite(first) and math.isfinite(last) and 0.<=first<=last<=path.length):
                raise ValueError('Source-cycle sampling interval leaves its original finite Path')
            if last<=first:continue
            if active and first<=active[-1][1]:active[-1]=(active[-1][0],max(last,active[-1][1]))
            else:active.append((float(first),float(last)))
        compulsory=[]
        for value in sorted({0.,path.length,*cuts.get(sid,()),*(v for interval in active for v in interval)}):
            if not compulsory or value-compulsory[-1]>1e-9:compulsory.append(value)
            elif value==path.length:compulsory[-1]=value
        if allocated+len(compulsory)>max_samples:
            raise PathSampleBudgetExceeded('Source-cycle sampling exceeds its allocation budget')
        values=[]
        for first,last in zip(compulsory,compulsory[1:]):
            if last<=first+1e-10:continue
            remaining=max_samples-allocated-len(values)
            if any(first<end and last>start for start,end in active):
                if remaining<3:raise PathSampleBudgetExceeded('Source-cycle sampling exceeds its allocation budget')
                values.extend(path.samples(first,last,min(.5,spacing),tolerance=tolerance,
                    max_samples=remaining, local_curvature=piecewise, piecewise=piecewise)[:-1])
            else:
                if remaining<2:raise PathSampleBudgetExceeded('Source-cycle sampling exceeds its allocation budget')
                values.append(first)
        values=np.r_[values,path.length]
        result[sid]=values[np.r_[True,np.diff(values)>1e-9]]
        allocated+=len(values)
    return result


def _port_footprint(owner, approach, station, radial, spacing, tolerance, partitions):
    at=[0.,station]
    for i,piece in enumerate(approach.path.pieces):
        begin=float(approach.path.stations[i]);finish=min(station,float(approach.path.stations[i+1]))
        if finish<=begin:break
        samples=partitions[piece.section.id]
        if not piece.forward:samples=piece.path.length-samples[::-1]
        at.extend(begin+float(s) for s in samples if 1e-9<s<finish-begin-1e-9)
    return _port_at_stations(owner,approach,station,radial,tolerance,np.array(sorted(set(at))))


def _port_at_stations(owner,approach,station,radial,tolerance,at):
    """Build the two banks and complete mouth from the same outward station grid."""
    at=np.asarray(at,dtype=float)
    if (at.ndim!=1 or len(at)<2 or not np.all(np.isfinite(at))
            or at[0]!=0. or at[-1]!=station or np.any(np.diff(at)<=0.)):
        raise ValueError('A port station grid must cover its complete finite outward interval')
    right=approach.point(at,approach.offset+.5*approach.width)
    left=approach.point(at,approach.offset-.5*approach.width)[::-1]
    u=np.linspace(1.,0.,2*radial+1)
    mouth=approach.point(station,approach.offset+approach.width*(u-.5))
    cap=np.vstack((left[-1],right[0]))
    curves=tuple(PolylineCurb(points) for points in (right,left,cap))
    from .source_curb_coordinates import CurbCoordinates,CurbCoordinatesUnavailable
    try:
        declarations=(
            CurbCoordinates(approach.path,at,np.full(len(at),approach.offset+.5*approach.width),
                            np.ones(len(at)),approach=approach,maximum_chord_error=tolerance,evaluation_grid=at),
            CurbCoordinates(approach.path,at[::-1],np.full(len(at),approach.offset-.5*approach.width),
                            -np.ones(len(at)),approach=approach,maximum_chord_error=tolerance,evaluation_grid=at),
            CurbCoordinates(approach.path,(0.,0.),
                            (approach.offset-.5*approach.width,approach.offset+.5*approach.width),
                            (-1,1),approach=approach,maximum_chord_error=tolerance,evaluation_grid=at))
    except CurbCoordinatesUnavailable as unavailable:
        for curve in curves:curve.source_coordinates_unavailable=str(unavailable)
    else:
        for curve,coordinates in zip(curves,declarations):curve.source_coordinates=coordinates
    points=np.vstack((right[:-1],mouth[:-1],left[:-1],cap[:-1]))
    labels=[]
    for index,curve in enumerate(curves):
        if index==1:
            labels.extend(RimSegment('mouth',0,float(a),float(b)) for a,b in zip(u,u[1:]))
        labels.extend(RimSegment('internal' if index==2 else 'curb',index,float(a/curve.length),float(b/curve.length))
                      for a,b in zip(curve.stations,curve.stations[1:]))
    layout=SimpleNamespace(approaches=(approach,),outward_stations=(station,),boundaries=curves,
        rows=len(points),radial=radial,mouths={},corners=(),
        original_root_caps={len(labels)-1:End(approach.source.id,approach.end)})
    return SampledFootprint(Region(owner,points),tuple(labels),layout)


def source_cycle_group(planner, members, consumed, internal, owners, seams, **retries):
    return drain_steps(source_cycle_group_steps(planner, members, consumed, internal, owners, seams, **retries))


def source_cycle_group_steps(planner,members,consumed,internal,owners,seams,*,retry_internal=(),retry_dimensional=(),retry_roots=(),retry_support=(),retry_median=(),retry_inner=(),retry_owned_connectors=()):
    """Union finite source pavement and local fillets before chart construction.

    Each original local contact is solved independently. A short internal road
    can be entirely owned while its two endpoint fillets overlap; their union
    clips that overlap instead of demanding an impossible retained curb span.
    """
    from .corner_radii import CornerContactError,CornerRadiusPolicy
    from .intersections import IntersectionProvider,JunctionLayout,Mouth,half_intervals
    from .compounds import CompoundGroup,SharedClaim
    settings=planner.context.network.settings
    tolerance=min(.001,.05*settings.curve_tolerance)
    radial=half_intervals(settings)
    corners=[];stations={};sources=[];patch=0;unexposed=[];unexposed_partial=[];dimensional={};local_corners={}
    approaches=[]
    from .profiles import memberships
    bundled=memberships(planner.context.network)
    for key in members:
        ordered=sorted(planner.layouts[key].approaches,
            key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
        for a in ordered:
            end=End(a.source.id,a.end)
            if abs(a.offset)>1e-8:
                raise ValueError('Source-cycle regions need centred original carriageways')
            stations.setdefault(end,0.)
            if end not in internal:approaches.append(a)
        for i,a in enumerate(ordered):
            yield
            b=ordered[(i+1)%len(ordered)]
            first,last=End(a.source.id,a.end),End(b.source.id,b.end)
            ma,mb=bundled.get(first.section_id),bundled.get(last.section_id)
            if (first in internal and last in internal and ma is not None and mb is not None
                    and ma[0].id==mb[0].id and ma[1]!=mb[1]
                    and shared_bundle_sites(planner.context.network,ma[0])):
                # Opposed members share this original terminal. Their actual
                # pavement union supplies the merged entrance; an inner gore
                # circle would add a fictitious boundary inside the overlap.
                continue
            radius,mode=planner.sites[key].corner_radius(first)
            try:
                from .corner_source_support import complete_source_roots,raise_source_root_refusal
                from .compound_contact_routes import regular_contact
                def accept_source_root(candidate):
                    if ((first,last) in retry_roots and not all(regular_contact(approach,station)
                            for approach,station in zip((a,b),candidate.stations))):
                        raise CornerContactError('Circular source contact lies on an internal original Route seam')
                    return complete_source_roots(candidate,settings.spacing,tolerance,
                        support_approaches=ordered,audit_all_contacts=(first,last) in retry_roots,
                        audit_facing_halves=(first,last) in retry_roots,on_refusal=raise_source_root_refusal)
                corner=yield from CornerRadiusPolicy.solve_steps(a,b,radius,mode,
                    math.radians(settings.continuation_angle_degrees),
                    accept=accept_source_root)
                if (first,last) in retry_roots:
                    from dataclasses import replace
                    corner=replace(corner,fit_policy='SOURCE_COMPLETE_FACING_ROOT_RETRY')
            except CornerContactError as failure:
                from .compound_corner_contacts import retry_steps
                # Preserve an already covered internal sector. A remote but
                # valid fillet must not add pavement where no curb is exposed.
                # Only a failed physical coverage proof below enables that
                # original pair's bounded alternative contact search.
                eligible=(first not in internal or last not in internal or (first,last) in retry_internal
                    or (first,last) in retry_roots)
                prepared_contact=next((record for (pair,endpoint),record in retry_support
                    if pair==(first,last) and record.report.get('policy')=='FINITE_ORIGINAL_SOURCE_CONTACT_WALK_ENVELOPE'),None)
                corner=prepared_contact.corner if prepared_contact is not None else None
                if prepared_contact is None and eligible:
                    corner=yield from retry_steps(planner.context,a,b,radius,mode,set(members),consumed,
                        ordered,force_facing=(first,last) in retry_roots)
                if corner is None:
                    corner=next((prepared for pair,prepared in retry_median if pair==(first,last)),None)
                if (corner is None and (first,last) in retry_roots and mode=='STYLE'
                        and abs(a.width-b.width)<=1e-6 and not a.profile.half and not b.profile.half):
                    from .compound_root_envelopes import prepare as root_envelope
                    from dataclasses import replace
                    try:
                        original=yield from CornerRadiusPolicy.solve_steps(a,b,radius,mode,
                            math.radians(settings.continuation_angle_degrees))
                        record=root_envelope(Key('cycle_corner',patch),original,
                            replace(settings,curve_tolerance=tolerance))
                    except ValueError:
                        pass  # The original failed contact remains authoritative.
                    else:
                        corner=record.corner;corner.root_envelope=record
                        record.source.layout.corners=(corner,)
                if corner is None:
                    from dataclasses import replace
                    from .native_zero_bank_root import prepare as prepare_endpoint_root
                    record = prepare_endpoint_root(planner, key, Key('cycle_corner', patch), a, b,
                                                   radius, mode, replace(settings, curve_tolerance=tolerance))
                    if record is not None:
                        corner = record.corner
                if corner is None:
                    from .cycle_overlap_expansion import required as overlap_expansion
                    expansion=overlap_expansion(planner,key,a,b,radius,mode,members,internal,consumed)
                    if expansion is not None:
                        from .compound_preparation import CompoundExpansionRequired
                        neighbor,proof=expansion
                        records=getattr(planner,'source_overlap_expansions',{})
                        records[key,neighbor,first,last]=proof;planner.source_overlap_expansions=records
                        raise CompoundExpansionRequired((neighbor,)) from failure
                    if first not in internal or last not in internal:
                        from .compound_cycle_hidden import partial_corner_candidate
                        if mode!='STYLE' or not partial_corner_candidate(first,last,(a,b),internal,consumed):raise
                        unexposed_partial.append((first,last,str(failure)))
                        continue
                    unexposed.append((first,last,str(failure)))
                    continue
            from .owned_prefixes import contact_walks as owned_contact_walks
            corner=owned_contact_walks(planner.context,key,corner)
            from .cycle_dimensional_retry import eligible as dimensional_eligible
            if (End(corner.first.source.id,corner.first.end),End(corner.second.source.id,corner.second.end))==(first,last):
                local_corners[first,last]=(key,corner)
            if dimensional_eligible(corner):dimensional[first,last]=(key,corner)
            if (first,last) in retry_dimensional:
                from dataclasses import replace
                from .compound_dimensional_corners import prepare as dimensional_corner,DimensionalLeadExhausted
                proxy=replace(corner,center=None,continuation=True,stations=(0.,0.),
                    fit_policy='FINITE_ORIGINAL_WIDTH_CONTINUATION')
                try:
                    from .ordinary_singleton_sources import dimensional_contact_candidate
                    record=dimensional_contact_candidate(planner,key,Key('cycle_corner',patch),proxy,
                        replace(settings,curve_tolerance=tolerance))
                    if record is None:
                        record=dimensional_corner(Key('cycle_corner',patch),proxy,
                            replace(settings,curve_tolerance=tolerance),
                            original_outer_half_profiles=bool(a.profile.half or b.profile.half))
                except DimensionalLeadExhausted as failure:
                    if a.profile.half or b.profile.half:raise
                    from .compound_root_envelopes import prepare_dimensional
                    record=prepare_dimensional(Key('cycle_corner',patch),proxy,
                        replace(settings,curve_tolerance=tolerance),failure)
                    corner=record.corner;corner.root_envelope=record
                else:
                    corner=replace(proxy,stations=record.stations)
                    corner.dimensional_corner=record
                record.source.layout.corners=(corner,)
            supported=next((record for (pair,endpoint),record in retry_support if pair==(first,last)),None)
            if supported is not None:
                from dataclasses import replace
                record=replace(supported,source=replace(supported.source,
                    region=Region(Key('cycle_corner',patch),supported.source.region.xy)))
                from .generated_junction_curbs import GeneratedJunctionArea
                if isinstance(record,GeneratedJunctionArea):record.source.layout.generated_junction_area=record
                corner=record.corner;corner.root_support_envelope=record
                record.source.layout.corners=(corner,)
            corners.append(corner)
            for approach,station in zip((corner.first,corner.second),corner.stations):
                contact_end=End(approach.source.id,approach.end)
                stations[contact_end]=max(stations.get(contact_end,0.),station)
            if hasattr(corner,'dimensional_corner'):
                source=corner.dimensional_corner.source
                patch+=1
                source.layout.original_junction=key
                sources.append(source)
                continue
            if hasattr(corner,'root_envelope'):
                source=corner.root_envelope.source
                patch+=1;source.layout.original_junction=key;sources.append(source)
                sources.extend(getattr(corner.root_envelope, 'retained_sources', ()))
                continue
            if hasattr(corner,'root_support_envelope'):
                source=corner.root_support_envelope.source
                patch+=1;source.layout.original_junction=key;sources.append(source)
                continue
            if corner.center is None:continue
            if hasattr(corner,'source_contact_walks'):
                from .compound_displaced_contacts import patch as contact_patch
            else:
                from .corner_patch import corner_patch as contact_patch
            source=contact_patch(Key('cycle_corner',patch),corner,settings.spacing,tolerance);patch+=1
            source.layout.original_junction=key
            sources.append(source)
    for token, original in retry_inner:
        from dataclasses import replace
        record=replace(original,source=replace(original.source,
            region=Region(Key('cycle_corner',patch),original.source.region.xy)))
        patch+=1
        record.source.layout.original_junction=token[2]
        record.corner.dimensional_corner=record
        record.source.layout.corners=(record.corner,)
        sources.append(record.source);corners.append(record.corner)
        for approach,station in zip((record.corner.first,record.corner.second),record.stations):
            end=End(approach.source.id,approach.end)
            stations[end]=max(stations[end],station)
    from .source_owned_connector_fillets import insert as insert_owned_connectors
    insert_owned_connectors(planner,members,consumed,sources,stations,retry_owned_connectors,settings.spacing,tolerance)
    paths={sid:planner.context[Key('path',sid)] for sid in consumed}
    paths.update((piece.section.id,piece.path) for a in approaches for piece in a.path.pieces)
    contact_cuts={};contact_ports={};direct_contacts=[]
    for corner in corners:
        for a,distance in zip((corner.first,corner.second),corner.stations):
            if distance<0.:continue
            direct_contacts.append((a,distance))
            for (sid,end),at in a.path.trims(distance).items():
                contact_cuts.setdefault(sid,set()).add(at)
                contact_ports.setdefault((sid,end),set()).add(at)
        envelope=getattr(corner,'source_contact_envelope',None)
        if envelope is not None:
            for sid,values in envelope.report['source_phase_parameters'].items():
                contact_cuts.setdefault(sid,set()).update(values)
    from .source_owned_connector_fillets import contact_cuts as owned_connector_cuts
    owned_connector_cuts(contact_cuts,retry_owned_connectors)
    for corner in corners:
        for retained in getattr(getattr(corner, 'root_envelope', None), 'retained_sources', ()):
            for original in retained.layout.corners:
                for arm, station in zip((original.first, original.second), original.stations):
                    for (sid, end), native in arm.path.trims(max(0., station)).items():
                        contact_cuts.setdefault(sid, set()).add(native)
    base_sources=tuple(sources)
    from .compound_contact_routes import outward_station
    outward=[IntersectionProvider.mouth_station(a,max(0.,stations[End(a.source.id,a.end)],outward_station(a,contact_ports,direct_contacts=direct_contacts))
        +min(.25,.05*a.width),settings) for a in approaches]
    from .ordinary_singleton_sources import OriginalSingletonRootBounds
    bounds = getattr(planner, 'ordinary_singleton_root_bounds', None)
    if type(bounds) is OriginalSingletonRootBounds and tuple(members) == (bounds.site,):
        outward = bounds.preserve_mouths(approaches, outward)
    from .cycle_terminal_ports import discover as terminal_pairs,prepare_requested,footprint as terminal_footprint
    paired=terminal_pairs(planner,members,consumed,approaches)
    separated_terminals=[]
    if paired:
        original_outward={End(a.source.id,a.end):s for a,s in zip(approaches,outward)}
        prepared=[];selected=[]
        for pair in paired:
            approach,station,gap=prepare_requested(pair,planner.context,contact_cuts,original_outward)
            if gap is not None:separated_terminals.append(gap.report())
            else:prepared.append((approach,station));selected.append(pair)
        paired=tuple(selected)
        replaced={end for pair in paired for end in pair.endpoints}
        kept=[(a,s) for a,s in zip(approaches,outward) if End(a.source.id,a.end) not in replaced]
        kept.extend(prepared)
        approaches=[a for a,s in kept];outward=[s for a,s in kept]
    from .compound_preparation import DeferredCompoundFootprint,CompoundExpansionRequired
    from .compounds import rim_label
    def check_extent(a,station):
        if IntersectionProvider.mouth_fits(a,station):return
        if hasattr(a,'original_approaches'):
            raise ValueError('Merged terminal profile exhausts its certified finite original prefixes')
        route=planner.topology.trace(End(a.source.id,a.end))
        if (route.terminal_site in planner.sites and route.terminal_site not in members
                and not route.terminal_scope and not route.terminal_selection
                and tuple(p.section.id for p in a.path.pieces)==tuple(v.section_id for v in route.traversals)):
            raise CompoundExpansionRequired((route.terminal_site,))
        raise ValueError('Source-cycle external fillets exhaust their finite approach')
    finite_profile_search=[];finite_profile_candidates=0;source_mouth_clearance=[]
    profile_lead_advances=[]
    for attempt in range(64):
        yield
        mouths={};ports=[];cuts={sid:set(values) for sid,values in contact_cuts.items()}
        for a,station in zip(approaches,outward):
            check_extent(a,station)
            actual=({(alias.endpoint.section_id,alias.endpoint.end):alias.station for alias in a.source_mouths(station)}
                    if hasattr(a,'original_approaches') else a.path.trims(station))
            for identity,at in actual.items():
                mouths[identity]=Mouth(float(at),2*radial)
                cuts.setdefault(identity[0],set()).add(at)
        # These same grids also certify future mouth candidates and the whole
        # finite remainder. Sparse owned-prefix rendering rows cannot support
        # those proofs, even when their first/last stations remain 0 and L.
        partitions=_canonical_samples(paths,cuts,settings.spacing,tolerance)
        strips=[]
        for sid in sorted(consumed):
            path=paths[sid];at=partitions[sid];section=path.section
            strip_points=np.vstack((
                path.point(at,section.offset+.5*section.width),
                path.point(at,section.offset-.5*section.width)[::-1]))
            from .source_curb_coordinates import CurbCoordinates
            coordinates=CurbCoordinates(path,np.r_[at,at[::-1]],
                np.r_[np.full(len(at),section.offset+.5*section.width),
                      np.full(len(at),section.offset-.5*section.width)],
                np.r_[np.ones(len(at)),-np.ones(len(at))],maximum_chord_error=tolerance,evaluation_grid=at)
            source=_footprint(Key('cycle_strip',sid),strip_points,
                internal_edges=(len(at)-1,2*len(at)-1),source_coordinates=coordinates)
            if source is None:raise ValueError('A source-cycle carriageway has no physical region')
            cap_ends={tuple(sorted((tuple(strip_points[i]),tuple(strip_points[(i+1)%len(strip_points)])))):End(sid,end)
                for i,end in ((len(at)-1,1),(2*len(at)-1,0))}
            source.layout.original_root_caps={i:cap_ends[token] for i,(a,b) in enumerate(zip(source.region.xy,np.roll(source.region.xy,-1,axis=0)))
                if (token:=tuple(sorted((tuple(a),tuple(b))))) in cap_ends}
            strips.append(source)
        for a,station in zip(approaches,outward):
            if hasattr(a,'original_approaches'):
                ports.append(terminal_footprint(Key('cycle_terminal',2*a.bundle_id+a.end_index),a,station,radial))
            else:
                source=_port_footprint(Key('cycle_port',2*a.source.id+a.end),a,
                    station,radial,settings.spacing,tolerance,partitions)
                if retry_support:
                    from .source_root_templates import refine_requested_port
                    source=refine_requested_port(planner.context,source,tuple(record for token,record in retry_support))
                ports.append(source)
        neighbors=set()
        for key,layout in planner.layouts.items():
            if key in members or isinstance(layout,DeferredCompoundFootprint):continue
            for (sid,end),mouth in mouths.items():
                opposite=layout.mouths.get((sid,1-end))
                if opposite is not None and (mouth.station>opposite.station+1e-7 if end==0
                                             else mouth.station<opposite.station-1e-7):neighbors.add(key)
        if neighbors:raise CompoundExpansionRequired(neighbors)
        sources=(*base_sources,*strips,*ports)
        if sum(len(source.region.xy) for source in sources)>20000:
            raise ValueError('Source-cycle physical regions exceed their sampling budget')
        # The source-cycle branch has an explicit 20k sampling allocation;
        # retain that same bound when the arrangement indexes those samples.
        # Pair-test, epsilon and geometric qualification limits stay shared.
        shared_bindings,shared_signature=(),()
        try:
            loops=yield from RegionUnion(max_segments=20000).build_steps(source.region for source in sources)
        except BoundaryIncidenceError:
            from .shared_bank_sampling import normalize
            sources,partitions,shared_bindings,shared_signature=normalize(planner.context,sources,partitions,settings)
            if not shared_signature:raise
            loops=yield from RegionUnion(max_segments=20000).build_steps(source.region for source in sources)
            inventory={s.region.owner:s for s in sources}
            ports=[inventory[s.region.owner] for s in ports]
            strips=[inventory[s.region.owner] for s in strips]
        if sum(not loop.hole for loop in loops)!=1:
            raise ValueError('Source-cycle pavement needs one connected physical exterior')
        covered_corners=[]
        if unexposed:
            from .compound_cycle_hidden import certify_hidden_corners,HiddenCornerExposureError
            by_end={End(a.source.id,a.end):a for key in members for a in planner.layouts[key].approaches}
            try:
                covered_corners=certify_hidden_corners(unexposed,by_end,partitions,loops,consumed)
            except HiddenCornerExposureError as failure:
                pair=(failure.first,failure.last)
                if pair in retry_internal:
                    if any(pair==old for old,corner in retry_median):raise
                    from .thin_median_corners import prepare as thin_median
                    try:
                        prepared=thin_median(planner,pair,consumed,members,exhausted=retry_internal)
                    except ValueError:
                        raise failure
                    if len(retry_median)>=16:
                        raise ValueError('Source-cycle thin median retry exceeds its finite pair budget') from failure
                    return (yield from source_cycle_group_steps(planner,members,consumed,internal,owners,seams,
                        retry_internal=retry_internal,retry_dimensional=retry_dimensional,retry_roots=retry_roots,
                        retry_support=retry_support,retry_median=(*retry_median,(pair,prepared)),retry_inner=retry_inner,retry_owned_connectors=retry_owned_connectors))
                if len(retry_internal)>=32:
                    raise ValueError('Source-cycle exposed-corner retry exceeds its bounded search budget') from failure
                return (yield from source_cycle_group_steps(planner,members,consumed,internal,owners,seams,
                    retry_internal=(*retry_internal,pair),retry_dimensional=retry_dimensional,retry_roots=retry_roots,
                    retry_support=retry_support,retry_median=retry_median,retry_inner=retry_inner,retry_owned_connectors=retry_owned_connectors))
        inventory={source.region.owner:source for source in sources}
        exposed={source.region.owner:0. for source in ports}
        for loop in loops:
            for span in loop.spans:
                label=rim_label(span,inventory)
                if label.kind=='mouth':exposed[span.owner]+=label.first-label.last
        missing=[i for i,source in enumerate(ports) if abs(exposed[source.region.owner]-1.)>1e-8]
        if not missing:
            from .source_mouth_clearance import advance_past_overlap
            advanced,records=advance_past_overlap(planner,approaches,outward,ports,sources,partitions,members,
                max_candidates=4096-finite_profile_candidates)
            if not records:break
            if attempt==63:
                raise ValueError('Source-cycle remaining-strip clearance exceeds its finite lead budget')
            finite_profile_candidates+=records[-1]['candidates_tested']
            for record in records:record['normal_union_iteration']=attempt+1
            source_mouth_clearance.extend(records);outward=advanced
            continue
        if attempt==63:
            raise ValueError('Source-cycle complete-profile search exceeds its finite lead budget')
        if attempt>=31:
            from .finite_profile_search import advance,FiniteProfileUnavailable
            if finite_profile_candidates>=4096:
                raise ValueError('Source-cycle finite profile search exceeds its bounded candidate budget')
            previous=outward
            try:
                outward,records=advance(planner,approaches,outward,missing,ports,sources,partitions,members,
                    max_candidates=4096-finite_profile_candidates,full_source_exclusion=attempt>=32)
            except FiniteProfileUnavailable:
                from .source_mouth_clearance import advance_past_overlap
                try:
                    # An actually owned prefix can cover every remaining cut
                    # on another original road. Only its proven finite-end
                    # neighbor may expand this pool; proof-only future strips
                    # and returned alternative stations cannot authorize it.
                    advance_past_overlap(planner,[approaches[i] for i in missing],
                        [outward[i] for i in missing],[ports[i] for i in missing],sources,
                        partitions,members,max_candidates=0)
                except CompoundExpansionRequired:
                    raise
                except ValueError:
                    pass
                raise
            if any(outward[i]!=previous[i] for i in range(len(outward)) if i not in missing):
                raise ValueError('Finite profile closure changed an already complete source mouth')
            finite_profile_candidates+=records[-1]['candidates_tested']
            for record,index in zip(records,missing):
                record.update(normal_union_iteration=attempt+1,
                    actual_exposed_profile_measure=exposed[ports[index].region.owner],
                    other_complete_mouths_unchanged=True)
            finite_profile_search.extend(records)
            continue
        # Extend only real original curb leads. Circle contacts, widths,
        # radii, source controls and transverse mouth partitions stay fixed.
        for i in missing:
            previous=outward[i]
            outward[i]=IntersectionProvider.mouth_station(approaches[i],
                previous+max(.25,settings.spacing),settings)
            from .incomplete_profile_leads import capture as capture_profile_lead
            evidence=capture_profile_lead(approaches[i],previous,outward[i],ports[i],sources,loops,attempt+1)
            if evidence is not None:profile_lead_advances.append(evidence)
    else:raise ValueError('Source-cycle complete-profile search exceeds its finite lead budget')
    for corner,record in shared_bindings:
        corner.dimensional_corner=record
    covered_partial=[]
    if unexposed_partial:
        from .compound_cycle_hidden import certify_partially_hidden_corners
        by_end={End(a.source.id,a.end):a for key in members for a in planner.layouts[key].approaches}
        covered_partial=certify_partially_hidden_corners(unexposed_partial,by_end,partitions,loops,consumed,mouths)
    from .bundles import CenterlineRange,certify_disjoint
    proofs={}
    for bundle in planner.context.network.bundles:
        if not set(bundle.sections)<=consumed or not shared_bundle_sites(planner.context.network,bundle):continue
        ranges=[]
        for member in bundle.members:
            active=[]
            for i,visit in enumerate(member.path):
                path=planner.context[Key('path',visit.section_id)]
                first,last=0.,path.length
                for terminal,end in ((i==0,0 if visit.forward else 1),
                                      (i==len(member.path)-1,1 if visit.forward else 0)):
                    if not terminal:continue
                    trim=max(0.,stations[End(visit.section_id,end)])+min(.25,.05*path.section.width)
                    if end==0:first=max(first,trim)
                    else:last=min(last,path.length-trim)
                if last<=first+1e-7:
                    raise ValueError('Shared terminal consumes a complete Bundle proof span')
                active.append(CenterlineRange.build(path,first,last))
            ranges.append(tuple(active))
        proofs[bundle.id]=certify_disjoint(*ranges)
    for sid in sorted(consumed):
        length=planner.context[Key('path',sid)].length
        mouths[sid,0]=Mouth(length,2*radial);mouths[sid,1]=Mouth(0.,2*radial)
    physical=JunctionLayout(tuple(approaches),tuple(corners),(),tuple(outward),mouths,
        max(4,len(loops[0].xy)),radial,np.zeros(2),[0.]*len(approaches),
        mesh_strategy='boundary_chart',absorbed_transitions=tuple(sorted(seams)))
    physical.source_cycle_proofs=proofs
    physical.shared_bank_signature=shared_signature
    physical.source_cycle_profile_lead_advances=tuple(profile_lead_advances)
    physical.source_cycle_unexposed_corners=unexposed
    physical.source_cycle_covered_corners=covered_corners
    physical.source_cycle_covered_partial_corners=covered_partial
    physical.source_cycle_displaced_contacts=[c.source_contact_support for c in corners if hasattr(c,'source_contact_support')]
    physical.source_cycle_dimensional_contacts=[c.dimensional_corner.report for c in corners if hasattr(c,'dimensional_corner')]
    physical.source_cycle_root_envelopes=[c.root_envelope.report for c in corners if hasattr(c,'root_envelope')]
    physical.source_cycle_root_support_envelopes=[c.root_support_envelope.report for c in corners if hasattr(c,'root_support_envelope')]
    physical.source_cycle_source_contact_envelopes=[c.source_contact_envelope.report for c in corners if hasattr(c,'source_contact_envelope')]
    physical.source_cycle_thin_median_gores=[c.thin_median_gore.report for c in corners if hasattr(c,'thin_median_gore')]
    physical.source_cycle_dimensional_retries=tuple(retry_dimensional)
    physical.source_cycle_inner_bundle_width_support=[record.report for token,record in retry_inner]
    physical.source_cycle_root_retries=tuple(retry_roots)
    physical.source_cycle_contact_retries=tuple(retry_internal)
    physical.source_cycle_exposed_root_contact_retries=[proof for (token,pair),proof in
        getattr(planner,'source_exposed_root_retries',{}).items() if token==tuple(members) and pair in retry_internal]
    physical.source_cycle_terminal_ports=[source.layout.terminal_merged_port for source in ports
        if hasattr(source.layout,'terminal_merged_port')]
    physical.source_cycle_separated_terminal_cuts=separated_terminals
    physical.source_cycle_finite_profile_search=finite_profile_search
    physical.source_cycle_mouth_clearance=source_mouth_clearance
    physical.source_cycle_mouth_expansions=[proof for proof in getattr(planner,'source_mouth_expansions',{}).values()
        if set(proof['previous_members'])<={key.text() for key in members}
        and proof['original_neighbor'] in {key.text() for key in members}]
    physical.source_cycle_overlap_expansions=[proof for proof in getattr(planner,'source_overlap_expansions',{}).values()
        if proof['junction'] in {key.text() for key in members} and proof['terminal_junction'] in {key.text() for key in members}]
    from .compound_ports import CompoundPortPlanner
    group=CompoundGroup(members,tuple(SharedClaim(sid,owners[sid][0],owners[sid][1],
        planner.context[Key('path',sid)].length,0.) for sid in sorted(consumed)),tuple(sources),loops,physical)
    from .native_zero_bank_root import validate_group as validate_endpoint_roots
    group = validate_endpoint_roots(planner.context, group)
    inventory={source.region.owner:source for source in sources}
    if any(rim_label(span,inventory).kind=='internal' for loop in loops for span in loop.spans):
        from .cycle_root_support import prepare_required_steps as root_support
        supported=yield from root_support(planner,loops,sources,corners,tuple(token for token,record in retry_support),consumed,members,settings,
            retry_roots=retry_roots,retry_dimensional=retry_dimensional,source_group=group)
        if supported is not None:
            if len(retry_support)>=16:raise ValueError('Source-cycle adjacent root support exceeds its finite pair budget')
            return (yield from source_cycle_group_steps(planner,members,consumed,internal,owners,seams,
                retry_internal=retry_internal,retry_dimensional=retry_dimensional,retry_roots=retry_roots,
                retry_support=(*retry_support,supported),retry_median=retry_median,retry_inner=retry_inner,retry_owned_connectors=retry_owned_connectors))
        from .cycle_dimensional_retry import required as dimensional_retry,circular_required
        root_pair=circular_required(planner,loops,sources,local_corners,retry_roots,settings)
        if root_pair is not None:
            if len(retry_roots)>=16:raise ValueError('Source-cycle original-root retry exceeds its finite pair budget')
            try:
                return (yield from source_cycle_group_steps(planner,members,consumed,internal,owners,seams,
                    retry_internal=retry_internal,retry_dimensional=retry_dimensional,retry_roots=(*retry_roots,root_pair),
                    retry_support=retry_support,retry_median=retry_median,retry_inner=retry_inner,retry_owned_connectors=retry_owned_connectors))
            except CornerContactError:
                from .source_contact_envelopes import prepare_required as contact_envelopes
                from dataclasses import replace
                prepared=contact_envelopes(planner,loops,sources,corners,unexposed,consumed,members,
                    replace(settings,curve_tolerance=tolerance))
                if prepared is None:raise
                if len(retry_support)+len(prepared)>16:
                    raise ValueError('Source-cycle source contact envelopes exceed their finite pair budget')
                return (yield from source_cycle_group_steps(planner,members,consumed,internal,owners,seams,
                    retry_internal=retry_internal,retry_dimensional=retry_dimensional,retry_roots=retry_roots,
                    retry_support=(*retry_support,*prepared),retry_median=retry_median,retry_inner=retry_inner,retry_owned_connectors=retry_owned_connectors))
        pair=dimensional_retry(planner,loops,sources,dimensional,retry_dimensional)
        if pair is not None:
            if len(retry_dimensional)>=16:raise ValueError('Source-cycle dimensional retry exceeds its finite pair budget')
            return (yield from source_cycle_group_steps(planner,members,consumed,internal,owners,seams,
                retry_internal=retry_internal,retry_dimensional=(*retry_dimensional,pair),retry_roots=retry_roots,
                retry_support=retry_support,retry_median=retry_median,retry_inner=retry_inner,retry_owned_connectors=retry_owned_connectors))
        from .cycle_exposed_root_contacts import required as exposed_root_contact
        required=exposed_root_contact(planner,loops,sources,unexposed,internal,consumed,retry_internal)
        if required is not None:
            pair,proof=required
            if len(retry_internal)>=32:
                raise ValueError('Source-cycle exposed-corner retry exceeds its bounded search budget')
            records=getattr(planner,'source_exposed_root_retries',{})
            records[tuple(members),pair]=proof;planner.source_exposed_root_retries=records
            return (yield from source_cycle_group_steps(planner,members,consumed,internal,owners,seams,
                retry_internal=(*retry_internal,pair),retry_dimensional=retry_dimensional,retry_roots=retry_roots,
                retry_support=retry_support,retry_median=retry_median,retry_inner=retry_inner,retry_owned_connectors=retry_owned_connectors))
        from .inner_bundle_width_support import required as inner_width_support
        from dataclasses import replace
        inner=inner_width_support(planner,loops,sources,tuple(token for token,record in retry_inner),
            consumed,members,replace(settings,curve_tolerance=tolerance))
        if inner is not None:
            if len(retry_inner)>=16:
                raise ValueError('Source-cycle inner Bundle support exceeds its finite pair budget')
            return (yield from source_cycle_group_steps(planner,members,consumed,internal,owners,seams,
                retry_internal=retry_internal,retry_dimensional=retry_dimensional,retry_roots=retry_roots,
                retry_support=retry_support,retry_median=retry_median,retry_inner=(*retry_inner,inner),retry_owned_connectors=retry_owned_connectors))
        from .native_owned_cap_chains import prepare as prepare_cap_chain
        repaired = prepare_cap_chain(planner, group)
        if repaired is None:
            raise ExposedPreparationSeam('Source-cycle union exposes an internal fillet preparation seam')
        group = repaired
        physical, loops, sources = group.physical_layout, group.loops, list(group.footprints)
        inventory = {source.region.owner: source for source in sources}
    for source in sources:
        if not any(getattr(getattr(c,'root_support_envelope',None),'report',{}).get('admission')
                =='EXPOSED_CYCLIC_ORIGINAL_ZERO_HEIGHT_ROOT' or
                getattr(getattr(c,'source_contact_envelope',None),'report',{}).get('policy')
                =='FINITE_ORIGINAL_SOURCE_CONTACT_WALK_ENVELOPE' for c in source.layout.corners):continue
        amount=sum(abs(rim_label(span,inventory).last-rim_label(span,inventory).first)
            for loop in loops for span in loop.spans if not loop.hole and span.owner==source.region.owner
            and rim_label(span,inventory).kind=='curb')
        if abs(amount-1.)>1e-8:
            raise ValueError('Full source union clips a cyclic root support curve before its original tangent endpoints')
    exposed={span.owner for loop in loops for span in loop.spans if rim_label(span,inventory).kind=='mouth'}
    expected={source.region.owner for source in sources if source.layout.approaches}
    if exposed!=expected:
        raise ValueError('Source-cycle union hides an original incoming profile')
    if retry_median:
        from .thin_median_corners import certify as certify_thin_median
        certify_thin_median(corners,sources,loops)
    if CompoundPortPlanner.incomplete_ports(group):
        raise ValueError('Source-cycle union clips an original incoming profile')
    CompoundPortPlanner(planner.context,tuple(planner.sites.values()),{}).check_support(consumed,loops[0])
    if paired or separated_terminals or covered_partial or retry_dimensional or retry_roots or retry_support or retry_median or retry_inner or finite_profile_search or source_mouth_clearance or any(hasattr(c,'source_contact_walks') for c in corners):
        from .cycle_dimensional_retry import certify_layout_steps
        physical.source_cycle_full_width_proof=yield from certify_layout_steps(planner.context,group,consumed)
    if paired and len(group.loops)==1:
        from .source_rim_normalization import normalize
        group=normalize(group)
        group.physical_layout.source_cycle_full_width_proof=yield from certify_layout_steps(planner.context,group,consumed)
    from .source_owned_connector_fillets import bind_contacts,required as required_owned_connector
    group=bind_contacts(planner.context,group)
    reservation=required_owned_connector(planner,group,retry_owned_connectors)
    if reservation is None:
        from .source_bank_preparation import UnfittedBankSwitch
        from .source_bank_reconciliation import reconcile_steps as reconcile_source_banks
        yield
        try:
            group=yield from reconcile_source_banks(planner,group,partitions)
        except UnfittedBankSwitch as failure:
            from .successor_bank_contacts import reserve
            reservation=reserve(planner,failure.group,failure.proposal)
            if reservation is None or any(r.token == reservation.token for r in retry_owned_connectors):
                raise
    if reservation is not None:
        if len(retry_owned_connectors)>=16:
            raise ValueError('Source-cycle owned connector repair exceeds its finite reservation budget')
        return (yield from source_cycle_group_steps(planner,members,consumed,internal,owners,seams,
            retry_internal=retry_internal,retry_dimensional=retry_dimensional,retry_roots=retry_roots,
            retry_support=retry_support,retry_median=retry_median,retry_inner=retry_inner,
            retry_owned_connectors=(*retry_owned_connectors,reservation)))
    from .source_parent_chain_fillets import reconcile as reconcile_parent_chains,validate as validate_parent_chains
    yield
    group=reconcile_parent_chains(planner,group)
    from .residual_road_gaps import reconcile as reconcile_residual_gaps
    yield
    group=reconcile_residual_gaps(planner,group)
    from .retained_arc_fillets import reconcile_steps as reconcile_retained_arcs
    yield
    group=yield from reconcile_retained_arcs(planner,group)
    from .arc_lead_joins import reconcile_steps as reconcile_arc_leads
    yield
    group=yield from reconcile_arc_leads(planner,group)
    from .exposed_bank_joins import reconcile_steps as reconcile_bank_joins
    yield
    group=yield from reconcile_bank_joins(planner,group)
    from .bank_circle_preparation import validate as validate_bank_circles
    validate_bank_circles(group)
    from .successor_bank_contacts import validate_exposure
    validate_exposure(group)
    return validate_parent_chains(validate_endpoint_roots(planner.context, group))
