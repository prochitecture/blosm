"""Source-grid refinement and unchanged-owner guards for finite parent chains."""
from copy import copy
from dataclasses import replace
from types import SimpleNamespace
import math
import numpy as np
from .model import Key,End
from .compounds import SampledFootprint,RimSegment,rim_label
from .regions import Region,RegionUnion,signed_area
from .finite_free_caps import PreparedCapCurb
from .source_curb_coordinates import CurbCoordinates
from .source_parent_chain_fillets import (SOURCE_PARENT_CHAIN_POLICY,ParentChainUnavailable,
    CertifiedCoreContact,OriginalCircleChain,PreparedParentChainFillet,discover,fit,validate)
from .source_owned_connector_fillets import NativeCurbContact


def _refine_core(source,index,parameter):
    curve=source.layout.boundaries[index];points=[];labels=[];added=False
    for at,label in enumerate(source.segments):
        points.append(source.region.xy[at])
        if label.kind=='curb' and label.index==index and label.first<parameter<label.last:
            if added:raise ParentChainUnavailable('A core contact has ambiguous original parameter ancestry')
            labels.append(replace(label,last=parameter));points.append(curve(parameter))
            labels.append(replace(label,first=parameter));added=True
        else:labels.append(label)
    if not added:raise ParentChainUnavailable('The finite core contact is not a new original curve cut')
    points=np.asarray(points);layout=copy(source.layout);curves=list(layout.boundaries);cloned=copy(curve)
    indices=[i for i,label in enumerate(labels) if label.kind=='curb' and label.index==index]
    refs=np.asarray([points[i] for i in indices]+[points[(indices[-1]+1)%len(points)]])
    refs.flags.writeable=False;cloned.contact_reference_points=refs;curves[index]=cloned
    layout.boundaries=tuple(curves);layout.rows=len(points)
    if not {tuple(p) for p in source.region.xy}<={tuple(p) for p in points}:
        raise ParentChainUnavailable('Parent-chain refinement moved an original core vertex')
    region=Region(source.region.owner,points);RegionUnion().build((region,))
    return SampledFootprint(region,tuple(labels),layout),cloned,curve(parameter)


def _refine_strip(source,path,station,side):
    from .compound_cycles import _footprint
    ledgers=[getattr(curve,'source_coordinates',None) for curve in source.layout.boundaries]
    if not ledgers or any(type(data) is not CurbCoordinates or data.path is not path or data.approach is not None for data in ledgers):
        raise ParentChainUnavailable('Parent-chain bank lacks its exact consumed native-strip coordinates')
    if any(data.maximum_chord_error is None for data in ledgers):
        raise ParentChainUnavailable('Parent-chain bank lacks its declared original sampling budget')
    grid=np.unique(np.concatenate([data.stations for data in ledgers]+[np.asarray((station,))]));section=path.section;count=len(grid)
    tolerance=min(data.maximum_chord_error for data in ledgers)
    right=path.point(grid,section.offset+.5*section.width);left=path.point(grid,section.offset-.5*section.width)
    points=np.vstack((right,left[::-1]))
    coordinates=CurbCoordinates(path,np.r_[grid,grid[::-1]],
        np.r_[np.full(count,section.offset+.5*section.width),np.full(count,section.offset-.5*section.width)],
        np.r_[np.ones(count),-np.ones(count)],maximum_chord_error=tolerance,evaluation_grid=grid)
    result=_footprint(source.region.owner,points,internal_edges=(count-1,2*count-1),source_coordinates=coordinates)
    caps={tuple(sorted((tuple(points[i]),tuple(points[(i+1)%len(points)])))):End(section.id,end)
        for i,end in ((count-1,1),(2*count-1,0))}
    result.layout.original_root_caps={i:caps[token] for i,(a,b) in enumerate(zip(result.region.xy,np.roll(result.region.xy,-1,axis=0)))
        if (token:=tuple(sorted((tuple(a),tuple(b))))) in caps}
    for old in source.layout.boundaries:
        if not {tuple(point) for point in old.points}<={tuple(point) for curve in result.layout.boundaries for point in curve.points}:
            raise ParentChainUnavailable('Parent-chain refinement moved an original native source vertex')
    at=int(np.flatnonzero(grid==station)[0]);point=right[at] if side==1 else left[at]
    return result,point,grid


def _contact_node(loop,sources,owner,index,parameter,point):
    matches=[]
    for at,span in enumerate(loop.spans):
        label=rim_label(span,sources)
        if span.owner==owner and label.kind=='curb' and label.index==index and abs(label.first-parameter)<1e-12:
            if np.linalg.norm(loop.xy[at]-point)>1e-10:
                raise ParentChainUnavailable('Declared source contact and its physical boundary disagree')
            matches.append(at)
    return matches


def ordered_curb_closure(loops, sources, start, finish, intermediate=()):
    """Trace the source rim between exact contacts, visiting each declared curb in order."""
    expected = (start[:2], *intermediate, finish[:2])
    matches = []
    for loop in loops:
        first = _contact_node(loop, sources, *start)
        last = _contact_node(loop, sources, *finish)
        if len(first) == len(last) == 1 and first != last:
            matches.append((loop, first[0], last[0]))
    if len(matches) != 1:
        raise ParentChainUnavailable('Finite contacts must bound one exposed source rim')
    loop, first, last = matches[0]
    indices = np.arange(first, first + (last-first) % len(loop.xy)) % len(loop.xy)
    ancestry = []
    actual = []
    for i in indices:
        span = loop.spans[i]
        label = rim_label(span, sources)
        if label.kind != 'curb':
            raise ParentChainUnavailable('A mouth or cap interrupts the declared curb chain')
        key = (span.owner, label.index)
        if not actual or actual[-1] != key:
            actual.append(key)
        ancestry.append((span.owner, label))
    if tuple(actual) != expected:
        raise ParentChainUnavailable('The exposed curbs do not follow the declared source chain')
    points = loop.xy[np.r_[indices, last]][::-1].copy()
    points[0], points[-1] = finish[-1], start[-1]
    return points, tuple(ancestry)


def _insert(context,candidate,solved,radius,attempts,regularity,owner):
    original=candidate.group;curve=candidate.core_curve
    parameter=float((curve.lengths[candidate.piece_index]+solved['poly_contact']*candidate.polynomial.length)/curve.length)
    core,refined_curve,P=_refine_core(candidate.core_source,candidate.core_index,parameter)
    bank=candidate.bank;arm=candidate.bank_arm
    native=solved['bank_outward_contact'] if arm.end==0 else arm.length-solved['bank_outward_contact']
    strip=next(area for area in original.footprints if area.region.owner==bank.owner)
    native_area,Q,grid=_refine_strip(strip,bank.path,native,bank.side)
    areas=tuple(core if area is candidate.core_source else native_area if area is strip else area for area in original.footprints)
    loops=RegionUnion(max_segments=20000).build(area.region for area in areas if not getattr(area.region,'boundary_only',False))
    sources={area.region.owner:area for area in areas};qindex=qparameter=None
    for i,bankcurve in enumerate(native_area.layout.boundaries):
        data=bankcurve.source_coordinates
        for at in np.flatnonzero(np.all(bankcurve.points==Q,axis=1)):
            if at<len(bankcurve.points)-1 and data.sides[at]==bank.side:
                if qindex is not None:raise ParentChainUnavailable('Native contact has ambiguous source curve incidence')
                qindex,qparameter=i,float(bankcurve.stations[at]/bankcurve.length)
    if qindex is None:raise ParentChainUnavailable('Native contact lost its exact prepared source vertex')
    closure,ancestry=ordered_curb_closure(loops,sources,
        (core.region.owner,candidate.core_index,parameter,P), (native_area.region.owner,qindex,qparameter,Q),
        ((candidate.circle_source.region.owner,candidate.circle_index),))
    tolerance=min(.001,.05*context.network.settings.curve_tolerance)
    count=max(2,int(math.ceil(radius*abs(solved['sweep'])/min(.25,context.network.settings.spacing,math.sqrt(8.*tolerance*radius)))))
    if count+len(closure)>8192:raise ParentChainUnavailable('Parent-chain operand exceeds its source sample budget')
    parameters=tuple(np.linspace(0.,1.,count+1))
    contact=CertifiedCoreContact(core.region.owner,candidate.core_index,refined_curve,parameter,
        candidate.core_visible_parameters,curve.corner,candidate.core_site)
    native_contact=NativeCurbContact(bank.path,bank.side,float(native),bank.intervals)
    parent=OriginalCircleChain(candidate.circle_source.region.owner,candidate.circle_index,
        candidate.circle_curve,candidate.circle_curve.original_circle_interval,candidate.circle,candidate.circle_site)
    record=PreparedParentChainFillet(owner,solved['center'],radius,solved['start_angle'],solved['sweep'],
        tuple(P),tuple(Q),contact,native_contact,parent,tuple(ancestry),candidate.radius,candidate.mode,
        attempts,parameters,tolerance,regularity)
    arc=record(np.asarray(parameters));internal=PreparedCapCurb(closure)
    polygon=np.vstack((arc[:-1],closure[:-1]))
    if signed_area(polygon)<=1e-10:raise ParentChainUnavailable('Parent-chain infill is not a positive physical area')
    labels=tuple([RimSegment('curb',0,float(a),float(b)) for a,b in zip(parameters,parameters[1:])]
        +[RimSegment('internal',1,float(a/internal.length),float(b/internal.length)) for a,b in zip(internal.stations,internal.stations[1:])])
    layout=SimpleNamespace(approaches=(),outward_stations=(),mouths={},corners=(),rows=len(polygon),
        boundaries=(record,internal),generated_junction_area=record,source_parent_chain_fillet=record)
    wedge=SampledFootprint(Region(owner,polygon),labels,layout);RegionUnion().build((wedge.region,))
    areas=(*areas,wedge)
    if sum(len(area.region.xy) for area in areas)>20000:raise ParentChainUnavailable('Parent-chain union exceeds its original segment allowance')
    physical=copy(original.physical_layout)
    result=replace(original,footprints=areas,physical_layout=physical,
        loops=RegionUnion(max_segments=20000).build(area.region for area in areas if not getattr(area.region,'boundary_only',False)))
    return result,record,grid


def _future_guard(planner,group,record=None):
    from .compound_cycles import _canonical_samples
    from .source_mouth_clearance import advance_past_overlap
    context=planner.context;paths={claim.section_id:context[Key('path',claim.section_id)] for claim in group.claims}
    for arm in group.physical_layout.approaches:
        for piece in arm.path.pieces:paths[piece.section.id]=piece.path
    cuts={sid:set() for sid in paths}
    for source in group.footprints:
        for curve in source.layout.boundaries:
            coordinates=getattr(curve,'source_coordinates',None)
            if type(coordinates) is not CurbCoordinates:continue
            for edge in range(len(coordinates.stations)-1):
                for piece in coordinates.edge_sources(edge):
                    sid=piece.path.section.id
                    if sid in cuts:cuts[sid].update((piece.first_station,piece.last_station))
    for (sid,end),mouth in group.physical_layout.mouths.items():cuts[sid].add(mouth.station)
    if record is not None:
        cuts[record.native_contact.path.section.id].add(record.native_contact.station)
    grids=_canonical_samples(paths,cuts,context.network.settings.spacing,min(.001,.05*context.network.settings.curve_tolerance))
    ports=[]
    for arm in group.physical_layout.approaches:
        matches=[area for area in group.footprints if len(getattr(area.layout,'approaches',()))==1 and area.layout.approaches[0] is arm]
        if len(matches)!=1:raise ParentChainUnavailable('Parent-chain infill lost an exact original complete port')
        ports.append(matches[0])
    outward=group.physical_layout.outward_stations
    unchanged,remote=advance_past_overlap(planner,group.physical_layout.approaches,outward,ports,
        group.area_footprints,grids,group.members,max_candidates=0)
    if tuple(unchanged)!=tuple(outward) or remote:raise ParentChainUnavailable('Parent-chain infill requires an unauthorized original mouth advance')
    for key,layout in planner.layouts.items():
        if key in group.members:continue
        for (sid,end),mouth in group.physical_layout.mouths.items():
            opposite=getattr(layout,'mouths',{}).get((sid,1-end))
            if opposite is not None and (mouth.station>opposite.station+1e-7 if end==0 else mouth.station<opposite.station-1e-7):
                raise ParentChainUnavailable('Parent-chain infill conflicts with a prepared original neighbor mouth')


def reconcile(planner,group):
    from .compound_ports import CompoundPortPlanner
    from .cycle_dimensional_retry import certify_layout
    from .sampling import AnalyticalBoundarySampler
    from .retained_arc_fillets import _topology,_mouth_points
    original=group;records=[];proofs=[];visited=set();settings=planner.context.network.settings
    if getattr(group.physical_layout,'source_parent_chain_reconciliation',None):return validate(group)
    for _ in range(16):
        candidates=discover(planner,group)
        if not candidates:
            if records:
                group.physical_layout.source_parent_chain_reconciliation=dict(policy=SOURCE_PARENT_CHAIN_POLICY,
                    records=tuple(record.report() for record in records),proofs=tuple(proofs),
                    original_logical_corners_unchanged=group.physical_layout.corners is original.physical_layout.corners,
                    prepared_owner_area_guard='PENDING_FINAL_PLAN_INVENTORY')
            return validate(group)
        candidate=candidates[0]
        token=(candidate.core_source.region.owner,candidate.core_index,candidate.circle_source.region.owner,candidate.circle_index,candidate.bank.token)
        if token in visited:raise ParentChainUnavailable('Parent-chain reconciliation repeats an original boundary chain')
        visited.add(token);solved,radius,attempts,regularity=fit(candidate)
        owner=Key('source_parent_chain_fillet',len(records))
        if owner in {area.region.owner for area in group.footprints}:raise ParentChainUnavailable('Parent-chain generated owner namespace collided')
        mouths=_mouth_points(group)
        result,record,grid=_insert(planner.context,candidate,solved,radius,attempts,regularity,owner)
        validate(result);_topology(group,result,record.source_chord_tolerance)
        before=certify_layout(planner.context,result,{claim.section_id for claim in result.claims})
        sampler=AnalyticalBoundarySampler(tolerance=settings.curve_tolerance,spacing=settings.spacing,max_samples=20000)
        result=sampler.build(result);validate(result);_topology(group,result,record.source_chord_tolerance)
        after=certify_layout(planner.context,result,{claim.section_id for claim in result.claims})
        actual=_mouth_points(result)
        if mouths.keys()!=actual.keys() or any(not np.array_equal(point,actual[key]) for key,point in mouths.items()):
            raise ParentChainUnavailable('Parent-chain infill changed an original complete profile')
        if CompoundPortPlanner.incomplete_ports(result):raise ParentChainUnavailable('Parent-chain infill clips an original complete profile')
        if result.members!=original.members or result.claims!=original.claims:
            raise ParentChainUnavailable('Parent-chain infill changed original logical ownership')
        _future_guard(planner,result,record)
        result.physical_layout.source_cycle_full_width_proof=after
        records.append(record);proofs.append(dict(raw_source_width=before,canonical_source_width=after,
            full_arc_and_hidden_closure=True,original_source_vertices_retained=True,
            complete_future_strip_guard='PASS_WITHOUT_MOUTH_ADVANCE',
            prepared_owner_area_guard='PENDING_FINAL_PLAN_INVENTORY'))
        group=result
    raise ParentChainUnavailable('Parent-chain reconciliation exceeds its finite 16-pair allowance')
