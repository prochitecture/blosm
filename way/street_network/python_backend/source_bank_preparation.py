"""Shared original sampling and admission of typed source-bank infills."""
from dataclasses import replace
from copy import copy
import numpy as np
from .model import Key,End
from .compound_cycles import _port_footprint,_footprint
from .source_curb_coordinates import CurbCoordinates


def _refine_areas(context,group,grids):
    sources=[];proof=[];tolerance=min(.001,.05*context.network.settings.curve_tolerance)
    for source in group.footprints:
        old=source
        if source.region.owner.kind=='cycle_port':
            arm,=source.layout.approaches;station,=source.layout.outward_stations
            source=_port_footprint(source.region.owner,arm,station,source.layout.radial,
                context.network.settings.spacing,tolerance,grids)
        elif source.region.owner.kind=='cycle_strip':
            sid=source.region.owner.id;path=context[Key('path',sid)];at=grids[sid];section=path.section
            points=np.vstack((path.point(at,section.offset+.5*section.width),
                              path.point(at,section.offset-.5*section.width)[::-1]))
            coordinates=CurbCoordinates(path,np.r_[at,at[::-1]],
                np.r_[np.full(len(at),section.offset+.5*section.width),
                      np.full(len(at),section.offset-.5*section.width)],
                np.r_[np.ones(len(at)),-np.ones(len(at))],maximum_chord_error=tolerance,evaluation_grid=at)
            source=_footprint(old.region.owner,points,internal_edges=(len(at)-1,2*len(at)-1),
                source_coordinates=coordinates)
            caps={tuple(sorted((tuple(points[i]),tuple(points[(i+1)%len(points)])))):End(sid,end)
                    for i,end in ((len(at)-1,1),(2*len(at)-1,0))}
            source.layout.original_root_caps={i:caps[token] for i,(a,b) in enumerate(zip(
                source.region.xy,np.roll(source.region.xy,-1,axis=0)))
                if (token:=tuple(sorted((tuple(a),tuple(b))))) in caps}
        if source is not old:
            for a,b in zip(old.layout.boundaries,source.layout.boundaries):
                if not {tuple(point) for point in a.points}<={tuple(point) for point in b.points}:
                    raise ValueError('Source-bank shared refinement lost an original physical vertex')
            proof.append(dict(owner=source.region.owner.text(),
                all_original_vertices_retained=True,source_tolerance=tolerance))
        sources.append(source)
    return replace(group,footprints=tuple(sources)),proof

from .compounds import rim_label
from .compound_ports import CompoundPortPlanner
from .cycle_dimensional_retry import certify_layout_steps
from .updates import drain_steps
from .source_curb_rounding import patch
from .source_mouth_clearance import advance_past_overlap
from .corner_radii import CornerContactError

from .source_bank_reconciliation import PreparedBankFillet,SourceBankDiscovery,fit_exposed_switch,SOURCE_BANK_POLICY,SourceBankReconciliationError,SourceBankExposureError


class UnfittedBankSwitch(CornerContactError):
    def __init__(self, group, proposal, reason):
        super().__init__(str(reason))
        self.group, self.proposal = group, proposal


def _apply_steps(planner,original,grids,proposals_and_fits):
    context=planner.context
    settings=context.network.settings;tolerance=min(.001,.05*settings.curve_tolerance)
    newgrids={sid:values.copy() for sid,values in grids.items()}
    for proposal,fitted in proposals_and_fits:
        for arm,contact,crossing in zip((fitted.corner.first,fitted.corner.second),fitted.corner.stations,fitted.intersection):
            native=[value if arm.end==0 else arm.length-value for value in (contact,crossing)]
            newgrids[arm.source.id]=np.unique(np.r_[newgrids[arm.source.id],native])
    if sum(len(values) for values in newgrids.values())>20000:
        raise ValueError('Joint source-bank refinement exceeds its unchanged original sample budget')
    group,oldvertices=_refine_areas(context,original,newgrids)
    sources=list(group.footprints);records=[];owners=[]
    for index,(proposal,fitted) in enumerate(proposals_and_fits):
        corner=fitted.corner;owner=Key('source_bank_fillet',index)
        if owner in {source.region.owner for source in sources}:raise ValueError('Generated owner namespace collided')
        partitions=tuple(newgrids[arm.source.id] if arm.end==0 else arm.length-newgrids[arm.source.id][::-1]
            for arm in (corner.first,corner.second))
        source=patch(owner,corner,fitted.intersection,settings.spacing,tolerance,partitions=partitions)
        source.layout.original_junction=proposal.site
        generated=PreparedBankFillet(proposal,fitted)
        source.layout.source_bank_fillet=generated
        source.layout.generated_junction_area=generated
        sources.append(source);records.append(generated);owners.append(owner)
    if sum(len(source.region.xy) for source in sources)>20000:
        raise ValueError('Source-bank area operands exceed their unchanged original sampling budget')
    physical=copy(group.physical_layout)
    physical.source_bank_fillets=tuple(records)
    assert physical.corners is original.physical_layout.corners
    assert physical.approaches is original.physical_layout.approaches
    assert physical.mouths==original.physical_layout.mouths
    from .shallow_bank_fillets import candidate_union_steps
    loops, shallow_arc = yield from candidate_union_steps(context, group, sources, proposals_and_fits)
    if shallow_arc:
        physical.shallow_bank_fillet = shallow_arc
    candidate=replace(group,footprints=tuple(sources),physical_layout=physical,loops=loops)
    inventory={source.region.owner:source for source in sources};arcproof=[]
    for owner,(proposal,fitted) in zip(owners,proposals_and_fits):
        spans=[rim_label(span,inventory) for loop in candidate.loops for span in loop.spans
            if span.owner==owner and rim_label(span,inventory).kind=='curb']
        indices={span.index for span in spans}
        if len(indices)!=1:raise SourceBankExposureError(owner,proposal,fitted,spans,'lost its unique physical arc')
        covered=0.;minimum=1.;maximum_gap=0.
        for span in sorted(spans,key=lambda label:(label.first,label.last)):
            minimum=min(minimum,span.first);maximum_gap=max(maximum_gap,span.first-covered)
            if span.first>covered+1e-9:raise SourceBankExposureError(owner,proposal,fitted,spans,'has an unexposed interior arc interval')
            covered=max(covered,span.last)
        if minimum>1e-9 or covered<1.-1e-9:raise SourceBankExposureError(owner,proposal,fitted,spans,'is not fully exposed')
        arcproof.append(dict(owner=owner.text(),coverage=(minimum,covered),maximum_interior_gap=maximum_gap))
    consumed={claim.section_id for claim in candidate.claims}
    before=yield from certify_layout_steps(context,candidate,consumed)
    from .original_bank_handoffs import sample as sample_original_banks
    candidate=sample_original_banks(planner,candidate,newgrids)
    after=yield from certify_layout_steps(context,candidate,consumed)
    candidate.physical_layout.source_cycle_full_width_proof=after
    if CompoundPortPlanner.incomplete_ports(candidate):raise ValueError('A new bank fillet clips an original incoming profile')
    if sum(loop.hole for loop in candidate.loops)!=sum(loop.hole for loop in original.loops):
        raise ValueError('Bank infill changed the original physical hole topology')
    if candidate.members!=original.members or candidate.claims!=original.claims:
        raise ValueError('Bank infill changed original logical ownership')
    # Use this planner's effective prepared neighbor inventory. A new fillet
    # cannot authorize a moved mouth or steal an opposite owner's interval.
    for key,layout in planner.layouts.items():
        if key in candidate.members:continue
        for (sid,end),mouth in candidate.physical_layout.mouths.items():
            opposite=getattr(layout,'mouths',{}).get((sid,1-end))
            if opposite is not None and (mouth.station>opposite.station+1e-7 if end==0
                                        else mouth.station<opposite.station-1e-7):
                raise SourceBankReconciliationError('Source bank fillet conflicts with a prepared original neighbor mouth')
    arms=candidate.physical_layout.approaches;outward=candidate.physical_layout.outward_stations
    ports=[]
    for arm in arms:
        matches=[source for source in candidate.footprints if len(getattr(source.layout,'approaches',()))==1
            and source.layout.approaches[0] is arm]
        if len(matches)!=1:raise ValueError('Bank infill lost its exact original port inventory')
        ports.append(matches[0])
    unchanged,remote=advance_past_overlap(planner,arms,outward,ports,
        tuple(source for source in candidate.footprints if not getattr(source.region,'boundary_only',False)),
        newgrids,candidate.members,max_candidates=0)
    if tuple(unchanged)!=tuple(outward) or remote:raise ValueError('Bank infill requires an unauthorized mouth advance')
    return candidate,newgrids,dict(all_original_area_vertices_retained=oldvertices,
        original_corners_identical=physical.corners is original.physical_layout.corners,
        original_approaches_identical=physical.approaches is original.physical_layout.approaches,
        original_mouths_exact=physical.mouths==original.physical_layout.mouths,
        original_members_claims_exact=True,arcs=arcproof,raw_source_width=before,canonical_source_width=after,
        complete_future_strip_guard='PASS_WITHOUT_MOUTH_ADVANCE',
        neighbor_mouth_inventory='CURRENT_PLANNER_EFFECTIVE_ORIGINAL_OWNERS',
        prepared_owner_area_guard='PENDING_FINAL_PLAN_INVENTORY')


def _record(proposal,fitted):
    banks=[]
    for bank,arm,contact,crossing in zip((proposal.first,proposal.last),
            (fitted.corner.first,fitted.corner.second),fitted.corner.stations,fitted.intersection):
        banks.append(dict(section=bank.path.section.id,side=bank.side,owned_intervals=bank.intervals,
            native_contact=contact if arm.end==0 else arm.length-contact,
            native_crossing=crossing if arm.end==0 else arm.length-crossing,
            original_origins=[dict(site=origin.site.text(),section=origin.endpoint.section_id,
                end=origin.endpoint.end,outward_side=origin.outward_side) for origin in bank.origins]))
    record = dict(kind=proposal.kind,original_site=proposal.site.text(),banks=banks,
        original_sectors=[dict(section=end.section_id,end=end.end,radius=radius,mode=mode)
                          for end,radius,mode in proposal.sectors],
        preferred_radius=proposal.radius,radius_mode=proposal.mode,radius=fitted.corner.radius,
        centre=tuple(map(float,fitted.corner.center)),start_angle=fitted.corner.start_angle,sweep=fitted.corner.sweep)
    if getattr(fitted.corner,'fit_policy',None) == 'FINITE_NATIVE_LOCAL_OFFSET':
        record['contact_fit'] = fitted.corner.fit_policy
    return record


class SourceBankReconciler:
    """Prepare local geometry while keeping every original ownership contract."""
    def __init__(self,planner,group,partitions,*,maximum_fillets=64):
        if type(maximum_fillets) is not int or not 1<=maximum_fillets<=64:
            raise ValueError('Source-bank reconciliation needs a bounded fillet allowance')
        self.planner,self.original,self.partitions=planner,group,partitions
        self.context=planner.context;self.maximum_fillets=maximum_fillets

    def build(self):
        return drain_steps(self.build_steps())

    def build_steps(self):
        original=self.original
        if (original.physical_layout is None
                or getattr(original.physical_layout,'source_bank_reconciliation',None)):
            return original
        if not SourceBankDiscovery(self.context,original).discover():return original
        from .original_bank_handoffs import sample as sample_original_banks
        current=sample_original_banks(self.planner,original,self.partitions)
        if len(getattr(current.physical_layout,'source_bank_handoff_canonicalizations',()))>len(
                getattr(original.physical_layout,'source_bank_handoff_canonicalizations',())):
            # A certified exact original exit changes only sampled ownership.
            # Keep that corrected source base if discovery finds no new fillet.
            original=current
        fitted=[];proofs=[];grids=self.partitions
        for _ in range(self.maximum_fillets):
            selected=SourceBankDiscovery(self.context,current).discover()
            if not selected:
                if not fitted:return original
                current.physical_layout.source_bank_reconciliation=dict(policy=SOURCE_BANK_POLICY,
                    records=[_record(proposal,fit) for proposal,fit in fitted],proofs=proofs,
                    original_corners_unchanged=True,original_profiles_and_mouths_unchanged=True,
                    prepared_owner_area_guard='PENDING_FINAL_PLAN_INVENTORY')
                return current
            if len(fitted)+len(selected)>self.maximum_fillets:
                raise SourceBankReconciliationError('Source bank repair exceeds its finite transition budget')
            for proposal in selected:
                try:
                    fit = fit_exposed_switch(self.context,current,proposal)
                except CornerContactError as error:
                    raise UnfittedBankSwitch(current,proposal,error) from error
                fitted.append((proposal,fit))
            current,grids,proof=yield from _apply_steps(self.planner,original,grids,fitted)
            proofs.append(proof)
        raise SourceBankReconciliationError('Source bank repair did not remove its typed boundary switches')
