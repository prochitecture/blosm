"""Failed-only ownership repair at an explicitly retained original arc exit.

A numerical union can assign a nearly coincident native bank to the wrong
original operand. This does not authorize another fillet, a moved contact, or
accepting a geometric gap. Retain the original complete circle and recover its
exact source-bank ledger only when every finite ownership certificate agrees.
"""
from copy import copy
from dataclasses import replace
import numpy as np
from .model import End,Key
from .profiles import Crown
from .compounds import rim_label
from .regions import BoundaryLoop,SourceSpan,Region,RegionUnion
from .source_curb_coordinates import CurbCoordinates
from .source_bank_reconciliation import SourceBankDiscovery,CircularCornerWitness
from .sampling import AnalyticalBoundarySampler,AnalyticalContactError


class OriginalBankHandoffUnavailable(ValueError):
    """The failed handoff has no complete exact original-source certificate."""


class _ContactSampler(AnalyticalBoundarySampler):
    def _contact(self,first,last,labels,layouts):
        try:return super()._contact(first,last,labels,layouts)
        except AnalyticalContactError:
            self.failed_contact=first,last,labels
            raise


def _candidate(context,group,first,last,labels,*,maximum_spans=8):
    sources={source.region.owner:source for source in group.footprints}
    discovery=SourceBankDiscovery(context,group)
    witness=discovery.describe(first,True)
    if (not isinstance(witness,CircularCornerWitness) or first.owner.kind!='cycle_corner'
            or getattr(sources[first.owner].layout,'generated_junction_area',None) is not None
            or labels[0].last!=1. or witness.entry_parameter!=0. or witness.exit_parameter!=1.):
        raise OriginalBankHandoffUnavailable('Failed handoff needs one complete original circle and its native banks')
    source=sources[first.owner];corner=witness.corner
    curve=source.layout.boundaries[labels[0].index]
    site=discovery.sites[witness.site]
    ends=(End(corner.first.source.id,corner.first.end),End(corner.second.source.id,corner.second.end))
    if (ends[0]==ends[1] or not set(ends)<=set(site.ends)
            or source.layout.original_junction!=witness.site
            or not any(value is corner for value in group.physical_layout.corners)):
        raise OriginalBankHandoffUnavailable('Original circle changed its actual original incidence inventory')
    profiles=[]
    for arm in (corner.first,corner.second):
        native=context[Key('profile',arm.source.id)]
        if type(native) is not Crown or type(arm.profile) is not Crown or arm.profile!=native.oriented(arm.end):
            raise OriginalBankHandoffUnavailable('Original circle changed its exact oriented source profile')
        profiles.append((arm.source.id,arm.end,native.rise,native.datum,native.half,arm.width,arm.offset))
    wanted_radius,wanted_mode=site.corner_radius(ends[0])
    if (corner.radius_mode!=wanted_mode or corner.requested_radius!=wanted_radius
            or wanted_mode=='STRICT' and corner.radius!=wanted_radius):
        raise OriginalBankHandoffUnavailable('Original circle changed its explicit radius contract')
    declared=getattr(curve,'original_circle_interval',None)
    if declared is not None:
        declared.validate(curve,corner)
        if declared.prepared_bounds!=(0.,1.) or declared.auxiliary_leads:
            raise OriginalBankHandoffUnavailable('Original exit repair cannot replace an auxiliary or partial circle')
    entry,exit=curve(0.),curve(1.)
    expected=corner.second.point(corner.stations[1],corner.second.offset+.5*corner.second.width)
    if not np.array_equal(exit,expected):
        raise OriginalBankHandoffUnavailable('Original exit is not its exact native bank contact')
    matching=[]
    for li,loop in enumerate(group.loops):
        for at,span in enumerate(loop.spans):
            previous=loop.spans[at-1]
            if (previous.owner==first.owner and span.owner==last.owner
                    and rim_label(previous,sources)==labels[0] and rim_label(span,sources)==labels[1]):
                matching.append((li,at))
    if len(matching)!=1:raise OriginalBankHandoffUnavailable('Failed original exit is not a unique physical handoff')
    li,start=matching[0];loop=group.loops[li]
    if not np.array_equal(loop.xy[start],exit):
        raise OriginalBankHandoffUnavailable('Failed boundary node is not the exact original exit')
    at=start-1;segments=[]
    while at>=0:
        label=rim_label(loop.spans[at],sources)
        if (loop.spans[at].owner,label.kind,label.index)!=(first.owner,'curb',labels[0].index):break
        segments.append(label);at-=1
    segments.reverse()
    if (not segments or segments[0].first!=0. or segments[-1].last!=1.
            or any(a.last!=b.first for a,b in zip(segments,segments[1:]))
            or not np.array_equal(loop.xy[at+1],entry)):
        raise OriginalBankHandoffUnavailable('Failed handoff does not expose the whole original arc')
    entry_bank=discovery.describe(loop.spans[at],True)
    if entry_bank is None or not discovery.same_bank(entry_bank,witness.first):
        raise OriginalBankHandoffUnavailable('Original entry does not retain its declared first bank')
    previous_label=rim_label(loop.spans[at],sources)
    previous_curve=sources[loop.spans[at].owner].layout.boundaries[previous_label.index]
    if not np.array_equal(previous_curve(previous_label.last),entry):
        raise OriginalBankHandoffUnavailable('Original entry has no exact physical bank handoff')
    target=None
    for index in range(start,min(len(loop.spans),start+maximum_spans)):
        span=loop.spans[index];operand=sources[span.owner]
        bank=discovery.describe(span,False)
        if index>start and bank is not None and discovery.same_bank(bank,witness.last):
            if (span.last==1. and np.array_equal(operand.region.xy[span.segment],exit)
                    and np.array_equal(operand.region.xy[(span.segment+1)%len(operand.region.xy)],
                                       loop.xy[(index+1)%len(loop.xy)])):
                target=index;break
        if bank is None or not discovery.same_bank(bank,witness.first):break
    if target is None:
        raise OriginalBankHandoffUnavailable('No exact finite declared exit-bank interval follows the wrong-owner run')
    old=loop.spans[target];replacement=SourceSpan(old.owner,old.segment,0.,old.last)
    operand=sources[old.owner];base=operand.region.xy[old.segment]
    endpoint=operand.region.xy[(old.segment+1)%len(operand.region.xy)]
    target_label=operand.segments[old.segment]
    target_curve=operand.layout.boundaries[target_label.index]
    if (type(getattr(target_curve,'source_coordinates',None)) is not CurbCoordinates
            or not np.array_equal(target_curve(target_label.first),base)
            or not np.array_equal(target_curve(target_label.last),endpoint)):
        raise OriginalBankHandoffUnavailable('Declared exit bank lost its exact finite native chord ledger')
    vector=endpoint-base;length2=float(vector@vector)
    reference=[loop.xy[start:target+2]]
    for span in loop.spans[start:target]:
        label=rim_label(span,sources);wrong=sources[span.owner].layout.boundaries[label.index]
        if type(getattr(wrong,'source_coordinates',None)) is not CurbCoordinates:
            raise OriginalBankHandoffUnavailable('Wrong-owner run has no explicit native bank ledger')
        parameters=[label.first,label.last]
        parameters.extend(float(s/wrong.length) for s in wrong.stations if label.first<s/wrong.length<label.last)
        reference.append(wrong(np.asarray(parameters)))
    reference=np.vstack(reference);fraction=((reference-base)@vector)/length2
    if np.any(fraction<0.) or np.any(fraction>1.):
        raise OriginalBankHandoffUnavailable('Wrong-owner run leaves the finite declared source chord')
    maximum=float(np.linalg.norm(reference-base-fraction[:,None]*vector,axis=1).max())
    if maximum>RegionUnion().epsilon:
        raise OriginalBankHandoffUnavailable('Wrong-owner run exceeds the unchanged numerical welding band')
    minimum=np.minimum(reference.min(axis=0),base);maximum_xy=np.maximum(reference.max(axis=0),endpoint)
    allowed={first.owner,replacement.owner,*(span.owner for span in loop.spans[start:target])}
    for other in group.footprints:
        if other.region.owner in allowed:continue
        if (np.all(other.region.xy.max(axis=0)>=minimum)
                and np.all(other.region.xy.min(axis=0)<=maximum_xy)):
            raise OriginalBankHandoffUnavailable('An unrelated third source can cover the candidate ownership run')
    xy=np.concatenate((loop.xy[:start],loop.xy[start:start+1],loop.xy[target+1:]))
    spans=(*loop.spans[:start],replacement,*loop.spans[target+1:])
    loops=list(group.loops);loops[li]=BoundaryLoop(xy,spans)
    RegionUnion().build((Region(group.members[0],xy[::-1] if loop.hole else xy),))
    physical=copy(group.physical_layout)
    candidate=replace(group,loops=tuple(loops),physical_layout=physical)
    construction=dict(kind='EXACT_ORIGINAL_CIRCLE_EXIT_BANK',site=witness.site.text(),
        ends=tuple((end.section_id,end.end) for end in ends),profiles=tuple(profiles),
        radius=corner.radius,requested_radius=corner.requested_radius,mode=corner.radius_mode,
        centre=tuple(map(float,corner.center)),start_angle=corner.start_angle,sweep=corner.sweep,
        contact_stations=tuple(corner.stations),first_bank=witness.first.token,last_bank=witness.last.token,
        first_owned_intervals=witness.first.intervals,last_owned_intervals=witness.last.intervals,
        original_owner=first.owner.text(),original_curve=labels[0].index,
        exit_bank_owner=old.owner.text(),exit_bank_curve=target_label.index,
        exit_interval=(target_label.first,target_label.last),
        exit_vertex=tuple(map(float,base)),other_end_vertex=tuple(map(float,endpoint)))
    return candidate,dict(construction=construction,whole_original_arc=(0.,1.),
        source_exit_exact=True,other_end_exact=True,third_source_candidates=0,
        whole_run_reference_distance=maximum,removed_wrong_owner_spans=target-start,
        original_source_operands_unchanged=True,original_circle_and_radius_unchanged=True)


def _validate(planner,original,result,partitions):
    from .cycle_dimensional_retry import certify_layout
    from .compound_ports import CompoundPortPlanner
    from .source_mouth_clearance import advance_past_overlap
    proof=certify_layout(planner.context,result,{claim.section_id for claim in result.claims})
    if CompoundPortPlanner.incomplete_ports(result):
        raise OriginalBankHandoffUnavailable('Original exit repair clips an incoming profile')
    if sum(loop.hole for loop in result.loops)!=sum(loop.hole for loop in original.loops):
        raise OriginalBankHandoffUnavailable('Original exit repair changed physical hole topology')
    old,new=original.physical_layout,result.physical_layout
    if (result.footprints is not original.footprints or result.claims is not original.claims
            or result.members is not original.members or new.corners is not old.corners
            or new.approaches is not old.approaches or new.mouths is not old.mouths):
        raise OriginalBankHandoffUnavailable('Original exit repair changed original source operands or profiles')
    arms,outward=new.approaches,new.outward_stations;ports=[]
    for arm in arms:
        matches=[source for source in result.footprints if len(getattr(source.layout,'approaches',()))==1
                 and source.layout.approaches[0] is arm]
        if len(matches)!=1:raise OriginalBankHandoffUnavailable('Original exit repair changed a complete incoming port')
        ports.append(matches[0])
    unchanged,remote=advance_past_overlap(planner,arms,outward,tuple(ports),
        tuple(source for source in result.footprints if not getattr(source.region,'boundary_only',False)),
        partitions,result.members,max_candidates=0)
    if tuple(unchanged)!=tuple(outward) or remote:
        raise OriginalBankHandoffUnavailable('Original exit repair needs an unauthorized mouth advance')
    new.source_cycle_full_width_proof=proof
    return dict(fresh_width_proof=proof,original_profiles_complete=True,
                future_strip_guard='PASS_WITHOUT_MOUTH_ADVANCE',original_holes_unchanged=True)


def sample(planner,group,partitions):
    """Keep normal successful sampling; repair only certified failed handoffs."""
    settings=planner.context.network.settings;current=group;records=[];first_failure=None
    for _ in range(17):
        sampler=_ContactSampler(tolerance=settings.curve_tolerance,spacing=settings.spacing,max_samples=20000)
        try:result=sampler.build(current)
        except AnalyticalContactError as failure:
            if first_failure is None:first_failure=failure
            if len(records)==16:raise first_failure
            try:
                current,record=_candidate(planner.context,current,*sampler.failed_contact)
            except (ValueError,KeyError,AttributeError) as unavailable:
                raise first_failure from unavailable
            records.append(record)
            continue
        if not records:return result
        validation=_validate(planner,group,result,partitions)
        previous=tuple(getattr(group.physical_layout,'source_bank_handoff_canonicalizations',()))
        result.physical_layout.source_bank_handoff_canonicalizations=(*previous,*(dict(record,**validation) for record in records))
        return result
    raise first_failure
