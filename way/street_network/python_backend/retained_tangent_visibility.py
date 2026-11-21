"""Recover an explicit analytical tangent hidden only by prepared chords.

Proof refinement never becomes an area operand or a source ownership claim.
The ordinary sampler and all complete-source guards still qualify the result.
"""
import math
from dataclasses import replace
import numpy as np
from .compounds import RimSpan,RimSegment,rim_label
from .finite_containment import segment_on_polygon_side,separated_segment
from .model import Key
from .regions import BoundaryLoop


def _fail(message):
    from .retained_arc_fillets import RetainedArcReconciliationError
    raise RetainedArcReconciliationError(message)


def _spend(budget,count):
    if budget['used']+count>budget['maximum']:
        _fail('Hidden retained tangent exceeds its proof-refinement budget')
    budget['used']+=count


def _native_strip_clear(context,group,source,points,error,coarse_error,budget):
    from .paths import Path
    from .source_curb_coordinates import CurbCoordinates
    sid=source.region.owner.id
    if context is None or source.region.owner.kind!='cycle_strip' or sid not in {c.section_id for c in group.claims}:
        _fail('Nearby hidden-tangent source needs its complete original native strip')
    path=context[Key('path',sid)]
    coordinates=[getattr(curve,'source_coordinates',None) for curve in source.layout.boundaries]
    if (type(path) is not Path or not coordinates
            or any(type(c) is not CurbCoordinates or c.path is not path or c.approach is not None
                   or c.maximum_chord_error is None or c.maximum_chord_error>coarse_error for c in coordinates)):
        _fail('Hidden tangent has no certified complete native source grid')
    for curve,c in zip(source.layout.boundaries,coordinates):c.validate(curve.points,context)
    grid=np.unique(np.concatenate([c.stations for c in coordinates]))
    if grid[0]!=0. or grid[-1]!=path.length or len(grid)>20001:
        _fail('Hidden tangent exclusion must retain the whole finite original source')
    bound=max(c.maximum_chord_error for c in coordinates)
    right=path.point(grid,path.section.offset+.5*path.section.width)
    left=path.point(grid,path.section.offset-.5*path.section.width)
    selected=[]
    for i,(first,last) in enumerate(zip(grid,grid[1:])):
        cell=np.array((right[i],right[i+1],left[i+1],left[i]))
        if np.any(cell.max(0)+bound+error<points.min(0)) or np.any(cell.min(0)-bound-error>points.max(0)):continue
        remaining=budget['maximum']-budget['used']
        if remaining<3:_fail('Hidden retained tangent exceeds its proof-refinement budget')
        try:stations=path.samples(first,last,.1,tolerance=error,max_samples=remaining)
        except ValueError as failure:
            if str(failure)=='Path sampling exceeds its allocation budget':
                _fail('Hidden retained tangent exceeds its proof-refinement budget')
            raise
        _spend(budget,len(stations))
        polygon=np.vstack((path.point(stations,path.section.offset+.5*path.section.width),
                           path.point(stations,path.section.offset-.5*path.section.width)[::-1]))
        if not all(separated_segment(a,b,polygon,clearance=2.*error+1e-7) for a,b in zip(points,points[1:])):
            _fail('A real complete native source covers the hidden retained tangent')
        selected.append((float(first),float(last),len(stations)))
    return dict(whole_native_interval=(0.,path.length),whole_grid_cells=len(grid)-1,
                refined_cells=tuple(selected),source_error=error,arc_error=error)


def restore_hidden(group,fillet,checker,*,context=None,budget=None):
    from .retained_arc_fillets import OriginalAnalyticalCircleCurve,shared_contact
    if budget is None:budget=dict(maximum=4096,used=0)
    if (set(budget)!={'maximum','used'} or type(budget['maximum']) is not int
            or not 3<=budget['maximum']<=4096 or type(budget['used']) is not int
            or not 0<=budget['used']<=budget['maximum']):
        _fail('Invalid bounded hidden-tangent proof allocation')
    sources={source.region.owner:source for source in group.footprints}
    layouts={owner:source.layout for owner,source in sources.items()}
    if any(span.owner==fillet.owner for loop in group.loops for span in loop.spans):
        _fail('Hidden tangent insertion cannot replace an already visible arc')
    parents=fillet.parent_arcs;contacts=[]
    for role,parent in enumerate(parents):
        curve=layouts[parent.owner].boundaries[parent.curve_index]
        if (not isinstance(curve,OriginalAnalyticalCircleCurve) or curve.corner is not parent.original_corner
                or curve.original_prepared_curve is not parent.original_curve):
            _fail('Hidden tangent changed its exact original parent ancestry')
        low,high=sorted(parent.exposed_circle_parameters)
        consumed=parent.consumed_circle_parameters
        finite=(parent.analytical_contact_parameter,parent.original_chord_contact_parameter,
                *parent.original_chord_contact_xy,*parent.analytical_contact_xy,*consumed,low,high)
        if (not np.isfinite(finite).all() or not low<=consumed[0]<consumed[1]<=high
                or consumed[0 if role==0 else 1]!=parent.analytical_contact_parameter):
            _fail('Hidden tangent leaves its finite original circular interval')
        contacts.append(float(curve.source_parameter(parent.analytical_contact_parameter)))
    matches=[]
    for li,loop in enumerate(group.loops):
        labels=[rim_label(span,sources) for span in loop.spans];choices=[]
        for role,parent in enumerate(parents):
            value=parent.original_chord_contact_parameter
            choices.append([i for i,(span,label) in enumerate(zip(loop.spans,labels))
                if span.owner==parent.owner and label.kind=='curb' and label.index==parent.curve_index
                and label.first-1e-12<=value<=label.last+1e-12
                and (label.first<contacts[0] if role==0 else contacts[1]<label.last)])
        if len(choices[0])==len(choices[1])==1:matches.append((li,choices[0][0],choices[1][0]))
    if len(matches)!=1:_fail('Hidden tangent lacks its unique recorded parent-chord incidence')
    li,begin,end=matches[0];loop=group.loops[li]
    labels=[rim_label(span,sources) for span in loop.spans]
    order=np.roll(np.arange(len(loop.spans)),-begin);stop=int(np.flatnonzero(order==end)[0])
    if not 0<stop<len(order)-1:_fail('Hidden tangent reverses its finite parent chain')
    ancestry=[]
    for position,index in enumerate(order[:stop+1]):
        span=loop.spans[index];label=labels[index]
        candidates=[parent for parent in parents if parent.owner==span.owner and parent.curve_index==label.index]
        if len(candidates)!=1 or label.kind!='curb':_fail('Hidden tangent cannot consume a third source boundary')
        parent=candidates[0]
        first=max(label.first,contacts[0]) if position==0 else label.first
        last=min(label.last,contacts[1]) if position==stop else label.last
        if not any(old.index==label.index and old.first-1e-10<=first<last<=old.last+1e-10 for old in parent.source_labels):
            _fail('Hidden tangent escapes its recorded consumed parent-chord intervals')
        ancestry.append((span.owner.text(),label.index,first,last))
    virtual=RimSpan(fillet.owner,RimSegment('curb',0,0.,1.))
    # These are contacts in the INTERIOR of existing prepared chords. Their
    # recorded chord parameters, rather than that chord's unrelated endpoint,
    # anchor the unchanged ordinary neighbourhood and normal-error checks.
    chord_labels=[]
    for role,(parent,index) in enumerate(zip(parents,(begin,end))):
        label=labels[index];value=parent.original_chord_contact_parameter
        u=(value-label.first)/(label.last-label.first)
        point=loop.xy[index]+u*(loop.xy[(index+1)%len(loop.xy)]-loop.xy[index])
        if np.linalg.norm(point-np.asarray(parent.original_chord_contact_xy))>1e-7:
            _fail('Hidden tangent changed its recorded prepared chord contact')
        chord_labels.append(replace(label,**({'last':value} if role==0 else {'first':value})))
    shared_contact(checker,loop.spans[begin],virtual,(chord_labels[0],virtual.label),layouts)
    shared_contact(checker,virtual,loop.spans[end],(virtual.label,chord_labels[1]),layouts)
    error=min(1e-8,.01*fillet.source_chord_tolerance)
    count=max(2,int(math.ceil(fillet.length/math.sqrt(8.*fillet.radius*error))))
    _spend(budget,count+1);points=fillet(np.linspace(0.,1.,count+1));proof=[]
    parent_by_owner={parent.owner:parent for parent in parents}
    for source in group.area_footprints:
        if source.region.owner==fillet.owner:continue
        polygon=source.region.xy;band=2.*fillet.source_chord_tolerance+error+1e-7
        if np.any(points.max(0)+band<polygon.min(0)) or np.any(points.min(0)-band>polygon.max(0)):continue
        if source.region.owner in parent_by_owner:
            original=parent_by_owner[source.region.owner].original_corner
            if (abs(np.linalg.norm(np.asarray(fillet.center)-original.center)+fillet.radius-original.radius)>1e-7
                    or not all(segment_on_polygon_side(a,b,polygon,inside=False,
                        tolerance=fillet.source_chord_tolerance-error) for a,b in zip(points,points[1:]))):
                _fail('Hidden tangent exceeds its complete original parent chord band')
            proof.append((source.region.owner.text(),'EXACT_PARENT_TANGENCY_AND_CONTINUOUS_CHORD_BAND'))
        elif all(separated_segment(a,b,polygon,clearance=band) for a,b in zip(points,points[1:])):
            proof.append((source.region.owner.text(),'CONTINUOUS_THIRD_OPERAND_SEPARATION',band))
        else:
            proof.append((source.region.owner.text(),'COMPLETE_NATIVE_STRIP_SEPARATION',
                _native_strip_clear(context,group,source,points,error,fillet.source_chord_tolerance,budget)))
    tail=[loop.spans[i] for i in order[stop+1:]]
    spans=(RimSpan(parents[0].owner,replace(labels[begin],last=contacts[0])),virtual,
           RimSpan(parents[1].owner,replace(labels[end],first=contacts[1])),*tail)
    xy=np.vstack((loop.xy[begin],fillet(0.),fillet(1.),loop.xy[order[stop+1:]]))
    if len(xy)!=len(spans):_fail('Hidden tangent lost its complete boundary incidence')
    loops=list(group.loops);loops[li]=BoundaryLoop(xy,tuple(spans))
    return replace(group,loops=tuple(loops)),dict(fillet=fillet.owner.text(),
        proof='EXPLICIT_FINITE_PARENT_CONTACTS_AND_CONTINUOUS_SOURCE_EXCLUSION',
        consumed_parent_chords=tuple(ancestry),continuous_source_proof=tuple(proof),
        proof_points=budget['used'],maximum_proof_points=budget['maximum'],arc_error=error)
