"""One bounded derived-exterior sampling retry after a typed chart failure.

Tiny generated intervals at a genuine original circular contact can force radial
cells far below the useful source chord phase. Only their already sampled
exterior is projected, within one micrometre; circles, original Paths, source
caps, complete profiles and genuine holes stay authoritative and unchanged.
"""
import math
from contextvars import ContextVar
from functools import wraps
import numpy as np
from .compounds import rim_label


_construction=ContextVar('source_exterior_sampling_construction',default=None)


def _pool_key(group,context):
    return (id(context),tuple(sorted(getattr(group,'members',()))),
            tuple(sorted(c.section_id for c in getattr(group,'claims',()))))


def chart_sampling_scope(function):
    """One temporary attempt budget across nested charts of this source pool."""
    @wraps(function)
    def scoped(self,group,context,*args,**kwargs):
        state=_construction.get();token=None
        if state is None:
            state=dict(attempted=set(),declined={});token=_construction.set(state)
        try:
            result=function(self,group,context,*args,**kwargs)
            decline=state['declined'].get(_pool_key(group,context))
            if token is not None and decline is not None and isinstance(getattr(result,'metadata',None),dict):
                result.metadata=dict(result.metadata,primary_generated_phase_sampling_declined=decline)
            return result
        finally:
            if token is not None:_construction.reset(token)
    return scoped


def chart_sampling_steps(function):
    """Keep one retry budget for this chart without exposing it between advances."""
    @wraps(function)
    def scoped(self, group, context, *args, **kwargs):
        inherited = _construction.get()
        state = inherited if inherited is not None else dict(attempted=set(), declined={})
        steps = function(self, group, context, *args, **kwargs)
        try:
            while True:
                token = _construction.set(state)
                try:
                    event = next(steps)
                except StopIteration as done:
                    result = done.value
                    break
                finally:
                    _construction.reset(token)
                yield event
        finally:
            token = _construction.set(state)
            try:
                steps.close()
            finally:
                _construction.reset(token)
        decline = state['declined'].get(_pool_key(group, context))
        if inherited is None and decline is not None and isinstance(getattr(result, 'metadata', None), dict):
            result.metadata = dict(result.metadata, primary_generated_phase_sampling_declined=decline)
        return result
    return scoped


class PrimarySamplingDeclined(ValueError):
    """The original source passed; one typed candidate construction did not."""
    def __init__(self,stage,failure,proof):
        super().__init__(f'Generated-phase candidate declined at {stage}: {failure}')
        self.report=dict(stage=stage,reason=str(failure),failure_type=type(failure).__name__,
            original_width_proof=proof,maximum_candidate_rebuilds=1)


def _facing_policy(corner):
    a,b=corner.first,corner.second
    zero_outer=(abs(float(a.profile.value(-1.)))<=1e-10
                and abs(float(b.profile.value(1.)))<=1e-10)
    shared_plateau=(a.profile.half==-1 and b.profile.half==1
        and abs(a.profile.datum-b.profile.datum)<=1e-10
        and abs(a.profile.rise-b.profile.rise)<=1e-10
        and abs(float(a.profile.derivative(-1.)))<=1e-10
        and abs(float(b.profile.derivative(1.)))<=1e-10)
    if (corner.center is None or not np.isfinite(corner.radius) or corner.radius<=0.
            or np.linalg.norm(a.origin-b.origin)>1e-8
            or abs(corner.gap-math.pi)>math.radians(15.)):
        return None
    return ('ORIGINAL_EQUAL_PLATEAU' if shared_plateau else
            'ORIGINAL_ZERO_OUTER' if zero_outer else None)


def _phase_knots(curve,label,approach,station,side):
    """The exposed interval must lie between stored root/contact knots."""
    points=getattr(curve,'points',None);stations=getattr(curve,'stations',None)
    if points is None or stations is None or not np.isfinite(curve.length) or curve.length<=0.:
        return None
    root=approach.point(0.,approach.offset+side*.5*approach.width)
    contact=approach.point(station,approach.offset+side*.5*approach.width)
    tolerance=max(1e-11,64.*np.spacing(max(1.,float(np.max(abs(points))))))
    found=[]
    for point in (root,contact):
        at=np.flatnonzero(np.linalg.norm(points-point,axis=1)<=tolerance)
        if len(at)!=1:return None
        found.append(int(at[0]))
    if found[0]==found[1]:return None
    low,high=sorted(float(stations[index]/curve.length) for index in found)
    allowance=tolerance/curve.length
    if min(label.first,label.last)<low-allowance or max(label.first,label.last)>high+allowance:
        return None
    return dict(root_knot=found[0],contact_knot=found[1],curb_side=side,
                source_contact_station=float(station),curve_parameter_interval=[low,high])


def _generated_original_phase(approach,station):
    """Do not confuse an authored control knot with generated contact phase."""
    from .paths import Path
    path=approach.path
    if (getattr(path,'absorbed_sites',()) or hasattr(path,'mouth_piece')
            or getattr(path,'deferred_profile_reserve',None) is not None
            or abs(approach.width-approach.source.width)>1e-8):return False
    if getattr(path,'outward',False):
        if not getattr(path,'pieces',()):return False
        piece=path.pieces[0];original=piece.path
        if not isinstance(original,Path):return False
        knots=original.station(original.knots)
        values=knots if piece.forward else original.length-knots
    else:
        original=path
        if not isinstance(original,Path):return False
        knots=original.station(original.knots)
        values=knots if approach.end==0 else original.length-knots
    tolerance=64.*np.spacing(max(1.,original.length))
    return station<=original.length and not np.any((values>tolerance)&(values<=station+tolerance))


def _positive_source_phase(source,label,approach,station,side):
    from .model import End
    approaches=getattr(source.layout,'approaches',())
    if len(approaches)!=1 or label.index not in (0,1):return None
    port=approaches[0]
    if (End(approach.source.id,approach.end) not in getattr(source.layout,'original_root_caps',{}).values()
            or label.index!=(1 if side==-1 else 0)
            or port.path is not approach.path or port.profile is not approach.profile
            or port.source is not approach.source or port.end!=approach.end
            or not _generated_original_phase(approach,station)):return None
    return _phase_knots(source.layout.boundaries[label.index],label,approach,station,side)


def tiny_auxiliary_edges(group,settings):
    layout=group.physical_layout
    if (layout is None or not group.claims
            or not getattr(layout,'source_cycle_full_width_proof',None)
            or getattr(layout,'source_rim_normalization',None)
            or getattr(layout,'source_exterior_sampling',None)):
        return ()
    if (not group.loops or group.loops[0].hole
            or any(not loop.hole for loop in group.loops[1:])):return ()
    # This bounds a generated contact interval, not geometric approximation.
    # Projection retains its separate, unchanged one-micrometre allowance.
    limit=.005
    inventory={source.region.owner:source for source in group.footprints}
    outer=group.loops[0];witnesses=[]
    exposed={span.owner for span in outer.spans if rim_label(span,inventory).kind=='curb'}
    corners=[(source.region.owner,source.layout.corners[0]) for source in group.footprints
        if source.region.owner in exposed and len(getattr(source.layout,'corners',()))==1
        and _facing_policy(source.layout.corners[0]) is not None]
    for index,span in enumerate(outer.spans):
        label=rim_label(span,inventory)
        if label.kind!='curb':continue
        source=inventory[span.owner]
        length=float(np.linalg.norm(outer.xy[(index+1)%len(outer.xy)]-outer.xy[index]))
        if not 0.<length<=limit:continue
        curve=source.layout.boundaries[label.index];chosen=None
        for owner,corner in corners:
            for a,station,side in zip((corner.first,corner.second),corner.stations,(-1,1)):
                if not 0.<abs(station)<=limit:continue
                if station<0.:
                    if owner!=span.owner:continue
                    policy='NEGATIVE_AUXILIARY_CONTACT'
                    proof=_phase_knots(curve,label,a,float(station),side)
                else:
                    policy='POSITIVE_ORIGINAL_CURB_CONTACT'
                    proof=_positive_source_phase(source,label,a,float(station),side)
                if proof is not None:chosen=(corner,policy,proof);break
            if chosen is not None:break
        if chosen is None:continue
        corner,policy,proof=chosen;a,b=corner.first,corner.second
        witnesses.append(dict(source_owner=span.owner.text(),boundary_segment=index,
            original_sources=[a.source.id,b.source.id],original_ends=[a.end,b.end],
            original_radius=corner.radius,original_radius_mode=corner.radius_mode,
            contact_stations=list(corner.stations),original_gap=corner.gap,
            facing_height_policy=_facing_policy(corner),generated_contact_phase=policy,phase_knots=proof,
            exposed_auxiliary_edge_length=length,maximum_auxiliary_length=limit))
    return tuple(witnesses)


def eligible(group,context):
    state=_construction.get()
    return (context is not None
        and (state is None or _pool_key(group,context) not in state['attempted'])
        and bool(tiny_auxiliary_edges(group,context.network.settings)))


def _chart_failure_types():
    from .partition import QuadChartGradingError
    from .quad_charts import QuadChartCellBudgetError,QuadChartLayerBudgetError
    from .decomposition import VisibilityCutFamilyError
    from .physical_hole_cuts import PhysicalHoleCutError
    return (QuadChartGradingError,QuadChartCellBudgetError,QuadChartLayerBudgetError,
            VisibilityCutFamilyError,PhysicalHoleCutError)


def _typed_chart_failure(failure):
    return isinstance(failure,_chart_failure_types())


def recover_sampling(parent,group,context,failure,*,decline_candidate=False):
    """Rebuild normal charts once over a source-proven, equally accurate rim."""
    if not _typed_chart_failure(failure) or not eligible(group,context):raise failure
    from .cycle_dimensional_retry import certify_layout
    from .source_rim_normalization import normalize_exterior
    settings=context.network.settings
    witnesses=tiny_auxiliary_edges(group,settings)
    consumed={claim.section_id for claim in group.claims}
    state=_construction.get()
    if state is not None:state['attempted'].add(_pool_key(group,context))
    # This original-input proof is deliberately OUTSIDE candidate handling.
    # A source failure here is fatal and cannot fall through to old charts.
    before=certify_layout(context,group,consumed)
    allowance=min(1e-6,.1*settings.curve_tolerance)
    chord_allowance=min(.001,.05*settings.curve_tolerance)
    if chord_allowance+allowance>settings.curve_tolerance:
        raise ValueError('Derived exterior sampling exhausts the authored curve budget')
    normalized=normalize_exterior(group,maximum_projection=allowance)
    from .source_strip_containment import SourceStripContainmentError
    stage='candidate width proof'
    try:
        after=certify_layout(context,normalized,consumed)
        stage='candidate chart'
        result=parent.build(normalized,context)
        stage='final candidate width proof'
        final=certify_layout(context,result.partition.group,consumed)
    except (*_chart_failure_types(),SourceStripContainmentError) as rejection:
        if not decline_candidate:raise
        raise PrimarySamplingDeclined(stage,rejection,before) from rejection
    maximum=max(row['maximum_projection_error']
                for row in normalized.physical_layout.source_rim_normalization)
    result.metadata=dict(result.metadata,derived_exterior_sampling_recovery=dict(
        policy='EXPOSED_GENERATED_ORIGINAL_CIRCULAR_CONTACT_PHASE',
        previous_chart_failure=str(failure),witnesses=witnesses,
        original_source_chord_allowance=chord_allowance,added_polyline_projection=maximum,
        maximum_added_projection=allowance,composed_curve_error_bound=chord_allowance+maximum,
        authored_curve_budget=settings.curve_tolerance,original_analytical_curves_unchanged=True,
        original_profiles_exact=True,original_holes_exact=True,original_area_operands_unchanged=True,
        original_width_proof=before,normalized_width_proof=after,final_width_proof=final,
        normalization=normalized.physical_layout.source_exterior_sampling,
        maximum_rebuilds=1,cell_budget=parent.max_cells))
    return result


def try_primary_sampling(parent,group,context,failure,stage):
    """Only a recognized failed primary chart can try the generated phase."""
    if not _typed_chart_failure(failure) or not eligible(group,context):return None
    try:
        result=recover_sampling(parent,group,context,failure,decline_candidate=True)
    except PrimarySamplingDeclined as rejection:
        state=_construction.get()
        if state is not None:state['declined'][_pool_key(group,context)]=dict(rejection.report,primary_stage=stage)
        return None
    result.metadata=dict(result.metadata,primary_generated_phase_sampling_recovery=dict(
        primary_stage=stage,primary_failure_type=type(failure).__name__,primary_failure=str(failure),
        original_primary_succeeded=False,maximum_candidate_rebuilds=1))
    return result
