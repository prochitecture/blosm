"""Continuous complete-width ownership proofs for original source strips."""
import math
import numpy as np
from .finite_containment import segments_on_polygon_side
from .mesh import cross
from .updates import drain_steps
from .model import Key
from .regions import _classify
from .prepared_source_classifier import SourceProofClassifier

SOURCE_RANGE_POLICY='3-certified-owned-prefix-curvature'


def _strip_spacing(k, reach, epsilon, spacing):
    if reach * k >= 1.:
        raise ValueError('Complete source-strip proof needs a positive offset curvature certificate')
    if k <= 1e-12:
        return spacing
    curvature = k / (1. - reach*k)
    speed = 1. + reach*k
    return min(spacing, math.sqrt(8.*epsilon/curvature)/speed, math.pi/(16.*curvature*speed))


class SourceStripContainmentError(ValueError):
    """A certified finite source cell is clipped by its physical owner."""
    def __init__(self,message,*,section_id=None):
        self.section_id=section_id
        super().__init__(message)


def certify_intervals(context,group,section_intervals,*,maximum_error=.002,maximum_cells=8192):
    return drain_steps(certify_intervals_steps(context, group, section_intervals,
        maximum_error=maximum_error, maximum_cells=maximum_cells))


def certify_intervals_steps(context,group,section_intervals,*,maximum_error=.002,maximum_cells=8192):
    """Certify exact owned ranges without filling an unowned gap between them.

    Each call retains the original Path and the existing continuous proof.
    Multiple ranges share one finite cell allowance; a returning Section's
    two endpoint prefixes cannot reset that allowance or claim its middle.
    """
    if (not math.isfinite(maximum_error) or not 0.<maximum_error<=.002
            or type(maximum_cells) is not int or not 1<=maximum_cells<=8192):
        raise ValueError('Invalid complete source-strip containment budget')
    normalized={};count=0
    for sid,ranges in sorted(section_intervals.items()):
        path=context[Key('path',sid)]
        ordered=[]
        for first,last in ranges:
            if not all(math.isfinite(v) for v in (first,last)) or not 0.<=first<last<=path.length:
                raise ValueError('Source-strip interval must remain inside its original finite Section')
            ordered.append((first,last));count+=1
            if count>maximum_cells:
                raise ValueError('Complete source-strip proof exceeds its finite interval budget')
        merged=[]
        for first,last in sorted(ordered):
            if merged and first<=merged[-1][1]:
                merged[-1]=(merged[-1][0],max(last,merged[-1][1]))
            else:
                merged.append((first,last))
        if not merged:
            raise ValueError('Source-strip ownership needs a positive original interval')
        normalized[sid]=tuple(merged)
    # Preserve the established call and report exactly for the ordinary case.
    if all(len(ranges)==1 for ranges in normalized.values()):
        intervals={sid:ranges[0] for sid,ranges in normalized.items()}
        return (yield from certify_steps(context,group,intervals,maximum_error=maximum_error,
                       maximum_cells=maximum_cells,intervals=intervals))
    reports=[];used=0;result=None
    for index in range(max(len(ranges) for ranges in normalized.values())):
        intervals={sid:ranges[index] for sid,ranges in normalized.items() if index<len(ranges)}
        if used>=maximum_cells:
            raise ValueError('Complete source-strip proof exceeds its finite cell budget')
        proof=yield from certify_steps(context,group,intervals,maximum_error=maximum_error,
                      maximum_cells=maximum_cells-used,intervals=intervals)
        reports.extend(proof['sections']);used+=proof['total_cells'];result=proof
    return dict(result,sections=reports,total_cells=used,
        interval_policy='exact_disjoint_original_source_ranges')


def _distance(points,segments):
    first=segments[:,0];edge=segments[:,1]-first
    squared=np.sum(edge*edge,axis=1)
    result=[]
    # Long retained strips can have thousands of sides. Bound all
    # point-by-side arrays while keeping each point's complete minimum.
    batch=max(1,min(128,32768//max(1,len(segments))))
    for start in range(0,len(points),batch):
        offset=points[start:start+batch,None]-first
        t=np.clip(np.sum(offset*edge,axis=-1)/squared,0.,1.)
        result.extend(np.linalg.norm(offset-t[...,None]*edge,axis=-1).min(axis=1))
    return np.asarray(result)


def certify(context,group,section_ids,*,maximum_error=.002,maximum_cells=8192,intervals=None):
    return drain_steps(certify_steps(context, group, section_ids, maximum_error=maximum_error,
        maximum_cells=maximum_cells, intervals=intervals))


def certify_steps(context,group,section_ids,*,maximum_error=.002,maximum_cells=8192,intervals=None):
    """Contain every original cross-section, including its entire half-profile.

    Certified source-curb chords form finite convex cells. Every cell edge is
    checked continuously against the outer boundary and every known hole;
    contained holes are checked separately from edge crossings. At most 2 mm
    accounts for source-curb chord error plus a physical boundary band.
    """
    if (not math.isfinite(maximum_error) or not 0.<maximum_error<=.002
            or type(maximum_cells) is not int or maximum_cells<1):
        raise ValueError('Invalid complete source-strip containment budget')
    outer=[loop for loop in group.loops if not loop.hole]
    holes=[loop.xy for loop in group.loops if loop.hole]
    if len(outer)!=1:raise ValueError('Complete source-strip ownership needs one physical exterior')
    exterior=outer[0].xy
    # This immutable polygon inventory belongs only to this proof invocation.
    classifier=SourceProofClassifier()
    epsilon=min(.00025,.125*maximum_error,.05*context.network.settings.curve_tolerance)
    rows=0;reports=[]
    identities=sorted(set(section_ids));intervals={} if intervals is None else dict(intervals)
    if set(intervals)-set(identities):raise ValueError('Source-strip intervals need their declared original Section IDs')
    for sid in identities:
        path=context[Key('path',sid)];section=path.section
        first,last=intervals.get(sid,(0.,path.length))
        if not all(math.isfinite(v) for v in (first,last)) or not 0.<=first<last<=path.length:
            raise ValueError('Source-strip interval must remain inside its original finite Section')
        reach=abs(section.offset)+.5*section.width
        if path.exact:
            # Exact centreline and miter-normal interpolation are linear on
            # each original segment, so these offset edges need no error band.
            stations=np.r_[first,np.asarray(path.knots)[(np.asarray(path.knots)>first)&(np.asarray(path.knots)<last)],last]
            chord_error=0.
        else:
            nominal = min(.5, context.network.settings.spacing)
            available = maximum_cells - rows
            spacing = _strip_spacing(path.curvature_bound, reach, epsilon, nominal)
            # A tight bend beyond the owned mouth can exhaust the whole budget.
            # Bound the same two offset curves on the actual interval in that case.
            local = 2 * math.ceil((last-first)/(2*spacing)) > available
            if local:
                bound = path.visible_curvature_bound(first, last)
                spacing = _strip_spacing(bound, reach, epsilon, nominal)
            if 2 * math.ceil((last-first)/(2*spacing)) > available:
                raise ValueError('Complete source-strip proof exceeds its finite sample budget')
            stations=path.samples(first,last,spacing,tolerance=epsilon,local_curvature=local,max_samples=available+1)
            chord_error=epsilon
        count=len(stations)-1;rows+=count
        if count<1 or rows>maximum_cells:
            raise ValueError('Complete source-strip proof exceeds its finite cell budget')
        high=path.point(stations,section.offset+.5*section.width)
        low=path.point(stations,section.offset-.5*section.width)
        cells=np.stack((high[:-1],high[1:],low[1:],low[:-1]),axis=1)
        edge=np.roll(cells,-1,axis=1)-cells
        if np.min(cross(edge,np.roll(edge,-1,axis=1)))<=1e-12:
            raise ValueError('Complete source-strip proof needs unfolded finite source cells')
        tolerance=maximum_error-chord_error
        for index in range(0, len(cells), 16):
            yield
            batch = cells[index:index+16]
            first = batch.reshape(-1, 2)
            last = np.roll(batch, -1, axis=1).reshape(-1, 2)
            # Preserve edge order when different parts fail on the same batch.
            failed_at, reason = len(first), None
            for loop, inside in ((exterior, True), *((hole, False) for hole in holes)):
                passed = segments_on_polygon_side(first[:failed_at], last[:failed_at], loop,
                    inside=inside, tolerance=tolerance, classifier=classifier)
                rejected = np.flatnonzero(~passed)
                if len(rejected):
                    failed_at = int(rejected[0])
                    reason = 'exterior clips' if inside else 'hole intersects'
                    if not failed_at:
                        break
            if reason is not None:
                raise SourceStripContainmentError(f'Compound {reason} the complete source strip of Section {sid}',section_id=sid)
        if len(exterior):
            # Only the two longitudinal source curbs bound the swept strip's
            # interior. Artificial sampling cross-sections must not hide a
            # small hole lying well inside the original carriageway.
            sides=np.concatenate((np.stack((high[:-1],high[1:]),axis=1),
                                  np.stack((low[:-1],low[1:]),axis=1)))
            caps=np.array(((high[0],low[0]),(high[-1],low[-1])))
            for boundary,is_hole in ((exterior,False),*((hole,True) for hole in holes)):
                occupied=np.zeros(len(boundary),bool)
                for index,cell in enumerate(cells):
                    if index % 128 == 0:
                        yield
                    selected=np.flatnonzero(np.all(boundary>=cell.min(axis=0)-1e-9,axis=1)
                        &np.all(boundary<=cell.max(axis=0)+1e-9,axis=1))
                    if not len(selected):continue
                    inside,hit=_classify(boundary[selected],cell,1e-9)
                    occupied[selected]|=inside | (hit>=0)
                interior=boundary[occupied]
                if not len(interior):
                    continue
                deep=_distance(interior,sides)>maximum_error
                terminal=_distance(interior,caps)<=maximum_error
                if np.any(deep & ~terminal):
                    kind='hole is contained' if is_hole else 'exterior enters'
                    raise SourceStripContainmentError(f'Compound {kind} inside the complete source strip of Section {sid}',section_id=sid)
        reports.append(dict(section=sid,interval=list(intervals.get(sid,(0.,path.length))),cells=count,source_chord_error=chord_error,
            boundary_band=tolerance,maximum_combined_error=maximum_error,
            proof='continuous_convex_source_cells_with_explicit_hole_exclusion'))
    return dict(sections=reports,total_cells=rows,maximum_error=maximum_error,
        outer_boundaries=1,holes=len(holes))
