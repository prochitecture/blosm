"""Bounded proof that one finite approach cannot separate from its neighbor.

This is a compound-planning condition, not permission to publish overlapping
roads. Only an original internal route to another physical junction may use
it; its eventual single exterior still needs complete source containment.
"""
import math
import numpy as np


def _samples(approach, *, spacing=.25, tolerance=.001, maximum=4096):
    from .routes import Route
    # Ordinary constant-dimensional routes only. An absorbed seam or owned
    # prefix has its own explicit interval and profile ownership contract.
    if not isinstance(approach.path,Route) or getattr(approach.path,'absorbed_sites',()):
        return None
    if approach.length/spacing>maximum-4:
        return None
    stations=[]
    error=0.
    for index,piece in enumerate(approach.path.pieces):
        path=piece.path
        if path.exact:
            return None  # Mitered corridor support needs a separate proof.
        curvature=path.curvature_bound
        offset=abs(approach.offset)
        if offset*curvature>=1.:
            return None
        speed=1.+offset*curvature
        bound=curvature/(1.-offset*curvature)
        step=spacing
        if bound>1e-12:
            step=min(step,math.sqrt(8.*tolerance/bound)/speed,
                     math.pi/(2.*bound*speed))
        if path.length/step+len(stations)>maximum-4:
            return None
        local=path.samples(0.,path.length,step,tolerance=tolerance)
        if not piece.forward:local=path.length-local[::-1]
        stations.extend((approach.path.stations[index]+local)[:-1])
        error=max(error,bound*(step*speed)**2/8.)
    stations.append(approach.length)
    if len(stations)>maximum:return None
    points=approach.point(np.asarray(stations),approach.offset)
    return points,error


def overlap_certificate(first,second,*,clearance=0.,refine=False):
    """Certify complete first-centerline proximity to the second corridor.

    Chords approximate continuously certified Smooth paths. Distance to a
    polyline is 1-Lipschitz, so endpoint distances plus half the first chord
    length bound every point of that chord. Both curve-approximation errors
    are subtracted from the combined physical half-width plus an optional
    construction clearance. A brief positive gap smaller than a required
    fillet diameter still cannot host that fillet. This admits planning only;
    it does not authorize filling grass or publishing any provisional faces.
    No extrapolated source geometry or isolated sampled coincidence establishes
    ownership.
    """
    if not math.isfinite(clearance) or clearance<0. or type(refine) is not bool:
        raise ValueError('Approach construction clearance must be finite and nonnegative')
    if np.linalg.norm(first.origin-second.origin)>1e-5:return None
    attempts=((.25,.001),(.0625,.00025),(.025,.0001)) if refine else ((.25,.001),)
    for attempt,(spacing,tolerance) in enumerate(attempts):
        a,b=(_samples(approach,spacing=spacing,tolerance=tolerance,maximum=4096) for approach in (first,second))
        if a is None or b is None:continue
        points,error=a;curve,other_error=b
        radius=.5*(first.width+second.width)+clearance-error-other_error-1e-7
        if radius<=0.:return None
        starts,edges=curve[:-1],np.diff(curve,axis=0)
        squared=np.sum(edges*edges,axis=1)
        if np.any(squared<=1e-20):return None
        distances=[]
        for start in range(0,len(points),32):
            delta=points[start:start+32,None]-starts[None]
            t=np.clip(np.sum(delta*edges[None],axis=2)/squared[None],0.,1.)
            distances.extend(np.linalg.norm(delta-t[...,None]*edges[None],axis=2).min(axis=1))
        distances=np.asarray(distances)
        segment_bound=np.maximum(distances[:-1],distances[1:])+.5*np.linalg.norm(np.diff(points,axis=0),axis=1)
        maximum=float(np.max(segment_bound))
        if maximum<radius:
            return dict(proof='continuous_chord_distance_lipschitz_bound',maximum_distance_bound=maximum,
                permitted_distance=radius,construction_clearance=clearance,source_curve_errors=[error,other_error],
                source_samples=[len(points),len(curve)],sampling_spacing=spacing,attempts=attempt+1,
                maximum_attempts=len(attempts),maximum_samples_per_source=4096,endpoint_extrapolation=False)
    return None


def persistent_overlap(first,second,*,clearance=0.,refine=False):
    """Boolean planning certificate; optional refinement keeps finite budgets."""
    return overlap_certificate(first,second,clearance=clearance,refine=refine) is not None
