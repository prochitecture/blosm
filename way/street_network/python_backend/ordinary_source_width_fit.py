"""Fit an existing preferred circle only after complete ordinary width fails."""
from dataclasses import replace
import math

from .corner_radii import CornerContactError,CornerRadiusPolicy
from .source_strip_containment import SourceStripContainmentError


def fit(context,site,layout,failure,certify_candidate,*,maximum_candidates=64):
    """Try one STYLE corner at a time with every original mouth kept fixed.

    Existing contact candidates and floors are authoritative. Each accepted
    candidate must retain every complete prefix on both its initial certified
    sampling and the normal analytical sampler's output. Other failures stay
    fatal; this is not a recovery from allocation or unrelated geometry errors.
    """
    if type(maximum_candidates) is not int or not 1<=maximum_candidates<=64:
        raise ValueError('Ordinary source-width fitting needs a bounded candidate budget')
    sid=getattr(failure,'section_id',None)
    if sid is None:
        return None
    from .intersections import _corner,BoundaryCurve,IntersectionProvider
    settings=context.network.settings
    incident=[i for i,c in enumerate(layout.corners) if c.center is not None
        and c.radius_mode=='STYLE' and (c.first.source.id==sid or c.second.source.id==sid)]
    incident.sort(key=lambda i:(layout.corners[i].first.source.id!=sid,i))
    tested=0;trials=[]
    for index in incident:
        original=layout.corners[index]
        if original.requested_radius is None:
            continue
        attempts=list(original.radius_attempts)
        for radius in CornerRadiusPolicy.contact_candidates(original.first,original.second,
                original.requested_radius,original.radius_mode):
            if radius>=original.radius:
                continue
            if tested>=maximum_candidates:
                raise ValueError('Ordinary source-width fitting exceeds its finite candidate budget')
            tested+=1;attempts.append(float(radius))
            trial=dict(corner=index,radius=float(radius));trials.append(trial)
            try:
                corner=_corner(original.first,original.second,radius,
                    math.radians(settings.continuation_angle_degrees))
            except CornerContactError:
                trial['result']='no_finite_circular_contact'
                continue
            if corner.center is None:
                trial['result']='not_an_original_circular_candidate'
                continue
            required=[IntersectionProvider.mouth_station(a,max(0.,s)+min(.25,a.width*.05),settings)
                for a,s in zip((corner.first,corner.second),corner.stations)]
            stations=(layout.outward_stations[index],
                      layout.outward_stations[(index+1)%len(layout.approaches)])
            if any(need>station+1e-9 or not IntersectionProvider.mouth_fits(a,station)
                    for a,need,station in zip((corner.first,corner.second),required,stations)):
                trial['result']='does_not_fit_original_complete_mouths'
                continue
            corner=replace(corner,requested_radius=original.requested_radius,
                radius_mode='STYLE',radius_attempts=tuple(attempts),
                fit_policy='STYLE_COMPLETE_SOURCE_WIDTH')
            corners=list(layout.corners);corners[index]=corner
            boundaries=list(layout.boundaries)
            boundaries[index]=BoundaryCurve(corner,stations,settings.curve_tolerance)
            candidate=replace(layout,corners=tuple(corners),boundaries=tuple(boundaries))
            try:
                prepared,proof=certify_candidate(context,site,candidate)
            except SourceStripContainmentError as rejection:
                trial.update(result='complete_source_width_rejected',section=rejection.section_id)
                continue
            trial['result']='complete_initial_and_canonical_source_width_accepted'
            record=dict(policy='STYLE_COMPLETE_SOURCE_WIDTH',corner=index,
                first=[corner.first.source.id,corner.first.end],
                second=[corner.second.source.id,corner.second.end],
                requested_radius=original.requested_radius,previous_radius=original.radius,
                actual_radius=corner.radius,radius_attempts=list(attempts),
                candidates_tested=tested,candidate_trials=trials,original_mouth_planes_preserved=True,
                other_original_corners_preserved=True)
            return prepared,proof,record
    return None
