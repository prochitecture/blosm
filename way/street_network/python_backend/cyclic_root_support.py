"""Finite C2 support of an exposed outer root in one authored cyclic junction.

The caller must first exhaust the failed facing corner's circular STYLE
alternatives. This helper replaces only an adjacent full-profile curb; every
other supplied circle and all original carriageway geometry stay authoritative.
"""
import math
import numpy as np
from .compound_root_envelopes import _search_steps
from .corner_patch import corner_patch,_stations
from .corner_source_support import _canonical_root_stations
from .intersections import IntersectionProvider,ApproachExtentError
from .model import End,Key
from .regions import BoundaryIncidenceError,Region,RegionUnion,signed_area
from .mesh import cross


def _prepare_steps(owner,donor,settings,*,site,original_approaches,support_corners,required_roots,
            source_regions=(),consumed=(),maximum_stations=None,crest=None):
    """Return an independently certified cyclic-root curve, not a full mesh.

    ``required_roots`` contains actual (original Approach, exposed outer side)
    pairs. Full finite local strips are proof operands only; the final source
    union must still retain every complete profile, source width and hole.
    """
    originals=tuple(original_approaches);support=tuple(support_corners);required=tuple(required_roots)
    source_regions=tuple(source_regions);consumed=frozenset(consumed)
    by_end={End(a.source.id,a.end):a for a in originals}
    if (not 3<=len(originals)<=8 or len(by_end)!=len(originals)
            or set(by_end)!=set(site.ends) or not 1<=len(required)<=len(originals)-2
            or not 1<=len(support)<=16):
        raise ValueError('Cyclic root support needs one complete declared original incidence set')
    origin=originals[0].origin
    if any(np.linalg.norm(a.origin-origin)>1e-8 or abs(a.offset)>1e-8
            or abs(a.width-a.source.width)>1e-8 or abs(a.length-a.path.length)>1e-8
            or getattr(a.path,'absorbed_sites',()) or getattr(a.path,'mouth_piece',None) is not None
            or getattr(a.path,'deferred_profile_reserve',None) is not None for a in originals):
        raise ValueError('Cyclic root support requires unchanged centred original approaches')
    ordered=sorted(originals,key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
    cyclic={(End(a.source.id,a.end),End(ordered[(i+1)%len(ordered)].source.id,ordered[(i+1)%len(ordered)].end))
            for i,a in enumerate(ordered)}
    def pair(c):return tuple(End(a.source.id,a.end) for a in (c.first,c.second))
    def original(c):return all(by_end.get(End(a.source.id,a.end)) is a for a in (c.first,c.second))
    own=pair(donor);preferred=donor.requested_radius if donor.requested_radius is not None else donor.radius
    requested,mode=site.corner_radius(own[0])
    if crest is None:
        valid_profile=(not donor.first.profile.half and not donor.second.profile.half
                       and abs(float(donor.first.profile.value(-1.)))<=1e-10
                       and abs(float(donor.second.profile.value(1.)))<=1e-10)
    else:
        valid_profile=crest.donor is donor
    if (not original(donor) or own not in cyclic or donor.center is None
            or mode!='STYLE' or donor.radius_mode!='STYLE' or abs(preferred-requested)>1e-8
            or not math.pi*.5<donor.gap<math.pi*1.5
            or not valid_profile):
        raise ValueError('Cyclic root support needs an original through-half-plane STYLE outer-curb donor')
    bounds=None
    if maximum_stations is not None:
        bounds=tuple(maximum_stations)
        if (len(bounds)!=2 or any(not math.isfinite(s) or not 0.<=s<=a.length
                for a,s in zip((donor.first,donor.second),bounds))):
            raise ValueError('Root support contact bounds need two finite original donor stations')
    roots=[]
    for a,side in required:
        endpoint=End(a.source.id,a.end)
        if (by_end.get(endpoint) is not a or endpoint in own or endpoint in roots
                or side not in (-1,1) or abs(float(a.profile.value(side)))>1e-10):
            raise ValueError('Cyclic root support requires actual exposed zero-height outer halves')
        roots.append(endpoint)
    seen={own};support_regions=[]
    for i,c in enumerate(support):
        token=pair(c)
        if not original(c) or c.center is None or token not in cyclic or token in seen:
            raise ValueError('Cyclic root support must retain distinct actual local circular corners')
        seen.add(token)
        support_regions.append(corner_patch(Key('cyclic_root_support_corner',i),c,
            settings.spacing,settings.curve_tolerance).region)

    def audit(candidate,source,settings,*,include_source_regions=False):
        regions=[source.region,*support_regions];caps=[]
        if include_source_regions:
            regions.extend(source_regions if crest is None else crest.regions(candidate,source_regions))
        contacts={end:[] for end in by_end}
        for c in (candidate,*support):
            for a,s in zip((c.first,c.second),c.stations):
                contacts[End(a.source.id,a.end)].append(max(0.,s))
        for i,a in enumerate(originals):
            fixed=contacts[End(a.source.id,a.end)]
            extent=min(a.length,max(2.*a.width,max(fixed,default=0.)+.25))
            at=_stations(a,extent,min(.5,settings.spacing),settings.curve_tolerance,8192)
            for s in fixed:at=_canonical_root_stations(at,extent,s)
            xy=np.vstack((a.point(at,.5*a.width),a.point(at,-.5*a.width)[::-1]))
            caps.append(xy[[-1,0]])
            if signed_area(xy)<0.:xy=xy[::-1]
            regions.append(Region(Key('cyclic_root_support_strip',i),xy))
        exposed=0.
        union=RegionUnion(max_segments=20000) if include_source_regions else RegionUnion()
        try:
            loops=yield from union.build_steps(regions)
        except BoundaryIncidenceError:
            # The candidate touches another operand at one point. Reject this
            # trial and let the bounded search test its remaining envelopes.
            return False
        for loop in loops:
            # A root cannot be accepted just because another owner supplied
            # the same cap edge. Reject every positive boundary overlap with
            # every actual original transverse root, regardless of its tag.
            for cap in caps:
                axis=cap[1]-cap[0];length=float(np.linalg.norm(axis))
                delta=loop.xy-cap[0]
                close=abs(cross(delta,axis))/length<=1e-8
                t=delta@axis/(axis@axis)
                low=np.maximum(0.,np.minimum(t,np.roll(t,-1)))
                high=np.minimum(1.,np.maximum(t,np.roll(t,-1)))
                if np.any(close&np.roll(close,-1)&((high-low)*length>1e-7)):return False
            for span in loop.spans:
                if loop.hole or span.owner!=source.region.owner:continue
                label=source.segments[span.segment]
                if label.kind=='curb':exposed+=abs((label.last-label.first)*(span.last-span.first))
        # Keeping the complete new curve exposed also prevents an accidental
        # non-tangent crossing of a third source curb after local clipping.
        return abs(exposed-1.)<=1e-8

    stations=[IntersectionProvider.mouth_station(a,max(0.,s)+min(.25,.05*a.width),settings)
              for a,s in zip((donor.first,donor.second),donor.stations)]
    from .continuations import ContinuationReachError
    support_width=max(a.width for a in originals)
    try:
        result=yield from _search_steps(owner,donor,settings,stations,root_audit_steps=audit,
            admission='EXPOSED_CYCLIC_ORIGINAL_ZERO_HEIGHT_ROOT',support_width=support_width)
    except ContinuationReachError as failure:
        previous=failure;result=None
    else:
        if source_regions and not (yield from audit(result.corner,result.source,settings,include_source_regions=True)):
            previous=ContinuationReachError('Other original source operands clip the complete cyclic support curve')
            result=None
    if result is None:
        if not source_regions:raise previous
        width=min(donor.first.width,donor.second.width)
        required_width=max(a.width for a,side in required)
        limit=max(.25*width,.5*max(0.,required_width-width))
        if limit<=.25*width:raise previous
        root_points=tuple(a.point(0.,side*.5*a.width) for a,side in required)+(
            donor.first.point(0.,-.5*donor.first.width),donor.second.point(0.,.5*donor.second.width))
        held=tuple(i for i,a in enumerate((donor.first,donor.second))
            if getattr(a.path,'pieces',()) and all(p.section.id in consumed for p in a.path.pieces))
        families=[]
        if len(held)==1:
            index=held[0]
            # These are bounded search candidates, not authored cuts or a
            # clearance proof. Every whole curb is tested against ALL other
            # prepared source operands before accepting one of these reaches.
            values=dict.fromkeys((.75*width,.5*width,.25*width,.125*width,stations[index]))
            families.extend((index,float(value)) for value in values)
        families.append((None,None))
        from .corner_radii import CornerRadiusPolicy
        trials=[];radii=CornerRadiusPolicy.contact_candidates(donor.first,donor.second,preferred,'STYLE')
        for radius in radii:
            for anchor,held_station in families:
                initial=list(stations)
                if anchor is not None:initial[anchor]=held_station
                try:
                    result=yield from _search_steps(owner,donor,settings,initial,anchor=anchor,
                        root_audit_steps=lambda c,s,st:audit(c,s,st,include_source_regions=True),
                        admission='EXPOSED_CYCLIC_ORIGINAL_ZERO_HEIGHT_ROOT',support_width=support_width,
                        minimum_radius=radius,root_points=root_points,outward_limit=limit)
                except (ContinuationReachError,ApproachExtentError) as failure:
                    trials.append(dict(minimum_radius=radius,held_end=anchor,
                        held_station=held_station,failure=str(failure)))
                    continue
                trials.append(dict(minimum_radius=radius,held_end=anchor,
                    held_station=held_station,qualified=True))
                result.report['cyclic_projection_recovery']=dict(
                    previous_failure=str(previous),actual_required_support_width=required_width,
                    donor_minimum_width=width,maximum_outward_bound=limit,
                    other_original_source_operands=len(source_regions),global_whole_curve_exposure=True,
                    candidate_station_policy='DONOR_WIDTH_FRACTIONS_CHECKED_AGAINST_ALL_ORIGINAL_SOURCE_OPERANDS',
                    bounded_families_per_radius=len(families),reach_steps_per_family=17,
                    amplitude_candidates_per_reach_at_most=8,trials=trials)
                break
            if result is not None:break
        if result is None:raise ContinuationReachError('Projected cyclic root support exhausts its bounded source families') from previous
    if bounds is not None:
        baseline=tuple(result.stations);trials=[]
        if any(s>limit for s,limit in zip(baseline,bounds)):
            previous=ContinuationReachError('Root support exceeds an existing original mouth lead')
            accepted=None
            # The caller supplies actual owned lead limits. Hold either
            # original donor endpoint with the established anchored search;
            # every result must satisfy BOTH limits and all physical audits.
            for anchor in (0,1):
                initial=list(stations);initial[anchor]=bounds[anchor]
                try:
                    candidate=yield from _search_steps(owner,donor,settings,initial,anchor=anchor,
                        root_audit_steps=lambda c,s,st:audit(c,s,st,include_source_regions=bool(source_regions)),
                        admission='EXPOSED_CYCLIC_ORIGINAL_ZERO_HEIGHT_ROOT',support_width=support_width)
                except (ContinuationReachError,ApproachExtentError) as failure:
                    previous=failure;trials.append(dict(anchor=anchor,failure=str(failure)))
                    continue
                if any(s>limit for s,limit in zip(candidate.stations,bounds)):
                    trials.append(dict(anchor=anchor,stations=list(candidate.stations),
                                       failure='exceeds_original_mouth_contact_bound'))
                    continue
                trials.append(dict(anchor=anchor,stations=list(candidate.stations),qualified=True))
                accepted=candidate;break
            if accepted is None:
                raise ContinuationReachError('Root support exhausts its original mouth lead bounds') from previous
            result=accepted
        result.report['original_mouth_contact_bounds']=dict(maximum_stations=list(bounds),
            baseline_stations=list(baseline),anchored_searches=len(trials),maximum_anchored_searches=2,
            trials=trials,both_original_lead_bounds_preserved=True)
    result.report.update(original_junction=site.key.text(),original_degree=len(originals),
        required_original_roots=[dict(section=a.source.id,end=a.end,outer_side=side,width=a.width,
                                     profile_half=a.profile.half) for a,side in required],
        complete_original_incidence=[dict(section=a.source.id,end=a.end,width=a.width,
                                         profile_half=a.profile.half) for a in originals],
        unchanged_circular_support=[dict(first=dict(section=pair(c)[0].section_id,end=pair(c)[0].end),
                                        second=dict(section=pair(c)[1].section_id,end=pair(c)[1].end),
                                        radius=c.radius,stations=c.stations) for c in support],
        complete_new_curb_exposed=True,complete_new_curb_on_exterior=True,
        all_original_root_caps_geometrically_internal=True)
    return result


def prepare(owner,donor,settings,**options):
    from .updates import drain_steps
    return drain_steps(prepare_steps(owner,donor,settings,**options))


def prepare_steps(owner,donor,settings,*,site,original_approaches,support_corners,required_roots,
            source_regions=(),consumed=(),maximum_stations=None,context=None,source_group=None,
            original_bounds=None):
    """Keep the established family first; admit typed native jets on failure."""
    from .continuations import ContinuationReachError
    try:
        return (yield from _prepare_steps(owner,donor,settings,site=site,original_approaches=original_approaches,
            support_corners=support_corners,required_roots=required_roots,source_regions=source_regions,
            consumed=consumed,maximum_stations=maximum_stations))
    except (ContinuationReachError,ApproachExtentError) as failure:
        if context is None or source_group is None:raise
        from .source_root_templates import prepare as native_template,OriginalRootTemplateUnavailable
        from .source_root_ports import OriginalRootPortUnavailable
        from .source_root_witnesses import StagedRecipientUnavailable
        try:
            return native_template(context,source_group,owner,donor,settings,site=site,
                original_approaches=original_approaches,required_roots=required_roots,
                maximum_stations=maximum_stations,original_bounds=original_bounds)
        except (OriginalRootTemplateUnavailable,OriginalRootPortUnavailable,StagedRecipientUnavailable) as unavailable:
            raise failure from unavailable
