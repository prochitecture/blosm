"""Admit a bounded curb support only for a real exposed same-node root."""
import math
import numpy as np
from .model import End
from .compounds import rim_label


def prepare_required(planner,loops,sources,corners,retried,consumed,members,settings,*,retry_roots=(),retry_dimensional=(),source_group=None):
    from .updates import drain_steps
    return drain_steps(prepare_required_steps(planner,loops,sources,corners,retried,consumed,members,settings,
                       retry_roots=retry_roots,retry_dimensional=retry_dimensional,source_group=source_group))


def prepare_required_steps(planner,loops,sources,corners,retried,consumed,members,settings,*,retry_roots=(),retry_dimensional=(),source_group=None):
    """Return one certified adjacent envelope after ordinary roots exhaust.

    A source tag only identifies the declared incidence. The actual exposed
    chord must also lie on that original finite root cross-section, on a
    zero-height outer half owned by one unchanged circular corner. The other
    incident dimensional curb is tried first. A larger declared cyclic node
    may use its separately certified adjacent outer-curb family only after
    those existing root and circular-contact families exhaust.
    """
    args=(planner,loops,sources,corners,retried,consumed,members,settings,
          retry_roots,retry_dimensional,source_group)
    ordinary=yield from _prepare_required(*args,crest=False)
    if ordinary is not None or getattr(planner,'ordinary_singleton_root_bounds',None) is None:
        return ordinary
    if any(c.first.profile.half==-1 and c.second.profile.half==1 for c in corners):
        return (yield from _prepare_required(*args,crest=True))
    return None


def _prepare_required(planner,loops,sources,corners,retried,consumed,members,settings,
                      retry_roots,retry_dimensional,source_group,crest):
    from .bundle_crest_roots import crest_bundle
    from .cycle_dimensional_retry import _exposed_facing_cap
    from .corner_radii import CornerRadiusPolicy,CornerContactError
    from .corner_source_support import complete_source_roots
    from .compound_contact_routes import regular_contact
    from .compound_corner_contacts import retry_steps
    from .compound_root_envelopes import prepare_adjacent
    from dataclasses import replace
    owners={e:key for key,site in planner.sites.items() for e in site.ends}
    inventory={s.region.owner:s for s in sources};targets={}
    for loop in loops:
        for span in loop.spans:
            if rim_label(span,inventory).kind!='internal' or not hasattr(span,'segment'):continue
            source=inventory[span.owner]
            endpoint=getattr(source.layout,'original_root_caps',{}).get(span.segment)
            site=owners.get(endpoint)
            if site is None or not 3<=len(planner.sites[site].ends)<=8:continue
            approach=next((a for a in planner.layouts[site].approaches if End(a.source.id,a.end)==endpoint),None)
            if approach is None:raise ValueError('Exposed root support lost its declared original approach')
            xy=source.region.xy;a=xy[span.segment];edge=xy[(span.segment+1)%len(xy)]-a
            points=a+np.array((span.first,span.last))[:,None]*edge
            low=approach.point(0.,approach.offset-.5*approach.width)
            high=approach.point(0.,approach.offset+.5*approach.width);axis=high-low
            u=(points-low)@axis/(axis@axis)
            if np.max(np.linalg.norm(points-(low+u[:,None]*axis),axis=1))>1e-7 or u.min()<-1e-8 or u.max()>1.+1e-8:
                raise ValueError('An exposed source tag does not describe its actual original root cross-section')
            eligible=((u.min()<.5 and abs(float(approach.profile.value(-1.)))<=1e-10)
                or (u.max()>.5 and abs(float(approach.profile.value(1.)))<=1e-10))
            if eligible:targets[site,endpoint]=(approach,dict(section=endpoint.section_id,end=endpoint.end,
                exposed_points=points.tolist(),original_root_parameters=u.tolist(),
                exposed_length=float(np.linalg.norm(points[1]-points[0]))))
    for (site,endpoint),(required,exposure) in sorted(targets.items()):
        local=[]
        for corner in corners:
            ends=tuple(End(a.source.id,a.end) for a in (corner.first,corner.second))
            if all(owners.get(end)==site for end in ends):local.append(corner)
        for corner in (() if crest else local):
            if not hasattr(corner,'dimensional_corner') or corner.radius_mode!='STYLE':continue
            pair=tuple(End(a.source.id,a.end) for a in (corner.first,corner.second))
            token=pair,endpoint
            if token in retried or endpoint in pair:continue
            support=tuple(c for c in local if c is not corner and c.center is not None
                and not hasattr(c,'root_envelope') and not hasattr(c,'root_support_envelope'))
            if len(support)!=2 or not any(_exposed_facing_cap(loops,sources,c) for c in support):continue
            active=next(c for c in support if _exposed_facing_cap(loops,sources,c))
            tolerance=min(.001,.05*settings.curve_tolerance)
            approaches=planner.layouts[site].approaches
            def accept(candidate):
                return (all(regular_contact(a,s) for a,s in zip((candidate.first,candidate.second),candidate.stations))
                    and complete_source_roots(candidate,settings.spacing,tolerance,support_approaches=approaches,
                        audit_all_contacts=True,audit_facing_halves=True))
            requested=active.requested_radius if active.requested_radius is not None else active.radius
            mode=active.radius_mode
            try:
                ordinary=yield from CornerRadiusPolicy.solve_steps(active.first,active.second,requested,mode,
                    math.radians(settings.continuation_angle_degrees),accept=accept)
            except CornerContactError:
                ordinary=yield from retry_steps(planner.context,active.first,active.second,requested,mode,set(members),
                    consumed,approaches,force_facing=True)
            if ordinary is not None:continue
            donor=next((s for s in sources if any(c is corner for c in getattr(s.layout,'corners',()))),None)
            if donor is None:raise ValueError('Adjacent root support lost its actual dimensional curb source')
            try:
                record=prepare_adjacent(donor.region.owner,corner,replace(settings,curve_tolerance=tolerance),
                    required_roots=(required,),support_corners=support)
            except ValueError:continue
            record.report.update(original_junction=site.text(),actual_exposed_original_root=exposure,
                ordinary_facing_source_roots_exhausted=True)
            return token,record
        maximum_contacts=None
        request=getattr(planner,'ordinary_singleton_root_bounds',None)
        if request is not None:
            from .ordinary_singleton_sources import OriginalSingletonRootBounds
            from .source_strip_containment import SourceStripContainmentError
            if (not isinstance(request,OriginalSingletonRootBounds)
                    or not isinstance(request.failure,SourceStripContainmentError)
                    or tuple(members)!=(site,) or consumed or request.site!=site
                    or request.layout is not planner.layouts[site]):
                raise ValueError('Ordinary root support lost its typed original singleton context')
            maximum_contacts=request.limits()
            if set(maximum_contacts)!=set(planner.sites[site].ends):
                raise ValueError('Ordinary root support lost its complete original End bounds')
        if len(planner.sites[site].ends)<4 and maximum_contacts is None:continue
        from .cycle_dimensional_retry import eligible as dimensional_eligible
        if any(dimensional_eligible(c) and not hasattr(c,'dimensional_corner')
                and tuple(End(a.source.id,a.end) for a in (c.first,c.second)) not in retry_dimensional
                for c in local):continue
        for active in local:
            if active.center is None or not _exposed_facing_cap(loops,sources,active):continue
            active_pair=tuple(End(a.source.id,a.end) for a in (active.first,active.second))
            if (not active.first.profile.half and not active.second.profile.half
                    and abs(active.first.width-active.second.width)<=1e-6
                    and abs(active.gap-math.pi)<=math.radians(15.) and active_pair not in retry_roots):
                continue  # Preserve the established standalone root family first.
            first_side=active_pair[0]==endpoint;side=-1 if first_side else 1
            if abs(float(required.profile.value(side)))>1e-10:continue
            candidates=[c for c in local if c is not active and c.center is not None
                and c.radius_mode=='STYLE'
                and (crest_bundle(planner.context.network,c) is not None if crest else
                     (not c.first.profile.half and not c.second.profile.half))
                and math.pi*.5<c.gap<math.pi*1.5
                and ((End(c.first.source.id,c.first.end)==active_pair[1]) if first_side else
                     (End(c.second.source.id,c.second.end)==active_pair[0]))
                and endpoint not in (End(c.first.source.id,c.first.end),End(c.second.source.id,c.second.end))]
            if not candidates:continue
            tolerance=min(.001,.05*settings.curve_tolerance);approaches=planner.layouts[site].approaches
            def accept(candidate):
                return (all(regular_contact(a,s) for a,s in zip((candidate.first,candidate.second),candidate.stations))
                    and complete_source_roots(candidate,settings.spacing,tolerance,support_approaches=approaches,
                        audit_all_contacts=True,audit_facing_halves=True))
            requested=active.requested_radius if active.requested_radius is not None else active.radius
            try:
                ordinary=yield from CornerRadiusPolicy.solve_steps(active.first,active.second,requested,active.radius_mode,
                    math.radians(settings.continuation_angle_degrees),accept=accept)
            except CornerContactError:
                ordinary=yield from retry_steps(planner.context,active.first,active.second,requested,active.radius_mode,
                    set(members),consumed,approaches,force_facing=True)
            if ordinary is not None:continue
            from .cyclic_root_support import prepare_steps as prepare_cyclic
            if crest:
                from .bundle_crest_roots import prepare_steps as prepare_cyclic
            for corner in candidates:
                pair=tuple(End(a.source.id,a.end) for a in (corner.first,corner.second));token=pair,endpoint
                if token in retried:continue
                donor=next((s for s in sources if any(c is corner for c in getattr(s.layout,'corners',()))),None)
                if donor is None:raise ValueError('Cyclic root support lost its actual original donor source')
                support=tuple(c for c in local if c is not corner and c.center is not None)
                try:
                    record=yield from prepare_cyclic(donor.region.owner,corner,replace(settings,curve_tolerance=tolerance),
                        site=planner.sites[site],original_approaches=approaches,support_corners=support,
                        required_roots=((required,side),),consumed=consumed,
                        context=planner.context,source_group=source_group,
                        original_bounds=getattr(planner,'ordinary_singleton_root_bounds',None),
                        source_regions=tuple(s.region for s in sources
                            if getattr(s.layout,'original_junction',None)!=site),
                        **({} if maximum_contacts is None else
                           dict(maximum_stations=tuple(maximum_contacts[e] for e in pair))))
                except ValueError as failure:
                    if maximum_contacts is not None:
                        from .continuations import ContinuationReachError
                        from .intersections import ApproachExtentError
                        if not isinstance(failure,(ContinuationReachError,ApproachExtentError)):raise
                    continue
                record.report.update(actual_exposed_original_root=exposure,
                    ordinary_facing_source_roots_exhausted=True,existing_dimensional_root_families_exhausted=True)
                return token,record
    return None
