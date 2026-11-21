"""Select a width transition only after an original root is exposed."""
import math
import numpy as np
from .compounds import rim_label
from .model import Key,End


def eligible(corner):
    a,b=corner.first,corner.second
    if (abs(a.width-b.width)<=1e-6 or abs(a.offset)>1e-8 or abs(b.offset)>1e-8
            or np.linalg.norm(a.origin-b.origin)>1e-8):return False
    if ((a.profile.half or b.profile.half) and (corner.radius_mode!='STYLE'
            or abs(corner.gap-math.pi)>math.radians(15.)
            or abs(float(a.profile.value(-1.)))>1e-10
            or abs(float(b.profile.value(1.)))>1e-10)):return False
    if corner.center is None and corner.continuation:return True
    return corner.radius_mode=='STYLE' and abs(corner.gap-math.pi)<=math.radians(15.)


def exposed_sites(planner,loops,sources):
    inventory={source.region.owner:source for source in sources}
    by_end={end:key for key,site in planner.sites.items() for end in site.ends}
    exposed=set()
    for loop in loops:
        for span in loop.spans:
            if rim_label(span,inventory).kind!='internal':continue
            source=inventory[span.owner];layout=source.layout
            original=getattr(layout,'original_junction',None)
            if original in planner.sites:exposed.add(original)
            endpoint=getattr(layout,'original_root_caps',{}).get(getattr(span,'segment',None))
            if endpoint in by_end:exposed.add(by_end[endpoint])
    return exposed


def required(planner,loops,sources,candidates,retried):
    """Return one declared width pair at a failed original cross-section."""
    exposed=exposed_sites(planner,loops,sources)
    for pair,(site,corner) in sorted(candidates.items()):
        if (site in exposed and pair not in retried and eligible(corner)
                and (not (corner.first.profile.half or corner.second.profile.half)
                     or _exposed_facing_cap(loops,sources,corner))):return pair
    return None


def circular_required(planner,loops,sources,candidates,retried,settings):
    """Force a circular facing audit only at an exposed original root.

    An acute unequal-width fork can erode a finite outer root too. Its
    existing STYLE circles must exhaust the same proof before any separately
    certified outward envelope is considered.
    """
    from .corner_source_support import complete_source_roots
    exposed=exposed_sites(planner,loops,sources)
    tolerance=min(.001,.05*settings.curve_tolerance)
    for pair,(site,corner) in sorted(candidates.items()):
        a,b=corner.first,corner.second
        if (site not in exposed or pair in retried or corner.center is None or corner.radius_mode!='STYLE'):continue
        if abs(a.width-b.width)>1e-6:
            if not 0.<corner.gap<.5*math.pi or not _exposed_facing_cap(loops,sources,corner):continue
            # First establish the site's actual dimensional curb. An acute
            # circle may need support from that adjacent certified curve;
            # forcing it away beforehand would remove a valid sector.
            if not any(getattr(source.layout,'original_junction',None)==site
                    and any(hasattr(c,'dimensional_corner') or
                        getattr(getattr(c,'root_envelope',None),'report',{}).get('admission')
                        =='TYPED_FINITE_NARROW_WIDTH_LEAD_EXHAUSTION' for c in getattr(source.layout,'corners',()))
                    for source in sources):continue
        support=planner.layouts[site].approaches
        if not complete_source_roots(corner,settings.spacing,tolerance,support_approaches=support,
                audit_all_contacts=True,audit_facing_halves=True):return pair
    return None


def _exposed_facing_cap(loops,sources,corner):
    """An acute retry must address this actual exposed original root half."""
    inventory={s.region.owner:s for s in sources}
    for loop in loops:
        for span in loop.spans:
            if rim_label(span,inventory).kind!='internal' or not hasattr(span,'segment'):continue
            source=inventory[span.owner]
            endpoint=getattr(source.layout,'original_root_caps',{}).get(span.segment)
            if endpoint is None:continue
            points=source.region.xy;a=points[span.segment];delta=points[(span.segment+1)%len(points)]-a
            exposed=np.array((a+span.first*delta,a+span.last*delta))
            for index,approach in enumerate((corner.first,corner.second)):
                if endpoint!=End(approach.source.id,approach.end):continue
                low=approach.point(0.,approach.offset-.5*approach.width)
                high=approach.point(0.,approach.offset+.5*approach.width)
                edge=high-low;length=float(np.linalg.norm(edge))
                u=(exposed-low)@edge/(edge@edge)
                if (np.max(np.linalg.norm(exposed-(low+u[:,None]*edge),axis=1))>1e-7
                        or u.min()<-1e-8 or u.max()>1.+1e-8):continue
                first,last=(0.,.5) if index==0 else (.5,1.)
                if (min(last,float(u.max()))-max(first,float(u.min())))*length>1e-7:return True
    return False


def certify_layout(context,group,consumed):
    from .updates import drain_steps
    return drain_steps(certify_layout_steps(context, group, consumed))


def certify_layout_steps(context,group,consumed):
    """Prove complete widths of consumed roads and every real finite prefix."""
    from .source_strip_containment import certify_intervals_steps
    intervals={sid:[(0.,context[Key('path',sid)].length)] for sid in consumed}
    for a,s in zip(group.physical_layout.approaches,group.physical_layout.outward_stations):
        if hasattr(a,'original_approaches'):continue
        for (sid,end),station in a.path.trims(s).items():
            length=context[Key('path',sid)].length
            first,last=(0.,station) if end==0 else (station,length)
            if last-first<=1e-9:continue
            intervals.setdefault(sid,[]).append((first,last))
    result=yield from certify_intervals_steps(context,group,intervals)
    from .cycle_terminal_ports import certify as certify_terminals
    result['terminal_prefixes']=certify_terminals(group)
    return result
