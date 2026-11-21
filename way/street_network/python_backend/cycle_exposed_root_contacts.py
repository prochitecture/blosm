"""Activate existing source contact walks only at an exposed original cap."""
from types import SimpleNamespace
from .cycle_dimensional_retry import _exposed_facing_cap
from .model import End


def required(planner,loops,sources,unresolved,internal,consumed,retried):
    incident={end:key for key,site in planner.sites.items() for end in site.ends}
    for first,last,failure in unresolved:
        pair=first,last
        if (pair in retried or first not in internal or last not in internal
                or first.section_id not in consumed or last.section_id not in consumed):continue
        site=incident.get(first)
        if site is None or incident.get(last)!=site:continue
        if planner.sites[site].corner_radius(first)[1]!='STYLE':continue
        approaches={End(a.source.id,a.end):a for a in planner.layouts[site].approaches}
        a,b=approaches.get(first),approaches.get(last)
        if a is None or b is None or a.profile.half or b.profile.half:continue
        if not _exposed_facing_cap(loops,sources,SimpleNamespace(first=a,second=b)):continue
        proof=dict(policy='EXPOSED_ORIGINAL_FACING_ROOT_CAP_CONTACT_RETRY',
            original_junction=site.text(),first=dict(section=first.section_id,end=first.end),
            second=dict(section=last.section_id,end=last.end),original_failure=failure,
            exact_facing_root_cap_exposed=True,source_sections_wholly_consumed=True,
            original_contact_walk_family=True,global_full_cap_and_source_width_guards_unchanged=True)
        return pair,proof
    return None
