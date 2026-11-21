"""Bounded alternative contacts on a compound's original finite road strips."""
from dataclasses import replace
import math
import numpy as np
from .corner_radii import CornerContactError,CornerRadiusPolicy
from .compound_contact_routes import extend,regular_contact


def retry(context, first, second, radius, mode, members, consumed, support,*,force_facing=False):
    from .updates import drain_steps
    return drain_steps(retry_steps(context, first, second, radius, mode, members, consumed, support,
                                   force_facing=force_facing))


def retry_steps(context, first, second, radius, mode, members, consumed, support, *, force_facing=False):
    """Try owned source continuations and a curved near-through normal side.

    The ordinary solver has already failed. Alternative roots preserve curb
    identities, analytical circles and the existing STYLE radius sequence.
    Every resulting patch must still pass the complete physical union, source
    coverage and exposed-seam guards. Explicit radius overrides are unchanged.
    """
    if mode!='STYLE':return None
    if type(force_facing) is not bool:raise ValueError('Compound contact root audit needs an explicit boolean policy')
    from .intersections import _corner_steps
    from .corner_patch import corner_patch
    from .corner_source_support import complete_source_roots
    from .model import Key
    settings=context.network.settings
    tolerance=min(.001,.05*settings.curve_tolerance)
    angle=math.radians(settings.continuation_angle_degrees)
    extended=(extend(context,first,members,consumed),extend(context,second,members,consumed))
    pairs=[(first,second)]
    if extended!=(first,second):pairs.append(extended)
    for index,(a,b) in enumerate(pairs):
        def accept(candidate):
            return (all(regular_contact(p,s) for p,s in zip((a,b),candidate.stations))
                and complete_source_roots(candidate,settings.spacing,tolerance,support_approaches=support,
                    audit_all_contacts=force_facing,audit_facing_halves=force_facing))
        if index:
            try:
                candidate=yield from CornerRadiusPolicy.solve_steps(a,b,radius,mode,angle,accept=accept)
                corner_patch(Key('contact_proof',0),candidate,settings.spacing,tolerance)
                return replace(candidate,fit_policy='SOURCE_INCIDENCE_CONTACT')
            except (CornerContactError,ValueError):pass
        gap=math.atan2(float(a.direction[0]*b.direction[1]-a.direction[1]*b.direction[0]),
                       float(a.direction@b.direction))%math.tau
        if (abs(gap-math.pi)>math.radians(15.) or
                np.linalg.norm(a.origin-b.origin)>1e-7):continue
        attempted=[]
        for value in CornerRadiusPolicy.candidates(radius,mode):
            attempted.append(float(value))
            try:
                candidate=yield from _corner_steps(a,b,value,angle,side_override=(-1. if gap<math.pi else 1.))
                if (max(abs(s) for s in candidate.stations)>2.*max(a.width,b.width)
                        or not all(regular_contact(p,s) for p,s in zip((a,b),candidate.stations))
                        or not complete_source_roots(candidate,settings.spacing,tolerance,
                            support_approaches=support,audit_all_contacts=True,audit_facing_halves=True)):continue
                corner_patch(Key('contact_proof',0),candidate,settings.spacing,tolerance)
                return replace(candidate,requested_radius=radius,radius_mode=mode,
                    radius_attempts=tuple(attempted),fit_policy='SOURCE_CURVED_CONTACT_SIDE')
            except (CornerContactError,ValueError):continue
    # Existing successful alternative roots keep their established behavior.
    # A remaining local sector owns first.LEFT and second.RIGHT; the other
    # root halves belong to other original corners. Force that sector audit
    # on the original finite pair, while the final complete union still has
    # to contain all original strips and hide every preparation cross-section.
    def facing(candidate):
        return (all(regular_contact(p,s) for p,s in zip((first,second),candidate.stations))
            and complete_source_roots(candidate,settings.spacing,tolerance,support_approaches=support,
                audit_all_contacts=True,audit_facing_halves=True))
    try:
        candidate=yield from CornerRadiusPolicy.solve_steps(first,second,radius,mode,angle,accept=facing)
        corner_patch(Key('contact_proof',0),candidate,settings.spacing,tolerance)
        return replace(candidate,fit_policy='SOURCE_CORNER_FACING_ROOTS')
    except (CornerContactError,ValueError):pass
    from .compound_displaced_contacts import retry_steps as displaced_retry
    result=yield from displaced_retry(context,first,second,radius,mode,members,consumed)
    if result is not None:return result
    from .bundle_terminal_contact_walks import retry as terminal_bank_retry
    return terminal_bank_retry(context,first,second,radius,mode,members,consumed)
