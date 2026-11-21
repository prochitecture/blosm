"""A short, rounded grade transition between a Bundle crest and an outer curb."""
from dataclasses import replace
import math
import numpy as np

from .compound_curbs import CompositeCurb
from .continuations import ContinuationReachError
from .intersections import BoundaryCurve
from .junction_bank_transitions import JunctionBankTransitionUnavailable
from .local_width_continuation import LocalBank, LocalTransition
from .source_strip_containment import SourceStripContainmentError


def _mixed_side(corner):
    arms = corner.first, corner.second
    if (corner.center is None or corner.radius_mode != 'STYLE'
            or not 0. < abs(corner.sweep) <= math.radians(15.)
            or any(s >= 0. for s in corner.stations)
            or sum(bool(a.profile.half) for a in arms) != 1):
        return None
    half = int(bool(arms[1].profile.half))
    side = (-1, 1)[half]
    profile = arms[half].profile
    if (profile.rise <= 0. or profile.peak != side or profile.derivative(side) != 0.
            or arms[1-half].profile.datum != profile.datum):
        return None
    return half


def fit(context, site, layout):
    from .compound_surfaces import CompoundSurfaceBuilder
    from .ordinary_mouth_clearance import _canonical_candidate, _group
    original = getattr(layout, 'ordinary_complete_mouth_clearance', None)
    if not original:
        return layout
    choices = [(i, h) for i, corner in enumerate(layout.corners) if (h := _mixed_side(corner)) is not None]
    if len(choices) != 1:
        return layout
    builder = CompoundSurfaceBuilder(context.patterns, source_heights=True)
    try:
        builder.boundary_conditions(_group(site, layout), context)
    except JunctionBankTransitionUnavailable:
        pass
    except ValueError:
        # An unmodified layout can still be repartitioned by the mesh builder.
        # Only a known bank-height refusal calls for this preparation step.
        return layout
    else:
        return layout
    i, half = choices[0]
    corner = layout.corners[i]
    arms = corner.first, corner.second
    mouths = (layout.outward_stations[i], layout.outward_stations[(i+1) % len(layout.approaches)])
    tolerance = context.network.settings.curve_tolerance
    full = 1-half
    for factor in (1., 2., 3., 5.):
        contacts = list(mouths)
        contacts[full] = mouths[half] * factor
        if not 0. < contacts[full] < mouths[full]:
            break
        contacts = tuple(contacts)
        trial = replace(corner, center=None, continuation=True, stations=contacts,
                        fit_policy='STYLE_MIXED_BANK_GRADE')
        for fraction in (0., .0025, .005, .01, .02):
            offset = fraction * min(a.width for a in arms)
            try:
                core = BoundaryCurve(trial, contacts, tolerance, outward_offset=offset)
            except ContinuationReachError:
                continue
            arm = arms[full]
            if full:
                lead = LocalBank(arm, contacts[1], mouths[1], arm.offset+.5*arm.width)
                pieces = (LocalTransition(core), lead)
            else:
                lead = LocalBank(arm, mouths[0], contacts[0], arm.offset-.5*arm.width)
                pieces = (lead, LocalTransition(core))
            curves = list(layout.boundaries)
            curves[i] = CompositeCurb(trial, pieces)
            corners = list(layout.corners)
            corners[i] = trial
            try:
                updated, proof = _canonical_candidate(context, site, replace(layout,
                    boundaries=tuple(curves), corners=tuple(corners)))
                _, _, limit = builder.boundary_conditions(_group(site, updated), context)
            except (SourceStripContainmentError, JunctionBankTransitionUnavailable):
                continue
            # TransitionHeight uses quintic smoothstep along the curb and
            # interpolates the inward crossfall between the endpoint values.
            # The derivative control hull bounds speed on the entire curve.
            axis = core.last-core.first
            axis /= np.linalg.norm(axis)
            speed_floor = float(min(core.continuation.first_derivative @ axis))
            if speed_floor <= 0.:
                continue
            crossfall = max(np.linalg.norm(a.physical_gradient(s, side))
                            for a, s, side in zip(arms, contacts, (-1., 1.)))
            grade_bound = math.hypot(1.875*arms[half].profile.rise/speed_floor, crossfall)
            if grade_bound > limit:
                continue
            updated.ordinary_complete_mouth_clearance = dict(original, complete_source_width=proof,
                original_corners_contacts_radii_preserved=False,
                mixed_bank_grade=dict(corner=i, sections=tuple(a.source.id for a in arms),
                    original_radius=corner.radius, requested_radius=corner.requested_radius,
                    curvature_bound=core.continuation.curvature_bound, outward_offset=offset,
                    contacts=contacts, mouths=mouths, grade_bound=grade_bound, slope_limit=limit))
            return updated
    return layout
