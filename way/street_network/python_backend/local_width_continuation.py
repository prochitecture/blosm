"""Fit a width transition inside the available ordinary approach leads."""
from dataclasses import replace
import numpy as np

from .compound_curbs import CompositeCurb, OffsetCurb
from .compound_dimensional_corners import _CertifiedContinuation
from .continuations import ContinuationReachError
from .intersections import BoundaryCurve, Mouth
from .model import End, Key
from .updates import drain_steps
from .paths import right
from .source_strip_containment import SourceStripContainmentError


class LocalBank(OffsetCurb):
    """Original carriageway bank between the local transition and the mouth."""


class LocalTransition(_CertifiedContinuation):
    """The quintic segment between the actual contacts on the facing banks."""


class TransitionHeight:
    def __init__(self, core, first, last):
        self.core, self.first, self.last = core, first, last
        heights, slopes = [], []
        corner = core.corner
        for t, arm, side, station in zip((0., 1.), (corner.first, corner.second),
                                         (-1., 1.), core.mouths):
            velocity = core.continuation.derivative(t)
            inward = -right(velocity / np.linalg.norm(velocity))
            gradient = arm.physical_gradient(station, side)
            slope = float(gradient @ inward)
            if slope < -1e-9 or np.linalg.norm(gradient - inward*slope) > 1e-6:
                raise ValueError('Local width transition cannot match its original bank gradient')
            heights.append(arm.profile.datum + arm.profile.value(side))
            slopes.append(slope)
        self.heights, self.slopes = heights, slopes

    def sample(self, t):
        u = self.first + np.asarray(t)*(self.last - self.first)
        weight = u**3 * (10. - 15.*u + 6.*u*u)
        velocity = self.core.continuation.derivative(u)
        speed = np.linalg.norm(velocity, axis=-1)
        tangent = velocity / speed[..., None]
        delta = self.heights[1] - self.heights[0]
        along = delta * 30.*u**2*(1.-u)**2 / speed
        inward = self.slopes[0] + weight*(self.slopes[1] - self.slopes[0])
        return self.heights[0] + delta*weight, along[..., None]*tangent - inward[..., None]*right(tangent)


def fit(context, site, layout, failing_section, certify, *, corner_index=None):
    return drain_steps(fit_steps(context, site, layout, failing_section, certify, corner_index=corner_index))


def _continuous_separator_crest(context, first, second):
    """Two consecutive Sections of one Bundle member can share its inner crest."""
    if (first.source.street_id != second.source.street_id or first.profile.rise != second.profile.rise
            or first.profile.datum != second.profile.datum):
        return False
    ends = {End(arm.source.id, arm.end) for arm in (first, second)}
    if len(ends) != 2 or any(arm.profile != context[Key('profile', arm.source.id)].oriented(arm.end)
                             for arm in (first, second)):
        return False
    return any(ends == {End(before.section_id, int(before.forward)), End(after.section_id, int(not after.forward))}
               for bundle in context.network.bundles for member in bundle.members
               for before, after in zip(member.path, member.path[1:]))


def fit_steps(context, site, layout, failing_section, certify, *, corner_index=None, extrapolated=False):
    tolerance = context.network.settings.curve_tolerance
    attempts = 0
    for i, corner in enumerate(layout.corners):
        if corner_index is not None and i != corner_index:
            continue
        arms = corner.first, corner.second
        if (any(a.offset != 0. for a in arms)
                or abs(arms[0].width - arms[1].width) <= 1e-6
                or failing_section not in {p.section.id for a in arms for p in a.path.pieces}):
            continue
        narrow = int(arms[1].width < arms[0].width)
        if extrapolated:
            # A circle touching the backward extension of a narrow street can
            # cut off part of the wide street at its root. Join the actual banks.
            if (corner.center is None or corner.radius_mode != 'STYLE' or corner.median_cap
                    or corner.stations[narrow] >= 0. or corner.stations[1-narrow] < 0.
                    or arms[1-narrow].source.id != failing_section
                    or arms[0].profile.datum != arms[1].profile.datum
                    or (arms[0].profile.half == -1 and arms[1].profile.half == 1
                        and not _continuous_separator_crest(context, *arms))
                    or arms[0].direction @ arms[1].direction >= 0.):
                continue
        elif not corner.continuation or corner.center is not None:
            continue
        mouth = layout.outward_stations[i], layout.outward_stations[(i+1) % len(layout.approaches)]
        increment = max(corner.radius, .5*arms[1-narrow].width)/8.
        for step in range(1, 18):
            contact = [0., 0.]
            contact[narrow] = step*increment
            if contact[narrow] >= mouth[narrow]:
                break
            attempts += 1
            yield
            trial_corner = (replace(corner, center=None, start_angle=0., sweep=0., continuation=True,
                                    stations=tuple(contact), fit_policy='STYLE_FINITE_WIDTH_TRANSITION')
                            if extrapolated else corner)
            try:
                core = BoundaryCurve(trial_corner, tuple(contact), tolerance)
            except ContinuationReachError:
                continue
            pieces = []
            if mouth[0] > contact[0]:
                pieces.append(LocalBank(arms[0], mouth[0], contact[0], -.5*arms[0].width))
            pieces.append(LocalTransition(core))
            if mouth[1] > contact[1]:
                pieces.append(LocalBank(arms[1], contact[1], mouth[1], .5*arms[1].width))
            boundaries = list(layout.boundaries)
            boundaries[i] = CompositeCurb(trial_corner, pieces)
            corners = list(layout.corners)
            corners[i] = trial_corner
            try:
                candidate, proof = certify(context, site, replace(layout, corners=tuple(corners), boundaries=tuple(boundaries)))
            except SourceStripContainmentError:
                continue
            return candidate, proof, dict(policy='LOCAL_SOURCE_BANK_WIDTH_CONTINUATION', corner=i,
                sections=[a.source.id for a in arms], contacts=contact, radius=corner.radius,
                **(dict(replaced_circle_stations=list(corner.stations)) if extrapolated else {}),
                curvature_bound=core.continuation.curvature_bound, candidates_tested=attempts,
                original_mouth_planes_preserved=True)
    return None


def fit_with_wider_approach_lead(context, site, layout, failing_section, certify):
    return drain_steps(wider_approach_steps(context, site, layout, failing_section, certify))


def wider_approach_steps(context, site, layout, failing_section, certify):
    """Give a width transition room on its wider approach after fixed mouths fail."""
    settings = context.network.settings
    incident = {end for node in (*context.network.connections, *context.network.intersections)
                for end in node.ends}
    for index, corner in enumerate(layout.corners):
        if not corner.continuation or corner.radius_mode != 'STYLE':
            continue
        arms = corner.first, corner.second
        wide = int(arms[1].width > arms[0].width)
        arm, narrow = arms[wide], arms[1-wide]
        if (arm.width <= narrow.width or narrow.source.id != failing_section
                or arm.offset or narrow.offset or len(arm.path.pieces) != 1
                or End(arm.source.id, 1-arm.end) in incident):
            continue
        slot = (index+wide) % len(layout.approaches)
        initial = layout.outward_stations[slot]
        for fraction in (.25, .5, .75, 1.):
            station = initial + fraction * arm.width
            if station > arm.length - settings.spacing:
                break
            stations = list(layout.outward_stations)
            stations[slot] = station
            curves = list(layout.boundaries)
            try:
                for i in ((slot-1) % len(curves), slot):
                    curves[i] = BoundaryCurve(layout.corners[i], (stations[i], stations[(i+1) % len(curves)]),
                                              settings.curve_tolerance)
            except ContinuationReachError:
                continue
            mouths = dict(layout.mouths)
            for token, value in arm.path.trims(station).items():
                mouths[token] = Mouth(float(value), 2*layout.radial)
            trial = replace(layout, boundaries=tuple(curves), mouths=mouths, outward_stations=tuple(stations))
            result = yield from fit_steps(context, site, trial, failing_section, certify, corner_index=index)
            if result is not None:
                candidate, proof, record = result
                record.update(original_mouth_planes_preserved=False, advanced_section=arm.source.id,
                              advanced_end=arm.end, original_station=initial, advanced_station=station)
                return candidate, proof, record
    return None
