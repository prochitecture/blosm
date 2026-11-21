"""Bounded contact fitting for preferred style radii; authored radii stay strict."""
from dataclasses import replace
import math


class CornerContactError(ValueError):
    """A radius has no accepted contact in the supported circular branch."""


class CornerRadiusPolicy:
    ratio = .8
    maximum_reductions = 7
    minimum_gore_radius = .5

    @classmethod
    def contact_candidates(cls, first, second, radius, mode, *, displaced=False):
        """Ordinary fitting first; an acute fork's nose has a separate scale.

        The nose between diverging carriageways is not a vehicle's turning
        radius. A category default such as 10 m may exceed their whole gap.
        Only long, acute, common-node approaches may try the additional range.
        Explicit radii and already feasible ordinary contacts stay unchanged.
        """
        ordinary = tuple(cls.candidates(radius, mode))
        yield from ordinary
        if mode != 'STYLE' or displaced or ordinary[-1] <= cls.minimum_gore_radius:
            return
        gap = math.atan2(float(first.direction[0]*second.direction[1]
                               - first.direction[1]*second.direction[0]),
                         float(first.direction @ second.direction)) % math.tau
        if (not 1e-8 < gap < math.pi/3.
                or min(first.length, second.length) <= 4.*max(first.width, second.width)
                or sum((first.origin-second.origin)**2) > 1e-10):
            return
        previous = ordinary[-1]
        # Eight additional attempts at most, with an absolute physical floor.
        for _ in range(cls.maximum_reductions):
            candidate = max(cls.minimum_gore_radius, previous*cls.ratio)
            if candidate < previous:
                yield candidate
            if candidate == cls.minimum_gore_radius:
                return
            previous = candidate
        if previous > cls.minimum_gore_radius:
            yield cls.minimum_gore_radius

    @classmethod
    def candidates(cls, radius, mode):
        if mode not in ('STRICT', 'STYLE'):
            raise ValueError('Unknown corner radius mode')
        yield radius
        if mode == 'STRICT':
            return
        # Do not replace an infeasible turn by an arbitrarily tiny fillet.
        # At most nine solves, including the original radius and the floor.
        minimum = min(radius, max(.5, .25 * radius))
        previous = radius
        for reduction in range(1, cls.maximum_reductions + 1):
            candidate = max(minimum, radius * cls.ratio ** reduction)
            if candidate < previous:
                yield candidate
            if candidate == minimum:
                return
            previous = candidate
        if minimum < previous:
            yield minimum

    @classmethod
    def solve(cls, first, second, radius, mode, continuation_angle=0., *, accept=None,
              displaced=False, median_cap=False):
        from .updates import drain_steps
        return drain_steps(cls.solve_steps(first, second, radius, mode, continuation_angle,
                                           accept=accept, displaced=displaced, median_cap=median_cap))

    @classmethod
    def solve_steps(cls, first, second, radius, mode, continuation_angle=0., *, accept=None,
                    displaced=False, median_cap=False):
        from .intersections import _corner_steps
        import math
        import numpy as np
        sides=(None,)
        if median_cap:
            if not displaced:
                raise ValueError('Paired median cap requires explicit displaced source-incidence ordering')
            sides=(1.,)
        elif displaced and np.linalg.norm(first.origin-second.origin)>1e-5:
            # A source-incidence walk establishes the cyclic order of mouths
            # at distinct origins. Their tangent angles alone do not identify
            # the inside of a rounded end cap between nearly parallel roads.
            gap=math.atan2(float(np.linalg.det(np.stack((first.direction,second.direction)))),
                           float(first.direction@second.direction))%math.tau
            sides=(None,-1. if gap<math.pi else 1.)
        attempted, last_error = [], None
        for candidate in cls.contact_candidates(first, second, radius, mode, displaced=displaced):
            if attempted:
                yield
            attempted.append(float(candidate))
            for side in sides:
                try:
                    corner = yield from _corner_steps(first, second, candidate, continuation_angle,
                                                      side_override=side, median_cap=median_cap)
                    if accept is not None and not accept(corner):
                        raise CornerContactError('Circular contact leaves insufficient approach for the full-profile entrance')
                    return replace(corner, requested_radius=radius, radius_mode=mode,
                                   radius_attempts=tuple(attempted))
                except CornerContactError as exc:
                    last_error = exc
        if mode == 'STRICT':
            raise last_error
        # A small through-angle and a width step can have no circular curb
        # contact at any permitted style radius. Use the existing certified
        # through-curb primitive as a bounded fallback, retaining the original
        # requested radius as its curvature lower bound. Boundary construction
        # must still certify its complete polynomial and finite source leads.
        mismatch = abs(first.offset + second.offset + .5*(second.width-first.width)) > 1e-6
        gap = math.atan2(float(first.direction[0]*second.direction[1]
                               - first.direction[1]*second.direction[0]),
                         float(first.direction @ second.direction)) % math.tau
        if (not displaced and not median_cap and mismatch
                and abs(gap-math.pi) <= math.radians(15.)
                and float(np.linalg.norm(first.origin-second.origin)) <= 1e-5):
            corner = yield from _corner_steps(first, second, radius, math.radians(15.))
            if accept is None or accept(corner):
                return replace(corner, requested_radius=radius, radius_mode=mode,
                               radius_attempts=tuple(attempted), fit_policy='STYLE_WIDTH_CONTINUATION')
        raise CornerContactError(f'Preferred style radius could not be fitted; tried '
            f'{", ".join(f"{r:g}" for r in attempted)} m: {last_error}') from last_error


def fitted_radius_report(owner, layout):
    """Persist actual physical changes without rewriting preferred input values."""
    return [dict(junction=owner.text(), first_section=c.first.source.id, first_end=c.first.end,
                 second_section=c.second.source.id, second_end=c.second.end,
                 requested_radius=c.requested_radius, actual_radius=c.radius,
                 mode=c.radius_mode, attempts=list(c.radius_attempts),
                 contact_source_routes=[dict(section=a.source.id,end=a.end,
                     sections=[piece.section.id for piece in a.path.pieces],
                     owned_sites=[key.text() for key in getattr(a.path,'contact_extension_through',())])
                     for a in (c.first,c.second) if getattr(a.path,'contact_extension_through',())],
                 fit_policy=(c.fit_policy if c.fit_policy != 'ORDINARY' else 'PAIRED_MEDIAN_CAP' if c.median_cap else 'ACUTE_FORK_NOSE' if c.radius < min(c.requested_radius,
                     max(.5, .25*c.requested_radius)) else 'ORDINARY'))
            for proxy in layout.corners for c in getattr(proxy, 'component_corners', (proxy,))
            if c.requested_radius is not None and (c.radius < c.requested_radius or c.fit_policy != 'ORDINARY')]
