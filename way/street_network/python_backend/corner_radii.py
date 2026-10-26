"""Bounded contact fitting for preferred style radii; authored radii stay strict."""
from dataclasses import replace


class CornerContactError(ValueError):
    """A radius has no accepted contact in the supported circular branch."""


class CornerRadiusPolicy:
    ratio = .8
    maximum_reductions = 7

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
    def solve(cls, first, second, radius, mode, continuation_angle=0., *, accept=None):
        from .intersections import _corner
        attempted, last_error = [], None
        for candidate in cls.candidates(radius, mode):
            attempted.append(float(candidate))
            try:
                corner = _corner(first, second, candidate, continuation_angle)
                if accept is not None and not accept(corner):
                    raise CornerContactError('Circular contact leaves insufficient approach for the full-profile entrance')
                return replace(corner, requested_radius=radius, radius_mode=mode,
                               radius_attempts=tuple(attempted))
            except CornerContactError as exc:
                last_error = exc
        if mode == 'STRICT':
            raise last_error
        raise CornerContactError(f'Preferred style radius could not be fitted; tried '
            f'{", ".join(f"{r:g}" for r in attempted)} m: {last_error}') from last_error


def fitted_radius_report(owner, layout):
    """Persist actual physical changes without rewriting preferred input values."""
    return [dict(junction=owner.text(), first_section=c.first.source.id, first_end=c.first.end,
                 second_section=c.second.source.id, second_end=c.second.end,
                 requested_radius=c.requested_radius, actual_radius=c.radius,
                 mode=c.radius_mode, attempts=list(c.radius_attempts))
            for c in layout.corners if c.requested_radius is not None and c.radius < c.requested_radius]
