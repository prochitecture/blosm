"""Retain original curb chains between displaced compound entrances.

A connected compound can have a long exposed side. Replacing that side with
one circle between its external roads need not have a finite solution. The
fallback here follows the original incidence walk instead: local prescribed
fillets, finite original internal curbs, then the next local fillet. The final
simple-region, source-containment and surface checks still decide acceptance.
"""
from dataclasses import replace
import math
import numpy as np

from .model import End,Key


class RetainedCurbOverlapError(ValueError):
    """Two real endpoint fillets overlap on an owned original curb interval."""


class RetainedSourceRouteError(ValueError):
    """The original curb walk cannot be represented by one retained lead."""


class RetainedSourceWidthError(ValueError):
    """A proposed exterior clips a complete original owned source strip."""
    def __init__(self, sections, failure):
        self.sections=tuple(sorted(sections))
        super().__init__(f'Retained exterior does not preserve its complete original source widths: {failure}')


class RetainedCornerContactError(ValueError):
    """One declared original local corner exhausts its circular contacts."""
    def __init__(self, junction, first, last, radius, mode, failure):
        self.junction,self.first,self.last=junction,first,last
        self.radius,self.mode=radius,mode
        super().__init__(f'Original retained corner {first} to {last} at {junction}: {failure}')


class SamplingBudgetCurb:
    """Keep a newly certified source-width boundary's final chord budget."""
    def __init__(self,curve,tolerance):self.curve,self.prepared_tolerance=curve,tolerance
    def __call__(self,values):return self.curve(values)
    def __getattr__(self,name):
        curve=self.__dict__.get('curve')
        if curve is None:raise AttributeError(name)
        return getattr(curve,name)
    def certified_parameters(self,sampler,first,last):
        from .sampling import CertifiedCurbSampler
        bounded=CertifiedCurbSampler(min(sampler.tolerance,self.prepared_tolerance),sampler.spacing,sampler.max_samples)
        return bounded.sample(self.curve,first,last)


class OffsetCurb:
    def __init__(self, approach, first, last, offset):
        self.approach, self.offset = approach, float(offset)
        self.parameters = np.linspace(first, last, max(2, int(abs(last-first)/.2)+2))
        points = approach.point(self.parameters, self.offset)
        self.cumulative = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
        self.length = float(self.cumulative[-1])
        if self.length <= 1e-10:
            raise ValueError('Retained original curb has no positive finite extent')

    def __call__(self, t):
        station = np.interp(np.asarray(t)*self.length, self.cumulative, self.parameters)
        return self.approach.point(station, self.offset)

    def sample(self, sampler, first, last):
        curvature, speed, fixed = sampler.lead_metrics(self.approach, self.parameters,
            self.cumulative, self.offset, 0., self.length)
        step = sampler.spacing/max(speed, 1e-12)
        if curvature > 1e-12:
            step = min(step, math.sqrt(8.*sampler.tolerance/curvature)/speed,
                       math.pi/(2.*curvature*speed))
        # This curb inverts its finite cumulative-length lookup table. Each
        # table interval is one distinct parameter map, even when its physical
        # bank is smooth. Boundary jets must not be fitted across an undeclared
        # lookup-map derivative break.
        lookup=self.cumulative/self.length
        knots = sorted({first, last, *(t for t in fixed if first < t < last),
                        *(float(t) for t in lookup if first < t < last)})
        values = []
        for a, b in zip(knots, knots[1:]):
            values.extend(sampler.regular(a, b, step/self.length)[:-1])
            if len(values)+1>sampler.max_samples:
                raise ValueError('Retained curb sampling exceeds its sample budget')
        return np.r_[values, last]


class CircularCurb:
    def __init__(self, corner):
        self.corner = corner
        self.length = abs(corner.sweep)*corner.radius
        if corner.center is None or self.length <= 1e-10:
            raise ValueError('Retained local corner requires a positive circular fillet')

    def __call__(self, t):
        angle = self.corner.start_angle + np.asarray(t)*self.corner.sweep
        return self.corner.center + self.corner.radius*np.stack((np.cos(angle), np.sin(angle)), axis=-1)

    def sample(self, sampler, first, last):
        r = self.corner.radius
        chord = min(sampler.spacing, math.sqrt(8.*sampler.tolerance*r), .5*math.pi*r)
        return sampler.regular(first, last, chord/self.length)


class CompositeCurb:
    """An arc-length-addressed chain with a certified sampling protocol."""
    continuation = None

    def __init__(self, corner, pieces):
        self.corner, self.pieces = corner, tuple(pieces)
        self.lengths = np.r_[0., np.cumsum([p.length for p in pieces])]
        self.length = float(self.lengths[-1])
        for first, last in zip(pieces, pieces[1:]):
            if np.linalg.norm(first(1.)-last(0.)) > 1e-6:
                raise ValueError('Retained compound curb has a disconnected source join')
        self.first, self.last = self.pieces[0](0.), self.pieces[-1](1.)

    def __call__(self, t):
        values = np.asarray(t)
        stations = np.clip(values, 0., 1.)*self.length
        result = np.empty(values.shape+(2,))
        indices = np.minimum(np.searchsorted(self.lengths, stations, side='right')-1, len(self.pieces)-1)
        for index, piece in enumerate(self.pieces):
            selected = indices == index
            if np.any(selected):
                result[selected] = piece((stations[selected]-self.lengths[index])/piece.length)
        return result

    def certified_parameters(self, sampler, first, last):
        values = []
        for i, piece in enumerate(self.pieces):
            lower = max(first*self.length, self.lengths[i])
            upper = min(last*self.length, self.lengths[i+1])
            if upper <= lower+1e-12:
                continue
            local = piece.sample(sampler, (lower-self.lengths[i])/piece.length,
                                 (upper-self.lengths[i])/piece.length)
            values.extend((self.lengths[i]+local[:-1]*piece.length)/self.length)
        return np.r_[values, last]


def _local_corners(planner, chain):
    from .corner_radii import CornerRadiusPolicy,CornerContactError
    by_end = {End(a.source.id, a.end): (key, a) for key, layout in planner.layouts.items()
              for a in layout.approaches}
    result = []
    for first, last in chain:
        key, a = by_end[first]
        other, b = by_end[last]
        if key != other:
            raise ValueError('Retained local curb pair lacks one original junction')
        radius, mode = planner.sites[key].corner_radius(first)
        try:
            corner = CornerRadiusPolicy.solve(a, b, radius, mode,
                math.radians(planner.context.network.settings.continuation_angle_degrees))
        except CornerContactError as failure:
            raise RetainedCornerContactError(key,first,last,radius,mode,failure) from failure
        if corner.center is None:
            if not corner.continuation:
                if np.linalg.norm(a.point(0.,a.offset-.5*a.width)-b.point(0.,b.offset+.5*b.width))>1e-8:
                    raise ValueError('Retained zero-turn source curbs do not share their original cross-section')
            else:
                from .compound_dimensional_corners import prepare
                settings=planner.context.network.settings
                record=prepare(Key('retained_dimensional',2*a.source.id+a.end),corner,
                    replace(settings,curve_tolerance=min(.001,.05*settings.curve_tolerance)))
                corner=replace(corner,stations=record.stations,fit_policy='FINITE_ORIGINAL_WIDTH_CONTINUATION')
                corner.dimensional_corner=record
        result.append(corner)
    return tuple(result)


def _composite(corner, local, mouths):
    pieces = []
    first, last = local[0], local[-1]
    def lead(approach, a, b, offset):
        if abs(a-b) > 1e-9:
            pieces.append(OffsetCurb(approach, a, b, offset))
    lead(first.first, mouths[0], first.stations[0], first.first.offset-.5*first.first.width)
    for i, c in enumerate(local):
        if hasattr(c,'dimensional_corner'):
            pieces.extend(c.dimensional_corner.curve.pieces)
        elif c.center is not None:
            pieces.append(CircularCurb(c))
        if i+1 == len(local):
            continue
        following = local[i+1]
        approach, reverse = c.second, following.first
        if (abs(approach.width-reverse.width) > 1e-6 or
                abs(approach.length-reverse.length) > 1e-6 or
                abs(float(approach.profile.value(1.))) > 1e-8 or
                abs(float(reverse.profile.value(-1.))) > 1e-8):
            raise RetainedSourceRouteError('Retained internal curb requires one compatible original route')
        end = approach.length-following.stations[0]
        if end < c.stations[1]-1e-7:
            raise RetainedCurbOverlapError('Retained local fillets overlap along their internal source curb')
        point = approach.point(end, approach.offset+.5*approach.width)
        target = reverse.point(following.stations[0], reverse.offset-.5*reverse.width)
        if np.linalg.norm(point-target) > 1e-6:
            raise RetainedSourceRouteError('Retained curb does not meet its opposite original incidence')
        lead(approach, c.stations[1], end, approach.offset+.5*approach.width)
    lead(last.second, last.stations[1], mouths[1], last.second.offset+.5*last.second.width)
    return CompositeCurb(corner, pieces)


def retained_layout(planner, specification, approaches, chains, *, retain_all=False):
    """Retry displaced circles using their original exposed curb walk.

    A finite displaced circle can still cut a returning original source road.
    retain_all uses the source walk for every multi-junction side in that case;
    local radii, profile mouths and all later qualification remain unchanged.
    """
    from .intersections import (Corner, CornerConstructionError, BoundaryCurve,
        IntersectionProvider, JunctionLayout, Mouth, ApproachExtentError, half_intervals)
    from .corner_radii import CornerContactError, CornerRadiusPolicy
    from .sampling import CertifiedCurbSampler
    settings = planner.context.network.settings
    by_end = {End(a.source.id, a.end): a for a in approaches}
    ordered = tuple(by_end[end] for end in specification.ordered_ends)
    corners, locals_by_side = [], {}
    for i, a in enumerate(ordered):
        b = ordered[(i+1) % len(ordered)]
        first, last = End(a.source.id, a.end), End(b.source.id, b.end)
        radius, mode = specification.corner_radius(first)
        chain = chains[first, last]
        retain=retain_all and len(chain)>1
        if not retain:
            try:
                corner = CornerRadiusPolicy.solve(a, b, radius, mode,
                    math.radians(settings.continuation_angle_degrees), displaced=True,
                    median_cap=(first, last) in specification.median_pairs)
            except CornerContactError as original:
                if len(chain) <= 1:
                    raise CornerConstructionError(str(original), specification, (a, b)) from original
                retain=True
        if retain:
            local = _local_corners(planner, chain)
            # The proxy is explicitly noncircular. Consumers that need one
            # real fillet must inspect component_corners or reject this case.
            corner = Corner(a, b, (local[0].stations[0], local[-1].stations[1]),
                            None, radius, 0., 0., 0.)
            corner.component_corners = local
            locals_by_side[i] = local
        corners.append(corner)
    stations = [max(0., corners[i].stations[0], corners[i-1].stations[1])
                for i in range(len(ordered))]
    stations = [IntersectionProvider.mouth_station(a, s+min(.25, a.width*.05), settings)
                for a, s in zip(ordered, stations)]
    exhausted = [(a, s) for a, s in zip(ordered, stations) if not IntersectionProvider.mouth_fits(a, s)]
    if exhausted:
        raise ApproachExtentError(specification, exhausted)
    boundaries = []
    for i, c in enumerate(corners):
        mouths = stations[i], stations[(i+1) % len(stations)]
        boundaries.append(_composite(c, locals_by_side[i], mouths) if i in locals_by_side
                          else BoundaryCurve(c, mouths, settings.curve_tolerance))
    has_dimensional=any(hasattr(local,'dimensional_corner') for c in corners for local in getattr(c,'component_corners',(c,)))
    if has_dimensional:
        boundaries=[SamplingBudgetCurb(curve,min(.001,.05*settings.curve_tolerance)) for curve in boundaries]
    sampler = CertifiedCurbSampler(min(.001,.05*settings.curve_tolerance) if has_dimensional else settings.curve_tolerance,
        settings.spacing, 16384)
    for curve in boundaries:
        if isinstance(curve, CompositeCurb) or has_dimensional:
            curve.initial_parameters = sampler.sample(curve, 0., 1.)
            curve.initial_parameters.flags.writeable = False
    radial = half_intervals(settings)
    rows = max(2*radial, *(len(sampler.sample(curve, 0., 1.))-1 for curve in boundaries))
    rows = 8*int(math.ceil(rows/8.))
    mouths = {}
    for a, distance in zip(ordered, stations):
        for identity, station in a.path.trims(distance).items():
            mouths[identity] = Mouth(float(station), 2*radial)
    return JunctionLayout(ordered, tuple(corners), tuple(boundaries), tuple(stations), mouths,
        rows, radial, np.zeros(2), [0.]*len(ordered), mesh_strategy='boundary_chart')
