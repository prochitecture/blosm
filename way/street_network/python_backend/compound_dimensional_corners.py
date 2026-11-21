"""Finite original-source curb transitions which retain the wider root.

A symmetric width taper can cut the original wider cross-section at the source
junction.  Anchor the certified width transition at that wider root instead,
and spend its finite reach on the narrower original approach.  Source strips,
complete incoming profiles and the final compound union remain authoritative.
"""
from dataclasses import dataclass
from types import SimpleNamespace
import math
import numpy as np

from .compound_curbs import CompositeCurb, OffsetCurb
from .compounds import RimSegment, SampledFootprint
from .continuations import ContinuationReachError
from .corner_patch import _PreparedCurve, _stations
from .intersections import BoundaryCurve, IntersectionProvider, ApproachExtentError
from .regions import Region, RegionUnion, signed_area
from .sampling import CertifiedCurbSampler


class _CertifiedContinuation:
    def __init__(self, curve):
        self.curve, self.length = curve, curve.length

    def __call__(self, t):
        return self.curve(t)

    def sample(self, sampler, first, last):
        count = self.curve.continuation.intervals(sampler.spacing, sampler.tolerance)
        return sampler.regular(first, last, 1./count)


class DimensionalCompositeCurb(CompositeCurb):
    """Keep the root/width proof's sampling budget during final redistribution."""
    def __init__(self, corner, pieces, tolerance):
        super().__init__(corner, pieces)
        self.prepared_tolerance = float(tolerance)

    def __call__(self, t):
        result = super().__call__(t)
        # The global-to-local arc-length conversion may round off a knot.
        # Reuse the port's exact point at each shared sample.
        for knots, points in getattr(self,'shared_bank_points',()):
            from .shared_bank_sampling import _snap_samples
            result = _snap_samples(t,knots,points,result)
        return result

    def certified_parameters(self, sampler, first, last):
        bounded = CertifiedCurbSampler(min(sampler.tolerance, self.prepared_tolerance),
                                       sampler.spacing, sampler.max_samples)
        return super().certified_parameters(bounded, first, last)


@dataclass(frozen=True)
class DimensionalCorner:
    source: SampledFootprint
    curve: CompositeCurb
    stations: tuple
    report: dict


class DimensionalLeadExhausted(ApproachExtentError):
    """The root-anchored search exhausted only its original narrower lead."""
    def __init__(self,owner,exhausted,*,narrow,wide,minimum_stations,stations,last_finite_stations):
        super().__init__(owner,exhausted)
        self.narrow,self.wide=int(narrow),int(wide)
        self.minimum_stations=tuple(minimum_stations)
        self.stations=tuple(stations)
        self.last_finite_stations=tuple(last_finite_stations)


def _patch(owner, corner, curve, stations, settings, maximum=8192):
    from .updates import drain_steps
    return drain_steps(_patch_steps(owner,corner,curve,stations,settings,maximum))


def _patch_steps(owner, corner, curve, stations, settings, maximum=8192):
    """Close the physical curb through the two real original centrelines."""
    spacing, tolerance = settings.spacing, settings.curve_tolerance
    first, second = corner.first, corner.second
    a = _stations(first, stations[0], spacing, tolerance, maximum)
    b = _stations(second, stations[1], spacing, tolerance, maximum)[::-1]
    sampler = CertifiedCurbSampler(tolerance, spacing, maximum)
    parameters = sampler.sample(curve, 0., 1.)
    curb = curve(parameters)
    curve.contact_reference_points = curb.copy()
    curve.contact_reference_points.flags.writeable = False
    before = np.vstack((first.frame(a)[0], curb[0]))
    after = np.vstack((curb[-1], second.frame(b)[0]))

    # Each construction seam keeps its exact prepared source points.  The
    # physical curb retains the actual certified polynomial/offset chain.
    curves, points, labels = [], [], []
    for chain, kind, actual, values in ((before, 'internal', None, None),
                                      (curb, 'curb', curve, parameters),
                                      (after, 'internal', None, None)):
        chain = np.asarray(chain)
        keep = np.r_[True, np.linalg.norm(np.diff(chain, axis=0), axis=1) > 1e-8]
        chain = chain[keep]
        if len(chain) < 2:
            continue
        actual = actual if actual is not None else _PreparedCurve(chain)
        values = values[keep] if values is not None else actual.stations/actual.length
        index = len(curves)
        curves.append(actual)
        points.extend(chain[:-1])
        labels.extend(RimSegment(kind, index, float(u), float(v)) for u, v in zip(values, values[1:]))
    points = np.asarray(points)
    if len(points) > maximum or len(points) != len(labels):
        raise ValueError('Dimensional source patch exceeds its finite sampling budget')
    if (np.linalg.norm(before[0]-after[-1]) > 1e-8 or signed_area(points) <= 0.):
        raise ValueError('Dimensional source patch lacks one simple oriented original junction')
    region = Region(owner, points)
    yield from RegionUnion().build_steps((region,))
    layout = SimpleNamespace(approaches=(), outward_stations=(), boundaries=tuple(curves),
                             rows=len(points), mouths={}, corners=(corner,))
    return SampledFootprint(region, tuple(labels), layout)


def prepare(owner, corner, settings, minimum_stations=(0., 0.), *, original_outer_half_profiles=False,
            original_inner_half_profiles=False, fixed_stations=False):
    """Return a certified finite curb, its source patch and real endpoints.

    The supplied corner already carries the accepted STYLE or STRICT radius.
    This operation never reduces it and never changes original source paths or
    dimensions.  Its bounded search only grows the narrower finite curb lead.
    The caller must retain both source strips and qualify the complete union.
    ``original_outer_half_profiles`` is reserved for an exposed original
    STYLE root with two zero-height facing curbs; it changes no half-profile
    values and leaves the other Bundle plateau edge authoritative.
    ``fixed_stations`` qualifies exactly one supplied original contact pair;
    it never searches farther along either source when that pair fails.
    """
    first, second = corner.first, corner.second
    approaches = first, second
    if (corner.center is not None or not corner.continuation
            or not math.isfinite(corner.radius) or corner.radius <= 0.
            or len(minimum_stations) != 2
            or any(not math.isfinite(s) or s < 0. for s in minimum_stations)
            or not all(math.isfinite(v) and v > 0. for v in (settings.spacing, settings.curve_tolerance))):
        raise ValueError('Invalid finite dimensional source-corner request')
    if np.linalg.norm(first.origin-second.origin) > 1e-8:
        raise ValueError('Dimensional curb needs one declared original source junction')
    if type(original_outer_half_profiles) is not bool or type(original_inner_half_profiles) is not bool:
        raise ValueError('Invalid original outer half-profile opt-in')
    if type(fixed_stations) is not bool:
        raise ValueError('Fixed dimensional contacts require an explicit boolean')
    if original_outer_half_profiles and original_inner_half_profiles:
        raise ValueError('A dimensional root cannot select both original inner and outer halves')
    has_half = any(a.profile.half for a in approaches)
    if original_inner_half_profiles:
        values=[float(first.profile.value(-1.)),float(second.profile.value(1.))]
        derivatives=[float(first.profile.derivative(-1.)),float(second.profile.derivative(1.))]
        if (first.profile.half!=-1 or second.profile.half!=1 or corner.radius_mode!='STYLE'
                or abs(corner.gap-math.pi)>math.radians(15.)
                or min(values)<=0. or abs(values[0]-values[1])>1e-10
                or max(map(abs,derivatives))>1e-10
                or abs(first.profile.datum-second.profile.datum)>1e-10):
            raise ValueError('Dimensional inner root needs two unchanged positive facing original half-profile peaks')
    elif has_half and (not original_outer_half_profiles or corner.radius_mode != 'STYLE'
            or abs(corner.gap-math.pi) > math.radians(15.)
            or abs(float(first.profile.value(-1.))) > 1e-10
            or abs(float(second.profile.value(1.))) > 1e-10):
        raise ValueError('Dimensional curb requires unmodified original full-profile approaches or opted zero-height outer halves')
    if any(abs(a.width-a.source.width) > 1e-8
           or (has_half and (abs(a.profile.rise-a.source.rise) > 1e-8
                            or abs(a.profile.datum-a.source.datum) > 1e-8))
           or getattr(a.path, 'absorbed_sites', ())
           or getattr(a.path, 'mouth_piece', None) is not None
           or getattr(a.path, 'deferred_profile_reserve', None) is not None for a in approaches):
        raise ValueError('Dimensional curb requires unmodified original full-profile approaches')
    # The wider original root contains the narrower one for centred roads.
    # Offset changes need a separate containment certificate, never a guessed
    # equivalent width or a fabricated uniform source route.
    if abs(first.offset) > 1e-8 or abs(second.offset) > 1e-8 or abs(first.width-second.width) <= 1e-6:
        raise ValueError('Dimensional root anchoring requires two distinct centred original widths')
    narrow = 0 if first.width < second.width else 1
    wide = 1-narrow
    stations = [IntersectionProvider.mouth_station(a, max(float(s), min(.25, .05*a.width)), settings)
                for a, s in zip(approaches, minimum_stations)]
    if fixed_stations and tuple(stations)!=tuple(minimum_stations):
        raise ContinuationReachError('Fixed dimensional contacts cannot move an original profile lead')
    reach = max(corner.radius, .5*max(first.width, second.width))
    increment = max(.125, reach/8.)
    last_error = None
    previous_stations=None
    for step in range(1 if fixed_stations else 17):
        exhausted = [(a, s) for a, s in zip(approaches, stations)
                     if not IntersectionProvider.mouth_fits(a, s)]
        if exhausted:
            if (step>0 and all(a is approaches[narrow] for a,s in exhausted)
                    and all(IntersectionProvider.mouth_fits(a,s) for a,s in zip(approaches,minimum_stations))):
                raise DimensionalLeadExhausted(owner,exhausted,narrow=narrow,wide=wide,
                    minimum_stations=minimum_stations,stations=stations,last_finite_stations=previous_stations)
            raise ApproachExtentError(owner, exhausted)
        endpoints = list(stations)
        endpoints[wide] = 0.
        try:
            core = BoundaryCurve(corner, endpoints, settings.curve_tolerance)
            pieces = []
            if wide == 0:
                pieces.append(OffsetCurb(first, stations[0], 0., -.5*first.width))
            pieces.append(_CertifiedContinuation(core))
            if wide == 1:
                pieces.append(OffsetCurb(second, 0., stations[1], .5*second.width))
            curve = DimensionalCompositeCurb(corner, pieces, settings.curve_tolerance)
            source = _patch(owner, corner, curve, tuple(stations), settings)
            report = dict(policy='FINITE_ORIGINAL_WIDTH_CONTINUATION',
                          first_section=first.source.id, first_end=first.end,
                          second_section=second.source.id, second_end=second.end,
                          first_width=first.width, second_width=second.width,
                          stations=list(stations), wider_root=wide,
                          minimum_radius=corner.radius,
                          curvature_bound=core.continuation.curvature_bound,
                          reach_steps=step, original_sources_unchanged=True)
            if fixed_stations:report['fixed_original_contact_candidate']=True
            if has_half:
                report.update(original_profile_halves=[a.profile.half for a in approaches],opposite_profiles_preserved=True)
                if original_inner_half_profiles:
                    report.update(original_inner_half_profiles=True,facing_inner_heights=values,
                                  facing_inner_derivatives=derivatives)
                else:
                    report.update(original_outer_half_profiles=True,
                                  facing_outer_heights=[float(first.profile.value(-1.)),float(second.profile.value(1.))])
            return DimensionalCorner(source, curve, tuple(stations), report)
        except ContinuationReachError as error:
            last_error = error
        previous_stations=tuple(stations)
        stations[narrow] = IntersectionProvider.mouth_station(
            approaches[narrow], stations[narrow]+increment, settings)
    raise ContinuationReachError('Dimensional source curb exhausts its bounded finite reach') from last_error
