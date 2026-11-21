"""Bounded STYLE alternatives when a finite circle erodes original root caps.

This helper is not a circular contact solver. Its caller must first establish
an exposed original root after the prescribed circular alternatives fail.
Original finite approach jets and the requested minimum radius remain fixed;
a certified zero-jet bell moves only the connecting curb outward.
"""
from dataclasses import dataclass, replace
import math
import numpy as np

from .compound_dimensional_corners import DimensionalCompositeCurb, _CertifiedContinuation, _patch_steps
from .continuations import WidthContinuation, ContinuationReachError, _evaluate
from .curb_envelope import CurbEnvelope
from .corner_patch import _stations
from .corner_source_support import _canonical_root_stations
from .free_pads import EnvelopeBoundary
from .intersections import IntersectionProvider, ApproachExtentError
from .model import Key
from .paths import right
from .regions import Region, RegionUnion, signed_area


@dataclass(frozen=True)
class RootEnvelope:
    source: object
    curve: object
    stations: tuple
    report: dict
    corner: object


def _complete_roots(corner, source, settings, *, extra_roots=(), support_regions=(), facing_halves=False):
    """Both complete actual zero-station cross-sections must be internal."""
    regions = [source.region, *support_regions]
    caps = {}
    targets = (*zip((corner.first, corner.second), corner.stations), *((a, 0.) for a in extra_roots))
    for index, (approach, contact) in enumerate(targets):
        extent = min(approach.length, max(2.*approach.width, contact+.25))
        values = _stations(approach, extent, min(.5, settings.spacing), settings.curve_tolerance, 8192)
        values = _canonical_root_stations(values, extent, contact)
        points = np.vstack((approach.point(values, approach.offset+.5*approach.width),
                            approach.point(values, approach.offset-.5*approach.width)[::-1]))
        cap = points[[-1, 0]]
        owner = Key('root_envelope_strip', index)
        if signed_area(points) < 0.:
            points = points[::-1]
        region = Region(owner, points)
        regions.append(region)
        caps[owner] = cap
    by_owner = {region.owner: region for region in regions}
    for loop in RegionUnion().build(regions):
        for span in loop.spans:
            if span.owner not in caps:
                continue
            points = by_owner[span.owner].xy
            a, b = points[span.segment], points[(span.segment+1) % len(points)]
            cap = caps[span.owner]
            if ((np.linalg.norm(a-cap[0]) <= 1e-8 and np.linalg.norm(b-cap[1]) <= 1e-8)
                    or (np.linalg.norm(a-cap[1]) <= 1e-8 and np.linalg.norm(b-cap[0]) <= 1e-8)):
                if np.linalg.norm(b-a)*(span.last-span.first) > 1e-7:
                    if not facing_halves or span.owner.id>=2:return False
                    direction=cap[1]-cap[0];squared=float(direction@direction)
                    endpoints=np.array((a+span.first*(b-a),a+span.last*(b-a)))
                    parameters=(endpoints-cap[0])@direction/squared
                    low,high=(0.,.5) if span.owner.id==0 else (.5,1.)
                    if (min(high,float(parameters.max()))-max(low,float(parameters.min())))*np.sqrt(squared)>1e-7:
                        return False
    return True


def _validate(corner, settings, minimum_stations, *, dimensional=False):
    approaches = corner.first, corner.second
    preferred = corner.requested_radius if corner.requested_radius is not None else corner.radius
    if (corner.radius_mode != 'STYLE' or (corner.center is None and not dimensional)
            or abs(corner.gap-math.pi) > math.radians(15.)
            or not math.isfinite(preferred) or preferred <= 0.
            or len(minimum_stations) != 2
            or any(not math.isfinite(s) or s < 0. for s in minimum_stations)
            or not all(math.isfinite(v) and v > 0. for v in (settings.spacing, settings.curve_tolerance))):
        raise ValueError('Finite original root envelopes require a near-through STYLE root failure')
    if np.linalg.norm(approaches[0].origin-approaches[1].origin) > 1e-8:
        raise ValueError('Root envelope needs one declared original source junction')
    if any(abs(a.offset) > 1e-8 or a.profile.half or abs(a.width-a.source.width) > 1e-8
           or getattr(a.path, 'absorbed_sites', ())
           or getattr(a.path, 'mouth_piece', None) is not None
           or getattr(a.path, 'deferred_profile_reserve', None) is not None for a in approaches):
        raise ValueError('Root envelope requires centred unmodified original full-profile approaches')
    if abs(approaches[0].profile.datum-approaches[1].profile.datum) > 1e-8:
        raise ValueError('Root envelope requires one original minimum datum')
    return preferred


def _search(owner, corner, settings, stations, **options):
    from .updates import drain_steps
    return drain_steps(_search_steps(owner, corner, settings, stations, **options))


def _search_steps(owner, corner, settings, stations, *, anchor=None, root_audit=None, root_audit_steps=None, admission=None,
            support_width=0., minimum_radius=None, root_points=(), outward_limit=None):
    approaches = corner.first, corner.second
    preferred = corner.requested_radius if corner.requested_radius is not None else corner.radius
    if minimum_radius is not None:preferred=float(minimum_radius)
    root_audit = _complete_roots if root_audit is None else root_audit
    stations = list(stations)
    if not math.isfinite(support_width) or support_width<0.:
        raise ValueError('Root envelope support width must be finite and nonnegative')
    reach_width=max(a.width for a in approaches)
    if support_width>0.:reach_width=max(reach_width,support_width)
    increment = max(.125, max(preferred, .5*reach_width)/8.)
    width = min(a.width for a in approaches)
    maximum_outward=.25*width
    if outward_limit is not None:
        bound=max(maximum_outward,.5*max(0.,support_width-width))
        if (not root_points or not math.isfinite(outward_limit)
                or outward_limit<0. or outward_limit>bound):
            raise ValueError('Explicit root envelope limit needs actual root projections and its source-width bound')
        maximum_outward=float(outward_limit)
    amounts = (0., .005, .01, .025, .05, .1, .2, .25)
    attempts = 0
    for step in range(17):
        exhausted = [(a, s) for a, s in zip(approaches, stations)
                     if not IntersectionProvider.mouth_fits(a, s)]
        if exhausted:
            raise ApproachExtentError(owner, exhausted)
        a, b = approaches
        first, va = a.jet(stations[0], -.5*a.width)
        last, vb = b.jet(stations[1], .5*b.width)
        try:
            base = WidthContinuation(first, last, -va, vb, preferred, refine_curvature=True)
        except ContinuationReachError:
            base = None
        if base is not None:
            outward = right(last-first)
            outward /= np.linalg.norm(outward)
            candidates=amounts
            required=0.
            if root_points:
                axis=(last-first)/np.linalg.norm(last-first)
                if np.min(base.first_derivative@axis)<=1e-8:
                    raise ContinuationReachError('Original width root needs a positive complete chord projection')
                for point in root_points:
                    target=float((point-first)@axis)
                    if not 0.<target<float((last-first)@axis):continue
                    low,high=0.,1.
                    for _ in range(60):
                        middle=.5*(low+high)
                        if float((_evaluate(base.controls,middle)-first)@axis)<target:low=middle
                        else:high=middle
                    t=.5*(low+high);bell=64.*t**3*(1.-t)**3
                    required=max(required,float((point-_evaluate(base.controls,t))@outward)/bell)
                margin=settings.curve_tolerance
                distances=(0.,required+margin,required+2.*margin,required+4.*margin,
                    required+8.*margin,required*1.01+margin,required*1.05+margin,required*1.1+margin)
                candidates=tuple(sorted({distance/width for distance in distances if 0.<=distance<=maximum_outward}))
            for fraction in candidates:
                attempts += 1
                yield
                try:
                    envelope = CurbEnvelope(base, outward, fraction*width, preferred)
                    proxy = replace(corner, stations=tuple(stations), center=None, radius=preferred,
                                    continuation=True, fit_policy='FINITE_ORIGINAL_ROOT_CURB_ENVELOPE')
                    curve = DimensionalCompositeCurb(proxy,
                        (_CertifiedContinuation(EnvelopeBoundary(proxy, envelope)),), settings.curve_tolerance)
                    source = yield from _patch_steps(owner, proxy, curve, tuple(stations), settings)
                    if root_audit_steps is None:
                        accepted = root_audit(proxy, source, settings)
                    else:
                        accepted = yield from root_audit_steps(proxy, source, settings)
                    if not accepted:
                        continue
                except ContinuationReachError:
                    continue
                report = dict(policy='FINITE_ORIGINAL_ROOT_CURB_ENVELOPE', mode='STYLE',
                    admission=admission or 'NEAR_THROUGH_ORIGINAL_SOURCE_ROOT',
                    first_section=a.source.id, first_end=a.end, second_section=b.source.id, second_end=b.end,
                    source_widths=[a.width, b.width], source_profile_halves=[a.profile.half, b.profile.half],
                    original_facing_curb_heights=[float(a.profile.value(-1.)), float(b.profile.value(1.))],
                    stations=list(stations),
                    preferred_minimum_radius=preferred, curvature_bound=envelope.curvature_bound,
                    outward_envelope_m=fraction*width, maximum_outward_envelope_m=maximum_outward,
                    reach_steps=step, candidate_count=attempts, complete_original_root_caps=True,
                    original_sources_unchanged=True, endpoint_extrapolation=False)
                if support_width>0.:
                    report.update(required_source_support_width=support_width,
                                  finite_reach_increment=increment)
                if root_points:
                    report.update(outward_amount_policy='EXACT_ORIGINAL_ROOT_PROJECTION_AND_BOUNDED_CHORD_MARGIN',
                        root_projection_required_outward_m=required,original_root_points=[p.tolist() for p in root_points],
                        candidates_per_reach_at_most=8)
                if outward_limit is not None:
                    report['outward_limit_policy']='BOUNDED_BY_ACTUAL_SUPPORT_WIDTH_MINUS_DONOR_WIDTH'
                return RootEnvelope(source, curve, tuple(stations), report, proxy)
        stations = [IntersectionProvider.mouth_station(a, s+(0. if i == anchor else increment), settings)
                    for i, (a, s) in enumerate(zip(approaches, stations))]
    raise ContinuationReachError('Original root curb envelope exhausts its bounded finite source search')


def prepare_dimensional(owner,corner,settings,failure):
    """Spend both finite original leads after a typed narrow-lead exhaustion.

    The old anchored width continuation remains first. Its accepted radius,
    source widths and endpoint jets stay fixed; the existing zero-jet outward
    envelope must additionally prove both complete original root caps.
    """
    from .compound_dimensional_corners import DimensionalLeadExhausted
    if not isinstance(failure,DimensionalLeadExhausted):
        raise ValueError('Finite width root support requires typed narrow-lead exhaustion')
    approaches=corner.first,corner.second
    if (corner.center is not None or not corner.continuation or corner.radius_mode not in ('STYLE','STRICT')
            or not math.isfinite(corner.radius) or corner.radius<=0.
            or abs(corner.gap-math.pi)>math.radians(15.)
            or abs(approaches[0].width-approaches[1].width)<=1e-6
            or np.linalg.norm(approaches[0].origin-approaches[1].origin)>1e-8
            or any(a.profile.half or abs(a.offset)>1e-8 or abs(a.width-a.source.width)>1e-8
                or getattr(a.path,'absorbed_sites',()) or getattr(a.path,'mouth_piece',None) is not None
                or getattr(a.path,'deferred_profile_reserve',None) is not None for a in approaches)):
        raise ValueError('Finite width root support needs unchanged near-through original full-profile sources')
    if abs(approaches[0].profile.datum-approaches[1].profile.datum)>1e-8:
        raise ValueError('Finite width root support needs its original common minimum datum')
    if failure.narrow != int(approaches[1].width<approaches[0].width) or failure.wide != 1-failure.narrow:
        raise ValueError('Finite width root support changed its exhausted original narrower source')
    stations=[IntersectionProvider.mouth_station(a,max(float(s),min(.25,.05*a.width)),settings)
              for a,s in zip(approaches,failure.minimum_stations)]
    stations[failure.narrow]=failure.last_finite_stations[failure.narrow]
    result=_search(owner,corner,settings,stations,anchor=failure.narrow,minimum_radius=corner.radius,
        admission='TYPED_FINITE_NARROW_WIDTH_LEAD_EXHAUSTION',
        root_points=(approaches[0].point(0.,-.5*approaches[0].width),approaches[1].point(0.,.5*approaches[1].width)),
        root_audit=lambda candidate,source,settings:_complete_roots(candidate,source,settings,facing_halves=True))
    result.report.update(mode=corner.radius_mode,accepted_minimum_radius=corner.radius,
        original_requested_radius=corner.requested_radius,previous_anchored_stations=list(failure.stations),
        exhausted_original_section=approaches[failure.narrow].source.id,
        both_finite_original_leads=True,held_finite_original_narrow_station=stations[failure.narrow])
    result.report.update(complete_original_root_caps=False,complete_original_facing_root_halves=True,
        opposite_root_sectors='mandatory unchanged final full-source union and complete-width proof')
    return result


def prepare(owner, corner, settings, minimum_stations=(0., 0.)):
    """Find a finite STYLE envelope after a typed circular root failure.

    STRICT circles remain strict. The source patch still needs the complete
    source union, mesh, drainage, stored-coordinate and overlap qualification.
    """
    _validate(corner, settings, minimum_stations)
    stations = [IntersectionProvider.mouth_station(a,
        max(float(s), max(0., contact)+min(.25, .05*a.width)), settings)
        for a, contact, s in zip((corner.first, corner.second), corner.stations, minimum_stations)]
    return _search(owner, corner, settings, stations)


def prepare_adjacent(owner, corner, settings, *, required_roots, support_corners):
    """Retain an exposed same-node root beside an accepted width transition.

    The caller must identify an actual exposed original cap at that declared
    source junction. Its adjacent circular corner stays unchanged. Only the
    already certified dimensional curve is replaced; the wider root anchor,
    actual endpoint jets and preferred minimum radius remain authoritative.
    """
    from .corner_patch import corner_patch
    preferred = _validate(corner, settings, (0., 0.), dimensional=True)
    record = getattr(corner, 'dimensional_corner', None)
    if (corner.center is not None or not corner.continuation or record is None
            or record.report.get('policy') != 'FINITE_ORIGINAL_WIDTH_CONTINUATION'
            or tuple(record.stations) != tuple(corner.stations)
            or abs(record.report['minimum_radius']-preferred) > 1e-8
            or abs(corner.first.width-corner.second.width) <= 1e-6):
        raise ValueError('Adjacent root support needs its already certified original width transition')
    required_roots, support_corners = tuple(required_roots), tuple(support_corners)
    if not 1 <= len(required_roots) <= 8 or not 1 <= len(support_corners) <= 16:
        raise ValueError('Adjacent root support needs bounded declared source roots and local corners')
    own = {(a.source.id, a.end) for a in (corner.first, corner.second)}
    roots = {(a.source.id, a.end) for a in required_roots}
    if len(roots) != len(required_roots) or own.intersection(roots):
        raise ValueError('Adjacent root support must identify distinct original neighboring approaches')
    for a in required_roots:
        if (np.linalg.norm(a.origin-corner.first.origin) > 1e-8 or abs(a.offset) > 1e-8
                or abs(a.width-a.source.width) > 1e-8 or a.profile.half not in (-1, 0, 1)
                or getattr(a.path, 'absorbed_sites', ()) or getattr(a.path, 'mouth_piece', None) is not None
                or getattr(a.path, 'deferred_profile_reserve', None) is not None):
            raise ValueError('Adjacent root support requires unchanged same-node original source approaches')
    local = []
    incidence = own | roots
    originals = {(a.source.id, a.end): a for a in (corner.first, corner.second, *required_roots)}
    for i, c in enumerate(support_corners):
        if (c.center is None or any((a.source.id, a.end) not in incidence
                or np.linalg.norm(a.origin-corner.first.origin) > 1e-8 for a in (c.first, c.second))):
            raise ValueError('Adjacent support corners must remain circular at the same original source node')
        if any(a.path is not originals[a.source.id, a.end].path
                or a.profile != originals[a.source.id, a.end].profile
                or a.width != originals[a.source.id, a.end].width for a in (c.first, c.second)):
            raise ValueError('Adjacent support corners lost their identical original source geometry and profiles')
        local.append(corner_patch(Key('root_envelope_support', i), c,
                                  settings.spacing, settings.curve_tolerance).region)
    def audit(candidate, source, settings):
        return _complete_roots(candidate, source, settings, extra_roots=required_roots, support_regions=local)
    anchor = int(corner.second.width > corner.first.width)
    result = _search(owner, corner, settings, corner.stations, anchor=anchor, root_audit=audit,
                     admission='EXPOSED_ADJACENT_ORIGINAL_SOURCE_ROOT')
    result.report.update(wider_root_anchor=anchor, wider_root_station=corner.stations[anchor],
        required_original_roots=[dict(section=a.source.id, end=a.end, width=a.width,
                                     profile_half=a.profile.half) for a in required_roots],
        unchanged_circular_support=[dict(first=(c.first.source.id, c.first.end),
                                        second=(c.second.source.id, c.second.end), radius=c.radius)
                                    for c in support_corners])
    return result
