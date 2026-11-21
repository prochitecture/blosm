"""Explicit junction ownership of a tiny, finite, true dead-end carriageway.

A pad already inside the through-road junction has no independent incoming
profile. Its complete XY width is contained in one bounded smooth through pool;
the two actual entrances retain their full profiles and minimum curb radius.
"""
from dataclasses import dataclass
import math
import numpy as np
from .continuations import WidthContinuation
from .curb_envelope import CurbEnvelope
from .model import End, Key
from .finite_containment import contains_polygon
from .regions import Region, RegionUnion, signed_area


@dataclass(frozen=True)
class OwnedFreePad:
    section: int
    end: int
    length: float
    proof_margin: float
    maximum_envelope: float

    def report(self):
        return dict(section=self.section, end=self.end, consumed_interval=[0., self.length],
            policy='FINITE_FREE_PAD_IN_THROUGH_POOL', containment_margin=self.proof_margin,
            maximum_curb_envelope=self.maximum_envelope, endpoint_extrapolation=False,
            independent_terminal_profile=False)


class EnvelopeBoundary:
    def __init__(self, corner, curve):
        self.corner, self.continuation = corner, curve
        self.length = curve.length
        self.first, self.last = curve(0.), curve(1.)

    def __call__(self, t):
        return self.continuation(t)


def _candidate(context, approaches, failure):
    from .intersections import CornerConstructionError, ApproachExtentError
    from .profiles import memberships
    if not isinstance(failure, (CornerConstructionError, ApproachExtentError)) or len(approaches) != 3:
        return None
    members = memberships(context.network)
    failed = {a.source.id for a in failure.approaches}
    result = []
    for pad in approaches:
        if (pad.source.id not in failed or not getattr(pad.path, 'terminal_free', False)
                or len(pad.path.pieces) != 1 or pad.source.scope_ends or pad.source.selection_ends
                or len(pad.source.points) != 2
                or (not pad.path.pieces[0].path.exact
                    and pad.path.pieces[0].path.curvature_bound > 1e-10)
                or pad.source.id in members or pad.profile.half or abs(pad.offset) > 1e-8
                or getattr(pad.path, 'absorbed_sites', ())):
            continue
        through = [a for a in approaches if a is not pad]
        a, b = through
        if (a.source.id in members or b.source.id in members or a.profile.half or b.profile.half
                or any(getattr(x.path, 'absorbed_sites', ())
                       or getattr(x.path, 'deferred_profile_reserve', None) is not None
                       or getattr(x.path, 'minimum_mouth', 0.) > 1e-9 for x in (a, b))
                or np.linalg.norm(a.origin-pad.origin) > 1e-8 or np.linalg.norm(b.origin-pad.origin) > 1e-8
                # Eligibility only bounds the local search. The complete
                # original pad must still fit inside certified smooth curbs.
                # A narrow stub can be longer than half its own width while
                # remaining entirely inside a wider, moderately bent road.
                or pad.width > min(a.width, b.width)+1e-8
                or pad.length > .5*min(a.width, b.width)
                or float(a.direction @ b.direction) > -math.cos(math.radians(60.))
                or max(abs(a.profile.datum-b.profile.datum), abs(a.profile.datum-pad.profile.datum)) > 1e-8):
            continue
        result.append((pad, a, b))
    return result[0] if len(result) == 1 else None


def from_failure(context, junction, approaches, failure):
    selected = _candidate(context, approaches, failure)
    if selected is None or any(junction.corner_radius(End(a.source.id, a.end))[1] != 'STYLE' for a in approaches):
        return None
    pad, a, b = selected
    from .intersections import Corner, JunctionLayout, Mouth, half_intervals
    from .transitions import _certify_orientation
    radius = max(junction.corner_radius(End(x.source.id, x.end))[0] for x in (a, b))
    settings = context.network.settings
    radial = half_intervals(settings)
    # Original finite pad footprint, including its complete transverse extent.
    source_rim = np.array((pad.point(0., -.5*pad.width), pad.point(pad.length, -.5*pad.width),
                           pad.point(pad.length, .5*pad.width), pad.point(0., .5*pad.width)))
    direction = b.direction-a.direction
    direction /= np.linalg.norm(direction)
    margin = max(4.*settings.curve_tolerance, .02)
    extent = float(np.max(abs((source_rim-pad.origin) @ direction)))+margin+.25
    maximum = min(.4*a.length, .4*b.length)
    for factor in (1., 1.5, 2., 3.):
        distance = extent*factor
        if distance >= maximum-1e-7:
            continue
        try:
            bases = []
            for q in (-1., 1.):
                first, va = a.jet(distance, a.offset-q*.5*a.width)
                last, vb = b.jet(distance, b.offset+q*.5*b.width)
                bases.append(WidthContinuation(first, last, -va, vb, radius, refine_curvature=True))
            mids = np.array([c(.5) for c in bases])
            middle = mids.mean(axis=0)
            outward = [(p-middle)/np.linalg.norm(p-middle) for p in mids]
            side = int(np.argmax([(pad.point(.5*pad.length, 0.)-middle) @ n for n in outward]))
            for amount in (0., .025, .05, .1, .2, .25):
                envelope = amount*min(a.width, b.width)
                curves = [CurbEnvelope(base, normal, envelope if i == side else 0., radius)
                          for i, (base, normal) in enumerate(zip(bases, outward))]
                _certify_orientation(*curves)
                count = max(c.intervals(min(settings.spacing, min(a.width,b.width)/radial),
                                        settings.curve_tolerance) for c in curves)
                if count > 2048:
                    raise ValueError('Finite free-pad curb exceeds its sampling budget')
                t = np.linspace(0., 1., count+1)
                rim = np.r_[curves[0](t), curves[1](t[::-1])]
                certified = RegionUnion().build((Region(Key('free_pad_candidate', junction.id),
                    rim if signed_area(rim) > 0. else rim[::-1]),))
                if len(certified) != 1 or certified[0].hole:
                    continue
                # Curbs are certified to the same chord tolerance. Requiring
                # an extra margin avoids accepting a pad outside the curve.
                if not contains_polygon(source_rim, rim, tolerance=1e-8):
                    continue
                clearance = min(_point_edge_distance(source_rim, rim),
                                _point_edge_distance(rim, source_rim))
                if clearance < margin:
                    continue
                corners = (Corner(a, b, (distance, distance), None, radius, 0., 0., math.pi, continuation=True),
                           Corner(b, a, (distance, distance), None, radius, 0., 0., math.pi, continuation=True))
                # Boundary zero leaves a's left edge and reaches b's right;
                # the second is the other curve traversed in reverse.
                boundaries = (EnvelopeBoundary(corners[0], curves[1]),
                              EnvelopeBoundary(corners[1], ReverseCurve(curves[0])))
                mouths = {}
                for approach in (a, b):
                    mouths.update({identity:Mouth(float(station), 2*radial)
                                   for identity, station in approach.path.trims(distance).items()})
                mouths[pad.source.id, pad.end] = Mouth(pad.path.length if pad.end == 0 else 0., 2*radial)
                claim = OwnedFreePad(pad.source.id, pad.end, pad.length, clearance, envelope)
                return JunctionLayout((a, b), corners, boundaries, (distance, distance), mouths,
                    max(8, 4*int(math.ceil(count/4.))), radial, middle.copy(), [0., 0.],
                    mesh_strategy='boundary_chart', finite_free_pads=(claim,))
        except ValueError:
            continue  # An invalid candidate never changes the original source.
    return None


class ReverseCurve:
    def __init__(self, curve):
        self.curve = curve
        self.length = curve.length

    def __call__(self, t):
        return self.curve(1.-np.asarray(t))

    def intervals(self, spacing, tolerance):
        return self.curve.intervals(spacing, tolerance)


def _point_edge_distance(points, polygon):
    """With nonintersecting boundaries, endpoint projections give exact clearance."""
    edges = np.roll(polygon, -1, axis=0)-polygon
    diffs = points[:, None]-polygon[None]
    u = np.clip(np.sum(diffs*edges, axis=-1)/np.sum(edges*edges, axis=-1), 0., 1.)
    return float(np.linalg.norm(diffs-u[..., None]*edges, axis=-1).min())
