"""Source conditions for a proven finite free-pad through-pool envelope."""
from dataclasses import dataclass
import math
import numpy as np
from .model import End, Intersection, Key
from .source_identity import matches_current_section


class FreePadCurbUnavailable(ValueError):
    """A generated pool boundary lost its original finite source certificate."""


@dataclass(frozen=True)
class FreePadCurbDomain:
    junction: object
    claim: object
    corner: object
    heights: tuple
    inward_slopes: tuple


def _require(value, message):
    if not value:
        raise FreePadCurbUnavailable(message)


def validate(context, group, sources, curve):
    """Prove the exact producer, complete owned pad and real incoming jets."""
    from .free_pads import EnvelopeBoundary, ReverseCurve, OwnedFreePad, _point_edge_distance
    from .curb_envelope import CurbEnvelope
    from .compound_curbs import SamplingBudgetCurb
    from .intersections import Approach
    from .routes import Route, RouteTopology
    from .ownership import footprint_key
    from .finite_containment import contains_polygon
    _require(type(curve) is EnvelopeBoundary and len(group.members) == 1 and not group.claims
             and len(group.footprints) == 1 and len(group.loops) == 1 and not group.loops[0].hole,
             'Free-pad source rule needs its single original through pool')
    owner = group.members[0]; source = sources.get(owner)
    node = context.network.records().get(owner)
    _require(source is not None and type(node) is Intersection
             and source.layout is context.products.get(footprint_key(owner))
             and not source.region.boundary_only,
             'Free-pad source rule lost its actual original layout')
    layout = source.layout
    records = getattr(layout, 'finite_free_pads', ())
    _require(len(records) == 1 and type(records[0]) is OwnedFreePad
             and len(layout.approaches) == 2 and len(layout.outward_stations) == 2
             and len(layout.boundaries) == 2 and len(layout.corners) == 2,
             'Free-pad source rule needs its typed whole-pad producer record')
    actual = []
    for boundary in layout.boundaries:
        for _ in range(16):
            if not isinstance(boundary, SamplingBudgetCurb):
                break
            boundary = boundary.curve
        actual.append(boundary)
    _require(sum(boundary is curve for boundary in actual) == 1
             and any(corner is curve.corner for corner in layout.corners),
             'Free-pad source rule borrowed another envelope or corner')
    claim = records[0]; endpoint = End(claim.section, claim.end)
    _require(len(node.ends) == 3 and endpoint in node.ends
             and all(node.corner_radius(end)[1] == 'STYLE' for end in node.ends),
             'Free-pad source rule cannot replace an authored strict corner')
    arms = tuple(layout.approaches)
    _require({End(arm.source.id, arm.end) for arm in arms} == set(node.ends)-{endpoint},
             'Free-pad source rule changed its two real through incidences')
    for arm in arms:
        _require(type(arm) is Approach and type(arm.path) is Route and arm.path.pieces
                 and matches_current_section(context, arm.source)
                 and all(piece.path is context[Key('path', piece.section.id)] for piece in arm.path.pieces)
                 and arm.profile == context[Key('profile', arm.source.id)].oriented(arm.end)
                 and not arm.profile.half and arm.width == arm.source.width
                 and arm.offset == arm.source.offset*(1 if arm.end == 0 else -1)
                 and not getattr(arm.path, 'absorbed_sites', ())
                 and not getattr(arm.path, 'minimum_mouth', 0.) and not hasattr(arm, 'field'),
                 'Free-pad source rule changed an original through profile or Path')
    radius = max(node.corner_radius(End(arm.source.id, arm.end))[0] for arm in arms)
    maxima = []
    for boundary in actual:
        _require(type(boundary) is EnvelopeBoundary,
                 'Free-pad source rule needs its two explicit envelope boundaries')
        corner = boundary.corner; a, b = corner.first, corner.second
        _require(any(a is arm for arm in arms) and any(b is arm for arm in arms) and a is not b
                 and corner.center is None and corner.continuation and corner.radius == radius,
                 'Free-pad source rule changed its original finite curb-radius policy')
        for arm, station in zip((a, b), corner.stations):
            at = next(value for candidate, value in zip(arms, layout.outward_stations) if candidate is arm)
            _require(math.isfinite(station) and station == at and 0. < station < .4*arm.length,
                     'Free-pad source rule moved an actual incoming profile plane')
        leaf = boundary.continuation; reverse = type(leaf) is ReverseCurve
        if reverse:
            leaf = leaf.curve
        _require(type(leaf) is CurbEnvelope and leaf.controls.shape == (7, 2)
                 and 0. <= leaf.curvature_bound <= 1./radius
                 and boundary.length == leaf.length,
                 'Free-pad source rule lost its certified polynomial envelope')
        maxima.append(leaf.envelope_distance)
        for i, (arm, station, side) in enumerate(zip((a, b), corner.stations, (-1, 1))):
            point, tangent = arm.jet(station, arm.offset+side*.5*arm.width)
            parameter = float(i)
            velocity = (-leaf.derivative(1.-parameter) if reverse else leaf.derivative(parameter))
            expected = (-tangent if i == 0 else tangent)
            _require(np.linalg.norm(boundary(parameter)-point) <= 1e-8
                     and np.linalg.norm(velocity) > 1e-10 and np.linalg.norm(expected) > 1e-10
                     and np.linalg.norm(velocity/np.linalg.norm(velocity)-expected/np.linalg.norm(expected)) <= 1e-8,
                     'Free-pad source rule changed a real endpoint or tangent')
    _require(max(maxima) == claim.maximum_envelope and min(maxima) == 0.
             and max(maxima) <= .25*min(arm.width for arm in arms),
             'Free-pad source rule changed its bounded envelope construction')
    pad_path = context[Key('path', claim.section)]; pad = pad_path.section
    trace = RouteTopology(context.network).trace(endpoint)
    participating = {pad.id, *(piece.section.id for arm in arms for piece in arm.path.pieces)}
    _require(matches_current_section(context, pad) and not pad.scope_ends and not pad.selection_ends
             and len(pad.points) == 2 and claim.length == pad_path.length
             and (pad_path.exact or pad_path.curvature_bound <= 1e-10)
             and len(trace.traversals) == 1 and trace.terminal_site is None
             and not trace.terminal_selection and not trace.terminal_scope
             and pad.width <= min(arm.width for arm in arms)+1e-8
             and claim.length <= .5*min(arm.width for arm in arms)
             and not any(participating.intersection(bundle.sections) for bundle in context.network.bundles),
             'Free-pad source rule lost its actual finite unbundled dead end')
    pad_arm = Approach(pad_path, claim.end, context[Key('profile', pad.id)])
    _require(not pad_arm.profile.half and abs(pad_arm.offset) <= 1e-8
             and all(np.linalg.norm(arm.origin-pad_arm.origin) <= 1e-8
                     and abs(arm.profile.datum-pad_arm.profile.datum) <= 1e-8 for arm in arms),
             'Free-pad source rule changed its original node or datum')
    points = np.array((pad_arm.point(0., -.5*pad.width), pad_arm.point(claim.length, -.5*pad.width),
                       pad_arm.point(claim.length, .5*pad.width), pad_arm.point(0., .5*pad.width)))
    mouth = layout.mouths.get((pad.id, claim.end))
    margin = max(4.*context.network.settings.curve_tolerance, .02)
    clearance = min(_point_edge_distance(points, source.region.xy), _point_edge_distance(source.region.xy, points))
    _require(mouth is not None and mouth.station == (pad_path.length if claim.end == 0 else 0.)
             and claim.proof_margin >= margin and contains_polygon(points, source.region.xy, tolerance=1e-8)
             and clearance >= margin-context.network.settings.curve_tolerance,
             'Free-pad source rule lost whole-width pad containment or its complete claim')
    corner = curve.corner
    heights = tuple(arm.profile.datum+float(arm.profile.value(side))
                    for arm, side in zip((corner.first, corner.second), (-1, 1)))
    slopes = tuple(float(np.linalg.norm(arm.physical_gradient(station, side)))
                   for arm, station, side in zip((corner.first, corner.second), corner.stations, (-1, 1)))
    return FreePadCurbDomain(owner, claim, corner, heights, slopes)
