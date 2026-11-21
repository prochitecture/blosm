"""Original radius authority for contacts reached through owned source roads.

The logical corner remains at its actual OSM junction. Its physical contact
may be on a later, differently sized Section. Validate both facts separately;
never pretend that the later Section has an End at the logical junction.
"""
from dataclasses import dataclass
import math
import numpy as np
from .model import End, Key
from .source_identity import matches_current_section


class DisplacedCornerAuthorityError(ValueError):
    """A recorded displaced contact lost its original finite source authority."""


@dataclass(frozen=True)
class DisplacedCornerAuthority:
    junction: object
    logical_ends: tuple
    contact_ends: tuple
    source_walks: tuple
    attempted_radii: tuple


def _require(condition, message):
    if not condition:
        raise DisplacedCornerAuthorityError(message)


def _native_arm(context, arm):
    from .compound_displaced_contacts import _ordinary
    from .intersections import Approach
    from .routes import Route
    _require(type(arm) is Approach and type(arm.path) is Route and _ordinary(arm),
             'Displaced authority needs an unchanged ordinary source Approach')
    _require(matches_current_section(context, arm.source)
             and arm.profile == context[Key('profile', arm.source.id)].oriented(arm.end)
             and arm.width == arm.source.width
             and arm.offset == arm.source.offset*(1 if arm.end == 0 else -1),
             'Displaced authority changed source dimensions or profile')
    _require(arm.path.pieces and all(piece.path is context[Key('path', piece.section.id)]
             for piece in arm.path.pieces), 'Displaced authority changed native Path ancestry')
    _require(arm.path.pieces[0].section is arm.source
             and arm.path.pieces[0].forward == (arm.end == 0),
             'Displaced authority changed the original directed incidence')


def _same_walk(a, b):
    from .junction_bank_transitions import same_source
    return (same_source(a.original, b.original) and same_source(a.approach, b.approach)
            and len(a.prefix) == len(b.prefix)
            and all(x.path is y.path and x.forward == y.forward
                    for x, y in zip(a.prefix, b.prefix)))


def validate_displaced_corner(context, group, node, corner):
    """Reprove the existing bounded source walk, policy and physical contacts.

    This accepts only the explicit original producer's dimensional-contact
    record. It neither searches for another radius nor alters its geometry.
    """
    from .compound_displaced_contacts import SourceContactWalk, _walks
    from .compound_contact_routes import regular_contact
    from .corner_radii import CornerRadiusPolicy
    from .ownership import footprint_key
    from .source_bank_reconciliation import owned_ranges
    _require(node.key in group.members and any(node is site for site in context.network.intersections),
             'Displaced authority lost its real original junction')
    _require(getattr(corner, 'fit_policy', None) == 'DISPLACED_ORIGINAL_DIMENSIONAL_SOURCE_CURBS'
             and corner.radius_mode == 'STYLE' and corner.center is not None,
             'Displaced authority needs its explicit original STYLE producer')
    walks = getattr(corner, 'source_contact_walks', ())
    _require(type(walks) is tuple and len(walks) == 2
             and all(type(walk) is SourceContactWalk for walk in walks),
             'Displaced authority needs two typed original source walks')
    _require(corner.first is walks[0].approach and corner.second is walks[1].approach
             and any(walk.prefix for walk in walks) and any(walk.changed_width for walk in walks),
             'Displaced authority does not describe its actual dimensional contacts')
    for walk in walks:
        for arm in (walk.original, walk.approach):
            _native_arm(context, arm)
    original = tuple(walk.original for walk in walks)
    logical = tuple(End(arm.source.id, arm.end) for arm in original)
    _require(logical[0] != logical[1] and all(end in node.ends for end in logical)
             and np.linalg.norm(original[0].origin-original[1].origin) <= 1e-8,
             'Displaced authority lost its logical original corner incidences')
    layout = context.products.get(footprint_key(node.key))
    _require(layout is not None, 'Displaced authority lacks the original source layout')
    from .junction_bank_transitions import same_source
    _require(all(any(same_source(arm, candidate) for candidate in layout.approaches)
                 for arm in original), 'Displaced authority substituted the original source arms')
    ordered = sorted(layout.approaches,
        key=lambda arm: math.atan2(arm.direction[1], arm.direction[0]) % math.tau)
    ends = tuple(End(arm.source.id, arm.end) for arm in ordered)
    _require(len(ends) == len(node.ends) and set(ends) == set(node.ends)
             and ends[(ends.index(logical[0])+1) % len(ends)] == logical[1],
             'Displaced authority changed the original counterclockwise corner sector')
    requested, mode = node.corner_radius(logical[0])
    _require(mode == 'STYLE' and corner.requested_radius == requested,
             'Displaced authority changed the original preferred radius policy')
    candidates = tuple(float(value) for value in
                       CornerRadiusPolicy.contact_candidates(*original, requested, mode))
    attempts = tuple(corner.radius_attempts)
    _require(attempts and attempts == candidates[:len(attempts)]
             and corner.radius == attempts[-1],
             'Displaced authority did not use the bounded original radius candidates')
    consumed = {claim.section_id for claim in group.claims}
    members = set(group.members)
    owned = owned_ranges(context, group)
    for walk, side, station in zip(walks, (-1, 1), corner.stations):
        _require(any(_same_walk(walk, candidate) for candidate in
                     _walks(context, walk.original, side, members, consumed)),
                 'Displaced authority lost its actual owned source-incidence walk')
        arm = walk.approach
        _require(math.isfinite(station) and 0. <= station <= arm.length
                 and regular_contact(arm, station),
                 'Displaced authority extrapolated or crossed a native contact seam')
        for (sid, end), at in arm.path.trims(station).items():
            length = context[Key('path', sid)].length
            first, last = (0., at) if end == 0 else (at, length)
            _require(any(a <= first <= last <= b for a, b in owned.get(sid, ())),
                     'Displaced authority left actual finite source ownership')
    u, v = original[0].direction, original[1].direction
    gap = math.atan2(float(u[0]*v[1]-u[1]*v[0]),
                     float(original[0].direction @ original[1].direction)) % math.tau
    side = 1. if gap < math.pi else -1.
    _require(math.isfinite(corner.sweep) and 0. < abs(corner.sweep) <= math.pi+1e-7
             and corner.sweep*side <= 1e-7,
             'Displaced authority selected another circular sector')
    angles = corner.start_angle+corner.sweep*np.array((0., 1.))
    arc = corner.center+corner.radius*np.column_stack((np.cos(angles), np.sin(angles)))
    for i, (arm, station, sign) in enumerate(zip((corner.first, corner.second),
                                               corner.stations, (-1, 1))):
        point, tangent = arm.jet(station, arm.offset+sign*.5*arm.width)
        radial = point-corner.center
        _require(np.isfinite(point).all() and np.isfinite(tangent).all()
                 and np.linalg.norm(tangent) > 1e-10
                 and np.linalg.norm(point-arc[i]) <= 1e-8
                 and abs(np.linalg.norm(radial)-corner.radius) <= 1e-5
                 and abs(float(radial @ tangent))/np.linalg.norm(tangent) <= 1e-7,
                 'Displaced authority lost its actual finite circular tangencies')
    return DisplacedCornerAuthority(node.key, logical,
        tuple(End(arm.source.id, arm.end) for arm in (corner.first, corner.second)),
        walks, attempts)


def original_parent_corner(resolver, node, corner):
    """Existing direct authority, or a fully proved displaced source record."""
    from .junction_bank_transitions import original_corner
    if original_corner(resolver, node, corner):
        return True
    if not hasattr(resolver, 'group'):
        return False
    support = getattr(corner, 'source_contact_support', None)
    if isinstance(support, dict) and support.get('policy') == 'EXPLICIT_WHOLE_OWNED_PREFIX_CONTACTS':
        from .owned_prefix_corner_authority import (
            validate_owned_prefix_corner, OwnedPrefixCornerAuthorityError)
        try:
            validate_owned_prefix_corner(resolver.context, resolver.group, node, corner)
        except OwnedPrefixCornerAuthorityError:
            return False
        return True
    try:
        validate_displaced_corner(resolver.context, resolver.group, node, corner)
    except DisplacedCornerAuthorityError:
        return False
    return True
