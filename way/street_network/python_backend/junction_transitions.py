"""A physical junction may own a nearby dimensional Section seam.

The original Section records and path polynomials survive. A displaced external
port uses the next Section's full dimensions; no isolated taper is constructed
under the junction. Incidence, not proximity, establishes that ownership.
"""
from dataclasses import dataclass, replace
import numpy as np
from .model import End, Key
from .routes import Piece, Route, RouteTopology


@dataclass(frozen=True)
class DeferredProfileReserve:
    """Provisional raw contact; ownership must prove this route is internal."""
    terminal_site: Key
    section_ids: tuple
    minimum_with_profile: float


class JunctionTransitions:
    maximum_expansions = 8

    @staticmethod
    def seam_signature(context, keys, neighbors):
        """Local input/owner signature despite a network-wide plan task."""
        from .ownership import PLAN
        plan = context[PLAN]
        owners = [(j.key.text(), plan.owner(j.key).text() if plan.owner(j.key) is not None else None)
                  for j in neighbors]
        return [(key.text(), context.digest(key)) for key in keys], owners

    @staticmethod
    def extend(context, approach):
        ordinary = JunctionTransitions._extend(context, approach)
        if ordinary is not None:
            return ordinary
        # An internal crown is not an exposed entrance condition. Only defer
        # its lead requirement provisionally; OwnershipPlan must subsequently
        # prove complete compound ownership before any surface can be built.
        return JunctionTransitions._extend(context, approach, compound_only=True)

    @staticmethod
    def _extend(context, approach, *, compound_only=False):
        """Cross one width/offset/rise change, retaining following pieces.

        A fully consumed crown-rise change is reconciled inside the junction;
        its exterior mouth uses the continuing Section's complete profile.
        Datum or profile-family changes, different Streets and different
        interpolation modes still need their own transition policy.
        """
        path = approach.path
        candidate = RouteTopology(context.network).trace(End(approach.source.id, approach.end))
        count = len(path.pieces)
        if count >= len(candidate.traversals):
            return None
        pieces = list(path.pieces)
        profile = approach.profile
        added = False
        sites = list(getattr(path, 'absorbed_sites', ()))
        minimum = getattr(path, 'minimum_mouth', 0.)
        profile_minimum = minimum
        width, offset = approach.width, approach.offset
        for index, visit in enumerate(candidate.traversals[count:], count):
            source_path = context[Key('path', visit.section_id)]
            source = source_path.section
            next_profile = context[Key('profile', source.id)].oriented(0 if visit.forward else 1)
            next_offset = source.offset * (1 if visit.forward else -1)
            if (source.street_id != approach.source.street_id
                    or source.interpolation != approach.source.interpolation
                    or abs(next_profile.datum - profile.datum) > 1e-9
                    or next_profile.half != profile.half):
                break
            rise_delta = abs(next_profile.rise - profile.rise)
            changed = (abs(source.width - width) > 1e-8 or abs(next_offset - offset) > 1e-8
                       or rise_delta > 1e-9)
            if added and changed:
                break
            piece = Piece(source_path, visit.forward)
            previous = pieces[-1]
            a, ta, na = previous.frame(previous.path.length)
            b, tb, nb = piece.frame(0.)
            reach = max(abs(offset) + width * .5, abs(next_offset) + source.width * .5)
            if (np.linalg.norm(a - b) > 1e-8 or np.linalg.norm(ta - tb) > 1e-8
                    or reach * np.linalg.norm(na - nb) > 1e-8):
                break
            if changed:
                geometric = max(2., 4. * abs(source.width - width) + 2. * abs(next_offset - offset))
                desired = max(geometric, 1.875 * rise_delta / .03)
                prefix = sum(p.path.length for p in pieces)
                margin = min(.25, .05 * source.width)
                profile_minimum = prefix + .5 * desired + margin
                minimum = prefix + .5 * (geometric if compound_only else desired) + margin
                sites.append(candidate.seams[index - 1])
                added = True
            pieces.append(piece)
            width, offset, profile = source.width, next_offset, next_profile
        if not added:
            return None
        # Constant terminal dimensions are used for the corner contacts, also
        # on consumed inner pieces. Certify that reach on each original cubic.
        reach = abs(offset) + width * .5
        if any(not p.path.exact and not p.path._certified(reach) for p in pieces):
            return None
        result = Route(pieces, terminal_free=(len(pieces) == len(candidate.traversals)
            and candidate.terminal_site is None and not candidate.terminal_selection))
        if minimum >= result.length - 1e-6:
            return None
        if compound_only:
            physical = {j.key for j in context.network.intersections if len(j.ends) > 2}
            ids = tuple(visit.section_id for visit in candidate.traversals)
            if (profile.half != 0 or candidate.terminal_site not in physical
                    or candidate.terminal_selection or result.terminal_free
                    or ids != tuple(piece.section.id for piece in pieces)):
                return None
            result.deferred_profile_reserve = DeferredProfileReserve(
                candidate.terminal_site, ids, profile_minimum)
        result.mouth_piece, result.mouth_profile = pieces[-1], profile
        result.minimum_mouth, result.absorbed_sites = minimum, tuple(sites)
        from .intersections import Approach
        return Approach(result, approach.end, context[Key('profile', approach.source.id)])

    @classmethod
    def footprint(cls, context, junction, approaches):
        from .intersections import IntersectionProvider, CornerConstructionError, ApproachExtentError
        current = list(approaches)
        for iteration in range(cls.maximum_expansions + 1):
            try:
                layout = IntersectionProvider._footprint(context, junction, list(current))
            except (CornerConstructionError, ApproachExtentError) as exc:
                needed = {a.source.id for a in exc.approaches}
                layout, failure = None, exc
            else:
                # Do not replace a qualified isolated taper. Expansion is
                # needed only when its complete source span is consumed or the
                # ordinary mandatory mouth clearance leaves no room for it.
                needed = {a.source.id for a, distance in zip(layout.approaches, layout.outward_stations)
                    if a.length - distance <= min(.25, .05 * a.width) + 1e-7}
                failure = None
            replacements = {}
            if iteration < cls.maximum_expansions:
                for a in current:
                    if a.source.id in needed:
                        extended = cls.extend(context, a)
                        if extended is not None:
                            replacements[a.source.id] = extended
            if not replacements:
                if failure is not None:
                    raise failure
                cls.check_claims(context, layout)
                return replace(layout, absorbed_transitions=tuple(sorted({s for a in current
                    for s in getattr(a.path, 'absorbed_sites', ())})))
            current = [replacements.get(a.source.id, a) for a in current]
        raise ValueError('Junction dimensional seam ownership exceeded its expansion budget')

    @staticmethod
    def check_claims(context, layout):
        """Every absorbed incidence must be strictly inside its physical rim.

        The old isolated taper is not a boundary constraint: it is replaced,
        not clipped. The actual exterior mouth and rounded curb remain the
        authoritative boundary. Whole-profile and surface guards run normally.
        """
        sites = {s.key: s for s in (*context.network.intersections, *context.network.connections)}
        claimed = {s for a in layout.approaches for s in getattr(a.path, 'absorbed_sites', ())}
        if not claimed:
            return
        from .compounds import sample_footprint
        from .regions import _classify
        polygon = sample_footprint(Key('claim', 0), layout).region.xy
        for key in sorted(claimed):
            site = sites[key]
            points = np.array([context.network.section_index[e.section_id].points[0 if e.end == 0 else -1][:2]
                               for e in site.ends])
            inside, boundary = _classify(points, polygon, 1e-7)
            if not np.all(inside & (boundary < 0)):
                raise ValueError('Absorbed dimensional seam is not inside its junction boundary')

    @staticmethod
    def absorbed(context, site, approaches, neighbors):
        """Retain the two-ended source site but publish no duplicate taper."""
        from .ownership import PLAN, footprint_key
        from .intersections import JunctionLayout
        claims = [(j, context[footprint_key(j.key)]) for j in neighbors
                  if site.key in context[footprint_key(j.key)].absorbed_transitions]
        if not claims:
            return None
        junction, footprint = claims[0]
        owner = junction.key
        if len(claims) != 1:
            owners = {context[PLAN].owner(j.key) for j, _ in claims}
            if None in owners or len(owners) != 1:
                raise ValueError('Multiple junctions claim a dimensional seam without one physical compound owner')
            owner = next(iter(owners))
        # A prefix ending before the seam cannot consume it. Both original
        # incidences must be present among the physical route trims.
        for _, footprint in claims:
            for end in site.ends:
                if not any(sid == end.section_id for sid, direction in footprint.mouths):
                    raise ValueError('Junction ownership does not reach both dimensional seam incidences')
        return JunctionLayout(tuple(approaches), (), (), (0., 0.), {}, 0, 0,
            approaches[0].origin.copy(), [0., 0.], absorbed_by=owner)
