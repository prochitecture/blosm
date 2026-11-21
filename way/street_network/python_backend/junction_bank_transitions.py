"""Finite height transitions on proven original inner banks inside a junction.

This opt-in condition policy changes no physical geometry, radius, source End,
profile plane or ownership. Original source coordinates, cyclic incidence and
the opposite bank's actual contact bound every selected native interval.
"""
from dataclasses import dataclass
import math
import numpy as np
from .conditions import BoundaryData, CurbEdge
from .model import End, Key
from .source_identity import matches_current_section


class JunctionBankTransitionUnavailable(ValueError):
    """An infeasible mixed chain lacks a proved finite original bank domain."""
    def __init__(self, message, *, scope=None):
        self.scope = dict(scope or {}, reason=message)
        super().__init__(message)


@dataclass(frozen=True)
class TransitionBoundaryData(BoundaryData):
    source_transitions: tuple = ()


@dataclass(frozen=True)
class OriginalJunctionBankDomain:
    junction: object
    approach: object
    side: int
    maximum_station: float
    mouth_station: float | None
    opposite_owner: object
    native_path: object
    native_forward: bool
    far_junction: object = None
    opposite_station: float | None = None
    far_contact_owner: object = None

    @classmethod
    def build(cls, resolver, junction, approach, side):
        from .intersections import Approach
        from .paths import Path
        from .routes import Route
        if type(approach) is not Approach or type(approach.path) not in (Path, Route):
            raise JunctionBankTransitionUnavailable('Inner-bank transition needs an original Approach')
        endpoint = End(approach.source.id, approach.end)
        if endpoint.end in approach.source.scope_ends or endpoint.end in approach.source.selection_ends:
            raise JunctionBankTransitionUnavailable('Scope or selection cuts cannot supply a junction-bank domain')
        node = next((item for item in resolver.context.network.intersections
                     if item.key == junction), None)
        if node is None or endpoint not in node.ends or junction not in resolver.group.members:
            raise JunctionBankTransitionUnavailable('Inner-bank transition lost its original junction incidence')
        profile = approach.profile
        if (not profile.half or side not in (-1, 1)
                or float(profile.value(side)) != profile.rise or profile.rise <= 0.
                or float(profile.derivative(side)) != 0.):
            raise JunctionBankTransitionUnavailable('Only a true flat inner half-profile bank can transition')
        if type(approach.path) is Route:
            piece = approach.path.pieces[0]
            native, forward = piece.path, piece.forward
            if (piece.section.id != endpoint.section_id or forward != (endpoint.end == 0)
                    or any(p.path is not resolver.context[Key('path', p.section.id)]
                           for p in approach.path.pieces)):
                raise JunctionBankTransitionUnavailable('Inner-bank transition changed original Route ancestry')
        else:
            native, forward = approach.path, approach.end == 0
        if (native is not resolver.context[Key('path', endpoint.section_id)]
                or profile != resolver.context[Key('profile', endpoint.section_id)].oriented(endpoint.end)
                or approach.width != native.section.width
                or hasattr(approach.path, 'mouth_piece') or getattr(approach, 'absorbed_sites', ())):
            raise JunctionBankTransitionUnavailable('Inner-bank transition changed native dimensions or profile')
        mouths, contacts = [], []
        for source in resolver.sources.values():
            layout = source.layout
            for candidate, station in zip(getattr(layout, 'approaches', ()),
                                          getattr(layout, 'outward_stations', ())):
                if (type(candidate) is Approach and candidate.source.id == endpoint.section_id
                        and candidate.end == endpoint.end and same_source(candidate, approach)):
                    mouths.append(float(station))
            if getattr(layout, 'original_junction', None) != junction:
                continue
            if source.region.owner.kind != 'cycle_corner' or getattr(layout, 'source_bank_fillet', None):
                continue
            for corner in getattr(layout, 'corners', ()):
                if not original_corner(resolver, node, corner):
                    continue
                for candidate, station, bank in ((corner.first, corner.stations[0], -1),
                                                  (corner.second, corner.stations[1], 1)):
                    if (bank == -side and station > 0. and type(candidate) is Approach
                            and candidate.source.id == endpoint.section_id
                            and candidate.end == endpoint.end and same_source(candidate, approach)):
                        other = corner.second if bank == -1 else corner.first
                        if End(other.source.id, other.end) not in node.ends:
                            continue
                        contacts.append((float(station), source.region.owner))
        if not mouths or not contacts or max(mouths)-min(mouths) > 1e-8:
            raise JunctionBankTransitionUnavailable('Inner bank needs its actual fixed mouth and opposite-bank contact')
        contact, owner = min(contacts, key=lambda row: (row[0], row[1]))
        mouth = min(mouths)
        maximum = min(contact, mouth, native.length)
        if not 0. < maximum <= mouth <= approach.length:
            raise JunctionBankTransitionUnavailable('Inner bank has no finite original junction-only interval')
        return cls(junction, approach, side, maximum, mouth, owner, native, forward)

    def native_interval(self, first, last):
        at = np.array((first, last))
        if not self.native_forward:
            at = self.native_path.length-at
        return tuple(sorted(map(float, at)))


def same_source(a, b):
    from .routes import Route
    if a.source.id != b.source.id or a.end != b.end or a.profile != b.profile:
        return False
    if type(a.path) is Route and type(b.path) is Route:
        return (len(a.path.pieces) == len(b.path.pieces)
                and all(x.path is y.path and x.forward == y.forward
                        for x, y in zip(a.path.pieces, b.path.pieces)))
    return a.path is b.path


def original_corner(resolver, node, corner):
    """Both incidences and every native control must belong to the real site."""
    from .intersections import Approach
    from .paths import Path
    from .routes import Route
    if getattr(corner, 'radius_mode', None) not in ('STYLE', 'STRICT'):
        return False
    endpoint = End(corner.first.source.id, corner.first.end)
    if endpoint not in node.ends:
        return False
    requested, mode = node.corner_radius(endpoint)
    declared = getattr(corner, 'requested_radius', None)
    if (mode != corner.radius_mode or declared is None
            or abs(float(declared)-requested) > 1e-9*max(1., requested)
            or mode == 'STRICT' and corner.radius < requested-1e-9*max(1., requested)):
        return False
    for a in (corner.first, corner.second):
        if (type(a) is not Approach or End(a.source.id, a.end) not in node.ends
                or not matches_current_section(resolver.context, a.source)
                or type(a.path) not in (Path, Route)
                or a.profile != resolver.context[Key('profile', a.source.id)].oriented(a.end)
                or a.width != a.source.width
                or a.offset != a.source.offset*(1 if a.end == 0 else -1)):
            return False
        paths = tuple(p.path for p in a.path.pieces) if type(a.path) is Route else (a.path,)
        if any(p is not resolver.context[Key('path', p.section.id)] for p in paths):
            return False
    return True


def terminal_contact_arc(resolver, spans, rules, join, closed):
    """Find the intact generated fillet meeting its own original paired bank."""
    from .compounds import rim_label
    from .routes import Route
    from .source_curb_conditions import _LinearTerminalRule, _SourceRule
    size=len(spans)
    for direction in (-1,1):
        arc=(join-1)%size if direction==-1 else join
        terminal=join if direction==-1 else (join-1)%size
        rule,field=rules[arc],rules[terminal]
        if not (isinstance(rule,_SourceRule) and rule.corner_fragment
                and isinstance(field,_LinearTerminalRule) and rule.z0==rule.z1):
            continue
        owner=spans[arc].owner
        layout=resolver.sources[owner].layout
        if owner.kind!='cycle_corner' or len(layout.corners)!=1:
            continue
        junction=getattr(layout,'original_junction',None)
        site=next((n for n in resolver.context.network.intersections if n.key==junction),None)
        corner=layout.corners[0]
        if (site is None or junction not in resolver.group.members or corner.center is None
                or not original_corner(resolver,site,corner)):
            continue
        _,merged,bank=field.original
        path=merged.paths[bank.role]
        t=field.first if direction==-1 else field.last
        coordinate=field.stations[0. if direction==-1 else 1.]
        if not field.common[0]<=coordinate<=field.common[1]:
            continue
        station=bank.first+t*(bank.last-bank.first)
        contact=False
        for arm,side,at in zip((corner.first,corner.second),(-1,1),corner.stations):
            native=arm.path
            if type(native) is Route:
                first_piece=native.pieces[0]
                if first_piece.forward!=(arm.end==0):
                    continue
                native=first_piece.path
            if (native is path and arm.source.id==bank.section and 0.<=at<=path.length
                    and side*(1 if arm.end==0 else -1)==bank.side
                    and abs(station-(at if arm.end==0 else path.length-at))<=1e-8):
                contact=True
                break
        if not contact:
            continue
        # Only the arc between its tangent contacts may reconcile the two
        # crossfalls. Lead-ins and all adjoining source banks retain their jets.
        label=rim_label(spans[arc],resolver.sources)
        indices=[]
        i=arc
        while len(indices)<size and spans[i].owner==owner:
            piece=rules[i]
            edge=rim_label(spans[i],resolver.sources)
            if (edge.index!=label.index or not isinstance(piece,_SourceRule)
                    or not piece.corner_fragment or piece.phase!=rule.phase
                    or piece.first<rule.phase[0]-1e-10 or piece.last>rule.phase[1]+1e-10):
                break
            indices.append(i)
            i+=direction
            if not closed and not 0<=i<size:
                break
            i%=size
        if indices and 0<=i<size and i not in indices and rules[i] is not None:
            first,last=sorted((indices[0],indices[-1]),key=lambda index:rules[index].first)
            z,gradient=resolver._endpoint(rules[i],1. if direction==-1 else 0.)
            if (abs(float(z)-rule.z0)<=1e-10 and np.linalg.norm(gradient)<=resolver.maximum_slope
                    and abs(rules[first].first-rule.phase[0])<=1e-10
                    and abs(rules[last].last-rule.phase[1])<=1e-10):
                return indices
    return ()


def gradient_maximum(rule):
    """All extrema of an existing CurbEdge polynomial, with no sampling gate."""
    if not isinstance(rule, CurbEdge):
        return 0.
    coefficients = rule.polynomial()[:, 1:]
    norm2 = np.zeros(11)
    for column in coefficients.T:
        values = np.polynomial.polynomial.polymul(column, column)
        norm2[:len(values)] += values
    roots = np.polynomial.polynomial.polyroots(np.polynomial.polynomial.polyder(norm2))
    at = [0., 1., *(float(root.real) for root in roots
                    if abs(root.imag) <= 1e-8 and 0. < root.real < 1.)]
    return math.sqrt(max(0., float(np.max(np.polynomial.polynomial.polyval(at, norm2)))))


def _generated(rule):
    if getattr(rule, 'corner_fragment', False):
        return True
    from .source_curb_conditions import _OffsetRule
    if isinstance(rule, _OffsetRule):
        at = np.interp(np.array((rule.first, rule.last))*rule.curve.length,
                       rule.curve.cumulative, rule.curve.parameters)
        return bool(max(at) <= 0. and min(at) < 0.)
    return False


def _bank_span(resolver, span, label, rule, junction):
    from .source_curb_conditions import _NativeCoordinateRule, _OffsetRule
    from .source_curb_coordinates import CurbCoordinates
    curve = resolver.sources[span.owner].layout.boundaries[label.index]
    if isinstance(rule, _OffsetRule):
        native = rule.curve
        approach = native.approach
        side = int(round(rule.q))
        at = np.interp(np.array((rule.first, rule.last))*native.length,
                       native.cumulative, native.parameters)
        knots = native.cumulative/native.length
        endpoint_knot = lambda t: bool(np.any(abs(knots-t) <= 1e-12))
        compulsory = (endpoint_knot(rule.first), endpoint_knot(rule.last))
    else:
        coordinates = getattr(curve, 'source_coordinates', None)
        if type(coordinates) is not CurbCoordinates:
            raise JunctionBankTransitionUnavailable('Inner bank has no exact original CurbCoordinates')
        coordinates.validate(curve.points, resolver.context)
        if coordinates.approach is None:
            if not isinstance(rule, _NativeCoordinateRule) or rule.curve is not curve:
                raise JunctionBankTransitionUnavailable('Consumed bank needs its native height and station ledger')
            return _claimed_bank_span(resolver, span, rule, junction, coordinates.path)
        if (len(set(coordinates.sides)) != 1
                or not (np.all(np.diff(coordinates.stations) > 0.)
                        or np.all(np.diff(coordinates.stations) < 0.))):
            raise JunctionBankTransitionUnavailable('Inner bank is not one monotone original bank')
        approach, side = coordinates.approach, int(coordinates.sides[0])
        at = np.interp(np.array((label.first, label.last))*curve.length,
                       curve.stations, coordinates.stations)
        knots = curve.stations/curve.length
        compulsory = tuple(bool(np.any(abs(knots-t) <= 1e-12)) for t in (label.first, label.last))
    domain = OriginalJunctionBankDomain.build(resolver, junction, approach, side)
    if min(at) < 0. or max(at) > domain.maximum_station:
        raise JunctionBankTransitionUnavailable('Inner bank interval crosses its actual original influence boundary')
    z, gradient = rule.sample(np.array((0., 1.)))
    if np.any(z != approach.profile.datum+approach.profile.rise) or np.any(gradient != 0.):
        raise JunctionBankTransitionUnavailable('Inner bank is not its true flat native peak')
    return domain, at, compulsory


def _claimed_bank_span(resolver, span, rule, junction, path):
    """Bound a fully consumed plateau bank by original contacts at both junctions."""
    sid = path.section.id
    claims = [claim for claim in resolver.group.claims if claim.section_id == sid]
    if (span.owner != Key('cycle_strip', sid) or len(claims) != 1
            or path is not resolver.context[Key('path', sid)]):
        raise JunctionBankTransitionUnavailable('Inner bank has no unique owned original Section strip')
    claim = claims[0]
    sites = (claim.first_owner, claim.last_owner)
    if (junction not in sites or sites[0] == sites[1] or not set(sites) <= set(resolver.group.members)
            or claim.first_trim != path.length or claim.last_trim != 0.
            or path.section.scope_ends or path.section.selection_ends):
        raise JunctionBankTransitionUnavailable('Inner bank needs a complete Section between two original junctions')
    local = sites.index(junction)
    nodes = {n.key:n for n in resolver.context.network.intersections if n.key in sites}
    if any(site not in nodes or End(sid, end) not in nodes[site].ends for end,site in enumerate(sites)):
        raise JunctionBankTransitionUnavailable('Consumed bank claim disagrees with original endpoint incidence')
    profile = resolver.context[Key('profile', sid)]
    native_side = rule.records[0][-1]
    if (not profile.half or profile.rise <= 0. or profile.value(native_side) != profile.rise
            or profile.derivative(native_side) != 0.
            or any(p is not path or prof != profile or side != native_side for _,_,p,prof,_,_,side in rule.records)):
        raise JunctionBankTransitionUnavailable('Consumed bank is not one original flat plateau edge')
    direction = np.sign(rule.records[0][5]-rule.records[0][4])
    if (not direction or any(direction*(r[5]-r[4]) <= 0. for r in rule.records)
            or any(abs(a[5]-b[4]) > 1e-9 for a,b in zip(rule.records,rule.records[1:]))):
        raise JunctionBankTransitionUnavailable('Consumed bank must follow one continuous native station interval')
    contacts = {}
    for source in resolver.sources.values():
        site = getattr(source.layout, 'original_junction', None)
        if (site not in nodes or source.region.owner.kind != 'cycle_corner'
                or getattr(source.layout, 'source_bank_fillet', None)):
            continue
        for corner in source.layout.corners:
            for arm, station, bank in ((corner.first, corner.stations[0], -1),
                                       (corner.second, corner.stations[1], 1)):
                if arm.source.id != sid:
                    continue
                if (arm.end != sites.index(site) or not original_corner(resolver, nodes[site], corner)
                        or not 0. <= station <= path.length or (arm.end, bank) in contacts):
                    raise JunctionBankTransitionUnavailable('Consumed bank has ambiguous original corner contacts')
                contacts[arm.end, bank] = (station, source.region.owner, arm)
    if len(contacts) != 4:
        raise JunctionBankTransitionUnavailable('Consumed bank needs both original corner contacts at each end')
    side = native_side if local == 0 else -native_side
    start, _, approach = contacts[local, side]
    opposite, owner, _ = contacts[local, -side]
    far = 1-local
    far_station, far_owner, _ = max((contacts[far, bank] for bank in (-1, 1)), key=lambda row:row[0])
    maximum = min(opposite, path.length-far_station)
    samples = np.array((rule.first, rule.last))*rule.curve.length
    stations = []
    for value, record in zip(samples, (rule.records[0], rule.records[-1])):
        lo, hi, _, _, a, b, _ = record
        stations.append(a + np.clip((value-lo)/(hi-lo), 0., 1.)*(b-a))
    at = np.asarray(stations) if local == 0 else path.length-np.asarray(stations)
    if min(at) < start-1e-8 or max(at) > maximum or not start < maximum:
        raise JunctionBankTransitionUnavailable('Consumed inner-bank interval crosses an original corner contact')
    knots = rule.curve.stations/rule.curve.length
    compulsory = tuple(bool(np.any(abs(knots-t) <= 1e-12)) for t in (rule.first, rule.last))
    domain = OriginalJunctionBankDomain(junction, approach, side, maximum, None, owner, path,
                                        local == 0, sites[far], opposite, far_owner)
    return domain, at, compulsory


def reconcile(resolver, spans, points, closed, *, original=None, before=None, after=None,
              baseline=None, original_error=None):
    """Retry only failed mixed generated chains; return rules and typed domains."""
    from .compounds import rim_label
    count = len(spans)
    if not closed or not 1 <= count <= 20000:
        if original_error is not None:
            raise original_error
        return baseline, ()
    labels = [rim_label(span, resolver.sources) for span in spans]
    raw = [original[i] if label.kind == 'mouth' else
           resolver.source_rule(span.owner, label.index, label.first, label.last)
           for i, (span, label) in enumerate(zip(spans, labels))]
    generated = [_generated(rule) for rule in raw]
    starts = [i for i in range(count) if generated[i] and not generated[i-1]]
    softened, reports, used = set(), [], set()
    for start in starts:
        ids, stop = [], start
        while generated[stop]:
            ids.append(stop); stop = (stop+1) % count
            if len(ids) > count:
                raise JunctionBankTransitionUnavailable('Generated transition exceeds its finite boundary budget')
        previous = (start-1) % count
        za, ga = raw[previous].sample(1.)
        zb, gb = raw[stop].sample(0.)
        if abs(float(za)-float(zb)) <= 1e-10:
            continue
        if baseline is not None and all(gradient_maximum(baseline[i]) <= resolver.maximum_slope for i in ids):
            continue
        owners = {spans[i].owner for i in ids}
        nodes = {getattr(resolver.sources[owner].layout, 'original_junction', None) for owner in owners}
        if None in nodes:
            raise JunctionBankTransitionUnavailable('Mixed generated chain needs declared original junctions')
        junction = next(iter(nodes)) if len(nodes) == 1 else None
        for owner in owners:
            source = resolver.sources[owner]
            original_junction = source.layout.original_junction
            node = next((item for item in resolver.context.network.intersections
                         if item.key == original_junction), None)
            if node is None or original_junction not in resolver.group.members:
                raise JunctionBankTransitionUnavailable('Mixed chain lost its original junction authority')
            corners = getattr(source.layout, 'corners', ())
            if (owner.kind != 'cycle_corner' or getattr(source.layout, 'source_bank_fillet', None)
                    or not corners or not all(original_corner(resolver, node, corner) for corner in corners)):
                raise JunctionBankTransitionUnavailable('Generated fragment lacks actual original corner incidence')
        active = list(ids)
        vectors = [points[(i+1) % count]-points[i] for i in active]
        length = sum(float(np.linalg.norm(vector)) for vector in vectors)
        k = max(float(np.linalg.norm(ga)), float(np.linalg.norm(gb)))
        limit = resolver.maximum_slope
        if not limit > k:
            raise JunctionBankTransitionUnavailable('Immutable source crossfall leaves no finite transition length')
        required = 1.875*abs(float(za)-float(zb))/math.sqrt(limit*limit-k*k)
        extension = None
        if length < required:
            if junction is None:
                raise JunctionBankTransitionUnavailable('A finite bank extension needs one shared original junction')
            for direction, anchor in ((-1, previous), (1, stop)):
                # Only a true nonzero flat inner bank may supply extra domain.
                if float(raw[anchor].sample(.5)[0]) <= min(float(za), float(zb))+1e-10:
                    continue
                added, total, domain, first_at, last_at = [], length, None, None, None
                index = anchor
                for _ in range(count):
                    if index in active or labels[index].kind != 'curb':
                        break
                    try:
                        current, at, compulsory = _bank_span(resolver, spans[index], labels[index], raw[index], junction)
                    except JunctionBankTransitionUnavailable:
                        break
                    if domain is not None and (current.side != domain.side or not same_source(current.approach, domain.approach)):
                        break
                    domain = current
                    added.append(index)
                    total += float(np.linalg.norm(points[(index+1) % count]-points[index]))
                    first_at = float(min(at)) if first_at is None else min(first_at, float(min(at)))
                    last_at = float(max(at)) if last_at is None else max(last_at, float(max(at)))
                    if total >= required and compulsory[0 if direction == -1 else 1]:
                        extension = (added, total, domain, first_at, last_at)
                        break
                    index = (index+direction) % count
                if extension is not None:
                    break
            if extension is None:
                raise JunctionBankTransitionUnavailable('No existing finite inner-bank knot meets the physical slope bound',
                    scope=dict(junction=junction.text(), required_length=required, generated_length=length))
            added, length, domain, first_at, last_at = extension
            active.extend(added)
            reports.append(dict(policy='FINITE_ORIGINAL_JUNCTION_INNER_BANK_TRANSITION',
                junction=junction.text(), section=domain.native_path.section.id,
                source_end=domain.approach.end, side=domain.side,
                native_interval=domain.native_interval(first_at, last_at),
                outward_interval=(first_at, last_at), opposite_contact=domain.maximum_station,
                opposite_owner=domain.opposite_owner.text(), unchanged_mouth=domain.mouth_station,
                required_length=required, total_transition_length=length,
                continuous_gradient_bound=math.hypot(1.875*abs(float(za)-float(zb))/length, k),
                physical_limit=limit, existing_native_knots=True, geometry_profiles_radii_unchanged=True))
            if domain.far_junction is not None:
                reports[-1].update(opposite_contact=domain.opposite_station,
                    maximum_outward_station=domain.maximum_station, far_junction=domain.far_junction.text(),
                    far_contact_owner=domain.far_contact_owner.text(), fully_consumed_section=True)
        if used.intersection(active):
            raise JunctionBankTransitionUnavailable('Adjacent finite transitions need a single explicit source-domain blend')
        used.update(active); softened.update(active)
    if not softened:
        if original_error is not None:
            raise original_error
        return baseline, ()
    result = resolver._chain_rules(spans, points, closed, before=before, after=after,
                                   original=original, soften=softened)
    if any(gradient_maximum(result[i]) > resolver.maximum_slope for i in softened):
        raise JunctionBankTransitionUnavailable('Finite source-bank transition lost its continuous slope certificate')
    return result, tuple(reports)


def certify_body_domains(guides, boundaries):
    """Never remove a real body constraint to make a junction taper fit."""
    domains = [row for boundary in boundaries for row in getattr(boundary, 'source_transitions', ())]
    for domain in domains:
        low, high = domain['native_interval']
        sources = [source for source in guides.ledger['source_influence']
                   if source['section'] == domain['section']
                   and source['owned'][0] <= low <= high <= source['owned'][1]]
        if not guides.ledger['complete_supported_scope'] or len(sources) != 1:
            raise JunctionBankTransitionUnavailable('Junction taper lacks a complete original source-body domain', scope=domain)
        for source in sources:
            if any(min(high, b) > max(low, a) for a, b in source['core_ranges']):
                raise JunctionBankTransitionUnavailable('Junction bank transition overlaps a hard original body core', scope=domain)
    return tuple(domains)
