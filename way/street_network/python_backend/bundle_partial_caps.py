"""Original curb continuations for a Bundle with one attached terminal.

The continuation stays outside Bundle membership and retains its own profile.
Only explicit two-port seams may be traversed; another junction is not skipped.
"""
import numpy as np
from .model import End, Key, Traversal

MAX_SECTIONS = 32
MAX_LENGTH = 2000.


class PhysicalJunctionStop(ValueError):
    """The bounded original walk reached an explicitly authored junction."""
    def __init__(self, path, terminal, site):
        super().__init__('Partial Bundle curb continuation reaches another physical junction')
        self.path, self.terminal, self.site = tuple(path), terminal, site


def terminal_chain(network, first, *, require_terminal=True):
    """Walk a bounded original Street to an actual free or grade-scope end."""
    from .bundle_regions import endpoint_sites
    sites, sections = endpoint_sites(network), network.section_index
    street_id = sections[first.section_id].street_id
    path, used, length, at = [], set(), 0., first
    while True:
        if at.section_id in used or len(path) >= MAX_SECTIONS:
            raise ValueError('Partial Bundle curb continuation exceeds its bounded Section walk')
        section = sections[at.section_id]
        if section.street_id != street_id:
            raise ValueError('Partial Bundle curb continuation changes its original Street')
        points = np.asarray(section.points)[:, :2]
        length += float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())
        if length > MAX_LENGTH:
            raise ValueError('Partial Bundle curb continuation exceeds its finite length budget')
        used.add(at.section_id)
        path.append(Traversal(at.section_id, at.end == 0))
        last = End(at.section_id, 1-at.end)
        site = sites.get(last)
        if site is None:
            if require_terminal and last.end in section.selection_ends:
                raise ValueError('Partial Bundle curb continuation cannot terminate at an import selection cut')
            return tuple(path), last
        following = [e for e in site.ends if e != last]
        if len(following) != 1 or sections[following[0].section_id].street_id != street_id:
            if require_terminal:
                raise PhysicalJunctionStop(path, last, site)
            return tuple(path), last
        at, = following


def dependencies(network, bundle):
    """Conservative source closure before the actual corner side is known."""
    from .bundle_ends import member_ends
    from .bundle_regions import endpoint_sites
    if len(bundle.members) != 2:
        return set()
    sites = endpoint_sites(network)
    ends = [member_ends(member) for member in bundle.members]
    result = set()
    for end in (0, 1):
        owners = [sites.get(member[end]) for member in ends]
        if sum(owner is not None for owner in owners) != 1:
            continue
        site = next(owner for owner in owners if owner is not None)
        for endpoint in site.ends:
            if endpoint.section_id in bundle.sections:
                continue
            path, _ = terminal_chain(network, endpoint, require_terminal=False)
            result.update(t.section_id for t in path)
    return result


def build(context, first, second, end, ports, *, rim_coordinates=None, rim_source=None):
    """Follow the prescribed inside corner and curb, then close the free gap."""
    from .bundles import Cap, _join, _road_rim
    from .intersections import IntersectionProvider
    from .ownership import CompoundProduct, surface_key
    members = (first, second)
    role = next(i for i, port in enumerate(ports) if port is not None)
    endpoint = members[role].ends[end]
    free_end = members[1-role].ends[end]
    if free_end.end in context.network.section_index[free_end.section_id].selection_ends:
        raise ValueError('Partial Bundle cap cannot close at an import selection cut')
    port = ports[role]
    layout = context[port.site.layout_key]
    index = port.approach_index
    approach = layout.approaches[index]
    q = approach.profile.peak if rim_coordinates is None else rim_coordinates[role] * (1 if endpoint.end == 0 else -1)
    if q not in (-1., 1.):
        raise ValueError('Partial Bundle cap needs an explicit physical inner curb')
    preceding = q == 1.
    side = (index-1) % len(layout.approaches) if preceding else index
    corner = layout.corners[side]
    other = corner.first if preceding else corner.second

    def actual_end(i):
        a = layout.approaches[i]
        return End(*list(a.path.trims(layout.outward_stations[i]).items())[-1][0])

    own_end = actual_end(index)
    other_end = actual_end(next(i for i, a in enumerate(layout.approaches) if a is other))
    stopped = None
    try:
        path, terminal = terminal_chain(context.network, other_end)
    except PhysicalJunctionStop as stop:
        path, terminal, stopped = stop.path, stop.terminal, stop
    source_ids = {t.section_id for t in path}
    if any(source_ids.intersection(bundle.sections) for bundle in context.network.bundles):
        raise ValueError('Partial Bundle curb continuation has competing Bundle ownership')
    product = context.products.get(surface_key(port.site.key))
    if isinstance(product, CompoundProduct):
        curve = product.edge_cap((own_end.section_id, own_end.end), (other_end.section_id, other_end.end),
                                 1. if preceding else 0., 0. if preceding else 1.)
    else:
        t = np.linspace(0., 1., layout.rows+1)
        xy, z = IntersectionProvider.band_surface(layout, side, t, np.zeros_like(t))
        ceiling = max(a.profile.rise for a in layout.approaches)
        curve = np.column_stack((xy, approach.profile.datum + context.network.settings.exaggeration*np.clip(z, 0., ceiling)))
        if preceding:
            curve = curve[::-1]
    # Both the rounded corner and the continuation run from the member outwards.
    oriented_side = -1 if preceding else 1
    rims = []
    for traversal in path:
        coordinate = oriented_side * (1 if traversal.forward else -1)
        rim = (rim_source or _road_rim)(context, traversal.section_id, coordinate)
        rims.append(rim if traversal.forward else rim[::-1])
    curb = _join(*rims)
    free = members[1-role].rim[0 if end == 0 else -1]
    closure = None
    if stopped is not None:
        from .finite_free_closure import prepare
        last = path[-1]
        coordinate = oriented_side*(1 if last.forward else -1)
        curve, curb, connector, closure = prepare(context, free_end, free, curve, curb,
            stopped.terminal, coordinate)
    else:
        count = max(4, 2*int(np.ceil(np.linalg.norm(curb[-1, :2]-free[:2])/(2*context.network.settings.spacing))))
        connector = np.linspace(curb[-1], free, count+1)
    record = dict(member=role, bundle_end=end, site=port.site.key.text(),
        path=[dict(section_id=t.section_id, forward=t.forward) for t in path],
        terminal=dict(section_id=terminal.section_id, end=terminal.end,
                      scope=terminal.end in context.network.section_index[terminal.section_id].scope_ends),
        opposite_terminal=dict(section_id=free_end.section_id, end=free_end.end,
            scope=free_end.end in context.network.section_index[free_end.section_id].scope_ends),
        membership_changed=False)
    if closure is not None:
        record['finite_free_closure'] = closure
        record['terminal_is_closure'] = False
    if closure is not None:
        from .finite_free_closure import FiniteFreeCap
        cap = FiniteFreeCap(curve, curb, connector, free, role)
    elif role == 0:
        cap = Cap(curve, _join(curb, connector), free[None, :], True)
    else:
        cap = Cap(free[None, :], _join(connector[::-1], curb[::-1]), curve[::-1], True)
    cap.continuation = record
    return cap
