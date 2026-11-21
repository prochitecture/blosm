"""A bounded grass closure on published curbs at a genuine free profile plane.

One additional physical-junction curb may be followed, stopping before its
next mouth. No next Street is traversed, source End created, or PV point added.
"""
import numpy as np
from .model import Key
from .ownership import PLAN, CompoundProduct, StandaloneProduct, surface_key
from .separator_ports import PORTS
from .source_mouths import source_mouths
from .bundles import Cap, _join


class FiniteFreeCap(Cap):
    """Only its computational exterior connector permits parity refinement."""
    def __init__(self, curve, curb, connector, free, role, *, refined=False):
        if role == 0:
            super().__init__(curve, _join(curb, connector), free[None, :], True)
        else:
            super().__init__(free[None, :], _join(connector[::-1], curb[::-1]), curve[::-1], True)
        self._curve, self._curb, self._connector = curve, curb, connector
        self._free, self._role, self.refined = free, role, bool(refined)

    def refine_exterior_connector(self):
        if self.refined: raise ValueError('Finite free connector parity refinement is bounded to one point')
        lengths = np.linalg.norm(np.diff(self._connector[:, :2], axis=0), axis=1)
        index = int(np.argmax(lengths))
        point = .5*(self._connector[index]+self._connector[index+1])
        connector = np.insert(self._connector, index+1, point, axis=0)
        result = FiniteFreeCap(self._curve, self._curb, connector, self._free, self._role, refined=True)
        result.continuation = dict(self.continuation,
            finite_free_closure=dict(self.continuation['finite_free_closure'],
                computational_connector_refinement=dict(added_collinear_samples=1,
                    physical_pavement_samples_added=0, original_connector_vertices_preserved=True)))
        return result


def inside_vertex(chain, origin, inward, *, tolerance=1e-8):
    """Choose the existing inside endpoint of one unambiguous plane crossing."""
    points = np.asarray(chain, float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 2 or not np.isfinite(points).all():
        raise ValueError('Finite free closure needs a complete published curb chain')
    distances = (points[:, :2]-origin)@inward
    signs = np.where(abs(distances) <= tolerance, 0, np.sign(distances)).astype(int)
    if signs[0] != 1: raise ValueError('Finite free closure must start inside its genuine free profile plane')
    if np.any((signs[:-1] == 0) & (signs[1:] == 0)):
        raise ValueError('Finite free closure has an ambiguous coplanar curb edge')
    events = []
    for i in range(len(points)-1):
        if signs[i]*signs[i+1] < 0: events.append((i, i+1, signs[i]))
    for i in np.flatnonzero(signs == 0):
        if i == len(points)-1:
            events.append((i-1, i, signs[i-1]))
        elif i == 0 or signs[i-1] == signs[i+1]:
            raise ValueError('Finite free closure has a tangent or ambiguous plane contact')
        else:
            events.append((i, i+1, signs[i-1]) if distances[i] >= 0.
                          else (i-1, i, signs[i-1]))
    if len(events) != 1 or events[0][2] != 1:
        raise ValueError('Finite free closure needs exactly one outward transverse plane crossing')
    node, following, _ = events[0]
    if distances[node] < 0.: node -= 1
    if node < 0: raise ValueError('Finite free closure has no original inside curb vertex')
    chord = float(np.linalg.norm(points[following, :2]-points[node, :2]))
    offset = float(distances[node])
    if offset < 0. or offset > chord+tolerance:
        raise ValueError('Finite free closure exceeds its one published chord plane reserve')
    return int(node), dict(existing_vertex=int(node), plane_offset=offset, source_chord_length=chord,
                          unique_transverse_crossing=True, physical_boundary_samples_added=False)


def _physical_corner(context, endpoint, coordinate, expected):
    """The declared bank selects one actual curb, never a new road connection."""
    port = context[PORTS].resolve(endpoint)
    if port is None or endpoint not in port.site.ends:
        raise ValueError('Finite free closure needs its visible original stopping incidence')
    plan = context[PLAN]
    owner = plan.owner(port.site.key) or port.site.key
    product = context.products.get(surface_key(owner))
    if not isinstance(product, (CompoundProduct, StandaloneProduct)) or product.pavement.mesh is None:
        raise ValueError('Finite free closure needs its completed published physical junction')
    layout = ((product.group.physical_layout or context[port.site.layout_key])
              if isinstance(product, CompoundProduct) else product.layout)
    approach = layout.approaches[port.approach_index]
    station = layout.outward_stations[port.approach_index]
    aliases = [a for a in source_mouths(approach, station) if a.original_incidence and a.endpoint == endpoint]
    if len(aliases) != 1:
        raise ValueError('Finite free closure lost its one declared original incoming mouth')
    first = endpoint.section_id, endpoint.end
    u = .5*(coordinate*(1 if endpoint.end == 0 else -1)+1.)
    if isinstance(product, CompoundProduct):
        li, node, at = product.port_points[(*first, u)]
        size = len(product.rims[li])
        steps = [step for step, edge in ((1, node), (-1, (node-1) % size)) if product.kinds[li][edge] == 'curb']
        if len(steps) != 1: raise ValueError('The declared stopping bank has no unique physical curb continuation')
        step = steps[0];nodes = [node];current = node
        while True:
            edge = current if step == 1 else (current-1) % size
            if product.kinds[li][edge] != 'curb': break
            current = (current+step) % size;nodes.append(current)
            if len(nodes) > min(size, 20000): raise ValueError('Finite physical curb trace exceeds its boundary budget')
        other = [(sid, end, value) for (sid, end, value), (loop, vertex, at) in product.port_points.items()
                 if loop == li and vertex == current and (sid, end) != first and value in (0., 1.)]
        if len(other) != 1: raise ValueError('The next physical mouth of the stopping curb is ambiguous')
        sid, end, value = other[0]
        curve = product.edge_cap(first, (sid, end), u, value)
        if not np.array_equal(curve, product.rims[li][nodes]):
            raise ValueError('Finite physical curb trace disagrees with its published whole-edge cap')
        next_end = sid, end
    else:
        index = port.approach_index
        preceding = u == 1.
        side = (index-1) % len(layout.approaches) if preceding else index
        corner = layout.corners[side]
        other = corner.first if preceding else corner.second
        other_index = next(i for i, a in enumerate(layout.approaches) if a is other)
        other_aliases = source_mouths(other, layout.outward_stations[other_index])
        if len(other_aliases) != 1 or not other_aliases[0].original_incidence:
            raise ValueError('Finite free closure cannot traverse an undeclared composite next mouth')
        next_end = other_aliases[0].endpoint.section_id, other_aliases[0].endpoint.end
        curve = product.edge_cap(first, next_end, u, 0. if preceding else 1.)
    if not np.array_equal(curve[0], expected):
        raise ValueError('The published stopping curb does not exactly meet its original road bank')
    return curve, dict(site=port.site.key.text(), owner=owner.text(), incoming=list(first),
                       next_mouth=list(next_end), next_street_traversed=False)


def prepare(context, free_end, free_point, curve, curb, stopping_end, coordinate):
    from .bundle_regions import endpoint_sites
    from .bundles import _join
    section = context.network.section_index[free_end.section_id]
    if (free_end in endpoint_sites(context.network) or free_end.end in section.selection_ends
            or free_end.end in section.scope_ends):
        raise ValueError('Finite physical-curb closure needs a genuinely free original end')
    path = context[Key('path', section.id)]
    profile = context[Key('profile', section.id)]
    if profile.peak not in (-1., 1.): raise ValueError('Finite free closure needs its original inward half-profile')
    station = 0. if free_end.end == 0 else path.length
    origin, tangent, normal = path.frame(station)
    expected = np.r_[path.point(station, section.offset+.5*section.width*profile.peak),
                     profile.datum+context.network.settings.exaggeration*profile.value(profile.peak)]
    if not np.array_equal(expected, free_point):
        raise ValueError('Finite free closure changed its complete original free profile endpoint')
    inward = tangent*(1. if free_end.end == 0 else -1.)
    combined = _join(curve, curb)
    extra = None
    if float(np.min((combined[:, :2]-origin)@inward)) > 1e-8:
        corner, extra = _physical_corner(context, stopping_end, coordinate, curb[-1])
        combined = _join(combined, corner)
        curb = _join(curb, corner)
    node, proof = inside_vertex(combined, origin, inward)
    if extra is not None and node >= len(combined)-1:
        raise ValueError('Finite free closure must stop strictly before the next physical mouth')
    if node < len(curve):
        curve = curve[:node+1]
        curb = curve[-1:]
    else:
        curb = curb[:node-len(curve)+2]
    last = curb[-1]
    count = max(4, 2*int(np.ceil(np.linalg.norm(last[:2]-free_point[:2]) /
                                (2.*context.network.settings.spacing))))
    connector = np.linspace(last, free_point, count+1)
    return curve, curb, connector, dict(proof, policy='genuine_free_profile_plane_on_published_curb',
        free_section=free_end.section_id, free_end=free_end.end, plane_origin=origin.tolist(),
        inward=inward.tolist(), extra_physical_curb=extra, new_source_end=False,
        new_physical_mouth=False, membership_changed=False)
