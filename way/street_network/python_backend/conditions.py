"""Profile conditions on combined junction boundaries, in physical metres.

These are inputs to the forthcoming compound surface provider, not an accepted
surface. Internal, consumed junction mouths impose no phantom height seam.
"""
from dataclasses import dataclass
import numpy as np
from .paths import right
from .compounds import rim_label


def _blend(t):
    return t ** 3 * (10. + t * (-15. + 6. * t))


def _profile_gradient(approach, station, q):
    # Derived entrances may have a longitudinally changing profile. Query the
    # approach's physical field rather than assuming a fixed Section profile.
    return approach.physical_gradient(station, q)


@dataclass(frozen=True)
class ProfileEdge:
    approach: object
    station: float
    first_u: float
    last_u: float

    def sample(self, t):
        u = self.first_u + np.asarray(t) * (self.last_u - self.first_u)
        q = 2. * u - 1.
        return self.approach.profile.datum + self.approach.profile.value(q), _profile_gradient(self.approach, self.station, q)


@dataclass(frozen=True)
class CurbEdge:
    first_height: float
    last_height: float
    first_slope: float
    last_slope: float
    length: float
    first_station: float
    last_station: float
    tangent: tuple

    def sample(self, t):
        fraction = (self.first_station + np.asarray(t) * (self.last_station - self.first_station)) / self.length
        weight = _blend(fraction)
        height = self.first_height + weight * (self.last_height - self.first_height)
        along = (self.last_height - self.first_height) * 30. * fraction ** 2 * (1. - fraction) ** 2 / self.length
        inward = self.first_slope + weight * (self.last_slope - self.first_slope)
        tangent = np.asarray(self.tangent)
        gradient = along[..., None] * tangent - inward[..., None] * right(tangent)
        return height, gradient


@dataclass(frozen=True)
class PlateauEdge:
    height: float

    def sample(self, t):
        t = np.asarray(t)
        return np.full(t.shape, self.height), np.zeros(t.shape + (2,))


@dataclass(frozen=True)
class Port:
    owner: object
    approach_index: int
    section_id: int
    end: int
    station: float
    vertices: tuple
    outward_u: tuple
    aliases: tuple = ()

    @property
    def source_u(self):
        if len(self.aliases) > 1:
            raise ValueError('A shared physical mouth has separate original source parameter maps')
        return self.outward_u if self.end == 0 else tuple(1. - u for u in self.outward_u)


@dataclass(frozen=True)
class BoundaryData:
    xy: np.ndarray
    heights: np.ndarray
    gradients: np.ndarray
    rules: tuple
    ports: tuple
    datum: float
    ceiling: float
    hole: bool

    def edge_samples(self, index, t):
        """Sample a retained chord and its unexaggerated boundary condition."""
        t = np.asarray(t)
        a, b = self.xy[index], self.xy[(index + 1) % len(self.xy)]
        height, gradient = self.rules[index].sample(t)
        return a + t[..., None] * (b - a), height, gradient


class BoundaryConditions:
    def build(self, group, loop, layouts):
        sources = {p.region.owner: p for p in group.footprints}
        # Physical boundary owners need not match the logical junction count.
        # Consumed internal profiles impose no phantom boundary height seam.
        # Only profiles exposed on the physical boundary constrain its datum
        # and crown ceiling. A caller may supply a larger layout dictionary;
        # fully consumed internal approaches also survive in raw footprints.
        exposed = set()
        for boundary in group.loops:
            for span in boundary.spans:
                label = rim_label(span, sources)
                if label.kind == 'mouth':
                    exposed.add((span.owner, label.index))
        source_profiles = [a.profile for owner in sources for a in layouts[owner].approaches]
        profiles = [a.profile for owner in sources for index, a in enumerate(layouts[owner].approaches)
                    if (owner, index) in exposed]
        if not profiles:
            raise ValueError('Compound boundary has no exposed incoming profiles; another boundary policy is required')
        datum, ceiling = profiles[0].datum, max(p.datum + p.rise for p in profiles)
        # Keep the existing common-minimum input invariant for this compound,
        # including consumed approaches. Only the crown ceiling is projected
        # onto exposed profiles; unrelated dictionary entries are excluded.
        if any(abs(p.datum - datum) > 1e-7 for p in source_profiles):
            raise ValueError('Compound boundary requires a common minimum height')
        tokens, parameters = [], []
        for span in loop.spans:
            label = rim_label(span, sources)
            tokens.append((span.owner, label.index) if label.kind == 'mouth' else None)
            parameters.append((label.first, label.last))
        count = len(loop.xy)
        rules, ports = [None] * count, []
        starts = [i for i in range(count) if tokens[i] is not None and tokens[i] != tokens[i - 1]]
        if not starts:
            if not loop.hole:
                raise ValueError('A closed exterior without incoming ports needs a separate boundary policy')
            rules = [PlateauEdge(ceiling)] * count
        else:
            for start in starts:
                owner, side = tokens[start]
                approach = layouts[owner].approaches[side]
                distance = layouts[owner].outward_stations[side]
                indices, values, index = [start], [parameters[start][0]], start
                while tokens[index] == (owner, side):
                    u, v = parameters[index]
                    if abs(u - values[-1]) > 1e-8 or v >= u:
                        raise ValueError('An incoming compound profile has a broken boundary partition')
                    rules[index] = ProfileEdge(approach, distance, float(u), float(v))
                    values.append(float(v))
                    index = (index + 1) % count
                    indices.append(index)
                    if index == start:
                        raise ValueError('A closed boundary cannot be one incoming profile')
                if abs(values[0] - 1.) > 1e-8 or abs(values[-1]) > 1e-8:
                    raise ValueError('Compound clipping retained only part of an incoming profile; extend its physical port')
                from .source_mouths import source_mouths
                aliases = source_mouths(approach, distance)
                first_alias = aliases[0]
                sid, end = first_alias.endpoint.section_id, first_alias.endpoint.end
                ports.append(Port(owner, side, sid, end, first_alias.station, tuple(indices), tuple(values), aliases))
            for current, next_port in zip(ports, (*ports[1:], ports[0])):
                first, last = current.vertices[-1], next_port.vertices[0]
                indices, index = [first], first
                while index != last:
                    if tokens[index] is not None:
                        raise ValueError('A curb chain crosses an unexpected incoming profile')
                    index = (index + 1) % count
                    indices.append(index)
                    if len(indices) > count:
                        raise ValueError('Compound curb chain does not close')
                if len(indices) < 2:
                    raise ValueError('Adjacent compound profiles need an explicit shared corner')
                xy = loop.xy[indices]
                vectors = np.diff(xy, axis=0)
                lengths = np.linalg.norm(vectors, axis=1)
                stations = np.r_[0., np.cumsum(lengths)]
                tangent = vectors / lengths[:, None]
                a = layouts[current.owner].approaches[current.approach_index]
                b = layouts[next_port.owner].approaches[next_port.approach_index]
                za, ga = ProfileEdge(a, layouts[current.owner].outward_stations[current.approach_index], 0., 0.).sample(0.)
                zb, gb = ProfileEdge(b, layouts[next_port.owner].outward_stations[next_port.approach_index], 1., 1.).sample(0.)
                ka, kb = float(ga @ -right(tangent[0])), float(gb @ -right(tangent[-1]))
                if min(ka, kb) < -1e-7:
                    raise ValueError('Compound curb does not face the incoming profile interior')
                for j, edge in enumerate(indices[:-1]):
                    rules[edge] = CurbEdge(float(za), float(zb), max(0., ka), max(0., kb), float(stations[-1]),
                        float(stations[j]), float(stations[j + 1]), tuple(tangent[j]))
        if any(rule is None for rule in rules):
            raise ValueError('Compound boundary has an unassigned height condition')
        heights, gradients = zip(*(rule.sample(0.) for rule in rules))
        heights, gradients = np.asarray(heights), np.asarray(gradients)
        # Port values and physical tangent planes remain authoritative, including
        # the final vertex where the outgoing edge belongs to the curb chain.
        for port in ports:
            approach = layouts[port.owner].approaches[port.approach_index]
            distance = layouts[port.owner].outward_stations[port.approach_index]
            q = 2. * np.asarray(port.outward_u) - 1.
            heights[list(port.vertices)] = approach.profile.datum + approach.profile.value(q)
            gradients[list(port.vertices)] = _profile_gradient(approach, distance, q)
        if not np.all(np.isfinite(gradients)) or np.min(heights) < datum - 1e-8 or np.max(heights) > ceiling + 1e-8:
            raise ValueError('Compound boundary conditions exceed their physical envelope')
        heights.flags.writeable = gradients.flags.writeable = False
        return BoundaryData(loop.xy, heights, gradients, tuple(rules), tuple(ports), datum, ceiling, loop.hole)
