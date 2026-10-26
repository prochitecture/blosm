"""Identity-preserving intersection approaches across compatible Section seams.

The route is derived geometry, not a replacement Street/Section or an edited
polyline. Every trim is mapped back to its original Section's directed end.
Only explicit two-ended incidence can continue an approach; proximity cannot.
"""
from dataclasses import dataclass
import numpy as np
from .model import End, Key, Traversal


def route_key(end):
    return Key('approach', 2 * end.section_id + end.end)


@dataclass(frozen=True)
class Candidate:
    traversals: tuple
    seams: tuple
    terminal_site: object = None
    terminal_selection: bool = False


class RouteTopology:
    def __init__(self, network):
        self.selection_ends = {End(s.id, end) for s in network.sections for end in s.selection_ends}
        self.seams = {}
        self.incidents = {end: site.key for site in (*network.intersections, *network.connections)
                          for end in site.ends}
        for site in (*network.intersections, *network.connections):
            if len(site.ends) == 2:
                for end, partner in (site.ends, site.ends[::-1]):
                    self.seams[end] = (partner, site.key)

    def trace(self, end):
        visits, seams, seen = [], [], set()
        while end.section_id not in seen:
            seen.add(end.section_id)
            visits.append(Traversal(end.section_id, end.end == 0))
            exit = End(end.section_id, 1 - end.end)
            next_end = self.seams.get(exit)
            if next_end is None or next_end[0].section_id in seen:
                break
            end, seam = next_end
            seams.append(seam)
        return Candidate(tuple(visits), tuple(seams), self.incidents.get(exit), exit in self.selection_ends)

    def attachments(self, network):
        """Potential dependency closure, including mouths beyond a seam."""
        result = {s.id: {} for s in network.sections}
        for site in (*network.intersections, *network.connections):
            for end in site.ends:
                visits = self.trace(end).traversals if len(site.ends) > 2 else (Traversal(end.section_id, end.end == 0),)
                for visit in visits:
                    at = 0 if visit.forward else 1
                    result[visit.section_id][(site.key, at)] = (site, at)
        return {sid: tuple(value[k] for k in sorted(value)) for sid, value in result.items()}


@dataclass(frozen=True)
class Piece:
    path: object
    forward: bool

    @property
    def section(self):
        return self.path.section

    def frame(self, station):
        p, t, n = self.path.frame(station if self.forward else self.path.length - station)
        return (p, t, n) if self.forward else (p, -t, -n)


class Route:
    outward = True

    def __init__(self, pieces, *, terminal_free=False):
        self.pieces = tuple(pieces)
        self.section = self.pieces[0].section
        self.stations = np.r_[0., np.cumsum([p.path.length for p in self.pieces])]
        self.length = float(self.stations[-1])
        self.terminal_free = terminal_free

    @classmethod
    def build(cls, context, candidate):
        pieces = []
        first_profile = None
        first_offset = None
        for visit in candidate.traversals:
            path = context[Key('path', visit.section_id)]
            profile = context[Key('profile', visit.section_id)].oriented(0 if visit.forward else 1)
            piece = Piece(path, visit.forward)
            offset = path.section.offset * (1 if visit.forward else -1)
            if pieces:
                previous = pieces[-1]
                if (abs(path.section.width - previous.section.width) > 1e-8
                        or abs(offset - first_offset) > 1e-8
                        or abs(profile.rise - first_profile.rise) > 1e-9
                        or abs(profile.datum - first_profile.datum) > 1e-9
                        or profile.half != first_profile.half):
                    break
                a, ta, na = previous.frame(previous.path.length)
                b, tb, nb = piece.frame(0.)
                # Preserve the complete source cross-section, including Exact
                # endpoint normals. A bent Exact seam needs its miter provider.
                reach = abs(offset) + path.section.width * .5
                if (np.linalg.norm(a - b) > 1e-8 or np.linalg.norm(ta - tb) > 1e-8
                        or reach * np.linalg.norm(na - nb) > 1e-8):
                    break
            pieces.append(piece)
            first_profile, first_offset = profile, offset
        return cls(pieces, terminal_free=(len(pieces) == len(candidate.traversals)
                                         and candidate.terminal_site is None and not candidate.terminal_selection))

    def frame(self, station):
        s = np.asarray(station)
        if np.any(s < -1e-8) or np.any(s > self.length + 1e-8):
            raise ValueError('Station falls outside the explicit approach route')
        s = np.clip(s, 0., self.length)
        if len(self.pieces) == 1:
            return self.pieces[0].frame(s)
        at = np.clip(np.searchsorted(self.stations, s, side='right') - 1, 0, len(self.pieces) - 1)
        point, tangent, normal = (np.empty(s.shape + (2,)) for _ in range(3))
        for index, piece in enumerate(self.pieces):
            mask = at == index
            if np.any(mask):
                p, t, n = piece.frame(np.clip(s[mask] - self.stations[index], 0., piece.path.length))
                point[mask], tangent[mask], normal[mask] = p, t, n
        return point, tangent, normal

    def trims(self, distance):
        """Map an outward mouth to all touched source Section end-stations."""
        result = {}
        for i, piece in enumerate(self.pieces):
            local = float(np.clip(distance - self.stations[i], 0., piece.path.length))
            end = 0 if piece.forward else 1
            result[(piece.section.id, end)] = local if piece.forward else piece.path.length - local
            if distance <= self.stations[i + 1] + 1e-9:
                break
        return result
