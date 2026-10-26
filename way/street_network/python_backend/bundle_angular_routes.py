"""Source-preserving metric routes for angular Bundle boundary construction.

The route concatenates already prepared Section paths, never their controls or
newly fitted cubics. Explicit two-port incidences establish source membership.
Compatible joins retain the complete cross-section; dimensional joins use the
original finite transition provider. Bent area-free Exact joins still need
their common miter construction.
"""
import numpy as np
from .bundle_axis import AxialChartUnavailable
from .model import End, Key
from .routes import Piece, Route


class AngularMemberRoute(Route):
    def __init__(self, context, member, role):
        pieces = tuple(Piece(context[Key('path', visit.section_id)], visit.forward)
                       for visit in member.path)
        sites = {end: site for site in (*context.network.intersections, *context.network.connections)
                 for end in site.ends}
        profiles = tuple(context[Key('profile', piece.section.id)].oriented(0 if piece.forward else 1)
                         for piece in pieces)
        first = pieces[0]
        self.transitions = {}
        self.maximum_reach = max(abs(piece.section.offset)+.5*piece.section.width for piece in pieces)
        self.offset = first.section.offset * (1 if first.forward else -1)
        self.width = first.section.width
        self.profile = profiles[0]
        # Blender stores source coordinates as float32. A collinear join can
        # consequently acquire different endpoint normals on its two pieces.
        # Reconcile only a bounded, at-most-ten-micrometre curb discrepancy; this is
        # not a license to bridge a real bent or dimensional Section seam.
        self.join_tolerance = min(1e-5, .01*context.network.settings.curve_tolerance)
        self.maximum_curb_discrepancy = 0.
        for i, (piece, profile) in enumerate(zip(pieces, profiles)):
            offset = piece.section.offset * (1 if piece.forward else -1)
            if (abs(profile.rise-self.profile.rise)>1e-9
                    or abs(profile.datum-self.profile.datum)>1e-9
                    or profile.half != (1 if role==0 else -1)):
                raise AxialChartUnavailable('Angular Section chains need compatible full cross-sections')
            if not i:
                continue
            previous = pieces[i-1]
            exit = End(previous.section.id, int(previous.forward))
            entry = End(piece.section.id, int(not piece.forward))
            site = sites.get(exit)
            if site is None or entry not in site.ends:
                raise ValueError('Angular Bundle chains need explicit two-ended Section seams')
            if len(site.ends)!=2:
                raise AxialChartUnavailable('Angular member crosses an internal physical junction')
            a, ta, na = previous.frame(previous.path.length)
            b, tb, nb = piece.frame(0.)
            if np.linalg.norm(a-b)>1e-8:
                raise ValueError('Angular Bundle Section endpoint copies do not meet')
            previous_offset = previous.section.offset*(1 if previous.forward else -1)
            if abs(piece.section.width-previous.section.width)>1e-8 or abs(offset-previous_offset)>1e-8:
                transition = context[site.layout_key].transition
                if transition is None:
                    raise AxialChartUnavailable('Angular dimensional seam needs a finite source-owned transition')
                self.transitions[site.key] = transition
                continue
            discrepancy = max(float(np.linalg.norm((a+d*na)-(b+d*nb)))
                              for d in (offset-.5*piece.section.width,offset+.5*piece.section.width))
            self.maximum_curb_discrepancy = max(self.maximum_curb_discrepancy,discrepancy)
            if np.linalg.norm(ta-tb)>1e-5 or discrepancy>self.join_tolerance:
                raise AxialChartUnavailable('Angular Section join needs its common miter or transition provider')
        super().__init__(pieces, terminal_free=True)
        self.points = np.concatenate([piece.path.points if piece.forward else piece.path.points[::-1]
                                      for piece in pieces[:1]] + [
            (piece.path.points if piece.forward else piece.path.points[::-1])[1:] for piece in pieces[1:]])
        self.exact = all(piece.path.exact for piece in pieces)
        self.smooth = all(not piece.path.exact for piece in pieces)

    def certify(self, frame):
        """Every original polynomial/offset gets a whole-source certificate."""
        return min(frame.certify(piece.path, forward=piece.forward) for piece in self.pieces)

    def curbs(self, settings):
        result = []
        for side in (-1., 1.):
            parts, total = [], 0
            for piece in self.pieces:
                path, source = piece.path, piece.section
                stations = path.samples(0., path.length, settings.spacing,
                                        tolerance=.1*settings.curve_tolerance)
                q = side if piece.forward else -side
                xy = path.point(stations, source.offset+q*.5*source.width)
                if not piece.forward:
                    xy = xy[::-1]
                if parts:
                    if np.linalg.norm(parts[-1][-1]-xy[0])>self.join_tolerance:
                        raise ValueError('Angular route physical curbs do not meet at their explicit seam')
                    xy = xy[1:]
                total += len(xy)
                if total>200000:
                    raise ValueError('Angular Section route exceeds its boundary sampling budget')
                parts.append(xy)
            result.append(np.concatenate(parts))
        return tuple(result)

    def source_spans(self):
        return [dict(section=piece.section.id, forward=piece.forward,
                     first=float(self.stations[i]), last=float(self.stations[i+1]))
                for i, piece in enumerate(self.pieces)]
