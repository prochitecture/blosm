"""Complete evaluated physical mouths, including non-incidental source aliases."""
from dataclasses import dataclass
import numpy as np
from .source_mouths import source_mouths


@dataclass(frozen=True)
class PreparedMergedMouth:
    owner: object
    approach_index: int
    approach: object
    station: float
    xyz: np.ndarray
    rim_references: tuple
    crest_column: int

    @classmethod
    def from_product(cls, product, approach_index):
        layout = product.group.physical_layout
        if layout is None or product.physical_rims is None:
            raise ValueError('A merged mouth needs its completed prepared physical product')
        if not 0 <= approach_index < len(layout.approaches):
            raise ValueError('A merged mouth needs an actual physical approach index')
        approach = layout.approaches[approach_index]
        station = layout.outward_stations[approach_index]
        aliases = source_mouths(approach, station)
        if len(aliases) != 2 or not hasattr(approach, 'end_index'):
            raise ValueError('A merged body interface needs two explicitly prepared source aliases')

        def sample(alias, u, *, peak=False):
            sid, end = alias.endpoint.section_id, alias.endpoint.end
            if alias.original_incidence:
                table = product.peaks if peak else product.port_points
                key = (sid, end) if peak else (sid, end, round(float(alias.outward_parameter(u)), 12))
            else:
                table = product.derived_peaks if peak else product.derived_port_points
                key = (product.owner, sid, end) if peak else (
                    product.owner, sid, end, round(float(alias.outward_parameter(u)), 12))
            if not table or key not in table:
                raise ValueError('A merged interface lost its explicitly owned original or derived source sample')
            loop, vertex, original_station = table[key]
            if abs(original_station-alias.station) > 1e-8:
                raise ValueError('A merged interface changed its original source station')
            return loop, vertex

        endpoints = []
        for u in (0., 1.):
            candidates = [alias for alias in aliases if alias.contains(u)]
            if len(candidates) != 1:
                raise ValueError('A merged outer curb needs one authoritative source alias')
            endpoints.append(sample(candidates[0], u))
        peaks = [sample(alias, alias.physical_peak_u, peak=True) for alias in aliases]
        if peaks[0] != peaks[1] or endpoints[0][0] != endpoints[1][0]:
            raise ValueError('A merged interface needs one shared crest on one physical boundary loop')
        loop, first = endpoints[0]; last = endpoints[1][1]
        rim, kinds = product.physical_rims[loop], product.kinds[loop]
        candidates = []
        for direction in (1, -1):
            indices, at = [first], first
            while at != last:
                edge = at if direction == 1 else (at-1) % len(rim)
                if kinds[edge] != 'mouth':break
                at = (at+direction) % len(rim);indices.append(at)
                if len(indices) > len(rim):break
            if at == last and peaks[0][1] in indices:candidates.append(indices)
        if len(candidates) != 1:
            raise ValueError('A merged interface must be one complete actual mouth chain without a curb')
        indices = candidates[0];xyz = rim[indices].copy()
        delta = xyz[-1, :2]-xyz[0, :2]
        squared = float(delta@delta)
        if squared <= 1e-12:raise ValueError('A merged physical mouth is collapsed')
        u = (xyz[:, :2]-xyz[0, :2])@delta/squared
        expected = approach.point(station, approach.offset+approach.width*(u-.5))
        if (np.any(np.diff(u) <= 1e-10) or np.max(np.linalg.norm(expected-xyz[:, :2], axis=1)) > 1e-7
                or np.max(abs(xyz[:, 2]-approach.field.values(xyz[:, :2]))) > 1e-8):
            raise ValueError('An evaluated mouth disagrees with its prepared whole physical profile')
        xyz.flags.writeable = False
        return cls(product.owner, approach_index, approach, float(station), xyz,
                   tuple((loop, i) for i in indices), indices.index(peaks[0][1]))

    def split_body_cut(self, points, *, crest_column=None):
        """Return outer-to-crest halves in the declared global member order."""
        points = np.asarray(points)
        crest = self.crest_column if crest_column is None else crest_column
        if (points.shape != self.xyz.shape or type(crest) is not int or not 0 < crest < len(points)-1
                or not np.isfinite(points).all()):
            raise ValueError('A transported body cut lost its complete whole-profile partition')
        first, last = points[:crest+1].copy(), points[crest:][::-1].copy()
        return (first, last) if self.approach.end_index == 0 else (last, first)
