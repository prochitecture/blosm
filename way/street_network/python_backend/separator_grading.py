"""Grade a compact run of short physical bank edges into a coarser grass disk."""
import numpy as np

from .fixed_rim_bands import BandedSeparatorChart, _qualified
from .hole_banded_chart import audit
from .mesh import Mesh, publication_coordinates
from .model import Key
from .quad_charts import ChartResult
from .validation import validate_overlaps


class LocalSeparatorBands(BandedSeparatorChart):
    """Nested quad bands share both endpoints of the dense boundary run."""

    def build(self, boundary):
        xy = np.asarray(boundary, dtype=float)
        count = len(xy)
        if xy.shape != (count, 2) or not 8 <= count <= 8192 or count % 2 or not np.isfinite(xy).all():
            raise ValueError('Local separator bands need a finite even rim')
        lengths = np.linalg.norm(np.roll(xy, -1, axis=0) - xy, axis=1)
        short = lengths < np.median(lengths) / 32
        starts = np.flatnonzero(short & ~np.roll(short, 1))
        size = int(np.sum(short))
        if len(starts) != 1 or not 2 <= size <= min(32, count - 4):
            raise ValueError('Local separator bands need one compact run of short edges')
        run = (int(starts[0]) + np.arange(size + 1)) % count
        left, right = int(run[0]), int(run[-1])
        before, after = (left - 1) % count, (right + 1) % count
        local = xy[np.r_[before, run, after]]
        chord = xy[right] - xy[left]
        half = .5 * np.linalg.norm(chord)
        if half <= 1e-8 or min(lengths[before], lengths[right]) < 8 * max(lengths[run[:-1]]):
            raise ValueError('Separator bank has no isolated spacing transition')
        tangent = chord / (2 * half)
        inward = np.array((-tangent[1], tangent[0]))
        center = (xy[right] + xy[left]) / 2
        if np.any(np.diff(local @ tangent) <= 0.) or np.max(abs((xy[run]-center) @ inward)) > .1 * half:
            raise ValueError('The fine bank run is not locally straight')
        endpoints = xy[[before, after]] - center
        radii = np.linalg.norm(endpoints, axis=1)
        angles = np.arctan2(endpoints @ inward, endpoints @ tangent)
        angles[0] %= 2 * np.pi
        # Try a small family of radial growth rates. All physical bank vertices,
        # including the two long adjacent edges, stay fixed in every candidate.
        levels = int(np.clip(round(np.log2(max(radii) / half)), 3, 8))
        for layers in dict.fromkeys((levels, min(8, levels+1), max(3, levels-1))):
            for inner_radius, aspect in ((1.4,.45),(1.3,.5),(1.2,.55),(1.15,.65),(1.,.85)):
                patch, interface = self._patch(local, center, tangent, inward, radii, angles,
                                               half * inner_radius, layers, aspect)
                if not _qualified(patch):
                    continue
                try:
                    return self._fill(xy, before, after, run, patch, interface)
                except ValueError:
                    continue
        raise ValueError('No local separator bands met the fixed-rim quality and coverage checks')

    @staticmethod
    def _patch(local, center, tangent, inward, radii, angles, inner_radius, layers, aspect):
        points = list(local)
        old = list(range(1, len(local)-1))
        faces = []
        for j, f in enumerate(np.linspace(0., 1., layers)):
            axes = radii**f * inner_radius**(1.-f)
            margin = 1.05 * (1.-f)
            theta = np.linspace(angles[0]-margin, angles[1]+margin, len(old))
            scale = np.linspace(*axes, len(old))[:, None]
            arc = center + scale*(np.cos(theta)[:, None]*tangent +
                                  np.sin(theta)[:, None]*inward*aspect**(1.-f))
            ids = []
            for i, p in enumerate(arc):
                if j == layers-1 and i in (0, len(arc)-1):
                    ids.append(0 if i == 0 else len(local)-1)
                else:
                    ids.append(len(points))
                    points.append(p)
            row = [old[0], *ids, old[-1]]
            # An extra interval at each end lets each band meet the fixed
            # endpoint with a full quad while its radial edges grow longer.
            faces.append((old[0], old[1], row[2], row[1]))
            faces.extend((old[i], old[i+1], row[i+2], row[i+1]) for i in range(1, len(old)-2))
            faces.append((old[-2], old[-1], row[-2], row[-3]))
            old = row
        xyz = np.column_stack((points, np.zeros(len(points))))
        return Mesh(xyz, np.asarray(faces)), np.asarray(old[1:-1])

    def _fill(self, xy, before, after, run, patch, interface):
        count = len(xy)
        sources = np.r_[before, run, after]
        added = patch.xyz[len(sources):, :2]
        remap = np.r_[sources, count + np.arange(len(added))]
        other = (after + 1 + np.arange((before-after-1) % count)) % count
        outer = np.r_[remap[interface], other]
        outer = np.roll(outer, -np.argmin(outer))
        points = np.vstack((xy, added))
        return _complete_disk(self, xy, points, remap[patch.quads], outer,
            construction='graded bands around a short boundary run', short_intervals=len(run)-1)


def _complete_disk(chart, xy, points, faces, outer, **metadata):
    """Weld the core by index and certify the entire disk against its original rim."""
    core = BandedSeparatorChart(chart.patterns, spacing=chart.spacing,
        interior_spacing=chart.interior_spacing, max_cells=chart.max_cells-len(faces),
        bands=chart.bands).build(points[outer]).mesh
    core_indices = np.r_[outer, len(points) + np.arange(len(core.xyz)-len(outer))]
    points = np.vstack((points, core.xyz[len(outer):, :2]))
    mesh = Mesh(np.column_stack((points, np.zeros(len(points)))),
        np.vstack((faces, core_indices[core.quads])), 'grass')
    mesh.seams = {0: np.arange(len(xy))}
    certificate = audit(mesh, xy)
    for candidate in (mesh, mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True))):
        validate_overlaps({Key('grass', 0): candidate})
    return ChartResult(mesh, dict(exact_whole_boundary_edges=True, rim_certificate=certificate,
        local_quads=len(faces), **metadata))


class PointedSeparatorChart(BandedSeparatorChart):
    """Stagger interior rows down a narrow tip between two fixed carriageway banks."""

    def build(self, boundary):
        xy = np.asarray(boundary, dtype=float)
        count = len(xy)
        if xy.shape != (count, 2) or not 16 <= count <= 8192 or count % 2 or not np.isfinite(xy).all():
            raise ValueError('Pointed separator needs a finite even rim')
        edges = np.roll(xy, -1, axis=0) - xy
        lengths = np.linalg.norm(edges, axis=1)
        if lengths.min() <= 1e-8:
            raise ValueError('Degenerate separator edge')
        unit = edges / lengths[:, None]
        prev = np.roll(unit, 1, axis=0)
        turn = prev[:, 0]*unit[:, 1] - prev[:, 1]*unit[:, 0]
        tips = np.flatnonzero((turn > 0.) & ((prev*unit).sum(axis=1) < -np.cos(np.deg2rad(10.))))
        if len(tips) != 1:
            raise ValueError('Pointed separator needs one acute convex tip')
        tip = int(tips[0])
        for rows in (12, 8, 16, 24, 5):
            if 2*(rows+2) >= count:
                continue
            indices = np.arange(rows+3)
            left, right = (tip-indices) % count, (tip+indices) % count
            axis = unit[tip] - prev[tip]
            axis /= np.linalg.norm(axis)
            banks = (xy[left], xy[right])
            depths = [(bank - xy[tip]) @ axis for bank in banks]
            if any(np.any(np.diff(t) <= 0.) for t in depths):
                continue
            # Two short intervals at each bank end lead into the ordinary spacing.
            if any(max(np.diff(t[:3])) >= .5*(t[3]-t[2]) for t in depths):
                continue
            for bias, stagger in ((.16,.45), (.18,.55), (.19,.7)):
                try:
                    points, faces, interface = self._tip(rows, left, right, xy, banks, depths, bias, stagger)
                    patch = Mesh(np.column_stack((points, np.zeros(len(points)))), faces)
                    if not _qualified(patch):
                        continue
                    other = (right[rows]+1+np.arange(count-2*rows-1)) % count
                    outer = np.r_[interface, other]
                    outer = np.roll(outer, -np.argmin(outer))
                    return _complete_disk(self, xy, points, faces, outer,
                        construction='staggered quad rows at an acute separator tip', tip_rows=rows)
                except ValueError:
                    continue
        raise ValueError('No bounded pointed separator chart passed its original rim checks')

    @staticmethod
    def _tip(rows, left, right, xy, banks, depths, bias, stagger):
        points = list(xy)
        a, b = depths

        def point(depth, fraction):
            if depth > min(a[-1], b[-1]):
                raise ValueError('Interior row leaves the common bank span')
            sides = [np.array([np.interp(depth, t, coord) for coord in bank.T])
                     for t, bank in zip(depths, banks)]
            points.append((1.-fraction)*sides[0] + fraction*sides[1])
            return len(points)-1

        step = min(a[3]-a[2], b[3]-b[2])
        middle = point(max(a[1], b[1]) + 2*min(a[1], b[1]), .5)
        near_left = point(a[2]+bias*step, .22)
        near_right = point(b[2]+bias*step, .78)
        end = point(max(a[2], b[2])+2*bias*step, .5)
        faces = [(left[1],left[0],right[1],middle), (left[2],left[1],middle,near_left),
                 (middle,right[1],right[2],near_right), (near_left,middle,near_right,end)]
        previous = [left[2], near_left, end, near_right, right[2]]
        for j in range(3, rows+1):
            step = min(a[j]-a[j-1], b[j]-b[j-1])
            row = [left[j]]
            for fraction in (.25, .5, .75):
                depth = (1.-fraction)*a[j] + fraction*b[j]
                row.append(point(depth + 2*min(fraction,1.-fraction)*stagger*step, fraction))
            row.append(right[j])
            faces.extend((row[i],previous[i],previous[i+1],row[i+1]) for i in range(4))
            previous = row
        return np.asarray(points), np.asarray(faces), previous
