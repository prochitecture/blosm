"""Bounded four-sided quad atlases for elongated, indented chart blocks.

Atlas corners are selected on existing rim vertices. All requested boundary
samples remain authoritative; the cached pattern terminates unequal boundary
partitions internally. Different local shape energies are construction choices,
not changes to source dimensions, heights or physical curb curves.
"""
import itertools
import math
import numpy as np
from .mesh import Mesh, cross
from .model import Key
from .quad_charts import ChartResult, FixedRimQuality
from .quad_rotations import QuadEdgeRotations
from .regions import Region, RegionUnion

_CORNER_CHUNK = 1024


def quality(mesh):
    result = mesh.validate()
    points = mesh.xyz[mesh.quads, :2]
    edge = np.roll(points, -1, axis=1) - points
    lengths = np.linalg.norm(edge, axis=-1)
    result['minimum_corner_sine'] = float(np.min(cross(edge, np.roll(edge, -1, axis=1))
        / (lengths * np.roll(lengths, -1, axis=1))))
    return result


class FourSidedQuadChart:
    def __init__(self, patterns, *, max_cells=50000, max_candidates=24,
                 quality_iterations=100, maximum_ratio=8.):
        if (type(max_cells) is not int or max_cells < 1
                or type(max_candidates) is not int or max_candidates < 1
                or type(quality_iterations) is not int or quality_iterations < 0
                or not math.isfinite(maximum_ratio) or maximum_ratio < 1.):
            raise ValueError('Invalid four-sided quad-chart budgets')
        self.patterns, self.max_cells, self.max_candidates = patterns, max_cells, max_candidates
        self.quality_iterations, self.maximum_ratio = quality_iterations, maximum_ratio

    @staticmethod
    def corners(boundary):
        n = len(boundary)
        if n == 4:
            return [(0., (0, 1, 2, 3))]
        edge = np.roll(boundary, -1, axis=0) - boundary
        lengths = np.linalg.norm(edge, axis=1)
        unit = edge / lengths[:, None]
        previous = np.roll(unit, 1, axis=0)
        turns = np.arctan2(cross(previous, unit), np.sum(previous * unit, axis=1))
        chosen = sorted(np.flatnonzero(turns > .05),
                        key=lambda i: (-round(float(turns[i]), 8), int(i)))[:16]
        choices = sorted({(int(i) + shift) % n for i in chosen for shift in (-1, 0, 1)})
        if len(choices) > 36:
            choices = sorted(sorted(choices, key=lambda i: (-round(float(turns[i]), 8), i))[:36])
        # The scalar corner angles are invariant across the candidate sets.
        # Keep math.atan2 and the existing operation order for exact scores.
        angles = np.zeros(n)
        for i in choices:
            angles[i] = math.atan2(float(-cross(-unit[i - 1], unit[i])),
                                  float(-unit[i - 1] @ unit[i])) % math.tau
        options = []
        # At most C(36,4)=58,905 preliminary sets; only max_candidates atlases
        # are allocated/mapped. The nearby choices permit opposite-side parity
        # without moving or inserting a physical boundary point.
        candidates = itertools.combinations(choices, 4)
        while True:
            block = list(itertools.islice(candidates, _CORNER_CHUNK))
            if not block:
                break
            corners = np.asarray(block, dtype=np.intp)
            counts = np.diff(np.column_stack((corners, corners[:, 0] + n)), axis=1)
            keep = ((counts.min(axis=1) >= 2) & ((counts[:, 0] - counts[:, 2]) % 2 == 0)
                    & ((counts[:, 1] - counts[:, 3]) % 2 == 0))
            corners = corners[keep]
            if not len(corners):
                continue
            points = boundary[corners]
            chord = np.roll(points, -1, axis=1) - points
            size = np.linalg.norm(chord, axis=2)
            sine = cross(chord, np.roll(chord, -1, axis=1)) / (size * np.roll(size, -1, axis=1))
            at = angles[corners]
            keep = (sine.min(axis=1) >= .15) & (at.min(axis=1) >= .20) & (at.max(axis=1) <= math.pi + .01)
            corners, size, sine, at = corners[keep], size[keep], sine[keep], at[keep]
            balance = np.sum(np.log(size[:, :2] / size[:, 2:]) ** 2, axis=1)
            score = (np.sum((at - math.pi / 2.) ** 2, axis=1)
                     + 1.5 * np.sum(1. - sine ** 2, axis=1) + .3 * balance)
            options.extend((round(float(value), 8), tuple(map(int, row))) for value, row in zip(score, corners))
        return sorted(options)

    def map(self, boundary, corners):
        n = len(boundary)
        order = np.r_[np.arange(corners[0], n), np.arange(corners[0])]
        origin = boundary[0]
        rim = boundary[order] - origin
        b, c, d = (i - corners[0] for i in corners[1:])
        sides = (rim[:b + 1], rim[b:c + 1], rim[c:d + 1][::-1],
                 np.vstack((rim[d:], rim[:1]))[::-1])
        counts = tuple(len(side) - 1 for side in sides)
        cells = (min(counts[0], counts[2]) * min(counts[1], counts[3])
                 + 3 * (abs(counts[0] - counts[2]) + abs(counts[1] - counts[3])) // 2)
        if cells > self.max_cells:
            raise ValueError('Four-sided atlas exceeds its cell budget')
        pattern = self.patterns.quad(*counts)
        u, v = pattern.uv.T
        samples, side_ids = [], []
        for side, points in enumerate(sides):
            ids = self.patterns.boundary(pattern, side)
            side_ids.append(ids)
            parameters = pattern.uv[ids, 0 if side in (0, 2) else 1]
            query = u if side in (0, 2) else v
            samples.append(np.column_stack([np.interp(query, parameters, points[:, axis]) for axis in range(2)]))
        bottom, right, top, left = samples
        a, bp, cp, dp = rim[[0, b, c, d]]
        bilinear = ((1. - u) * (1. - v))[:, None] * a + (u * (1. - v))[:, None] * bp
        bilinear += (u * v)[:, None] * cp + ((1. - u) * v)[:, None] * dp
        xy = (1. - v[:, None]) * bottom + v[:, None] * top + (1. - u[:, None]) * left + u[:, None] * right - bilinear
        for ids, points in zip(side_ids, sides):
            xy[ids] = points
        trace = np.r_[side_ids[0][:-1], side_ids[1][:-1], side_ids[2][::-1][:-1], side_ids[3][::-1][:-1]]
        if len(trace) != n or len(np.unique(trace)) != n:
            raise ValueError('Four-sided atlas lost its ordered physical rim')
        internal = np.flatnonzero(~np.isin(np.arange(len(xy)), trace))
        remap = np.empty(len(xy), dtype=np.int32)
        remap[trace], remap[internal] = order, np.arange(n, len(xy))
        points = np.zeros((len(xy), 3))
        points[remap, :2] = xy + origin
        points[:n, :2] = boundary
        return Mesh(points, remap[pattern.quads], seams={0: np.arange(n)}), counts

    def accepted(self, report):
        return report['maximum_edge_ratio'] <= self.maximum_ratio and report['minimum_corner_sine'] >= .03

    def build(self, boundary):
        boundary = np.asarray(boundary, dtype=float)
        if boundary.ndim != 2 or boundary.shape[1] != 2 or len(boundary) < 4 or len(boundary) % 2:
            raise ValueError('Four-sided chart needs an even XY rim')
        if len(boundary) // 2 - 1 > self.max_cells:
            raise ValueError('Four-sided rim exceeds its cell budget')
        RegionUnion().build((Region(Key('chart', 0), boundary),))
        options, attempted = self.corners(boundary), 0
        for score, corners in options[:self.max_candidates]:
            attempted += 1
            try:
                initial, counts = self.map(boundary, corners)
                report = quality(initial)
            except ValueError:
                continue
            variants = [(initial, report, dict(power=0, steps=0))]
            if not self.accepted(report) and report['maximum_edge_ratio'] <= 2. * self.maximum_ratio:
                for power in (2, 1):
                    candidate, improvement = FixedRimQuality(iterations=self.quality_iterations, power=power).improve(initial)
                    candidate_report = quality(candidate)
                    variants.append((candidate, candidate_report, dict(power=power, **improvement)))
                    if self.accepted(candidate_report):
                        break
                if not any(self.accepted(item[1]) for item in variants):
                    candidate, candidate_report, previous = min(variants, key=lambda item: item[1]['maximum_edge_ratio'])
                    candidate, rotations = QuadEdgeRotations(maximum_ratio=self.maximum_ratio,
                        max_cells=self.max_cells).improve(candidate, FixedRimQuality(iterations=self.quality_iterations, power=2))
                    variants.append((candidate, quality(candidate), dict(previous=previous, edge_rotations=rotations)))
            accepted = [item for item in variants if self.accepted(item[1])]
            if accepted:
                mesh, report, improvement = accepted[0]
                return ChartResult(mesh, dict(quality=report, initial=quality(initial), improvement=improvement,
                    atlas_corners=corners, rim_intervals=counts, candidates=attempted,
                    construction='four-sided direct quad atlas with fixed rim'))
        raise ValueError('Four-sided quad chart needs another compatible corner atlas')
