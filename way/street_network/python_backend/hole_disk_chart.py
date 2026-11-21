"""Fixed-rim direct quad atlases for narrow separator disks.

Grass keeps its established contract: positive corners and edge ratio <=8,
including float32 publication. A sharply tapered physical grass rim may force
small angles; their sines remain explicit metadata. Pavement's separate 0.03
sine criterion and the ordinary four-sided atlas policy are unchanged.
"""
import numpy as np
from .four_side_charts import FourSidedQuadChart, quality
from .mesh import publication_coordinates
from .model import Key
from .quad_charts import ChartResult, FixedRimQuality
from .regions import Region, RegionUnion


class _NarrowAtlas(FourSidedQuadChart):
    def accepted(self, report):
        return report['maximum_edge_ratio'] <= self.maximum_ratio

    def build(self, boundary):
        points = np.asarray(boundary, dtype=float)
        if len(points) < 4 or len(points) % 2:
            raise ValueError('Separator atlas needs an even prescribed rim')
        RegionUnion().build((Region(Key('separator_atlas', 0), points),))
        options = self.corners(points, minimum_corner_angle=self.minimum_corner_angle,
            minimum_polygon_sine=self.minimum_polygon_sine, include_smooth_corners=True)
        for score, corners in options[:self.max_candidates]:
            try:
                initial, counts = self.map(points, corners)
                report = quality(initial)
            except ValueError:
                continue
            variants = [(initial, dict(power=0, steps=0))]
            if report['maximum_edge_ratio'] <= 2.*self.maximum_ratio:
                for power in (2, 1):
                    candidate, improvement = FixedRimQuality(iterations=self.quality_iterations, power=power).improve(initial)
                    variants.append((candidate, dict(power=power, **improvement)))
            for candidate, improvement in variants:
                try:
                    physical = quality(candidate)
                    stored = quality(candidate.with_positions(publication_coordinates(candidate.xyz, 1., physical=True)))
                except ValueError:
                    continue
                if not self.accepted(physical) or not self.accepted(stored):
                    continue
                if not np.array_equal(candidate.xyz[candidate.seams[0], :2], points):
                    raise ValueError('Separator atlas changed its physical seam')
                edges = np.sort(np.stack((candidate.quads, np.roll(candidate.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
                unique, incidence = np.unique(edges, axis=0, return_counts=True)
                if (len(candidate.xyz)-len(unique)+len(candidate.quads) != 1
                        or not set(incidence).issubset({1, 2}) or np.sum(incidence == 1) != len(points)):
                    raise ValueError('Separator atlas is not a conforming quad disk')
                return ChartResult(candidate, dict(quality=physical, stored_quality=stored,
                    initial=report, improvement=improvement, atlas_corners=corners, rim_intervals=counts,
                    grass_quality_contract='positive corners and edge ratio <=8',
                    construction='fixed-rim narrow separator quad atlas'))
        raise ValueError('Separator atlas has no qualified stored-coordinate candidate')


def build_hole_disk(patterns, boundary, spacing, max_cells, quality_iterations):
    return _NarrowAtlas(patterns, max_cells=max_cells, max_candidates=128, quality_iterations=quality_iterations,
        minimum_corner_angle=.03, minimum_polygon_sine=.01, include_smooth_corners=True).build(boundary)
