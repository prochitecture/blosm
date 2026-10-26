"""Direct quad charts in angular coordinates of a prescribed physical rim.

The fitted centre chooses coordinates only. Road/junction samples, their
incidence, positions and heights remain authoritative in the physical plane.
"""
import math
import numpy as np
from .mesh import Mesh
from .model import Key
from .regions import Region, RegionUnion
from .quad_charts import ChartResult, LayeredQuadChart, FixedRimQuality
from .four_side_charts import FourSidedQuadChart


class AngularSeparatorChart:
    def __init__(self, patterns, *, spacing=.7, max_cells=50000):
        if not math.isfinite(spacing) or spacing <= 0. or type(max_cells) is not int or max_cells < 4:
            raise ValueError('Invalid angular separator chart limits')
        self.patterns, self.spacing, self.max_cells = patterns, spacing, max_cells

    def build(self, boundary):
        boundary = np.asarray(boundary, dtype=float)
        if (boundary.ndim != 2 or boundary.shape[1] != 2 or len(boundary) < 4
                or len(boundary) % 2 or not np.all(np.isfinite(boundary))):
            raise ValueError('Angular separator needs an even simple physical rim')
        RegionUnion().build((Region(Key('separator_chart', 0), boundary),))
        origin = boundary.mean(axis=0)
        relative = boundary-origin
        fitted, _, rank, _ = np.linalg.lstsq(np.c_[2.*relative, np.ones(len(relative))],
                                            np.sum(relative*relative, axis=1), rcond=None)
        if rank != 3:
            raise ValueError('Separator rim has no determinate angular centre')
        centre = origin+fitted[:2]
        delta = boundary-centre
        radius = np.linalg.norm(delta, axis=1)
        scale = float(radius.mean())
        if not math.isfinite(scale) or radius.min() <= 1e-6:
            raise ValueError('Angular separator has insufficient radial clearance')
        theta = np.arctan2(delta[:, 1], delta[:, 0])
        ordered = np.sort(theta)
        gaps = np.diff(np.r_[ordered, ordered[0]+2.*math.pi])
        index = int(np.argmax(gaps))
        cut = ordered[index]+.5*gaps[index]
        angle = cut+math.pi
        local_theta = np.arctan2(np.sin(theta-angle), np.cos(theta-angle))
        if np.max(abs(np.diff(np.r_[local_theta, local_theta[0]]))) >= math.pi:
            raise ValueError('A physical rim edge crosses the angular chart cut')
        local = np.c_[scale*local_theta, scale-radius]
        # This chart changes coordinates only. Its complete boundary and
        # final Jacobians are checked again in the physical XY plane.
        spacing = self.spacing*min(1., scale/float(radius.max()))
        failures = []
        from .separator_surfaces import PartitionedSeparatorChart
        partition_spacing=max(spacing, float(np.linalg.norm(np.roll(local,-1,axis=0)-local,axis=1).max())*(1.+1e-8))
        builders = (
            ('layered', LayeredQuadChart(self.patterns, spacing=spacing,
                max_cells=self.max_cells, quality_iterations=100)),
            ('four-sided', FourSidedQuadChart(self.patterns, max_cells=self.max_cells, quality_iterations=100)),
            ('partitioned', PartitionedSeparatorChart(self.patterns, spacing=partition_spacing,
                interior_spacing=spacing, max_cells=self.max_cells)))
        extent = np.ptp(local, axis=0)
        aspect = float(extent.max()/max(float(extent.min()), 1e-12))
        # Long thin runs need longitudinal blocks rather than many radial
        # layers around one centre. Keep the compact-chart alternatives as
        # fallbacks; all choices retain the same complete physical rim.
        if aspect > 4.:
            builders = (builders[2], *builders[:2])
        for label, builder in builders:
            try:
                chart = builder.build(local)
                parameters = chart.mesh.xyz[:, :2]
                r = scale-parameters[:, 1]
                if r.min() <= 1e-6:
                    raise ValueError('Angular separator chart reaches its origin')
                a = angle+parameters[:, 0]/scale
                xy = centre+r[:, None]*np.c_[np.cos(a), np.sin(a)]
                xy[:len(boundary)] = boundary
                mesh = Mesh(np.c_[xy, np.zeros(len(xy))], chart.mesh.quads.copy(), 'grass',
                            {0:np.arange(len(boundary))})
                quality = mesh.validate()
                improvement = None
                if quality['maximum_edge_ratio'] > 8.:
                    mesh, improvement = FixedRimQuality(iterations=100,target_ratio=8.).improve(mesh)
                    quality = mesh.validate()
                if quality['maximum_edge_ratio'] > 8.:
                    raise ValueError('Mapped angular separator exceeds physical edge ratio 8')
                assert np.array_equal(mesh.xyz[:len(boundary), :2], boundary)
                return ChartResult(mesh, dict(construction='angular coordinates of the prescribed physical separator rim',
                    centre=centre.tolist(), scale=scale, branch_angle=angle,
                    angular_span=float(np.ptp(local_theta)), inner_chart=label,
                    local_aspect_ratio=aspect, partitioned_first=aspect>4.,
                    inner_metadata=chart.metadata, physical_quality=quality,
                    physical_improvement=improvement, failed_charts=failures,
                    physical_boundary_exact=True, all_direct_quads=True))
            except ValueError as error:
                failures.append(label+': '+str(error))
        raise ValueError('No accepted angular separator chart: '+'; '.join(failures))
