"""Local physical patches with radial quad grading and fixed profiles."""
import numpy as np
from .mesh import Mesh
from .quad_charts import LayeredQuadChart
from .regions import signed_area
from .surface_quality import ProfiledQuadCorrection, ProfiledSlopeCorrection


class LayeredAngularPatch:
    """Retain concave cap/curb steps without forcing one four-side atlas."""
    def __init__(self, builder, frame, field):
        self.builder, self.frame, self.field = builder, frame, field

    def build(self, outer, seam, name):
        physical, kinds, boundary_keys = outer
        shared, shared_keys = seam
        if (np.linalg.norm(physical[0]-shared[0]) > 1e-8
                or np.linalg.norm(physical[-1]-shared[-1]) > 1e-8):
            raise ValueError('A layered angular patch has unmatched crown endpoints')
        rim = np.vstack((physical, shared[-2:0:-1]))
        if len(rim) % 2:
            raise ValueError('A layered angular patch needs coordinated even boundary parity')
        boundary_keys = list(boundary_keys)
        boundary_keys[0], boundary_keys[-1] = shared_keys[0], shared_keys[-1]
        keys = boundary_keys+shared_keys[-2:0:-1]
        values = self.field.values(rim)
        for index, kind in enumerate(kinds):
            if kind[0] == 'outer':
                values[index:index+2] = self.field.datum
        values[len(physical):] = self.field.ceiling
        values[0] = values[len(physical)-1] = self.field.ceiling
        if signed_area(rim) < 0.:
            rim, keys, values = rim[::-1].copy(), keys[::-1], values[::-1].copy()
        chart = LayeredQuadChart(self.builder.patterns, spacing=self.builder.spacing,
            max_cells=self.builder.max_cells, quality_iterations=100).build(
                self.frame.local(rim), require_quality=False)
        xy = self.frame.world(chart.mesh.xyz[:, :2])
        xy[:len(rim)] = rim
        heights = self.field.values(xy)
        heights[:len(rim)] = values
        mesh = Mesh(np.c_[xy, heights], chart.mesh.quads)
        mesh, grading = ProfiledQuadCorrection(max_cells=self.builder.max_cells,
            max_iterations=500).improve(mesh, self.field.values, datum=self.field.datum,
                ceiling=self.field.ceiling, exaggeration=1., maximum_slope=self.builder.maximum_slope)
        mesh, slopes = ProfiledSlopeCorrection(max_cells=self.builder.max_cells).improve(
            mesh, self.field.values, datum=self.field.datum, ceiling=self.field.ceiling,
            maximum_slope=self.builder.maximum_slope, longitudinal_bound=self.frame.longitudinal_bound,
            maximum_longitudinal_slope=self.builder.maximum_longitudinal_slope)
        if not np.array_equal(mesh.xyz[:len(rim)], np.c_[rim, values]):
            raise ValueError('Layered angular grading changed a prescribed physical boundary')
        keys.extend((name, 'layered inside', index) for index in range(len(rim), len(mesh.xyz)))
        return mesh, keys, dict(construction='local physical radial quad layers',
            chart=chart.metadata, world_quality=mesh.validate(), boundary_vertices=len(rim),
            exact_prescribed_xyz=True, grading=grading, profile_slope_correction=slopes,
            all_direct_quads=True)
