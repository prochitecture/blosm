"""Keep physical slope bounds after Blender's float32 vertex storage.

Only an actual rounded-buffer violation triggers a correction. XY, topology,
boundary profiles and curb heights remain fixed; the numerical query uses the
same corrected physical heights as the mesh later sent for publication.
"""
import numpy as np
from .drainage import BoundaryDrainage
from .mesh import Mesh
from .slope_projection import SlopeProjection, corner_matrices
from .surface_query import QuadHeightField


class PublicationSlopeGuard:
    tolerance = 1e-7
    maximum_rounds = 8

    @staticmethod
    def physical_buffer(mesh, datum, exaggeration):
        xyz = mesh.xyz.copy()
        xyz[:, 2] = datum + exaggeration * (xyz[:, 2] - datum)
        xyz = xyz.astype(np.float32).astype(np.float64)
        xyz[:, 2] = datum + (xyz[:, 2] - datum) / exaggeration
        return Mesh(xyz, mesh.quads, mesh.material, dict(mesh.seams))

    def correct(self, mesh, field, fixed, *, exaggeration, limit, max_cells, check_joins,
                allow_projection=True):
        if exaggeration == 0.:
            return mesh, field, None
        rounded = self.physical_buffer(mesh, field.datum, exaggeration)
        initial = float(rounded.quad_slopes().max())
        if initial <= limit + self.tolerance:
            return mesh, field, None
        if not allow_projection:
            raise ValueError('Publication slope correction is required but height projection is disabled')
        original = mesh
        floors = np.full(len(mesh.xyz), field.datum)
        rounds = []
        for _ in range(self.maximum_rounds):
            def gradients(candidate):
                z = candidate.xyz[candidate.quads, 2] - field.datum
                return np.einsum('fcij,fj->fci', corner_matrices(candidate), z)
            error = float(np.linalg.norm(gradients(rounded) - gradients(mesh), axis=-1).max())
            margin = max(1e-6, 4. * error)
            target = limit - margin
            if target <= 0.:
                raise ValueError('Published coordinate precision cannot represent the physical slope limit')
            mesh, projection = SlopeProjection(max_cells=max_cells).solve(mesh, fixed,
                datum=field.datum, ceiling=field.ceiling, maximum_slope=target,
                minimum_heights=floors)
            filled, drainage = BoundaryDrainage(max_cells=max_cells).fill(mesh, fixed)
            rounded = self.physical_buffer(mesh, field.datum, exaggeration)
            maximum = float(rounded.quad_slopes().max())
            rounds.append(dict(roundoff_gradient_error=error, target=target,
                               published_maximum=maximum, projection=projection, drainage=drainage))
            if maximum <= limit + self.tolerance and drainage['maximum_fill'] <= 1e-6:
                check_joins(mesh)
                if (not np.array_equal(mesh.xyz[fixed], original.xyz[fixed])
                        or not np.array_equal(mesh.xyz[:, :2], original.xyz[:, :2])
                        or not np.array_equal(mesh.quads, original.quads)):
                    raise ValueError('Publication slope correction changed its fixed geometry')
                details = dict(initial_published_maximum=initial, published_maximum=maximum,
                    maximum_height_change=float(abs(mesh.xyz[:, 2] - original.xyz[:, 2]).max()),
                    exact_boundary=True, exact_xy=True, exact_quads=True, rounds=rounds)
                query = QuadHeightField(mesh, datum=field.datum, ceiling=field.ceiling,
                    metadata=dict(field.metadata, publication_roundoff=details), max_cells=max_cells)
                return mesh, query, details
            raised = filled > mesh.xyz[:, 2] + 1e-10
            floors[raised] = np.maximum(floors[raised], filled[raised])
        raise ValueError('Publication slope correction exceeded its iteration budget')
