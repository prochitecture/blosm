"""Strict original-boundary terrain sampling for the Python geometry backend."""
from dataclasses import dataclass

from .errors import PadlError


@dataclass(frozen=True)
class TerrainSample:
    minimum: float = 0
    maximum: float = 0
    applied: bool = False
    warning: str | None = None


class BlenderTerrainSampler:
    def __init__(self, terrain, *, policy='error', depsgraph=None):
        import bpy
        from mathutils.bvhtree import BVHTree
        if policy not in ('error', 'unprojected_with_warning'):
            raise ValueError('Unknown terrain policy')
        self.policy = policy
        self.bvh = None
        if terrain is None:
            self.requested = False
            return
        self.requested = True
        depsgraph = depsgraph or bpy.context.evaluated_depsgraph_get()
        obj = terrain.evaluated_get(depsgraph)
        mesh = obj.to_mesh()
        try:
            if mesh and mesh.polygons:
                vertices = [obj.matrix_world @ v.co for v in mesh.vertices]
                self.minimum = min(v.z for v in vertices)
                self.maximum = max(v.z for v in vertices)
                # Use Blender's evaluated tessellation, as GN Raycast does.
                # Re-tessellating nonplanar quads in BVHTree.FromPolygons can
                # choose another diagonal and give a different T_min.
                mesh.calc_loop_triangles()
                self.bvh = BVHTree.FromPolygons(vertices, [tuple(t.vertices) for t in mesh.loop_triangles], all_triangles=True)
        finally:
            obj.to_mesh_clear()

    def sample(self, footprint):
        from mathutils import Vector
        if not self.requested:
            return TerrainSample()
        hits = []
        if self.bvh is not None:
            margin = max(1, max(abs(self.minimum), abs(self.maximum)) * 1e-5)
            length = self.maximum - self.minimum + 2 * margin
            for ring in (footprint.outer, *footprint.holes):
                for vertex in ring:
                    hit, *_ = self.bvh.ray_cast(Vector((vertex[0], vertex[1], self.maximum + margin)), Vector((0, 0, -1)), length)
                    if hit is not None:
                        hits.append(hit.z)
        if len(hits) != sum(len(ring) for ring in (footprint.outer, *footprint.holes)):
            if self.policy == 'unprojected_with_warning':
                return TerrainSample(warning='terrain_coverage')
            raise PadlError(f'Terrain does not cover every boundary vertex of {footprint.identity}', code='terrain_coverage')
        return TerrainSample(min(hits), max(hits), True)
