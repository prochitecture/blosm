"""Main-thread terrain sampling without changing terrain geometry or transforms."""
import math

from mathutils import Vector


class GroundSurface:
    def __init__(self, context):
        self.context, self.scene = context, context.scene
        name = self.scene.blosm.terrainObject.strip()
        self.terrain = self.scene.objects.get(name) if name else None
        if self.terrain is not None and self.terrain.type != 'MESH':
            raise ValueError('The selected Mapillary terrain must be a mesh')
        self.depsgraph = None

    def validate(self):
        if self.terrain is not None:
            try:
                if self.scene.objects.get(self.terrain.name) == self.terrain:
                    return
            except ReferenceError:
                pass
            raise RuntimeError('The Mapillary terrain was deleted or unlinked from the import scene')

    def height(self, x, y):
        self.validate()
        if self.terrain is None:
            return 0., 'scene_plane'
        # Acquire evaluated geometry lazily within a cooperative camera turn.
        # Object.ray_cast uses Blender's cached mesh acceleration structure;
        # do not rebuild a Python BVH or force a depsgraph update per camera.
        if self.depsgraph is None:
            self.depsgraph = self.context.evaluated_depsgraph_get()
        terrain = self.terrain.evaluated_get(self.depsgraph)
        matrix = terrain.matrix_world
        try:
            inverse = matrix.inverted()
        except ValueError:
            raise RuntimeError('The Mapillary terrain has a singular transform') from None
        top = max((matrix @ Vector(corner)).z for corner in terrain.bound_box) + 1.
        origin = inverse @ Vector((x, y, top))
        direction = (inverse.to_3x3() @ Vector((0., 0., -1.))).normalized()
        hit, location, _, _ = terrain.ray_cast(origin, direction)
        if hit:
            z = (matrix @ location).z
            if not math.isfinite(z):
                raise RuntimeError('The Mapillary terrain returned an invalid elevation')
            return z, 'terrain_surface'
        return 0., 'scene_plane_outside_terrain'
