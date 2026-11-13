"""Main-thread terrain sampling without changing terrain geometry or transforms."""
import math

import bpy
from mathutils import Vector
from mathutils.bvhtree import BVHTree


class GroundSurface:
    def __init__(self, context):
        self.context, self.scene = context, context.scene
        name = self.scene.blosm.terrainObject.strip()
        self.terrain = self.scene.objects.get(name) if name else None
        if self.terrain is not None and self.terrain.type != 'MESH':
            raise ValueError('The selected Mapillary terrain must be a mesh')
        self.tree = None
        self.matrix = self.inverse = None
        self.project_z = None

    def snapshot(self):
        """Evaluate a temporary visible copy, then keep only Python/math data.

        A hidden/excluded source can have a Mesh datablock but no evaluated mesh
        for Object.ray_cast. A private copy also evaluates its modifiers without
        changing the user's visibility, mesh, selection, or object transform.
        This scope never yields: all temporary IDs are removed before an undo
        checkpoint or viewport redraw can observe them.
        """
        collection = proxy = None
        try:
            collection = bpy.data.collections.new('.Mapillary terrain evaluation')
            self.scene.collection.children.link(collection)
            proxy = self.terrain.copy()
            collection.objects.link(proxy)
            proxy.hide_viewport = False
            proxy.hide_render = True
            proxy.hide_set(False, view_layer=self.context.view_layer)
            depsgraph = self.context.evaluated_depsgraph_get()
            depsgraph.update()
            evaluated = proxy.evaluated_get(depsgraph)
            matrix = evaluated.matrix_world.copy()
            try:
                inverse = matrix.inverted()
            except ValueError:
                raise RuntimeError('The Mapillary terrain has a singular transform') from None
            try:
                tree = BVHTree.FromObject(proxy, depsgraph)
            except (ValueError, RuntimeError):
                raise RuntimeError('The selected Mapillary terrain could not be evaluated as a mesh') from None
            if tree is None:
                raise RuntimeError('The selected Mapillary terrain could not be evaluated as a mesh')
            self.project_z = max((matrix @ Vector(corner)).z for corner in evaluated.bound_box) + 1.
            self.tree, self.matrix, self.inverse = tree, matrix, inverse
        finally:
            if proxy is not None:
                bpy.data.objects.remove(proxy, do_unlink=True)
            if collection is not None:
                bpy.data.collections.remove(collection)

    def clear(self):
        self.tree = self.matrix = self.inverse = self.project_z = None

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
        # Snapshot once, lazily inside the first camera's main-thread turn.
        # The independent BVH survives depsgraph updates and undo checkpoints;
        # subsequent cameras need neither reevaluation nor temporary objects.
        if self.tree is None:
            self.snapshot()
        origin = self.inverse @ Vector((x, y, self.project_z))
        direction = (self.inverse.to_3x3() @ Vector((0., 0., -1.))).normalized()
        location, _, _, _ = self.tree.ray_cast(origin, direction)
        if location is not None:
            z = (self.matrix @ location).z
            if not math.isfinite(z):
                raise RuntimeError('The Mapillary terrain returned an invalid elevation')
            return z, 'terrain_surface'
        return 0., 'scene_plane_outside_terrain'
