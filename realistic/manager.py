"""
This file is a part of Blosm addon for Blender.
Copyright (C) 2014-2018 Vladimir Elistratov
prokitektura+support@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from .. import parse
from ..manager import BaseManager
from .renderer import AreaRenderer#, TerrainRenderer
from ..renderer.layer import MeshLayer
from contextlib import contextmanager


@contextmanager
def terrain_changes(terrain):
    """Rollback only the existing terrain slots touched by a synchronous unit."""
    import bpy
    from ..util.blender import makeActive
    modifiers = {m.as_pointer() for m in terrain.modifiers}
    particles = {p.as_pointer() for p in terrain.particle_systems}
    groups = {g.as_pointer() for g in terrain.vertex_groups}
    counts = [(p.settings, p.settings.count, p.settings.display_percentage) for p in terrain.particle_systems]
    surfaces = {m.as_pointer(): {s.as_pointer() for s in m.canvas_settings.canvas_surfaces}
                for m in terrain.modifiers if m.type == 'DYNAMIC_PAINT' and m.canvas_settings}
    try:
        yield
    except Exception:
        makeActive(terrain)
        for index in reversed(range(len(terrain.particle_systems))):
            if terrain.particle_systems[index].as_pointer() not in particles:
                terrain.particle_systems.active_index = index
                bpy.ops.object.particle_system_remove()
        for modifier in list(terrain.modifiers):
            if modifier.as_pointer() not in modifiers:
                terrain.modifiers.remove(modifier)
            elif modifier.as_pointer() in surfaces:
                collection = modifier.canvas_settings.canvas_surfaces
                for index in reversed(range(len(collection))):
                    if collection[index].as_pointer() not in surfaces[modifier.as_pointer()]:
                        collection.active_index = index
                        bpy.ops.dpaint.surface_slot_remove()
        for group in list(terrain.vertex_groups):
            if group.as_pointer() not in groups:
                terrain.vertex_groups.remove(group)
        for settings, count, display in counts:
            settings.count, settings.display_percentage = count, display
        raise


class ProgressiveAreaLayer(MeshLayer):
    def finalizeBlenderObject(self, obj):
        # AreaRenderer applies intersection/slicing/shrinkwrap once, before the
        # area and its scattering changes are committed together.
        pass


class AreaManager(BaseManager):
    
    def __init__(self, osm, app, defaultRenderer, **renderers):
        """
        Args:
            osm (parse.Osm): Parsed OSM data
        """
        super().__init__(osm)
        self.app = app
        self.defaultRenderer = defaultRenderer
        self.renderers = renderers
        #self.terrainRenderer = TerrainRenderer()

    def createLayer(self, layerId, app, **kwargs):
        if getattr(app, 'osm_session', None) and app.osm_session.requested:
            return app.createLayer(layerId, ProgressiveAreaLayer, **kwargs)
        if not app.singleObject:
            # patching kwargs: forcing single object mode
            kwargs["singleObject"] = True
        return super().createLayer(layerId, app, **kwargs)

    def render_steps(self, session=None):
        if session is None or not session.requested:
            yield from super().render_steps(session)
            return
        from copy import copy
        from math import ceil, floor
        import bpy
        from ..util.blender import getModifier
        from ..renderer import Renderer
        from .renderer import ForestRenderer
        self._area_states = {}
        osm = self.osm
        for elements in (osm.relations.values(), osm.ways.values()):
            for element in elements:
                if not element.valid or not element.r:
                    continue
                def render_area(element=element):
                    layer = element.l
                    renderer = element.rr or self.renderer
                    renderer.preRender(element)
                    if element.t is parse.multipolygon:
                        renderer.renderMultiPolygon(element, osm)
                    else:
                        renderer.renderPolygon(element, osm)
                    renderer.postRender(element)
                    area = copy(layer)
                    area.obj = renderer.obj
                    area_renderer = self.renderers.get(layer.id, self.defaultRenderer)
                    area_renderer.finalizeBlenderObject(area, self.app)
                    terrain = self.app.terrain.terrain
                    state = self._area_states.get(layer.id)
                    if state is None:
                        # Names are only labels; retain actual new identities.
                        area.id = '%s_%s' % (Renderer.collection.name, layer.id)
                        AreaRenderer.addSubsurfModifier(terrain)
                        AreaRenderer.beginDynamicPaintCanvas(terrain)
                        area_renderer.renderTerrain(area, terrain)
                        AreaRenderer.endDynamicPaintCanvas(terrain)
                        area.copyParticleSettings = True
                        area_renderer.renderArea(area, self.app)
                        area.brush_collection = getModifier(terrain, 'DYNAMIC_PAINT').canvas_settings.canvas_surfaces[-1].brush_collection
                        if isinstance(area_renderer, ForestRenderer):
                            area.particle_system = terrain.particle_systems[-1]
                        self._area_states[layer.id] = area
                    else:
                        state.brush_collection.objects.link(area.obj)
                        AreaRenderer.renderArea(area_renderer, area, self.app)
                        state.surfaceArea += getattr(area, 'surfaceArea', 0.)
                        if isinstance(area_renderer, ForestRenderer):
                            count = ceil(self.app.treeDensity / 10000 * state.surfaceArea)
                            state.particle_system.settings.count = count
                            state.particle_system.settings.display_percentage = min(100, floor(area_renderer.maxDisplayCount/max(count, 1)*100))
                    # This unit intentionally adds its owned scattering state;
                    # subsequent external terrain edits still invalidate input.
                    session.terrain_signature = session._terrain_key()
                    session.terrain_dirty = False
                def render():
                    with terrain_changes(self.app.terrain.terrain):
                        render_area()
                yield from session.unit(render)
        # Tree nodes retain their independent publication units.
        for node in osm.rNodes:
            yield from session.node_unit(node, lambda node=node: (node.rr or self.nodeRenderer).renderNode(node, osm))

    def parseWay(self, element, elementId):
        # exactly the same code as in manager.Polygon.parseWay(..)
        if element.closed:
            element.t = parse.polygon
            # render it in <BaseManager.render(..)>
            element.r = True
        else:
            element.valid = False
    
    def renderExtra(self):
        if getattr(self.app, 'osm_session', None) and self.app.osm_session.requested:
            return
        app = self.app
        # Blender object for the terrain
        terrain = app.terrain.terrain
        
        # a counter for the number of valid area layers
        numLayers = 0
        # set renderer for each layer
        for layer in app.layers:
            layerId = layer.id
            if layer.area and layer.obj:
                numLayers += 1
                # find a render for <layer>
                for _layerId in self.renderers:
                    if _layerId.startswith(layerId):
                        renderer = self.renderers[layerId]
                        break
                else:
                    renderer = self.defaultRenderer
                # save <renderer> in <layer> for future references
                layer.renderer = renderer
        if not numLayers:
            return
        
        for layer in app.layers:
            if layer.area and layer.obj and layer.renderer:
                layer.renderer.finalizeBlenderObject(layer, app)
        
        AreaRenderer.addSubsurfModifier(terrain)
        AreaRenderer.beginDynamicPaintCanvas(terrain)
                
        for layer in app.layers:
            if layer.area and layer.obj and layer.renderer:
                layer.renderer.renderTerrain(layer, terrain)
        
        AreaRenderer.endDynamicPaintCanvas(terrain)
        
        for layer in app.layers:
            if layer.area and layer.obj and layer.renderer:
                layer.renderer.renderArea(layer, app)
        
        # set material for the terrain
        #self.terrainRenderer.render(app)
        
        terrain.select_set(False)
