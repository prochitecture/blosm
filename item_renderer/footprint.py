import os

from blosm.util.blender import appendNodeGroupFromFile


class Footprint:

    def __init__(self):
        pass

    def init(self, itemRenderers, globalRenderer):
        self.itemRenderers = itemRenderers
        self.r = globalRenderer
        self.app = globalRenderer.app

        # Load the node group for the case when buildings will be generated as a single Blender object
        self.gnBldgSingleObject = appendNodeGroupFromFile(
            os.path.join(os.path.dirname(self.app.baseAssetPath), "prochitecture_buildings.blend"),
            "Blosm Building Single object"
        )

        # Load the node group for the case when each building will be generated as a separate Blender object
        self.gnBldgSeparateObject = appendNodeGroupFromFile(
            os.path.join(os.path.dirname(self.app.baseAssetPath), "prochitecture_buildings.blend"),
            "Blosm Building Separate Object"
        )
    
    def render(self, footprint):
        # get the first vertex index from the first facade
        firstVertIndex = footprint.facades[0].indices[0]

        face = self.r.createFace(
            footprint,
            range(firstVertIndex, firstVertIndex + footprint.polygon.numEdges)
        )
        # set the attribute "roof_shape" to <face>
        # 1: flat
        # 4: skillion
        layer = footprint.building.element.l.bm.faces.layers.int["roof_shape"]
        face[layer] = 4 if footprint.getStyleBlockAttr("roofShape") == "skillion" else 1

        # set the attribute "facade_class" to each edge of <face>
        bmLayer = footprint.building.element.l.bm.edges.layers.int["facade_class"]
        loop = face.loops[0]
        for facade in footprint.facades:
            loop.edge[bmLayer] = facade.cl
            loop = loop.link_loop_next
    
    def finalize(layer):
        pass
        