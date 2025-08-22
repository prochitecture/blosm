
class Footprint:

    def __init__(self):
        pass

    def init(self, itemRenderers, globalRenderer):
        self.itemRenderers = itemRenderers
        self.r = globalRenderer
        self.app = globalRenderer.app
    
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
        layer = footprint.building.element.l.bm.edges.layers.int["facade_class"]
        loop = face.loops[0]
        for facade in footprint.facades:
            loop.edge[layer] = facade.cl
            loop = loop.link_loop_next
            