
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

        self.r.createFace(
            footprint,
            range(firstVertIndex, firstVertIndex + footprint.polygon.numEdges)
        )