

class Facade:

    def init(self, itemRenderers, globalRenderer):
        return
        #self.itemRenderers = itemRenderers
        #self.r = globalRenderer
        #self.app = globalRenderer.app
    
    def render(self, footprint):
        for facade in footprint.facades:
            footprint.element.l.genVolumes.createFace(footprint, facade.indices)


class RoofBase:

    def init(self, itemRenderers, globalRenderer):
        return

    def render(self, roofItem):
        for roofSide in roofItem.roofSides:
            roofItem.footprint.element.l.genVolumes.createFace(
                roofItem.footprint,
                roofSide.indices
            )