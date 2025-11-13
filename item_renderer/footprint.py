import os

from util.blender import appendNodeGroupFromFile, addGeometryNodesModifier


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
            "Blosm Building Single Object"
        )

        # Load the node group for the case when each building will be generated as a separate Blender object
        self.gnBldgSeparateObject = appendNodeGroupFromFile(
            os.path.join(os.path.dirname(self.app.baseAssetPath), "prochitecture_buildings.blend"),
            "Blosm Building Separate Object"
        )
    
    def render(self, footprint, data):
        footprint.building.element.l.genFootprints.createFootprint(footprint, data, True)
        if footprint.innerPolygons:
            footprint.building.element.l.genHoles.createFootprintsForHoles(footprint)
    
    def renderPart(self, footprint):
        # same as <render>, but for building parts
        footprint.building.element.l.genPartFootprints.createFootprint(footprint, None, False)
        if footprint.innerPolygons:
            footprint.building.element.l.genHoles.createFootprintsForHoles(footprint)
    
    def finalize(self, layer):
        layer.genFootprints.finalize()
        addGeometryNodesModifier(
            layer.genFootprints.obj,
            self.gnBldgSingleObject if layer.singleObject else self.gnBldgSeparateObject,
            "Building Init"
        )