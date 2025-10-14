from mathutils import Vector
from .roof_flat import RoofFlat
from item.facade import Facade
from item.roof_flat_multi import RoofFlatMulti as ItemRoofFlatMulti
from blosm.building import BldgPolygonCW


class RoofMulti:
    
    def do(self, footprint):
        roofItem = self.init(footprint)
        if footprint.valid:
            if roofItem.innerPolygons:
                if self.renderAfterExtrude:
                    self.render(footprint, roofItem)
                else:
                    self.extrude(footprint, roofItem)
                    footprint.roofItem = roofItem
                    footprint.roofRenderer = self.roofRenderer
                footprint.bldgPart.footprint = footprint
            else:
                footprint.element.makePolygon()
                self.volumeAction.volumeGenerators[footprint.getStyleBlockAttr("roofShape")].do(
                    footprint
                )
    
    def extrudeInnerPolygons(self, footprint, roofItem):
        #
        # deal with the inner polygons below
        #
        facades = footprint.facades
        verts = footprint.element.l.genVolumes.verts
        indexOffset = len(verts)
        
        for innerPolygon in roofItem.innerPolygons:
            numVerts = innerPolygon.numEdges

            self._generateVerts(footprint, innerPolygon, roofItem)

            vectors = innerPolygon.getVectors()
            
            # the starting side
            _in = indexOffset+numVerts
            facades.append(
                Facade(
                    footprint,
                    (_in-1, indexOffset, _in, _in+numVerts-1),
                    next(vectors),
                    self,
                )
            )
            # the rest of the sides
            facades.extend(
                Facade(
                    footprint,
                    (indexOffset+i-1, indexOffset+i, _in+i, _in+i-1),
                    vector,
                    self
                ) for i,vector in zip(range(1, numVerts), vectors)
            )
            # mark the created facades as inner
            for i in range(-numVerts, 0):
                facades[i].outer = False
            
            indexOffset += 2*numVerts
    
    def init(self, footprint):
        data = self.data
        roofItem = super().init(footprint)
        if not footprint.valid:
            return
        z1 = footprint.minHeight
        element = footprint.element
        innerPolygons = roofItem.innerPolygons
        
        for _l in element.ls:
            if _l.role is data.outer:
                continue
            # create an inner polygon located at <minHeight>
            innerPolygon = BldgPolygonCW(_l, self.volumeAction.manager, footprint.building)
            if innerPolygon.numEdges < 3:
                continue
            innerPolygons.append(innerPolygon)
        return roofItem


class RoofFlatMulti(RoofMulti, RoofFlat):
    
    def __init__(self, data, volumeAction, itemRenderers):
        super().__init__("RoofFlatMulti", data, volumeAction, itemRenderers)

    def extrude(self, footprint, roofItem):
        super().extrude(footprint, roofItem)
        self.extrudeInnerPolygons(footprint, roofItem)
    
    def getRoofItem(self, footprint):
        return ItemRoofFlatMulti(
            footprint,
            self.getRoofFirstVertIndex(footprint)
        )