import math
from itertools import accumulate
from operator import add
from .roof_flat_multi import RoofMulti
from .roof_hipped import RoofHipped
from .roof_flat import RoofLeveled
from item.roof_hipped_multi import RoofHippedMulti as ItemRoofHippedMulti
from mathutils import Vector

from lib.bpypolyskel.bpypolyskel import polygonize

from util import zAxis

#from util.debug import dumpInputHippedRoof


class RoofHippedMulti(RoofMulti, RoofHipped):
    
    def __init__(self, data, volumeAction, itemRenderers):
        super().__init__(data, volumeAction, itemRenderers)
        
        # Python dictionary used for mapping:
        # the index of the first face vertex -> counter index;
        # the face is formed by the straight skeleton;
        # counter indices: 0 for the outer polygon, 1 for the first hole and so on;
        # the dictionary is used if there two or more holes in the polygon
        self.faceToContourIndex = {}
    
    def getRoofItem(self, footprint):
        return ItemRoofHippedMulti(
            footprint,
            self.getRoofFirstVertIndex(footprint)
        )
    
    def extrudeInnerPolygons(self, footprint, roofItem):
        if footprint.noWalls:
            z = footprint.roofVerticalPosition
            # the basement of the roof
            footprint.building.verts.extend(
                Vector((v.x, v.y, z)) for innerPolygon in roofItem.innerPolygons for v in innerPolygon.verts
            )
            return
        super().extrudeInnerPolygons(footprint, roofItem)
        
    def extrude(self, footprint, roofItem):
        # <firstVertIndex> is the index of the first vertex of the polygon that defines the roof base
        firstVertIndex = self.getRoofFirstVertIndex(footprint)
        
        # extrude the outer polygon
        RoofLeveled.extrude(self, footprint, roofItem)
        self.extrudeInnerPolygons(footprint, roofItem)
        
        # now generate the roof
        ok = self.generateRoof(footprint, roofItem, firstVertIndex)
        
        if not ok:
            # Unable to generate the roof.
            # Generate a flat roof as a fallback solution
            self.volumeAction.volumeGeneratorMultiFlat.do(footprint)
        
        return ok
    
    def render(self, footprint, roofItem):
        if self.extrude(footprint, roofItem):
            self.facadeRenderer.render(footprint)
            self.roofRenderer.render(roofItem)
    
    def generateRoof(self, footprint, roofItem, firstVertIndex):
        verts = footprint.element.l.genVolumes.verts
        numPolygonVerts = footprint.polygon.numEdges
        innerPolygons = roofItem.innerPolygons
        numHoles = len(innerPolygons)
        
        roofSideIndices = []
        
        holesInfo = []
        
        # the outer contour
        unitVectors = [
            vector.unitVector3d for vector in footprint.polygon.getVectors()
        ]
        lengths = [
            vector.length for vector in footprint.polygon.getVectors()
        ]
        
        if footprint.noWalls:
            _offset = firstVertIndex + numPolygonVerts
            holesInfo.append((_offset, innerPolygons[0].numEdges))
            holesInfo.extend(
                zip(
                    (_offset + v for v in accumulate( (innerPolygons[i].numEdges for i in range(numHoles-1)), add)), (innerPolygons[i].numEdges for i in range(1, numHoles))
                )
            )
        else:
            _offset = firstVertIndex + numPolygonVerts + innerPolygons[0].numEdges # FIXME
            holesInfo.append((_offset, innerPolygons[0].numEdges))
            holesInfo.extend(
                zip(
                    (_offset + v for v in accumulate( (innerPolygons[i].numEdges + innerPolygons[i+1].numEdges for i in range(numHoles-1)), add)), (innerPolygons[i].numEdges for i in range(1, numHoles))
                )
            )
        
        # the holes
        for innerPolygon in innerPolygons:
            unitVectors.extend(
                vector.unitVector3d for vector in innerPolygon.getVectors()
            )
            lengths.extend(
                vector.length for vector in innerPolygon.getVectors()
            )

        if numHoles > 1:
            faceToContourIndex = self.faceToContourIndex
            faceToContourIndex.clear()
            for i in range(firstVertIndex, firstVertIndex+numPolygonVerts):
                faceToContourIndex[i] = firstVertIndex
            _offset = numPolygonVerts
            for firstVertIndexHole,numVertsHole in holesInfo:
                for i in range(firstVertIndexHole, firstVertIndexHole+numVertsHole):
                    faceToContourIndex[i] = firstVertIndexHole - _offset
                _offset += numVertsHole
        
        #dumpInputHippedRoof(verts, firstVertIndex, numPolygonVerts, holesInfo, unitVector)
        #return
        
        # calculate polygons formed by the straight skeleton
        polygonize(
            verts,
            firstVertIndex,
            numPolygonVerts,
            holesInfo,
            footprint.roofHeight,
            0,
            roofSideIndices,
            unitVectors
        )
        
        if not self.validatePolygonizeOutput(roofSideIndices):
            return False
        
        roofVerticalPosition = verts[firstVertIndex][2]
        
        # calculate tangent of the roof pitch angle
        tan = ( verts[ roofSideIndices[0][2] ][2] - roofVerticalPosition ) / \
        (verts[ roofSideIndices[0][2] ] - verts[ roofSideIndices[0][1] ]).dot( zAxis.cross(unitVectors[0]) )
        factor = math.sqrt(1. + tan*tan)
        
        for indices in roofSideIndices:
            if numHoles == 1:
                edgeIndex = indices[0] - holesInfo[0][0] + numPolygonVerts\
                    if indices[0] >= holesInfo[0][0] else\
                    indices[0] - firstVertIndex
            else:
                edgeIndex = indices[0] - faceToContourIndex[indices[0]]
            roofItem.addRoofSide(
                indices,
                # UV-coordinates
                ( (0., 0.), (lengths[edgeIndex], 0.) ) + tuple(
                    (
                        (verts[ indices[_index] ] - verts[ indices[0] ]).dot(unitVectors[edgeIndex]),
                        (verts[ indices[_index] ][2] - roofVerticalPosition) * factor
                    ) for _index in range(2, len(indices))
                ),
                edgeIndex
            )
        
        return True