import math
from mathutils import Vector
from .roof_flat import RoofLeveled
from item.roof_hipped import RoofHipped as ItemRoofHipped

from lib.bpypolyskel.bpypolyskel import polygonize

from util import zAxis

#from util.debug import dumpInputHippedRoof

# auxiliary indices to deal with quadrangles
_prevIndices = (3, 0, 1, 2)
_indices = (0, 1, 2, 3)
_nextIndices = (1, 2, 3, 0)
_oppositeIndices = (2, 3, 0, 1)


class RoofHipped(RoofLeveled):
    
    height = 4.
    
    def __init__(self, data, volumeAction, itemRenderers):
        super().__init__("RoofHipped", data, volumeAction, itemRenderers)
        
        self.extrudeTillRoof = True
        self.gnForFootprintAvailable = False
    
    def getRoofItem(self, footprint):
        return ItemRoofHipped(footprint)

    def validate(self, footprint):
        """
        Additional validation
        """
        if not footprint.roofHeight:
            footprint.valid = False
    
    def validatePolygonizeOutput(self, roofSideIndices):
        for indices in roofSideIndices:
            if not ( len(indices) >= 3 and len(indices) == len(set(indices)) ):
                return False
        return True
    
    def extrude(self, footprint, roofItem):
        # <firstVertIndex> is the index of the first vertex of the polygon that defines the roof base
        firstVertIndex = self.getRoofFirstVertIndex(footprint)
        
        super().extrude(footprint, roofItem)
        
        # now generate the roof
        if footprint.polygon.numEdges == 4:
            ok = self.generateRoofQuadrangle(footprint, roofItem, firstVertIndex)
        else:
            ok = self.generateRoof(footprint, roofItem, firstVertIndex)
        
        if not ok:
            # Unable to generate the hipped roof.
            # Generate a flat roof as a fallback solution
            self.volumeAction.volumeGenerators["flat"].do(
                footprint
            )
        
        return ok

    def render(self, footprint, roofItem):
        if self.extrude(self, footprint, roofItem):
            self.facadeRenderer.render(footprint, self.data)
            self.roofRenderer.render(roofItem)
    
    def generateRoofQuadrangle(self, footprint, roofItem, firstVertIndex):
        def _getRoofVert(vert, i, tan):
            """
            Assuming the edge event occured for the polygon edge with the index <i>,
            get the vertex of the edge event (i.e. where the bisectors originating
            from the polygon edge vertices meet).
            <tan> defines the tangent of the roof pitch angle, so an additional offset
            <self.distance[i] * tan * zAxis> is added to the resulting vertex.
            """
            factor = (1. + cos[i]) / sin[i]
            return\
                vert + \
                distances[i]/lengths[i] * \
                Vector((-vectors[i][1] + factor*vectors[i][0], vectors[i][0] + factor*vectors[i][1], 0.)) + \
                + \
                distances[i] * tan * zAxis

        verts = footprint.building.element.l.genVolumes.verts

        # Vectors along the edges of the footprint polygon
        vectors = [
            v.vector for v in footprint.polygon.getVectors()
        ]
        
        # Lengths of the edges of the footprint polygon
        lengths = [
            v.length for v in footprint.polygon.getVectors()
        ]
        
        # Cosines of the angles of the footprint polygon
        cos = [
            - v.unitVector.dot(v.prev.unitVector) for v in footprint.polygon.getVectors()
        ]
        
        # Sines of the angles of the footprint polygon
        sin = [
            - v.unitVector.cross(v.prev.unitVector) for v in footprint.polygon.getVectors()
        ]
        
        # Distances between the edge of the footprint polygon and
        # the point of the edge event (i.e. where the bisectors originating
        # from the edge vertices meet)
        distances = [
            lengths[i]/( (1.+cos[i])/sin[i] + (1.+cos[_next])/sin[_next] ) \
                for i, _next in zip(_indices, _nextIndices)
        ]
        
        if distances[0] == distances[1]:
            # the special case of the square footprint
            return
        else:
            # The first of the two newly created vertices of the hipped roof:
            # we find the very first occurance of the edge event,
            # namely the index of the polygon edge with the minimum <distance>
            minDistanceIndex1 = min(_indices, key = lambda i: distances[i])
            minDistanceIndex1Next = _nextIndices[minDistanceIndex1]
            # The second of the two newly created vertices of the hipped roof.
            # Note that it's an assumption that the other edge event for the quadrangle
            # occurs for the polygon edge opposite to the edge with the index <minDistanceIndex1>
            minDistanceIndex2 = _oppositeIndices[minDistanceIndex1]
            minDistanceIndex2Next = _nextIndices[minDistanceIndex2]
            
            # tangent of the roof pitch angle
            tan = footprint.roofHeight / max(distances[minDistanceIndex1], distances[minDistanceIndex2])
            factor = math.sqrt(1. + tan*tan)
            
            # add two new vertices to the Python list <verts>
            vertIndex1 = len(verts)
            verts.append( _getRoofVert(verts[firstVertIndex + minDistanceIndex1], minDistanceIndex1, tan) )
            
            vertIndex2 = vertIndex1 + 1
            verts.append( _getRoofVert(verts[firstVertIndex + minDistanceIndex2], minDistanceIndex2, tan) )
            
            # variable below are used for the assignment of the UV-coordinates
            if self.setUvs:
                u1 = distances[minDistanceIndex1] * (1. + cos[minDistanceIndex1]) / sin[minDistanceIndex1]
                v1 = distances[minDistanceIndex1] * factor
                
                u2 = distances[minDistanceIndex2] * (1. + cos[minDistanceIndex2]) / sin[minDistanceIndex2]
                v2 = distances[minDistanceIndex2] * factor
            
            # Triangle of the hipped roof originating from the polygon edge
            # with the index <minDistanceIndex1>
            roofItem.addRoofSide(
                (firstVertIndex + minDistanceIndex1, firstVertIndex + minDistanceIndex1Next, vertIndex1),
                ( (0., 0.), (lengths[minDistanceIndex1], 0.), (u1, v1) ) if self.setUvs else None,
                minDistanceIndex1
            )
            
            # Quadrangle of the hipped roof originating from the polygon edge
            # with the index <minDistanceIndex1Next>
            roofItem.addRoofSide(
                (firstVertIndex + _nextIndices[minDistanceIndex1], firstVertIndex + minDistanceIndex2, vertIndex2, vertIndex1),
                (
                    (0., 0.),
                    (lengths[minDistanceIndex1Next], 0.),
                    (lengths[minDistanceIndex1Next] - u2, v2),
                    (lengths[minDistanceIndex1] - u1, v1)
                ) if self.setUvs else None,
                minDistanceIndex1Next
            )
            
            # Triangle of the hipped roof originating from the polygon edge
            # with the index <minDistanceIndex2>
            roofItem.addRoofSide(
                (firstVertIndex + minDistanceIndex2, firstVertIndex + minDistanceIndex2Next, vertIndex2),
                ( (0., 0.), (lengths[minDistanceIndex2], 0.), (u2, v2) ) if self.setUvs else None,
                minDistanceIndex2
            )
            
            # Quadrangle of the hipped roof originating from the polygon edge
            # with the index <minDistanceIndex2Next>
            roofItem.addRoofSide(
                (firstVertIndex + minDistanceIndex2Next, firstVertIndex + minDistanceIndex1, vertIndex1, vertIndex2),
                (
                    (0., 0.),
                    (lengths[minDistanceIndex2Next], 0.),
                    (lengths[minDistanceIndex2Next] - u1, v1),
                    (lengths[minDistanceIndex2] - u2, v2)
                ) if self.setUvs else None,
                minDistanceIndex2Next
            )
        return True
    
    def generateRoof(self, footprint, roofItem, firstVertIndex):
        verts = footprint.building.element.l.genVolumes.verts
        numPolygonVerts = footprint.polygon.numEdges

        roofSideIndices = []
        
        unitVectors = [
            vector.unitVector3d for vector in footprint.polygon.getVectors()
        ]

        lengths = [
            vector.length for vector in footprint.polygon.getVectors()
        ]
        
        #dumpInputHippedRoof(verts, firstVertIndex, numPolygonVerts, None, unitVector)
        #return
        
        # calculate polygons formed by the straight skeleton
        polygonize(
            verts,
            firstVertIndex,
            numPolygonVerts,
            None,
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
            edgeIndex = indices[0] - firstVertIndex
            if edgeIndex < numPolygonVerts:
                # The normal case:
                # all faces of the roof have a common edge with the original footprint
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
            else:
                # A special exotic case:
                # the face of the roof doesn't have a common edge with the original footprint
                origin = verts[indices[0]]
                # 1) get a normal to the face
                n = (verts[indices[1]] - origin).cross(verts[indices[2]]-verts[indices[1]])
                n.normalize()
                # 2) get the unit vector for the U-axis located in a horizontal plane
                u = zAxis.cross(n)
                u.normalize()
                # 3) get the unit vector for the X-axis
                v = n.cross(u)
                roofItem.addRoofSide(
                    indices,
                    # UV-coordinates
                    tuple(
                        (
                            (verts[indices[_index]] - origin).dot(u),
                            (verts[indices[_index]] - origin).dot(v)
                        ) for _index in range(len(indices))
                    ),
                    -1
                )
        return True