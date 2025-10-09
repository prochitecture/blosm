import bmesh
from mathutils import Vector
from . import ItemRendererTexture
from util import zAxis


def _getEdge(bmVert):
    return (
        bmVert.link_loops[0]\
        if bmVert.co[2] > bmVert.link_loops[1].link_loop_next.vert.co[2] else\
        bmVert.link_loops[1]
    ).edge


class RoofFlatMulti(ItemRendererTexture):
    
    def render(self, roofItem):
        # all needed BMesh verts have been already created during facade generation
        bmVerts = roofItem.footprint.element.l.genVolumes.bmVerts
        
        # create a Python list of BMesh edges for the outer and inner polygons
        indexOffset = roofItem.firstVertIndex
        # treat the outer polygon
        polygon = roofItem.footprint.polygon
        edges = [_getEdge(bmVerts[i]) for i in range(indexOffset, indexOffset + polygon.numEdges) ]
        
        # treat the inner polygons
        indexOffset += polygon.numEdges
        for polygon in roofItem.innerPolygons:
            # skipping the verts for the lower cap
            indexOffset += polygon.numEdges
            edges.extend(_getEdge(bmVerts[i]) for i in range(indexOffset, indexOffset + polygon.numEdges) )
            # skipping the verts for the upper cap
            indexOffset += polygon.numEdges
        
        # <bmesh.ops.triangle_fill(..)> a magic function that does everything
        self.renderCladding(
            roofItem,
            [
                face for face in bmesh.ops.triangle_fill(roofItem.footprint.element.l.genVolumes.bm, use_beauty=False, use_dissolve=False, edges=edges)\
                ["geom"] if isinstance(face, bmesh.types.BMFace)
            ],
            None
        )
    
    def setCladdingUvs(self, roofItem, faces, claddingTextureInfo, uvs):
        textureWidthM = claddingTextureInfo["textureWidthM"]
        textureHeightM = textureWidthM * claddingTextureInfo["textureSize"][1] / claddingTextureInfo["textureSize"][0]
        
        polygon = roofItem.footprint.polygon
        
        # Arrange the texture along the longest edge of <polygon>,
        # so the longest edges surves as u-axis for the texture
        bldgVector = polygon.getLongestVector()
        offset = bldgVector.v1
        uVec = bldgVector.unitVector
        # <vVec> is perpendicular to <uVec>, <uVec.cross(vVec)> is pointing up (Z+)
        vVec = Vector((-uVec[1], uVec[0]))
        
        for face in faces:
            self.r.setUvs(
                face,
                # a generator!
                (
                    (
                        (Vector((vert.co[0], vert.co[1]))-offset).dot(uVec)/textureWidthM,
                        (Vector((vert.co[0], vert.co[1]))-offset).dot(vVec)/textureHeightM
                    ) for vert in face.verts
                ),
                roofItem.footprint.element.l,
                roofItem.footprint.element.l.uvLayerNameCladding
            )
    
    def setMaterial(self, item, faces, materialId):
        for face in faces:
            self.r.setMaterial(item.building.element.l, face, materialId)