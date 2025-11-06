from util.blender import getBmesh, setBmesh


class MeshGen:

    def __init__(self, obj):
        self.obj = obj
        self.bm = getBmesh(obj)
        self.verts = []
        # counterparts for <self.verts> in the BMesh
        self.bmVerts = []


class MeshGenByIndices(MeshGen):

    def createFace(self, footprint, indices):
        bm = self.bm
        renderInfo = footprint.building.renderInfo
        verts = self.verts
        bmVerts = self.bmVerts
        
        # extend <bmVerts> to have the same number of vertices as in <verts>
        bmVerts.extend(None for _ in range(len(verts)-len(bmVerts)))
        
        # check if we have BMVerts for for all <indices>
        for index in indices:
            if not bmVerts[index]:
                bmVerts[index] = bm.verts.new(
                    (verts[index] + renderInfo.offsetVertex) if renderInfo.offsetVertex else verts[index]
                )
        
        return bm.faces.new(bmVerts[index] for index in indices)

    def createFootprint(self, footprint):
        face = self.bm.faces.new(
            self.bm.verts.new(vector.v1_3d) for vector in footprint.bldgPart.polygon.vectors
        )

        # set the attribute "roof_shape" to <face>
        # 1: flat
        # 4: skillion
        bmLayer = self.bm.faces.layers.int["roof_shape"]
        face[bmLayer] = 4 if footprint.getStyleBlockAttr("roofShape") == "skillion" else 1

        # set the attribute "facade_class" to each edge of <face>
        bmLayer = self.bm.edges.layers.int["facade_class"]
        loop = face.loops[0]
        for facade in footprint.facades:
            loop.edge[bmLayer] = facade.cl
            loop = loop.link_loop_next

    def createFootprintsForHoles(self, footprint):
        for innerPolygon in footprint.innerPolygons:
            face = self.bm.faces.new(
                self.bm.verts.new(vector.v1_3d) for vector in innerPolygon.vectors
            )
    
    def finalize(self):
        setBmesh(self.obj, self.bm)


class MeshGenByCoords(MeshGen):
    pass