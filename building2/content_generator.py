from util.blender import getBmesh, setBmesh


def _getCachedBMVert(vector, data, bm):
    node = data.nodes[vector.id1]
    if node.cache:
        return node.cache
    else:
        # check if <vector.edge> is shared by two building footprints
        if len(vector.edge.vectors) == 2:
            # cache BMVerts for both ends of <vector>
            if not node.cache:
                # cache BMVert for the 1st vertex
                node.cache = bm.verts.new(vector.v1_3d)
            if not data.nodes[vector.id2].cache:
                # cache BMVert for the 2nd vertex
                data.nodes[vector.id2].cache = bm.verts.new(vector.v2_3d)
            return node.cache
        else:
            return bm.verts.new(vector.v1_3d)


class MeshGen:

    def __init__(self, obj=None):
        self.verts = []
        if obj:
            self.init(obj)
    
    def init(self, obj):
        self.obj = obj
        self.bm = getBmesh(obj)
        # counterparts for <self.verts> in the BMesh
        self.bmVerts = []


class MeshGenByIndices(MeshGen):

    def initFootprintAttributes(self):
        self.bm.faces.layers.int.new("roof_shape")
        self.bm.edges.layers.int.new("facade_class")

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

    def createFootprint(self, footprint, data, shareBMVerts):
        face =\
            self.bm.faces.new(
                _getCachedBMVert(vector, data, self.bm) for vector in footprint.bldgPart.polygon.getVectors()
            )\
            if shareBMVerts else\
            (
                self.bm.faces.new(
                    self.bm.verts.new(vector.v1_3d + footprint.building.renderInfo.offsetVertex)\
                        for vector in footprint.bldgPart.polygon.getVectors()
                )\
                if footprint.building.renderInfo.offsetVertex else\
                self.bm.faces.new(
                    self.bm.verts.new(vector.v1_3d) for vector in footprint.bldgPart.polygon.getVectors()
                )
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
        if footprint.building.renderInfo.offsetVertex:
            for innerPolygon in footprint.innerPolygons:
                self.bm.faces.new(
                    self.bm.verts.new(vector.v1_3d + footprint.building.renderInfo.offsetVertex)\
                        for vector in innerPolygon.getVectors()
                )
        else:
            for innerPolygon in footprint.innerPolygons:
                self.bm.faces.new(
                    self.bm.verts.new(vector.v1_3d) for vector in innerPolygon.getVectors()
                )
    
    def finalize(self):
        setBmesh(self.obj, self.bm)


class MeshGenByCoords(MeshGen):
    pass