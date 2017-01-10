import bmesh
from mathutils import Vector
from util.polygon import Polygon
from . import Roof


class RoofFlat(Roof):
    """
    A class to deal with building or building parts with the flat roof
    for a simple polygon (not a multipolygon)
    """
    
    defaultHeight = 0.01
    
    def make(self, bldgMaxHeight, roofMinHeight, bldgMinHeight, osm):
        polygon = self.polygon
        n = len(self.verts)
        
        # Create sides for the prism with the height <bldgMaxHeight - bldgMinHeight>,
        # that is based on the <polygon>
        polygon.sidesPrism(bldgMaxHeight, self.wallIndices)
        # vertices of the top part of the prism serve for the flat roof
        self.roofIndices.append( tuple(range(n, n+polygon.n)) )
        return True


class RoofFlatMulti(RoofFlat):
    """
    A class to deal with building or building parts with the flat roof
    for a multipolygon
    """
    
    def __init__(self):
        self.verts = []
        self.wallIndices = []
        self.polygons = []
    
    def init(self, element, minHeight, osm):
        self.verts.clear()
        self.wallIndices.clear()
        self.polygons.clear()
        
        self.element = element
    
    def make(self, bldgMaxHeight, roofMinHeight, bldgMinHeight, osm):
        if bldgMinHeight is None:
            bldgMinHeight = roofMinHeight
        element = self.element
        
        # create vertices for all linestrings of the multipolygon
        verts = self.verts
        polygons = self.polygons
        indexOffset = 0
        for _l in element.l:
            verts.extend(
                Vector((v[0], v[1], bldgMaxHeight)) for v in element.getLinestringData(_l, osm)
            )
            n = len(verts)
            polygons.append(
                Polygon(
                    tuple(range(indexOffset, n)),
                    verts
                )
            )
            indexOffset = n
        # vertices for the bottom
        verts.extend(Vector((verts[i].x, verts[i].y, bldgMinHeight)) for i in range(n))
        return True
    
    def render(self, r):
        element = self.element
        polygons = self.polygons
        bm = r.bm
        # create BMesh vertices
        verts = tuple( bm.verts.new(v) for v in self.verts )
        # create BMesh edges out of <verts>
        edges = tuple(
            bm.edges.new( (verts[polygon.indices[i-1]], verts[polygon.indices[i]]) )\
            for polygon in polygons\
            for i in range(polygon.n)
        )
        
        # a magic function that does everything
        geom = bmesh.ops.triangle_fill(bm, use_beauty=True, use_dissolve=True, edges=edges)
        materialIndex = r.getRoofMaterialIndex(element)
        # check the normal direction of the created faces and assign material to all BMFace
        for f in geom["geom"]:
            if isinstance(f, bmesh.types.BMFace):
                if f.normal.z < 0.:
                    f.normal_flip()
                f.material_index = materialIndex
        
        # create BMesh faces for building sides
        indexOffset1 = 0
        indexOffset2 = len(verts)//2
        wallIndices = self.wallIndices
        for polygon in polygons:
            n = polygon.n
            # the first edge of the polygon
            edge = edges[indexOffset1]
            if not edge.link_loops:
                # something wrong with the topology of the related OSM multipolygon
                # update index offsets to switch to the next polygon (i.e. a closed linestring)
                indexOffset1 += n
                indexOffset2 += n
                # skip that polygon
                continue
            # a BMLoop for <edge>
            l = edge.link_loops[0]
            # check if the direction of <polygon> needs to be reverted
            keepDirection = l.link_loop_next.vert == verts[indexOffset1]
            
            wallIndices.extend(
                (
                    indexOffset1 - 1 + (i if i else n),
                    indexOffset2 - 1 + (i if i else n),
                    indexOffset2 + i, 
                    indexOffset1 + i
                )\
                if keepDirection else\
                (
                    indexOffset1 - 1 + (i if i else n),
                    indexOffset1 + i,
                    indexOffset2 + i,
                    indexOffset2 - 1 + (i if i else n)
                )\
                for i in range(n)
            )
            
            indexOffset1 += n
            indexOffset2 += n
        
        materialIndex = r.getWallMaterialIndex(element)
        # actual code to create BMesh faces for the building walls out of <verts> and <wallIndices>
        for f in (bm.faces.new(verts[i] for i in indices) for indices in wallIndices):
            f.material_index = materialIndex