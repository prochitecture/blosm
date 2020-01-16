import math
from .. import ItemRenderer

from util import zAxis


# Generatrix for a dome roof.
# generatix for dome is circle equation in the parameterized form x=x(t), y=y(t)
def generatrix_dome(rows):
    profile=[]
    for j in range(rows):
        x = math.cos(j/rows*math.pi/2)
        z = math.sin(j/rows*math.pi/2)
        profile.append((x,z))
    profile.append((0., 1.))
    return profile

# Generatrix for an onion roof.
# There is no explicit formula for an onion, that's why we will just write it as a set of vertex cordinates. 
# Note, that there are actually more than 5 typical onion types even for orthodox churches only.
# Other forms will be added later.
# (r, z) 
generatrix_onion =(
    (1.0000,    0.0000),
    (1.2971,    0.0999),
    (1.2971,    0.2462),
    (1.1273,    0.3608),
    (0.6219,    0.4785),
    (0.2131,    0.5984),
    (0.1003,    0.7243),
    (0.0000,    1.0000)
)


class RoofGeneratrix(ItemRenderer):
    
    def __init__(self, generatrix):
        super().__init__(False)
        self.generatrix = generatrix
        # The variable below indicates if the last point of the generatrix is located at zero,
        # i.e. in the center of the underlying polygon
        self.hasCenter = not self.generatrix[-1][0]
    
    def render(self, roofItem):
        gen = self.generatrix
        footprint = roofItem.footprint
        polygon = footprint.polygon
        building = roofItem.building
        verts = building.verts
        indices = roofItem.indices
        
        roofHeight = footprint.roofHeight
        roofVerticalPosition = footprint.roofVerticalPosition
        
        center = polygon.centerBB(roofVerticalPosition)
        
        n = polygon.n
        numRows = len(self.generatrix)
        
        vertIndexOffset = len(verts)
        
        verts.extend(
            center + gen[gi][0]*(verts[indices[vi]]-center) + gen[gi][1]*roofHeight*zAxis\
            for gi in range(1, numRows-1 if self.hasCenter else numRows) for vi in range(n)
        )
        # Also create a vertex at the center if the last point of the generatrix is located at zero,
        # i.e. in the center of the underlying polygon
        if self.hasCenter:
            verts.append(center + gen[-1][1]*roofHeight*zAxis)
                
        # the first row
        for vi in range(n-1):
            self.createFace(
                building,
                (indices[vi], indices[vi+1], vertIndexOffset+vi+1, vertIndexOffset+vi)
            )
        # and the closing quad for the ring
        self.createFace(building, (indices[-1], indices[0], vertIndexOffset, vertIndexOffset+n-1))
        
        # The rest of rows except the last row made of triangles ending at the center of
        # the underlying triangle
        for gi in range(1, numRows-2 if self.hasCenter else numRows-1):
            for vi in range(vertIndexOffset, vertIndexOffset+n-1):
                self.createFace(building, (vi, vi+1, vi+n+1, vi+n))
            # and the closing quad for the ring
            self.createFace(
                building,
                (vertIndexOffset+n-1, vertIndexOffset, vertIndexOffset+n, vertIndexOffset+2*n-1)
            )
            vertIndexOffset += n
        
        if self.hasCenter:
            # The last row made of triangles ending at the center of the underlying triangle
            for vi in range(vertIndexOffset, vertIndexOffset+n-1):
                self.createFace(building, (vi, vi+1, -1))
            # and the closing triangle for the ring
            self.createFace(building, (vertIndexOffset+n-1, vertIndexOffset, -1))
    
    def createFace(self, building, indices):
        face = self.r.createFace(building, indices)
        face.smooth = True