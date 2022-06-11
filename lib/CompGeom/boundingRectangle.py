from mathutils import Vector, Matrix
from math import inf

from lib.pygeos.shared import Coordinate

def boundingRectangle(verts):
     
    if isinstance(verts[0],Coordinate):
        points = [Vector((v.x,v.y)) for v in verts]
    else:
        points = [Vector((v[0],v[1])) for v in verts]
    if points[0] == points[-1]:
        points.pop()
    points.sort()

    # Compute convex hull by monotone chain algorithm
    # first lower hull 
    lower = []
    for p in points:
        while len(lower) >= 2 and (lower[-1]-lower[-2]).cross(p-lower[-2]) <= 0:
            lower.pop()
        lower.append(p)

    # then upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and (upper[-1]-upper[-2]).cross(p-upper[-2]) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list. 
    hull = lower[:-1] + upper[:-1]

    # unit vectors of hull edges
    edges = [(v1-v0).normalized() for v0,v1 in zip(hull, hull[1:] + hull[:1]) ]

    # start rotating caliper
    minArea = [inf, Vector((0.,0.)), 0., 0., 0., 0.]
    for e in edges:
        fwd = Matrix( [ (e[0], -e[1]), (e[1], e[0]) ] )
        rotV = [ (v @ fwd) for v in hull ]

        vertsX = [ p[0] for p in rotV ]
        vertsY = [ p[1] for p in rotV ]
        area = (max(vertsX)-min(vertsX)) * (max(vertsY)-min(vertsY))
        if area < minArea[0]:
            minArea = [area, e, min(vertsX), max(vertsX), min(vertsY), max(vertsY)]

    a = abs(minArea[2]-minArea[3])
    b = abs(minArea[4]-minArea[5])
    length,width = (a,b) if a > b else (b,a)

    # construct found minimum bounding rectangle        
    minRect = [
        Vector((minArea[2],minArea[4])),
        Vector((minArea[3],minArea[4])),
        Vector((minArea[3],minArea[5])),
        Vector((minArea[2],minArea[5])),
    ]

    # rotate it back
    back = Matrix( [ (minArea[1][0], minArea[1][1]), (-minArea[1][1], minArea[1][0]) ] )
    rect = [ (v @ back) for v in minRect ]

    return length, width, rect
