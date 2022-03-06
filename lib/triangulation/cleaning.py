from mathutils import Vector
from collections import defaultdict
import matplotlib.pyplot as plt
from math import atan2, pi

from .Vertex import Vertex

def removeStraightAngleVerts(poly):
    # poly is a polygon formed by vertices of type mathutils.Vector.
    vectors = [v2-v1 for v1,v2 in zip([poly[-1]]+poly[:-1],poly)]
    angles = [abs(atan2(v1.cross(v2),v1.dot(v2))) for v1,v2 in zip(vectors,vectors[1:]+[vectors[0]])]
    cleanPoly = [point for point,angle in zip(poly,angles) if angle > 3*pi/180]
    return cleanPoly


def cleaningForTriangulation(geosPoly):
    # get polygons from pyGEOS Polygon
    poly = [Vector((v.x,v.y)).freeze() for v in geosPoly.exterior.coords[:-1]]
    if not geosPoly.exterior.is_ccw:    # polygon is ordered counter-clockwise
        poly.reverse()
    holes = []
    for geosHole in geosPoly.interiors:
        hole = [Vector((v.x,v.y)).freeze() for v in geosHole.coords[:-1]]
        if geosHole.is_ccw:             # hole is ordered clockwise
            hole.reverse()
        holes.append(hole)

    if holes:
        # check for every hole
        for hole in holes:
            sharedVertices = set(poly).intersection(set(hole))
            if sharedVertices:
                # We correct the shared vertex of the hole by moving it 1mm along
                # the bisector into the hole.
                for shared in sharedVertices:
                    # find the shared vertex and its neighbors in hole 
                    indx =  hole.index(shared)
                    p1 = hole[indx-1]
                    p2 = hole[(indx+1)%len(hole)]
                    # unity vectors of the edges
                    u1 = (shared-p1)/(shared-p1).length
                    u2 = (p2-shared)/(p2-shared).length
                    # move this point 1mm along bisector
                    if u1.cross(u2) < 0:
                        u1 = -u1
                    else:
                        u2 = -u2
                    bisector = (u1 + u2)/(u1+u2).length
                    hole[indx] = shared + bisector * 0.001

    # cleanPoly = removeStraightAngleVerts(poly)
    polyV = [Vertex(v) for v in poly]
    holesV = []
    for hole in holes:
        # cleanHole = removeStraightAngleVerts(hole)
        holeV = [Vertex(v) for v in hole]
        holesV.append(holeV)
    return polyV, holesV

def plotPoly(polygon,vertsOrder,color='k',width=1.,order=100):
    count = 0
    for v1,v2 in zip(polygon[:-1],polygon[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
        if vertsOrder:
            plt.text(v1[0],v1[1],str(count),fontsize=12)
        count += 1
        # plt.plot(v1[0],v1[1],'kx')
    v1, v2 = polygon[-1], polygon[0]
    plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
    if vertsOrder:
        plt.text(v1[0],v1[1],str(count),fontsize=12)

def plotEnd():
    plt.gca().axis('equal')
    plt.show()

    # polygon1 = [Vector((v.x,v.y)).freeze() for v in poly.exterior.coords[:-1]]
    # if not poly.exterior.is_ccw:    # polygon is ordered counter-clockwise
    #     polygon1.reverse()
    # holes1 = []
    # for geosHole in poly.interiors:
    #     hole = [Vector((v.x,v.y)).freeze() for v in geosHole.coords[:-1]]
    #     if geosHole.is_ccw:             # hole is ordered clockwise
    #         hole.reverse()
    #     holes.append(hole)
