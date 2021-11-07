import numpy as np
from math import sqrt

import matplotlib.pyplot as plt

from way.way_network import WayNetwork, NetSegment
from way.way_algorithms import createSectionNetwork
from defs.way import allRoadwayCategoriesRank
from mathutils import Vector as mutlsVector
from lib.py2D.Polygon import Polygon
from lib.py2D.Vector import Vector


def polylineBisectors(polyline, dist=1.):
    bisectors = []
    # unit vectors
    v = [(p2-p1).normalized() for p1,p2 in zip(polyline[:-1],polyline[1:])]
    v = [v[0]]+v+[v[-1]]
    for v1,v2 in zip(v[:-1],v[1:]):
        bisect = (v1-v2).normalized()
        if bisect == mutlsVector((0,0)):
            bisectors.append( mutlsVector((-v1[1],v1[0]))*dist )
        else:
            bisectors.append( bisect*dist if v1.cross(v2)<0. else -bisect*dist )
    return bisectors

def polygonBisectors(polygon, dist=1.):
    bisectors = []
    # unit vectors
    v = [(p2-p1).normalized() for p1,p2 in zip(polygon, polygon[1:] + polygon[:1])]
    for v1,v2 in zip(v, v[1:] + v[:1]):
        bisect = (v1-v2).normalized()
        bisectors.append( bisect*plt.disconnect if v1.cross(v2)<0. else -bisect*dist )
    return bisectors

def mergeBuildingsToBlocks(buildings):
    segList = dict()
    for building in buildings:
        verts = [v for v in building.polygon.verts]
        for _, v1, v2 in building.polygon.edgeInfo(verts, 0, skipShared=True):
            segList[v1.freeze()] = v2

    polyList = []
    while segList:
        firstVert = next(iter(segList))
        vertList = [Vector(firstVert[0],firstVert[1])]
        nextVert = segList.pop(firstVert, None).freeze()
        while nextVert != firstVert and nextVert is not None:
            vertList.append(Vector(nextVert[0],nextVert[1]))
            nextVert = segList.pop(nextVert, None).freeze()
        polyList.append( Polygon.from_pointlist(vertList) )
    return polyList

class RoadPolygons:

    def __init__(self):
        self.allNetwork = None
        self.sectionNetwork = None
        self.polylines = None
        self.polyVerts = None
        self.vertIndexToPolyIndex = []
        self.kdTree = None
        self.searchHeight = 80.
        self.searchHeight_2 = self.searchHeight*self.searchHeight

    def do(self, manager):
        self.createAllWaySectionNetwork()
        self.createPolygons(manager)
        self.createRoads()

    def cleanup(self):
        self.kdTree = None
        self.bldgVerts = None
        self.vertIndexToPolyIndex.clear()

    def createPolygons(self,manager):
        # create list <polyList> of polygon objects from polylines. Currently, only closed polylines
        # are accepted. Later, lines like fences have to get an area, for instnce a width of 0.5m
        # test = [edge.v1 for polyline in manager.polylines for edge in polyline.edges if polyline.element.closed]
        # pass
        self.polyList = []
        for polyline in manager.polylines:
            if polyline.element.closed:
                pointList = [Vector(edge.v1[0],edge.v1[1]) for edge in polyline.edges]
                self.polyList.append( Polygon.from_pointlist(pointList))

        # add building polylines as polygon objects to <polyList>
        buildingPolynoms = mergeBuildingsToBlocks(manager.buildings)
        self.polyList.extend(buildingPolynoms)

         # create mapping between the index of the vertex and index of the polygon in <polyList>
        self.vertIndexToPolyIndex.extend(
            polyIndex for polyIndex, polygon in enumerate(self.polyList) for _ in range(len(polygon))
        )

        # the total number of vertices
        totalNumVerts = sum(len(polygon) for polygon in self.polyList )

        # allocate the memory for an empty numpy array
        self.polyVerts = np.zeros((totalNumVerts, 2))

        # fill vertices in <polyVerts>
        index = 0
        for polygon in self.polyList:
            # store the index of the first vertex of <polygon> in <self.polyVerts>
            # polygon.auxIndx = index
            for vert in polygon:
                self.polyVerts[index] = (vert[0],vert[1])
                index += 1

        self.createKdTree()
        pass

    def createRoads(self):
        searchHeight = 80.      # to be defined later in defs, !!!! needs to be adapted to way category !!!!

        for waySeg in self.sectionNetwork.iterAllSegments():
            if allRoadwayCategoriesRank[waySeg.category] > allRoadwayCategoriesRank['secondary']:
                continue

            # c = (waySeg.s+waySeg.t)/2.
            # plt.text(c[0],c[1],str(waySeg.ID))
            if waySeg.ID not in [3788]:
                continue

            # create polygon (counter-clockwise) of search range for the whole segment
            bisectors = polylineBisectors(waySeg.path,searchHeight)
            searchPoly = Polygon()
            for v,b in zip(waySeg.path,bisectors):
                p = v-b
                searchPoly.add_point(Vector(p[0],p[1]))
            for v,b in zip(waySeg.path[::-1],bisectors[::-1]):
                p = v+b
                searchPoly.add_point(Vector(p[0],p[1]))

            plotWaySeg(waySeg,'b',2.,101)

            processedPolys = set()
            holes = []
            for v1,v2 in zip(waySeg.path[:-1],waySeg.path[1:]):
                v1, v2 = np.array(v1), np.array(v2) # <sectionNetwork> delivers mathutils.Vector's
                way = (v2-v1)
                wayLength = np.linalg.norm(way)
                wayCenter= (v1 + v2)/2.
                searchWidth = wayLength/2.

                # Get all polygon vertices for the polygons whose vertices
                # are returned by the query to the KD tree
                queryPolyIndices = set(
                    self.makeKdQuery(wayCenter, sqrt(searchWidth*searchWidth + searchHeight*searchHeight))
                )

                queryPolyIndices -= processedPolys
                processedPolys = processedPolys | queryPolyIndices

                if queryPolyIndices:
                    for indx in queryPolyIndices:
                        try:
                            difference = Polygon.subtract(searchPoly,self.polyList[indx])
                        except:
                            # plotPoly(searchPoly,'b',1,100)
                            # plotPoly(self.polyList[indx],'r',1,100)
                            # plt.title('EXCEPTION')
                            # plotEnd()
                            break

                        N = len(difference)
                        if N == 1:  # object polygon outside <searchPoly>
                            searchPoly = difference[0]
                            continue
                        for s in difference:
                            if s: 
                                if s.is_clockwise() and s.contains_point(Vector(wayCenter[0],wayCenter[1])):
                                    # object polygon intersected <searchPoly>, remaining is new <searchPoly>
                                    searchPoly = s
                                elif not s.is_clockwise():
                                    # object polygon is a hole
                                    holes.append(s)
                                else:
                                    # residual polygon cutted from <searchPoly>
                                    pass

            plotResult = False
            if plotResult:
                plt.subplot(2,1,1)
                for h in holes:
                    plotPoly(h,'b',1,100)
                plotPoly(searchPoly,'r',1,100)
                plt.gca().axis('equal')
                plt.subplot(2,1,2)
                plotRange(searchPoly,holes)
                plotEnd()
            
            # plot to final image
            plotRange(searchPoly,holes)

    def createKdTree(self):
        from scipy.spatial import KDTree
        self.kdTree = KDTree(self.polyVerts)

    def makeKdQuery(self, searchCenter, searchRadius):
        return (
            self.vertIndexToPolyIndex[vertIndex] for vertIndex in self.kdTree.query_ball_point(searchCenter, searchRadius)
        )

    def createAllWaySectionNetwork(self):
        wayManager = self.app.managersById["ways"]

         # create full way network
        self.allNetwork = WayNetwork()
        for way in wayManager.getAllWays():
            if allRoadwayCategoriesRank[way.category] <= allRoadwayCategoriesRank['secondary']:
                for segment in way.segments:
                    netSeg = NetSegment(segment.v1,segment.v2,way.category,segment.length)
                    self.allNetwork.addSegment(netSeg)

        # create way-section network
        self.sectionNetwork = createSectionNetwork(self.allNetwork)

# plotting functions for development
#-----------------------------------------------------------------------------------
def plotPoly(polygon,color='k',width=1.,order=100):
    for v1,v2 in zip(polygon[:-1],polygon[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
    v1, v2 = polygon[-1], polygon[0]
    plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)

def plotWaySeg(wayseg,color='k',width=1.,order=100):
    for v1,v2 in zip(wayseg.path[:-1],wayseg.path[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)

def plotRange(poly,holes,color='#ff0000',alpha = 0.7,zorder=2):
    from lib.py2D.patchify import patchify
    polyList = []
    poly_np = np.array((
        tuple(p[0] for p in poly),
        tuple(p[1] for p in poly),
    ))
    polyList.append(poly_np)
    holes_np = []
    for hole in holes:
        hole_np = np.array((
            tuple(p[0] for p in hole),
            tuple(p[1] for p in hole),
        ))
        polyList.append(hole_np)
    patch = patchify(polyList,color,2,alpha,zorder)
    plt.gca().add_patch(patch)

def plotEnd():
    plt.gca().axis('equal')
    plt.show()
