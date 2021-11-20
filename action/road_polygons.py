import numpy as np
from math import sqrt

import matplotlib.pyplot as plt

from way.way_network import WayNetwork, NetSegment
from way.way_algorithms import createSectionNetwork
from defs.way import allRoadwayCategoriesRank
from mathutils import Vector
from lib.CompGeom.polygon import Polygon
from lib.CompGeom.properties import isClockwise, pointInPolygon
from lib.CompGeom.algorithms import circumCircle, SCClipper
from lib.CompGeom.splitPolygonHoles import splitPolygonHoles


def polylineBisectors(polyline, dist=1.):
    bisectors = []
    # unit vectors
    v = [(p2-p1).normalized() for p1,p2 in zip(polyline[:-1],polyline[1:])]
    v = [v[0]]+v+[v[-1]]
    for v1,v2 in zip(v[:-1],v[1:]):
        bisect = (v1-v2).normalized()
        if bisect == Vector((0,0)):
            bisectors.append( Vector((-v1[1],v1[0]))*dist )
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
        vertList = [Vector((firstVert[0],firstVert[1]))]
        nextVert = segList.pop(firstVert, None).freeze()
        while nextVert != firstVert and nextVert is not None:
            vertList.append(Vector((nextVert[0],nextVert[1])))
            nextVert = segList.pop(nextVert, None).freeze()
        polyList.append( Polygon.fromPointlist(vertList) )
    return polyList

class RoadPolygons:

    def __init__(self):
        self.allNetwork = None
        self.sectionNetwork = None
        self.polylines = None
        self.polyVerts = None
        self.vertIndexToPolyIndex = []
        self.kdTree = None
        self.cycles = None
        self.searchHeight = 80.
        self.searchHeight_2 = self.searchHeight*self.searchHeight

    def do(self, manager):
        self.createWaySectionNetwork()
        self.createPolylinePolygons(manager)
        self.cycles = self.sectionNetwork.iterCycles()
        self.createGraphCyclePolygons()

    def cleanup(self):
        self.kdTree = None
        self.bldgVerts = None
        self.vertIndexToPolyIndex.clear()

    def createPolylinePolygons(self,manager):
        # create list <polyList> of polygon objects from polylines. Currently, only closed polylines
        # are accepted. Later, lines like fences have to get an area, for instnce a width of 0.5m
        self.polyList = []
        for polyline in manager.polylines:
            if polyline.element.closed:
                pointList = [Vector((edge.v1[0],edge.v1[1])) for edge in polyline.edges]
                self.polyList.append(Polygon.fromPointlist(pointList))

        # add building polylines as polygon objects to <polyList>
        self.buildingPolynoms = mergeBuildingsToBlocks(manager.buildings)
        self.polyList.extend(self.buildingPolynoms)

        # assert that all polyline objects are oriented clockwise
        for poly in self.polyList:
            if not poly.isClockwise():
                poly.points.reverse()

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
            for vert in polygon:
                self.polyVerts[index] = (vert[0],vert[1])
                index += 1

        self.createKdTree()

    def createGraphCyclePolygons(self):

        processedCycles = []
        # plotNetwork(self.sectionNetwork)
        for cycle in self.cycles:
            # if any(seg.category=='scene_border' for seg in cycle):
            #     continue # cycles at scene order require special handling

            cyclePoly = Polygon()
            points = []
            for c in cycle:
                for p in c.path[:1]:                   
                    indx = next((i for i, v in enumerate(points) if np.all(np.isclose(v,p)) ), -1)
                    if indx >= 0:
                        del points[indx:] # remove spikes by dead-ends.
                    points.append(p)

            # points = [p for s in cycle for p in s.path[:-1]]
            cyclePoly.addPoints( points )

            # only the polygon on the out-side of the scene can be oriented clockwise
            if cyclePoly.isClockwise():
                continue

            # plotPoly(cyclePoly)
            # plt.title('CW' if cyclePoly.isClockwise() else 'CCW')
            # plotEnd()
            # continue

            center, radius = circumCircle(points)

            # Get all polygon vertices for the polygons whose vertices
            # are returned by the query to the KD tree
            queryIndices = set(
                self.makeKdQuery(center, radius)
            )

            # # keep only those that are in the cycle's polygon
            # queryCycleIndicesTmp = [self.vertIndexToPolyIndex[indx] for indx in queryIndices if \
            #        pointInPolygon(cyclePoly.points,Vector((self.polyVerts[indx])))==2]

            # for indx in queryCycleIndicesTmp:
            #     v = self.polyVerts[self.vertIndexToPolyIndex[indx]]
            #     plt.plot(v[0],v[1], 'kx',zorder=600)


            queryCycleIndices = set(self.vertIndexToPolyIndex[indx] for indx in queryIndices)

            debug = True
            if queryCycleIndices:
                if debug:
                    ax1 = plt.subplot(1,2,1)
                # plotNetwork(self.sectionNetwork)
                plotPoly(cyclePoly,'b',2,100)
                for indx in queryCycleIndices:
                    objPoly = self.polyList[indx]
                    plotPoly(objPoly,'g:',1,100)
                plt.gca().axis('equal')
                if not debug:
                    plotEnd()

            origCyyclePoly = cyclePoly

            if queryCycleIndices:
                possibleHoles = []
                for indx in queryCycleIndices:
                    objPoly = self.polyList[indx]

                    if not objPoly.isClockwise():
                        print('objPoly was counter-clockwise')
                        objPoly.points.reverse()
                        print('checked:',objPoly.isClockwise())

                    if cyclePoly.isClockwise():
                        print('cyclePoly was clockwise')
                        cyclePoly.points.reverse()
                        print('checked:',not cyclePoly.isClockwise())

                    wasException = False
                    try:
                        differences = Polygon.subtract(cyclePoly,objPoly)
                    except Exception as e:
                        import sys, os, traceback
                        print(traceback.format_exc())
                        plt.close()
                        plotNetwork(self.sectionNetwork)                        
                        plotPoly(cyclePoly,'g',1,100)
                        plotPoly(origCyyclePoly,'b',1,100)
                        plotPoly(objPoly,'r',1,100)
                        plt.title('EXCEPTION')
                        plotEnd()
                        wasException = True
                        break

                    if differences:
                        if  differences[0].points:
                            cyclePoly = differences[0]
                            for possibleHole in differences[1:]:
                                if possibleHole and possibleHole.isClockwise():
                                    possibleHoles.append(possibleHole)
                        else:
                            break
                    if cyclePoly.isClockwise():
                        cyclePoly.points.reverse()

                if possibleHoles:
                    for hole in possibleHoles:
                        pointsInCycle = [pointInPolygon(cyclePoly.points,p) for p in hole.points]
                        if not all(pointsInCycle):
                            possibleHoles.remove(hole)
                            print('hole removed')
                if possibleHoles:
                    splittedPoly = splitPolygonHoles(cyclePoly,possibleHoles)
                    cyclePoly.clear()
                    cyclePoly.addPoints(splittedPoly)


                if not wasException:
                    if debug:
                        ax2 = plt.subplot(1,2,2)
                    plotPoly(origCyyclePoly,'b:',1,100)
                    x = [n[0] for n in cyclePoly.points]
                    y = [n[1] for n in cyclePoly.points]
                    plt.fill(x,y,'#0000ff',alpha = 0.1,zorder = 100)
                    plotPoly(cyclePoly,'r',2,200)
                    plotEnd()


    def createKdTree(self):
        from scipy.spatial import KDTree
        self.kdTree = KDTree(self.polyVerts)

    def makeKdQuery(self, searchCenter, searchRadius):
        return (
            vertIndex for vertIndex in self.kdTree.query_ball_point(searchCenter, searchRadius)
        )

    def createWaySectionNetwork(self):
        # get border polygon (ounter-clockwise) of scene frame
        minX, minY = self.app.projection.fromGeographic(self.app.minLat, self.app.minLon)
        maxX, maxY = self.app.projection.fromGeographic(self.app.maxLat, self.app.maxLon)

        # prepare clipper for this frame
        clipper = SCClipper(minX,maxX,minY,maxY)

        wayManager = self.app.managersById["ways"]

         # create full way network
        self.allNetwork = WayNetwork()
        excludedTags = ['tunnel','steps','covered']
        for way in wayManager.getAllWays():
            if [tag for tag in excludedTags if tag in way.element.tags]:
                continue
            for segment in way.segments:
                way.element.tags
                v1, v2 = Vector(segment.v1),Vector(segment.v2)
                accepted, v1, v2 = clipper.clip(v1,v2)
                if accepted:
                    netSeg = NetSegment(v1,v2,way.category,(v2-v1).length)
                    self.allNetwork.addSegment(netSeg)
        borderPolygon = clipper.getPolygon()
        for v1,v2 in zip(borderPolygon[:-1],borderPolygon[1:]):
            netSeg = NetSegment(v1,v2,'scene_border',(v2-v1).length)
            self.allNetwork.addSegment(netSeg)

        # create way-section network
        self.sectionNetwork = createSectionNetwork(self.allNetwork)

# plotting functions for development
#-----------------------------------------------------------------------------------
def plotPoly(polygon,color='k',width=1.,order=100):
    count = 0
    for v1,v2 in zip(polygon[:-1],polygon[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
        # plt.text(v1[0],v1[1],str(count))
        count += 1
        # plt.plot(v1[0],v1[1],'kx')
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

def plotNetwork(network):
    for seg in network.iterAllSegments():
        plotWaySeg(seg,'k',0.5)


def plotEnd(xlim=None,ylim=None):
    plt.gca().axis('equal')
    if xlim and ylim:
        plt.gca().set_xlim(xlim)
        plt.gca().set_ylim(ylim)
    # plt.savefig("C:/Users/Roger/Desktop/test.svg")
    plt.show()
