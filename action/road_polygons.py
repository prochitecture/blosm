import numpy as np
from math import sqrt
from bisect import bisect_left

import matplotlib.pyplot as plt

from way.way_network import WayNetwork, NetSegment
from way.way_algorithms import createSectionNetwork

class PolyLine():
    def __init__(self,polyline):
        self.tags = polyline.element.tags
        self.closed = polyline.element.closed
        self.edges = [edge for edge in polyline.edges]
        if self.closed:
            self.verts = [edge.v1 for edge in polyline.edges]
            self.numVerts = len(self.verts)
            self.forceCcwDirection()
        else:
            self.verts = [edge.v1 for edge in polyline.edges] + [self.edges[-1].v2]
            # can't define counter-clockwise direction
            self.numVerts = len(self.verts)
        self.auxIndx = 0

    def forceCcwDirection(self):
        # Get index of the vertex with the minimum Y-coordinate and maximum X-coordinate,
        # i.e. the index of the rightmost lowest vertex
        index = min(
            range(self.numVerts),
            key = lambda index: ( self.verts[1], -self.verts[0] )
        )
        # <vectorIn>: the vector entering the vertex with the <index>
        vectorIn = self.verts[index]-self.verts[index-1]
        # <vectorOut>: the vector leaving the vertex with the <index>
        vectorOut = self.verts[(index+1)%self.numVerts]-self.verts[index]
        # Check if the <vectorOut> is to the left from the vector <vectorIn>;
        # in that case the direction of vertices is counterclockwise,
        # it's clockwise in the opposite case.
        if vectorIn[0] * vectorOut[1] - vectorIn[1] * vectorOut[0] < 0.:
            self.verts.reverse()

    def edgeInfo(self, queryPolyVerts, firstVertIndex):
        edges = iter(self.edges)
        verts = queryPolyVerts[firstVertIndex:firstVertIndex+self.numVerts]
        test=1
        for v1,v2 in zip(verts[:-1],verts[1:]):
            yield v1, v2

        if self.closed:
            yield verts[-1], verts[0]


class Event():
    def __init__(self, relStart, v1, v2, positiveY):
        if v1[0] > v2[0]:
            v1,v2 = v2,v1
        self.v1 = v1
        self.v2 = v2
        dv = v2-v1
        self.slope = dv[1]/dv[0] if dv[0] else None
        # None if start event, else start event that corresponds
        # to this end event
        self.relStart = relStart

        # values used for comparison, depending on positiveY
        self.sX = v1[0]
        self.sY = v1[1] if positiveY else -v1[1]    
        self.eX = v2[0]
        self.eY = v2[1] if positiveY else -v2[1]
        self.key = (self.eX,self.eY) if relStart else (self.sX,self.sY)  

    def yAt(self, x):
        return self.v1[1] + (x-self.v1[0])*self.slope if self.slope is not None else self.v2[1]

    def isInFrontOf(self,activeEvent):
        isInFront = False
        if self.sX == activeEvent.sX: # when both start at same x, y at end decides
            return  self.sY < activeEvent.eY
        if self.sY < max(activeEvent.sY,activeEvent.eY):
            # if this event starts behind the active event's minimum y
            if self.sY > min(activeEvent.sY,activeEvent.eY):
                # the decision has to be based on the true y-coord at this start x
                activeY = -activeEvent.yAt(self.sX) if activeEvent.relStart else activeEvent.yAt(self.sX)
                return self.sY < activeY
            else:
                # else, the new edge is in front for sure
                return True

    def __lt__(self,other):
        return self.key < other.key

class PriorityQueue():
    def __init__(self):
        self.keys = []
        self.queue = []
        self.ids = []

    def push(self, key, item):
        i = bisect_left(self.keys, key) # Determine where to insert item.
        self.keys.insert(i, key)        # Insert key of item to keys list.
        self.queue.insert(i, item)      # Insert the item itself in the corresponding place.
        self.ids.insert(i, id(item))

    def pop(self):
        self.keys.pop(0)
        self.ids.pop(0)
        return self.queue.pop(0)

    def remove(self, item):
        itemIndex = self.ids.index(id(item))
        del self.ids[itemIndex]
        del self.queue[itemIndex]
        del self.keys[itemIndex]

    def empty(self):
        return not self.queue

    def cleanup(self):
        self.keys.clear()
        self.queue.clear()
        self.ids.clear()

class CyrusBeckClipper():
    # polygon in expected as numpy array with dimensions (N,2).
    # it must consist of N > 2 2D vectors in counter-clockwise order
    def __init__(self, polygon):
        self.polyVertices = polygon
        self.eV = np.roll(polygon,-1,axis=0) - polygon  # polygon vectors
        self.eN = self.eV[:,[1,0]]*np.array((1,-1))     # vetor normal to outside of polygon

    def clipSegment(self, p1, p2):
        assert any(p1 != p2), 'degenerated segment'
        dS = p2 - p1
        tE, tL = 0., 1.
        for eNi, Pi in zip(self.eN,self.polyVertices):
            dotP = np.dot(eNi, dS)
            if dotP:
                t = np.dot(eNi, p1 - Pi) / (-dotP)
                tE, tL = (max(tE, t), tL) if dotP < 0 else (tE, min(tL, t))
        if tE <= tL:
            return p1 + dS*tE, p1 + dS*tL
        else:
            return None, None

class Skyline():
    def __init__(self,positiveY,width,height):
        self.positiveY = positiveY
        self.events = []
        self.skyline = []

        # add far baseline
        h = height if positiveY else -height
        s1 = np.array((-width,h))
        s2 = np.array((width,h))
        self.events.append( Event(None,s1,s2,self.positiveY))
        self.events.append( Event(self.events[-1],s1,s2,self.positiveY))

        # prepare clipper
        clipPolygon = np.array([
            (-width,-height),
            ( width,-height),
            ( width, height),
            (-width, height)
        ])
        self.clipper = CyrusBeckClipper(clipPolygon)

    def addSegment(self,v1,v2):
        if v1[0] != v2[0]: # avoid segments that are perpendicular to way segment
            p1, p2 = self.clipper.clipSegment(v1, v2)
            if p1 is not None:
                self.events.append( Event(None,p1,p2,self.positiveY))
                self.events.append( Event(self.events[-1],p1,p2,self.positiveY))
        pass

    def addToSkyline(self,v):
        # index of a numpy array in a list
        indx = next((i for i, val in enumerate(self.skyline) if np.all(np.isclose(val,v))), -1)
        if indx >= 0:
            del self.skyline[indx:] # remove back-spike loop
        self.skyline.append(v)
       
    def process(self):
        self.events.sort()
        pendingEvents= PriorityQueue()
        activeEvent = None

        for event in self.events:
            if not activeEvent:          # start of scan
                activeEvent = event
                self.addToSkyline(event.v1)
                continue
            
            elif not event.relStart: # a new event starts here
                if event.isInFrontOf(activeEvent):
                    self.addToSkyline( np.array((event.sX,activeEvent.yAt(event.sX))) )
                    self.addToSkyline(event.v1)
                    pendingEvents.push(activeEvent.eY,activeEvent)
                    activeEvent = event
                else:
                    pendingEvents.push(event.eY,event)

            else:
                if activeEvent is event.relStart:
                    self.addToSkyline(activeEvent.v2)
                    if not pendingEvents.empty():
                        activeEvent = pendingEvents.pop()
                        self.addToSkyline( np.array((event.eX,activeEvent.yAt(event.eX))) )
                    else:
                        activeEvent = None
                else:
                    pendingEvents.remove(event.relStart)

        # remove tail at end
        while abs(self.skyline[-1][0]-self.skyline[-2][0]) < 10.e-3:
            del self.skyline[-1]

        return self.skyline


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
        self.priorityQueue = PriorityQueue()

    def do(self, manager):
        self.createAllWaySectionNetwork()
        self.createPolylineVerts(manager)
        self.computePolylineVisibility()

    def cleanup(self):
        self.kdTree = None
        self.bldgVerts = None
        self.vertIndexToPolyIndex.clear()

    def createPolylineVerts(self,manager):
        # create polyline set with first vertex as key
        self.polylines = [PolyLine(polyline) for polyline in manager.polylines]

        # the total number of vertices
        totalNumVerts = sum(polyline.numVerts for polyline in self.polylines )

         # create mapping between the index of the vertex and index of the polyline in <polylines>
        self.vertIndexToPolyIndex.extend(
            polyIndex for polyIndex, polyline in enumerate(self.polylines) for _ in range(polyline.numVerts)
        )

        # allocate the memory for an empty numpy array
        self.polyVerts = np.zeros((totalNumVerts, 2))

        # fill vertices in <polyVerts>
        index = 0
        for polyline in self.polylines:
            # store the index of the first vertex of <polyline> in <self.polyVerts>
            polyline.auxIndx = index
            for vert in polyline.verts:
                self.polyVerts[index] = vert
                index += 1

        self.createKdTree()

    def computePolylineVisibility(self):
        from defs.way import allRoadwayCategoriesRank

        for crx in self.sectionNetwork.iterAllIntersectionNodes():
            plt.plot(crx[0],crx[1],'ko',zorder=300)

        for waySeg in self.sectionNetwork.iterAllSegments():
            finalPosSkyline = []
            finalNegSkyline = []
            if allRoadwayCategoriesRank[waySeg.category] > allRoadwayCategoriesRank['secondary']:
                continue
            if waySeg.ID != 5833:
                continue

            v1,v2 = waySeg.s, waySeg.t
            if True:
            # for v1,v2 in zip(waySeg.path[:-1],waySeg.path[1:]):
                v1, v2 = np.array(v1), np.array(v2) # <sectionNetwork> delivers mathutils.Vector's
                way = (v2-v1)
                wayLength = np.linalg.norm(way)
                wayCenter, wayUnitVector = (v1 + v2)/2., way/wayLength
                searchWidth = wayLength/2.

                # Get all polygon vertices for the polygons whose vertices
                # are returned by the query to the KD tree
                queryPolyIndices = set(
                    self.makeKdQuery(wayCenter, sqrt(searchWidth*searchWidth + self.searchHeight_2))
                )

                if queryPolyIndices:
                    queryPolyVertIndices = tuple(
                        vertIndex for polyIndex in queryPolyIndices for vertIndex in \
                            range(self.polylines[polyIndex].auxIndx, \
                            self.polylines[polyIndex].auxIndx + self.polylines[polyIndex].numVerts)
                    )

                    # slice <self.polyVerts> using <queryPolyIndices>
                    queryPolyVerts = self.polyVerts[queryPolyVertIndices,:]

                    # transform <queryPolygVerts> to the system of reference of <way> by ..
                    # ... computing the rotation matrix
                    matrix = np.array(( (wayUnitVector[0], -wayUnitVector[1]), (wayUnitVector[1], wayUnitVector[0]) ))
                    # ... shifting the origin to the segment center
                    queryPolyVerts -= wayCenter
                    # ... rotating shifted <queryPolygVerts>, output back to <queryPolygVerts>
                    np.matmul(queryPolyVerts, matrix, out=queryPolyVerts)

                    # # plot polylines
                    # firstVertIndex = 0
                    # for polyIndex in queryPolyIndices:
                    #     polyline = self.polylines[polyIndex]
                    #     for edge, edgeVert1, edgeVert2 in polyline.edgeInfo(queryPolyVerts, firstVertIndex):
                    #         v1, v2 = edgeVert1, edgeVert2
                    #         plt.plot([v1[0],v2[0]],[v1[1],v2[1]],'m:',linewidth=2)
                    #         plt.plot(v1[0],v1[1],'k.')
                    #     firstVertIndex += polyline.numVerts

                    posSkyliner = Skyline(True,searchWidth,self.searchHeight)
                    negSkyliner = Skyline(False,searchWidth,self.searchHeight)

                    firstVertIndex = 0
                    for polyIndex in queryPolyIndices:
                        polyline = self.polylines[polyIndex]
                        for edgeVert1, edgeVert2 in polyline.edgeInfo(queryPolyVerts, firstVertIndex):
                            if edgeVert1[1] > 0. or edgeVert2[1] > 0.:
                                posSkyliner.addSegment(edgeVert1, edgeVert2)
                            if edgeVert1[1] < 0. or edgeVert2[1] < 0.:
                                negSkyliner.addSegment(edgeVert1, edgeVert2)

                        firstVertIndex += polyline.numVerts

                    posSkyline = posSkyliner.process()
                    negSkyline = negSkyliner.process()

                    # back transform skylines
                    revMatrix = np.array(( (wayUnitVector[0], wayUnitVector[1]), (-wayUnitVector[1], wayUnitVector[0]) ))
                    for v in posSkyline:
                        np.matmul(v, revMatrix, out=v)
                        v += wayCenter 
                    for v in negSkyline:
                        np.matmul(v, revMatrix, out=v)
                        v += wayCenter

                    finalPosSkyline.extend(posSkyline) 
                    finalNegSkyline.extend(negSkyline) 

                    w = np.array(( (-searchWidth,0.),(searchWidth,0.) ))
                    np.matmul(w, revMatrix, out=w)
                    w += wayCenter 
                    w1, w2 = w[0],w[1]
                    plt.plot([w1[0],w2[0]],[w1[1],w2[1]],'b',linewidth = 3,zorder=300)

                polygon = finalPosSkyline[::-1] + finalNegSkyline
                plt.fill(
                    tuple(v[0] for v in polygon),
                    tuple(v[1] for v in polygon),
                    '#ff0000',
                    alpha=0.3,
                    zorder=300
                )

                for v1,v2 in zip(polygon[:-1],polygon[1:]):
                    plt.plot([v1[0],v2[0]],[v1[1],v2[1]],'r',linewidth=2,zorder=300)
                v1, v2 = polygon[-1], polygon[0]
                plt.plot([v1[0],v2[0]],[v1[1],v2[1]],'r',linewidth=2,zorder=300)


                # plt.gca().axis('equal')
                # plt.show()

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
            for segment in way.segments:
                netSeg = NetSegment(segment.v1,segment.v2,way.category,segment.length)
                self.allNetwork.addSegment(netSeg)

        # create way-section network
        self.sectionNetwork = createSectionNetwork(self.allNetwork)

    # def plotPolys(self):
    #     for v in self.polyVerts:
    #         plt.plot(v[0],v[1],'k.')
    # def plotEnd(self):
    #     plt.gca().axis('equal')
    #     plt.show()