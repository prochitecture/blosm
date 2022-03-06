from collections import defaultdict
from mathutils import Vector
from math import inf
import matplotlib.pyplot as plt


from lib.pygeos.geom import GeometryFactory
from lib.pygeos.shared import TopologyException

class GraphCycle():
    ID = 0
    def __init__(self,segList):
        self.id = GraphCycle.ID
        GraphCycle.ID += 1
        self.subPolys = []
        self.triangles = []

        try:
            boundary, holes, spurs = GraphCycle.cleanCycle(segList)
        except:
            nodes = [n for s in segList for n in s.path[:-1]]
            plotPoly(nodes,True,'b',3)
            print('Problem on cleaning graph-cycle')
            return

        geosF = GeometryFactory()
        coords = [ geosF.createCoordinate(v) for v in boundary+[boundary[0]] ] 
        self.cyclePoly = geosF.createPolygon(geosF.createLinearRing(coords))
        if holes:
            for hole in holes:
                coords = [ geosF.createCoordinate(v) for v in hole+[hole[0]] ] 
                holePoly = geosF.createPolygon(geosF.createLinearRing(coords))
                try:
                    self.cyclePoly = self.cyclePoly.difference(holePoly)
                except (TopologyException,ValueError) as e:
                    import traceback
                    traceback.print_exception(type(e), e, e.__traceback__)
        self.spurs = spurs

    @staticmethod
    def cleanCycle(segList):
        boundary = []   # boundary of the cycle
        holes = []      # holes in the cycle
        spurs = []      # removed spur segments

        # remove coincident segments (spurs)
        nodes = [n for s in segList for n in s.path[:-1]]
        boundStart = Vector((inf,0.))     # start point for boundary
        segments = defaultdict(list)
        for source, target in zip(nodes,nodes[1:]+[nodes[0]]):
            if source in segments.get(target,[]):
                spurs.append((source,target))
                # edge and its twin detected -> coincident
                segments[target].remove(source)
                if not segments[target]:
                    del segments[target]
            else:
                segments[source].append(target)
                if source[0] < boundStart[0]:
                    boundStart = source

        sharedVertices = [ source for source,targets in segments.items() if len(targets)>1 ]

        # follow boundary, begin at <boundStart>
        boundary.append(boundStart)
        nextNode = boundStart
        while segments:
            if len(segments[nextNode]) > 1:
                # Graph-cycles are delivered in counter-clockwise order. To keep the outer
                # boundary in case of an intersection, we have to take the rightmost vertex.
                ccw = (segments[nextNode][0]-nextNode).cross(segments[nextNode][1]-nextNode)
                thisNode = segments[nextNode].pop(0) if ccw >= 0 else segments[nextNode].pop(1)
            else:
                thisNode = segments[nextNode].pop(0)
            if not segments[nextNode]: del segments[nextNode]
            if thisNode == boundStart:
                break
            boundary.append(thisNode)
            nextNode = thisNode

        # collect remaining holes, if any
        while segments:
            firstNode = next(iter(segments))
            hole = [firstNode]
            nextNode = firstNode
            while True:
                thisNode = segments[nextNode].pop(0)
                if not segments[nextNode]: del segments[nextNode]
                if thisNode == firstNode:
                    break
                hole.append(thisNode)
                nextNode = thisNode
            holes.append(hole)

        if sharedVertices:
            # Here we treat the case of shared vertices between the boundary and
            # the holes, which can't finally be triangulated by our algorithm.
            # such cases my look like this:
            #          -------------
            #         |             |
            #     ----O-------      |
            #     |   |       |     |
            #     |    -------      |
            #     |                 |
            #      -----------------
            # We correct the shared vertex of the hole by moving it 1mm along
            # the bisector into the hole.
            for shared in sharedVertices:
                for hole in holes:
                    if shared in hole:
                        # find the shared vertex and its neighbors in hole 
                        indx =  hole.index(shared)
                        p1 = hole[indx-1]
                        p2 = hole[(indx+1)%len(hole)]
                        # unity vectors to the neigbors
                        u1 = (p1-shared)/(p1-shared).length
                        u2 = (p2-shared)/(p2-shared).length
                        # move this point 1mm along bisector
                        hole[indx] = shared + (u1+u2)/(u1+u2).length * 0.001

        return boundary, holes, spurs

    @staticmethod
    def createHoleCycles(segList):
        geosF = GeometryFactory()
        holePolys = []      # holes in the cycle
        spurs = []      # removed spur segments

        nodes = [n for s in segList for n in s.path[:-1]]

        # remove coincident segments (spurs)
        segments = defaultdict(list)
        for source, target in zip(nodes,nodes[1:]+[nodes[0]]):
            if source in segments.get(target,[]):
                spurs.append((source,target))
                # edge and its twin detected -> coincident
                segments[target].remove(source)
                if not segments[target]:
                    del segments[target]
            else:
                segments[source].append(target)

        # collect closed contours, if any
        while segments:
            firstNode = next(iter(segments))
            contour = [firstNode]
            nextNode = firstNode
            while True:
                thisNode = segments[nextNode][0]
                contour.append(thisNode)
                del segments[nextNode]
                if thisNode == firstNode:
                    break
                nextNode = thisNode

            coords = [ geosF.createCoordinate(v) for v in contour ] 
            holePoly = geosF.createPolygon(geosF.createLinearRing(coords))
            holePolys.append(holePoly)

        return holePolys

def plotPoly(polygon,vertsOrder,color='k',width=1.,order=100):
    count = 0
    for v1,v2 in zip(polygon[:-1],polygon[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
        if vertsOrder:
            plt.text(v1[0],v1[1],str(count),fontsize=12)
        count += 1
        # plt.plot(v1[0],v1[1],'kx')
    # v1, v2 = polygon[-1], polygon[0]
    # plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
    if vertsOrder:
        plt.text(v1[0],v1[1],str(count),fontsize=12)

def plotSegments(adj):
    for v1,v2s in adj.items():
        for v2 in v2s:
            plt.plot([v1[0],v2[0]],[v1[1],v2[1]],'k')

def plotEnd():
    plt.gca().axis('equal')
    plt.show()