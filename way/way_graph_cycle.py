from collections import defaultdict
from mathutils import Vector
from math import inf

from lib.pygeos.geom import GeometryFactory
from lib.pygeos.shared import TopologyException

class GraphCycle():
    ID = 0
    def __init__(self,segList):
        self.id = GraphCycle.ID
        GraphCycle.ID += 1
        self.subPolys = []
        self.triangles = []
 
        boundary, holes, spurs = GraphCycle.cleanCycle(segList)
        geosF = GeometryFactory()
        coords = [ geosF.createCoordinate(v) for v in boundary ] 
        self.cyclePoly = geosF.createPolygon(geosF.createLinearRing(coords))
        if holes:
            for hole in holes:
                coords = [ geosF.createCoordinate(v) for v in hole ] 
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

        # follow boundary, begin at <boundStart>
        boundary.append(boundStart)
        nextNode = boundStart
        while True:
            thisNode = segments[nextNode][0]
            boundary.append(thisNode)
            del segments[nextNode]
            if thisNode == boundStart:
                break
            nextNode = thisNode

        # collect remaining holes, if any
        while segments:
            firstNode = next(iter(segments))
            hole = [firstNode]
            nextNode = firstNode
            while True:
                thisNode = segments[nextNode][0]
                hole.append(thisNode)
                del segments[nextNode]
                if thisNode == firstNode:
                    break
                nextNode = thisNode
            holes.append(hole)

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

