from collections import deque
from way.way_network import WayNetwork, NetSegment
from way.way_algorithms import createSectionNetwork, findWayJunctionsFor
from defs.way import allRoadwayCategories, mainRoads, smallRoads
from lib.CompGeom.algorithms import SCClipper
from mathutils import Vector

class WayClustering:
    
    def __init__(self):
        self.networkGraph = None
    
    def do(self, wayManager):
        # get border polygon (ounter-clockwise) of scene frame
        minX, minY = self.app.projection.fromGeographic(self.app.minLat, self.app.minLon)
        maxX, maxY = self.app.projection.fromGeographic(self.app.maxLat, self.app.maxLon)

        # prepare clipper for this frame
        clipper = SCClipper(minX,maxX,minY,maxY)
        # create full way network
        wayManager.networkGraph = WayNetwork()
        for way in wayManager.getAllWays():
            for segment in way.segments:
                v1, v2 = Vector(segment.v1),Vector(segment.v2)
                accepted, v1, v2 = clipper.clip(v1,v2)
                if accepted:
                    netSeg = NetSegment(v1,v2,way.category,(v2-v1).length)
                    wayManager.networkGraph.addSegment(netSeg)
        borderPolygon = clipper.getPolygon()
        for v1,v2 in zip(borderPolygon[:-1],borderPolygon[1:]):
            netSeg = NetSegment(v1,v2,'scene_border',(v2-v1).length)
            wayManager.networkGraph.addSegment(netSeg)
        # create way-section network
        graph = wayManager.waySectionGraph = createSectionNetwork(wayManager.networkGraph)

        # uncomment to display cycles
        # cycles = graph.iterCycles()

        # find way-junctions for principal roads
        allCrossings = graph.getCrossingsThatContain(allRoadwayCategories)
        mainJunctions = findWayJunctionsFor(graph, allCrossings, mainRoads, 20.)

        # expand them with near crossings of small roads
        for cluster in mainJunctions:
            for side_cluster in findWayJunctionsFor(graph, cluster, smallRoads, 15.):
                cluster |= side_cluster

        # remove these crossings from <allCrossings>
        remainingCrossings = list({crossing for crossing in allCrossings} -\
                        {crossing for cluster in mainJunctions for crossing in cluster })

        # find way-junctions for small roads in <remainingCrossings>
        smallJunctions = findWayJunctionsFor(graph, remainingCrossings, smallRoads, 15.)

        wayManager.junctions = (
            [],#mainJunctions,
            []#smallJunctions
        )
    
    def cleanup(self):
        pass

