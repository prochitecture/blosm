from collections import defaultdict
from mathutils import Vector
import matplotlib.pyplot as plt
from  action.plotUtilities import *
from itertools import combinations, chain, product

from lib.SweepIntersectorLib.SweepIntersector import SweepIntersector
from lib.CompGeom.algorithms import SCClipper

from way.way_network import WayNetwork, NetSection
from way.way_algorithms import createSectionNetwork
from way.way_intersections import Section, Intersection

from defs.road_polygons import ExcludedWayTags

class RoadIntersections:

    def __init__(self):
        self.networkGraph = None
        self.sectionNetwork = None
        self.sections = dict()
        self.intersections = []

    def do(self, manager):
        self.findSelfIntersections(manager)
        self.createWaySectionNetwork()
        self.createSections()
        self.createIntersections()
        pass

    def findSelfIntersections(self, manager):
        wayManager = self.app.managersById["ways"]

        # some way tags to exclude, used also in createWaySectionNetwork(),
        # should be moved to defs.
        uniqueSegments = defaultdict(set)
        for way in wayManager.getAllWays():            
            if [tag for tag in ExcludedWayTags if tag in way.element.tags]:
                continue
            for segment in way.segments:
                v1, v2 = (segment.v1[0],segment.v1[1]),  (segment.v2[0],segment.v2[1])
                if v1 not in uniqueSegments.get(v2,[]):
                    uniqueSegments[v1].add(v2)
        cleanedSegs = [(v1,v2) for v1 in uniqueSegments for v2 in uniqueSegments[v1]]

        intersector = SweepIntersector()
        self.intersectingSegments = intersector.findIntersections(cleanedSegs)

    # Creates the network graph <self.sectionNetwork> for way-sctions (ways between crossings)
    def createWaySectionNetwork(self):
        # get border polygon (ounter-clockwise) of scene frame
        minX, minY = self.app.projection.fromGeographic(self.app.minLat, self.app.minLon)
        maxX, maxY = self.app.projection.fromGeographic(self.app.maxLat, self.app.maxLon)

        # prepare clipper for this frame
        clipper = SCClipper(minX,maxX,minY,maxY)

        wayManager = self.app.managersById["ways"]
        # Not really used. This is a relict from way_clustering.py
        wayManager.junctions = (
            [],#mainJunctions,
            []#smallJunctions
        )

         # create full way network
        wayManager.networkGraph = self.networkGraph = WayNetwork()

        # some way tags to exclude, used also in findSelfIntersections(),
        # should be moved to defs.
        for way in wayManager.getAllWays():
            # Exclude ways with unwanted tags
            if [tag for tag in ExcludedWayTags if tag in way.element.tags]:
                continue

            for waySegment in way.segments:
                # Check for segments splitted by self-intersections
                segments = []
                newSegments = self.intersectingSegments.get( (tuple(waySegment.v1),tuple(waySegment.v2)), None)
                if newSegments:
                    for v1,v2 in zip(newSegments[:-1],newSegments[1:]):
                        segments.append((v1,v2))
                else:
                    segments.append((waySegment.v1,waySegment.v2))

                for segment in segments:
                    way.element.tags
                    v1, v2 = Vector(segment[0]),Vector(segment[1])
                    accepted, v1, v2 = clipper.clip(v1,v2)
                    if accepted:
                        netSeg = NetSection(v1,v2,way.category,way.element.tags,(v2-v1).length)
                        wayManager.networkGraph.addSegment(netSeg,False)

        borderPolygon = clipper.getPolygon()
        for v1,v2 in zip(borderPolygon[:-1],borderPolygon[1:]):
            netSeg = NetSection(v1,v2,'scene_border',None, (v2-v1).length) 
            wayManager.networkGraph.addSegment(netSeg)

        # create way-section network
        wayManager.sectionNetwork = self.sectionNetwork = createSectionNetwork(wayManager.networkGraph)

    def createSections(self):
        for net_section in self.sectionNetwork.iterAllSegments():
            if net_section.category != 'scene_border':
                section = Section(net_section)
                self.sections[net_section.sectionId] = section
                line = section.polyline
                plotPolyLine(line,section.halfWidth,'#1f85ba')
        test=1

    def createIntersections(self):
        processed_nodes = []
        for nodeNr,node in enumerate(self.sectionNetwork):
            # plt.text(node[0],node[1],str(nodeNr),fontsize=12,zorder=900)
            # plt.close()
            # if nodeNr == 227:
            #     test=1
            #     continue
            if node not in processed_nodes:
                # processed_nodes.append(node)
                intersection = Intersection(node,self.sections)
                nodes_to_merge = intersection.mergeNode(node,self.sections,self.sectionNetwork)
                # while nodes_to_merge:
                #     node = nodes_to_merge.pop()
                #     if node not in processed_nodes:
                #         processed_nodes.append(node)
                #         nodes_to_merge.extend(intersection.mergeNode(node,self.sections,self.sectionNetwork))

                if len(intersection.outPolyLines) > 2:
                    intersection.sortSections()
                    intersection.findIntersectionPoly()
                    # intersection.findCollisions()

                    # for line in intersection.outPolyLines:
                    #     plotPolyLine(line['line'],line['halfWidth'],'#0000ff')
                    #     plotPolyLine(line['line'],0.,'#0000ff')
                    #     c = intersection.position
                    #     plt.plot(c.x,c.y,'ko',markersize=2)
                    # for id in intersection.shortSectionIDs:
                    #     line = self.sections[id].polyline
                    #     plotPolyLine(line,self.sections[id].halfWidth,'#00ff00')

                    # plt.title(str(nodeNr))
                    # plotEnd()

        test=1

                # width = 4.#getRoadWidth(section.tags)
                # color = '#ff0000' if width == 4 * 2.6 else '#0000ff'
                # geosCoords = [geosF.createCoordinate(v) for v in section.path]
                # geosString = geosF.createLineString(geosCoords)
                # buffer = geosString.buffer(width/2.,resolution=3,cap_style=CAP_STYLE.flat)
                # plotGeosPolyFill(buffer,color,1,0.4,600)

        # for node in self.sectionNetwork:
        #     isect = Intersection(node)
        #     isect.addSections([section for section in self.sectionNetwork.iterOutSegments(node) if section.category != 'scene_border'])
        #     if len(isect.sections) < 3:
        #         continue
        #     # isect.boundaryCollisions()
        #     isect.unionIntersections()

        #     for sect in isect.sections:
        #         section = sect.originalSection
        #         width = sect.halfWidth
        #         geosCoords = [geosF.createCoordinate(v) for v in section.path]
        #         geosString = geosF.createLineString(geosCoords)
        #         buffer = geosString.buffer(width,resolution=3,cap_style=CAP_STYLE.flat)
        #         plotGeosPolyFill(buffer,'b',1,0.2,600)
        #         plotGeomPoly(buffer,False,'b')
        #         for v1,v2 in zip(section.path[:-1],section.path[1:]):
        #             plt.plot([v1[0],v2[0]],[v1[1],v2[1]],'k')

        #         for c in sect.collisions:
        #             plt.plot(c.x,c.y,'rx',markersize=8)

        #     plotEnd()

        #     # if len(isect.ways) < 2:
        #     #     continue

        #     # isect.ways = sorted(isect.ways, key = sort_key)
        #     # isect.centerLines = []
        #     # isect.leftLines = []
        #     # isect.rightLines = []
        #     # isect.widths = []

        #     # # find intersections
        #     # geosF = GeometryFactory()
        #     # plt.plot(node.x,node.y,'ro')
        #     # for way in isect.ways:
        #     #     width = getRoadWidth(way.tags)
        #     #     centerLine = geosF.createLineString([geosF.createCoordinate(v) for v in way.path])
        #     #     leftLine = centerLine.parallel_offset(width/2,3,JOIN_STYLE.round,5)
        #     #     rightLine = centerLine.parallel_offset(-width/2,3,JOIN_STYLE.round,5)
        #     #     isect.centerLines.append(centerLine)
        #     #     isect.leftLines.append(leftLine)
        #     #     isect.rightLines.append(rightLine)
        #     #     isect.widths.append(width)
        #     #     plotGeosString(centerLine,False,'k')
        #     #     plotGeosString(leftLine,False,'b')
        #     #     plotGeosString(rightLine,False,'r')
        #     # # plotEnd()

        #     # for leftLine,rightLine in zip(isect.leftLines,isect.rightLines[1:]+[isect.rightLines[0]]):
        #     #     leftL = leftLine.coords
        #     #     rightL = rightLine.coords
        #     #     leftsegs = [((v1.x,v1.y),(v2.x,v2.y)) for v1,v2 in zip(leftL[:-1],leftL[1:])]
        #     #     rightsegs = [((v1.x,v1.y),(v2.x,v2.y)) for v1,v2 in zip(rightL[:-1],rightL[1:])]
        #     #     segs = leftsegs + rightsegs
        #     # #     # plotGeom(leftLine,False,'b')
        #     #     # plotGeom(rightLine,False,'r')
        #     #     # plotEnd()
        #     #     intersector = SweepIntersector()
        #     #     intersectingSegments = intersector.findIntersections(segs)
        #     #     # print(len(intersectingSegments))
        #     #     for seg,isects in intersectingSegments.items():
        #     #         for p in isects[1:-1]:
        #     #             plt.plot(p[0],p[1],'rx',markersize=8)
        #     #         break
        #     # plotEnd()
        #     # test=1


