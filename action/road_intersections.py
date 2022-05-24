from collections import defaultdict
from mathutils import Vector
import matplotlib.pyplot as plt

from lib.SweepIntersectorLib.SweepIntersector import SweepIntersector
from lib.CompGeom.algorithms import SCClipper

from way.way_network import WayNetwork, NetSection
from way.way_algorithms import createSectionNetwork
from way.way_intersections import Intersection
from way.way_section import WaySection

from defs.road_polygons import ExcludedWayTags

# from plotUtilities import *

class RoadIntersections:

    def __init__(self):
        self.networkGraph = None
        self.sectionNetwork = None
        self.waySections = dict()
        self.intersections = []

    def do(self, manager):
        self.findSelfIntersections(manager)
        self.createWaySectionNetwork()
        self.createSections()
        self.createIntersections()
        self.plotSections()
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
                section = WaySection(net_section)
                self.waySections[net_section.sectionId] = section
                # line = section.polyline
                # plotPolyLine(line,section.leftWidth,'#1f85ba')
        test=1

    def createIntersections(self):
        for nodeNr,node in enumerate(self.sectionNetwork):
            # print(nodeNr)
            # plt.text(node[0],node[1],str(nodeNr),fontsize=12,zorder=900)
            # plt.close()
            # if nodeNr != 213:
            #     test=1
            #     continue
            intersection = Intersection(node, self.sectionNetwork, self.waySections)

            if intersection.order > 2:
                polygon = intersection.findIntersectionPoly() 
                plotPolygon(polygon,False,'k','r',1.,True,0.2,100)

                # plt.title(str(nodeNr))
                # plotEnd()

    def plotSections(self):
        for nr,section in enumerate(self.waySections.values()):
            totalT = section.trimS +section.trimT
            if 0. <= totalT < len(section.polyline)-1:
                # center = sum(section.polyline.verts,Vector((0,0)))/len(section.polyline.verts)
                # plt.text(center[0],center[1],str(section.id),zorder=120)
                # print(nr)
                waySlice = section.polyline.slice(section.trimS,section.trimT)
                waySegPoly = waySlice.buffer(section.leftWidth, section.rightWidth)
                plotPolygon(waySegPoly,False,'b','#aaaaff',1.,True,0.2,100)
            # else:
            #     center = sum(section.polyline.verts,Vector((0,0)))/len(section.polyline.verts)
            #     trimText = ' %4.2f %4.2f'%(section.trimS,section.trimT)
                # plt.text(center[0],center[1],str(section.id)+trimText,zorder=120)
                # print(nr)
                # waySegPoly = section.polyline.buffer(section.leftWidth, section.rightWidth)
                # plotPolygon(waySegPoly,False,'g','#00ff00',1.,True,0.8,100)

        test=1

def plotPolygon(poly,vertsOrder,lineColor='k',fillColor='k',width=1.,fill=True,alpha = 0.2,order=100):
    x = [n[0] for n in poly] + [poly[0][0]]
    y = [n[1] for n in poly] + [poly[0][1]]
    if fill:
        plt.fill(x[:-1],y[:-1],color=fillColor,alpha=alpha,zorder = order)
    plt.plot(x,y,lineColor,linewidth=width,zorder=order)
    if vertsOrder:
        for i,(xx,yy) in enumerate(zip(x[:-1],y[:-1])):
            plt.text(xx,yy,str(i))
