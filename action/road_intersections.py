from collections import defaultdict, deque
from mathutils import Vector
import matplotlib.pyplot as plt

from lib.SweepIntersectorLib.SweepIntersector import SweepIntersector
from lib.CompGeom.algorithms import SCClipper, polygonAdjacencyRook
from lib.CompGeom.boundingRectangle import boundingRectangle

from way.way_network import WayNetwork, NetSection
from way.way_algorithms import createSectionNetwork
from way.way_intersections import Intersection
from way.way_section import WaySection
from way.way_graph_cycle import GraphCycle

from defs.road_polygons import ExcludedWayTags

class RoadIntersections:

    def __init__(self):
        self.networkGraph = None
        self.sectionNetwork = None
        self.waySections = dict()
        self.intersections = dict()
        self.transitionPolys = dict()
        self.intersectionPolys = dict()
        self.wayPolys = dict()

    def do(self, manager):
        evalType = 'overView'
        self.useFillet = True

        if evalType == 'overView':
            self.findSelfIntersections(manager)
            self.createWaySectionNetwork()
            self.createSections()
            self.createIntersectionPolys()
            self.createWayPolys()
            self.plotPolygonAreas()

        if evalType == 'cycleDetails':
            self.findSelfIntersections(manager)
            self.createWaySectionNetwork()
            self.createSections()
            self.createIntersectionPolys()
            self.typeDetection()

        elif evalType == 'clusterWays':
            self.findSelfIntersections(manager)
            self.createWaySectionNetwork()
            self.createSections()
            self.clusterElongatedWays()

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
        for way in wayManager.getAllIntersectionWays():
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
        for net_section in self.sectionNetwork.iterAllForwardSegments():
            if net_section.category != 'scene_border':
                section = WaySection(net_section)
                self.waySections[net_section.sectionId] = section
                # line = section.polyline
                # plotPolyLine(line,section.leftWidth,'#1f85ba')
        test=1

    def createIntersectionPolys(self):       
        for nodeNr,node in enumerate(self.sectionNetwork):
            # plt.text(node[0],node[1],str(isectNr),fontsize=12,zorder=900)
            intersection = Intersection(node, self.sectionNetwork, self.waySections)
            self.intersections[node] = intersection

        # Transitions have to be processed first, because way widths may be altered.
        for node,intersection in self.intersections.items():
            if intersection.order == 2:
                polygon = intersection.findTransitionPoly()
                if polygon:
                    self.transitionPolys[node] = polygon

        for node,intersection in self.intersections.items():
            if intersection.order > 2:
                debug = False
                if debug:
                    plt.close()
                if self.useFillet:
                    polygon = intersection.intersectionPoly()
                else:
                    polygon = intersection.intersectionPoly_noFillet()
                if debug:
                    plotPolygon(polygon,False,'k','r',1.,True,0.4,100)
                    plotEnd()
                if polygon:
                    self.intersectionPolys[node] = polygon

    def createWayPolys(self):
        for sectionNr,section in enumerate(self.waySections.values()):
            # center = sum(section.polyline.verts,Vector((0,0)))/len(section.polyline.verts)
            # plt.text(center[0],center[1],str(nr),zorder=120)
            waySlice = None
            wayLength = 0.
            if section.trimT > section.trimS:
                waySlice = section.polyline.trimmed(section.trimS,section.trimT)
                if section.turnParams:
                    waySegPoly = waySlice.parabolicBuffer(section.leftWidth, section.rightWidth,section.turnParams[0],section.turnParams[1])
                else:
                    waySegPoly = waySlice.buffer(section.leftWidth, section.rightWidth)
                self.wayPolys[sectionNr] = waySegPoly
            else:
                section.isValid = False

    def plotPolygonAreas(self):
        for node,polygon in self.transitionPolys.items():
            plotPolygon(polygon,False,'m','m',1.,True,0.4,100)
        for node,polygon in self.intersectionPolys.items():
            plotPolygon(polygon,False,'r','r',1.,True,0.4,200)
        for nr,polygon in self.wayPolys.items():
            plotPolygon(polygon,False,'b','b',1.,True,0.4,200)

    def clusterElongatedWays(self):
        hiThresh = 0.8
        loThresh = 0.7
        widthThresh = 20.
        cycles = []
        elongations = []
        widths = []

        cycleSegs, _, _ = self.sectionNetwork.getCycles()
        polyNr = 0
        for sectList in cycleSegs:
            graphCycle = GraphCycle(sectList)
            length, width, _, dx, dy = boundingRectangle(graphCycle.cyclePoly.coords)
            elongation = 1. - width/length

            # Split polygons with high elongation but large width into tiles
            if elongation > hiThresh and width > widthThresh:
                for slice in GraphCycle.sliceLargeCycles(graphCycle.cyclePoly):
                    length,width,_,dx,dy = boundingRectangle(slice.coords)
                    elongation = 1. - width/length
                    if elongation > hiThresh:
                        cycles.append(slice)
                        elongations.append(elongation)
                        widths.append(width)
            else:
                cycles.append(graphCycle.cyclePoly)
                elongations.append(elongation)
                widths.append(width)

        # hysteresis thresholding
        adj = polygonAdjacencyRook(cycles)
        queue = deque( i for i in range(len(elongations)) if elongations[i] > hiThresh )
        selected = set()
        while queue:
            indx = queue.popleft()
            if indx not in selected and elongations[indx] > loThresh:
                selected.add(indx)
                for i in adj[indx]:
                    if i != indx:
                        queue.append(i)

        for i in range(len(cycles)):
            cycVerts = [(v.x,v.y) for v in cycles[i].coords]
            if elongations[i] > hiThresh and widths[i] < widthThresh:
                plotPolygon(cycVerts,False,'k','g',0.5,True,0.4,999)
            elif i in selected and widths[i] < widthThresh:
                plotPolygon(cycVerts,False,'c','c',0.5,True,0.4,999)


    def typeDetection(self):
        # local functions --------------------------------------------------
        from math import atan2, pi
        def argmax(x):
            return max(range(len(x)), key=lambda i: x[i])
        def argmin(x):
            return min(range(len(x)), key=lambda i: x[i])
        def cyclicOrientationDiff(a,b):
            d = (pi - abs(pi - (2.*a-2.*b)%(2.*pi))) / 2.
            return d
        def divideAtEndpoints(poly,uV):
            # poly is list of mathutils Vectors
            # uV is unit Vector along elongation

            # project all polygon vertices onto elongation vector
            proj = [v.dot(uV) for v in poly]

            # find indices of endpoints
            mi, ma = argmin(proj), argmax(proj)
            mi, ma = (mi,ma) if mi<ma else (ma,mi)

            # show elongation direction
            # v1 = poly[mi]
            # v2 = poly[mi] + (poly[ma]-poly[mi]).length * uV
            # plt.plot( (v1[0], v2[0]), (v1[1], v2[1]),'r:')

            # create lines between these endpoints
            line1 = [poly[i] for i in range(mi,ma+1)]
            line2 = [poly[i] for i in range(mi,ma-len(poly)-1,-1)]

            return line1, line2
        def removePerpendEdges(line,dx,dy):
            # Removes perpendicular edges at the beginning and end of the line.
            # dx, dy: components of unit vector along elongation
            angleLim = 35./180.*pi  # 35°
            elongAngle = atan2(dy,dx)
            # remove at the beginning and at the end
            for s in [slice(None,None,1), slice(None,None,-1)]:
                while len(line)>2:
                    edge = line[s][1]-line[s][0]
                    angle = atan2(edge[1],edge[0])
                    if cyclicOrientationDiff(angle,elongAngle) < angleLim:
                        break
                    line = line[s][1:][s]
            return line

        def iniPlot():
            plt.close()
            plt.figure(figsize=(10,8))
        def plotCycle(sLists,nr):
            plt.subplot(2,2,1)
            for section in sLists[nr]:
                plotFull(section)
            plt.gca().axis('equal')
        def plotWaysAndIntersections(sLists,nr):
            plt.subplot(2,2,2)
            for section in sLists[nr]:
                if section.category == 'scene_border':
                    plotLine(section.path,'r:',2,100)
                else:
                    order = sum(1 for v in self.sectionNetwork.iterOutSegments(section.s))
                    if order > 1:
                        if section.s in self.intersectionPolys:
                            poly = self.intersectionPolys[section.s]
                            plotPolygon(poly,False,'k','r',1.,True,0.4,100)
                    section = self.waySections[section.sectionId]
                    if section.trimT > section.trimS:
                        waySlice = section.polyline.trimmed(section.trimS,section.trimT)
                        if section.turnParams:
                            waySegPoly = waySlice.parabolicBuffer(section.leftWidth, section.rightWidth,section.turnParams[0],section.turnParams[1])
                        else:
                            waySegPoly = waySlice.buffer(section.leftWidth, section.rightWidth)
                        plotPolygon(waySegPoly,False,'b','#aaaaff',1.,True,0.4,100)
            plt.gca().axis('equal')

        def plotRemainingPolygon(sLists,nr):
            from lib.pygeos.geom import GeometryFactory
            from lib.pygeos.shared import TopologyException, PrecisionModel
            from lib.pygeos.shared import CAP_STYLE,JOIN_STYLE
            # precModel = PrecisionModel(scale=10)
            # geosF = GeometryFactory(precisionModel=precModel)
            geosF = GeometryFactory()

            plt.subplot(2,2,4)

            polyList = []
            for section in sLists[nr]:
                if section.category == 'scene_border':
                    verts = section.path
                    way_coords = [ geosF.createCoordinate(v) for v in verts ]
                    way_line = geosF.createLineString(way_coords)
                    poly = way_line.buffer(0.001,cap_style=CAP_STYLE.round,join_style=JOIN_STYLE.round )
                    polyList.append(poly)
                else:
                    order = sum(1 for v in self.sectionNetwork.iterOutSegments(section.s))
                    if order > 1:
                        if section.s in self.intersectionPolys:
                            polyV = self.intersectionPolys[section.s]
                            poly = geosF.createPolygon(geosF.createLinearRing([geosF.createCoordinate(v) for v in polyV] + [geosF.createCoordinate(polyV[0])]))
                            polyList.append(poly)
                    waySection = self.waySections[section.sectionId]
                    totalT = waySection.trimS + waySection.trimT
                    if 0. <= totalT < len(waySection.polyline)-1:
                        waySlice = waySection.polyline.slice(waySection.trimS,waySection.trimT)
                        way_coords = [ geosF.createCoordinate(v) for v in waySlice.verts ]
                        way_line = geosF.createLineString(way_coords)
                        width = (waySection.leftWidth + waySection.rightWidth)/2.
                        poly = way_line.buffer(width,cap_style=CAP_STYLE.round,join_style=JOIN_STYLE.round )

                        # if waySection.turnParams:
                        #     polyV = waySlice.parabolicBuffer(waySection.leftWidth, waySection.rightWidth,waySection.turnParams[0],waySection.turnParams[1])
                        # else:
                        #     polyV = waySlice.buffer(waySection.leftWidth, waySection.rightWidth)
                        # poly = geosF.createPolygon(geosF.createLinearRing([geosF.createCoordinate(v) for v in polyV] + [geosF.createCoordinate(polyV[0])]))

                        polyList.append(poly)

                for v1,v2 in zip(section.path[:-1],section.path[1:]):
                    plt.plot( (v1[0], v2[0]), (v1[1], v2[1]), 'k:', zorder=50 )
                    # plt.plot(v1[0],v1[1],'k.',zorder=55)
                    # plt.plot(v2[0],v2[1],'k.',zorder=55)

            multiPoly = geosF.createMultiPolygon(polyList)
            try:
                merged = multiPoly.union()
            except (TopologyException) as e:
                pass
            else:
                if merged.geom_type == 'Polygon':
                    for ring in merged.interiors:
                        p = [(c.x,c.y) for c in ring.coords]
                        plotPoly(p,False,'r',1,200)
                else:
                    for geom in merged.geoms:
                        for ring in geom.interiors:
                            p = [(c.x,c.y) for c in ring.coords]
                            plotPoly(p,False,'r',1,200)
                
                plt.gca().axis('equal')

        def plotLine(line,color='k',width = 1, order=10):
            x = [n[0] for n in line]
            y = [n[1] for n in line]
            plt.plot(x,y,color,linewidth=width,zorder=order)
        def plotFull(section):
            from mpl.renderer.road_polygons import RoadPolygonsRenderer
            for v1,v2 in zip(section.path[:-1],section.path[1:]):
                plt.plot( (v1[0], v2[0]), (v1[1], v2[1]), **RoadPolygonsRenderer.styles[section.category], zorder=50 )
                # plt.plot(v1[0],v1[1],'k.',zorder=55)
                # plt.plot(v2[0],v2[1],'k.',zorder=55)
            plt.plot(section.path[0][0],section.path[0][1],'bo',markersize=6,zorder=56)
            plt.plot(section.path[-1][0],section.path[-1][1],'bo',markersize=6,zorder=56)

        # end of local functions -------------------------------------------

        cycleSegs, _, _ = self.sectionNetwork.getCycles()

        cycles = []
        sLists = []
        for sectList in cycleSegs:
            cycles.append(GraphCycle(sectList).cyclePoly)
            sLists.append(sectList)
        for nr,cycle in enumerate(cycles):
            length, width, _, dx, dy = boundingRectangle(cycle.coords)
            elongation = 1. - width/length
            # if nr != 54:
            #     continue

            if elongation > 0.8:
                iniPlot()
                plotCycle(sLists,nr)
                plt.title('Map')
                plotWaysAndIntersections(sLists,nr)
                plt.title('Ways & Intersection Areas')

                poly = [Vector((v.x,v.y)) for v in cycle.coords[:-1]]
                line1, line2 = divideAtEndpoints(poly,Vector((dx,dy)))
                line1 = removePerpendEdges(line1,dx,dy)
                line2 = removePerpendEdges(line2,dx,dy)

                plt.subplot(2,2,3)
                plotLine(line1,'r',2)
                plotLine(line2,'b',2)
                plotPolygon(poly,False,'k:','k',1.,False,0.2,55)
                plt.title('Long Line Pair')

                plt.gca().axis('equal')

                # plotRemainingPolygon(sLists,nr)
                # plt.title('Inner Area')

                plt.subplot(2,2,4)
                plotNetwork(self.sectionNetwork)
                p = [(c.x,c.y) for c in cycle.coords]
                plotPoly(p,False,'r',4,200)
                plt.title('Position')
                plt.gca().axis('equal')

                plt.suptitle('Cycle '+str(nr)+'     (close figure to see next cycle)')
                plt.tight_layout()
                plt.show()




from itertools import *
def _iterEdges(poly):
        p1, p2= tee(poly)
        p2 = islice(cycle(p2), 1, None)
        return zip(p1,p2)


def plotPolygon(poly,vertsOrder,lineColor='k',fillColor='k',width=1.,fill=True,alpha = 0.2,order=100):
    x = [n[0] for n in poly] + [poly[0][0]]
    y = [n[1] for n in poly] + [poly[0][1]]
    if fill:
        plt.fill(x[:-1],y[:-1],color=fillColor,alpha=alpha,zorder = order)
    plt.plot(x,y,lineColor,linewidth=width,zorder=order)
    if vertsOrder:
        for i,(xx,yy) in enumerate(zip(x[:-1],y[:-1])):
            plt.text(xx,yy,str(i))

def plotNetwork(network):
    from mpl.renderer.road_polygons import RoadPolygonsRenderer
    for count,seg in enumerate(network.iterAllSegments()):
        # plt.plot(seg.s[0],seg.s[1],'k.')
        # plt.plot(seg.t[0],seg.t[1],'k.')
        color = 'r' if seg.category=='scene_border' else 'y'

        for v1,v2 in zip(seg.path[:-1],seg.path[1:]):
            plt.plot( (v1[0], v2[0]), (v1[1], v2[1]), **RoadPolygonsRenderer.styles[seg.category], zorder=50 )
        # plotWaySeg(seg,color,1)

def plotWaySeg(wayseg,color='k',width=1.,order=100):
    for v1,v2 in zip(wayseg.path[:-1],wayseg.path[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
        # plt.plot(v1[0],v1[1],'k.')
        # plt.plot(v2[0],v2[1],'k.')
        x = (v1[0]+v2[0])/2
        y = (v1[1]+v2[1])/2
        # plt.text(x,y,str(wayseg.ID),fontsize=12)

def plotGeosPoly(geosPoly,vertsOrder,color='k',width=1.,order=100):
    poly = [(c.x,c.y) for c in geosPoly.coords[:-1]]
    plotPoly(poly,vertsOrder,color,width,order)

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

def plotGeosWithHoles(geosPoly,vertsOrder,color='k',width=1.,order=100):
    poly = [(c.x,c.y) for c in geosPoly.exterior.coords[:-1]]
    plotPoly(poly,vertsOrder,color,width,order)
    for ring in geosPoly.interiors:
        p = [(c.x,c.y) for c in ring.coords[:-1]]
        plotPoly(p,vertsOrder,'r',width,order)

def plotGeneralPoly(geosPoly,vertsOrder,color='k',width=1.,order=100):
    if geosPoly.geom_type == 'Polygon':
         plotGeosWithHoles(geosPoly,vertsOrder,color,width,order)
    else:
        for geom in geosPoly.geoms:
            plotGeosWithHoles(geom,vertsOrder,color,width,order)


def plotPolygon(poly,vertsOrder,lineColor='k',fillColor='k',width=1.,fill=True,alpha = 0.2,order=100):
    x = [n[0] for n in poly] + [poly[0][0]]
    y = [n[1] for n in poly] + [poly[0][1]]
    if fill:
        plt.fill(x[:-1],y[:-1],color=fillColor,alpha=alpha,zorder = order)
    plt.plot(x,y,lineColor,linewidth=width,zorder=order)
    if vertsOrder:
        for i,(xx,yy) in enumerate(zip(x[:-1],y[:-1])):
            plt.text(xx,yy,str(i),fontsize=12)

def plotEnd():
    plt.gca().axis('equal')
    plt.show()