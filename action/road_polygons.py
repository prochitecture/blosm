from mathutils import Vector
import numpy as np
from itertools import *
from collections import defaultdict

import matplotlib.pyplot as plt

from way.way_network import WayNetwork, NetSegment
from way.way_algorithms import createSectionNetwork

from lib.pygeos.geom import GeometryFactory
from lib.pygeos.shared import CAP_STYLE
from lib.pygeos.polygonDecomposition import polygonDecomposition
from lib.CompGeom.algorithms import circumCircle, SCClipper
from lib.CompGeom.splitPolygonHoles import splitPolygonHoles
from lib.CompGeom.poly_point_isect import isect_segments_include_segments

# cyclic iterate over polygons vertices
def _iterEdges(poly):
        p1, p2= tee(poly)
        p2 = islice(cycle(p2), 1, None)
        return zip(p1,p2)

# Gets all buildings in scene and combines those with shared edges to
# building blocks. Returns them as a list <geosPolyList> of PyGeos <Polynom>.
def mergeBuildingsToBlocks(buildings):
    segList = dict()
    for building in buildings:
        verts = [v for v in building.polygon.verts]
        for _, v1, v2 in building.polygon.edgeInfo(verts, 0, skipShared=True):
            segList[v1.freeze()] = v2

    geosF = GeometryFactory()
    geosPolyList = []
    while segList:
        firstVert = next(iter(segList))
        vertList = [Vector((firstVert[0],firstVert[1]))]
        nextVert = segList.pop(firstVert, None).freeze()
        while nextVert != firstVert and nextVert is not None:
            vertList.append(Vector((nextVert[0],nextVert[1])))
            nextVert = segList.pop(nextVert, None)
            if nextVert:
                nextVert.freeze()
        if len(vertList) > 2:
            geosCoords = [geosF.createCoordinate(v) for v in vertList+[vertList[0]]]
            geosRing = geosF.createLinearRing(geosCoords)
            geosPoly = geosF.createPolygon(geosRing)

            geosPolyList.append( geosPoly )
    return geosPolyList

class RoadPolygons:

    def __init__(self):
        self.sectionNetwork = None
        self.geosF = GeometryFactory()
        self.intersectingSegments = defaultdict(list)
        self.geosPolyList = None
        self.polyVerts = None
        self.vertIndexToPolyIndex = []
        self.kdTree = None
        self.cyclePolys = None

    def do(self, manager):
        self.findSelfIntersections(manager)
        self.createWaySectionNetwork()
        self.createPolylinePolygons(manager)
        self.createCyclePolygons()       
        self.createWayEnvironmentPolygons()

    def cleanup(self):
        self.kdTree = None
        self.bldgVerts = None
        self.vertIndexToPolyIndex.clear()

    def findSelfIntersections(self, manager):
        wayManager = self.app.managersById["ways"]

        # some way tags to exclude, used also in createWaySectionNetwork(),
        # should be moved to defs.
        excludedTags = ['steps']

        segList = []
        for way in wayManager.getAllWays():            
            if [tag for tag in excludedTags if tag in way.element.tags]:
                continue
            for segment in way.segments:
                v1, v2 = (segment.v1[0],segment.v1[1]),  (segment.v2[0],segment.v2[1])
                segList.append((v1,v2))

        # Unfortunatel, polylines can not be added to the intersection check,
        # as they may be collinear to ways, which is not supported by the
        # Bentley-Ottmann sweep-line algorithm.
        # for polyline in manager.polylines:
        #     for edge in polyline.edges:
        #         v1, v2 = (edge.v1[0],edge.v1[1]),  (edge.v2[0],edge.v2[1])
        #         segList.append((v1,v2))

        # Find self-intersections using Bentley-Ottmann sweep-line algorithm
        isects = isect_segments_include_segments(segList)

        # The algorithm delivers edge-ends, that are only nearby to the real edges,
        # and even sometimes reversed. nearbySeg() finds the true end-points.
        def nearbySeg(isectSeg,segLIst):
            s1,s2 = Vector(isectSeg[0]), Vector(isectSeg[1])
            for listSeg in segList:
                v1,v2 = Vector(listSeg[0]), Vector(listSeg[1])
                # nearby?
                if (v1-s1).length + (v2-s2).length < 0.001:
                    return (listSeg[0],listSeg[1]) 
                # nearby but reversed?
                if (v1-s2).length + (v2-s1).length < 0.001:
                    return (listSeg[0],listSeg[1]) 
            return None

        # Find the nearby edge-ends and correct them
        if isects:
            for node,segments in isects:
                for segment in segments:
                    if segment not in segList:
                        nearSeg = nearbySeg(segment,segList)
                        if nearSeg:
                            self.intersectingSegments[nearSeg].append(node)
                            v1,v2 = nearSeg
                    else:
                        self.intersectingSegments[segment].append(node)

        def _inorderExtend(segment, v1, v2, points):
            # Extend a segment <segment> by <points> that are on
            # between the points v1, v2
            k, r = None, False
            if v1[0] < v2[0]:   k, r = lambda i: i[0], True
            elif v1[0] > v2[0]: k, r = lambda i: i[0], False
            elif v1[1] < v2[1]: k, r = lambda i: i[1], True
            else:               k, r = lambda i: i[1], False
            l = [ p for p in sorted(points, key=k, reverse=r) ]
            i = next((i for i, p in enumerate(segment) if p == v2), -1)
            assert(i>=0)
            for e in l:
                # a vertex can appear only once in a segment
                if not e in segment:
                    segment.insert(i, e)

        # Extend the segments by the intersection points, if any.
        for segment,isectNodes in self.intersectingSegments.items():
            pathInSegment = list(segment)
            _inorderExtend(pathInSegment, segment[0], segment[1], isectNodes)
            self.intersectingSegments[segment] = pathInSegment

    # define way-widths from way tags. This part should later be replaced in the generation
    # of the way-segments by the way manager.
    def temporaryWayWidth(self,tags):
        wayWidths = {
            "motorway": 4.,
            "motorway_link": 3.,
            "trunk": 3.,
            "trunk_link": 3.,
            "primary": 3.,
            "primary_link": 3.,
            "secondary": 2.5,
            "secondary_link": 2.5,
            "tertiary": 2.5,
            "tertiary_link": 2.5,
            "unclassified": 2.,
            "residential": 2.5,
            "living_street": 3.,
            "service": 2.,
            "pedestrian": 1.5,
            "track": 1.,
            "escape": 2.,
            "raceway": 2.,
            "other": 1.,
            # "road", # other
            "steps": 2.,
            "footway": 1.5,
            "path": 1.,
            "cycleway": 1.5,
            "bridleway": 1.5           
        }

        if 'highway' in tags:
            if tags['highway'] in wayWidths:
                width = wayWidths[tags['highway']]
                if 'lanes' in tags:
                    width *= int(tags['lanes'])
                return width
        else:
            return 1.   # used as default

    # Creates the network graph <self.sectionNetwork> for way-sctions (ways between crossings)
    def createWaySectionNetwork(self):
        # get border polygon (ounter-clockwise) of scene frame
        minX, minY = self.app.projection.fromGeographic(self.app.minLat, self.app.minLon)
        maxX, maxY = self.app.projection.fromGeographic(self.app.maxLat, self.app.maxLon)

        # prepare clipper for this frame
        clipper = SCClipper(minX,maxX,minY,maxY)

        wayManager = self.app.managersById["ways"]

         # create full way network
        fullNetwork = WayNetwork()

        # some way tags to exclude, used also in findSelfIntersections(),
        # should be moved to defs.
        excludedTags = ['steps']

        for way in wayManager.getAllWays():
            # Exclude ways with unwanted tags
            if [tag for tag in excludedTags if tag in way.element.tags]:
                continue

            # Get the width of the way segment. Later this should be delivered from the
            # way-segement and be different for every category.
            width = self.temporaryWayWidth(way.element.tags)

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
                        netSeg = NetSegment(v1,v2,way.category,(v2-v1).length, width)
                        fullNetwork.addSegment(netSeg)

        borderPolygon = clipper.getPolygon()
        for v1,v2 in zip(borderPolygon[:-1],borderPolygon[1:]):
            netSeg = NetSegment(v1,v2,'scene_border',(v2-v1).length, 0.) # has no width
            fullNetwork.addSegment(netSeg)

        # create way-section network
        self.sectionNetwork = createSectionNetwork(fullNetwork)

    # Gets all polylines and buildings in scene (buildings with shared edges combined
    # to building blocks) and collects them in a list <self.geosPolyList> of PyGeos
    # polynoms. All their vertices are inserted into the KD-tree <self.kdTree>
    # Create mapping <self.vertIndexToPolyIndex> between the index of the vertex and
    # index of the polygon in <self.geosPolyList>
    def createPolylinePolygons(self,manager):

        # eventually polyline tags have to be excluded
        excludedTags = []
        self.geosPolyList = []
        for polyline in manager.polylines:
            if [tag for tag in excludedTags if tag in polyline.element.tags]:
                 continue

            edges = []
            for edge in polyline.edges:
                # Check for edges splitted by self-intersections
                newSegments = self.intersectingSegments.get( (tuple(edge.v1),tuple(edge.v2)), None)
                if newSegments:
                    for v1,v2 in zip(newSegments[:-1],newSegments[1:]):
                        edges.append((v1,v2))
                else:
                    edges.append((edge.v1,edge.v2))
 
            vertList = [Vector(edge[0]) for edge in edges] + [edges[-1][1]]
            if polyline.element.closed and not 'barrier' in polyline.element.tags:
                geosCoords = [self.geosF.createCoordinate(v) for v in vertList]
                geosRing = self.geosF.createLinearRing(geosCoords)
                geosPoly = self.geosF.createPolygon(geosRing)
            else: 
                # Linear objects, like fences, get here a small width to be a polygon
                width = 0.3
                geosCoords = [self.geosF.createCoordinate(v) for v in vertList]
                geosString = self.geosF.createLineString(geosCoords)
                geosPoly = geosString.buffer(width,cap_style=CAP_STYLE.flat)
            self.geosPolyList.append(geosPoly)

        # add building polylines as polygon objects to <polyList>
        self.buildingPolynoms = mergeBuildingsToBlocks(manager.buildings)
        self.geosPolyList.extend(self.buildingPolynoms)

        # create mapping between the index of the vertex and index of the polygon in <geosPolyList>
        self.vertIndexToPolyIndex.extend(
            polyIndex for polyIndex, polygon in enumerate(self.geosPolyList) for _ in range(polygon.numpoints)
        )

        # the total number of vertices
        totalNumVerts = sum(polygon.numpoints for polygon in self.geosPolyList )

        # allocate the memory for an empty numpy array
        self.polyVerts = np.zeros((totalNumVerts, 2))

        # fill vertices in <self.polyVerts>
        index = 0
        for polygon in self.geosPolyList:
            for vert in polygon.coords:
                self.polyVerts[index] = (vert.x,vert.y)
                index += 1

        self.createKdTree()

    # Creates polygons from graph cycles of the section network 
    def createCyclePolygons(self):
        cycleSegs = self.sectionNetwork.getCycles()

        self.cyclePolys = []
        for segList in cycleSegs:
            # Find dead-end ways by searching for segment pairs, where (src,dst) == (dst,src)
            # The cycle has been constructed correctly for these cases, but contains
            # antiparallel edges, which can't be processed by PyGEOS.
            src = [s.s for s in segList]
            dst = [s.t for s in segList]
            deadEndSegs = []
            toRemove = []
            rev_pairs = [pair for pair in zip(dst,src)]
            for indx, pair in enumerate(zip(src,dst)):
                if pair in rev_pairs:
                    rev_indx = rev_pairs.index(pair)
                    if indx < rev_indx:
                        revPath = segList[rev_indx].path.copy()
                        revPath.reverse()
                        if segList[indx].path == revPath:
                            deadEndSegs.append(segList[indx])
                            toRemove.extend([segList[indx],segList[rev_indx]])

            # Remove all eventual dead-ends from segment list and get the boundary of cycle polygon
            boundaryList = [seg for seg in segList if seg not in toRemove]

            if not boundaryList:
                continue

            # Create the cycle boundary polygon
            boundary = [v for s in boundaryList for v in s.path[:-1] ] + [boundaryList[0].s]
            geosCoords = [ self.geosF.createCoordinate(v) for v in boundary ] 
            geosRing = self.geosF.createLinearRing(geosCoords)
            boundaryPoly = self.geosF.createPolygon(geosRing)

            # The boundary of the cycle polygon consists of way segments. Subtract them
            # as polygon with the way as center line and a buffer polygon with a width
            # given by the way-category. Additionally, one of the two antiparallel dead-end
            # way-segments has to be subtracted.
            if deadEndSegs:
                boundaryList.extend(deadEndSegs)

           # Exclude segments from scene border, which have a zero width
            wayPolys = [seg.getBuffer() for seg in deadEndSegs if seg.width > 0.]
            for wayPoly in wayPolys:
                try:
                    boundaryPoly = boundaryPoly.difference(wayPoly)
                except:
                    plotNetwork(self.sectionNetwork)
                    plotGeosPoly(boundaryPoly,False,'b',2)
                    plotGeosPoly(wayPoly,False,'r',2)
                    plotEnd()

            # if polygon has been broken in parts, separate them
            if boundaryPoly.geom_type == 'Polygon':
                self.cyclePolys.append(boundaryPoly)
            else: # Multipolygon
                for geom in boundaryPoly.geoms:
                    self.cyclePolys.append(geom)

    def createWayEnvironmentPolygons(self):
        debug = False
        environmentPolys = []
        for polyNr,cyclePoly in enumerate(self.cyclePolys):

            # Construct a circumscribed circle around the polygon vertices
            # used as search range in KD-tree of polyline objects.
            cycleVerts = [Vector((v.x,v.y)) for v in cyclePoly.exterior.coords]
            center, radius = circumCircle(cycleVerts)

            # Query the KD-tree for indices of polyline objects that are at least
            # partly within the search range.
            queryCycleIndices = set(
                self.makeKdQuery(center, radius)
            )

            # if polyline objects found, subtract them from the cycle polygon <cyclePoly>
            if queryCycleIndices:
                # Get the polyline objects in <objPolys>.
                objPolys = [self.geosPolyList[indx] for indx in queryCycleIndices]

                # sort by size to subtract largest objects first
                objPolys.sort(key=lambda x:x.area,reverse=True)

                for objPoly in objPolys:
                    cyclePoly = cyclePoly.difference(objPoly)

                # def makeBridgedPolygon(geom):
                #     holes = geom.interiors
                #     if holes:
                #         # Split polygon holes by bridges (does not work with PyGeos).
                #         # The PyGeos coordinates have to be converted back to tuples.
                #         polyVerts = [(v.x,v.y) for v in geom.exterior.coords]
                #         holeVerts = []
                #         for hole in holes:
                #             holeVerts.append([(v.x,v.y) for v in hole.coords])
                #         bridgedPolyVerts = splitPolygonHoles(polyVerts,holeVerts)
                #         geosCoords = [self.geosF.createCoordinate(v) for v in bridgedPolyVerts+[bridgedPolyVerts[0]]]
                #         geosRing = self.geosF.createLinearRing(geosCoords)
                #         bridgedPoly = self.geosF.createPolygon(geosRing)
                #         environmentPolys.append(bridgedPoly)
                #     else:
                #         environmentPolys.append(geom)
 
                if cyclePoly.geom_type == 'Polygon':
                    environmentPolys.append(cyclePoly)
                    # makeBridgedPolygon(cyclePoly)
                else: # Multipolygon
                    for geom in cyclePoly.geoms:
                        environmentPolys.append(geom)
                        # makeBridgedPolygon(geom)
            else:
                environmentPolys.append(cyclePoly)

        colorCycler = iterColorCycle()
        for poly in environmentPolys:
            try:
                L = polygonDecomposition(poly)
            except Exception as ex:
                # import traceback
                # traceback.print_exception(type(ex), ex, ex.__traceback__)
                # plotGeosWithHoles(poly,True)
                # plt.title('exception')
                # plotEnd()
                continue
            plotFillMutliPolyList(L,colorCycler)
            # plotEnd()

    def createKdTree(self):
        from scipy.spatial import KDTree
        self.kdTree = KDTree(self.polyVerts)

    def makeKdQuery(self, searchCenter, searchRadius):
        return (
            self.vertIndexToPolyIndex[vertIndex] for vertIndex in self.kdTree.query_ball_point(searchCenter, searchRadius)
        )


# Plotting functions used during development
#-----------------------------------------------------------------------------------


def plotPoly(polygon,vertsOrder,color='k',width=1.,order=100):
    count = 0
    for v1,v2 in zip(polygon[:-1],polygon[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
        if vertsOrder:
            plt.text(v1[0],v1[1],str(count))
        count += 1
        # plt.plot(v1[0],v1[1],'kx')
    v1, v2 = polygon[-1], polygon[0]
    plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
    if vertsOrder:
        plt.text(v1[0],v1[1],str(count))

def plotGeosPoly(geosPoly,vertsOrder,color='k',width=1.,order=100):
    poly = [(c.x,c.y) for c in geosPoly.coords]
    plotPoly(poly,vertsOrder,color,width,order)

def plotGeosWithHoles(geosPoly,vertsOrder,color='k',width=1.,order=100):
    poly = [(c.x,c.y) for c in geosPoly.exterior.coords]
    plotPoly(poly,vertsOrder,color,width,order)
    for ring in geosPoly.interiors:
        p = [(c.x,c.y) for c in ring.coords]
        plotPoly(p,vertsOrder,color,width,order)


def plotGeosPolyFill(geosPoly,vertsOrder,color='k',width=1.,order=100):
    poly = [(v.x,v.y) for v in geosPoly.coords]
    x = [n[0] for n in poly]
    y = [n[1] for n in poly]
    plt.fill(x,y,'#ff0000',alpha = 0.2,zorder = 500)

def plotPolyFill(poly):
    x = [n[0] for n in poly]
    y = [n[1] for n in poly]
    plt.fill(x,y,'#ff0000',alpha = 0.2,zorder = 500)

def plotWaySeg(wayseg,color='k',width=1.,order=100):
    for v1,v2 in zip(wayseg.path[:-1],wayseg.path[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
        plt.plot(v1[0],v1[1],'k.')
        plt.plot(v2[0],v2[1],'k.')
        x = (v1[0]+v2[0])/2
        y = (v1[1]+v2[1])/2
        plt.text(x,y,str(wayseg.ID))

def iterColorCycle():
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    return cycle(colors)
    # import matplotlib.cm as cm
    # import numpy as np
    # colors = cm.prism(np.linspace(0, 1, 50))
    # return cycle(colors)

def plotFillMutliPolyList(polyList,colorCycler):
    for poly in polyList:
        color = next(colorCycler)
        coords = [(v.x,v.y) for v in poly]
        plotPoly(coords,False,color)
        x = [v.x for v in poly]
        y = [v.y for v in poly]
        plt.fill(x,y,color,alpha=1.0,zorder=10)

# def plotRange(poly,holes,color='#ff0000',alpha = 0.7,zorder=2):
#     from lib.CompGeom import patchify
#     polyList = []
#     poly_np = np.array((
#         tuple(p[0] for p in poly),
#         tuple(p[1] for p in poly),
#     ))
#     polyList.append(poly_np)
#     holes_np = []
#     for hole in holes:
#         hole_np = np.array((
#             tuple(p[0] for p in hole),
#             tuple(p[1] for p in hole),
#         ))
#         polyList.append(hole_np)
#     patch = patchify(polyList,color,2,alpha,zorder)
#     plt.gca().add_patch(patch)

def plotNetwork(network):
    for count,seg in enumerate(network.iterAllSegments()):
        plotWaySeg(seg,'k',0.5)

def plotCycle(cycle):
    nodes = [n for s in cycle for n in s.path[:-1]]
    plotPoly(nodes,True,'m',4)

def plotEnd():
    plt.gca().axis('equal')
    plt.show()
