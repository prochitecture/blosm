from mathutils import Vector
import numpy as np
from itertools import *

import matplotlib.pyplot as plt

from way.way_network import WayNetwork, NetSegment
from way.way_algorithms import createSectionNetwork

from lib.pygeos.geom import GeometryFactory
from lib.pygeos.shared import CAP_STYLE
from lib.CompGeom.algorithms import circumCircle, SCClipper
from lib.CompGeom.splitPolygonHoles import splitPolygonHoles

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
        self.geosPolyList = None
        self.polyVerts = None
        self.vertIndexToPolyIndex = []
        self.kdTree = None
        self.cyclePolys = None

    def do(self, manager):
        self.createWaySectionNetwork()
        self.createPolylinePolygons(manager)
        self.createCyclePolygons()
        self.createWayEnvironmentPolygons()

    def cleanup(self):
        self.kdTree = None
        self.bldgVerts = None
        self.vertIndexToPolyIndex.clear()

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

        # some way tags to exclude
        excludedTags = ['steps']

        for way in wayManager.getAllWays():
            if [tag for tag in excludedTags if tag in way.element.tags]:
                continue
            # eliminate crossings to avoid self intersecting cycles in graph
            if 'footway' in way.element.tags and way.element.tags['footway']=='crossing':
                continue

            # width of the way segment. Later this should be delivered from the
            # way-segement and be different for every catgory.
            width = self.temporaryWayWidth(way.element.tags)

            for segment in way.segments:
                way.element.tags
                v1, v2 = Vector(segment.v1),Vector(segment.v2)
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

        def triangleCapPoly(v0,v1,width):
            uVec = v1-v0
            uVec.normalize()
            nVec = Vector((-uVec[1],uVec[0]))
            vm = v0 - uVec*width*0.1
            poly = [
                vm + nVec*width*1.1,
                vm + nVec*width*1.1 + uVec*width*2.5,
                vm + uVec*width*0.2,
                vm - nVec*width*1.1 + uVec*width*3.,
                vm - nVec*width*1.1,
                vm + nVec*width*1.1
            ]
            geosCoords = [self.geosF.createCoordinate(v) for v in poly]
            geosRing = self.geosF.createLinearRing(geosCoords)
            capPoly = self.geosF.createPolygon(geosRing)
            return capPoly

        # eventually polyline tags have to be excluded
        excludedTags = []
        self.geosPolyList = []
        for polyline in manager.polylines:
            if [tag for tag in excludedTags if tag in polyline.element.tags]:
                 continue

            vertList = [Vector((edge.v1[0],edge.v1[1])) for edge in polyline.edges]
            vertList.append(Vector((polyline.edges[-1].v2[0],polyline.edges[-1].v2[1])))
            if polyline.element.closed and not 'barrier' in polyline.element.tags:
                geosCoords = [self.geosF.createCoordinate(v) for v in vertList]
                geosRing = self.geosF.createLinearRing(geosCoords)
                geosPoly = self.geosF.createPolygon(geosRing)
            else: # Linear objects, like fences, get here a small width to be a polygon
                # The end caps are made triangular to avoid unwanted intersections
                width = 0.3
                geosCoords = [self.geosF.createCoordinate(v) for v in vertList]
                geosString = self.geosF.createLineString(geosCoords)
                geosPoly = geosString.buffer(width,cap_style=CAP_STYLE.flat)
                cap = triangleCapPoly(vertList[0],vertList[1], width)
                geosPoly = geosPoly.difference(cap)
                cap = triangleCapPoly(vertList[-1],vertList[-2], width)
                geosPoly = geosPoly.difference(cap)
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
            # antiparallel edges, which can't be processed by PyGeos
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

            # Remove all eventual dead-ends from segment list and get boundary of cycle polygon
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
            wayPolys = [seg.getBuffer() for seg in boundaryList if seg.width > 0.]
            polysToSubtractFrom = [ boundaryPoly ]
            while wayPolys:
                wayPoly = wayPolys.pop()
                partsOfPoly = []
                for currentPoly in polysToSubtractFrom:
                    diffPoly = currentPoly.difference(wayPoly)
                    # if the difference is a single polygon ...
                    if diffPoly.geom_type == 'Polygon':
                        if currentPoly.contains(wayPoly): 
                            # segPoly is a hole in current polygon, try again at end
                            wayPolys.insert(0,wayPoly)
                            partsOfPoly.append(currentPoly)
                        else:
                            # keep the remaining polygon part
                            partsOfPoly.append(diffPoly)
                   # ... else, the polygon has been broken in parts
                    elif diffPoly.geom_type == 'MultiLinearRing':
                        # then, keep them all
                        for poly in diffPoly.geoms:
                            partsOfPoly.append(poly)

                # new list for next way segment to subtract
                polysToSubtractFrom = partsOfPoly

            # the remaing parts are collected as processed cycle polygons
            if polysToSubtractFrom:
                for poly in polysToSubtractFrom:
                    self.cyclePolys.append(poly)


    def createWayEnvironmentPolygons(self):
        debug = False
        environmentPolys = []
        for polyNr,cyclePoly in enumerate(self.cyclePolys):

            # Construct a circumscribed circle around the polygon vertices
            # used as search range in KD-tree of polyline objects.
            cycleVerts = [Vector((v.x,v.y)) for v in cyclePoly.coords]
            center, radius = circumCircle(cycleVerts)
 
            debug = False
            if debug:
                plt.text(center[0],center[1],str(polyNr),zorder=600)
                if polyNr != 158:# 31, 107, 15 fence, forests: 173, 177
                    continue

            # Query the KD-tree for indices of polyline objects that are at least
            # partly within the search range.
            queryCycleIndices = set(
                self.makeKdQuery(center, radius)
            )

            # if polyline objects found, subtract them from the cycle polygon <cyclePoly>
            if queryCycleIndices:
                # Get the polyline objects in <objPolys>.
                objPolys = [self.geosPolyList[indx] for indx in queryCycleIndices]

                # sort by size to subtract largets objects first
                objPolys.sort(key=lambda x:x.area)

                polysToProcess = [ {'poly':cyclePoly,'holes':[]} ]
                while objPolys:
                    objPoly = objPolys.pop()
                    partsOfCyclePoly = []
                    for currentPoly in polysToProcess:
                        if debug:
                            plotGeosPoly(currentPoly['poly'],False,'k')
                            plotGeosPoly(objPoly,False,'r')
                            plotEnd()

                        if objPoly.contains(currentPoly['poly']):               # object covers polynom
                            continue
                        elif currentPoly['poly'].disjoint(objPoly):             # object is completely outside
                            if debug:
                                print('-----> outside')
                            partsOfCyclePoly.append(currentPoly)
                            continue
                        elif currentPoly['poly'].relate(objPoly,'T**FF*FF*'):   # object is completely inside -> hole
                            if debug:
                                print('-----> hole')
                            partsOfCyclePoly.append(currentPoly)
                            currentPoly['holes'].append(objPoly)
                        elif currentPoly['poly'].intersects(objPoly):           # polynoms intersect -> subtraction
                            if debug:
                                print('-----> subtract')
                            diffPoly = currentPoly['poly'].difference(objPoly)
                            # the difference is a single polygon
                            if diffPoly.geom_type == 'Polygon':
                                newPoly = {'poly':diffPoly,'holes':[]}
                                if currentPoly['holes']:
                                    for hole in currentPoly['holes']:
                                        if newPoly['poly'].contains(hole):
                                            newPoly['holes'].append(hole)
                                        elif newPoly['poly'].intersects(hole):
                                            newPoly['poly'] = newPoly['poly'].difference(hole)
                                partsOfCyclePoly.append(newPoly)
                            # or has been broken in parts
                            elif diffPoly.geom_type == 'MultiLinearRing':
                                # we have to distribute the holes to the polygon parts
                                holes = currentPoly['holes']
                                for poly in diffPoly.geoms:
                                    newPoly = {'poly':poly,'holes':[]}
                                    if holes:
                                        newPoly['holes'] = [hole for hole in holes if poly.contains(hole)]
                                    partsOfCyclePoly.append(newPoly)
                    polysToProcess = partsOfCyclePoly

                    if debug:
                        N = len(polysToProcess)
                        plt.close()
                        for i, poly in enumerate(polysToProcess):
                            plt.subplot(1,N,i+1)
                            plotGeosPoly(cyclePoly,False,'k:')
                            plotGeosPoly(objPoly,False,'r:')
                            plotGeosPoly(poly['poly'],False,'b')
                            for hole in poly['holes']:
                                plotGeosPoly(hole,False,'r')
                        plotEnd()

                if polysToProcess:
                    for processedPoly in polysToProcess:

                        # Process now the holes in the polygons, if any
                        if processedPoly['holes']:
                            if len(processedPoly['holes'])>1:
                                # Holes in holes have to be removed first
                                # To find them, we have to check every combination
                                toRemove = []
                                combs = combinations(processedPoly['holes'],2)
                                for polyA,polyB in combs:
                                    if   polyA.contains(polyB): toRemove.append(polyB)
                                    elif polyB.contains(polyA): toRemove.append(polyA)
                                for poly in toRemove:
                                    processedPoly['holes'].remove(poly)

                            # Split polygon holes by bridges (does not work with PyGeos).
                            # The PyGeos coordinates have to be converted back to tuples.
                            polyVerts = [(v.x,v.y) for v in processedPoly['poly'].coords]
                            holeVerts = []
                            for hole in processedPoly['holes']:
                                holeVerts.append([(v.x,v.y) for v in hole.coords])
                            bridgedPolyVerts = splitPolygonHoles(polyVerts,holeVerts)
                            geosCoords = [self.geosF.createCoordinate(v) for v in bridgedPolyVerts+[bridgedPolyVerts[0]]]
                            geosRing = self.geosF.createLinearRing(geosCoords)
                            bridgedPoly = self.geosF.createPolygon(geosRing)
                            environmentPolys.append(bridgedPoly)
                        else:
                            environmentPolys.append(processedPoly['poly'])

            else:
                environmentPolys.append(cyclePoly)
        for poly in environmentPolys:
            plotGeosPolyFill(poly,False)


    def createKdTree(self):
        from scipy.spatial import KDTree
        self.kdTree = KDTree(self.polyVerts)

    def makeKdQuery(self, searchCenter, searchRadius):
        return (
            self.vertIndexToPolyIndex[vertIndex] for vertIndex in self.kdTree.query_ball_point(searchCenter, searchRadius)
        )


# plotting functions for development
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
        x = (v1[0]+v2[0])/2
        y = (v1[1]+v2[1])/2
        # plt.text(x,y,str(wayseg.ID))

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

def plotEnd():
    plt.gca().axis('equal')
    plt.show()
