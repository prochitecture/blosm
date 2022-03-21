from mathutils import Vector
import numpy as np
from itertools import *
from collections import defaultdict

import matplotlib.pyplot as plt

from way.way_network import WayNetwork, NetSection
from way.way_algorithms import createSectionNetwork
from way.way_graph_cycle import GraphCycle

from lib.pygeos.geom import GeometryFactory
from lib.pygeos.shared import TopologyException
from lib.CompGeom.algorithms import circumCircle, SCClipper, repairSimpleSelfIntersection, progress
from lib.CompGeom.GraphBasedAlgos import DisjointSets, ArticulationPoints
from lib.SweepIntersectorLib.SweepIntersector import SweepIntersector
from lib.triangulation.PolygonTriangulation import PolygonTriangulation
from lib.triangulation.Vertex import Vertex
from lib.triangulation.cleaning import cleaningForTriangulation

from defs.road_polygons import ExcludedWayTags, MinDetectionSize, DetectionGridWidth, SharedVertDist

# cyclic iterate over polygons vertices
def _iterEdges(poly):
        p1, p2= tee(poly)
        p2 = islice(cycle(p2), 1, None)
        return zip(p1,p2)

class RoadPolygons:

    def __init__(self):
        self.networkGraph = None
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
        self.checkAndRepairObjectPolys(manager)
        self.fillObjectsInKDTree()
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
                        netSeg = NetSection(v1,v2,way.category,(v2-v1).length)
                        wayManager.networkGraph.addSegment(netSeg,False)

        borderPolygon = clipper.getPolygon()
        for v1,v2 in zip(borderPolygon[:-1],borderPolygon[1:]):
            netSeg = NetSection(v1,v2,'scene_border',(v2-v1).length) 
            wayManager.networkGraph.addSegment(netSeg)

        # create way-section network
        wayManager.sectionNetwork = self.sectionNetwork = createSectionNetwork(wayManager.networkGraph)

    def checkAndRepairObjectPolys(self,manager):
        geosF = GeometryFactory()
        self.geosPolyList = []

        # some helper functions --------------------------------------------------------
        def vertsToGeosPoly(verts):
            # verts is axpected as a list of vertices of a polygon without end point
            coords = [ geosF.createCoordinate(v) for v in verts+[verts[0]] ]
            return geosF.createPolygon( geosF.createLinearRing(coords) )

        def mergePolygonsToGeos(polygons,polyIDs):
            geosList = []
            for polyID in polyIDs:
                geosList.append( vertsToGeosPoly(polygons[polyID]) )
            multiPoly = geosF.createMultiPolygon(geosList)
            try:
                merged = multiPoly.union()
            except (TopologyException) as e:
                print('Problem with merging polys in mergePolygonsToGeos()',polyIDs)
            return merged

        def mergeIntersectingPolygons(polygons, intersectingPolygonIDs, groupIDs, vertsToIds):
            # find connected groups of intersecting polygons
            disjointSets = DisjointSets()
            for IDs in intersectingPolygonIDs.values():
                IDlist = list(IDs)
                for iD in IDlist[1:]:
                    disjointSets.addSegment(IDlist[0],iD)

            # Merge connected group
            geosPolyGroup = []
            for group in disjointSets:
                for polyID in group:
                    wasRepaired,verts = repairSimpleSelfIntersection(polygons[polyID])
                    if wasRepaired:
                        print('Mapping Conflict repaired: Self intersection in object %d'%(polyID))
                    geosPolyGroup.append( vertsToGeosPoly(verts) )
                multiPoly = geosF.createMultiPolygon(geosPolyGroup)
                try:
                    merged = multiPoly.union()
                except (TopologyException) as e:
                    print('Problem with merging polys in mergeIntersectingPolygons()',group)

                for polyID in group:
                    for v in polygons[polyID]:
                        vertsToIds[v].remove(polyID)

                # result comes in first polygon
                verts = [Vector((v.x,v.y)) for v in merged.exterior.coords[:-1]]
                for v in verts:
                    v.freeze()
                    vertsToIds[v].add(group[0])
                polygons[group[0]] = verts
                for poly in group[1:]:
                    groupIDs.remove(poly)
                    polygons[poly] = []
                print('Mapping Conflict repaired: Intersecting objects: ', group)

            return polygons, groupIDs, vertsToIds

        def separateConnectedPolys(polygons,vertsToIds,artPoint,connectedPolys):
            # find the unit vectors of the edges at the corner srtPoint
            # for both polygons
            cornerData = []
            for i,poly in enumerate(connectedPolys):
                verts = polygons[poly]
                indx =  verts.index(artPoint)
                p1 = verts[indx-1]
                p2 = verts[(indx+1)%len(verts)]
                u1 = (p1-artPoint)/(p1-artPoint).length
                u2 = (p2-artPoint)/(p2-artPoint).length
                hypot = (u1-u2).length
                cornerData.append( (hypot,indx,u1,u2) )

            # select sharper corner, which is the one with minimal hypot
            # and shift its artPoint by <SharedVertDist> along bisector
            iP = 0 if cornerData[0][0] < cornerData[1][0] else 1
            _, indx, u1, u2 = cornerData[iP]
            sharedArtPoint = (artPoint + (u1+u2)/(u1+u2).length * SharedVertDist).freeze()
            print('Mapping Conflict repaired: Single shared vertex between objects: ', list(connectedPolys))

            # bookkeeping
            polyID = list(connectedPolys)[iP]
            polygons[polyID][indx] = sharedArtPoint
            vertsToIds[artPoint].remove(polyID)
            vertsToIds[sharedArtPoint].add(polyID)
            return polygons, vertsToIds

        # end of helper functions ------------------------------------------------------

        # Create dictionary polygon id -> vertices from all objects
        polygons = {} 
        for building in manager.buildings:
            verts = [v for v in building.polygon.verts]
            for v in verts: v.freeze()  # make vertices immutable
            polygons[int(building.element.tags["id"])] = [v for v in building.polygon.verts]
        for polyline in manager.polylines:
            if polyline.element.closed and not 'barrier' in polyline.element.tags:
                verts = [edge.v1 for edge in polyline.edges]
                ccw = sum( (v2[0]-v1[0])*(v2[1]+v1[1]) for v1,v2 in _iterEdges(verts)) < 0.
                if not ccw: verts.reverse()
                for v in verts: v.freeze()  # make vertices immutable
                polygons[int(polyline.element.tags["id"])] = verts
                

        # Create a mapping <vertsToIds> from vertex to polygon OSM IDs (there
        # may be multiple IDs for shared vertices).
        # Load class ConnectedComponents with segments
        disjointGroups = DisjointSets()
        vertsToIds = defaultdict(set)
        for osmID,verts in polygons.items():
            for v1,v2 in _iterEdges(verts):
                disjointGroups.addSegment(v1,v2)
                vertsToIds[v1].add(int(osmID))

        # All vertices in <connectedNodes> belong to a connected group of objects.
        for connectedNodes in disjointGroups:
            # find the OSM IDs of the connected objects
            groupIDs = set(osmID for v in connectedNodes for osmID in vertsToIds[v])

            if len(groupIDs) == 1: # Only a single building or polyline object
                verts = [v for v in polygons[next(iter(groupIDs))]]
                wasRepaired,verts = repairSimpleSelfIntersection(verts)     # eventually not required ???
                if wasRepaired:
                    print('Mapping Conflict repaired: Self intersection in object %d'%(next(iter(groupIDs))))
                self.geosPolyList.append( vertsToGeosPoly(verts) )

            else: # we have a group of connected objects.
                # Do we have segment intersections from mapping errors?
                # The class <SweepIntersector> requires segments without duplicates
                # stored in <uniqueSegments>
                uniqueSegments = defaultdict(set)
                for osmID in groupIDs:
                    for v1,v2 in _iterEdges(v for v in polygons[osmID]):
                        v1,v2 = tuple(v1),tuple(v2)
                        if v1 not in uniqueSegments.get(v2,[]):
                            uniqueSegments[v1].add(v2)
                cleanedSegs = [(v1,v2) for v1 in uniqueSegments for v2 in uniqueSegments[v1]]
                intersector = SweepIntersector()
                dictOfIntersections = intersector.findIntersections(cleanedSegs)

                # If there are intersections, try to merge the intersecting polygons
                if dictOfIntersections:
                    intersectingPolygonIDs = defaultdict(set)
                    for seg,isects in dictOfIntersections.items():
                        v = Vector(seg[0]).freeze()
                        isectID = vertsToIds[v]
                        for v in isects[1:-1]:
                            v = Vector(v).freeze()
                            intersectingPolygonIDs[v] = intersectingPolygonIDs[v] | isectID

                    # The merge result is a change of the polygons in <polygons>
                    polygons, groupIDs, vertsToIds = mergeIntersectingPolygons(polygons,intersectingPolygonIDs, groupIDs, vertsToIds)

                articulationPoints = ArticulationPoints()   
                for osmID in groupIDs:
                    for v1,v2 in _iterEdges(v for v in polygons[osmID]):
                        articulationPoints.addEdge(v1,v2)

                # Do we have single shared vertices (articulation points in graph notation)
                for artPoint in articulationPoints:
                    connectedPolys = vertsToIds[artPoint]
                    polygons,vertsToIds = separateConnectedPolys(polygons,vertsToIds,artPoint,connectedPolys)

 
                # At this position, all issues of this group (self-intersections, intersections
                # between objects and single shared vertices) should be repaired.

                # Now remove all shared segments using an adjacency list of the polygons
                # in this group.
                adj = defaultdict(list)
                for osmID in groupIDs:
                    for v1,v2 in _iterEdges(polygons[osmID]):
                        if v1 in adj.get(v2,[]):
                            # segment and its twin detected -> shared segment
                            adj[v2].remove(v1)
                            if not adj[v2]:
                                del adj[v2]
                        else:
                            adj[v1].append(v2)

                if any(len(v)!=1 for v in adj.values()):
                    # This is a very rare case. Due to error propapgation after the
                    # correction of mapping issues, new issues may have been generated.
                    # In this section, we try to repair them. Not very nice, but it works.
                    geosPoly = mergePolygonsToGeos(polygons, groupIDs)
                    if geosPoly.geom_type == 'Polygon':
                        self.geosPolyList.append(geosF.createPolygon(geosPoly.exterior))
                    elif geosPoly.geom_type == 'MultiLinearString' or geosPoly.geom_type == 'MultiLinearRing':
                        for geom in geosPoly.geoms:
                            if geom.geom_type == 'Polygon':
                                self.geosPolyList.append(geosF.createPolygon(geom.exterior))
                            else:
                                # for v1,v2s in adj.items():
                                #     if len(v2s)>1:
                                #         print('inner edge')
                                #     for v2 in v2s:
                                #         plt.plot([v1[0],v2[0]],[v1[1],v2[1]],'b',linewidth = 3, zorder=900)
                                print('There is a problem in following the segment chain with multiple successors')
                    else:
                        print('There is a problem in following the segment chain with multiple successors')
                        # for v1,v2s in adj.items():
                        #     for v2 in v2s:
                        #         plt.plot([v1[0],v2[0]],[v1[1],v2[1]],'b',linewidth = 3, zorder=900)
                    continue


                # The adjacency list <adj> gets now parsed for polygons. Take a starting vertex <firstVert>
                # and follow the segments in the list until this vertex is reached again and create a polygon
                # from its vertices.
                while adj:
                    firstVert = next(iter(adj))
                    verts = [firstVert]
                    thisVert = adj[firstVert].pop().freeze()
                    # maybe the list in <adj[firstVert]> is now empty, so remove this entry
                    if not adj[firstVert]: del adj[firstVert]

                    # now loop until first vertex is found again
                    while thisVert != firstVert and thisVert is not None:
                        verts.append(thisVert)
                        nextVert = adj[thisVert].pop().freeze()
                        if not nextVert:
                            print('There is a problem in following the segment chain')
                        if not adj[thisVert]: del adj[thisVert]
                        thisVert = nextVert

                    coords = [ geosF.createCoordinate(v) for v in verts+[verts[0]] ]
                    geosPoly = geosF.createPolygon(geosF.createLinearRing(coords))
                    self.geosPolyList.append(geosPoly)

    def fillObjectsInKDTree(self):
        # create mapping between the index of the vertex and index of the polygon in <geosPolyList>
        self.vertIndexToPolyIndex.extend(
            polyIndex for polyIndex, polygon in enumerate(self.geosPolyList) for _ in range(polygon.numpoints)
        )

        # the total number of vertices
        ordinaryNumVerts = sum(polygon.numpoints for polygon in self.geosPolyList )

        # Create detection vertices for large object polygons
        detectionPoints = []
        for polyIndex, polygon in enumerate(self.geosPolyList):
            bbox = polygon.envelope
            if bbox.width > MinDetectionSize or bbox.height > MinDetectionSize:
                nW = int(np.round(bbox.width/DetectionGridWidth))
                nH = int(np.round(bbox.height/DetectionGridWidth))
                s = DetectionGridWidth
                grid = [(bbox.minx + (s/2) + (x * s), bbox.miny + (s/2) + (y * s)) for x in range(nW) for y in range(nH)]
                polyV = [(v.x,v.y) for v in polygon.exterior.coords]
                
                # Point in polygon test
                thisPoly = []
                for p in grid:
                    odd = False
                    for v1,v2 in zip(polyV,polyV[1:]):
                        if ( ( (v1[1] > p[1]) != (v2[1] > p[1]) ) and \
                            (p[0] < ( (v2[0] - v1[0]) * (p[1] - v1[1]) / (v2[1] - v1[1]) ) + v1[0]) ):
                            odd = not odd
                    if odd:
                        thisPoly.append(p)
                detectionPoints.extend(thisPoly)
                self.vertIndexToPolyIndex.extend( [polyIndex]*len(thisPoly) )

        detectNumVerts = len(detectionPoints)
        # allocate the memory for an empty numpy array
        self.polyVerts = np.zeros((ordinaryNumVerts+detectNumVerts, 2))

        # fill vertices into <self.polyVerts>
        index = 0
        for polygon in self.geosPolyList:
            for vert in polygon.coords:
                self.polyVerts[index] = (vert.x,vert.y)
                index += 1

        # fill detection vertices into <self.polyVerts>
        for vert in detectionPoints:
            self.polyVerts[index] = vert
            index += 1

        self.createKdTree()

    # Creates polygons from graph cycles of the section network 
    def createCyclePolygons(self):
        cycleSegs, islandSegs, solitaireSegs = self.sectionNetwork.getCycles()
        self.graphCycles = []
        for segList in cycleSegs:
            self.graphCycles.append( GraphCycle(segList) )
        islandPolys = []
        islandSpurs = []
        islandIndxs = []
        islandIndx = 0
        for sectNr,segList in enumerate(islandSegs):
            polys, spurs = GraphCycle.createHoleCycles(segList)
            for _ in polys:
                islandIndxs.append(sectNr)
                islandSpurs.append(spurs)
                islandIndx += 1
            islandPolys.extend( polys )

        # Island graph-cycles are detected because they are in clockwise order.
        # They become holes in their enclosing polygons. This computation is not
        # very effective (O(n*m), where n is the number of graph-cycles and m is 
        # the number of island  graph cycles), but they are rare, so m is small.
        for islandIndx, islandPoly in zip(islandIndxs,islandPolys):
            enclosingPolys = [cycle for cycle in self.graphCycles if cycle.cyclePoly.relate(islandPoly,'T**FF*FF*')]
            if enclosingPolys:
                smallestEnclosingPoly = min(enclosingPolys,key=lambda x:x.cyclePoly.area)
                smallestEnclosingPoly.islandSegs.update(islandSegs[islandIndx])
                smallestEnclosingPoly.spurs.extend(islandSpurs[islandIndx])
                try:
                    smallestEnclosingPoly.cyclePoly = smallestEnclosingPoly.cyclePoly.difference(islandPoly)
                except (TopologyException,ValueError) as e:
                    import traceback
                    traceback.print_exception(type(e), e, e.__traceback__)

    def createWayEnvironmentPolygons(self):
        environmentPolys = []
        print(' ')
        for polyNr,graphCycle in enumerate(self.graphCycles):
            progress(polyNr+1,len(self.graphCycles),'polyline subtraction')
            # print('%d/%d polyline subtraction started'%(polyNr+1,len(self.graphCycles)))
            newSlices = []
            for cycleSlice in GraphCycle.sliceLargeCycles(graphCycle.cyclePoly):
 
                # Construct a circumscribed circle around the polygon vertices
                # used as search range in KD-tree of polyline objects.
                cycleVerts = [Vector((v.x,v.y)) for v in cycleSlice.exterior.coords]
                center, radius = circumCircle(cycleVerts)

                # Query the KD-tree for indices of polyline objects that are at least
                # partly within the search range.
                queryCycleIndices = set(
                    self.makeKdQuery(center, radius)
                )

                # if polyline objects found, subtract them from the cycle polygon <cyclePoly>
                subSlices = [cycleSlice]
                if queryCycleIndices:
                    # Get the polyline objects in <objPolys>.
                    objPolys = [self.geosPolyList[indx] for indx in queryCycleIndices]

                    # sort by size to subtract largest objects first
                    objPolys.sort(key=lambda x:x.area,reverse=True)

                    for i, objPoly in enumerate(objPolys):
                        subPolys = subSlices
                        subSlices = []
                        for subPoly in subPolys:
                            if not subPoly.area:
                                break
                            try:
                                subPoly = subPoly.difference(objPoly)
                            except (TopologyException,ValueError) as e:
                                import traceback
                                traceback.print_exception(type(e), e, e.__traceback__)
                                print('For cyclePoly Nr.: ', polyNr)
                                # plt.subplot(1,2,1)
                                # plotGeosWithHoles(subPoly,True)
                                # plt.gca().axis('equal')
                                # plt.subplot(1,2,2)
                                # plotGeosWithHoles(subPoly,False)
                                # plotGeosWithHoles(objPoly,True,'g',2)
                                # plt.title('exception')
                                # plotEnd()

                            if subPoly.geom_type == 'Polygon':
                                subSlices.append(subPoly)
                            else: # Multipolygon
                                for geom in subPoly.geoms:
                                    subSlices.append(geom)

                newSlices.extend(subSlices)
            graphCycle.tiles = newSlices

        # from patchify import plotGeosPatch
        # for polyNr,graphCycle in enumerate(self.graphCycles):
        #     fig = plt.figure()
        #     ax = plt.Axes(fig, [0., 0., 1., 1.])
        #     ax.set_axis_off()
        #     fig.add_axes(ax)
        #     for cycleSlice in graphCycle.tiles:
        #         plotGeosPatch(cycleSlice,'g','b',2,0.3,500)
        #     plotEnd()



        print(' ')
        triangulation = PolygonTriangulation()
        for polyNr,graphCycle in enumerate(self.graphCycles):
            progress(polyNr+1,len(self.graphCycles),'triangulation')
            for poly in graphCycle.tiles:
                polyVerts, holeVerts = cleaningForTriangulation(poly)
                try:
                    triangles = triangulation.triangulate(polyVerts,holeVerts)
                except Exception as e:
                    # saveData(polyNr,polyVerts,holeVerts)
                    import traceback
                    traceback.print_exception(type(e), e, e.__traceback__)
                    print('Exception cyclePoly Nr.: ', polyNr)
                    # plotPoly(polyVerts,True,'k',1)
                    # for hole in holeVerts:
                    #     plotPoly(hole,True,'r')
                    # printMultiPolyData(poly)
                    # plotEnd()
                    continue

                graphCycle.triangles.extend(triangles)

        print(' ')
        for polyNr,graphCycle in enumerate(self.graphCycles):
            progress(polyNr+1,len(self.graphCycles),'plotting triangles')
            for triangle in graphCycle.triangles:
                plotPoly(triangle,False,'r',0.5)

        print('\nplt.show() called  -->  Result')

    def createKdTree(self):
        from scipy.spatial import KDTree
        self.kdTree = KDTree(self.polyVerts)

    def makeKdQuery(self, searchCenter, searchRadius):
        return (
            self.vertIndexToPolyIndex[vertIndex] for vertIndex in self.kdTree.query_ball_point(searchCenter, searchRadius)
        )


# Plotting functions used during development
#-----------------------------------------------------------------------------------

def saveData(polyNr,poly,holes):
    fid = open('C:/Users/Roger/Desktop/data/poly'+str(polyNr)+'.py','w')
    fid.write('from mathutils import Vector\n\n')
    fid.write('polygon = [')
    for v in poly:
        fid.write('%s, '%(str(v)))
    fid.write(']\n')

    fid.write('holes = [\n')
    for hole in holes:
        fid.write('    [')
        for v in hole:
            fid.write('%s, '%(str(v)))
        fid.write('],\n')
    fid.write(']\n')
    fid.close()



def plotEdge(v1,v2,color='k',arrow=False,width=1,order=100):
    plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
    if arrow:
        v = (v1+v2)/2.
        plt.annotate(
            '',
            xytext = (v1[0],v1[1]),
            xy=(v[0],v[1]),
            arrowprops=dict(color=color, width = 0.25, shrink=0., headwidth=3, headlength=8)
        )

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

def plotLine(line,vertsOrder,color='k',width=1.,order=100):
    count = 0
    for v1,v2 in zip(line[:-1],line[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
        if vertsOrder:
            plt.text(v1[0],v1[1],str(count),fontsize=12)
        count += 1
        # plt.plot(v1[0],v1[1],'kx')
    if vertsOrder:
        plt.text(v1[0],v1[1],str(count),fontsize=12)


def plotGeosPoly(geosPoly,vertsOrder,color='k',width=1.,order=100):
    poly = [(c.x,c.y) for c in geosPoly.coords[:-1]]
    plotPoly(poly,vertsOrder,color,width,order)

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

def plotGeosPolyFill(geosPoly,vertsOrder,color='k',width=1.,order=100):
    poly = [(v.x,v.y) for v in geosPoly.coords[:-1]]
    x = [n[0] for n in poly]
    y = [n[1] for n in poly]
    plt.fill(x,y,'#ff0000',alpha = 0.5,zorder = 500)

def plotGeosPolyFillExterior(geosPoly,vertsOrder,color='k',width=1.,order=100):
    poly = [(v.x,v.y) for v in geosPoly.exterior.coords[:-1]]
    x = [n[0] for n in poly]
    y = [n[1] for n in poly]
    plt.plot(x,y,'b',linewidth=0.5)
    plt.fill(x,y,'#ff0000',alpha = 0.2,zorder = 500)

def polyCenter(geosPoly):
    poly = [(v.x,v.y) for v in geosPoly.coords[:-1]]
    x = sum([n[0] for n in poly])/len(poly)
    y = sum([n[1] for n in poly])/len(poly)
    return (x,y)

def plotPolyFill(poly):
    x = [n[0] for n in poly]
    y = [n[1] for n in poly]
    plt.fill(x,y,'#ff0000',alpha = 0.2,zorder = 500)

def plotWaySeg(wayseg,color='k',width=1.,order=100):
    for v1,v2 in zip(wayseg.path[:-1],wayseg.path[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
        # plt.plot(v1[0],v1[1],'k.')
        # plt.plot(v2[0],v2[1],'k.')
        x = (v1[0]+v2[0])/2
        y = (v1[1]+v2[1])/2
        # plt.text(x,y,str(wayseg.ID),fontsize=12)

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
        # color = next(colorCycler)
        color = "red"
        coords = [(v.x,v.y) for v in poly]
        plotPoly(coords,False,color,width=0.3)
        # x = [v.x for v in poly]
        # y = [v.y for v in poly]
        # plt.fill(x,y,color,alpha=1.0,zorder=10)

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
#     patch = patchify(polyList,color,color,2,alpha,zorder)
#     plt.gca().add_patch(patch)


def plotNetwork(network):
    for count,seg in enumerate(network.iterAllSegments()):
        plt.plot(seg.s[0],seg.s[1],'k.')
        plt.plot(seg.t[0],seg.t[1],'k.')

        plotWaySeg(seg,'k',0.5)

def plotCycle(cycle):
    nodes = [n for s in cycle for n in s.path[:-1]]
    plotPoly(nodes,True,'m',4)

def plotEnd():
    plt.gca().axis('equal')
    plt.show()

def printMultiPolyData(multipoly):
    fid = open('C:/Users/Roger/Desktop/polybug.py','w')
    poly = [(c.x,c.y) for c in multipoly.exterior.coords[:-1]]
    if not multipoly.exterior.is_ccw:
        poly.reverse()
    fid.write('polygon = [')
    for v in poly:
        fid.write('%s, '%(str(v)))
    fid.write(']\n')
    for i,ring in enumerate(multipoly.interiors):
        p = [(c.x,c.y) for c in ring.coords[:-1]]
        if ring.is_ccw:
            p.reverse()
        fid.write('hole%d = ['%(i))
        for v in p:
            fid.write('%s, '%(str(v)))
        fid.write(']\n')
    fid.write('holes = [')
    for i in range(len(multipoly.interiors)):
        fid.write('hole%d[:-1],'%(i))
    fid.write(']\n')    
    fid.close()

def plotGeosPolyFillID(poly,vertsOrder,polyNr,color='k',width=1.,order=100):
    c = polyCenter(poly)
    plt.text(c[0],c[1],str(polyNr))
    plotGeosPolyFill(poly,vertsOrder,color,width,order)
