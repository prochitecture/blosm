from mathutils import Vector
import numpy as np
from itertools import *
from collections import defaultdict

import matplotlib.pyplot as plt

from way.way_network import WayNetwork, NetSegment
from way.way_algorithms import createSectionNetwork
from way.way_graph_cycle import GraphCycle

from lib.pygeos.geom import GeometryFactory
from lib.pygeos.shared import CAP_STYLE, TopologyException
from lib.CompGeom.algorithms import circumCircle, SCClipper, repairSimpleSelfIntersection
from lib.CompGeom.GraphBasedAlgos import DisjointSets, ArticulationPoints
from lib.SweepIntersectorLib.SweepIntersector import SweepIntersector
from lib.triangulation.PolygonTriangulation import PolygonTriangulation
from lib.triangulation.Vertex import Vertex

# cyclic iterate over polygons vertices
def _iterEdges(poly):
        p1, p2= tee(poly)
        p2 = islice(cycle(p2), 1, None)
        return zip(p1,p2)

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
        # self.checkAndRepairObjectPolys_III(manager)
        # self.fillObjectsInKDTree()
        # self.createCyclePolygons()
        # self.createWayEnvironmentPolygons()

    def cleanup(self):
        self.kdTree = None
        self.bldgVerts = None
        self.vertIndexToPolyIndex.clear()

    def findSelfIntersections(self, manager):
        wayManager = self.app.managersById["ways"]
        wayManager.junctions = (
            [],#mainJunctions,
            []#smallJunctions
        )

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

        # Unfortunately, polylines can not be added to the intersection check,
        # as they may be collinear to ways, which is not supported by the
        # Bentley-Ottmann sweep-line algorithm.
        # for polyline in manager.polylines:
        #     for edge in polyline.edges:
        #         v1, v2 = (edge.v1[0],edge.v1[1]),  (edge.v2[0],edge.v2[1])
        #         segList.append((v1,v2))

        intersector = SweepIntersector()
        self.intersectingSegments = intersector.findIntersections(segList)
        # for s,isects in self.intersectingSegments.items():
        #     for p in isects[1:-1]:
        #         plt.plot(p[0],p[1],'r.',zorder=900)

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
        wayManager.networkGraph = WayNetwork()
        # fullNetwork = WayNetwork()

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
                        wayManager.networkGraph.addSegment(netSeg,False)

        borderPolygon = clipper.getPolygon()
        for v1,v2 in zip(borderPolygon[:-1],borderPolygon[1:]):
            netSeg = NetSegment(v1,v2,'scene_border',(v2-v1).length, 0.) # has no width
            wayManager.networkGraph.addSegment(netSeg)

        # create way-section network
        self.sectionNetwork = createSectionNetwork(wayManager.networkGraph)

    def checkAndRepairObjectPolys_III(self,manager):
        geosF = GeometryFactory()
        self.geosPolyList = []


        # some helper functions --------------------------------------------------------
        def vertsToGeosPoly(verts):
            # verts is axpected as a list of vertices of a polygon without end point
            coords = [ geosF.createCoordinate(v) for v in verts+[verts[0]] ]
            return geosF.createPolygon( geosF.createLinearRing(coords) )

        def mergeIntersectingPolygons(polygons, intersectingPolygonIDs):
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
                    print('Problem with merging polys ',group)

                # result comes in first polygon
                verts = [Vector((v.x,v.y)) for v in merged.exterior.coords[:-1]]
                for v in verts: v.freeze()
                polygons[group[0]] = verts
                for poly in group[1:]:
                    polygons[poly] = []
                print('Mapping Conflict repaired: Intersectiong objects: ', group)

            return polygons

        def separateConnectedPolys(polygons,artPoint,connectedPolys):
            # The shared vertex <artPoint> of the first polygon is shifted along 
            # the bisector by 1mm into the polygon ...
            firstPoly = next(iter(connectedPolys))
            verts = polygons[firstPoly]               # vertices of first polygon
            sharedIndx = verts.index(artPoint)        # index of shared vertex
            # find predecessor and successor of shared vertex
            N = len(verts)
            pred, succ = verts[(sharedIndx-1)%N], verts[(sharedIndx+1)%N]
            # compute bisector and move shared vertex for 1mm into the polygon
            u1 = (pred-artPoint)/(pred-artPoint).length
            u2 = (succ-artPoint)/(succ-artPoint).length
            shiftedShared = artPoint + (u1 +u2)/(u1 +u2).length * 0.001
            # freeze and put it back into polygon
            shiftedShared.freeze()
            polygons[firstPoly][sharedIndx] = shiftedShared
            return polygons
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

        # All vertices in connectedNodes belong to a connected group of objects.
        for connectedNodes in disjointGroups:
            # find the OSM IDs of the connected objects
            groupIDs = set(osmID for v in connectedNodes for osmID in vertsToIds[v])

            if len(groupIDs) == 1: # a single building or polyline object
                verts = [v for v in polygons[next(iter(groupIDs))]]
                wasRepaired,verts = repairSimpleSelfIntersection(verts)     # eventually not required ???
                if wasRepaired:
                    print('Mapping Conflict repaired: Self intersection in object %d'%(next(iter(groupIDs))))
                self.geosPolyList.append( vertsToGeosPoly(verts) )

            else: # we have a group of connected objects.
                plt.close()
                plt.subplot(1,3,1)
                for osmID in groupIDs:
                    plotPoly(polygons[osmID],False)
                plt.title('original')
                plt.gca().axis('equal')
                problem = False
                # plotEnd()

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
                    polygons = mergeIntersectingPolygons(polygons,intersectingPolygonIDs)
                    plt.subplot(1,3,2)
                    for osmID in groupIDs:
                        if polygons[osmID]:
                            plotPoly(polygons[osmID],False)
                    plt.title('intersection')
                    plt.gca().axis('equal')
                    problem = True
                    # plotEnd()

                articulationPoints = ArticulationPoints()   
                for osmID in groupIDs:
                    for v1,v2 in _iterEdges(v for v in polygons[osmID]):
                        articulationPoints.addEdge(v1,v2)

                # Do we have single shared vertices (articulation points in graph notation)
                for artPoint in articulationPoints:
                    connectedPolys = vertsToIds[artPoint]
                    polygons = separateConnectedPolys(polygons,artPoint,connectedPolys)

                    plt.subplot(1,3,3)
                    for osmID in groupIDs:
                        if polygons[osmID]:
                            plotPoly(polygons[osmID],False)
                    plt.title('articulation')
                    plt.gca().axis('equal')
                    problem = True

                if problem:
                    plotEnd()
 
                # import pickle
                # with open('D:/LEDA/SweepIntersector/problem.pkl', 'wb') as f:
                #     pickleList = []
                #     for seg in segList:
                #         v1, v2 = seg[0], seg[1]
                #         newseg = ((v1[0],v1[1]),(v2[0],v2[1]))
                #         pickleList.append(newseg)
                #     pickle.dump(pickleList, f)

                # if problem:
                #     for v1,v2 in subSet:
                #         plt.plot([v1[0],v2[0]],[v1[1],v2[1]],'k',linewidth=1,zorder=100)
                #     if ap:
                #         for p in ap:
                #             plt.plot(p[0],p[1],'bo',markersize=10)
                #     plotEnd()  
                # test=1

        test = 1
                # geosList = []
                # if len(polyIDs) > 1:
                #     segList = []
                #     intersector = SweepIntersector()
                #     for polyID in polyIDs:
                #         verts = [v for v in polygons[polyID]]
                #         # plotPoly(verts,False)
                #         for v1,v2 in _iterEdges(verts):
                #             segList.append((tuple(v1),tuple(v2)))
                #     intersectingSegments = intersector.findIntersections(segList)
                #     if intersectingSegments:
                #         for v1, v2 in segList:
                #             plt.plot([v1[0],v2[0]],[v1[1],v2[1]],'k',linewidth=1,zorder=100)
                #             plt.plot(v1[0],v1[1],'b.')
                #             plt.plot(v2[0],v2[1],'b.')
                #         for s,isects in intersectingSegments.items():
                #             for p in isects[1:-1]:
                #                 plt.plot(p[0],p[1],'ro',markersize=8)

                #         plotEnd()

        test=1

    def checkAndRepairObjectPolys_II(self,manager):
        import queue
        def BFSUtil(u, adj, visited):
            subnet = []       
            # Create a queue for BFS 
            q = queue.Queue()            
            # Mark the current node as visited
            # and enqueue it 
            visited[u] = True
            q.put(u)            
            # 'i' will be used to get all adjacent 
            # vertices 4 of a vertex list<int>::iterator i         
            while(not q.empty()):              
                # Dequeue a vertex from queue 
                # and print it 
                u = q.queue[0] 
                subnet.append(u)
                q.get() 
                # Get all adjacent vertices of the 
                # dequeued vertex s. If an adjacent 
                # has not been visited, then mark 
                # it visited and enqueue it 
                i = 0
                while i != len(adj[u]):
                    if (not visited[adj[u][i]]):
                            visited[adj[u][i]] = True
                            q.put(adj[u][i])
                    i += 1
            return subnet

        geosF = GeometryFactory()

        # Create dictionary polygon id -> vertices
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

        # for p,verts in polygons.items():
        #     plotPoly(verts,False)

        # Create an adjacency list <edges> of all edges. Edges that have a twin
        # are shared edges and get removed.
        # Create also a mapping <vertsToIds> from vertex to polygon IDs (there
        # may be multiple ids for shared vertices)
        self.geosPolyList = []

        adj = defaultdict(list)
        vertsToIds = defaultdict(set)
        for polyID,verts in polygons.items():
            for source,target in _iterEdges(verts):
                if source in adj.get(target,[]):
                    # edge and its twin detected -> shared edge
                    adj[target].remove(source)
                    if not adj[target]:
                        del adj[target]
                else:
                    adj[source].append(target)
                    vertsToIds[source].add(int(polyID))

        # for source,targets in adj.items():
        #     if len(set(targets)) > 1 and  max(len(vertsToIds[t]) for t in targets) == 1:
        #         for t in targets:
        #             print(vertsToIds[t])
        #             plt.plot([source[0],t[0]],[source[1],t[1]],'r',linewidth=4)
        #         for polyID in vertsToIds[source]:
        #             plotPoly(polygons[polyID],False)
        #         plt.plot(source[0],source[1],'ro',markersize=8)
        #         # plt.text(source[0],source[1],len(test),fontsize=12)
        #         plotEnd()

        for source,targets in adj.items():
            if len(set(targets)) > 1 and  max(len(vertsToIds[t]) for t in targets) == 1:
                id = list(vertsToIds[source])[0]
                verts = polygons[id]    
                indx = verts.index(source) 
                N = len(verts)
                pred, succ = verts[(indx-1)%N], verts[(indx+1)%N]
                u1 = (pred-source)/(pred-source).length
                u2 = (succ-source)/(succ-source).length
                shiftedShared = source + (u1 +u2)/(u1 +u2).length * 0.001
                shiftedShared.freeze()
                verts[indx] = shiftedShared
                # replace new shared vertex in polygon, maybe it will be considered again
                polygons[id][indx] = shiftedShared

        adj = defaultdict(list)
        vertsToIds = defaultdict(set)
        for polyID,verts in polygons.items():
            for source,target in _iterEdges(verts):
                adj[source].append(target)
                vertsToIds[source].add(int(polyID))

        # shared = [k for k,v in adj.items() for t in v if k in adj.get(t)]
        # for p in shared:
        #     polyIDs = set()
        #     for v in adj[p]:
        #         polyIDs = polyIDs.union(vertsToIds[v])
        #     # if len(polyIDs) == 2:
        #     plt.plot(p[0],p[1],'ro',markersize=8)
        #     plt.text(p[0],p[1],len(polyIDs),fontsize=12)
        #     print(len(polyIDs))
        # plotEnd()

        visited = {k:False for k,v in adj.items()} 
        for u in adj.keys():
            if (visited[u] == False): 
                subnet = BFSUtil(u, adj, visited)
                polyIDs = set()
                for v in subnet:
                    polyIDs = polyIDs.union(vertsToIds[v])
                geosList = []
                for polyID in polyIDs:
                    verts = [v for v in polygons[polyID]]
                    coord = [ geosF.createCoordinate(v) for v in verts+[verts[0]] ]
                    geosList.append( geosF.createPolygon( geosF.createLinearRing(coord)) )
                multiPoly = geosF.createMultiPolygon(geosList)
                try:
                    unified = multiPoly.union()
                except (TopologyException) as e:
                    plotPoly(subnet,False)
                    plt.title('Exception')
                    plotEnd()

                if unified.geom_type == 'Polygon':
                    self.geosPolyList.append( geosF.createPolygon(unified.exterior) )
                    exterior = [Vector((v.x,v.y)) for v in unified.exterior.coords[:-1]]
                    # plotPoly(exterior,False,'b',1)
                    # if unified.interiors:
                    #     test=1
                else:
                    plt.subplot(1,2,1)
                    shared = []
                    for i1,i2 in combinations(range(len(unified.geoms)), 2):
                        g1 = unified.geoms[i1].exterior.coords[:-1]
                        g2 = unified.geoms[i2].exterior.coords[:-1]
                        sharedVertList = list(set(g1).intersection(g2))
                        if sharedVertList:
                            shared.append([i1,i2,sharedVertList[0]])
                    if shared:
                        test = shared[0][2]
                        plt.plot(test.x,test.y,'ro')

                    for geom in unified.geoms:
                        self.geosPolyList.append( geosF.createPolygon(geom.exterior) )
                        exterior = [Vector((v.x,v.y)) for v in geom.exterior.coords[:-1]]
                        plotPoly(exterior,False,'k',1)
                    plt.gca().set_title(str(len(unified.geoms)))
                    plt.gca().axis('equal')
                    plt.subplot(1,2,2)
                    for polyID in polyIDs:
                        verts = [v for v in polygons[polyID]]
                        plotPoly(verts,False,'k')
                    plotEnd()
                # test = 1

        for poly in self.geosPolyList:
            plotGeosWithHoles(poly,False,'b',1)
        # plotEnd()



    # Gets all polylines and buildings in scene (those with shared edges combined
    # to blocks) and collects them in a list <self.geosPolyList> of PyGeos polynoms. 
    def checkAndRepairObjectPolys_I(self,manager):

        # Create dictionary polygon id -> vertices
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

        # Create an adjacency list <edges> of all edges. Edges that have a twin
        # are shared edges and get removed.
        # Create also a mapping <vertsToIds> from vertex to polygon IDs (there
        # may be multiple ids for shared vertices)
        edges = defaultdict(list)
        vertsToIds = defaultdict(set)
        for polyID,verts in polygons.items():
            for source,target in _iterEdges(verts):
                if source in edges.get(target,[]):
                    # edge and its twin detected -> shared edge
                    edges[target].remove(source)
                    if not edges[target]:
                        del edges[target]
                else:
                    edges[source].append(target)
                    vertsToIds[source].add(int(polyID))

        # # DEBUG-PLOT, to be removed when no more required
        # # Shows edges after shared edge removal.
        # for source,targets in edges.items():
        #     for target in targets:
        #         plotEdge(source,target)

        # Detect IDs of conflicting polygons due to shared vertices
        conflictingPolys = defaultdict(list)
        for source,targets in edges.items():
            if len(targets) > 1:
                conflictingPolys[tuple(vertsToIds[source])].append(source)

        # # DEBUG-PLOT, to be removed when no more required
        # # Shows conflicting polygons and shared vertices.
        # for polys,sources in conflictingPolys.items():
        #     for source in sources:
        #         plt.plot(source[0],source[1],'co',markersize=10,zorder=500)
        #     for poly in polys:
        #         plotPoly(polygons[poly],False,'r',3,500)
        # plotEnd()

        # Now try to repair the conflicting polygons
        geosF = GeometryFactory()
        for polyIDs,sharedVerts in conflictingPolys.items():
            # Remove all edges of these polygons from <edges>, but save those
            # that have been already detected as shared edges and are no more
            # in the adjacency list <edges> in <hiddenSharedEdges>.
            hiddenSharedEdges = []
            for polyID in polyIDs:
                for source,target in _iterEdges(polygons[polyID]):
                    if source in edges and target in edges[source]:
                        # remove those that exist and are not shared
                        edges[source].remove(target)
                        if not edges[source]: del edges[source]
                    else:
                        # store edges detected as shared
                        hiddenSharedEdges.append((source,target))

            # We will analyse this conflict using pyGEOS,so create pyGEOS polygons           
            polysToBeAnalyzed = []
            for polyID in polyIDs:
                verts = [v for v in polygons[polyID]]
                coords = [ geosF.createCoordinate(v) for v in verts+[verts[0]] ]
                polysToBeAnalyzed.append( geosF.createPolygon( geosF.createLinearRing(coords)) )

            # Try to detect the reasons of the conflict
            conflictReasons = []
            for p1,p2 in combinations(range(len(polysToBeAnalyzed)), 2):
                intersection = polysToBeAnalyzed[p1].intersection(polysToBeAnalyzed[p2]).area > 0.
                relation = polysToBeAnalyzed[p1].relate(polysToBeAnalyzed[p2])
                singleVertex = relation.matrix[1][1] == 0
                touchEdge = relation.matrix[1][1] == 1
                if intersection or touchEdge:
                    conflictReasons.append('intersect')
                elif singleVertex:
                    conflictReasons.append('singleVert')
                else:
                    conflictReasons.append('unknown')

            # Try to repair the conflict
            if len(conflictReasons) == 1:   # Only two polygons involved in conflict
                id0, id1 = list(polyIDs)
                if conflictReasons[0] == 'intersect':
                    # # DEBUG-PLOT, to be removed when no more required
                    # plt.close()
                    # plotGeosWithHoles(polysToBeAnalyzed[0],False,'r',2)
                    # plotGeosWithHoles(polysToBeAnalyzed[1],False,'g',2)
                    # plt.title('Intersection')
                    # plotEnd()
                    # Merge the intersecting polygons
                    joinedPoly = polysToBeAnalyzed[0].union(polysToBeAnalyzed[1])
                    if joinedPoly.geom_type=='Polygon' and joinedPoly.exterior.is_valid:
                        verts = [Vector((v.x,v.y)) for v in joinedPoly.exterior.coords[:-1]]
                        if not joinedPoly.exterior.is_ccw:
                            verts.reverse()
                        for v in verts: v.freeze()  # make vertices immutable
                        # put edges of joined polygon back to <edges>, when not shared.
                        for source,target in _iterEdges(verts):
                            if (source,target) not in hiddenSharedEdges:
                                edges[source].append(target)
                        # replace new vertices in both polygons, maybe they will be considered again
                        polygons[id0] = verts
                        polygons[id1] = verts
                        print('Mapping CONFLICT repaired: Overlapping or touching polygons %s and %s'%(id0,id1))
                    else:
                        print('Mapping CONFLICT unresolved: Polygons %s and %s removed'%(id0,id1))
                elif conflictReasons[0] == 'singleVert':
                    # # DEBUG-PLOT, to be removed when no more required
                    # plt.close()
                    # plotGeosWithHoles(polysToBeAnalyzed[0],False,'r',2)
                    # plotGeosWithHoles(polysToBeAnalyzed[1],False,'g',2)
                    # plt.title('Single Vertex')
                    # plotEnd()
                    # The first polygon remains as it was
                    for source,target in _iterEdges(polygons[id0]):
                        source.freeze(), target.freeze()
                        if (source,target) not in hiddenSharedEdges:
                            edges[source].append(target)
                    # The shared vertex of the second polygon is shifted along 
                    # the bisector by 1mm into the polygon ...
                    sharedVert = sharedVerts[0]
                    verts = polygons[id1]    
                    indx = verts.index(sharedVert) 
                    N = len(verts)
                    pred, succ = verts[(indx-1)%N], verts[(indx+1)%N]
                    u1 = (pred-sharedVert)/(pred-sharedVert).length
                    u2 = (succ-sharedVert)/(succ-sharedVert).length
                    shiftedShared = sharedVert + (u1 +u2)/(u1 +u2).length * 0.001
                    shiftedShared.freeze()
                    verts[indx] = shiftedShared
                    # replace new shared vertex in polygon, maybe it will be considered again
                    polygons[id1][indx] = shiftedShared

                    # ... and then put back to <edges>, when not shared.
                    for source,target in _iterEdges(verts):
                        if (source,target) not in hiddenSharedEdges:
                            edges[source].append(target)
                    print('Mapping CONFLICT repaired: Single shared vertex between polygons %s and %s'%(id0,id1))
                else:
                    id0, id1 = list(polyIDs)
                    print('Mapping CONFLICT unresolved: Polygons %s and %s removed'%(id0,id1))

            else: # More than two polygons involved in conflict
                # # DEBUG-PLOT, to be removed when no more required
                # plt.close()
                # for poly in polysToBeAnalyzed:
                #     plotGeosWithHoles(poly,False,'r',2)
                # plt.title('Multiple poly conflict')
                # plotEnd()

                if all(reason=='intersect' for reason in conflictReasons):
                    # Merge polygons involved in this conflict
                    polysToBeMerged = []
                    for polyID in polyIDs:
                        verts = [v for v in polygons[polyID]]
                        coords = [ geosF.createCoordinate(v) for v in verts+[verts[0]] ]
                        polysToBeMerged.append( geosF.createPolygon( geosF.createLinearRing(coords)) )
                    joinedPoly = polysToBeMerged[0]
                    for poly in polysToBeMerged[1:]:
                        joinedPoly = joinedPoly.union(poly)
                    if joinedPoly.geom_type=='Polygon' and joinedPoly.exterior.is_valid:
                        verts = [Vector((v.x,v.y)) for v in joinedPoly.exterior.coords[:-1]]
                        if not joinedPoly.exterior.is_ccw:
                            verts.reverse()
                        for v in verts: v.freeze()  # make vertices immutable
                        # put edges of joined polygon back to <edges>, when not shared.
                        for source,target in _iterEdges(verts):
                            if (source,target) not in hiddenSharedEdges:
                                edges[source].append(target)
                        print('Mapping CONFLICT repaired: Overlapping or touching polygons ', polyIDs, ' merged.')

        # # DEBUG-PLOT, to be removed when no more required
        # # Shows edges after conflict repair.
        # for source,targets in edges.items():
        #     for target in targets:
        #         plotEdge(source,target)
        # plotEnd()

        # All conflicts due to shared vertices should be repaired now. The adjacency list <edges>
        # is now parsed for polygons. Take a starting vertex <firstVert> and follow the edges in
        # the list until this edge is reached again and create a polygon from its vertices.
        self.geosPolyList = []
        while edges:
            firstVert = next(iter(edges))
            vertList = [firstVert]
            thisVert = edges[firstVert].pop().freeze() #segDict.get(firstVert).popleft().freeze()
            if not edges[firstVert]: del edges[firstVert]

            while thisVert != firstVert and thisVert is not None:
                vertList.append(thisVert)
                nextQueue = edges.get(thisVert)
                if nextQueue is None:
                    # Merge polygons involved in this conflict
                    involvedPolygonsIDs = set()
                    for v in vertList:
                        involvedPolygonsIDs = involvedPolygonsIDs.union(vertsToIds[v])
                    # remove their edges from edge list
                    for polyID in involvedPolygonsIDs:
                        for source,target in _iterEdges(polygons[polyID]):
                            if source in edges and target in edges[source]:
                                # remove those that exist
                                edges[source].remove(target)
                                if not edges[source]: del edges[source]

                    polysToBeMerged = []
                    for polyID in involvedPolygonsIDs:
                        verts = [v for v in polygons[polyID]]
                        plotPoly(verts,False,'r')
                        coords = [ geosF.createCoordinate(v) for v in verts+[verts[0]] ]
                        polysToBeMerged.append( geosF.createPolygon( geosF.createLinearRing(coords)) )
                    joinedPoly = polysToBeMerged[0]
                    for poly in polysToBeMerged[1:]:
                        joinedPoly = joinedPoly.union(poly)
                    if joinedPoly.geom_type=='Polygon' and joinedPoly.exterior.is_valid:
                        self.geosPolyList.append( joinedPoly )
                    print('Mapping CONFLICT repaired:: Polygons ', involvedPolygonsIDs, ' merged')
                    thisVert = None
                    break
                nextVert = nextQueue.pop().freeze()
                if not edges[thisVert]: del edges[thisVert]
                thisVert = nextVert

            if len(vertList) > 2:
                # simple check for self-intersections
                selfIsect, isectVerts = simpleSelfIntersection(vertList)
                if selfIsect:
                    # try to find the involved polygons
                    polyIDs = set()
                    for v in vertList:
                        polyIDs = polyIDs.union(vertsToIds[v])
                    conflictIDs = set()
                    for v in isectVerts:
                        conflictIDs = polyIDs.union(vertsToIds[v])
                    if len(polyIDs) == 2:
                        # we're lucky, two polygons. Do they intersect?
                        geosPolys = []
                        for polyID in polyIDs:
                            verts = [v for v in polygons[polyID]]
                            coord = [ geosF.createCoordinate(v) for v in verts+[verts[0]] ]
                            geosPolys.append( geosF.createPolygon( geosF.createLinearRing(coord)) )
                        # haveIntersection
                        doIntersect = geosPolys[0].intersection(geosPolys[1]).area > 0.
                        doTouch = geosPolys[0].relate(geosPolys[1],'212111212')
                        id0, id1 = list(polyIDs)
                        if doIntersect or doTouch:
                            joinedPoly = geosPolys[0].union(geosPolys[1])
                            if joinedPoly.geom_type=='Polygon' and joinedPoly.exterior.is_valid:
                                self.geosPolyList.append( joinedPoly )
                                
                                print('Mapping CONFLICT repaired: Overlapping or touching polygons %s and %s'%(id0,id1))
                            else:
                                print('Mapping CONFLICT unresolved: Polygons %s and %s removed'%(id0,id1))
                        else:
                            print('Mapping CONFLICT unresolved: Polygons %s and %s removed'%(id0,id1))
                    else:
                        print('Mapping CONFLICT repaired: Conflicting objects', conflictIDs)
                        geosList = []
                        for polyID in polyIDs:
                            verts = [v for v in polygons[polyID]]
                            coord = [ geosF.createCoordinate(v) for v in verts+[verts[0]] ]
                            geosList.append( geosF.createPolygon( geosF.createLinearRing(coord)) )
                        multiPoly = geosF.createMultiPolygon(geosList)
                        try:
                            unified = multiPoly.union()
                        except (TopologyException) as e:
                            plotPoly(vertList,False)
                            plt.title('Exception')
                            plotEnd()


                        # # DEBUG-PLOT, to be removed when no more required
                        # plt.close()
                        # plt.subplot(1,2,1)
                        # plotPoly(vertList,False)
                        # plt.gca().axis('equal')
                        # plt.subplot(1,2,2)
                        # # for polyID in polyIDs:
                        # #     verts = [v for v in polygons[polyID]]
                        # #     coord = [ geosF.createCoordinate(v) for v in verts+[verts[0]] ]
                        # #     debugPoly = geosF.createPolygon( geosF.createLinearRing(coord)) 
                        # plotGeosWithHoles(unified,False,'k',1)
                        # plt.title('Multiple poly conflict')
                        # plotEnd()
                else:
                    geosCoords = [geosF.createCoordinate(v) for v in vertList+[vertList[0]]]
                    geosRing = geosF.createLinearRing(geosCoords)
                    geosPoly = geosF.createPolygon(geosRing)
                    self.geosPolyList.append( geosPoly )

        # # DEBUG-PLOT, to be removed when no more required
        # # Show the final result of checkAndRepairObjectPolys()
        for poly in self.geosPolyList:
            plotGeosWithHoles(poly,False,'b',2)
        plotEnd()

    def fillObjectsInKDTree(self):
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
        cycleSegs, holeSegs, _ = self.sectionNetwork.getCycles()
        self.graphCycles = []
        for segList in cycleSegs:
            self.graphCycles.append( GraphCycle(segList) )
        holePolys = []
        for segList in holeSegs:
            holePolys.extend( GraphCycle.createHoleCycles(segList) )

        # Hole graph-cycles are detected because they are in clockwise order.
        # They become holes in their enclosing polygons. This computation is not
        # very effective (O(n*m), where n is the nmber of graph-cycles and m is 
        # the number of hole graph cycles), but they are rare, so m is small.
        for holePoly in holePolys:
            plotGeosWithHoles(holePoly,False,'r')                             # contains properly -> 'T**FF*FF*'
            enclosingPolys = [cycle for cycle in self.graphCycles if cycle.cyclePoly.relate(holePoly,'T**FF*FF*')]
            smallestEnclosingPoly = min(enclosingPolys,key=lambda x:x.cyclePoly.area)
            try:
                smallestEnclosingPoly.cyclePoly = smallestEnclosingPoly.cyclePoly.difference(holePoly)
            except (TopologyException,ValueError) as e:
                import traceback
                traceback.print_exception(type(e), e, e.__traceback__)

    def createWayEnvironmentPolygons(self):
        environmentPolys = []
        for polyNr,graphCycle in enumerate(self.graphCycles):
            print('%d/%d polyline subtraction started'%(polyNr+1,len(self.graphCycles)))
            cyclePoly = graphCycle.cyclePoly

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
                    if not cyclePoly.area:
                        break
                    try:
                        cyclePoly = cyclePoly.difference(objPoly)
                    except (TopologyException,ValueError) as e:
                        import traceback
                        traceback.print_exception(type(e), e, e.__traceback__)
                        print('For cyclePoly Nr.: ', polyNr)
                        # plt.subplot(1,2,1)
                        # plotGeosWithHoles(cyclePoly,True)
                        # plt.gca().axis('equal')
                        # plt.subplot(1,2,2)
                        # plotGeosWithHoles(cyclePoly,False)
                        # plotGeosWithHoles(objPoly,True,'g',2)
                        # plt.title('exception')
                        # plotEnd()

                if cyclePoly.geom_type == 'Polygon':
                    graphCycle.subPolys.append(cyclePoly)
                else: # Multipolygon
                    for geom in cyclePoly.geoms:
                        graphCycle.subPolys.append(geom)
            else:
                graphCycle.subPolys.append(cyclePoly)

        triangulation = PolygonTriangulation()
        for polyNr,graphCycle in enumerate(self.graphCycles):
            for polyNr,poly in enumerate(graphCycle.subPolys):
                # if polyNr==40:
                # c = polyCenter(poly)
                # plt.text(c[0],c[1],str(polyNr))
                # plotGeosWithHoles(poly,False,'k',1,900)
                # plotGeosPolyFillExterior(poly,False)
                holes = poly.interiors
                # the exterior polygon must be counter-clockwise
                polyVerts = [Vertex((v.x,v.y)) for v in poly.exterior.coords[:-1]]
                if not poly.exterior.is_ccw:
                    polyVerts.reverse()
                holeVerts = []
                for hole in holes:
                    holeV = [Vertex((v.x,v.y)) for v in hole.coords[:-1]]
                    if hole.is_ccw:
                        holeV.reverse()
                    holeVerts.append(holeV)

                try:
                    triangles = triangulation.triangulate(polyVerts,holeVerts)
                except Exception as e:
                    import traceback
                    traceback.print_exception(type(e), e, e.__traceback__)
                    print('For cyclePoly Nr.: ', polyNr)
                    plotGeosWithHoles(poly,True,'b',2)
                    # printMultiPolyData(poly)
                    # plotEnd()

                graphCycle.triangles.extend(triangles)

        # for polyNr,graphCycle in enumerate(self.graphCycles):
        #     for triangle in graphCycle.triangles:
        #         # plotPolyFill(triangle)
        #         plotPoly(triangle,False,'r',0.5)
        #     print('%d triangulated'%(polyNr+1))
        # # plotEnd()

    def createKdTree(self):
        from scipy.spatial import KDTree
        self.kdTree = KDTree(self.polyVerts)

    def makeKdQuery(self, searchCenter, searchRadius):
        return (
            self.vertIndexToPolyIndex[vertIndex] for vertIndex in self.kdTree.query_ball_point(searchCenter, searchRadius)
        )


# Plotting functions used during development
#-----------------------------------------------------------------------------------

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
        plt.plot(v1[0],v1[1],'kx')
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
        p = [(c.x,c.y) for c in ring.coords]
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
    plt.fill(x,y,'#ff0000',alpha = 0.2,zorder = 500)

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
    plt.fill(x,y,'#0000ff',alpha = 0.2,zorder = 500)

def plotWaySeg(wayseg,color='k',width=1.,order=100):
    for v1,v2 in zip(wayseg.path[:-1],wayseg.path[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
        plt.plot(v1[0],v1[1],'k.')
        plt.plot(v2[0],v2[1],'k.')
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
