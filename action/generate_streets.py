from collections import defaultdict
from itertools import tee, permutations
from mathutils import Vector
import re

from app import AppType

from defs.road_polygons import ExcludedWayTags
from defs.way_cluster_params import minTemplateLength, minNeighborLength, searchDist, dbScanDist

from way.item import Intersection, IntConnector, Section, Street, SideLane, SymLane
from way.item.bundle import Bundle, StreetGroup, mergePseudoMinors, removeSplittingStreets, orderHeadTail, \
                                    findInnerStreets, canBeMerged, joinBundles, mergeBundles,intersectBundles, \
                                    endBundleIntersection, parallelToSample
from way.way_network import WayNetwork, NetSection
from way.way_algorithms import createSectionNetwork
from way.way_properties import lanePattern
from way.streetRouter import streetsOnRoute

from lib.SweepIntersectorLib.SweepIntersector import SweepIntersector

from lib.CompGeom.StaticSpatialIndex import StaticSpatialIndex, BBox
from lib.CompGeom.algorithms import SCClipper
from lib.CompGeom.GraphBasedAlgos import DisjointSets
from lib.CompGeom.PolyLine import PolyLine
from lib.CompGeom.LinePolygonClipper import LinePolygonClipper
from lib.CompGeom.dbscan import dbClusterScan
from lib.CompGeom.chains import find_all_lines
from lib.CompGeom.ConvexHull2D import ConvexHull2D
from lib.CompGeom.centerline import pointInPolygon


# helper functions -----------------------------------------------
def pairs(iterable):
    # iterable -> (p0,p1), (p1,p2), (p2, p3), ...
    p1, p2 = tee(iterable)
    next(p2, None)
    return zip(p1,p2)

def isEdgy(polyline):
    vu = polyline.unitVectors()
    for v1,v2 in pairs(vu):
        if abs(v1.cross(v2)) > 0.6:
            return True
    return False

def _pseudoangle(d):
    p = d[0]/(abs(d[0])+abs(d[1])) # -1 .. 1 increasing with x
    return 3 + p if d[1] < 0 else 1 - p
# ----------------------------------------------------------------


class StreetGenerator():

    def __init__(self, styleStore, getStyle, leftHandTraffic=True,
            debugCircularStreets = False, debugParallelStreets = False, debugBundle=False, plotStreetID=None
        ):
        self.styleStore = styleStore
        self.getStyle = getStyle
        self.leftHandTraffic = leftHandTraffic
        
        #
        # Debug configuration
        #
        
        # debug for self.circularStreets
        self.debugCircularStreets = debugCircularStreets
        # If self.debugParallel is True, all groups of streets considered as parallel are plotted at once.
        self.debugParallelStreets = debugParallelStreets
        # If self.debugBundle is True, the heads and tails and the inner street of every Bundle are plotted.
        self.debugBundle = debugBundle
        # If ID is given: plots the street[ID] with identifier ID in createStreets()
        # This may help to localize a buggy street.
        self.plotStreetID = plotStreetID
        #
        # End of debug configuration
        #

        self.networkGraph = None
        self.sectionNetwork = None
        self.parallelStreets = None

        self.internalTransitionSideLanes = dict()
        self.internalTransitionSymLanes = dict()
        self.intersections = dict()
        self.processedNodes = set()

        self.prohibitedHulls = []
        self.curvedTramLines = []

        # Determine whether all ways are to be calculated or only those intended
        # for vehicles. If True: wayManager.getAllWays() are used, else those from
        # wayManager.getAllVehicleWays()
        self.allWays = True

    def do(self, manager):

        # Finally, the results of StreetGenerator.do() are stored by the WayManager.
        # Link the local attributes to them.
        self.wayManager = manager
        self.waymap = manager.waymap
        self.majorIntersections = manager.majorIntersections
        self.minorIntersections = manager.minorIntersections
        self.transitionSideLanes = manager.transitionSideLanes
        self.transitionSymLanes = manager.transitionSymLanes
        self.streets = manager.streets
        self.bundles = manager.bundles
        self.wayClusters = manager.wayClusters
        self.waySectionLines = manager.waySectionLines

        self.allSplittingStreets = []

        # This class variable of NetSection is not reset with a new instance
        # of StreetGenerator! So we do it here.
        NetSection.ID = 0

        # Find the way segments, that intersect each other.
        self.findSelfIntersections()

        # Split evetually crossing segments.
        # Join the OSM segments to sections between the intersection points.
        # Store them in the network graph <self.sectionNetwork>.
        self.createWaySectionNetwork()

        # Create instances of class Section sections, including style block parameters
        # and attributes. Cretae instances of class Street, containing only a single section.
        # Add these Street instances as edges to the waymap (class WayMap), construct Instances of class
        # Intersection and add them as nodes to the waymap..
        self.createWaymap()

        # Check all instances of Street for SymLane or SideLane intersections. Construct them,
        # create longer Street objects with both streets and the SymSideLane object. Replace
        # these two streets and the intersection by the new instance of Street.
        self.createSymSideLanes()

        # Construct the connectors off the intersections
        self.updateIntersections()

        # Create the dictionary <manager.streets> of Street instances.
        self.createStreets()

        # Find and plot possibly circular streets like roundabouts (curently not active)
        # self.circularStreets()

        # Create an instance self.parallelStreets of <DisjointSets> with
        # groups of streets, that are considered as parallel.
        self.createParallelStreets()

        # Constructs the Bundles and includes their inner streets
        self.createBundles()

        # If two bundles meet and there are no major inner streets at the meeting
        # intersection, they can be merged into one bundle.
        self.joinBundles()

        # Finally, the intersections of bundles are constructed.
        self.createBundleIntersections()


    # Find the way segments, that intersect each other.
    def findSelfIntersections(self):
        uniqueSegments = defaultdict(set)

        # Determine whether all ways are to be calculated or
        # only those intended for vehicles.
        if self.allWays:
            getWays = self.wayManager.getAllWays()
        else:
            getWays = self.wayManager.getAllVehicleWays()

        for way in getWays:
            # ExcludedWayTags is defined in <defs>.
            # It is used also in createWaySectionNetwork().
            if [tag for tag in ExcludedWayTags if tag in way.category]:
                continue
            for segment in way.segments:
                v1, v2 = (segment.v1[0],segment.v1[1]),  (segment.v2[0],segment.v2[1])
                if v1 not in uniqueSegments.get(v2,[]):
                    uniqueSegments[v1].add(v2)
        cleanedSegs = [(v1,v2) for v1 in uniqueSegments for v2 in uniqueSegments[v1]]

        # SweepIntersector provides a fast implementation of a sweep intersection algorithm.
        intersector = SweepIntersector()
        self.intersectingSegments = intersector.findIntersections(cleanedSegs)


    # Creates the network graph <self.sectionNetwork> for way-sections (ways between crossings).
    # If there are intersections, detected by findSelfIntersections(), split them at the
    # intersection point and add the parts to the network.
    # Join the OSM segments to sections between the intersection points.
    def createWaySectionNetwork(self):
        wayManager = self.wayManager

        # prepare clipper for this frame of the scene
        clipper = SCClipper(self.app.minX, self.app.maxX, self.app.minY, self.app.maxY)

        # create full way network
        wayManager.networkGraph = self.networkGraph = WayNetwork(self.leftHandTraffic)

        # Determine whether all ways are to be calculated or
        # only those intended for vehicles.
        if self.allWays:
            getWays = self.wayManager.getAllWays()
        else:
            getWays = self.wayManager.getAllVehicleWays()

        # some way tags to exclude, used also in createWaySectionNetwork(),
        # ExcludedWayTags is defined in <defs>.
        for way in getWays:
            # Exclude ways with unwanted tags
            if [tag for tag in ExcludedWayTags if tag in way.category]:
                continue

            for waySegment in way.segments:
                # If waySSegment has an intersection, detected by findSelfIntersections(), split it at the
                # intersection point and add the parts to the network.
                segments = []
                newSegments = self.intersectingSegments.get( (tuple(waySegment.v1),tuple(waySegment.v2)), None)
                if newSegments:
                    for v1,v2 in zip(newSegments[:-1],newSegments[1:]):
                        segments.append((v1,v2))
                else:
                    segments.append((waySegment.v1,waySegment.v2))

                # If a segment is clipped by the scene's frame, add only the inner part.
                for segment in segments:
                    v1, v2 = Vector(segment[0]),Vector(segment[1])
                    accepted, v1, v2 = clipper.clip(v1,v2)
                    if accepted:
                        netSeg = NetSection(v1,v2,way.category,way.element.tags,(v2-v1).length)
                        wayManager.networkGraph.addSegment(netSeg,False)

        # As a helper for further algorithms, add virtual segments of the category
        # 'scene_border' along the border polygon of the scene.
        borderPolygon = clipper.getPolygon()
        for v1,v2 in zip(borderPolygon[:-1],borderPolygon[1:]):
            netSeg = NetSection(v1,v2,'scene_border',None, (v2-v1).length)
            wayManager.networkGraph.addSegment(netSeg)

        # Join the OSM segments to sections between the intersection points.
        wayManager.sectionNetwork = self.sectionNetwork = createSectionNetwork(wayManager.networkGraph,self.leftHandTraffic)

    # Create instances of class Section sections, including style block parameters
    # and attributes. Cretae instances of class Street, containing only a single section.
    # Add these Street instances as edges to the waymap (class WayMap), construct Instances of class
    # Intersection and add them as nodes to the waymap..
    def createWaymap(self):
        for net_section in self.sectionNetwork.iterAllForwardSegments():
            if net_section.category != 'scene_border':

                # Create Section from net-section, including style block parameters
                section = Section(net_section,PolyLine(net_section.path),self.sectionNetwork)
                oneway = 'oneway' in section.tags and section.tags['oneway'] != 'no'
                street = Street(section.src, section.dst)
                street.insertEnd(section)
                streetStyle = self.styleStore.get( self.getStyle(street) )
                street.style = streetStyle
                street.setStyleBlockFromTop(streetStyle)
                section.street = street

                # Derive Section attributes
                if oneway:
                    totalNumLanesOneway = street.getStyleBlockAttr("totalNumLanesOneway")
                    nrLanes = totalNumLanesOneway if totalNumLanesOneway else street.getStyleBlockAttr("totalNumLanes")
                else:
                    nrLanes = street.getStyleBlockAttr("totalNumLanes")
                props = {
                    'nrLanes' : nrLanes,
                    'laneWidth' : street.getStyleBlockAttr("laneWidth")
                }
                _,fwdPattern,bwdPattern,bothLanes = lanePattern(section.category,section.tags,self.leftHandTraffic,props)
                section.setSectionAttributes(oneway, fwdPattern, bwdPattern, bothLanes, props)

                self.waymap.addNode(Intersection(section.src))
                self.waymap.addNode(Intersection(section.dst))
                street.street = street # Fill superclass Item
                self.waymap.addEdge(street)

        # Add ways to intersections
        for location, intersection in self.waymap.iterNodes(Intersection):
            inStreets, outStreets = self.waymap.getInOutEdges(location)
            intersection.update(inStreets, outStreets)


    # Check all instances of Street for SymLane or SideLane intersections. Construct them,
    # create longer Street objects with both streets and the SymSideLane object. Replace
    # these two streets and the intersection by the new instance of Street.
    def createSymSideLanes(self):
        # Local method to determine, if one of the ends or both of the street
        # is possibly a SymLane or a SideLane. Find their directions and
        # if there are turn lanes.
        def findSymSideIsects(street):
            srcIsectObj = self.waymap.getNode(street.src)
            if srcIsectObj and isinstance(srcIsectObj['object'], Intersection):
                srcIsect = srcIsectObj['object']
                if srcIsect:
                    if srcIsect.order==2:   # Only intersections with two ways can be SymLanes odr SideLanes
                        street0 = srcIsect.leaveWays[0].section.street
                        street1 = srcIsect.leaveWays[1].section.street
                        arriving, leaving = (street0, street1) if street0.dst == street1.src else (street1, street0)
                        # Check directions
                        if arriving.dst == srcIsect.location and leaving.src == srcIsect.location:
                            hasTurns = bool( re.search(r'[^N]', leaving.head.lanePatterns[0] ) )
                            srcIsect = {'isect':srcIsect, 'arriving':arriving, 'leaving':leaving, 'hasTurns':hasTurns}
                        else:
                            srcIsect = None
                    else:
                        srcIsect = None
            else:
                srcIsect = None

            dstIsectObj = self.waymap.getNode(street.dst)
            if dstIsectObj and isinstance(dstIsectObj['object'], Intersection):
                dstIsect = dstIsectObj['object']
                if dstIsect:
                    if dstIsect.order==2:   # Only intersections with two ways can be SymLanes odr SideLanes
                        street0 = dstIsect.leaveWays[0].section.street
                        street1 = dstIsect.leaveWays[1].section.street
                        arriving, leaving = (street0, street1) if street0.dst == street1.src else (street1, street0)
                        # Check directions
                        if arriving.dst == dstIsect.location and leaving.src == dstIsect.location:
                            hasTurns = bool( re.search(r'[^N]', leaving.head.lanePatterns[0] ) )
                            dstIsect = {'isect':dstIsect, 'arriving':arriving, 'leaving':leaving, 'hasTurns':hasTurns}
                        else:
                            dstIsect = None
                    else:
                        dstIsect = None
            else:
                dstIsect = None

            return srcIsect, dstIsect

        processedStreets = set()
        nodesToRemove = []
        longStreets = []
        for _, _, _, street in self.waymap.edges(data='object',keys=True):
            if street in processedStreets:
                continue

            # Find possible SymLane or a SideLane at source or destination
            # of this street.
            srcIsectInit, dstIsectInit = findSymSideIsects(street)
            hasSymSide = srcIsectInit or dstIsectInit

            if hasSymSide:
                # Prepare a new instance of Street, which will store both steets of the SymSide intersection
                longStreet = Street(street.src,street.dst)
                longStreet.insertStreetEnd(street)
                processedStreets.add(street)
                longStreet.pred = street.pred

                if srcIsectInit: # SymSide intersection at the front of this street
                    # If there are turns, it's a SideLane, else it's a SymLane
                    # Create it and insert it into the corresponding list.
                    if srcIsectInit['hasTurns']:
                        newLane = SideLane(srcIsectInit['isect'].location, srcIsectInit['arriving'].head, srcIsectInit['leaving'].head)
                        self.transitionSideLanes.append(newLane)
                    else:
                        newLane = SymLane(srcIsectInit['isect'].location, srcIsectInit['arriving'].head, srcIsectInit['leaving'].head)
                        self.transitionSymLanes.append(newLane)
                    newLane.street = longStreet
                    longStreet.insertFront(newLane)   # insert new SymSide object
                    nodesToRemove.append(srcIsectInit['isect'])
                    srcIsectCurr = srcIsectInit
                    while True: # Continue, if there are more SymSide intersections t the arriving end of this intersection
                        prevStreet = srcIsectCurr['arriving']
                        if prevStreet in processedStreets:
                            break
                        longStreet.insertStreetFront(prevStreet)
                        processedStreets.add(prevStreet)
                        srcIsectCurr, _ = findSymSideIsects(prevStreet)
                        if not srcIsectCurr:
                            break
                        if srcIsectCurr['hasTurns']:
                            newLane = SideLane(srcIsectCurr['isect'].location, srcIsectCurr['arriving'].head, srcIsectCurr['leaving'].head)
                            self.transitionSideLanes.append(newLane)
                        else:
                            newLane = SymLane(srcIsectCurr['isect'].location, srcIsectCurr['arriving'].head, srcIsectCurr['leaving'].head)
                            self.transitionSymLanes.append(newLane)
                        newLane.street = longStreet
                        longStreet.insertFront(newLane)   # insert minor intersection object
                        nodesToRemove.append(srcIsectInit['isect'])

                if dstIsectInit: # SymSide intersection at the end of this street
                    # If there are turns, it's a SideLane, else it's a SymLane
                    # Create it and insert it into the corresponding list.
                    if dstIsectInit['hasTurns']:
                        newLane = SideLane(dstIsectInit['isect'].location, dstIsectInit['arriving'].head, dstIsectInit['leaving'].head)
                        self.transitionSideLanes.append(newLane)
                    else:
                        newLane = SymLane(dstIsectInit['isect'].location, dstIsectInit['arriving'].head, dstIsectInit['leaving'].head)
                        self.transitionSymLanes.append(newLane)
                    newLane.street = longStreet
                    longStreet.insertEnd(newLane)   # insert new lane object
                    nodesToRemove.append(dstIsectInit['isect'])
                    dstIsectCurr = dstIsectInit
                    while True: # Continue, if there are more SymSide intersections at the leaving end of this intersection
                        nextStreet = dstIsectCurr['leaving']
                        if nextStreet in processedStreets:
                            break
                        longStreet.insertStreetEnd(nextStreet)
                        processedStreets.add(nextStreet)

                        dstIsectCurr, _ = findSymSideIsects(nextStreet)
                        if not dstIsectCurr:
                            break
                        if dstIsectCurr['hasTurns']:
                            newLane = SideLane(dstIsectCurr['isect'].location, dstIsectCurr['arriving'].head, dstIsectCurr['leaving'].head)
                            self.transitionSideLanes.append(newLane)
                        else:
                            newLane = SymLane(dstIsectCurr['isect'].location, dstIsectCurr['arriving'].head, dstIsectCurr['leaving'].head)
                            self.transitionSymLanes.append(newLane)
                        newLane.street = longStreet
                        longStreet.insertEnd(newLane)   # insert minor intersection object
                        nodesToRemove.append(dstIsectInit['isect'])

                longStreets.append(longStreet)

        # At this stage, intersections do not yet have connectors. They will get them
        # when processIntersection() is called (which occurs at self.updateIntersections).
        # But the leaving ways structure of the intersections at the end of the new
        # long Street needs to be updated.
        streetIDs = [s.id for s in processedStreets]
        for longStreet in longStreets:
            node = self.waymap.getNode(longStreet.src)
            if node:
                srcIsect = node['object']
                for way in srcIsect.leaveWays:
                    if way.street.id in streetIDs:
                        way.street = longStreet
            node = self.waymap.getNode(longStreet.dst)
            if node:
                dstIsect = node['object']
                for way in dstIsect.leaveWays:
                    if way.street.id in streetIDs:
                        way.street = longStreet

        for node in nodesToRemove:
            self.waymap.removeNode(node.location)

        for longStreet in longStreets:
            self.waymap.addEdge(longStreet)

    # Construct the connectors of the intersections
    def updateIntersections(self):
        # At this stage, it is assumed, that SideLanes, SymLanes are already built
        # and that their nodes are stored in <self.processedNodes>.
        for location, intersection in self.waymap.iterNodes(Intersection):
            if location in self.processedNodes:
                continue

            # Construct the connectors
            intersection.processIntersection()
            if intersection.isMinorIntersection():
                intersection.transformToMinor()
                # see https://github.com/prochitecture/blosm/issues/106#issuecomment-2305297075
                if intersection.isMinor:
                    self.minorIntersections[intersection.location] = intersection
                else:
                    if intersection.order > 1:
                        self.majorIntersections[intersection.location] = intersection
            else:
                if intersection.order > 1:
                    self.majorIntersections[intersection.location] = intersection

            self.processedNodes.add(location)

    # Create the dictionary <manager.streets> of Street instances.
    def createStreets(self):
        # Create instances of Street from the current waymap, by concatenated
        # Streets, internally separated by minor intersections. These are intersections
        # with two major ways and one or more minor ways (of categories 'footway',
        # 'cycleway' or'service').
        # The concatenation is done in WayManager
        for street in self.wayManager.iterStreetsFromWaymap():
            self.streets[street.id] = street

        # Debug plot of street[self.plotStreetID].
        if self.plotStreetID and self.app.type == AppType.commandLine:
            from debug import plt, plotQualifiedNetwork, randomColor, plotStreet, plotEnd

            plotStreet(self.streets[self.plotStreetID],'red')
            p = self.streets[self.plotStreetID].dst
            plt.plot(p[0],p[1],'co',markersize=12,alpha=0.4)
            plotQualifiedNetwork(self.sectionNetwork)
            plotEnd()

    # Find and plot possibly circular streets like roundabouts
    def circularStreets(self):
        # Get tags of the street
        def tagsOfStreet(street):
            for item in street.iterItems():
                if isinstance(item, Section):
                    return item.tags

        # Get the centerline of the wholestreet
        def centerlineOfStreet(street):
            # Find the centerline of the whole street.
            centerlineVerts = []
            for item in street.iterItems():
                if isinstance(item, Section):
                    centerlineVerts.extend( item.centerline)

            # Remove duplicate verices and create polyLine
            centerlineVerts = list(dict.fromkeys(centerlineVerts))
            centerline = PolyLine(centerlineVerts)
            return centerline, centerlineVerts

        # Using tags, create list of streets that are possibly circular.
        circularStreets = []
        for street in self.wayManager.iterStreets():
            tags = tagsOfStreet(street)
            possibleCircular = 'junction' in tags and tags['junction']=='circular'
            possibleCircular = possibleCircular or ('junction' in tags and tags['junction']=='roundabout')

            if possibleCircular:
                circularStreets.append(street)

        # Debug plot of possibly circular streets
        if self.debugCircularStreets and self.app.type == AppType.commandLine:
            from debug import plt, plotQualifiedNetwork, randomColor, plotEnd

            plotQualifiedNetwork(self.sectionNetwork,False)
            for street in circularStreets:
                centerline,verts = centerlineOfStreet(street)
                centerline.plot('red',2,'solid')
                centerline.plotWithArrows('red',1,0.5,'solid',False,950)
            plotEnd()

    # Create an instance self.parallelStreets of <DisjointSets> with
    # groups of streets, that are considered as parallel.
    def createParallelStreets(self):

        def categoryOfStreet(street):
            for item in street.iterItems():
                if isinstance(item, Section):
                    break
            return item.category

        def centerlineOfStreet(street):
            # Find the centerline of the whole street.
            centerlineVerts = []
            for item in street.iterItems():
                if isinstance(item, Section):
                    centerlineVerts.extend( item.centerline)

            # Remove duplicates and create polyLine
            centerlineVerts = list(dict.fromkeys(centerlineVerts))
            centerline = PolyLine(centerlineVerts)
            return centerline, centerlineVerts
        
        debugParallelStreets = self.debugParallelStreets and self.app.type == AppType.commandLine

        # Spatial index (R-tree) of candidate Streets
        candidateIndex = StaticSpatialIndex()

        # Dictionary from index in candidateIndex to street.
        index2Street = dict()

        # The bounding boxes of the streets. The dictionary key is <dictKey>.
        boxes = dict()

        # Some computed attributes of the streets. The dictionary key is <dictKey>.
        attributes = dict()

        # Add the bounding boxes of all streets to the index.
        for street in self.wayManager.iterStreets():
            # Some categories are excluded.
            category =  categoryOfStreet(street)
            if category in ('steps', 'footway', 'cycleway', 'path', 'service', 'runway', 'taxiway'):
                continue

            # Find the centerline of the whole street.
            centerline, centerlineVerts = centerlineOfStreet(street)

            # Exclude if too curvy
            ds = (centerline[0]-centerline[-1]).length / centerline.length()
            if isEdgy(centerline) and ds < 0.9:
                continue

            if category=='tram' and centerline.localCurvature(1) > 0.01:
                self.curvedTramLines.append(street)
                continue

            # Exclude if too short
            if centerline.length() < min(minTemplateLength,minNeighborLength):
                continue

            # Find bounding box and fill in index
            min_x = min(v[0] for v in centerlineVerts)
            min_y = min(v[1] for v in centerlineVerts)
            max_x = max(v[0] for v in centerlineVerts)
            max_y = max(v[1] for v in centerlineVerts)
            bbox = BBox(None,min_x,min_y,max_x,max_y)
            index = candidateIndex.add(min_x,min_y,max_x,max_y)
            index2Street[index] = street
            bbox.index = index
            boxes[street] = (min_x,min_y,max_x,max_y)
            attributes[street] = ( category, centerline, centerlineVerts )

        # Finalize the index for usage.
        candidateIndex.finish()

        # Sometimes, there ar gaps left between the curved tram lines.
        # Find them and add them to the curved tram lines.
        if self.curvedTramLines:
            tramLineSegs = []
            d_max = 5.0
            for tramLine0, tramLine1 in permutations(self.curvedTramLines, 2):
                if tramLine0 != tramLine1  and 0. < (tramLine0.dst-tramLine1.src).length < d_max:
                    tramLineSegs += streetsOnRoute(tramLine0, tramLine1, 2.*d_max)
            self.curvedTramLines.extend(tramLineSegs)

        # This is the structure we use to collect the parallel streets
        self.parallelStreets = DisjointSets()

        # Every street that was inserted into the spatial index becomes now
        # as sample. We expand it to a buffer area around it.
        for sampleStreet in self.wayManager.iterStreets():
            # Use only accepted streets
            if sampleStreet in boxes:
                if sampleStreet in self.curvedTramLines:
                    continue

                sampleCategory, sampleCenterline, _ = attributes[sampleStreet]
                # Create buffer polygon around the sample street with a width according
                # to the category of the sample street.
                bufferWidth = searchDist[sampleCategory]
                bufferPoly = sampleCenterline.buffer(bufferWidth,bufferWidth)

                # Create a line clipper using this polygon.
                clipper = LinePolygonClipper(bufferPoly.verts)

                # Get neighbors of this sample street from the static spatial index, using its
                # bounding box, expanded by the buffer width as additional search range.
                min_x,min_y,max_x,max_y = boxes[sampleStreet]
                results = stack = []
                neighborIndices = candidateIndex.query(min_x-bufferWidth,min_y-bufferWidth,
                                               max_x+bufferWidth,max_y+bufferWidth,results,stack)

                # Now test all these neighbors of the sample street for parallelism.
                for neigborIndex in neighborIndices:
                    neighborStreet = index2Street[neigborIndex]
                    if neighborStreet == sampleStreet:
                        continue # Skip, the sample street is its own neighbor.
                    if neighborStreet in self.curvedTramLines:
                        continue

                    if parallelToSample(sampleStreet, neighborStreet, clipper):
                        self.parallelStreets.addSegment(sampleStreet,neighborStreet)

        # DEBUG: Show clusters of parallel way-sections.
        if debugParallelStreets:
            from debug import plt, plotQualifiedNetwork, randomColor, plotEnd

            inBundles = False

            if not inBundles:
                plotQualifiedNetwork(self.sectionNetwork,False)
            colorIter = randomColor(19)
            for bIndx,streets in enumerate(self.parallelStreets):
                # if bIndx not in [0,1]:
                #     continue
                if inBundles:
                    plotQualifiedNetwork(self.sectionNetwork,False)
                    plt.title("Bundle "+str(bIndx))
                color = next(colorIter)
                allVerts = []
                for street in streets:
                    srcVec, _ = street.endVectors()
                    vu = srcVec/srcVec.length * 3
                    p = street.src
                    plt.text(p[0]+vu[0],p[1]+vu[1],'   S'+str(street.id),fontsize=10,color=color,zorder=960)
                    width = 2
                    if inBundles:
                        color = "red"
                        width = 3
                    centerline,verts = centerlineOfStreet(street)
                    allVerts.extend(verts)
                    centerline.plot(color,width,'solid')
                    centerline.plotWithArrows(color,1,0.5,'solid',False,950)
                    if inBundles:
                        plt.scatter(centerline[0][0], centerline[0][1], s=80, facecolors='none', edgecolors='g',zorder=999)
                        plt.scatter(centerline[-1][0], centerline[-1][1], s=80, facecolors='none', edgecolors='g',zorder=999)
                        # plt.plot(polyline[0][0], polyline[0][1], 'go', markersize=8,zorder=999)
                        # plt.plot(polyline[-1][0], polyline[-1][1], 'go', markersize=8,zorder=999)
                center = sum(allVerts,Vector((0,0)))/len(allVerts)
                plt.text(center[0],center[1],str(bIndx),fontsize=10)
                if inBundles:
                    plt.title('Parallel streets')
                    plotEnd()

            if not inBundles:
                plt.title('Parallel streets')
                plotEnd()
            # END DEBUG

    # Constructs the Bundles and includes their inner streets
    def createBundles(self):
        debugBundle = self.debugBundle and self.app.type == AppType.commandLine

        if debugBundle:
            from debug import plt, plotQualifiedNetwork, randomColor, plotEnd
            colorIter = randomColor(19)


        # The disjoint set <self.parallelStreets> is not a convenient container for
        # the streets group. Mainly, modified streets cannot be stored there.
        streetGroups = []
        for streetGroup in self.parallelStreets:
            streetGroups.append(StreetGroup(streetGroup))

#----------------------------------------------------------------------------------------------
# This part is rarely used. A use case is the scene osm_extracts/streets/toronto_queen_west.osm

        # If there is a cluster of intersections within a streetGroup, the group has to be split
        # into multiple groups and the cluster will later become a bundle intersection.
        # newGroups = []
        # toRemoveGroups = []
        # for streetGroup in streetGroups:
        #     allDistinctIntersections = set()
        #     for street in streetGroup:
        #         allDistinctIntersections |= set([(street.src,street), (street.dst,street)]) 
        #     # Find clusters of intersections
        #     clusters = dbClusterScan(list(allDistinctIntersections), dbScanDist, 2)

        #     # If a cluster is large, eventually split the group into multiple groups
        #     for cluster in clusters:
        #         if len(cluster) > 15:
        #             toRemoveGroups.append(streetGroup)

        #             # Remove all streets that are inside the convex hull of the
        #             # intersections in the cluster
        #             clusterIntersections = [p[0].freeze() for p in cluster]
        #             hullModel = ConvexHull2D()
        #             hull = hullModel(clusterIntersections)
        #             self.prohibitedHulls.append([h.freeze() for h in hull])
        #             toRemoveStreets = [ street for street in streetGroup
        #                                     if pointInPolygon(hull, street.src) in ["IN", "ON"] 
        #                                     and pointInPolygon(hull, street.dst) in ["IN", "ON"] ]
        #             for street in toRemoveStreets:
        #                 streetGroup.remove(street)

        #             # Any other streets and intersections within this hull have to be removed too
        #             for group in streetGroups:
        #                 if group != streetGroup:
        #                     toRemoveStreets = [ street for street in group
        #                                         if pointInPolygon(hull, street.src) in ["IN", "ON"] 
        #                                         and pointInPolygon(hull, street.dst) in ["IN", "ON"] ]
        #                     if toRemoveStreets:
        #                         toRemoveGroups.append(group)

        #             # Now collect all remaining streets into new groups, that remained parallel
        #             processed = set()
        #             streetGroup.sort() # Start wit long streets to get good templates
        #             while len(processed) != len(streetGroup):
        #                 for i, sampleStreet in enumerate(streetGroup):
        #                     if sampleStreet not in processed:
        #                         processed.add(sampleStreet)

        #                         # Create a new street group, with all streets parallel to the sample street
        #                         newGroup = StreetGroup([sampleStreet])
        #                         for neighborStreet in streetGroup[(i+1):]:
        #                             if neighborStreet not in processed: 
        #                                 if parallelToSample(sampleStreet, neighborStreet):
        #                                     processed.add(neighborStreet)
        #                                     newGroup.append(neighborStreet)

        #                         # Collect this group for later appending to streetGroups
        #                         newGroups.append( newGroup )

        # # Bookkeeping for cluster of intersections within a streetGroup
        # if toRemoveGroups:
        #     for group in toRemoveGroups:
        #         if group in streetGroups:
        #             streetGroups.remove(group)

        # if newGroups:
        #     streetGroups.extend(newGroups)
#----------------------------------------------------------------------------------------------

        # Find intersections between streets of different groups.
        # The dictionary key is the location of the intersection
        # and its lists contains the indices of the groups that
        # intersect there.
        groupIntersections = defaultdict(list)
        for streetGroup in streetGroups:
            for street in streetGroup:
                for p in [street.src, street.dst]:
                    groupIntersections[p].append(streetGroup.id)

        for streetGroup in streetGroups:
            # see https://github.com/prochitecture/blosm/issues/104#issuecomment-2322836476
            # Major intersections in bundles, with only one side street, are merged into a long street,
            # similar to minor intersections.
            mergePseudoMinors(self, streetGroup, groupIntersections)


        newGroups = []
        delGroups = []
        for streetGroup in streetGroups:
            # Sometimes, createParallelStreets() is not stopped by intersections with
            # other bundles and delivers streets, that pass these intersections. In
            # such a situation, its  proposed group has to be split into two groups, while
            # the splitting streets at the intersection have to be removed.
            wasSplit, splittedGroups, splittingStreets = removeSplittingStreets(self,streetGroup.id,streetGroup,groupIntersections)
            self.allSplittingStreets.extend(splittingStreets)
            if debugBundle:
                if wasSplit:
                    for group in splittedGroups:
                        streetGroup.plot('blue', 1, False)
                        splittingStreets.plot('whitesmoke', 1, False)
                        group.innerPlot('green', 1, False)
                        group.plot('red',1, True)

            if wasSplit:
                newGroups.extend(splittedGroups)
                delGroups.append(streetGroup)
                for street in splittingStreets:
                    if street.id in self.streets:
                        del self.streets[street.id]

        streetGroups.extend(newGroups)
        for group in delGroups:
            streetGroups.remove(group)

        for gIndex,streetGroup in enumerate(streetGroups):
            if not streetGroup:
                continue
            head, tail = orderHeadTail(streetGroup)

            if debugBundle:
                plotQualifiedNetwork(self.sectionNetwork)
                streetGroup.plot()
                for indx in range(len(head)):
                    item = head[indx]
                    p = item['firstVert']
                    # plt.text(p[0],p[1],'  '+str(item['i']),fontsize=12)
                    plt.plot(p[0],p[1],'coral',marker='o',markersize=14,zorder=998)
                    plt.text(p[0],p[1],'H'+str(indx),fontsize=10,zorder=999,horizontalalignment='center',verticalalignment='center')
                    plt.text(p[0]+2,p[1]-2,str(item['street'].id))

                for indx in range(len(tail)):
                    item = tail[indx]
                    p = item['firstVert']
                    # plt.text(p[0],p[1],'  '+str(item['i']),fontsize=12)
                    plt.plot(p[0],p[1],'skyblue',marker='o',markersize=14,zorder=998)
                    plt.text(p[0],p[1],'T'+str(indx),fontsize=10,zorder=999,horizontalalignment='center',verticalalignment='center')
                    plt.text(p[0]+2,p[1]-2,str(item['street'].id))

                plt.title('Heads (H) and tails (T) in streetGroup %d'%gIndex)
                plotEnd()

            innerStreets = findInnerStreets(streetGroup,self.leftHandTraffic)

            if debugBundle:

                for street in innerStreets:
                    if street not in streetGroup.group and street not in self.curvedTramLines:
                        allVertices = []
                        for item in street.iterItems():
                            if isinstance(item, Section):
                                item.polyline.plot('red',2,'solid',False,999)
                                allVertices.extend(item.centerline)
                        if len(allVertices):
                            c = sum(allVertices,Vector((0,0))) / len(allVertices)
                            plt.text(c[0],c[1],'S '+str(street.id),color='k',fontsize=8,zorder=130,ha='left', va='top', clip_on=True)

                if innerStreets:
                    plotQualifiedNetwork(self.sectionNetwork)
                    streetGroup.plot()
                    plt.title('Inner streets of street group %d'%gIndex)
                else:
                    plt.title('No inner streets in street group %d'%gIndex)
                plotEnd()

            bundle = Bundle()
            for item in head:
                street = item['street']
                street.bundle = bundle
                bundle.streetsHead.append(street)
                bundle.headLocs.append(item['firstVert'])
            for item in tail:
                street = item['street']
                street.bundle = bundle
                bundle.streetsTail.append(street)
                bundle.tailLocs.append(item['firstVert'])
            self.bundles[bundle.id] = bundle
            for street in innerStreets:
                if street not in streetGroup.group and street not in self.curvedTramLines:
                    street.bundle = bundle

            for street in streetGroup:
                street.bundle = bundle


    # If two bundles meet and there are no inner streets at the meeting intersection, they
    # can be merged into one bundle.
    def joinBundles(self):

        # An edge is a tuple  of the form (headKey, tailKey, bundle), where headKey is the set of locations
        # in self.headLocs and tailKey is the set of locations in self.tailLocs of a the bundle.
        # directConnected() checks, if two bundles, given as a sequence egde0, edge1 are connected only by
        # there key locations and connect only streets of the bundles.
        def directConnected(edge0,edge1):
            mergeKey0, bundle0 = edge0[0], edge0[2]
            bundleEndType0 = keyDict[(mergeKey0, bundle0)]
            mergeKey1, bundle1 = edge1[0], edge1[2]
            bundleStartType1 = keyDict[(mergeKey1, bundle1)]
            allBundleStreets = []
            allBundleStreets.extend(bundle0.streetsTail if bundleEndType0=='head' else bundle0.streetsHead)
            allBundleStreets.extend(bundle1.streetsHead if bundleStartType1=='head' else bundle1.streetsTail)

            # Collect the intersections at the connecting locations. All of them must be major intersections
            intersections = set(self.majorIntersections.get(loc,None) for loc in (bundle0.tailLocs if bundleEndType0=='head' else bundle0.headLocs) )
            if not all(intersections):
                return False

            # The streets of the intersections must be streets of the bundles.
            for intersection in  intersections:
                if intersection:
                    for conn in intersection:
                        if conn.item not in allBundleStreets:
                            return False                      
            return True

        # Create a list of all bundles as edges of a graph. An edge is a tuple  of the form
        # (headKey, tailKey, bundle), where headKey is the set of locations in self.headLocs
        # and tailKey is the set of locations in self.tailLocs of a the bundle. headKey and
        # tailKey are the nodes of the graph. The dictionary keyDict relates the combination
        # of the nodes and the bundle to the end type of the bundle.
        edgeList = []
        keyDict = dict()
        for id,bundle in self.bundles.items():
            headKey = tuple( sorted(set(bundle.headLocs)) )
            tailKey = tuple( sorted(set(bundle.tailLocs)) )
            keyDict[(headKey,bundle)] = 'head'
            keyDict[(tailKey,bundle)] = 'tail'
            edgeList.append( (headKey, tailKey, bundle) )

        # Decomposes the edges into maximally long, connected paths. Each edge may only appear
        # in a single path. Preference is given to starting at nodes with only one connection.
        # See description of find_all_lines in lib/CompGeom/chains.py.
        paths = find_all_lines(edgeList)

        # Only those edge pairs may be joined, that apply a direct connection (see description
        # of directConnected() above). The result is found in <cleanedPaths>.
        cleanedPaths = []    
        for path in paths:
            tempPath = []            
            for i, (edge0, edge1) in enumerate(pairs(path)):
                if not tempPath:
                    tempPath.append(edge0)  # Start a new path             
                if directConnected(edge0,edge1):
                    tempPath.append(edge1)
                else:
                    if len(tempPath) > 1:
                        cleanedPaths.append(tempPath)  # Save valid path
                    tempPath = []  # Reset path
            
            if tempPath and len(tempPath) > 1:
                cleanedPaths.append(tempPath)  # Save remaining path
 
        paths = cleanedPaths

        # If there are paths, their edges may be merged into a new bundle.
        for path in paths:
            newBundle = Bundle()

            # Fill this new bundle with the values of the first bundle in <path>.
            # Depending on the direction of ths edge, the bundle has to be reversed.
            firstBundle = path[0][2]
            bundleEndType = keyDict[(path[0][1], firstBundle)]
            if bundleEndType == 'tail':
                newBundle._pred = firstBundle._pred
                newBundle._succ = firstBundle._succ
                newBundle.streetsHead = firstBundle.streetsHead
                newBundle.streetsTail = firstBundle.streetsTail
                newBundle.headLocs = firstBundle.headLocs
                newBundle.tailLocs = firstBundle.tailLocs
            else:
                newBundle._pred = firstBundle._succ
                newBundle._succ = firstBundle._pred
                newBundle.streetsHead = firstBundle.streetsTail
                newBundle.streetsTail = firstBundle.streetsHead
                newBundle.headLocs = firstBundle.tailLocs
                newBundle.tailLocs = firstBundle.headLocs

            # Link all streets of <firstBundle> to the new bundle
            for id,street in self.streets.items():
                if street.bundle == firstBundle:
                    street.bundle = newBundle

            # Now, we start to join the remaining bundles to the new bundle.
            for edge in path[1:]:
                mergeKey = edge[0]
                bundle = edge[2]
                bundleEndType = keyDict[(mergeKey, bundle)]

                newBundle = joinBundles(self, newBundle, bundle, bundleEndType=='head')

            for edge in path:
                bundleID = edge[2].id
                if bundleID in self.bundles:
                    del self.bundles[bundleID]

            self.bundles[newBundle.id] = newBundle

        self.mergeBundles()

    # If two bundles meet and there are no inner streets at the meeting intersection, they
    # can be merged into one bundle.
    def mergeBundles(self):
        # Find all ends of streets of all bundles and cluster them to groups,
        # that are potential intersections.
        endPoints = []
        for id, bundle in self.bundles.items():
            for end, street in zip(bundle.headLocs,bundle.streetsHead):
                endPoints.append( (end,{'end':end, 'type':'head', 'street':street, 'bundle':bundle}) )
            for end, street in zip(bundle.tailLocs,bundle.streetsTail):
                endPoints.append( (end,{'end':end, 'type':'tail', 'street':street, 'bundle':bundle}) )
        isectCandidates = dbClusterScan(endPoints, dbScanDist, 2)

        toBeMerged = []
        for candidates in isectCandidates:
            involvedBundles = defaultdict(list)
            for _,cand in candidates:
                data = {'end':cand['end'], 'type':cand['type'], 'street':cand['street'], 'bundle':cand['bundle']}
                involvedBundles[cand['bundle']].append(data)

            if len(involvedBundles) == 2:
                if len(candidates) == 2:
                    continue
                if canBeMerged(self, involvedBundles):
                    toBeMerged.append(involvedBundles)

        for involvedBundles in toBeMerged:
            mergeBundles(self,involvedBundles)

    # Finally, the intersections of bundles are constructed.
    def createBundleIntersections(self):
        toBeIntersected = []

        endPoints = []
        for id, bundle in self.bundles.items():
            for end, street in zip(bundle.headLocs,bundle.streetsHead):
                endPoints.append( (end,{'end':end, 'type':'head', 'street':street, 'bundle':bundle}) )
            for end, street in zip(bundle.tailLocs,bundle.streetsTail):
                endPoints.append( (end,{'end':end, 'type':'tail', 'street':street, 'bundle':bundle}) )
        isectCandidates = dbClusterScan(endPoints, dbScanDist, 2)

        # <isectCandidates> is a list of potential intersection candidates.
        # These candidates are lists of dictionaries, that hold the location of the
        # streets end, the street's end type ('head' or 'tail'), the street instance
        # itself and the bundle, they belong to.
        for candidates in isectCandidates:
            # Check if there is a prohibited area. Skip, if any.
            cluster = [end.freeze() for end,_ in candidates]
            if any([len(set(cluster) & set(prohibited) ) for prohibited in self.prohibitedHulls]):
                continue

            involvedBundles = defaultdict(list)
            for _,cand in candidates:
                data = {'end':cand['end'], 'type':cand['type'], 'street':cand['street'], 'bundle':cand['bundle']}
                involvedBundles[cand['bundle']].append(data)

            nrOfBundles = len(involvedBundles)
            involvedBundleIDs = set( b.id for b,_ in involvedBundles.items())
            involvedBundleTypes = []
            for _,data in involvedBundles.items():
                involvedBundleTypes.extend( [d['type'] for d in data] )

            if nrOfBundles==1:
                # These are ends of bundles. Because these are not detected reliably
                # here, they are processed at the end of this method.
                pass

            if nrOfBundles==2:
                # These are bundles, that touch each other. If there is no street
                # to the inner side, they can be merged using pseudo minors.
                # Else, an intersection needs to be created.
                if len(candidates) == 2:
                    bundleIDs = [b.id for b,_ in involvedBundles.items()]
                    print('Single common bundle end between bundles', bundleIDs)
                    continue

                # If there are no common endpoints between these two bundles,
                # these are ends of bundles. They are processed at the end of this method.
                bundleInfo = list(involvedBundles.items())
                bundleEnds0 = {e for e in (bundleInfo[0][0].headLocs if bundleInfo[0][1][0]['type']=='head' else bundleInfo[0][0].tailLocs) }
                bundleEnds1 = {e for e in (bundleInfo[1][0].headLocs if bundleInfo[1][1][0]['type']=='head' else bundleInfo[1][0].tailLocs) }
                commonEnds = bundleEnds0.intersection(bundleEnds1)
                if not commonEnds:
                    continue


                toBeIntersected.append(involvedBundles)

            if nrOfBundles>2:
                toBeIntersected.append(involvedBundles)

        for involvedBundles in toBeIntersected:
            intersectBundles(self, involvedBundles)

        # Finally process the open Bundle ends.
        for id,bundle in self.bundles.items():
            if not bundle.pred or not bundle.succ:
                endBundleIntersection(self, bundle)




