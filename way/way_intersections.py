
from lib.pygeos.geom import GeometryFactory
from lib.pygeos.shared import Coordinate, CAP_STYLE, JOIN_STYLE
from lib.SweepIntersectorLib.SweepIntersector import SweepIntersector
from way.PolyLine import PolyLine

from math import sin, cos, atan2, pi, sqrt, fmod
from mathutils import Vector
import numpy as np

import matplotlib.pyplot as plt
from action.plotUtilities import *

roadTypes = {
    "motorway" : 'driving',
    "motorway_link" : 'driving',
    "trunk" : 'driving',
    "trunk_link" : 'driving',
    "primary" : 'driving',
    "primary_link" : 'driving',
    "secondary" : 'driving',
    "secondary_link" : 'driving',
    "tertiary" : 'driving',
    "tertiary_link" : 'driving',
    "unclassified" : 'default',
    "residential" : 'driving',
    "living_street" : 'driving',
    "service" : 'driving',
    "pedestrian" : 'walking',
    "track" : 'walking',
    "escape" : 'driving',
    "raceway" : 'driving',
    "steps" : 'walking',
    "footway" : 'walking',
    "path" : 'walking',
    "cycleway" : 'cycling',
    "bridleway" : 'walking'
}

def getRoadWidth(tags):
    rType = tags.get('highway',None)
    if rType:
        rType = roadTypes[rType]
        if rType == 'driving':
            lanes = tags.get('lanes',None)
            if lanes:
                nrOfLanes = int(lanes)
            else:
                nrOfLanes = 1
            if not tags.get('oneway','no') == 'yes':
                nrOfLanes *= 2
            return nrOfLanes * 2.7#2.5
        elif rType == 'cycling':
            return 1.8
        elif rType == 'walking':
            return 3.0
        else:
            return 1.0  # ????
    else:
        return 1.0

    pass

class Section():
    def __init__(self,section):
        self.originalSection = section
        self.halfWidth = getRoadWidth(section.tags)/2.
        geosF = GeometryFactory()
        self.polyline = PolyLine(section.path)

# sort key for angles
def sort_key(section):
    vec = section.firstV
    length = vec.length
    angle = atan2(vec[1],vec[0])
    if angle < 0:
        return 2*pi+angle, length
    else:
        return angle, length

class Intersection():
    def __init__(self,position,sections):
        self.position = position
        self.sections = sections
        self.outPolyLines = []

    def addSection(self,section):
        polyline = { 'halfWidth':section.halfWidth,'line':section.polyline}
        if section.originalSection.s != self.position:
            polyline = { 'halfWidth':section.halfWidth,'line':section.polyline.reversed()}
        self.outPolyLines.append(polyline)

    def mergeNode(self,node,sections,network):
        nodes_to_merge = []
        for net_section in network.iterOutSegments(node):
            if net_section.category != 'scene_border':
                if net_section.length < 10.:
                    pass
                # #     # nodes_to_merge.append(net_section.t)
                # #     # self.shortSectionIDs.append(net_section.sectionId)
                # #     # if net_section.sectionId in sectionList:
                # #     #     del sectionList[net_section.sectionId]
                else:
                    self.addSection(sections[net_section.sectionId])
        return nodes_to_merge

    def sortSections(self):
    # sort key for angles
        def sort_key(line):
            vec = line['line'].coords[1]-line['line'].coords[0]
            length = vec.length
            angle = atan2(vec.y,vec.x)
            if angle < 0:
                return 2*pi+angle, length
            else:
                return angle, length
        self.outPolyLines = sorted(self.outPolyLines, key = sort_key)

    def findIntersectionPoly(self):
        N = len(self.outPolyLines)
        percLeft = [0.]*N
        percRight = [0.]*N
        allCircleVerts = [None]*N
        leftLines = [None]*N
        rightLines = [None]*N
        for l1, l2, line1,line2 in zip(range(N),range(1,N+1),self.outPolyLines,self.outPolyLines[1:]+[self.outPolyLines[0]]):
            try:
                l2 = l2%N
                leftLine = line1['line'].createOffsetPolyLine(line1['halfWidth'])
                rightLine = line2['line'].createOffsetPolyLine(-line2['halfWidth'])
                leftLines[l1] = leftLine
                rightLines[l2] = rightLine
                # plotPolyLine(leftLine,0.,'#00ff00')
                # plotPolyLine(rightLine,0.,'#ff0000')

                # find intersection between left border of line 1 and right border of line 2
                # segLength is the sum of the segment number and the percentage of this segment
                isectP, segLength1, segLength2, unitV1, unitV2 = PolyLine.intersection(leftLine,rightLine)
                if not isectP:
                    continue
                # plt.plot(isectP.x,isectP.y,'mo',markersize=6)

                # get first estimate of fillet between these segments
                maxVel = 50.
                radius = maxVel/9.81/0.7/2. #???
                while True:
                    origin, tangentP1, tangentP2, legLength = filletCircle(isectP, unitV1, unitV2, radius)
                    if origin is None:
                        break
                    # check whether the leg ends we are on the same segment as the intersection point
                    percentage1 = fmod(segLength1,1) + legLength /leftLine.segmentLength(int(segLength1))
                    percentage2 = fmod(segLength2,1) + legLength / rightLine.segmentLength(int(segLength2))
                    if percentage1 < 1. and percentage2 < 1.:
                        break
                    # if not, reduce radius
                    ax = plt.gca()
                    ax.add_artist(plt.Circle(
                        (isectP.x,isectP.y),
                        10,
                        alpha=0.1,
                        color='g',
                        zorder=100
                    )) 

                    radius = 0.9 * radius

                if origin:
                    circleVerts = filletLine(origin, tangentP1, tangentP2, radius)
                    percLeft[l1] = percentage1
                    percRight[l2] = percentage2
                    allCircleVerts[l1] = circleVerts
                    for v1,v2 in zip(circleVerts[:-1],circleVerts[1:]):
                        plt.plot([v1.x,v2.x],[v1.y,v2.y],'r')
                    # if circleVerts:
                    #     p = circleVerts[0]
                    #     plt.text(p.x,p.y,' %f5.2'%(percentage1))
                    #     p = circleVerts[-1]
                    #     plt.text(p.x,p.y,' %f5.2'%(percentage2))
                        # p,t = line2['line'].trimLength(circleVerts[-1],Coordinate(unitV2.y,-unitV2.x))
                        # if t > trimLength[l2%N]:
                        #     maxLine[l2%N] = 'right'
                        # trimLength[l2%N] = max(trimLength[l2%N], t)
                        # p,t = line1['line'].trimLength(circleVerts[0],Coordinate(unitV1.y,-unitV1.x))
                        # if t > trimLength[l1]:
                        #     maxLine[l1] = 'left'
                        # trimLength[l1] = max(trimLength[l1], t)
                else:
                    continue

                    # plt.plot(isectP.x,isectP.y,'mo',markersize=6,zorder=700)
            except Exception as e:
                # plotEnd()
                print('Exception 1 ')# + str(nodeNr))

        # for i,t in enumerate(trimLength):
        #     v1,v2 = self.outPolyLines[i]['line'].segment(int(t))
        #     d = v2-v1
        #     p = v1 + (v2-v1) * fmod(t,1)
        #     plt.plot(p.x,p.y,'mo',markersize=6,zorder=700)
        #     plt.text(p.x,p.y,str(i))

        try:
            polygon = []
            for i in range(N):
                if percLeft[i] > percRight[i]:
                    p = rightLines[i].trimPoint(percLeft[i])
                    # plt.text(p.x,p.y,' %f5.2'%(percLeft[i]))
                    polygon.append(p)
                else:
                    p = leftLines[i].trimPoint(percRight[i])
                    # plt.text(p.x,p.y,' %f5.2'%(percRight[i]))
                    polygon.append(p)
                # plt.text(p.x,p.y,str(i)+' x')
                if allCircleVerts[i]:
                    p = allCircleVerts[i][0]
                    # plt.text(p.x,p.y,str(i))
                if allCircleVerts[i] is not None:
                    polygon.extend(allCircleVerts[i])
                test=1
            #     if maxLine[i] == 'left':
            #         polygon.extend(allCircleVerts[i])
            #         p = rightLine.trimPoint(trimLength[i])
            #         plt.plot(p.x,p.y,'co',markersize=6,zorder=700)
            #         plt.text(p.x,p.y,str(k))
            #         k += 1
            #         polygon.append(p)
            #     else:
            #         p = leftLine.trimPoint(trimLength[i])
            #         plt.text(p.x,p.y,str(k))
            #         k += 1
            #         plt.plot(p.x,p.y,'co',markersize=6,zorder=700)
            #         polygon.append(p)
            #         polygon.extend(allCircleVerts[i])
            geosF = GeometryFactory()
            poly = geosF.createPolygon( geosF.createLinearRing(polygon+[polygon[0]]) )
            plotGeosPolyFill(poly,'#cf4b23',1,1.0,700)
        except Exception as e:
            # plotEnd()
            print('Exception 2 ')# + str(nodeNr))
            

    # def findCollisions(self):
    #     N = len(self.outPolyLines)
    #     trimLength = [0.]*N
    #     for l1, l2, line1,line2 in zip(range(N),range(1,N+1),self.outPolyLines,self.outPolyLines[1:]+[self.outPolyLines[0]]):
    #         try:
    #             leftLine = line1['line'].createOffsetPolyLine(line1['halfWidth'])
    #             rightLine = line2['line'].createOffsetPolyLine(-line2['halfWidth'])
    #             # plotPolyLine(leftLine,0.,'#00ff00')
    #             # plotPolyLine(rightLine,0.,'#ff0000')
    #             c, i1, i2, u1, u2 = PolyLine.intersection(leftLine,rightLine)
    #             if c:
    #                 # plt.plot(c.x,c.y,'kx',markersize=10)
    #                 r = 2.#50./9.81/0.7/2
    #                 segs = filletCircle(c, u1, u2, r)
    #                 if segs:
    #                     for v1,v2 in zip(segs[:-1],segs[1:]):
    #                         plt.plot([v1.x,v2.x],[v1.y,v2.y],'r')

    #                     p,t = line1['line'].trimLength(segs[-1],Coordinate(u1.y,-u1.x))
    #                     trimLength[l1] = max(trimLength[l1], t)
    #                     # if p:
    #                     #     plt.plot(p.x,p.y,'ko',markersize=6,zorder=700)
    #                     p,t = line2['line'].trimLength(segs[0],Coordinate(u2.y,-u2.x))
    #                     trimLength[l2] = max(trimLength[l2], t)
    #                     # if p:
    #                     #     plt.plot(p.x,p.y,'mo',markersize=6,zorder=700)
    #                         # test = 1

                
    #         except:
    #             print('exception')
    #     for i,t in enumerate(trimLength):
    #         v1,v2 = self.outPolyLines[i]['line'].segment(int(t))
    #         d = v2-v1
    #         p = v1 + (v2-v1) * fmod(t,1)
    #         plt.plot(p.x,p.y,'mo',markersize=6,zorder=700)

def intersectPolyLines(listOfLines):
    segs = []
    for line in listOfLines:
        segs.extend([((v1.x,v1.y),(v2.x,v2.y)) for v1,v2 in zip(line.coords[:-1],line.coords[1:])])
    intersector = SweepIntersector()
    intersector.findIntersections(segs)
    return set(item for sublist in intersector.isectDict.values() for item in sublist)

def filletCircle(p, uv1, uv2, radius):
    # p: corner point (class Coordinate)
    # uv1, uv2: unit vectors from p to legs direction
    # cos(a) for angle a between uv1 and uv2
    cos_a = uv1.x*uv2.x + uv1.y*uv2.y # dot product
    if abs(cos_a) >= 1.-1.e-4:
        return None, None, None, None
    # tan(a/2) = sqrt((1 - cos(a)) / (1 + cos(a))
    tan_a2 = sqrt( (1 - cos_a) / (1 + cos_a) )
    # length of legs to tangent points
    length = radius / tan_a2
    # tangent points
    tp1 = p + uv1 * length
    tp2 = p + uv2 * length
    # origin of circle
    o = tp2 + Coordinate(uv2.y,-uv2.x) * radius
    # plt.plot(tp1.x,tp1.y,'b.')
    # plt.plot(tp2.x,tp2.y,'r.')
    # plt.plot(o.x,o.y,'g.')

    return  o, tp1, tp2, length

def filletLine(o, tp1, tp2, radius):
    vertList = [tp2]
    a2 = atan2(tp1.y-o.y, tp1.x-o.x)
    a1 = atan2(tp2.y-o.y, tp2.x-o.x)

    if a1 > a2:
        a2 += 2*pi
    atot = a2 - a1

    QUADRANT_SEGMENTS = 8
    filletAngleQuantum = pi / 2.0 / QUADRANT_SEGMENTS
    nSegs = int(atot / filletAngleQuantum + 0.5)

    # no segments because angle is less than increment-nothing to do!
    if nSegs < 1:
        return []

    # choose angle increment so that each segment has equal length
    ainc = atot / nSegs

    acurr = ainc
    pt = Coordinate()
    while acurr < atot:
        a = a1 + acurr
        pt = o + Coordinate(radius * cos(a),radius * sin(a))
        vertList.append(pt)
        acurr += ainc

    # pt = o - Vector( (radius * cos(startAngle-currAngleInc),radius * sin(startAngle-currAngleInc) ) )
    vertList.append(tp1)

    # segList.append(p2)
    vertList.reverse()
    return vertList


# def filletCircle(p, uv1, uv2, radius):
#     # p: corner point (class Coordinate)
#     # uv1, uv2: unit vectors from p to legs direction
#     # cos(a) for angle a between uv1 and uv2
#     cos_a = uv1.x*uv2.x + uv1.y*uv2.y # dot product
#     if abs(cos_a) >= 1.:
#         return []
#     # tan(a/2) = sqrt((1 - cos(a)) / (1 + cos(a))
#     tan_a2 = sqrt( (1 - cos_a) / (1 + cos_a) )
#     # length of legs to tangent points
#     length = radius / tan_a2
#     p1 = p + uv1 * length
#     p2 = p + uv2 * length
#     # sign of cross product
#     # sign = 1. if uv1[0]*uv2[1] + uv1[1]*uv2[0]>0. else -1.
#     o = p2 + Coordinate(uv2.y,-uv2.x) * radius
#     plt.plot(p1.x,p1.y,'b.')
#     plt.plot(p2.x,p2.y,'r.')
#     plt.plot(o.x,o.y,'g.')

#     segList = [p2]
#     a2 = atan2(p1.y-o.y, p1.x-o.x)
#     a1 = atan2(p2.y-o.y, p2.x-o.x)

#     if a1 > a2:
#         a2 += 2*pi
#     atot = a2 - a1

#     QUADRANT_SEGMENTS = 8
#     filletAngleQuantum = pi / 2.0 / QUADRANT_SEGMENTS
#     nSegs = int(atot / filletAngleQuantum + 0.5)

#     # no segments because angle is less than increment-nothing to do!
#     if nSegs < 1:
#         return []

#     # choose angle increment so that each segment has equal length
#     ainc = atot / nSegs

#     acurr = ainc
#     pt = Coordinate()
#     while acurr < atot:
#         a = a1 + acurr
#         pt = o + Coordinate(radius * cos(a),radius * sin(a))
#         segList.append(pt)
#         acurr += ainc

#     # pt = o - Vector( (radius * cos(startAngle-currAngleInc),radius * sin(startAngle-currAngleInc) ) )
#     segList.append(p1)

#     # segList.append(p2)
#     return segList



# class Intersection():
#     def __init__(self,position):
#         self.position = position
#         self.sections = []

#     def addSections(self,sectionList):
#         for section in sectionList:
#             self.sections.append(Section(section))

#         # sort them by angle
#         self.sections = sorted(self.sections, key = sort_key)

#     def computeCollisions(self):
#         for section1,section2 in zip(self.originalSections,self.originalSections[1:] + [self.originalSections[0]]):
#             # from https://stackoverflow.com/a/565282
#             p = section1.leftStart
#             r = section1.firstV
#             q = section2.rightStart
#             s = section2.firstV

#             r_cross_s = r.cross(s)
#             if r_cross_s != 0.:
#                 t = (q-p).cross(s / r_cross_s)
#                 c = p + t * r
#                 section1.collisions.append(c)
#                 section2.collisions.append(c)
#         test=1

#     def boundaryCollisions(self):
#         geosF = GeometryFactory()
#         for section1,section2 in zip(self.sections,self.sections[1:] + [self.sections[0]]):
#             poly1 = PolyLine(section1.leftPolyline.coords)
#             poly2 = PolyLine(section2.rightPolyline.coords).reversed()
#             c,i1,i2,u1,u2 = PolyLine.intersection(poly1,poly2)
#             if c:
#                 # v1 = poly1[i1]
#                 # v2 = poly1[i1+1]
#                 # plt.plot([v1.x,v2.x],[v1.y,v2.y],'k:')
#                 # # plt.plot(v2.x,v2.y,'bo')
#                 # v1 = poly2[i2]
#                 # v2 = poly2[i2+1]
#                 # plt.plot([v1.x,v2.x],[v1.y,v2.y],'k:')
#                 # plt.plot(v2.x,v2.y,'ro')
#                 collisions = intersectPolyLines([section1.leftPolyline,section2.rightPolyline])
#                 segs = filletCircle(c, u1, u2, 1)
#                 if segs:
#                     for v1,v2 in zip(segs[:-1],segs[1:]):
#                         plt.plot([v1.x,v2.x],[v1.y,v2.y],'r')
#                 section1.collisions.append(c)
#                 section2.collisions.append(c)

#         # for section1,section2 in zip(reversed(self.sections),reversed(self.sections[1:] + [self.sections[0]])) :
#         #     collisions = intersectPolyLines([section1.rightPolyline,section2.leftPolyline,section2.rightPolyline])
#         #     for c in collisions:
#         #         section1.collisions.append(c)
#         #         section2.collisions.append(c)

#         print([len(section.collisions) for section in self.sections])


#             # geosCoords1 = [geosF.createCoordinate(v) for v in section1.originalSection.path]
#             # geosString = geosF.createLineString(geosCoords1)
#             # buffer1 = geosString.buffer(section1.halfWidth,resolution=3,cap_style=CAP_STYLE.flat)
#             # boundary1 = buffer1.boundary

#             # geosCoords2 = [geosF.createCoordinate(v) for v in section2.originalSection.path]
#             # geosString = geosF.createLineString(geosCoords2)
#             # buffer2 = geosString.buffer(section2.halfWidth,resolution=3,cap_style=CAP_STYLE.square)
#             # boundary2 = buffer2.boundary

#             # collisionPoly = buffer1.intersection(buffer2)
#             # # coll = np.array(boundary1.intersection(boundary2))
#             # for cv in collisionPoly.exterior.coords[:-1]:
#             #     if cv not in geosCoords1 and cv not in geosCoords2:
#             #         c = Vector((cv.x,cv.y))
#             #         section1.collisions.append(c)
#             #         section2.collisions.append(c)

#             # test=1

#     def unionIntersections(self):
#         geosF = GeometryFactory()
#         lines = []
#         for section in self.sections:
#             lines.append(section.centerPolyline.buffer(0.001,resolution=3,cap_style=CAP_STYLE.square))

#         line = geosF.createMultiPolygon(lines)
#         uni = line.union()
#         f = uni.buffer(3,resolution=3,cap_style=CAP_STYLE.square)
#         plotUtilities.plotGeosWithHoles(f,False)
#         plotUtilities.plotEnd()
#         test=1


