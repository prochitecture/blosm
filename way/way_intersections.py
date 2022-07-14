from math import sin, cos, atan2, pi, sqrt, floor
from mathutils import Vector
from itertools import tee,islice, cycle
from itertools import cycle
from lib.CompGeom.PolyLine import PolyLine
from lib.CompGeom.algorithms import repairSimpleSelfIntersection
from way.way_properties import estFilletRadius
import matplotlib.pyplot as plt

# helper functions -----------------------------------------------
def cyclePair(lst):
    prevs, nexts = tee(lst)
    prevs = islice(cycle(prevs), len(lst) - 1, None)
    return zip(prevs,nexts)

def spline_4p( t, v0, p0, p1, v1 ):
    # (Catmull-Rom) Cubic curve goes from p0 to p1, and outer points v0 and v1
    # determine the slopes at p0 and p1.
    # assert 0 <= t <= 1
    return (
          t*((2-t)*t - 1)     * v0
        + (t*t*(3*t - 5) + 2) * p0
        + t*((4 - 3*t)*t + 1) * p1
        + (t-1)*t*t           * v1 ) / 2
# ----------------------------------------------------------------

class OutgoingWay():
    # This class holds a way-section that leaves an intersection. The direction
    # of its polyline, the sides 'left' and 'right', and the start and target
    # trim positions are redirected correctly to the original instance of WaySection.
    def __init__(self,section,fwd):
        # section:  Instance of a WaySection, that holds the way-section in its
        #           original direction.
        # fwd:      Direction of the outgoing section. True, if same as original.
        self.section = section
        self.fwd = fwd
        self.polyline = section.polyline.clone()
        self.polyline.setView( PolyLine.fwd if fwd else PolyLine.rev)
        self.isLoop = self.section.isLoop

    @property
    def leftW(self):
        return self.section.leftWidth if self.fwd else self.section.rightWidth

    @property
    def rightW(self):
        # negative, so that the offset in PolyLine's offset methods starts at the
        # beginning (where the intersection is) of the polyline
        return -self.section.rightWidth if self.fwd else -self.section.leftWidth

    def setTrim(self,trim):
        if self.fwd:
            self.section.trimS = trim
        else:
            self.section.trimT = len(self.polyline)-1 - trim

class Intersection():
    def __init__(self, position, network, waySections):
        self.position = position
        self.outWays = []
        self.order = 0

        for net_section in network.iterOutSegments(self.position):
            if net_section.category != 'scene_border':
                if waySections[net_section.sectionId].isValid:
                    self.addSection( position, waySections[net_section.sectionId] )

        self.sortSections()

    def addSection(self, position, waySection):
        # The polyline of the outgoing way of the intersection in <self.outWays> starts
        # at the intersection <position>. The instance of the class OutgoingWay controls the
        # relation to the original way-section given in <waySection>. <self.order> is the
        # number of outgoing ways and determines the type of the intersection. Loop ways have
        # their two ends at the intersection position and are marked in <waySection> with 
        # waySection.isLoop == True.
        if waySection.isLoop:
            self.outWays.append( OutgoingWay(waySection,True))
            self.outWays.append( OutgoingWay(waySection,False))
            self.order += 2
        else:
            if position == waySection.originalSection.s:
                self.outWays.append( OutgoingWay(waySection,True))
            else:
                self.outWays.append( OutgoingWay(waySection,False))
            self.order += 1

    def sortSections(self):
    # Sort the outgoing ways in <self.outWays> first by their angle in counter-clockwise order
    # around the cener of gravity of their start positions and then by the distance of their
    # end-point from this center.
        def sort_key(outway):
            vec = outway.polyline[-1] - self.position
            length = vec.length
            angle = atan2(vec.y,vec.x)
            if angle < 0:
                return 2*pi+angle, length
            else:
                return angle, length
        self.outWays = sorted(self.outWays, key = sort_key)
        # for i,outway in enumerate(self.outWays):
        #     p = outway.polyline[-1]
        #     plt.text(p[0],p[1],str(i),fontsize=12)

    def intersectionPoly_noFillet(self):
        tmax = [0.]*self.order  # List of maximum line parameters of the projections of the intersection
                                # points onto the center-lines of the ways.
        isectPoints = [None]*self.order

        doDebug = False

        for (indx1, indx2), (way1,way2) in zip(cyclePair(range(self.order)),cyclePair(self.outWays)):
            # <way1> and <way2> are consecutive outgoing polylines in counter-clockwise order of the
            # intersection. These are the centerlines of the way-sections, <way1> is to the right
            # (clockwise) of <way2>.

            # Compute the left border <borderL> of the way-section given by way1 and the 
            # right border <borderR> of the way-section given by way2.
            # For way-loops, the forward and backward polyline is in expected consecutive order
            # in <self.outWayLines>. For way-loops, a parallel offset can't be computed
            # correctly, therefore, only their first segments are used.
            if way1.section.isLoop and way2.section.isLoop:
                trimmed1 = PolyLine(way1.polyline[:2])
                trimmed2 = PolyLine(way2.polyline[:2])
                borderL = trimmed1.parallelOffset(way1.leftW)
                borderR = trimmed2.parallelOffset(way2.rightW)
            else:
                borderL = way1.polyline.parallelOffset(way1.leftW)
                borderR = way2.polyline.parallelOffset(way2.rightW)

            # Find the intersection between the left border <borderL> of <way1> 
            # and the right border <borderR> of <way2>.
            # In <isectParams>, <iP> is the intersection point, <tL> and <tR> are the
            # line parameters of the intersection point (see PolyLine.py) for the 
            # left and the right border. <uVL> and <uVR> are the unit vectors of the intersecting
            # segments.
            isectParams = PolyLine.intersection(borderL,borderR)
            if doDebug:
                borderL.plot('g')
                borderR.plot('r')
            if isectParams is None:
                if doDebug:
                    borderL.plot('c')
                    borderR.plot('m')
                # No intersection found
                continue
            iP,_,_,_,_ = isectParams
            if doDebug:
                plt.plot(iP[0],iP[1],'mx',markersize=10)

            # The line parameters of the intersection may be different from the line parameter of
            # the center-line. To get consistent and equal vertices for the way-polygon and the
            # intersection area, the projection of the intersection onto the center-line is used.
            _,tPL = way1.polyline.orthoProj(iP)
            _,tPR = way2.polyline.orthoProj(iP)

            # Find the maximal line parameter for every way
            tmax[indx1] = max(tmax[indx1],tPL)
            tmax[indx2] = max(tmax[indx2],tPR)
            isectPoints[indx1] = way1.polyline.offsetPointAt(tPL,way1.leftW)

        polygon = []
        # Construct the intersection area polygon
        for (indx1, indx2), (way1,way2) in zip(cyclePair(range(self.order)),cyclePair(self.outWays)):
            # <way1> and <way2> are consecutive outgoing polylines in counter-clockwise order of the
            # intersection. These are the centerlines of the way-sections, <way1> is to the right
            # (clockwise) of <way2>.

            isectP =  isectPoints[indx1]
            t1 = tmax[indx1]  # Line parameters of trim points.
            t2 = tmax[indx2]

            # The first way section will have to be trimmed to the length not occupied by
            # the polygon of the intersection area.
            way1.setTrim(t1)

            # Starting with the left border of the centerline <way1>.
            # Given the line parameter of the trim point on the centerline, we can 
            # compute the point <pL> perpendicularly offset to the left.
            pL = way1.polyline.offsetPointAt(t1,way1.leftW)
            if isectP:
                # If we had an intersection between the borders, by a simple check of
                # the coordinate difference to the trim point we may decide if
                # <pL> is the same as this intersection.
                dL = sum(pL-isectP)
                if abs(dL) > 1.e-4:
                    # If not on the intersection, add to polygon
                    polygon.append(pL)
                    # then extend polygon by the intersection.
                polygon.append(isectP)
            else:
                # Use only the trim point.
                polygon.append(pL)

            # The same procedure is then done for the right border of the centerline <way2>
            pR = way2.polyline.offsetPointAt(t2,way2.rightW)
            if isectP:
                dR = sum(pR-isectP)
                if abs(dR) > 1.e-4:
                    polygon.append(pR)
            else:
                polygon.append(pR)

        return polygon

    def intersectionPoly(self):
        tmax = [0.]*self.order  # List of maximum line parameters of the projections of the intersection
                                # points onto the center-lines of the ways.
        isectPoints = [None]*self.order
        filletVertsList = [None]*self.order

        doDebug = False

        for (indx1, indx2), (way1,way2) in zip(cyclePair(range(self.order)),cyclePair(self.outWays)):
            # <way1> and <way2> are consecutive outgoing polylines in counter-clockwise order of the
            # intersection. These are the centerlines of the way-sections, <way1> is to the right
            # (clockwise) of <way2>.

            # Compute the left border <borderL> of the way-section given by way1 and the 
            # right border <borderR> of the way-section given by way2.
            # For way-loops, the forward and backward polyline is in expected consecutive order
            # in <self.outWayLines>. For way-loops, a parallel offset can't be computed
            # correctly, therefore, only their first segments are used.
            if way1.section.isLoop and way2.section.isLoop:
                trimmed1 = PolyLine(way1.polyline[:2])
                trimmed2 = PolyLine(way2.polyline[:2])
                borderL = trimmed1.parallelOffset(way1.leftW)
                borderR = trimmed2.parallelOffset(way2.rightW)
            else:
                borderL = way1.polyline.parallelOffset(way1.leftW)
                borderR = way2.polyline.parallelOffset(way2.rightW)

            # Find the intersection between the left border <borderL> of <way1> 
            # and the right border <borderR> of <way2>.
            # In <isectParams>, <iP> is the intersection point, <tL> and <tR> are the
            # line parameters of the intersection point (see PolyLine.py) for the 
            # left and the right border. <uVL> and <uVR> are the unit vectors of the intersecting
            # segments.
            isectParams = PolyLine.intersection(borderL,borderR)
            if doDebug:
                borderL.plot('g')
                borderR.plot('r')
            if isectParams is None:
                if doDebug:
                    borderL.plot('c')
                    borderR.plot('m')
                # No intersection found
                continue
            iP,tL,tR,uVL,uVR = isectParams
            if doDebug:
                plt.plot(iP[0],iP[1],'mx',markersize=10)

            # Now, the fillet of this intersection gets constructed.
            # Get the fillet radius between these ways based on the ways' categories.
            # The smaller of them is used.
            radius1 = estFilletRadius(way1.section.originalSection.category,way1.section.originalSection.tags)
            radius2 = estFilletRadius(way2.section.originalSection.category,way2.section.originalSection.tags)
            radius = min(radius1,radius2)

            # For small angles, reduce the fillet radius continuously
            cosAngle = uVL.dot(uVR)
            if cosAngle > 0.3:  # ~28°
                radius *= cosAngle

            # <lengthL> and <lengthR> are the lengths on the intersected segments between
            # their start point and the intersection
            lengthL = borderL.lengthOfSegment(floor(tL))
            lengthR = borderR.lengthOfSegment(floor(tR))

            # Starting with the radius found above, a fillet is created within the intersected
            # segments. Its legs must not be longer than the way-segments. If required,
            # the raidus is decreased iteratively in the follwoing loop.
            while True:
                # Try to find an arc with radius <radius> and using the segments starting at <iP>
                # and suing the unit vectors <uVL> and <uVR> as tangents. In the result, <origin>
                # is the circle origin of the arc, <pLegEndL> and <pLegEndR> are the tangent vertices
                # and <legLength> is the length of the legs between the intersection point and the
                # tangent vertices.
                origin, pLegEndL, pLegEndR, legLength = filletArc(iP, uVL, uVR, radius)
                if origin is None: # No arc possible, when segments are almost parallel.
                    break

                # To check whether the leg ends are on the intersecting segments, the line
                # parameters of the tangent vertices are computed (see PolyLine.py). For a 
                # valid fillet, they must be both less than 1.
                tLegL = tL%1. + legLength / lengthL
                tLegR = tR%1. + legLength / lengthR
                if tLegL < 1. and tLegR < 1.:
                    # OK, fillet fits between these segments using this radius
                    break
                # if not, reduce radius
                radius = 0.9 * radius

            # If a fillet arc with legs that fit onto the segments has been found.
            if origin:
                if doDebug:
                    plt.plot(pLegEndL[0],pLegEndL[1],'go',markersize=6)
                    plt.plot(pLegEndR[0],pLegEndR[1],'ro',markersize=6)
                # Create the vertices for the fillet.
                filletVerts = filletLine(origin, pLegEndL, pLegEndR, radius)

                # The line parameters of the intersection may be different from the line parameter of
                # the center-line. To get consistent and equal vertices for the way-polygon and the
                # intersection area, the projection of the intersection onto the center-line is used.
                _,tPL = way1.polyline.orthoProj(pLegEndL)
                _,tPR = way2.polyline.orthoProj(pLegEndR)

                # Find the maximal line parameter for every way
                tmax[indx1] = max(tmax[indx1],tPL)
                tmax[indx2] = max(tmax[indx2],tPR)
                isectPoints[indx1] = way1.polyline.offsetPointAt(tPL,way1.leftW)

                # Store the fillet vertices for <way1>.
                filletVertsList[indx1] = filletVerts

            else:
                # No possible fillet found
                print('No possible fillet found.')

        polygon = []
        # Construct the intersection area polygon
        for (indx1, indx2), (way1,way2) in zip(cyclePair(range(self.order)),cyclePair(self.outWays)):
            # <way1> and <way2> are consecutive outgoing polylines in counter-clockwise order of the
            # intersection. These are the centerlines of the way-sections, <way1> is to the right
            # (clockwise) of <way2>.

            filletVerts = filletVertsList[indx1]
            isectP =  isectPoints[indx1]
            t1 = tmax[indx1]  # Line parameters of trim points.
            t2 = tmax[indx2]

            # The first way section will have to be trimmed to the length not occupied by
            # the polygon of the intersection area.
            way1.setTrim(t1)

            # Starting with the left border of the centerline <way1>.
            # Given the line parameter of the trim point on the centerline, we can 
            # compute the point <pL> perpendicularly offset to the left.
            pL = way1.polyline.offsetPointAt(t1,way1.leftW)
            if filletVerts:
                # If we had an fillet between the borders, by a simple check of
                # the coordinate difference of the first vertex of the fillet 
                # to the trim point we may decide if <pL> is the same as this vertex.
                dL = sum(pL-filletVerts[0])
                if abs(dL) > 1.e-4:
                    # If not on the fillet, add to polygon
                    polygon.append(pL)
                # then extend polygon by fillet verts.
                polygon.extend(filletVerts)
            else:
                # Use only the trim point.
                polygon.append(pL)


            # The same procedure is then done for the right border of the centerline <way2>
            pR = way2.polyline.offsetPointAt(t2,way2.rightW)
            if filletVerts:
                dR = sum(pR-filletVerts[-1])
                if abs(dR) > 1.e-4:
                    polygon.append(pR)
            else:
                polygon.append(pR)

        return polygon

    def findTransitionPoly(self):
        # Transitions are intersections of order==2 and have an incoming and
        # an outgoing way. They have to be processed before the other intersections,
        # because way widths may be altered due to turning lanes.
        way1, way2 = self.outWays
        outWay, inWay = (way1, way2) if way1.fwd else (way2, way1)
        inTags, outTags = inWay.section.originalSection.tags, outWay.section.originalSection.tags

        # Do we have turn lanes? They are only possible in the outLine.
        if 'turn:lanes' in outTags:
            # There is no transition polygon required. The outgoing way section
            # becomes eventually a turning lane.
            laneDescs = outTags['turn:lanes'].split('|')
            leftTurnLanes = sum(1 for tag in laneDescs if 'left' in tag)
            rightTurnLanes = sum(1 for tag in laneDescs if 'right' in tag)
            if leftTurnLanes or rightTurnLanes:
                leftWidthDifference = outWay.section.leftWidth - inWay.section.leftWidth
                rightWidthDifference = outWay.section.rightWidth - inWay.section.rightWidth
                if leftWidthDifference or rightWidthDifference:
                    outWay.section.turnParams = [leftWidthDifference,rightWidthDifference]

        # Prepare transition polygon 
        if outWay.section.turnParams:
            transitionWidth = 0.5
        else:
            inWidth = inWay.section.leftWidth + inWay.section.leftWidth
            outWidth = outWay.section.leftWidth + outWay.section.leftWidth
            widthDiff = abs(inWidth-outWidth)
            transitionWidth = min(10.,max(1.,2 * widthDiff))
        inT = inWay.polyline.d2t(transitionWidth)
        outT = outWay.polyline.d2t(transitionWidth)
        if inT is None or outT is None:
            pass
            # plt.plot(self.position[0],self.position[1],'rx',markersize=20,zorder=950)
        else:
            assert inT is not None
            assert outT is not None

            # compute cornerpoints
            # inLeft = inLine.polyline.offsetPointAt(inT,inLine.section.leftWidth)
            # inRight = inLine.polyline.offsetPointAt(inT,-inLine.section.rightWidth)
            inLeft = inWay.polyline.offsetPointAt(inT,inWay.section.rightWidth)
            inRight = inWay.polyline.offsetPointAt(inT,-inWay.section.leftWidth)
            outLeftWidth = inWay.section.leftWidth if outWay.section.turnParams else outWay.section.leftWidth
            outRightWidth = inWay.section.rightWidth if outWay.section.turnParams else outWay.section.rightWidth
            outLeft = outWay.polyline.offsetPointAt(outT,outLeftWidth)
            outRight = outWay.polyline.offsetPointAt(outT,-outRightWidth)

            poly = []
            import numpy as np
            t0 = np.linspace(0,1,10)
            p0, p1 = inRight, outLeft
            v0, v1 = p0 - inWay.polyline.unitEndVec(False)*transitionWidth, p1 - outWay.polyline.unitEndVec(False)*transitionWidth
            for t in t0:
                sp = spline_4p( t, v0, p0, p1, v1 )
                poly.append( sp )
            p0, p1 = outRight, inLeft
            v0, v1 = p0 - outWay.polyline.unitEndVec(False)*transitionWidth, p1 - inWay.polyline.unitEndVec(False)*transitionWidth
            for t in t0:
                sp = spline_4p( t, v0, p0, p1, v1 )
                poly.append( sp )


            if inWay.fwd:
                inWay.section.trimS = inT
            else:
                inWay.section.trimT = len(inWay.polyline)-1 - inT
            if outWay.fwd:
                outWay.section.trimS = outT
            else:
                outWay.section.trimT = len(outWay.polyline)-1 - outT

            return poly


def filletArc(p, uv1, uv2, radius):
    # p: corner point (class Vector)
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
    o = tp2 + Vector((uv2.y,-uv2.x)) * radius
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

    acurr = 0.
    # pt = Coordinate()
    while acurr < atot:
        a = a1 + acurr
        pt = o + Vector((radius * cos(a),radius * sin(a)))
        vertList.append(pt)
        acurr += ainc

    vertList.append(tp1)
    vertList.reverse()
    return vertList
