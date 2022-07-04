from math import sin, cos, atan2, pi, sqrt, floor, ceil
from mathutils import Vector
from itertools import tee,islice, cycle
from itertools import cycle
from lib.CompGeom.PolyLine import PolyLine
from way.way_properties import estFilletRadius
from way.way_properties import estimateWayWidth, getLanes
import matplotlib.pyplot as plt

# helper functions and classes -----------------------------------
def cyclePair(lst):
    prevs, nexts = tee(lst)
    prevs = islice(cycle(prevs), len(lst) - 1, None)
    return zip(prevs,nexts)

class OutLine():
    def __init__(self,section,polyline,fwd):
        self.section = section
        self.polyline = polyline
        self.fwd = fwd

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

class Intersection():
    def __init__(self,position,network,waySections,shortestSection=0.):
        self.position = position
        self.outWayLines = []
        self.order = 0

        for net_section in network.iterOutSegments(self.position):
            if net_section.category != 'scene_border':
                if net_section.length < shortestSection:
                    pass
                else:
                    self.addSection(waySections[net_section.sectionId])

        self.sortSections()

    def addSection(self,waySection):
        # The polylines of the outgoing lines of the intersection in <self.outWayLines> start all
        # at the intersection point. If this is the same direction as in the corresponding way
        # section <waySection>, the attribute <fwd> in <Outline> is set to True, and False,
        # if reversed. <self.order> is the number of outgoing ways and determines the type of the
        # intersection. Loop ways have their two ends at the intersection and are marked in 
        # <waySection> with waySection.isLoop == True.
        if waySection.isLoop:
            self.outWayLines.append( OutLine(waySection,waySection.polyline.clone(),True))
            self.outWayLines.append( OutLine(waySection,waySection.polyline.reversed(),False))
            self.order += 2
        else:
            if self.position == waySection.originalSection.s:
                self.outWayLines.append( OutLine(waySection,waySection.polyline.clone(),True))
            else:
                self.outWayLines.append( OutLine(waySection,waySection.polyline.reversed(),False))
            self.order += 1

    def sortSections(self):
    # Sort the outlines of the intersection first by their angle in counter-clockwise order
    # and secondary by their distance from their start point.
        def sort_key(outline):
            vec = outline.polyline[1]-outline.polyline[0]
            length = vec.length
            angle = atan2(vec.y,vec.x)
            if angle < 0:
                return 2*pi+angle, length
            else:
                return angle, length
        self.outWayLines = sorted(self.outWayLines, key = sort_key)

        # indx = 0
        # plt.close()
        # for line in self.outWayLines:
        #     line.polyline.plot('k')
        #     e = line.polyline.verts[1]
        #     plt.text(e[0],e[1],str(indx))
        #     # if indx==1:
        #     #     borderL = line.polyline.parallelOffset( line['halfWidth'])
        #     #     borderL.plot('#00ff00')
        #     indx += 1
        # plotEnd()

    def findIntersectionPoly(self):
        # The outgoing polylines in <self.outWayLines> of the intersection 
        # are of class Outline. They are expected here to be sorted in
        # counter-clockwise order.
        N = self.order              # Number of outgoing polylines
        tSmax = [0.]*N              # List of maximum line parameter along polyline (meaning see PolyLine.py)
        filletVertsList = [None]*N      # List of fillet vertices per border line intersection
        bordersL = [None]*N         # List of left border polylines of outgoing ways
        bordersR = [None]*N         # list of right border polylines of outgoing ways

        doDebug = False

        for (lineIndx1, lineIndx2), (line1,line2) in zip(cyclePair(range(N)),cyclePair(self.outWayLines)):
            # <line1> and <line2> are consecutive outgpoing polylines in counter-clockwise order of the
            # intersection. These are the centerlines of the way-sections, <line1> is to the right of <line2>.
            try:
                # Compute the left border <borderL> of the way-section given by line1 and the 
                # right border <borderR> of the way-section given by line2.
                # For way-loops, the forward and backward polyline is in expected consecutive order
                # in <self.outWayLines>. For loop-ways, a parallel offset can't be computed
                # correctly, therefore, only their first segments are used.
                w1 = line1.section.leftWidth if line1.fwd else line1.section.rightWidth
                w2 = line2.section.rightWidth if line2.fwd else line2.section.leftWidth
                if line1.section.isLoop and line2.section.isLoop:
                    trimmed1 = PolyLine(line1.polyline[:2])
                    trimmed2 = PolyLine(line2.polyline[:2])
                    borderL = trimmed1.parallelOffset(w1 )
                    borderR = trimmed2.parallelOffset(-w2 )
                else:
                    borderL = line1.polyline.parallelOffset(w1 )
                    borderR = line2.polyline.parallelOffset(-w2 )

                # Find intersection between left border <borderL> of <line1> 
                # and right border <borderR> of <line2>
                # In <isectParams>, <iP> is the intersection point, <tL> and <tR> are the
                # line parameters of the intersection point (see PolyLine.py) for left and
                # right border. <uVL> and <uVR> are the unit vectors of the intersecting
                # segments
                isectParams = PolyLine.intersection(borderL,borderR)
                if isectParams is None:
                    if doDebug:
                        borderL.plot('c')
                        borderR.plot('m')
                    continue
                iP,tL,tR,uVL,uVR = isectParams

                if doDebug:
                    plt.plot(iP[0],iP[1],'mx',markersize=10)

                # If an intersection parameter is negative, then the intersection point
                # is on the extended infinite line of the first segment, before the segment.
                # This segment is then extended here to the intersection point to simplify
                # the following calculations.
                if tL < 0.:
                    borderL.verts[0] = iP
                    tL = 0.
                if tR < 0.:
                    borderR.verts[0] = iP
                    tR = 0.

                # Store the polylines of the borders for later use.
                bordersL[lineIndx1] = borderL
                bordersR[lineIndx2] = borderR

                if doDebug:
                    borderL.plot('#00ff00')
                    borderR.plot('#ff0000')

                # Get an estimate of fillet radius between these segments
                radius1 = estFilletRadius(line1.section.originalSection.category,line1.section.originalSection.tags)
                radius2 = estFilletRadius(line2.section.originalSection.category,line2.section.originalSection.tags)
                radius = min(radius1,radius2)

                # Reduce fillet radius continuously when angle gets smaller
                cosAngle = uVL.dot(uVR)
                if cosAngle > 0.1:
                    radius *= cosAngle

                # Starting with this radius, we try to create a fillet within the intersected
                # segments. If this is not possible, the radius is decreased.
                # <lengthL> and <lengthR> are the lengths on the intersected segments between
                # their start point and the intersection
                lengthL = borderL.segmentLength(floor(tL))
                lengthR = borderR.segmentLength(floor(tR))

                while True: # Loop to find a valid radius
                    # Try to find an arc with radius <radius> and using the segments <iP>-<uVLY and
                    # <iP>-<uVR> as tangents. In the result, <origin> is the circle origin of the arc,
                    # <legEndL> and <legEndR> are the tangent vertices and <legLength> is the distance
                    # between the intersection point and the tangent vertices.
                    origin, tLegEndL, tLegEndR, legLength = filletArc(iP, uVL, uVR, radius)
                    if origin is None: # No arc possible, segments are almost parallel.
                        break

                    # To check whether the leg ends we are on the intersecting segments, the line
                    # parameters of the tangent vertices are computed (see PolyLine.py). For a 
                    # valid fillet, they must be both less than 1.
                    tLegL = tL%1. + legLength / lengthL
                    tLegR = tR%1. + legLength / lengthR
                    if tLegL < 1. and tLegR < 1.:
                        # OK, fillet fits between these segments using this radius
                        break

                    # if not, reduce radius
                    radius = 0.9 * radius

                # If a fillet arc has been found.
                if origin:
                    if doDebug:
                        plt.plot(tLegEndL[0],tLegEndL[1],'go',markersize=6)
                        plt.plot(tLegEndR[0],tLegEndR[1],'ro',markersize=6)

                    # Create the vertices for the fillet.
                    filletVerts = filletLine(origin, tLegEndL, tLegEndR, radius)

                    # Maybe, the borders and the centerline of a way-segment do not have
                    # the same shape. To get reliable line parameters, the tangent vertices
                    # are projected onto the centerline and the line parameter of this 
                    # projection is computed (<t1> for <line1> and <t2> for <line2>).
                    t1 = line1.polyline.projectOrthogonal(tLegEndL,True)
                    if doDebug:
                        q = line1.polyline.verts[-1]
                        plt.text(q[0],q[1],str(lineIndx1))
                    t2 = line2.polyline.projectOrthogonal(tLegEndR,True)
                    if doDebug:
                        q = line2.polyline.verts[-1]
                        plt.text(q[0],q[1],str(lineIndx2))

                    # store the maximum of these intersection parameters,
                    # indexed by the line index.
                    tSmax[lineIndx1] = max(tSmax[lineIndx1],t1)
                    tSmax[lineIndx2] = max(tSmax[lineIndx2],t2)

                    # Store the fillet vertices for <line1>.
                    filletVertsList[lineIndx1] = filletVerts

                    if doDebug:
                        for v1,v2 in zip(filletVerts[:-1],filletVerts[1:]):
                            plt.plot([v1.x,v2.x],[v1.y,v2.y],'c:')
                else:
                    pass    # No possible fillet found
                    print('No possible fillet found')

            except (Exception,ValueError) as e:
                import traceback
                traceback.print_exception(type(e), e, e.__traceback__)
                print('Exception on fillet computation')
        try:
            polygon = []
            # Construct the intersection polygon
            for (lineIndx1, lineIndx2), (line1,line2) in zip(cyclePair(range(N)),cyclePair(self.outWayLines)):
            # <line1> and <line2> are consecutive outgpoing polylines in counter-clockwise order of the
            # intersection. These are the centerlines of the way-sections, <line1> is to the right of <line2>.
                filletVerts = filletVertsList[lineIndx1]
                tS1 = tSmax[lineIndx1]  # Line parameters of trim points.
                tS2 = tSmax[lineIndx2]

                # The way sections will have to be trimmed to the length not occupied by
                # the polygon of the intersection area. Depending on the direction of
                # the way section as outgoing polyline (forward => fwd==True), these
                # values have to be stored in the corresponding way section object.
                if line1.fwd:
                    line1.section.trimS = tS1
                    # pS = line1.section.polyline.t2v(tS1,True)
                    # plt.plot(pS[0],pS[1],'gx',markersize=15)
                    # plt.text(pS[0],pS[1],' %5.2f'%(tS1))
                else:
                    line1.section.trimT = tS1
                    # pS = line1.section.polyline.t2v(tS1,False)
                    # plt.plot(pS[0],pS[1],'rx',markersize=15)
                    # plt.text(pS[0],pS[1],' %5.2f'%(tS1))

                # Starting with the left border of the centerline <line1>.
                # Given the line parameter of the trim point on the centerline, we can 
                # compute the point <pL> perpendicularly offset to the left.
                w1 = line1.section.leftWidth if line1.fwd else line1.section.rightWidth
                w2 = line2.section.rightWidth if line2.fwd else line2.section.leftWidth
                pL = line1.polyline.offsetPointAt(tS1,w1)
                if filletVerts:
                    # By a simple check of coordinate differences, we may decide if
                    # <pL> is a vertex of the fillet.
                    dL = sum(pL-filletVerts[0])
                    if abs(dL) > 1.e-4:
                        # If not on the fillet, add to polygon
                        polygon.append(pL)
                    # then extend polygon by fillet verts.
                    polygon.extend(filletVerts)
                else:
                    polygon.append(pL)

                # The same procedure is then done for the right border of the centerline <line2>
                pR = line2.polyline.offsetPointAt(tS2,-w2)
                if filletVerts:
                    dR = sum(pR-filletVerts[-1])
                    if abs(dR) > 1.e-4:
                        polygon.append(pR)
                else:
                    polygon.append(pR)

        except (Exception,ValueError) as e:
            import traceback
            traceback.print_exception(type(e), e, e.__traceback__)
            print('Exception on construction of intersection polygon')

        return polygon

    def findTransitionPoly(self):
        # Transitions are intersections of order==2 and have an incoming and
        # an outgoing way. They have to be processed before the others, as
        # way widths may be altered due to turning lanes.
        line1, line2 = self.outWayLines
        outLine, inLine = (line1, line2) if line1.fwd else (line2, line1)
        inTags, outTags = inLine.section.originalSection.tags, outLine.section.originalSection.tags

        # Do we have turn lanes? They are only possible in the outLine.
        if 'turn:lanes' in outTags:
            # There is no transition polygon required. The outgoing way section
            # becomes eventually a turning lane.
            laneDescs = outTags['turn:lanes'].split('|')
            leftTurnLanes = sum(1 for tag in laneDescs if 'left' in tag)
            rightTurnLanes = sum(1 for tag in laneDescs if 'right' in tag)
            if leftTurnLanes or rightTurnLanes:
                leftWidthDifference = outLine.section.leftWidth - inLine.section.leftWidth
                rightWidthDifference = outLine.section.rightWidth - inLine.section.rightWidth
                if leftWidthDifference or rightWidthDifference:
                    outLine.section.turnParams = [leftWidthDifference,rightWidthDifference]

        # Prepare transition polygon 
        if outLine.section.turnParams:
            transitionWidth = 0.5
        else:
            inWidth = inLine.section.leftWidth + inLine.section.leftWidth
            outWidth = outLine.section.leftWidth + outLine.section.leftWidth
            widthDiff = abs(inWidth-outWidth)
            transitionWidth = min(10.,max(1.,2 * widthDiff))
        inT = inLine.polyline.length2t(transitionWidth,inLine.fwd)
        outT = outLine.polyline.length2t(transitionWidth,outLine.fwd)
        if inT is None or outT is None:
            pass
            # plt.plot(self.position[0],self.position[1],'rx',markersize=20,zorder=950)
        else:
            assert inT is not None
            assert outT is not None

            # compute cornerpoints
            # inLeft = inLine.polyline.offsetPointAt(inT,inLine.section.leftWidth)
            # inRight = inLine.polyline.offsetPointAt(inT,-inLine.section.rightWidth)
            inLeft = inLine.polyline.offsetPointAt(inT,inLine.section.rightWidth)
            inRight = inLine.polyline.offsetPointAt(inT,-inLine.section.leftWidth)
            outLeftWidth = inLine.section.leftWidth if outLine.section.turnParams else outLine.section.leftWidth
            outRightWidth = inLine.section.rightWidth if outLine.section.turnParams else outLine.section.rightWidth
            outLeft = outLine.polyline.offsetPointAt(outT,outLeftWidth)
            outRight = outLine.polyline.offsetPointAt(outT,-outRightWidth)

            poly = []#[inLeft,inRight,outLeft,outRight]
            import numpy as np
            t0 = np.linspace(0,1,10)#[0.0,0.2,0.4,0.6,0.8,1.0]
            p0, p1 = inRight, outLeft
            v0, v1 = p0 - inLine.polyline.unitEndVec(False)*transitionWidth, p1 - outLine.polyline.unitEndVec(False)*transitionWidth
            for t in t0:
                sp = spline_4p( t, v0, p0, p1, v1 )
                poly.append( sp )
            p0, p1 = outRight, inLeft
            v0, v1 = p0 - outLine.polyline.unitEndVec(False)*transitionWidth, p1 - inLine.polyline.unitEndVec(False)*transitionWidth
            for t in t0:
                sp = spline_4p( t, v0, p0, p1, v1 )
                poly.append( sp )


            if inLine.fwd:
                inLine.section.trimS = inT
            else:
                inLine.section.trimT = inT
            if outLine.fwd:
                outLine.section.trimS = outT
            else:
                outLine.section.trimT = outT

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

    # segList.append(p2)
    vertList.reverse()
    return vertList


