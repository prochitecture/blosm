# from lib.pygeos.shared import Coordinate
from way.PolyLine import PolyLine

from math import sin, cos, atan2, pi, sqrt, floor
from mathutils import Vector
from itertools import tee,islice, cycle
from itertools import cycle

def cyclePair(lst):
    prevs, nexts = tee(lst)
    prevs = islice(cycle(prevs), len(lst) - 1, None)
    return zip(prevs,nexts)

class OutLine():
    def __init__(self,section,polyline,fwd):
        self.section = section
        self.polyline = polyline
        self.fwd = fwd

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
        # section <waySection>, the attribuet <fwd> in <Outline> is set to True, and False,
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
        # The outgoing polylines in <self.outWayLines> of class Outline are expected 
        # to be sorted in counter-clockwise order.
        N = self.order              # number of outgoing polylines
        tParams = [0.]*N            # Parameter along a polyline, sum of the index of the first
                                    # vertex of a segmentand the percentage of the length of this segment.
        allFilletVerts = [None]*N
        bordersL = [None]*N         # list of left borders of outgoing ways
        bordersR = [None]*N         # list of right borders of outgoing ways
        # plt.close()

        for (lineIndx1, lineIndx2), (line1,line2) in zip(cyclePair(range(N)),cyclePair(self.outWayLines)):
        # for lineIndx1, lineIndx2, line1,line2 in zip(range(N),range(1,N+1),self.outWayLines,self.outWayLines[1:]+[self.outWayLines[0]]):
            try:
                # plt.close()
                lineIndx2 %= N  # circular index for outgoing 

                # line1['line'].plot('k')
                # e = line1['line'].verts[1]
                # plt.text(e[0],e[1],str(indx1))
                # line2['line'].plot('k')
                # e = line2['line'].verts[1]
                # plt.text(e[0],e[1],str(indx2))

                # Compute the left and right border of the outgoing polylines.
                # For way-loops, the forward and backward polyline is in consecutive order
                # in <self.outWayLines>. For loop-ways, a parallel offset can't be computed
                # correctly, therefore, only the first segments are used.
                if line1.section.isLoop and line2.section.isLoop:
                    trimmed1 = PolyLine(line1.polyline[:2])
                    trimmed2 = PolyLine(line2.polyline[:2])
                    borderL = trimmed1.parallelOffset( line1.section.leftWidth)
                    borderR = trimmed2.parallelOffset(-line2.section.rightWidth)
                else:
                    borderL = line1.polyline.parallelOffset( line1.section.leftWidth)
                    borderR = line2.polyline.parallelOffset(-line2.section.rightWidth)

                # Find intersection between left border <borderL> of <line1> 
                # and right border <borderR> of <line2>
                # In <isectParams>, <iP> is the intersection point, <tL> and <tR> are the
                # intersection parameters for left and right line. They consist of the 
                # sum of the index of the first vertex of the intersected segments
                # and the percentage (0 ... 1) within the length of this segment.
                # <uVL> and <uVR> are the unit vectors of the intersected segments
                isectParams = PolyLine.intersection(borderL,borderR)
                if isectParams is None:
                    # print('no intersection')
                    continue
                iP,tL,tR,uVL,uVR = isectParams

                # plt.plot(iP[0],iP[1],'mo',markersize=6)

                # If the intersection parameter is negative, then the intersection point
                # is on the extended infinite line of the first segment. This segment is
                # extended to the intersection point to simplify the following calculations.
                if tL < 0.:
                    borderL.verts[0] = iP
                    tL = 0.
                if tR < 0.:
                    borderR.verts[0] = iP
                    tR = 0.

                # Store them for later use.
                bordersL[lineIndx1] = borderL
                bordersR[lineIndx2] = borderR

                # borderL.plot('#00ff00')
                # borderR.plot('#ff0000')

                # get first estimate of fillet radius between these segments
                # TODO: Replace this by better radius values.
                maxVel = 50.
                radius = maxVel/9.81/0.7/2. #???

                # Starting with this radius, we try to create a fillet within the intersected
                # segments. If this is not possible, the radius is decreased.
                # <lengthL> and <lengthR> are the lengths on the intersected segments between
                # their start point and the intersection
                lengthL = borderL.segmentLength(floor(tL))
                lengthR = borderR.segmentLength(floor(tR))

                while True: # Loop to find possible radius
                    # Try to find an arc using the segments given by <iP>, <uVLY, <uVR> as tangents.
                    # <origin> is the circle origin of the arc, <legEndL> and <legEndR> are the 
                    # tangent points and <legLength> is the distance between the intersection point
                    # and the tangent points.
                    origin, tLegEndL, tLegEndR, legLength = filletArc(iP, uVL, uVR, radius)
                    if origin is None: # No arc possible, segments are almost parallel.
                        break

                    # Check whether the leg ends we are on the intersecting segments.
                    tLegL = tL%1. + legLength / lengthL
                    tLegR = tR%1. + legLength / lengthR
                    if tLegL < 1. and tLegR < 1.:
                        # OK, fillet fits with these segments and this radius
                        break

                    # if not, reduce radius
                    radius = 0.9 * radius

                    # ax = plt.gca()
                    # ax.add_artist(plt.Circle(
                    #     (isectP.x,isectP.y),
                    #     10,
                    #     alpha=0.1,
                    #     color='g',
                    #     zorder=100
                    # )) 

                # If a fillet arc has been found.
                if origin:
                    # plt.plot(tLegEndL[0],tLegEndL[1],'go',markersize=6)
                    # plt.plot(tLegEndR[0],tLegEndR[1],'ro',markersize=6)

                    # Create the vertices for the fillet.
                    filletVerts = filletLine(origin, tLegEndL, tLegEndR, radius)

                    t1 = tL + legLength / lengthL
                    t2 = tR + legLength / lengthR
                    # t1 = line1.polyline.projectOrthogonal(tLegEndL,True)
                    # t2 = line2.polyline.projectOrthogonal(tLegEndR,True)
                    test=1

                    # Project the tangent points onto the center line and
                    # get their intersection parameters.
                    # t1 = line1.polyline.projectOrthogonal(tLegEndL,True)
                    # t2 = line2.polyline.projectOrthogonal(tLegEndR,True)

                    # p1 = line1.polyline.t2v(t1,line1.fwd)
                    # p2 = line2.polyline.t2v(t2,line2.fwd)
                    # plt.plot(tLegEndL[0],tLegEndL[1],'co',markersize=6)
                    # plt.plot(tLegEndR[0],tLegEndR[1],'mo',markersize=6)
                    # plt.plot(p1[0],p1[1],'go',markersize=6)
                    # plt.plot(p2[0],p2[1],'ro',markersize=6)
                    # line1.polyline.plot('k')
                    # line2.polyline.plot('k')

                    # store the maximum of the intersection parameters
                    tParams[lineIndx1] = max(tParams[lineIndx1],t1)
                    tParams[lineIndx2] = max(tParams[lineIndx2],t2)

                    # Store the fillet vertices for <line1>.
                    allFilletVerts[lineIndx1] = filletVerts

                    # for v1,v2 in zip(filletVerts[:-1],filletVerts[1:]):
                    #     plt.plot([v1.x,v2.x],[v1.y,v2.y],'c:')
                else:
                    pass
                    # print('no fillet')

            except (Exception,ValueError) as e:
                import traceback
                traceback.print_exception(type(e), e, e.__traceback__)
                print('exception 1')
        try:
        # if True:
            polygon = []
            # Construct the intersection polygon
            for (indx1, indx2), (line1,line2) in zip(cyclePair(range(N)),cyclePair(self.outWayLines)):
                filletVerts = allFilletVerts[indx1]
                tParamL = tParams[indx1]
                tParamR = tParams[indx2]

                if line1.fwd:
                    line1.section.trimS = max(tParamL,line1.section.trimS)
                else:
                    line1.section.trimT = max(tParamL,line1.section.trimT)

                # we start with the offset point of the left border of line1
                pL = line1.polyline.offsetPointAt(tParamL,line1.section.leftWidth)
                if filletVerts:
                    dL = sum(pL-filletVerts[0])
                    if abs(dL) > 1.e-4:
                        polygon.append(pL)
                    polygon.extend(filletVerts)
                else:
                    polygon.append(pL)

                

                pR = line2.polyline.offsetPointAt(tParamR,-line2.section.rightWidth)
                if filletVerts:
                    dR = sum(pR-filletVerts[-1])
                    if abs(dR) > 1.e-4:
                        polygon.append(pR)
                else:
                    polygon.append(pR)
                # bordersL[indx1].plot('g')
                # bordersR[indx2].plot('r')
                
            # x = [n[0] for n in polygon]
            # y = [n[1] for n in polygon]
            # plt.fill(x,y,'#cf4b23',alpha = 0.4,zorder = 900)
            # plt.plot(x,y,'r')

        except (Exception,ValueError) as e:
            import traceback
            traceback.print_exception(type(e), e, e.__traceback__)
            print('exception 2')

        return polygon

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


