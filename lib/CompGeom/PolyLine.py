from ctypes.wintypes import tagSIZE
from itertools import product
from lib.pygeos.shared import Coordinate
from lib.CompGeom.OffsetGenerator import OffsetGenerator
from itertools import tee
import numpy as np

from math import floor, ceil
from mathutils import Vector

import matplotlib.pyplot as plt

# ----------------------------------------------------------------
# PolyLine holds an ordered list of vertices that describe a line.
# The vertices are of the class mathutils.Vector. Line parameters
# tS and tT describe position onto the line. They are a compound of
# the last vertex number and the percentage of the following segment.
# tS measures the distance from the start vertex, while tE measures
# it from the end vertex, as for example:
#
#          tS=1.4                       tE= 2.8
#        ------------>         <--------------------------
#       |             |       |                           |
#       |             V       V                           |
#       o---------o---x-----o-x-------o---------o---------o
#       0         1         2         3         4         5.
#
# ----------------------------------------------------------------

# helper functions -----------------------------------------------
def pairs(iterable):
	# s -> (s0,s1), (s1,s2), (s2, s3), ...
    p1, p2 = tee(iterable)
    next(p2, None)
    return zip(p1,p2)

def unitDoubleParabola(nrPoints):
    x = np.linspace(0.,1.,nrPoints)
    return np.where(x <= 0.5, 2.*x*x, 1.-2.*(1.-x)*(1.-x))

def unitHalfParabola(nrPoints):
    x = np.linspace(0.,1.,nrPoints)
    return np.where(x <= 0.5, 1.-4*(0.5-x)*(0.5-x), 1.)

# ----------------------------------------------------------------

class PolyLine():
    def __init__(self, vertList):
        # Initilaizes an instance from a list of vertices of the
        # classes mathutils.Vector, pyGEOS Coordinate or of
        # tuples of two floats for x and y values.
        self.verts = None
        if isinstance(vertList[0],Vector):
            self.verts = vertList
        elif isinstance(vertList[0],Coordinate):
            self.verts = [Vector((v.x,v.x)) for v in vertList]
        elif isinstance(vertList[0],tuple):
            self.verts = [Vector(v) for v in vertList]
        else:
            raise TypeError('PolyLine: Unexpected input type.')

    def __len__(self):
        # Return the length of the line
        return len(self.verts)

    def clone(self):
        # Create a copy of the PolyLine
        return PolyLine([v for v in self.verts])

    def reversed(self):
        # Create a reversedcopy of the PolyLine
        return PolyLine([v for v in reversed(self.verts)])

    def __getitem__(self, key):
        # Acces a vertex by index
        return self.verts[key]

    def __iter__(self):
        # Iterator of the vertices list
        return iter(self.verts)

    def segments(self):
        # Returns a list of vertex tuples for all segments of the line
        return [(v1,v2) for v1, v2 in pairs(self.verts)]

    def segment(self,i):
        # Returns the ith segment of the line
        if i>len(self.verts)-2:
            i = len(self.verts)-2
        return (self.verts[i],self.verts[i+1])

    def segmentLength(self,i):
        # Returns the length of the ith segment of the line
        assert i<len(self.verts)+1
        return (self.verts[i+1]-self.verts[i]).length

    def unitVectors(self):
        # Returns a list of the unit vectors of the line's segments
        return [(v2-v1)/(v2-v1).length for v1, v2 in pairs(self.verts)]

    def projectOrthogonal(self,p,fwd):
        # Projects the point <p> onto the line. The return value
        # is the line parameter tS, when <fwd> is true, or tE, 
        # when <fwd> is False.
        s = slice(None,None,1) if fwd else slice(None,None,-1)
        pNr = 0
        for p1,p2 in pairs(self.verts[s]):
            segV = p2-p1
            t = (p-p1).dot(segV)/segV.dot(segV)
            if  t < 1.:
                p0 = p1 + t*segV
                break
            pNr += 1
        t += pNr
        return t if t<=len(self.verts) else None

    def offsetPointAt(self,tS,dist):
        # Computes the position of the point on the line, given by
        # the line parameter <tS> and offfsets this point perpendicular
        # to the segment it lies on by the distance <dist>. The result
        # is on the left of the line when <dist> is positive and on
        # their right else.
        t0 = floor(tS)
        if t0 >= len(self.verts)-1:
            v0 = self.verts[-2]
            v1 = self.verts[-1]
            p0 = v1
        else:
            v0 = self.verts[t0]
            v1 = self.verts[t0+1]
            p0 = v0
        p = p0 + (v1 - v0) * (tS % 1)
        d = v1-v0
        u = d/d.length
        pOffs = p + Vector((-u[1],u[0])) * dist
        return pOffs

    def unitEndVec(self,fwd):
        s = slice(None,None,1) if fwd else slice(None,None,-1)
        d = self.verts[s][1]-self.verts[s][0]
        return d/d.length

    def length2t(self,length,fwd):
        t = 0.
        sumLength = 0.
        s = slice(None,None,1) if fwd else slice(None,None,-1)
        for i in range(len(self.verts)-1):
            segLength = (self.verts[s][i+1] - self.verts[s][i]).length
            if sumLength + segLength < length:
                t += 1.
                sumLength += segLength
            else:
                t1 = (length-sumLength)/segLength
                t += (length-sumLength)/segLength
                return t
        return None
 

    def t2v(self,t,fwd):
       # Computes the vertex given by the line parameter <t>, where 
        # <t> is the line parameter tS, when <fwd> is true, or tE, 
        # when <fwd> is False.
        s = slice(None,None,1) if fwd else slice(None,None,-1)
        v0 = self.verts[s][floor(t)]
        v1 = self.verts[s][ceil(t)]
        return v0 + (v1 - v0) * (t % 1)

    def slice(self,tS,tE):
        # Returns the slice of the line between and including the
        # vertices given by the line parameters <tS> and <tE>
        N = len(self.verts)
        if tS==0.:
            pS = []
            iS = -1
        else:
            iS = floor(tS)
            v0,v1 = self.verts[iS], self.verts[iS+1]
            pS = [v0 + (v1 - v0) * (tS % 1)]

        if tE==0.:
            pE = []
            iE = N-1
        else:
            t = N-1-tE
            iE = floor(t)
            v0,v1 = self.verts[iE], self.verts[iE+1]
            pE = [v0 + (v1 - v0) * (t % 1)]

        inter = [v for v in self.verts[iS+1:iE+1]]
        return PolyLine(pS+inter+pE)
        
    def parallelOffset(self, dist):
        # Returns the PolyLine that is perpendicular offset by
        # the distance <dist>. The result s on the left of
        # the line when <dist> is positive and on their right else.
        offsetter = OffsetGenerator()
        offsetVerts = offsetter.parallel_offset(self.verts,dist)
        return PolyLine(offsetVerts)

    def buffer(self,leftW,rightW):
        # Expands the line to a polygon, to the left by the distance
        # <leftW> and to the right by <rightW>.
        offsetter = OffsetGenerator()
        offsetL = offsetter.parallel_offset(self.verts,leftW)
        offsetR = offsetter.parallel_offset(self.verts,-rightW)
        offsetR.reverse()
        offsetVerts = offsetL + offsetR
        offsetPoly = PolyLine(offsetVerts)
        return offsetPoly

    def parabolicOffset(self,dist,delta,nrPoints = 30):
        # Prepare transform of segment lengths along polyline to line parameters <t>
        segLengths = np.array([0]+[(v1-v2).length for v1,v2 in pairs(self.verts)])
        sumLengths = np.cumsum(segLengths)
        t = np.array([i for i in range(len(self.verts))])

        # Create <samples> and transform to line parameter samples <tSamples>
        samples = np.linspace(sumLengths[0],sumLengths[-1],nrPoints)
        tSamples = np.interp(samples,sumLengths,t)

        # Create offset samples offSamples, shifted by dist and added double parabola
        # with width of delta.
        offVerts = []
        para = unitHalfParabola(nrPoints)#unitDoubleParabola(nrPoints)
        # plt.close()
        # plt.plot(para)
        # plt.show()
        for i in range(nrPoints):
            o = self.offsetPointAt(tSamples[i],dist+para[i]*delta)
            offVerts.append(o)
        return offVerts

    def parabolicBuffer(self,leftW,rightW,leftD,rightD ):
        p = self.verts[0]
        # plt.plot(p[0],p[1],'ro',markersize=10,zorder=950)
        offsetter = OffsetGenerator()
        if leftD != 0.:
            offsetL = self.parabolicOffset(leftW-leftD,leftD)
        else:
            offsetL = offsetter.parallel_offset(self.verts,leftW)
        if rightD != 0.:
            offsetR = self.parabolicOffset(-(rightW-rightD),-rightD)
        else:
            offsetR = offsetter.parallel_offset(self.verts,-rightW)

        offsetR.reverse()
        offsetVerts = offsetL + offsetR
        offsetPoly = PolyLine(offsetVerts)
        return offsetPoly


    @staticmethod
    def intersectInfinite(p1,p2,p3,p4):
        # Finds the intersection of two infinite lines, line 1 given by
        # the vertices <p1> and <p2> and line 2 the vertices <p3> and <p4>.
        # Returns the intersection point <p> ant the line parameters <tS1>
        # and <ts2> relative to the segments given. When the lines are almost
        # paralllel (angle between the less tha 1°), None is returned.
        d1, d2 = p2-p1, p4-p3
        cross = d1.cross(d2)
        dot = d1.dot(d2)
        if abs(cross/dot) < 0.0175: # tangent of ~1°, almost parallel
            return None
        d3 = p1-p3
        tS1 = (d2[0]*d3[1] - d2[1]*d3[0])/cross
        tS2 = (d1[0]*d3[1] - d1[1]*d3[0])/cross
        p = p1 + d1*tS1
        return p, tS1, tS2

    @staticmethod
    def intersection(poly1,poly2):
        # Finds the intersection between two PolyLines <poly1> and <poly2>.
        # Returns the intersection point <iP>, the line parameters <tS1> and
        # <tS2> of the intersection on the lines and the unit vectors <uV1>
        # and <uV2> of the intersected segments.
        # When there is no intersection between the PolyLine segments, the
        # intersection may be on the infinite lines through the first segments.
        # In this case, <tS1> and/or <tS2> may be negative. Whene these infinte
        # lines are almost paralllel (angle between the less tha 1°), or when
        # they intersect after the segments of the polylines, None is returned.
        segs1, segs2 = poly1.segments(), poly2.segments()
        found = False

        # Intersections for the first few segments are most probable.
        # The combinations of the indices are constructed so that these
        # are first tested.
        combinations = list(product(range(len(segs1)),range(len(segs2))))
        combinations = sorted(combinations,key=lambda x: sum(x))
        for i1,i2 in combinations:
            p1, p2 = segs1[i1]
            p3, p4 = segs2[i2]
            isect = PolyLine.intersectInfinite(p1,p2,p3,p4)
            if isect is None: continue

            iP, t1, t2 = isect
            # Accept intersection of infinite lines (negative line parameters)
            # only between first segments.
            if (i1,i2) == (0,0): 
                if t1 > 1. or t2 > 1.:
                    continue # out of segment
            else:
                if t1 < 0. or t1 > 1. or t2 < 0. or t2 > 1.:
                    continue # out of segment
            found = True
            break # valid intersection between PolyLines found

        if found:
            d1, d2 = p2-p1, p4-p3
            return iP, i1+t1, i2+t2, d1/d1.length, d2/d2.length
        else:
            return None

    def plot(self,color):
        for v1,v2 in self.segments():
            plt.plot([v1.x,v2.x],[v1.y,v2.y],color)
            plt.plot(v1.x,v1.y,'k.')
            plt.plot(v2.x,v2.y,'k.')
