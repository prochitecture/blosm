from itertools import product
from lib.pygeos.shared import Coordinate
from lib.CompGeom.WayOffsetGenerator import WayOffsetGenerator
from itertools import tee

from math import floor, ceil
from mathutils import Vector
import matplotlib.pyplot as plt

# helper functions -----------------------------------------------
def pairs(iterable):
	# s -> (s0,s1), (s1,s2), (s2, s3), ...
    p1, p2 = tee(iterable)
    next(p2, None)
    return zip(p1,p2)
# ----------------------------------------------------------------

class PolyLine():
    def __init__(self, vertList):
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
        return len(self.verts)

    def clone(self):
         return PolyLine([v for v in self.verts])

    def reversed(self):
        return PolyLine([v for v in reversed(self.verts)])

    def __getitem__(self, key):
        return self.verts[key]

    def __iter__(self):
        return iter(self.verts)

    def segments(self):
        return [(v1,v2) for v1, v2 in zip(self.verts[:-1],self.verts[1:])]

    def segment(self,i):
        if i>len(self.verts)-2:
            i = len(self.verts)-2
        return (self.verts[i],self.verts[i+1])

    def segmentLength(self,i):
        assert i<len(self.verts)+1
        return (self.verts[i+1]-self.verts[i]).length

    def unitVectors(self):
        return [(v2-v1)/(v2-v1).length for v1, v2 in zip(self.verts[:-1],self.verts[1:])]

    def projectOrthogonal(self,p,fwd):
        s = slice(None,None,1) if fwd else slice(None,None,-1)
        pNr = 0
        for p1,p2 in pairs(self.verts[s]):
            segV = p2-p1
            t = (p-p1).dot(segV)/segV.dot(segV)
            if 0. <= t < 1.:
                p0 = p1 + t*segV
                break
            pNr += 1
        t += pNr
        return t if t<=len(self.verts) else None

    def offsetPointAt(self,t,width):
        if t == 0.:
            v0 = self.verts[0]
            v1 = self.verts[1]
        else:
            v0 = self.verts[floor(t)]
            v1 = self.verts[ceil(t)]
        p = v0 + (v1 - v0) * (t % 1)
        d = v1-v0
        u = d/d.length
        pOffs = p + Vector((-u[1],u[0])) * width
        return pOffs


    def t2v(self,t,fwd):
        # computes the vertex given by the intersection parameter t, where t is
        # the sum of the index of the start vertex and the percentage
        # of its following segment.
        s = slice(None,None,1) if fwd else slice(None,None,-1)
        v0 = self.verts[s][floor(t)]
        v1 = self.verts[s][ceil(t)]
        return v0 + (v1 - v0) * (t % 1)

    def trimS(self,t):
        if t == 0.: return
        n = ceil(t)
        assert(n<len(self.verts))
        p = self.t2v(t,True)
        self.verts = [p] + self.verts[n:]

    def slice(self,tS,tT):
        N = len(self.verts)
        if tS==0.:
            pS = []
            iS = -1
        else:
            iS = floor(tS)
            v0,v1 = self.verts[iS], self.verts[iS+1]
            pS = [v0 + (v1 - v0) * (tS % 1)]

        if tT==0.:
            pT = []
            iT = N-1
        else:
            t = N-1-tT
            iT = floor(t)
            v0,v1 = self.verts[iT], self.verts[iT+1]
            pT = [v0 + (v1 - v0) * (t % 1)]

        inter = [v for v in self.verts[iS+1:iT+1]]
        return PolyLine(pS+inter+pT)
        

    def clone(self):
        return PolyLine(self.verts)

    def parallelOffset(self, distance):
        offsetter = WayOffsetGenerator()
        offsetVerts = offsetter.parallel_offset(self.verts,distance)
        return PolyLine(offsetVerts)

    def buffer(self,leftW,rightW):
        offsetter = WayOffsetGenerator()
        offsetL = offsetter.parallel_offset(self.verts,leftW)
        offsetR = offsetter.parallel_offset(self.verts,-rightW)
        offsetR.reverse()
        offsetVerts = offsetL + offsetR
        offsetPoly = PolyLine(offsetVerts)
        return offsetPoly

    @staticmethod
    def intersectInfinite(p1,p2,p3,p4):
        d1, d2 = p2-p1, p4-p3
        cross = d1.cross(d2)
        dot = d1.dot(d2)
        if abs(cross/dot) < 0.0175: # tangent of ~1°, almost parallel
            return None
        d3 = p1-p3
        t1 = (d2[0]*d3[1] - d2[1]*d3[0])/cross
        t2 = (d1[0]*d3[1] - d1[1]*d3[0])/cross
        p = p1 + d1*t1
        return p, t1, t2

    @staticmethod
    def intersect(seg1,seg2,i1,i2):
        p1, p2 = seg1.segment(i1)
        p3, p4 = seg2.segment(i2)
        d1, d2 = p2-p1, p4-p3 
        denom = d2.y*d1.x - d2.x*d1.y
        if denom == 0.: # parallel
            return None
        d3 = p1 - p3
        ua = (d2.x*d3.y - d2.y*d3.x)/denom 
        ub = (d1.x*d3.y - d1.y*d3.x)/denom              
        c = p1 + d1*ua
        p5 = p2 + d1*5
        p6 = p4 + d2*5
        # plt.plot([p1.x,p5.x],[p1.y,p5.y],'k:')
        # plt.plot([p3.x,p6.x],[p3.y,p6.y],'k:')
        return c, i1+ua, i2+ub, d1/d1.length, d2/d2.length


    @staticmethod
    def intersection(poly1,poly2):
        segs1, segs2 = poly1.segments(), poly2.segments()
        found = False

        # Intersections for the first few segments are most probable.
        # Construct indices for combinations so that these are first tested.
        combinations = list(product(range(len(segs1)),range(len(segs2))))
        combinations = sorted(combinations,key=lambda x: sum(x))
        for i1,i2 in combinations:
            p1, p2 = segs1[i1]
            p3, p4 = segs2[i2]
            isect = PolyLine.intersectInfinite(p1,p2,p3,p4)
            if isect is None: continue

            iP, t1, t2 = isect
            if t1 < 0. or t1 > 1.: continue # out of segment
            if t2 < 0. or t2 > 1.: continue # out of segment
            found = True
            # p5 = p2 + d1*5
            # p6 = p4 + d2*5
            # plt.plot([p1.x,p5.x],[p1.y,p5.y],'k:')
            # plt.plot([p3.x,p6.x],[p3.y,p6.y],'k:')
            break # valid intersection found

        if not found:
            # try to intersect the infinite lines of the first segments
            p1, p2 = segs1[0]
            p3, p4 = segs2[0]
            isect = PolyLine.intersectInfinite(p1,p2,p3,p4)
            if isect is None:
                return None
            iP, t1, t2 = isect
            if t1>1. or t2>1.:  # want to have only intersections before start of segments
                return None
            d1, d2 = p2-p1, p4-p3
            return iP, t1, t2, d1/d1.length, d2/d2.length

        d1, d2 = p2-p1, p4-p3
        return iP, i1+t1, i2+t2, d1/d1.length, d2/d2.length

    def plot(self,color):
        for v1,v2 in self.segments():
            plt.plot([v1.x,v2.x],[v1.y,v2.y],color)
