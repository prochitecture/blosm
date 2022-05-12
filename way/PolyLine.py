from itertools import product
from lib.pygeos.shared import Coordinate, CAP_STYLE, JOIN_STYLE
from lib.pygeos.geom import GeometryFactory
from lib.CompGeom.WayOffsetGenerator import WayOffsetGenerator

from math import floor, ceil, pi, asin, sin, tan
from mathutils import Vector
import matplotlib.pyplot as plt

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

    def t2v(self,t):
        # computes the vertex given by the intersection parameter t, where t is
        # the sum of the index of the start vertex and the percentage
        # of its following segment.
        v1 = self.verts[floor(t)]
        v2 = self.verts[ceil(t)]
        return v1 + (v2 - v1) * (t % 1)

    def clone(self):
        return PolyLine(self.verts)

    def parallelOffset(self, distance):
        geosF = GeometryFactory()
        offsetter = WayOffsetGenerator()
        offsetVerts = offsetter.parallel_offset(self.verts,distance)
        return PolyLine(offsetVerts)

    def trimLength(self,p,nv):
        segs = self.segments()
        p3,p4 = p,p+nv
        c = None
        for i1 in range(len(segs)):
            p1, p2 = segs[i1]
            d1, d2 = p2-p1, p4-p3 
            denom = d2.y*d1.x - d2.x*d1.y
            if denom == 0.: continue # parallel
            d5 = p1 - p3
            ua = (d2.x*d5.y - d2.y*d5.x)/denom               
            if ua < 0. or ua > 1.: continue # out of range
            c = p1 + d1*ua
            break
        return c, i1+ua

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
