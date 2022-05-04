from itertools import product
from lib.pygeos.shared import Coordinate, CAP_STYLE, JOIN_STYLE
from lib.pygeos.geom import GeometryFactory
from math import fmod
from mathutils import Vector
import matplotlib.pyplot as plt

class PolyLine():
    def __init__(self, coords):
        self.coords = None
        if isinstance(coords[0],Coordinate):
            self.coords = coords
        elif isinstance(coords[0],Vector):
            self.coords = [Coordinate(v[0],v[1]) for v in coords]
        else: # # tuples?
            self.coords = [Coordinate(v[0],v[1]) for v in coords]

    def reversed(self):
        return PolyLine([v for v in reversed(self.coords)])

    def __getitem__(self, key):
        return self.coords[key]

    def __iter__(self):
        return iter(self.coords)

    def segments(self):
        return [(v1,v2) for v1, v2 in zip(self.coords[:-1],self.coords[1:])]

    def segment(self,i):
        if i>len(self.coords)-2:
            i = len(self.coords)-2
        return (self.coords[i],self.coords[i+1])

    def segmentLength(self,i):
        assert i<len(self.coords)+1
        return (self.coords[i+1]-self.coords[i]).length

    def unitVectors(self):
        return [(v2-v1)/(v2-v1).length for v1, v2 in zip(self.coords[:-1],self.coords[1:])]

    def clone(self):
        return PolyLine(self.coords)

    def createOffsetPolyLine(self, distance, resolution=3, join_style=JOIN_STYLE.round, mitre_limit=5):
        geosF = GeometryFactory()
        offset = geosF.createLineString(self.coords).parallel_offset(distance, resolution, join_style, mitre_limit)
        if distance < 0.:
            # pyGEOS reverses order for negative offset
            return PolyLine([v for v in reversed(offset.coords)])
        else:
            return PolyLine(offset.coords)

    def trimPoint(self,trimLength):
        seg = self.segment(int(trimLength))
        p = seg[0] + (seg[1]-seg[0]) * fmod(trimLength,1)
        return p

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
        c = None

        # intersections for the first few segments are most probable.
        # construct indices for combinations so that these are first tested.
        combinations = list(product(range(len(segs1)),range(len(segs2))))
        combinations = sorted(combinations,key=lambda x: sum(x))
        for i1,i2 in combinations:
            p1, p2 = segs1[i1]
            p3, p4 = segs2[i2]
            d1, d2 = p2-p1, p4-p3 
            denom = d2.y*d1.x - d2.x*d1.y
            if denom == 0.: continue # parallel
            d3 = p1 - p3
            ua = (d2.x*d3.y - d2.y*d3.x)/denom               
            if ua < 0. or ua > 1.: continue # out of range
            ub = (d1.x*d3.y - d1.y*d3.x)/denom
            if ub < 0. or ub > 1.: continue # out of range
            c = p1 + d1*ua
            p5 = p2 + d1*20
            p6 = p4 + d2*20
            # plt.plot([p1.x,p5.x],[p1.y,p5.y],'k:')
            # plt.plot([p3.x,p6.x],[p3.y,p6.y],'k:')
            break

        if c is None:
            # try to intersect the infinite lines of the first segments
                p1, p2 = segs1[0]
                p3, p4 = segs2[0]
                d1, d2 = p2-p1, p4-p3 
                denom = d2.y*d1.x - d2.x*d1.y
                if abs(denom) < 1.e-9:
                    return None,None,None,None,None
                d3 = p1 - p3
                ua = (d2.x*d3.y - d2.y*d3.x)/denom
                ub = (d1.x*d3.y - d1.y*d3.x)/denom
                c = p1 + d1*ua
                return c, i1, i2, d1/d1.length, d2/d2.length
        else:
            return c, i1+ua, i2+ub, d1/d1.length, d2/d2.length
