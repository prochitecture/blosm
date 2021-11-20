from mathutils import Vector
from collections import defaultdict

# import matplotlib.pyplot as plt

from .properties import *
from .line import *

EPSILON = 0.0001

class Polygon():
    # 2D polygon. A list of points of type mathutils.Vector that form a closed polygon

    def __init__(self):
        self.points = []

    def __len__(self):
        return len(self.points)

    def __iter__(self):
        return iter(self.points)

    def __getitem__(self, key):
        return self.points[key]

    def __setitem__(self, key, value):
        self.points[key] = value

    def clear(self):
        self.points = []

    def flip(self):
		# reverses the orientation of the polygon
        self.points.reverse()
        return self

    def clone(self):
		# return a shallow copy of the polygon (points are not cloned)
        poly = Polygon()
        poly.points = [ p for p in self.points ]
        return poly

    @staticmethod
    def fromPointlist(points):
        # points: List of Vectors that make up the polygon
        p = Polygon()
        p.points = [p.freeze() for p in points]
        return p

    @staticmethod
    def fromTuples(tuples):
    # tuples: List of 2-tuples that make up the polygon
        p = Polygon()
        p.points = [ Vector((t[0], t[1])).freeze() for t in tuples ]
        return p

    def addPoint(self, point):
        self.points.append(point.freeze())

    def addPoints(self, points):
        self.points.extend([p.freeze() for p in points])

    def isClockwise(self):
        return Polygon.isClockwiseS(self.points)

    @staticmethod
    def isClockwiseS(points):
        vectors = [(p2-p1) for p1,p2 in zip(points,points[1:]+[points[0]])]
        area = 0.5 * sum(v1.cross(v2) for v1,v2 in zip(vectors[:-1],vectors[1:]) )
        return area < 0.

        # original code did not work correctly:

        # # get index of point with minimal x value
        # iMin = min(range(len(points)), key=lambda i: points[i][0])
        # # get previous, current and next points
        # a = points[iMin-1]
        # b = points[iMin]
        # c = points[(iMin+1) % len(points)]
        # return orientation(a,b,c)

    def isClockwise(self):
        return Polygon.isClockwiseS(self.points)

    @staticmethod
    def containsPointS(points, p) :
        # Checks if the polygon defined by the point list <points> contains the point <p>

		# see if we find a line segment that p is on
        for a,b in list(zip(points[0:], points[1:])) + [(points[-1], points[0])]:
            d = distPointSegment(p, a, b)
            if d < EPSILON: return 2

        # p is not on the boundary, cast ray and intersect to see if we are inside
        intersections = set(intersectPolygonRay(points, p, p + Vector((1,0))))


        # filter intersection points that are boundary points
        for isectPoint in set(filter(lambda x: x in points, intersections)):
            i = points.index(isectPoint)
            prv = points[i-1]
            nxt = points[(i+1) % len(points)]

            if orientation(p, isectPoint, nxt) == orientation(p, isectPoint, prv):
                intersections.remove(isectPoint)

        # we are inside if we have an odd amount of polygon intersections
        return 1 if len(intersections) % 2 == 1 else 0

    def containsPoint(self, p):
            return Polygon.containsPointS(self.points, p)

    def booleanOperation(polygonA, polygonB, operation):
		# for union and intersection, we want the same orientation on both polygons. 
        # for difference, we want different orientation.
        matchingOrientation = polygonA.isClockwise() == polygonB.isClockwise()
        if matchingOrientation != (operation != 'd'):
            polygonB = polygonB.clone()
            polygonB.flip()

        def inorderExtend(v, v1, v2, ints):
            # Extend a sequence v by points ints that are on the segment v1, v2
            k, r = None, False
            if v1.x < v2.x:
                k, r = lambda i: i.x, True
            if v1.x > v2.x:
                k, r = lambda i: i.x, False
            if v1.y < v2.y:
                k, r = lambda i: i.y, True
            else:
                k, r = lambda i: i.y, False
            l = [ (p, 2) for p in sorted(ints, key=k, reverse=r) ]
            i = next((i for i, p in enumerate(v) if p[0] == v2), -1)
            assert(i>=0)
            for e in l:
                v.insert(i, e)

        # initialize vector rings
        vA = [(p, polygonB.containsPoint(p)) for p in polygonA.points]
        vB = [(p, polygonA.containsPoint(p)) for p in polygonB.points]

        if sum(v[1] for v in vA) + sum(v[1] for v in vB) == 0:
           # no intersections at all, no subtraction
        #    print('no intersection')
           return [polygonA]

        # find all intersections
        intersectionsA = defaultdict(list)
        intersectionsB = defaultdict(list)
        for a1, a2 in list(zip(vA, vA[1:])) + [(vA[-1], vA[0])]:
            for b1, b2 in list(zip(vB, vB[1:])) + [(vB[-1], vB[0])]:
                i = intersectSegmentSegment(a1[0],a2[0],b1[0],b2[0])
                if i:
                    intersectionsA[(a1[0],a2[0])].append(i)
                    intersectionsB[(b1[0],b2[0])].append(i)

        # extend vector rings by intersections
        for k, v in intersectionsA.items():
            inorderExtend(vA, k[0], k[1], v)

        for k, v in intersectionsB.items():
            inorderExtend(vB, k[0], k[1], v)

        edgeFragments = defaultdict(list)
        def extendFragments(v, poly, fragmentType):
            for v1, v2 in list(zip(v, v[1:])) + [(v[-1], v[0])]:
                if v1[1] == fragmentType or v2[1] == fragmentType:
                    # one of the vertices is of the required type
                    edgeFragments[v1[0]].append( v2[0] )

                elif v1[1] == 2 and v2[1] == 2:
                    # we have two boundary vertices
                    m = (v1[0] + v2[0]) / 2.0
                    t = poly.containsPoint(m)
                    if t == fragmentType or t == 2:
                        edgeFragments[v1[0]].append( v2[0] )

        fragmentTypeA = 1 if operation == 'i' else 0
        fragmentTypeB = 1 if operation != 'u' else 0

        extendFragments(vA, polygonB, fragmentTypeA)
        extendFragments(vB, polygonA, fragmentTypeB)

        output = []
        while edgeFragments:
            start = list(edgeFragments.keys())[0]
            current = edgeFragments[start][0]
            sequence = [start]

            # follow along the edge fragments sequence
            while not current in sequence:
                sequence.append(current)
                current = edgeFragments[current][0]
 

            # get only the cyclic part of the sequence
            sequence = sequence[sequence.index(current):]

            for c,n in list(zip(sequence, sequence[1:])) + [(sequence[-1], sequence[0])]:
                edgeFragments[c].remove(n)

                if not edgeFragments[c]:
                    del edgeFragments[c]
            
            output.append(Polygon.fromPointlist(Polygon.simplifySequence(sequence)))

        return output

    @staticmethod
    def simplifySequence(seq):
		# simplify a point sequence so that no subsequent points are on the same line
        i = 0
        while i < len(seq):
            p, c, n = seq[i-1], seq[i], seq[(i + 1) % len(seq)]
            if p == c or c == n or p == n or distPointSegment(c, p, n) < EPSILON:
                del seq[i]
            else:
                i+=1
        return seq

    @staticmethod
    def union(polygon_a, polygon_b):
        return Polygon.booleanOperation(polygon_a, polygon_b, 'u')

    @staticmethod
    def intersect(polygon_a, polygon_b):
        return Polygon.booleanOperation(polygon_a, polygon_b, 'i')

    @staticmethod
    def subtract(polygon_a, polygon_b):
        return Polygon.booleanOperation(polygon_a, polygon_b, 'd')

