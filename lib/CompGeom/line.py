from mathutils import Vector

def intersectLineLine(p1, p2, q1, q2):
    d = (q2[1] - q1[1]) * (p2[0] - p1[0]) - (q2[0] - q1[0]) * (p2[1] - p1[1])
    if d == 0: return None
    n1 = (q2[0] - q1[0]) * (p1[1] - q1[1]) - (q2[1] - q1[1]) * (p1[0] - q1[0])
    n2 = (p2[0] - p1[0]) * (p1[1] - q1[1]) - (p2[1] - p1[1]) * (p1[0] - q1[0])
    u1 = n1 / d
    u2 = n2 / d
    return Vector((u1, u2)).freeze()

def intersectSegmentSegment(p1, p2, q1, q2):
    # <p1,p2> and <q1,q2> are the segments
    if max(q1[0], q2[0]) < min(p1[0], p2[0]): return None
    if min(q1[0], q2[0]) > max(p1[0], p2[0]): return None
    if max(q1[1], q2[1]) < min(p1[1], p2[1]): return None
    if min(q1[1], q2[1]) > max(p1[1], p2[1]): return None

    ll = intersectLineLine(p1, p2, q1, q2)

    if ll == None: return None
    if ll[0] < 0 or ll[0] > 1: return None
    if ll[1] < 0 or ll[1] > 1: return None

    return Vector( (p1[0] + ll[0] * (p2[0] - p1[0]) , p1[1] + ll[0] * (p2[1] - p1[1]) ) ).freeze()

def intersectSegmentRay(p1, p2, q1, q2):
    # <p1,p2> is segment, <q1,q2> are points on ray
    l = intersectLineLine(p1, p2, q1, q2)
    if l == None: return None
    if l[0] < 0 or l[0] > 1: return None
    if l[1] < 0: return None
    return Vector( ((p1[0] + l[0] * (p2[0] - p1[0]) , p1[1] + l[0] * (p2[1] - p1[1]))) ).freeze()

def intersectSegmentsRay(segs, q1, q2):
    # <segs> is a list of point pairs of segments
    # <q1,q2> are points on ray
    intersectPoints = []
    for seg in segs:
        iSect = intersectSegmentRay(seg[0], seg[1], q1, q2)
        if iSect:
            intersectPoints += [iSect]
    return intersectPoints

def intersectPolygonRay(polyPoints, p1, p2):
    # <polyPoints> is list of points in the polygon
    # <p1,p2> are points on ray
    return intersectSegmentsRay(list(zip(polyPoints[0:], polyPoints[1:])) + [(polyPoints[-1], polyPoints[0])], p1, p2)

def distPointSegment(p, p1, p2):
    # shortest distance from a point <p> to a line segment <p1,p2>
    pVec = (p2-p1)
    L = pVec.length
    t = ((p[0] - p1[0]) * (p2[0] - p1[0]) + (p[1] - p1[1]) * (p2[1] - p1[1])) / L
    t = max(0., min(1., t))
    intSect = p1 + pVec*t
    return (intSect-p).length

