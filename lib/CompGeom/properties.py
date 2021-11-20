from mathutils import Vector

def orientation(a,b,c):
    # Return True if a,b,c are oriented clock-wise.
    return (b-a).cross(c-a) > 0.

def isClockwise(points):
    # get index of point with minimal x value
    iMin = min(range(len(points)), key=lambda i: points[i][0])
    # get previous, current and next points
    a = points[iMin-1]
    b = points[iMin]
    c = points[(iMin+1) % len(points)]
    return orientation(a,b,c)

def pointInPolygon(points, p):
    """
    Location of a point vis a vis a polygon.
    Returns either 2 if in, 1 if on, or 0 if out of polygom.
    Based on algorithm 7 from:
        Kai Horman and Alexander Agathos,
        "The point in polygon problem for arbitrary polygons".
        Computational Geometry: Theory and Applications, 
        Volume 20 Issue 3, November 2001
    See: https://www.sciencedirect.com/science/article/pii/S0925772101000128
    """
    # shift polygon so that p is (0,0)
    pi = [v-p for v in points]
    
    if pi[0] == Vector((0,0)): return 1 # on vertex

    wNr = 0
    for p1,p2 in zip(pi,pi[1:]+[pi[0]]):
        if p2[1] == 0.:
            if p2[0] == 0.:
                return 1 # on vertex
            else:
                if p1[1] == 0. and (p2[0] > 0.) == (p1[0] < 0.):
                    return 1 # on edge
        # if crossing horizontal line
        if (p1[1] < 0. and p2[1] >= 0.) or (p1[1] >= 0. and p2[1] < 0.):
            if p1[0] >= 0.:
                if p2[0] > 0.:
                    # modify w
                    wNr += 1 if p2[1] > p1[1] else -1
                else:
                    det = p1[0] * p2[1] - p2[0] * p1[1]
                    if det == 0: return 1 # on edge
                    # if right crossing
                    if (det > 0 and p2[1] > p1[1]) or (det < 0 and p2[1] < p1[1]):
                        # modify w
                        wNr += 1 if p2[1] > p1[1] else -1
            else:
                if p2[0] > 0.:
                    det = p1[0] * p2[1] - p2[0] * p1[1]
                    if det == 0: return 1 # on edge
                    # if right crossing
                    if (det > 0 and p2[1] > p1[1]) or (det < 0 and p2[1] < p1[1]):
                        # modify w
                        wNr += 1 if p2[1] > p1[1] else -1
    if (wNr % 2) != 0:
        return 2    # in polygon
    else:
        return 0    # out of polygon


    # if polygon.first.y == point.y and polygon.first.x == point.x:
    #     return "on" # vertex
    # w =0
    # for v in polygon.iter():
    #     if v.next.y == point.y:
    #         if v.next.x == point.x:
    #             return "on" # vertex
    #         else:
    #             if v.y == point.y and (v.next.x > point.x) == (v.x < point.x):
    #                 return "on" # edge
    #     # if crossing horizontal line
    #     if (v.y < point.y and v.next.y >= point.y)\
    #            or (v.y >= point.y and v.next.y < point.y):
    #         if v.x >= point.x:
    #             if v.next.x > point.x:
    #                 # modify w
    #                 if v.next.y > v.y: w += 1
    #                 else: w -= 1
    #             else:
    #                 det = (v.x - point.x) * (v.next.y - point.y) \
    #                     - (v.next.x - point.x) * (v.y - point.y)
    #                 if det == 0: return "on" # edge
    #                 # if right crossing
    #                 if (det > 0 and v.next.y > v.y)\
    #                    or (det < 0 and v.next.y < v.y):
    #                     # modify w
    #                     if v.next.y > v.y: w += 1
    #                     else: w -= 1
    #         else:
    #             if v.next.x > point.x:
    #                 det = (v.x - point.x) * (v.next.y - point.y) \
    #                     - (v.next.x - point.x) * (v.y - point.y)
    #                 if det == 0: return "on" # edge
    #                 # if right crossing
    #                 if (det > 0 and v.next.y > v.y)\
    #                    or (det < 0 and v.next.y < v.y):
    #                     # modify w
    #                     if v.next.y > v.y: w += 1
    #                     else: w -= 1
    # if (w % 2) != 0:
    #     return "in"
    # else:
    #     return "out"