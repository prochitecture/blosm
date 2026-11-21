"""Finite, exact-rational polynomial exclusion from original polygon interiors.

The stored binary64 controls and polygon vertices are exact input constants.
All subdivision, Bernstein signs, finite projections, and winding arithmetic
use Fraction; a positive result has no floating point rounding assumption.
"""
from fractions import Fraction


class NativeC2ExposureError(ValueError):
    """The declared finite proof could not establish continuous exposure."""


def certify(controls, areas, *, maximum_depth=32, maximum_cells=20000):
    if not 1 <= maximum_depth <= 32 or not 1 <= maximum_cells <= 20000:
        raise NativeC2ExposureError('C2 exposure needs its existing finite proof bounds')
    q=lambda value: Fraction.from_float(float(value))
    controls=tuple(tuple(map(q,row)) for row in controls)
    if len(controls)<3 or any(len(row)!=2 for row in controls):
        raise NativeC2ExposureError('C2 exposure needs finite planar polynomial controls')
    sub=lambda a,b:(a[0]-b[0],a[1]-b[1])
    cross=lambda a,b:a[0]*b[1]-a[1]*b[0]
    dot=lambda a,b:a[0]*b[0]+a[1]*b[1]
    def bbox(points):
        return tuple(min(p[j] for p in points) for j in (0,1)),tuple(max(p[j] for p in points) for j in (0,1))
    def split(values):
        work=list(values);left=[work[0]];right=[work[-1]]
        for _ in range(len(work)-1):
            work=[tuple((a+b)/2 for a,b in zip(x,y)) for x,y in zip(work,work[1:])]
            left.append(work[0]);right.append(work[-1])
        return tuple(left),tuple(reversed(right))
    edges=[];polygons=[]
    for area in areas:
        points=tuple(tuple(map(q,row)) for row in area.region.xy)
        polygons.append((area.region.owner.text(),points))
        for i,(a,b) in enumerate(zip(points,(*points[1:],points[0]))):
            low,high=bbox((a,b));edges.append((area.region.owner.text(),i,a,b,low,high))
    def exclude(values,edge,start,end):
        _,_,a,b,_,_=edge;vector=sub(b,a);length2=dot(vector,vector)
        projection=[dot(sub(p,a),vector) for p in values]
        if max(projection)<0 or min(projection)>length2:return 'exact_finite_segment_projection'
        coefficients=[cross(sub(p,a),vector) for p in values];factors=0
        # Only global endpoints can be declared native contact exceptions.
        while len(coefficients)>1 and start==0 and coefficients[0]==0:
            n=len(coefficients)-1
            coefficients=[Fraction(n,i)*coefficients[i] for i in range(1,n+1)];factors+=1
        while len(coefficients)>1 and end==1 and coefficients[-1]==0:
            n=len(coefficients)-1
            coefficients=[Fraction(n,n-i)*coefficients[i] for i in range(n)];factors+=1
        if min(coefficients)>0 or max(coefficients)<0:
            return 'exact_global_endpoint_factor_sign' if factors else 'exact_Bernstein_strict_sign'
        return None
    pending=[(controls,Fraction(0),Fraction(1),tuple(range(len(edges))),0)]
    cells=leaves=used_depth=0;methods={}
    while pending:
        values,start,end,indices,depth=pending.pop();cells+=1;used_depth=max(used_depth,depth)
        if cells>maximum_cells:raise NativeC2ExposureError('C2 exposure exhausts its finite proof cells')
        low,high=bbox(values);remaining=[]
        for i in indices:
            edge=edges[i]
            if any(high[j]<edge[4][j] or low[j]>edge[5][j] for j in (0,1)):
                method='exact_control_box'
            else:method=exclude(values,edge,start,end)
            if method:methods[method]=methods.get(method,0)+1
            else:remaining.append(i)
        if not remaining:leaves+=1;continue
        if depth>=maximum_depth:
            raise NativeC2ExposureError('C2 exposure has unresolved original polygon contacts')
        left,right=split(values);middle=(start+end)/2
        pending.extend(((left,start,middle,tuple(remaining),depth+1),(right,middle,end,tuple(remaining),depth+1)))
    middle=split(controls)[0][-1];witness={}
    for owner,polygon in polygons:
        crossings=0
        for a,b in zip(polygon,(*polygon[1:],polygon[0])):
            v=sub(b,a);w=sub(middle,a)
            if cross(w,v)==0 and 0<=dot(w,v)<=dot(v,v):
                raise NativeC2ExposureError('C2 midpoint remains on an original polygon boundary')
            if (a[1]>middle[1])!=(b[1]>middle[1]):
                left=(middle[0]-a[0])*v[1];right=(middle[1]-a[1])*v[0]
                if (v[1]>0 and left<right) or (v[1]<0 and left>right):crossings+=1
        if crossings%2:raise NativeC2ExposureError('C2 curve is inside an original polygon')
        witness[owner]='outside'
    return dict(policy='EXACT_BINARY64_INPUT_RATIONAL_BERNSTEIN_POLYGON_EXCLUSION',
        original_areas=len(polygons),original_edges=len(edges),proof_cells=cells,
        certified_leaves=leaves,maximum_used_depth=used_depth,exclusion_methods=methods,
        exact_midpoint_classifications=witness,continuous_outside_original_polygon_interiors=True,
        underlying_native_offset_separation_claimed=False)
