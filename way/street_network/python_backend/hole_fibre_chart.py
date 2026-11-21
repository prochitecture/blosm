"""Protected-point capped quad strips for unequal sampled separator banks."""
from types import SimpleNamespace
from dataclasses import dataclass
import math
import numpy as np
from .bundle_axis import Boundary, Island
from .fibre_strips import CappedFibreFrame, CappedFibreStrips
from .four_side_charts import quality
from .mesh import publication_coordinates
from .quad_charts import ChartResult
from .regions import signed_area


def _pair_rows(first, second, origin, axis, left, *, max_pairs=500000):
    xa, xb = (first-origin) @ axis, (second-origin) @ axis
    if np.any(np.diff(xa) <= 1e-8) or np.any(np.diff(xb) <= 1e-8):
        raise ValueError('Separator banks are not resolved monotone source curves')
    na, nb = len(first), len(second)
    if na*nb > max_pairs:
        raise ValueError('Separator source-row pairing exceeds its pair budget')
    ya, yb = (first-origin) @ left, (second-origin) @ left
    costs = np.full((na, nb), np.inf)
    previous = {}
    costs[0, 0] = 0.

    def update(i, j, ni, nj, extra):
        value = costs[i, j]+extra
        if value < costs[ni, nj]:
            costs[ni, nj] = value
            previous[ni, nj] = (i, j)

    for i in range(na):
        for j in range(nb):
            if not np.isfinite(costs[i, j]):
                continue
            if i+1 < na and j+1 < nb and ((i+1 == na-1) == (j+1 == nb-1)):
                at = .5*(xa[i+1]+xb[j+1])
                width = np.interp(at, xa, ya)-np.interp(at, xb, yb)
                update(i, j, i+1, j+1, ((xa[i+1]-xb[j+1])/max(width, 1e-4))**2)
            if i+1 < na-1 and j+1 < nb and xb[j]+1e-8 < xa[i+1] < xb[j+1]-1e-8:
                update(i, j, i+1, j, .1)
            if j+1 < nb-1 and i+1 < na and xa[i]+1e-8 < xb[j+1] < xa[i+1]-1e-8:
                update(i, j, i, j+1, .1)
    if not np.isfinite(costs[-1, -1]):
        raise ValueError('Separator banks have no bounded monotone row pairing')
    steps, at = [], (na-1, nb-1)
    while at != (0, 0):
        before = previous[at]
        steps.append((before, at))
        at = before
    a, b = [first[0]], [second[0]]
    for (i, j), (ni, nj) in steps[::-1]:
        if ni == i+1 and nj == j+1:
            a.append(first[ni]); b.append(second[nj])
        elif ni == i+1:
            a.append(first[ni])
            b.append(second[j]+(xa[ni]-xb[j])/(xb[j+1]-xb[j])*(second[j+1]-second[j]))
        else:
            b.append(second[nj])
            a.append(first[i]+(xb[nj]-xa[i])/(xa[i+1]-xa[i])*(first[i+1]-first[i]))
    return np.asarray(a), np.asarray(b), float(costs[-1, -1])


class _SkewFrame(CappedFibreFrame):
    def __init__(self, plan, island, number, spacing):
        self.plan, self.island, self.number, self.spacing = plan, island, number, spacing
        a, b = island.first.points-plan.origin, island.second.points-plan.origin
        xa, xb = a @ plan.axis, b @ plan.axis
        self.original = .5*(xa+xb)
        self.upper, self.lower = a @ plan.left, b @ plan.left
        if np.any(np.diff(self.original) <= 0.):
            raise ValueError('Separator logical rows are not strictly monotone')
        rows, split = [self.original[0]], {1, 2, len(a)-4, len(a)-3}
        for i, (x, y) in enumerate(zip(self.original, self.original[1:])):
            if i in split:
                rows.extend((x+(y-x)/3., x+2.*(y-x)/3.))
            rows.append(y)
        self.rows = np.asarray(rows)
        self.extra_start = max(xa[1], xb[1])-self.original[1]+max(
            np.linalg.norm(a[1]-a[0]), np.linalg.norm(b[1]-b[0]))/7.8
        self.extra_end = self.original[-2]-min(xa[-2], xb[-2])+max(
            np.linalg.norm(a[-2]-a[-1]), np.linalg.norm(b[-2]-b[-1]))/7.8
        self._update_columns()

    def _update_columns(self):
        delta = np.diff(self.rows)
        step = np.minimum(np.r_[delta[0], delta], np.r_[delta, delta[-1]])
        width = np.interp(self.rows, self.original, self.upper-self.lower)
        self.across = max(4, int(np.max(2.*np.ceil(width/(2.*np.minimum(self.spacing, 6.*step))))))
        if self.across > 128 or len(self.rows) > 8192:
            raise ValueError('Separator fibre chart exceeds its row/column budget')

    def refine(self):
        super().refine()
        self._update_columns()

    def columns(self, row):
        if row in (1, len(self.rows)-2):
            return 2
        if row in (0, len(self.rows)-1):
            return self.across
        step = min(self.rows[row]-self.rows[row-1], self.rows[row+1]-self.rows[row])
        width = float(np.interp(self.rows[row], self.original, self.upper-self.lower))
        return min(self.across, max(2, 2*int(math.ceil(width/(2.*min(self.spacing, 6.*step))))))

    def point(self, x, v):
        x, v = np.broadcast_arrays(np.asarray(x), np.asarray(v))
        source_x = x.copy()
        width = np.interp(x, self.original, self.upper-self.lower)
        fade = np.maximum(0., 1.-(width/(.6*self.spacing))**2)
        fade = np.where((x == self.rows[1]) | (x == self.rows[-2]), np.maximum(fade, .5), fade)
        length = self.rows[-1]-self.rows[0]
        delta = np.maximum(np.minimum((x-self.rows[0])/7.8, .28*self.spacing), self.extra_start)*(
            1.-self.smooth((x-self.rows[0]-.25*length)/(.25*length)))
        delta -= np.maximum(np.minimum((self.rows[-1]-x)/7.8, .28*self.spacing), self.extra_end)*(
            self.smooth((x-self.rows[0]-.5*length)/(.25*length)))
        x = x+delta*fade*(1.-abs(2.*v-1.))
        upper, lower = np.interp(x, self.original, self.upper), np.interp(x, self.original, self.lower)
        result = self.plan.world(x, (1.-v)*upper+v*lower)
        # Original world coordinates remain authoritative on both banks,
        # including unequal-X matched rows and the two collapsed tip points.
        for side, points in ((0., self.island.first.points), (1., self.island.second.points)):
            selected = v == side
            if np.any(selected):
                result[selected] = np.column_stack([np.interp(source_x[selected], self.original, points[:, k]) for k in (0, 1)])
        for row in (0, -1):
            selected = source_x == self.original[row]
            if np.any(selected):
                result[selected] = ((1.-v[selected, None])*self.island.first.points[row]
                                    +v[selected, None]*self.island.second.points[row])
        return result


def _existing_fibre_hole(patterns, boundary, spacing, max_cells):
    points = np.asarray(boundary, dtype=float)
    if (points.ndim != 2 or points.shape[1] != 2 or len(points) < 4
            or not np.isfinite(points).all() or signed_area(points) <= 0.):
        raise ValueError('Separator fibre chart requires a finite CCW physical rim')
    if not math.isfinite(spacing) or spacing <= 0. or type(max_cells) is not int or max_cells < 1:
        raise ValueError('Invalid separator fibre-chart spacing or cell budget')
    origin = points.mean(axis=0)
    _, vectors = np.linalg.eigh((points-origin).T @ (points-origin))
    axis = vectors[:, -1]
    left = np.array((-axis[1], axis[0]))
    along = (points-origin) @ axis
    a, b, n = int(np.argmin(along)), int(np.argmax(along)), len(points)
    first = points[[(a+j) % n for j in range((b-a) % n+1)]]
    second = points[[(a-j) % n for j in range((a-b) % n+1)]]
    if np.mean((first-origin) @ left) < np.mean((second-origin) @ left):
        first, second = second, first
    first, second, pairing_cost = _pair_rows(first, second, origin, axis, left)
    settings = SimpleNamespace(spacing=spacing, exaggeration=1.)
    plan = SimpleNamespace(origin=origin, axis=axis, left=left, datum=0., rise=0., settings=settings,
        islands=(Island(Boundary(first, []), Boundary(second, []), (True, True)),))
    plan.world = lambda x, y: origin+np.asarray(x)[..., None]*axis+np.asarray(y)[..., None]*left
    context = SimpleNamespace(patterns=patterns, network=SimpleNamespace(settings=settings))
    frame = _SkewFrame(plan, plan.islands[0], 0, spacing)
    builder = CappedFibreStrips(max_cells=max_cells)
    failures = []
    for attempt in range(3):
        try:
            mesh, tips = builder.assemble((frame,), context, 0.)
            report = quality(mesh)
            stored = quality(mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True)))
            if max(report['maximum_edge_ratio'], stored['maximum_edge_ratio']) > 8.:
                raise ValueError('Separator fibre chart exceeds its stored edge-ratio limit')
            edges = np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2)
            unique, inverse, counts = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True)
            if len(mesh.xyz)-len(unique)+len(mesh.quads) != 1 or not set(counts).issubset({1, 2}):
                raise ValueError('Separator fibres do not form a conforming disk')
            following = {int(a): int(b) for a, b in edges[counts[inverse] == 1]}
            starts = [i for i in following if np.array_equal(mesh.xyz[i, :2], points[0])]
            if len(starts) != 1:
                raise ValueError('Separator fibre chart lost an original physical tip')
            ids, current = [], starts[0]
            while current not in ids:
                ids.append(current)
                current = following[current]
            if current != ids[0] or len(ids) != len(following):
                raise ValueError('Separator fibre boundary is not one simple cycle')
            rim = mesh.xyz[ids, :2]
            positions = []
            for point in points:
                matches = np.flatnonzero(np.all(rim == point, axis=1))
                if len(matches) != 1:
                    raise ValueError('Separator fibre chart changed an original boundary vertex')
                positions.append(int(matches[0]))
            if positions != sorted(positions):
                raise ValueError('Separator fibre chart changed boundary order')
            for i, start in enumerate(positions):
                end = positions[(i+1) % n]+(len(rim) if i+1 == n else 0)
                p, q = points[i], points[(i+1) % n]
                d = q-p
                sample = rim[np.arange(start, end+1) % len(rim)]
                t = (sample-p) @ d/(d @ d)
                if (np.max(np.linalg.norm(sample-(p+t[:, None]*d), axis=1)) > 1e-8
                        or np.any(np.diff(t) <= 0.)):
                    raise ValueError('Separator fibre samples leave an original physical chord')
            if abs(signed_area(rim)-signed_area(points)) > max(1e-8, signed_area(points)*1e-10):
                raise ValueError('Separator fibre chart changed its physical area')
            mesh.seams = {0: np.array(ids, dtype=np.int32)}
            return ChartResult(mesh, dict(quality=report, stored_quality=stored, refinements=attempt,
                previous_failures=failures, pairing_cost=pairing_cost, protected_vertices=len(points),
                added_boundary_samples=len(rim)-len(points), exact_physical_polygon=True,
                construction='protected-point capped separator fibre strips'))
        except ValueError as exc:
            failures.append(str(exc))
            if attempt == 2:
                raise ValueError('Separator fibre chart remains unqualified: '+str(exc)) from exc
            frame.refine()


@dataclass(frozen=True)
class RuledHoleBanks:
    """Two finite monotone banks with actual adjacent edges at both tips."""
    boundary: object
    origin: object
    axis: object
    left: object
    first: object
    second: object
    pairing_cost: float
    intrinsic_tip_sine: float
    tips: tuple

    @classmethod
    def build(cls,boundary):
        points=np.asarray(boundary,float).copy()
        if (points.ndim!=2 or points.shape[1]!=2 or not 6<=len(points)<=20000
                or not np.isfinite(points).all() or signed_area(points)<=0.):
            raise ValueError('Ruled closed-hole banks need a bounded original CCW polygon')
        origin=points.mean(0)
        _,vectors=np.linalg.eigh((points-origin).T@(points-origin))
        axis=vectors[:,-1];left=np.array((-axis[1],axis[0]));along=(points-origin)@axis
        a,b,n=int(np.argmin(along)),int(np.argmax(along)),len(points)
        first=points[[(a+j)%n for j in range((b-a)%n+1)]]
        second=points[[(a-j)%n for j in range((a-b)%n+1)]]
        if min(len(first),len(second))<3:
            raise ValueError('Ruled closed-hole tips need two actual adjacent bank edges')
        if np.mean((first-origin)@left)<np.mean((second-origin)@left):first,second=second,first
        for bank in (first,second):
            if np.any(np.diff((bank-origin)@axis)<=1e-8):
                raise ValueError('Ruled closed-hole banks lack finite ordered source progress')
        tips=[];intrinsic=1.
        for end,row,neighbor in ((0,0,1),(1,-1,-2)):
            if not np.array_equal(first[row],second[row]):
                raise ValueError('Ruled closed-hole banks do not share their original physical tip')
            vectors=np.array((first[neighbor]-first[row],second[neighbor]-second[row]))
            lengths=np.linalg.norm(vectors,axis=1)
            sine=float(abs(np.linalg.det(vectors))/np.prod(lengths))
            if not math.isfinite(sine) or sine<=0.:
                raise ValueError('Ruled closed-hole tip has no positive original angular sector')
            intrinsic=min(intrinsic,sine)
            tips.append(dict(end=end,original_tip=first[row].tolist(),
                adjacent_original_bank_vertices=(first[neighbor].tolist(),second[neighbor].tolist()),
                intrinsic_sine=sine))
        # Preserve the real two adjacent tip vertices. The row-pairing DP may
        # interpolate only between the interior protected bank rows.
        aa,bb,cost=_pair_rows(first[1:-1],second[1:-1],origin,axis,left)
        aa=np.vstack((first[0],aa,first[-1]));bb=np.vstack((second[0],bb,second[-1]))
        if np.any((aa[1:-1]-bb[1:-1])@left<=0.):
            raise ValueError('Ruled closed-hole banks change their original lateral order')
        for value in (points,aa,bb):value.flags.writeable=False
        return cls(points,origin,axis,left,aa,bb,cost,intrinsic,tuple(tips))


class _RuledHoleFrame(_SkewFrame):
    def __init__(self,plan,island,number,spacing):
        super().__init__(plan,island,number,spacing)
        self.rows=self.original.copy()
        self.closed,self.index=island.collapsed,number
        self._update_columns()

    def columns(self,row):
        count=super().columns(row)
        # Each half then has an odd rim count, compatible with the existing
        # 1:3 templates. Keep the existing 128-column ceiling.
        return min(126,max(2,4*int(math.ceil((count-2)/4.))+2))

    def point(self,x,v):
        x,v=np.broadcast_arrays(np.asarray(x),np.asarray(v));source_x=x.copy()
        width=np.interp(x,self.original,self.upper-self.lower)
        fade=np.maximum(0.,1.-(width/(.6*self.spacing))**2)
        fade=np.where((x==self.rows[1])|(x==self.rows[-2]),np.maximum(fade,.5),fade)
        length=self.rows[-1]-self.rows[0]
        delta=np.maximum(np.minimum((x-self.rows[0])/7.8,.28*self.spacing),self.extra_start)*(
            1.-self.smooth((x-self.rows[0]-.25*length)/(.25*length)))
        delta-=np.maximum(np.minimum((self.rows[-1]-x)/7.8,.28*self.spacing),self.extra_end)*(
            self.smooth((x-self.rows[0]-.5*length)/(.25*length)))
        shifted=x+delta*fade*(1.-abs(2.*v-1.))
        if (np.any(source_x<self.original[0]) or np.any(source_x>self.original[-1])
                or np.any(shifted<self.original[0]) or np.any(shifted>self.original[-1])):
            raise ValueError('Ruled hole coordinates cannot extrapolate their finite protected banks')
        banks=[np.stack([np.interp(shifted,self.original,bank[:,k]) for k in (0,1)],axis=-1)
               for bank in (self.island.first.points,self.island.second.points)]
        result=(1.-v[...,None])*banks[0]+v[...,None]*banks[1]
        for row in (0,-1):
            selected=source_x==self.original[row]
            if np.any(selected):
                result[selected]=(1.-v[selected,None])*self.island.first.points[row]+v[selected,None]*self.island.second.points[row]
        return result


def _ruled_paired_hole(patterns,boundary,spacing,max_cells,previous):
    from .bundle_open_fibre_islands import _assemble
    from .hole_banded_chart import audit
    from .quality_correction import QuadShapeCorrection
    pair=RuledHoleBanks.build(boundary)
    settings=SimpleNamespace(spacing=spacing,exaggeration=1.)
    plan=SimpleNamespace(origin=pair.origin,axis=pair.axis,left=pair.left,
        islands=(Island(Boundary(pair.first,[]),Boundary(pair.second,[]),(True,True)),))
    context=SimpleNamespace(patterns=patterns,network=SimpleNamespace(settings=settings))
    frame=_RuledHoleFrame(plan,plan.islands[0],0,spacing)
    mesh,_,_=_assemble(frame,context,0.,max_cells)
    mesh,grading=QuadShapeCorrection(max_cells=max_cells,max_iterations=1000,
        minimum_sine=min(.03,.5*pair.intrinsic_tip_sine),
        target_sine=min(.05,.75*pair.intrinsic_tip_sine)).improve(mesh)
    edges=np.stack((mesh.quads,np.roll(mesh.quads,-1,axis=1)),axis=-1).reshape(-1,2)
    unique,inverse,counts=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    boundary_edges=edges[counts[inverse]==1]
    following={int(a):int(b) for a,b in boundary_edges}
    if len(following)!=len(boundary_edges) or set(following)!=set(following.values()):
        raise ValueError('Ruled closed-hole boundary is not one manifold cycle')
    starts=[i for i in following if np.array_equal(mesh.xyz[i,:2],pair.boundary[0])]
    if len(starts)!=1:raise ValueError('Ruled closed-hole chart lost its prescribed boundary start')
    ids=[];current=starts[0]
    while current not in ids:
        ids.append(current);current=following[current]
    if current!=ids[0] or len(ids)!=len(following):
        raise ValueError('Ruled closed-hole chart has more than one boundary cycle')
    mesh.seams={0:np.asarray(ids,np.int32)}
    report=audit(mesh,pair.boundary)
    return ChartResult(mesh,dict(report,pairing_cost=pair.pairing_cost,
        ruled_paired_hole=dict(original_tips=pair.tips,intrinsic_tip_sine=pair.intrinsic_tip_sine,
            continuous_cartesian_banks=True,protected_tip_adjacent_rows=True,half_strip_templates=True,
            maximum_rows=8192,maximum_columns=128,actual_rows=len(frame.rows),
            maximum_actual_columns=max(frame.columns(row) for row in range(len(frame.rows))),
            original_failure=str(previous),grading=grading,cell_budget=max_cells),
        construction='protected original closed-hole banks with continuous ruled half-strip quads'))


def build_fibre_hole(patterns,boundary,spacing,max_cells,*,ruled_paired=False):
    """Keep the established chart first; a typed caller may try ruled banks."""
    if type(ruled_paired) is not bool:
        raise ValueError('Ruled paired-hole recovery requires an explicit Boolean policy')
    if not ruled_paired:return _existing_fibre_hole(patterns,boundary,spacing,max_cells)
    try:return _existing_fibre_hole(patterns,boundary,spacing,max_cells)
    except ValueError as previous:
        try:return _ruled_paired_hole(patterns,boundary,spacing,max_cells,previous)
        except ValueError as failure:
            raise ValueError('Original and ruled paired-hole charts remain unqualified: '+str(failure)) from previous
