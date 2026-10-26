"""Common-fibre Bundle pavement union, dynamic crown and separator islands.

The physical offset boundaries are retained as bounded-error polylines. A
monotone axial chart finds their events and outer envelope without polygon
triangulation or a numerical surface solve. Curved member interiors are allowed;
this first chart branch requires compatible, aligned free end caps. Junction
caps continue through the separated boundary provider until compound charts
are available. No eligible geometry is silently trimmed to a common rectangle.
"""
from dataclasses import dataclass
import math
import numpy as np
from .mesh import Assembler, Mesh
from .paths import unit
from .providers import Part
from .model import End, Key


class AxialChartUnavailable(ValueError):
    """Try another chart without treating unsupported geometry as valid."""


class SkewCaps(AxialChartUnavailable):
    """The monotone body needs finite entrance strips at its physical ends."""


def _simplify(points, tolerance, protected=()):
    """Polyline error bound in arbitrary coordinate dimension, bounded stack."""
    if len(points) <= 2:
        return points
    keep = {0, len(points) - 1, *protected}
    ordered = sorted(keep)
    pending = list(zip(ordered, ordered[1:]))
    while pending:
        a, b = pending.pop()
        if b <= a + 1:
            continue
        delta = points[b] - points[a]
        t = np.clip((points[a + 1:b] - points[a]) @ delta / max(float(delta @ delta), 1e-30), 0., 1.)
        error = np.linalg.norm(points[a + 1:b] - points[a] - t[:, None] * delta, axis=1)
        k = int(np.argmax(error)) + a + 1
        if error[k - a - 1] > tolerance:
            keep.add(k)
            pending.extend(((a, k), (k, b)))
    return points[sorted(keep)]


def _regular_rows(points, spacing, multiple=2):
    pieces = []
    for a, b in zip(points, points[1:]):
        count = max(1, int(math.ceil(np.linalg.norm(b - a) / spacing)))
        if len(pieces) + count + 1 > 200_000:
            raise ValueError('Bundle boundary exceeds its allocation budget; increase cell size')
        pieces.extend(np.linspace(a, b, count + 1)[:-1])
    pieces.append(points[-1])
    result = np.asarray(pieces)
    while (len(result) - 1) % multiple:
        at = int(np.argmax(np.linalg.norm(np.diff(result, axis=0), axis=1)))
        result = np.insert(result, at + 1, .5 * (result[at] + result[at + 1]), axis=0)
    return result


@dataclass
class Boundary:
    points: np.ndarray
    keys: list


@dataclass
class Island:
    first: Boundary
    second: Boundary
    collapsed: tuple


@dataclass
class SeparatedPlan:
    mode: str = 'separated'
    reason: str | None = None


class AxialPlan:
    mode = 'union'
    coordinated_fibres = False

    def __init__(self, context, bundle, *, allow_separated=False):
        self.bundle, self.settings = bundle, context.network.settings
        self.members = bundle.members
        self._separator_cache = None
        self.separator_correction = None
        self.origin = context[Key('path', self.members[0].path[0].section_id)].points[
            0 if self.members[0].path[0].forward else -1].copy()
        starts, directions = [], []
        for member in self.members:
            a = member.path[0]
            path = context[Key('path', a.section_id)]
            station = 0. if a.forward else path.length
            point, tangent, _ = path.frame(station)
            starts.append(point)
            directions.append(tangent * (1. if a.forward else -1.))
        self.axis = self.chart_axis(context,directions)
        self.left = np.array((-self.axis[1], self.axis[0]))
        if float((starts[0] - starts[1]) @ self.left) <= 1e-5:
            raise ValueError('Bundle members disagree with their declared left-to-right order')
        self.curves = []
        self.profiles = []
        self.source_ids = tuple(t.section_id for m in self.members for t in m.path)
        self.smooth = all(not context[Key('path', sid)].exact for sid in self.source_ids)
        from .bundle_boundaries import MemberBoundary
        owners, exact_corners = set(), []
        for role, member in enumerate(self.members):
            boundary = MemberBoundary(context, member, self.axis)
            owners.update(boundary.owners)
            for side in (-1, 1):
                xy = boundary.points(side) - self.origin
                curve = np.column_stack((xy @ self.axis, xy @ self.left))
                if np.any(np.diff(curve[:, 0]) <= 1e-9):
                    raise AxialChartUnavailable('Bundle physical boundary is not monotone in this axial chart')
                self.curves.append(curve)
            exact_corners.extend(boundary.exact_corners)
            self.profiles.extend(context[Key('profile', t.section_id)] for t in member.path)
        self.transition_owners = tuple(sorted(owners))
        fixed = np.sort((np.asarray(exact_corners)-self.origin) @ self.axis) if exact_corners else np.empty(0)
        # Mirrored curb corners can have the same axial coordinate. Rigid
        # transforms introduce round-off; do not turn one computational cut
        # into a microscopic extra interval. Physical curbs stay untouched.
        self.fixed_peak_rows = fixed[np.r_[True,np.diff(fixed)>1e-8]] if len(fixed) else fixed
        # Stored order is A outer, A inner, B inner, B outer.
        self.prepare_curves(context)
        self.first, self.last = float(self.curves[0][0, 0]), float(self.curves[0][-1, 0])
        if any(abs(p[0, 0] - self.first) > 1e-6 or abs(p[-1, 0] - self.last) > 1e-6 for p in self.curves):
            raise SkewCaps('Axial Bundle free caps are skew: a finite cap chart is required')
        self.datum, self.rise = self.profiles[0].datum, self.profiles[0].rise
        if any(abs(p.datum - self.datum) > 1e-7 or abs(p.rise - self.rise) > 1e-7 for p in self.profiles):
            raise ValueError('Axial Bundle members require one common datum and rise')
        base = np.unique(np.concatenate([p[:, 0] for p in self.curves]))
        values = self.edges(base)
        if np.min(values[0] - values[1]) <= 1e-8 or np.min(values[2] - values[3]) <= 1e-8:
            raise AxialChartUnavailable('A member has no positive physical width in this axial chart')
        events = list(base)
        for i in range(4):
            for j in range(i):
                difference = values[i] - values[j]
                hit = difference[:-1] * difference[1:] < 0.
                for k in np.flatnonzero(hit):
                    events.append(float(base[k] + (base[k + 1] - base[k]) * difference[k] / (difference[k] - difference[k + 1])))
        events = np.asarray(sorted(events))
        self.nodes = events[np.r_[True, np.diff(events) > 1e-9]]
        event_values = self.edges(self.nodes)
        gap = np.maximum(event_values[1], event_values[3]) - np.minimum(event_values[0], event_values[2])
        self.minimum_gap, self.maximum_gap = float(gap.min()), float(gap.max())
        self.physical_roles_swapped = bool(np.any(event_values[0] + event_values[1] < event_values[2] + event_values[3]))
        self.fill_micro_gap = self.minimum_gap >= 0. and self.maximum_gap <= min(.001, self.settings.curve_tolerance * .1)
        self.islands, self.peak_a, self.peak_b = [], [], []
        if (self.minimum_gap > 1e-8 and not self.fill_micro_gap and not self.physical_roles_swapped
                and not allow_separated):
            return
        self.regions = self._regions()
        self._make_peaks()
        self.outer_a = self._outer(upper=True)
        self.outer_b = self._outer(upper=False)
        if self.islands:
            try:
                preview = self.separator_chart(context)
                quality = preview.validate()
            except ValueError:
                preview = None
                quality = dict(maximum_edge_ratio=math.inf)
            if quality['maximum_edge_ratio'] > 8.:
                from .fibre_strips import CappedFibreStrips
                replacement = CappedFibreStrips().build(context, self)
                for name, extras in (('peak_a', replacement.extra_first), ('peak_b', replacement.extra_second)):
                    boundary = getattr(self, name)
                    points = np.vstack((boundary.points, [p for p, key in extras])) if extras else boundary.points
                    keys = boundary.keys + [key for p, key in extras]
                    order = np.argsort((points - self.origin) @ self.axis, kind='stable')
                    setattr(self, name, Boundary(points[order], [keys[i] for i in order]))
                self.separator_correction = replacement.metadata
                preview = replacement.mesh
            self._separator_cache = preview

    def chart_axis(self, context, directions):
        return unit(np.asarray(directions).sum(axis=0))

    def prepare_curves(self, context):
        """A specialized end chart may restrict these original boundaries."""

    def edges(self, x):
        return np.asarray([np.interp(x, p[:, 0], p[:, 1]) for p in self.curves])

    def peaks(self, x):
        ao, ai, bi, bo = self.edges(x)
        high, low = np.maximum(ao, bi), np.minimum(ai, bo)
        inner_high, inner_low = np.maximum(ai, bo), np.minimum(ao, bi)
        gap = inner_high - inner_low
        if self.fill_micro_gap:
            gap = np.zeros_like(gap)
        overlap = np.maximum(0., -gap)
        fraction = np.clip(overlap / np.maximum(1e-10, np.minimum(ao - ai, bi - bo)), 0., 1.)
        blend = fraction ** 3 * (10. + fraction * (-15. + 6. * fraction))
        peak = (1. - blend) * .5 * (inner_high + inner_low) + blend * .5 * (high + low)
        return np.where(gap > 0., inner_high, peak), np.where(gap > 0., inner_low, peak)

    def gap(self, x):
        ao, ai, bi, bo = self.edges(x)
        return np.maximum(ai, bo) - np.minimum(ao, bi)

    def height(self, xy):
        local = np.asarray(xy) - self.origin
        x, y = local @ self.axis, local @ self.left
        ao, ai, bi, bo = self.edges(x)
        high, low = np.maximum(ao, bi), np.minimum(ai, bo)
        a, b = self.peaks(x)
        distance = np.where(y >= a, (y - a) / np.maximum(high - a, 1e-12),
                    np.where(y <= b, (y - b) / np.maximum(b - low, 1e-12), 0.))
        return self.datum + self.rise * np.clip(1. - distance * distance, 0., 1.) * self.settings.exaggeration

    def world(self, x, y):
        return self.origin + np.asarray(x)[..., None] * self.axis + np.asarray(y)[..., None] * self.left

    def _regions(self):
        if self.fill_micro_gap:
            return [(self.first, self.last, False)]
        # Only actual inner-edge contacts split gap topology. Other boundary
        # crossings are retained for outer-envelope ownership, not new islands.
        gap = self.gap(self.nodes)
        roots = [self.first, self.last]
        roots.extend(float(self.nodes[i]) for i in range(1, len(self.nodes) - 1) if abs(gap[i]) < 1e-8)
        roots = np.unique(roots)
        regions = []
        for a, b in zip(roots, roots[1:]):
            mid = .5 * (a + b)
            positive = bool(self.gap(mid) > 1e-8)
            if regions and not positive and not regions[-1][2]:
                regions[-1] = (regions[-1][0], b, False)
            else:
                regions.append((a, b, positive))
        return regions

    def _make_peaks(self):
        for region, (a, b, positive) in enumerate(self.regions):
            xs = np.r_[a, self.nodes[(self.nodes > a + 1e-9) & (self.nodes < b - 1e-9)], b]
            refined = [a]
            for start, end in zip(xs, xs[1:]):
                pending = [(start, end, 0)]
                while pending:
                    lo, hi, depth = pending.pop()
                    t = np.linspace(lo, hi, 5)
                    ya, yb = self.peaks(t)
                    p = np.column_stack((t, ya, yb))
                    line = p[0] + np.linspace(0., 1., 5)[:, None] * (p[-1] - p[0])
                    if np.max(np.linalg.norm(p - line, axis=1)) > self.settings.curve_tolerance * .2:
                        if depth >= 14:
                            raise ValueError('Bundle peak approximation exceeded its refinement budget')
                        middle = .5 * (lo + hi)
                        pending.extend(((middle, hi, depth + 1), (lo, middle, depth + 1)))
                    else:
                        refined.append(hi)
            t = np.asarray(refined)
            ya, yb = self.peaks(t)
            reference = np.column_stack((t, ya, yb))
            p = _simplify(reference, 1e-9 if positive else self.settings.curve_tolerance * .2)
            multiple = 8 if positive else 2
            p = _regular_rows(p, self.settings.spacing, multiple=multiple)
            if self.transition_owners or self.coordinated_fibres:
                sampled = self._peak_interval(a, b, t, multiple, limit=512)
                if sampled is None:
                    if self.coordinated_fibres:
                        query = self._partition_peak_rows(a,b,t)
                    else:
                        # Retain the established partition for dimensional
                        # charts; their fixed-corner parity is independently
                        # qualified, including rotated Exact-mode junctions.
                        cuts=np.r_[a,self.fixed_peak_rows[(self.fixed_peak_rows>a+1e-8)
                                                        &(self.fixed_peak_rows<b-1e-8)],b]
                        parts=[self._peak_interval(lo,hi,t,multiple,limit=4096)
                               for lo,hi in zip(cuts,cuts[1:])]
                        if any(p is None for p in parts):
                            raise ValueError('Bundle peak rows exceed the approximation budget')
                        query=np.r_[np.concatenate([p[0][:-1] for p in parts]),b]
                        if len(query)>65536:
                            raise ValueError('Bundle peak boundary exceeds its allocation budget')
                    qa,qb = self.peaks(query)
                else:
                    query,qa,qb = sampled
                p = np.column_stack((query, qa, qb))
            elif self.smooth:
                # Smooth lookup knots are not mandatory mesh rows. Nearly
                # coincident knots from opposite offsets distort a collapsed
                # quad atlas. Uniform axial rows retain the physical curves to
                # a separately checked error, rather than copying the lookup.
                count = max(multiple, multiple * int(math.ceil((b - a) / (multiple * self.settings.spacing))))
                while True:
                    query = np.linspace(a, b, count + 1)
                    qa, qb = self.peaks(query)
                    error = max(float(np.max(abs(np.interp(t, query, q) - y)))
                                for q, y in ((qa, ya), (qb, yb)))
                    if error <= self.settings.curve_tolerance * .8:
                        break
                    count *= 2
                    if count > 65536:
                        raise ValueError('Bundle peak rows exceed the approximation budget')
                p = np.column_stack((query, qa, qb))
            x = p[:, 0]
            ya, yb = self.peaks(x)
            boundary = []
            for role, y in enumerate((ya, yb)):
                xy = self.world(x, y)
                keys = []
                for row, xi in enumerate(x):
                    if row == 0 or row == len(x) - 1:
                        endpoint = region if row == 0 else region + 1
                        shared = not positive or abs(self.gap(xi)) < 1e-8
                        keys.append(('axial event', endpoint, None if shared else role))
                    else:
                        keys.append(('axial peak', region, row, role if positive else None))
                boundary.append(Boundary(xy, keys))
            if positive:
                collapsed = (np.linalg.norm(boundary[0].points[0] - boundary[1].points[0]) < 1e-7,
                             np.linalg.norm(boundary[0].points[-1] - boundary[1].points[-1]) < 1e-7)
                self.islands.append(Island(*boundary, collapsed))
            self.peak_a.append(boundary[0])
            self.peak_b.append(boundary[1])
            if sum(len(p.points) for p in self.peak_a) > 200_000:
                raise ValueError('Bundle peak boundary exceeds its allocation budget; increase cell size')
        self.peak_a = self._combine(self.peak_a)
        self.peak_b = self._combine(self.peak_b)

    def _peak_interval(self, a, b, reference, multiple, *, limit):
        """Bounded rows against authoritative lookup data; no curve is moved."""
        t = np.r_[a, reference[(reference > a) & (reference < b)], b]
        ya,yb = self.peaks(t)
        edges = self.edges(t)
        outer = np.maximum(edges[0],edges[2]), np.minimum(edges[1],edges[3])
        count = max(multiple, multiple * int(math.ceil((b-a)/(multiple*self.settings.spacing))))
        if self.coordinated_fibres:
            # Prefer splitting a sharp corner to making every wide separator
            # row much shorter than its transverse cell width.
            limit=min(limit,2*count)
        while count <= limit:
            query = np.linspace(a, b, count+1)
            qa,qb = self.peaks(query)
            error = max(float(np.max(abs(np.interp(t,query,q)-y))) for q,y in ((qa,ya),(qb,yb)))
            qe = self.edges(query)
            error = max(error, *(float(np.max(abs(np.interp(t,query,q)-y))) for q,y in
                zip((np.maximum(qe[0],qe[2]),np.minimum(qe[1],qe[3])),outer)))
            if error <= self.settings.curve_tolerance * .8:
                return query,qa,qb
            count *= 2
        return None

    def _partition_peak_rows(self, first, last, reference):
        """Split only where a bounded uniform approximation actually fails.

        Forcing every outer/inner Exact kink into the common row partition can
        introduce microscopic intervals for a harmless small-angle bend. The
        worst measured lookup error instead selects necessary cuts; every leaf
        still passes the same exterior/crest approximation check.
        """
        pending=[(first,last,0)]
        parts=[];checks=0;rows=0
        while pending:
            a,b,depth=pending.pop();checks+=1
            if checks>1024:
                raise ValueError('Bundle peak row partition exceeded its search budget')
            sampled=self._peak_interval(a,b,reference,1,limit=512)
            if sampled is not None:
                parts.append(sampled[0][:-1]);rows+=len(sampled[0])-1
                if rows>65536:
                    raise ValueError('Bundle peak boundary exceeds its allocation budget')
                continue
            t=np.r_[a,reference[(reference>a)&(reference<b)],b]
            count=(min(512,2*max(1,int(math.ceil((b-a)/self.settings.spacing))))
                   if self.coordinated_fibres else 512)
            query=np.linspace(a,b,count+1)
            errors=np.zeros(len(t))
            for x,put in ((t,False),(query,True)):
                peaks=self.peaks(x);edges=self.edges(x)
                values=(*peaks,np.maximum(edges[0],edges[2]),np.minimum(edges[1],edges[3]))
                if not put:
                    expected=values
                else:
                    for value,actual in zip(values,expected):
                        errors=np.maximum(errors,abs(np.interp(t,query,value)-actual))
            cut=float(t[int(np.argmax(errors))])
            if float(errors.max()) <= self.settings.curve_tolerance*.8:
                # The 512-row budget can also be exceeded by a long straight
                # interval whose shape is already exact. Split its length.
                cut=.5*(a+b)
            if depth>=16 or min(cut-a,b-cut)<=1e-8:
                raise ValueError('Bundle peak rows exceed the approximation budget')
            pending.extend(((cut,b,depth+1),(a,cut,depth+1)))
        return np.r_[np.concatenate(parts),last]

    @staticmethod
    def _combine(parts):
        return Boundary(np.concatenate([p.points[:-1] for p in parts] + [parts[-1].points[-1:]]),
                        sum((p.keys[:-1] for p in parts), []) + parts[-1].keys[-1:])

    def _outer(self, upper):
        v = self.edges(self.nodes)
        y = np.maximum(v[0], v[2]) if upper else np.minimum(v[1], v[3])
        points = _simplify(np.column_stack((self.nodes, y)), 1e-9)
        p = _regular_rows(points, self.settings.spacing)
        return self.world(p[:, 0], p[:, 1])

    def pavement_chart(self, context):
        """First direct chart and semantic crown guides, before grading."""
        if self.transition_owners or self.coordinated_fibres:
            return self.transition_pavement_chart(context)
        assembly = Assembler()
        for role, outer, peak in ((0, self.outer_a, self.peak_a), (1, self.outer_b, self.peak_b)):
            # Even longitudinal boundaries allow local direct-quad refinement.
            points, keys = peak.points, list(peak.keys)
            if (len(points) - 1) % 2:
                raise ValueError('Bundle peak partition lost its even boundary parity')
            # Use actual paired physical half-widths, not distance to chart origin.
            x = (points - self.origin) @ self.axis
            edges = self.edges(x)
            y = (points - self.origin) @ self.left
            width = np.max((np.maximum(edges[0], edges[2]) - y) if role == 0 else (y - np.minimum(edges[1], edges[3])))
            across = max(2, 2 * int(math.ceil(float(width) / (2. * self.settings.spacing))))
            across = getattr(self, 'port_intervals', across)
            pattern = context.patterns.quad(len(outer) - 1, across, len(points) - 1, across)
            u, v = pattern.uv.T
            outer_ids, peak_ids = (context.patterns.boundary(pattern, side) for side in (0, 2))
            U0, U1 = pattern.uv[outer_ids, 0], pattern.uv[peak_ids, 0]
            p0 = np.column_stack([np.interp(u, U0, outer[:, c]) for c in range(2)])
            p1 = np.column_stack([np.interp(u, U1, points[:, c]) for c in range(2)])
            # Remap in the physical fibre chart. Interpolating the XY points
            # directly can cut outside a curved boundary when opposite rim
            # partitions have different axial coordinates at the same atlas u.
            x0, x1 = (p0 - self.origin) @ self.axis, (p1 - self.origin) @ self.axis
            axial = (1. - v) * x0 + v * x1
            outer_local, peak_local = outer - self.origin, points - self.origin
            y0 = np.interp(axial, outer_local @ self.axis, outer_local @ self.left)
            y1 = np.interp(axial, peak_local @ self.axis, peak_local @ self.left)
            xy = self.world(axial, (1. - v) * y0 + v * y1)
            xyz = np.column_stack((xy, self.height(xy)))
            xyz[peak_ids, :2] = points
            xyz[peak_ids, 2] = self.datum + self.rise * self.settings.exaggeration
            xyz[outer_ids, :2] = outer
            xyz[outer_ids, 2] = self.datum
            addresses = [None] * len(xyz)
            for i, identity in zip(peak_ids, keys):
                addresses[int(i)] = identity
            faces = pattern.quads[:, ::-1] if role == 0 else pattern.quads
            assembly.add(xyz, faces, addresses)
        mesh = assembly.mesh()
        guides = np.array(sorted(set(assembly.lookup.values())), dtype=np.int32)
        return mesh, guides

    def transition_pavement_chart(self, context):
        """Coordinated axial rows keep both sides of every quad positive.

        The selected rows bound both exterior-curb and crest approximation.
        Separator refinement supplies any extra crest rows before this build.
        Every cross-row edge lies on one physical fibre, so a sharp Exact
        corner cannot invert a cell through unequal outer/crest atlas indices.
        """
        assembly=Assembler()
        guide_keys=set()
        for role,peak in enumerate((self.peak_a,self.peak_b)):
            x=(peak.points-self.origin)@self.axis
            y=(peak.points-self.origin)@self.left
            edges=self.edges(x)
            outer=np.maximum(edges[0],edges[2]) if role==0 else np.minimum(edges[1],edges[3])
            width=float(np.max(abs(outer-y)))
            across=getattr(self,'port_intervals',max(2,2*int(math.ceil(width/(2*self.settings.spacing)))))
            if 2*(len(x)-1)*across > 50000:
                counts=np.maximum(2,2*np.ceil(abs(outer-y)/(2*self.settings.spacing))).astype(np.int32)
                if hasattr(self,'port_intervals'):
                    counts[[0,-1]]=self.port_intervals
                for row in range(len(x)-1):
                    pattern=context.patterns.strip(1,int(counts[row]),int(counts[row+1]))
                    u,v=pattern.uv.T
                    axial=(1-u)*x[row]+u*x[row+1]
                    high=(1-u)*outer[row]+u*outer[row+1]
                    low=(1-u)*y[row]+u*y[row+1]
                    xy=self.world(axial,(1-v)*high+v*low)
                    xyz=np.column_stack((xy,self.height(xy)))
                    xyz[v==0.,2]=self.datum
                    xyz[v==1.,2]=self.datum+self.rise*self.settings.exaggeration
                    keys=[]
                    for i,(a,b) in enumerate(pattern.uv):
                        if a not in (0.,1.):
                            keys.append(None)
                            continue
                        at=row+int(a)
                        if b==1.:
                            keys.append(peak.keys[at])
                            xyz[i,:2]=peak.points[at]
                        else:
                            keys.append(('dimensional fibre',role,at,round(float(b),12)))
                    assembly.add(xyz,pattern.quads[:,::-1] if role==0 else pattern.quads,keys)
                guide_keys.update(peak.keys)
                continue
            pattern=context.patterns.grid(len(x)-1,across)
            v=np.linspace(0.,1.,across+1)
            xy=self.world(x[:,None],(1-v)*outer[:,None]+v*y[:,None])
            heights=self.height(xy.reshape(-1,2)).reshape(xy.shape[:2])
            xyz=np.dstack((xy,heights))
            xyz[:,-1,:2]=peak.points
            xyz[:,-1,2]=self.datum+self.rise*self.settings.exaggeration
            xyz[:,0,2]=self.datum
            keys=[None]*(xyz.shape[0]*xyz.shape[1])
            for row,key in enumerate(peak.keys):
                keys[row*(across+1)+across]=key
            faces=pattern.quads[:,::-1] if role==0 else pattern.quads
            assembly.add(xyz.reshape(-1,3),faces,keys)
            guide_keys.update(peak.keys)
        mesh=assembly.mesh()
        if len(mesh.quads)>50000:
            raise ValueError('Dimensional Bundle pavement exceeds its cell budget')
        return mesh,np.array(sorted({assembly.lookup[key] for key in guide_keys}),dtype=np.int32)

    def pavement(self, context):
        mesh, guides = self.pavement_chart(context)
        quality = mesh.validate()
        correction = None
        if quality['maximum_edge_ratio'] > 8.:
            from .surface_quality import ProfiledQuadCorrection
            incoming = max(2. * float(np.max(abs(context[Key('profile', sid)].derivative(np.array((-1., 1.)))))) /
                           context.network.section_index[sid].width for sid in self.source_ids)
            mesh, correction = ProfiledQuadCorrection().improve(mesh, self.height, fixed_vertices=guides,
                datum=self.datum, ceiling=self.datum + self.rise * self.settings.exaggeration,
                exaggeration=self.settings.exaggeration, maximum_slope=max(.10, 1.1 * incoming))
            quality = mesh.validate()
        return Part(self.bundle.key, mesh, dict(quality=quality, construction='common-fibre pavement union',
            dimensional_transition_owners=[key.text() for key in self.transition_owners],
            quad_correction=correction,
            profile='dynamic half crowns; no plateau in overlap', islands=len(self.islands),
            minimum_gap=self.minimum_gap, maximum_gap=self.maximum_gap,
            submillimetre_gap_filled=self.fill_micro_gap,
            physical_roles_swapped=self.physical_roles_swapped,
            physical_slopes=mesh.slope_statistics(self.settings.exaggeration)))

    def separator_chart(self, context):
        assembly = Assembler()
        for island in self.islands:
            n = len(island.first.points) - 1
            width = float(np.linalg.norm(island.first.points - island.second.points, axis=1).max())
            across = max(4, 4 * int(math.ceil(width / (4. * self.settings.spacing))))
            first, last = (0 if collapsed else across for collapsed in island.collapsed)
            pattern = context.patterns.band(n, first, last, middle=across)
            u, v = pattern.uv.T
            row = np.linspace(0., 1., n + 1)
            a = np.column_stack([np.interp(u, row, island.first.points[:, c]) for c in range(2)])
            b = np.column_stack([np.interp(u, row, island.second.points[:, c]) for c in range(2)])
            xy = (1. - v[:, None]) * a + v[:, None] * b
            xyz = np.column_stack((xy, np.full(len(xy), self.datum + self.rise * self.settings.exaggeration)))
            assembly.add(xyz, pattern.quads[:, ::-1])
        mesh = assembly.mesh('grass')
        return mesh

    def separator(self, context):
        if not self.islands:
            return None
        if self._separator_cache is None:
            if self.separator_correction is None:
                self._separator_cache = self.separator_chart(context)
            else:
                from .fibre_strips import CappedFibreStrips
                self._separator_cache = CappedFibreStrips().build(context, self).mesh
        mesh = self._separator_cache
        return Part(self.bundle.key, mesh, dict(quality=mesh.validate(), islands=len(self.islands),
                    fibre_correction=self.separator_correction,
                    construction='common-fibre direct quad separator islands'))

    def release_geometry(self):
        self._separator_cache = None


def plan(context, bundle):
    if len(bundle.members) != 2 or not all(m.complete for m in bundle.members):
        return SeparatedPlan()
    if bundle.interior_section_ids:
        from .bundle_merged_entrances import preparation_key
        merged = context.products.get(preparation_key(bundle))
        if merged is not None:
            from .bundle_merged_plan import MergedBundlePlan
            return MergedBundlePlan(context, bundle, merged)
        from .bundle_angular_internal import InteriorAngularPlan
        try:
            return InteriorAngularPlan.build(context, bundle) or SeparatedPlan()
        except AxialChartUnavailable as error:
            return SeparatedPlan(reason=str(error))
    external = set()
    for member in bundle.members:
        a, b = member.path[0], member.path[-1]
        external.update((End(a.section_id, 0 if a.forward else 1), End(b.section_id, 1 if b.forward else 0)))
    if any(e in external for j in (*context.network.intersections, *context.network.connections) for e in j.ends):
        # The separated provider certifies the actual visible ranges after
        # junction trims. Authored terminal centerlines may meet at a declared
        # merge junction, which belongs to that junction rather than the Bundle.
        from .bundle_ends import AttachedAxialPlan
        try:
            candidate = AttachedAxialPlan.build(context, bundle)
        except SkewCaps:
            from .bundle_caps import finite_plan
            try:
                candidate = finite_plan(context, bundle)
            except AxialChartUnavailable as exc:
                return attached_angular_fallback(context, bundle, str(exc))
        except AxialChartUnavailable as exc:
            return attached_angular_fallback(context, bundle, str(exc))
        return candidate if candidate is not None else SeparatedPlan()
    from .bundles import certify_disjoint
    certify_disjoint(tuple(context[Key('path', t.section_id)] for t in bundle.members[0].path),
                     tuple(context[Key('path', t.section_id)] for t in bundle.members[1].path))
    # Prove the authored centerline constraint even when union removes an
    # otherwise invalid overlap between member pavement surfaces.
    try:
        candidate = AxialPlan(context, bundle)
    except SkewCaps:
        from .bundle_caps import finite_plan
        try:
            candidate = finite_plan(context, bundle)
        except AxialChartUnavailable as exc:
            return angular_fallback(context, bundle, str(exc))
        if candidate is None:
            return SeparatedPlan()
    except AxialChartUnavailable as exc:
        return angular_fallback(context, bundle, str(exc))
    return candidate if candidate.minimum_gap <= 1e-8 or candidate.fill_micro_gap or candidate.physical_roles_swapped else SeparatedPlan()


def attached_angular_fallback(context, bundle, reason):
    from .bundle_angular_attached import AttachedAngularPlan
    try:
        return AttachedAngularPlan.build(context, bundle)
    except AxialChartUnavailable as exc:
        return SeparatedPlan(reason=f'{reason}; {exc}')


def angular_fallback(context, bundle, reason):
    """Try a certified curvilinear union after established free-end branches.

    Eligibility failures retain the separated provider and its own geometry
    checks. A failed eligible angular surface is an explicit construction error,
    never permission to suppress the member pavements without a replacement.
    """
    from .bundle_angular import AngularBundlePlan
    try:
        return AngularBundlePlan(context, bundle)
    except AxialChartUnavailable as exc:
        return SeparatedPlan(reason=f'{reason}; {exc}')
