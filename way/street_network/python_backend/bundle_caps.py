"""Finite, profiled quad strips joining skew mouths to a common-fibre body.

Original paths, mouth profiles and rounded junction caps remain authoritative.
The internal cut planes are computational boundaries, not edits to source data.
"""
from dataclasses import replace
import numpy as np
from . import model
from .bundle_axis import AxialPlan, AxialChartUnavailable, _regular_rows
from .bundle_ends import AttachedAxialPlan
from .mesh import Assembler, Mesh
from .providers import Part
from .routes import Piece, Route
from .surface_query import QuadHeightField
from .entrance_slopes import EntranceSlopeGuard

_ENTRANCE_SLOPE = EntranceSlopeGuard()


def interp(curve, x):
    x = np.asarray(x)
    y = np.interp(x, curve[:, 0], curve[:, 1])
    lo = (curve[1, 1]-curve[0, 1])/(curve[1, 0]-curve[0, 0])
    hi = (curve[-1, 1]-curve[-2, 1])/(curve[-1, 0]-curve[-2, 0])
    y = np.where(x < curve[0, 0], curve[0, 1] + lo*(x-curve[0, 0]), y)
    y = np.where(x > curve[-1, 0], curve[-1, 1] + hi*(x-curve[-1, 0]), y)
    at = np.clip(np.searchsorted(curve[:, 0], x, side='right')-1, 0, len(curve)-2)
    slope = (curve[at+1, 1]-curve[at, 1])/(curve[at+1, 0]-curve[at, 0])
    return y, slope


def clip(curve, first, last):
    x = np.r_[first, curve[(curve[:, 0] > first+1e-9) & (curve[:, 0] < last-1e-9), 0], last]
    return np.column_stack((x, interp(curve, x)[0]))


class FiniteAxialPlan(AttachedAxialPlan):
    """Select a bounded finite end chart, retaining every prescribed rim."""
    margin_factors = (2., 2.5, 3., 1.5, 1., .75, .5, 4.)

    @classmethod
    def build(cls, context, bundle, **chart):
        from .bundle_compounds import compound_owners
        if compound_owners(context, bundle):
            # A compound may consume or crop an original mouth. Its visible
            # region owner must supply the cap; a full-profile finite strip
            # cannot restore the consumed port or overlap that pavement.
            raise AxialChartUnavailable('Compound-owned Bundle mouths retain their visible-region construction')
        failures = []
        for margin in cls.margin_factors:
            try:
                plan = cls(context, bundle, margin_factor=margin, **chart)
                ends=np.asarray((plan.first,plan.last))[list(plan.strip_ends())]
                if plan.physical_roles_swapped or np.min(plan.gap(ends)) <= 1e-7:
                    raise ValueError('Finite caps need positively separated, consistently ordered body ends')
                if plan.minimum_gap > 1e-8 and not plan.fill_micro_gap:
                    straight = []
                    for points in plan.visible_points:
                        delta, offsets = points[-1]-points[0], points-points[0]
                        error = abs(offsets[:, 0]*delta[1]-offsets[:, 1]*delta[0])/np.linalg.norm(delta)
                        straight.append(bool(np.max(error) <= 1e-8))
                    if all(straight):
                        return None
                plan._prefer_fibre_strips(context)
                part = plan.pavement(context)
                from .validation import validate_overlaps
                validate_overlaps({bundle.key: part.mesh}, cell_size=max(1., context.network.settings.spacing*2.))
                plan.cap_selection = dict(margin_factor=margin, prior_failures=failures)
                part.metadata['finite_cap_selection'] = plan.cap_selection
                return plan
            except ValueError as exc:
                failures.append(dict(margin_factor=margin, reason=str(exc)))
        raise AxialChartUnavailable('No qualified finite cap chart: ' + '; '.join(
            f"{r['margin_factor']:g}: {r['reason']}" for r in failures))

    def __init__(self, context, bundle, *, margin_factor):
        self.cap_meshes, self.cap_queries, self.cap_rims = {}, [], {}
        self.cap_keys = []
        self.cap_entrance_bounds = {}
        self.combined = None
        self.margin_factor, self.cap_selection = margin_factor, None
        super().__init__(context, bundle)
        for role in (0, 1):
            for end in self.strip_ends():
                self.make_strip(context, role, end)
        # Construction-only paths are not retained behind surface queries.
        del self.routes

    def prepare_curves(self, context):
        self.original_curves = tuple(p.copy() for p in self.curves)
        ends=self.strip_ends()
        self.cut_first = max(p[0, 0] for p in self.curves) + (self.margin_factor*self.settings.spacing if 0 in ends else 0.)
        self.cut_last = min(p[-1, 0] for p in self.curves) - (self.margin_factor*self.settings.spacing if 1 in ends else 0.)
        if self.cut_last <= self.cut_first + 2*self.settings.spacing:
            raise AxialChartUnavailable('Finite caps leave no positive common body')
        self.curves = [clip(p, self.cut_first, self.cut_last) for p in self.curves]
        self.routes = []
        for member in self.members:
            pieces = [Piece(context[model.Key('path', t.section_id)], t.forward) for t in member.path]
            self.routes.append(Route(pieces))

    def strip_ends(self):
        return (0,1)

    def make_strip(self, context, role, end):
        """Keep the source profile on a bounded quad strip with fixed rims.

        Axial coordinates interpolate the two prescribed side partitions.
        A cubic Hermite correction in the transverse coordinate matches the
        actual mouth and its tangent, fading to the common-fibre chart. It is
        zero on both lateral curbs. Profile height depends only on transverse
        q, so no longitudinal height overshoot is introduced by the blend.
        """
        route = self.routes[role]
        traversal = self.members[role].path[0 if end == 0 else -1]
        source = context.network.section_index[traversal.section_id]
        profile = context[model.Key('profile', source.id)].oriented(0 if traversal.forward else 1)
        width, offset = source.width, source.offset * (1 if traversal.forward else -1)
        station = 0. if end == 0 else route.length
        plane = self.first if end == 0 else self.last
        curves = self.original_curves[2*role:2*role+2]
        sides = []
        for curve in curves:
            local = clip(curve, curve[0, 0], plane) if end == 0 else clip(curve, plane, curve[-1, 0])[::-1]
            sides.append(_regular_rows(self.world(local[:, 0], local[:, 1]),
                                      min(self.settings.spacing, 2*width/self.port_intervals)))
        pattern = context.patterns.quad(len(sides[0])-1, self.port_intervals,
                                        len(sides[1])-1, self.port_intervals)
        if len(pattern.quads) > 4096:
            raise ValueError('Finite entrance strip exceeds its cell budget')
        u, v = pattern.uv.T
        ids0, ids1 = (context.patterns.boundary(pattern, side) for side in (0, 2))
        x0 = np.interp(u, pattern.uv[ids0, 0], (sides[0]-self.origin)@self.axis)
        x1 = np.interp(u, pattern.uv[ids1, 0], (sides[1]-self.origin)@self.axis)
        x = (1-v)*x0 + v*x1
        q = 2*v-1
        p, tangent, normal = route.frame(station)
        front = p + (offset+.5*width*q)[:, None]*normal
        front_x, front_y = (front-self.origin)@self.axis, (front-self.origin)@self.left
        y0, _ = interp(curves[0], x)
        y1, _ = interp(curves[1], x)
        front0, fs0 = interp(curves[0], front_x)
        front1, fs1 = interp(curves[1], front_x)
        shift = front_y - ((1-v)*front0+v*front1)
        step = min(1e-4, .01*route.length)
        a, b = max(0., station-step), min(route.length, station+step)
        pa, _, na = route.frame(a); pb, _, nb = route.frame(b)
        jet = (pb + (offset+.5*width*q)[:, None]*nb - pa - (offset+.5*width*q)[:, None]*na)/(b-a)
        source_slope = (jet@self.left)/(jet@self.axis)
        derivative = (source_slope-((1-v)*fs0+v*fs1))*(plane-front_x)
        shift[(v == 0.) | (v == 1.)] = 0.
        derivative[(v == 0.) | (v == 1.)] = 0.
        t = np.clip((x-front_x)/(plane-front_x), 0., 1.)
        y = (1-v)*y0+v*y1 + (1-3*t*t+2*t*t*t)*shift + (t-2*t*t+t*t*t)*derivative
        physical_heights = profile.value(q)
        xyz = np.column_stack((self.world(x, y), profile.datum+physical_heights*self.settings.exaggeration))
        for side, points in ((0, sides[0]), (2, sides[1])):
            ids = context.patterns.boundary(pattern, side)
            xyz[ids, :2] = points
        front_ids = context.patterns.boundary(pattern, 3)
        fv = pattern.uv[front_ids, 1]
        xyz[front_ids, :2] = p + (offset+.5*width*(2*fv-1))[:, None]*normal
        back_ids = context.patterns.boundary(pattern, 1)
        bv = pattern.uv[back_ids, 1]
        back_y = (1-bv)*interp(curves[0], plane)[0]+bv*interp(curves[1], plane)[0]
        xyz[back_ids, :2] = self.world(np.full(len(bv), plane), back_y)
        mesh = Mesh(xyz, pattern.quads[:, ::-1] if end == 0 else pattern.quads)
        quality = mesh.validate()
        if quality['maximum_edge_ratio'] > 8.:
            raise ValueError(f'Finite cap {role}/{end} exceeds ratio 8: {quality}')
        incoming = 2*float(np.max(abs(profile.derivative(np.array((-1.,1.))))))/width
        if mesh.slope_statistics(self.settings.exaggeration)['maximum'] > max(.1,1.1*incoming)+1e-8:
            raise ValueError('Finite cap exceeds physical slope limit')
        try:
            entrance_bound = _ENTRANCE_SLOPE.validate(mesh, front_ids, jet, physical_heights)
        except ValueError as exc:
            raise ValueError(f'Finite cap {role}/{end}: {exc}') from exc
        self.cap_entrance_bounds[role, end] = entrance_bound
        self.cap_meshes[role, end] = mesh
        self.cap_keys.append((role,end))
        self.cap_queries.append(QuadHeightField(mesh, datum=self.datum,
            ceiling=self.datum+self.rise*self.settings.exaggeration))
        rim = xyz[context.patterns.boundary(pattern, 2 if role == 0 else 0)].copy()
        self.cap_rims[role, end] = rim if end == 0 else rim[::-1]

    def height(self, xy):
        result = np.asarray(super().height(xy))
        rows = np.atleast_2d(xy)
        output = result.reshape(-1).copy()
        sentinel = self.datum-1.
        for query in self.cap_queries:
            lo, hi = query.origin+query.minimum.min(axis=0), query.origin+query.maximum.max(axis=0)
            ids = np.flatnonzero(np.all((rows >= lo) & (rows <= hi), axis=1))
            if len(ids):
                values = query.evaluate(rows[ids], default=sentinel)
                valid = values != sentinel
                output[ids[valid]] = values[valid]
        return output.reshape(result.shape)

    def keys(self, mesh):
        points = (mesh.xyz[:, :2]-self.origin)
        x, y = points@self.axis, points@self.left
        keys = [None]*len(points)
        for role in (0, 1):
            for end in self.strip_ends():
                plane=(self.first,self.last)[end]
                upper, lower = self.edges(plane)[2*role:2*role+2]
                t = (y-upper)/(lower-upper)
                ids = np.flatnonzero((abs(x-plane)<1e-8) & (t>=-1e-8) & (t<=1+1e-8))
                for i in ids:
                    k = int(round(float(t[i])*self.port_intervals))
                    if abs(t[i]-k/self.port_intervals) > 1e-8: raise ValueError('Nonmatching cap partition')
                    keys[i] = ('finite cap/body', role, end, k)
        return keys

    def pavement(self, context):
        if self.combined is None:
            body = super().pavement(context)
            if not self.cap_meshes:
                for key,query in zip(self.cap_keys,self.cap_queries):
                    xyz = np.column_stack((query.xy+query.origin,query.heights))
                    self.cap_meshes[key] = Mesh(xyz,query.quads.copy())
            assembly = Assembler()
            assembly.add(body.mesh.xyz, body.mesh.quads, self.keys(body.mesh))
            for cap in self.cap_meshes.values(): assembly.add(cap.xyz, cap.quads, self.keys(cap))
            mesh = assembly.mesh()
            if len(mesh.quads) > 50000:
                raise ValueError('Finite Bundle pavement exceeds its cell budget')
            quality = mesh.validate()
            if quality['maximum_edge_ratio'] > 8.:
                raise ValueError('Finite Bundle pavement exceeds its unchanged quad-quality limit')
            incoming = max(2*float(np.max(abs(context[model.Key('profile',sid)].derivative(np.array((-1.,1.)))))) /
                context.network.section_index[sid].width for sid in self.source_ids)
            slopes = mesh.slope_statistics(self.settings.exaggeration)
            if slopes['maximum'] > max(.1,1.1*incoming)+1e-8:
                raise ValueError('Finite Bundle pavement exceeds its physical slope limit')
            self.combined = Part(self.bundle.key, mesh, dict(body.metadata, quality=mesh.validate(),
                finite_caps=len(self.cap_meshes), cut_planes=[self.first,self.last],
                finite_cap_selection=self.cap_selection,
                maximum_entrance_longitudinal_slope_bound=max(self.cap_entrance_bounds.values()),
                entrance_longitudinal_slope_limit=_ENTRANCE_SLOPE.limit,
                cap_entrance_longitudinal_slope_bounds={f'{role}/{end}': bound
                    for (role,end),bound in self.cap_entrance_bounds.items()},
                physical_slopes=slopes, measured_slope_limit=max(.1,1.1*incoming),
                cap_quads=sum(len(c.quads) for c in self.cap_meshes.values())))
        return self.combined

    def separator(self, context):
        from .bundles import BundleProvider, Member, Cap
        body = AxialPlan.separator(self, context)
        members = tuple(Member.build(context, member, role) for role, member in enumerate(self.visible_members))
        caps, joins, attached = [], [None,None], [False,False]
        for end in self.strip_ends():
            cap = BundleProvider.cap(context, *members, end)
            attached[end]=cap is not None
            island = self.islands[0 if end == 0 else -1]
            across = max(4, 4*int(np.ceil(np.linalg.norm(island.first.points-island.second.points,axis=1).max()/(4*self.settings.spacing))))
            pieces = tuple(replace(m, rim=self.cap_rims[role, end]) for role,m in enumerate(members))
            flat = Cap.free(pieces[0].rim[-1 if end == 0 else 0],pieces[1].rim[-1 if end == 0 else 0],across)
            joins[end]=(flat.rim[0], flat.rim[-1], across)
            caps.append(BundleProvider.region(context,self.bundle,pieces,[cap,flat] if end == 0 else [flat,cap],self.proof_checks))
        assembly = Assembler()
        for part in (body,*caps): assembly.add(part.mesh.xyz,part.mesh.quads,self.seam_keys(part.mesh,joins))
        mesh = assembly.mesh('grass')
        edges = np.sort(np.stack((mesh.quads,np.roll(mesh.quads,-1,axis=1)),axis=-1).reshape(-1,2),axis=1)
        unique, counts = np.unique(edges,axis=0,return_counts=True)
        if counts.max() > 2 or len(mesh.xyz)-len(unique)+len(mesh.quads) != len(self.islands):
            raise ValueError('Finite separator caps do not form their declared islands')
        quality = mesh.validate()
        if quality['maximum_edge_ratio'] > 8.:
            raise ValueError('Finite Bundle separator exceeds its unchanged quad-quality limit')
        return Part(self.bundle.key,mesh,dict(quality=quality,islands=len(self.islands),finite_caps=True,
            attached_caps=attached, finite_cap_selection=self.cap_selection,
            physical_slopes=mesh.slope_statistics(self.settings.exaggeration)))

    def release_geometry(self):
        super().release_geometry()
        self.combined = None
        self.cap_meshes.clear()


def finite_plan(context, bundle):
    """Keep the established branch first, then try the merged free-end chart."""
    try:
        return FiniteAxialPlan.build(context,bundle)
    except AxialChartUnavailable as original:
        from .bundle_free_caps import MergedFreeCapPlan
        try:
            alternative=MergedFreeCapPlan.build(context,bundle)
            if alternative is not None:
                return alternative
        except AxialChartUnavailable:
            pass
        from .bundle_terminal_caps import SkewMergedCapPlan
        try:
            return SkewMergedCapPlan.build(context,bundle)
        except AxialChartUnavailable as exc:
            raise AxialChartUnavailable(str(original)+'; '+str(exc)) from exc

