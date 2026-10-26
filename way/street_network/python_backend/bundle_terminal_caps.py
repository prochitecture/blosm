"""Finite union caps for nonparallel free ends of a monotone Bundle.

The shared body ends at an internal cross-section. Actual member curbs and
terminal cuts bound the remaining union; no source end is rotated or shortened.
Direct four-sided quad patterns retain that union and the body's whole edges.
"""
from dataclasses import replace
import numpy as np
from .bundle_axis import AxialPlan, AxialChartUnavailable, _simplify
from .bundle_caps import FiniteAxialPlan, clip, interp
from .bundle_ends import AttachedAxialPlan, member_ends
from .entrance_slopes import EntranceSlopeGuard
from .four_side_charts import FourSidedQuadChart
from .mesh import Mesh
from .model import Key
from .paths import unit
from .quality_correction import QuadShapeCorrection
from .regions import Region, RegionUnion
from .routes import Piece, Route
from .surface_query import QuadHeightField


class ExtendedCrown:
    """Continue the body's dynamic crown over its retained physical end cuts."""
    def __init__(self, plan, *, physical=False):
        self.origin, self.axis, self.left = plan.origin, plan.axis, plan.left
        self.curves = plan.original_curves
        self.datum, self.rise = plan.datum, plan.rise
        self.fill_micro_gap = plan.fill_micro_gap
        self.settings = replace(plan.settings, exaggeration=1.) if physical else plan.settings

    def edges(self, x):
        return np.asarray([interp(p, x)[0] for p in self.curves])

    peaks = AxialPlan.peaks
    __call__ = AxialPlan.height


def sample_chain(points, count, tolerance):
    """Uniform arc rows with a checked bound against the supplied polyline."""
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    stations = np.r_[0., np.cumsum(lengths)]
    t = np.linspace(0., stations[-1], count + 1)
    result = np.column_stack([np.interp(t, stations, points[:, i]) for i in range(2)])
    result[[0, -1]] = points[[0, -1]]
    # Every source segment between retained samples lies within its endpoint
    # error from the corresponding new chord. Include both station partitions
    # so a new chord never conceals an original polyline corner.
    query = np.unique(np.r_[t, stations])
    original = np.column_stack([np.interp(query, stations, points[:, i]) for i in range(2)])
    actual = np.column_stack([np.interp(query, t, result[:, i]) for i in range(2)])
    if np.max(np.linalg.norm(actual-original, axis=1)) > tolerance:
        raise ValueError('Terminal curb rows exceed the remaining curve-accuracy budget')
    return result


def graded_segment(length, spacing, first, last, extra=0):
    """Grow boundary cells locally, without copying a tiny step everywhere.

    The first/last cell has the supplied upper length; interior lengths grow
    by at most two per cell. Normalizing the weights preserves the exact line
    segment. Only the number of direct boundary intervals is adaptive.
    """
    if not np.all(np.isfinite((length, spacing, first, last))) or min(length, spacing, first, last) <= 0.:
        raise ValueError('Graded terminal segment needs positive finite lengths')
    count = max(1, int(np.ceil(length/spacing)))
    while count <= 8192:
        i = np.arange(count, dtype=float)
        weights = np.minimum(1., np.minimum(first/spacing*2.**np.minimum(i, 32.),
                                           last/spacing*2.**np.minimum(count-1-i, 32.)))
        cells = length*weights/weights.sum()
        if cells.max() <= spacing*(1.+1e-12):
            break
        count += 1
    if count > 8192 or not 1 <= count+extra <= 8192:
        raise ValueError('Graded terminal segment exceeds its allocation budget')
    count += extra
    i = np.arange(count, dtype=float)
    weights = np.minimum(1., np.minimum(first/spacing*2.**np.minimum(i, 32.),
                                       last/spacing*2.**np.minimum(count-1-i, 32.)))
    result = np.r_[0., np.cumsum(weights)/weights.sum()]
    result[-1] = 1.
    return result


def terminal_samples(points, spacing, body_count, bias, *, graded=False):
    """Keep all cut/curb corners; choose a compatible boundary partition.

    A compatible partition of the long terminal cuts can move an interior
    spoke away from a reflex corner. Only sample counts change, never the
    physical boundary. Every original cut/curb corner stays present.
    """
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    counts = np.maximum(1, np.ceil(lengths/spacing).astype(np.int64))
    extras = np.zeros(len(lengths), dtype=np.int64)
    for at, extra in zip(np.argsort(-lengths, kind='stable')[:2], bias):
        extras[at] = extra
    fractions = []
    for i, length in enumerate(lengths):
        if graded:
            first = min(spacing, 16.*lengths[i-1]) if i else spacing
            last = min(spacing, 16.*lengths[i+1]) if i+1<len(lengths) else spacing
            t = graded_segment(length, spacing, first, last, int(extras[i]))
            counts[i] = len(t)-1
        else:
            counts[i] += extras[i]
            if not 1 <= counts[i] <= 8192:
                raise ValueError('Terminal cut samples exceed their allocation budget')
            t = np.linspace(0., 1., int(counts[i])+1)
        fractions.append(t)
    if counts.min() < 1 or counts.sum() > 8192:
        raise ValueError('Terminal cut samples exceed their allocation budget')
    parts = [((a+(b-a)*t[:,None]) if graded else np.linspace(a,b,int(n)+1))[:-1]
             for a,b,n,t in zip(points, points[1:], counts, fractions)]
    result = np.vstack((*parts, points[-1]))
    if (len(result)-body_count) % 2:
        at = int(np.argmax(np.linalg.norm(np.diff(result, axis=0), axis=1)))
        result = np.insert(result, at+1, .5*(result[at]+result[at+1]), axis=0)
    return result


def bounded_terminal_front(points, accuracy):
    """Approximate sub-accuracy steps without changing the authored geometry.

    Smooth curbs use 0.1*accuracy for sampling and another 0.1 for polyline
    simplification. Reserve the remaining 0.8 for this front approximation,
    keeping both outer terminal corners exact.
    Each retained chord bounds its entire continuous source polyline, not just
    the final mesh vertices. Larger features retain their original corners.
    """
    points = np.asarray(points, dtype=float)
    budget = .8*accuracy
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    result = _simplify(points, budget) if lengths.min() <= budget else points
    changed = len(result) != len(points)
    return result, dict(reference_vertices=len(points), retained_vertices=len(result),
        approximated=changed, maximum_front_chord_error_bound=budget if changed else 0.,
        source_curb_error_budget=.2*accuracy, outer_terminal_corners_exact=bool(
            np.array_equal(result[[0,-1]], points[[0,-1]])))


def seed_terminal_chart(mesh, rim_count):
    """Place incident free vertices on the inner side of each physical edge.

    A short reflex step can invert a Coons cell even though its topology is
    usable. This bounded geometric seed changes only interior coordinates;
    the normal convexity, shape, overlap and surface guards still follow.
    """
    try:
        mesh.validate()
        return mesh, dict(needed=False, passes=0, moved_vertices=0)
    except ValueError:
        pass
    xyz = mesh.xyz.copy()
    constraints = []
    for face in mesh.quads:
        for j in range(4):
            a, b = face[j], face[(j+1)%4]
            if a >= rim_count or b >= rim_count or (int(b)-int(a))%rim_count != 1:
                continue
            edge = xyz[b, :2]-xyz[a, :2]
            length = np.linalg.norm(edge)
            normal = np.array((-edge[1], edge[0]))/length
            for k in face:
                if k >= rim_count:
                    constraints.append((k, a, normal, .25*length))
    moved = set()
    for iteration in range(8):
        for k, a, normal, margin in constraints:
            distance = float((xyz[k, :2]-xyz[a, :2])@normal)
            if distance < margin:
                xyz[k, :2] += (margin-distance)*normal
                moved.add(int(k))
    result = Mesh(xyz, mesh.quads, mesh.material, dict(mesh.seams))
    result.validate()
    return result, dict(needed=True, passes=iteration+1, moved_vertices=len(moved),
                       exact_fixed_rim=bool(np.array_equal(result.xyz[:rim_count], mesh.xyz[:rim_count])))


class SkewMergedCapPlan(FiniteAxialPlan):
    coordinated_fibres = True
    terminal_margins = (.5, .75, 1., 1.5, 2., 2.5, 3., 4.)

    @classmethod
    def build(cls, context, bundle):
        from .bundle_compounds import compound_owners
        from .validation import validate_overlaps
        if compound_owners(context, bundle):
            return None
        failures = []
        for margin in cls.terminal_margins:
            try:
                plan = cls(context, bundle, margin_factor=margin)
                if plan.physical_roles_swapped:
                    raise ValueError('Finite terminal caps retain the declared physical lateral order')
                plan._prefer_fibre_strips(context)
                part = plan.pavement(context)
                validate_overlaps({bundle.key: part.mesh}, cell_size=max(1., plan.settings.spacing*2.))
                plan.cap_selection = dict(margin_factor=margin, prior_failures=failures)
                part.metadata['finite_cap_selection'] = plan.cap_selection
                return plan
            except ValueError as exc:
                failures.append(dict(margin_factor=margin, reason=str(exc)))
        raise AxialChartUnavailable('No qualified merged terminal cap: ' + '; '.join(
            f"{f['margin_factor']:g}: {f['reason']}" for f in failures))

    def __init__(self, context, bundle, *, margin_factor):
        self.union_profiles = {}
        self.union_cap_metadata = {}
        self.union_caps_prepared = False
        super().__init__(context, bundle, margin_factor=margin_factor)

    def chart_axis(self, context, directions):
        return unit(sum((p[-1]-p[0] for p in self.visible_points), np.zeros(2)))

    def prepare_curves(self, context):
        self.original_curves = tuple(p.copy() for p in self.curves)
        first = max(p[0, 0] for p in self.curves)
        last = min(p[-1, 0] for p in self.curves)
        planes = np.array((first, last)) + np.array((1., -1.))*self.margin_factor*self.settings.spacing
        values = np.asarray([interp(p, planes)[0] for p in self.curves])
        self.separated_ends = tuple(int(e) for e in np.flatnonzero(values[1]-values[2] > 1e-7))
        self.merged_ends = tuple(e for e in (0, 1) if e not in self.separated_ends)
        if not self.merged_ends:
            raise AxialChartUnavailable('This terminal branch needs an actual merged end')
        sites = {e for j in (*context.network.intersections, *context.network.connections) for e in j.ends}
        for end in self.merged_ends:
            if any(member_ends(m)[end] in sites for m in self.members):
                raise AxialChartUnavailable('Merged terminal caps need genuine free member ends')
            planes[end] = (first, last)[end] + (1. if end == 0 else -1.)*max(2., self.margin_factor)*self.settings.spacing
        if planes[1] <= planes[0] + 2*self.settings.spacing:
            raise AxialChartUnavailable('Terminal caps leave no positive common body')
        self.curves = [clip(p, *planes) for p in self.curves]
        values = np.asarray([interp(p, planes)[0] for p in self.curves])
        if any(values[1, e]-values[2, e] >= -1e-7 for e in self.merged_ends):
            raise AxialChartUnavailable('A merged cap needs a strictly overlapping body cross-section')
        self.routes = [Route([Piece(context[Key('path', t.section_id)], t.forward) for t in m.path])
                       for m in self.members]

    def strip_ends(self):
        return self.separated_ends

    def terminal_boundary(self, body, end):
        plane = (self.first, self.last)[end]
        regions = []
        for role in (0, 1):
            upper, lower = self.original_curves[2*role:2*role+2]
            upper = clip(upper, plane, upper[-1, 0]) if end else clip(upper, upper[0, 0], plane)
            lower = clip(lower, plane, lower[-1, 0]) if end else clip(lower, lower[0, 0], plane)
            regions.append(Region(Key('terminal_member', role), np.vstack((lower, upper[::-1]))))
        loops = RegionUnion().build(regions)
        if len(loops) != 1 or loops[0].hole:
            raise ValueError('Merged terminal cap needs one positive exterior union')
        boundary = loops[0].xy
        hits = np.flatnonzero(abs(boundary[:, 0]-plane) < 1e-8)
        start = hits[np.argmin(boundary[hits, 1])] if end else hits[np.argmax(boundary[hits, 1])]
        boundary = np.roll(boundary, -int(start), axis=0)
        outside = np.flatnonzero(abs(boundary[:, 0]-plane) > 1e-8)
        if not len(outside) or outside[-1]+1 >= len(boundary):
            raise ValueError('Terminal union lost its complete body-cut boundary')
        front = self.world(boundary[:outside[-1]+2, 0], boundary[:outside[-1]+2, 1])
        points = body.xyz[:, :2]-self.origin
        ids = np.flatnonzero(abs(points@self.axis-plane) < 1e-8)
        ids = ids[np.argsort(points[ids]@self.left)]
        if end:
            ids = ids[::-1]
        back = body.xyz[ids].copy()
        if (np.linalg.norm(front[-1]-back[0, :2]) > 1e-7
                or np.linalg.norm(front[0]-back[-1, :2]) > 1e-7):
            raise ValueError('Terminal union and body exterior corners disagree')
        # Source terminal cuts define the two semantic atlas corners, rather
        # than searching unrelated sets of four points on the polygon.
        ends = [self.original_curves[i][-1 if end else 0] for i in (3, 0)]
        corners = sorted(int(np.argmin(np.linalg.norm(front-self.world(p[0], p[1]), axis=1))) for p in ends)
        a, b = corners
        if not 0 < a < b < len(front)-1:
            raise ValueError('Terminal union does not expose both outer terminal corners')
        return front[:a+1], front[a:b+1], front[b:], back

    def make_union_cap(self, context, body, end):
        from .validation import validate_overlaps
        bottom, front, top, back = self.terminal_boundary(body, end)
        original_front, approximation = bounded_terminal_front(front, self.settings.curve_tolerance)
        shortest = float(np.linalg.norm(np.diff(original_front, axis=0), axis=1).min())
        # At a boundary vertex with one interior spoke, adjacent edge lengths
        # cannot differ by more than ratio_limit**2. Grade the retained short
        # feature's neighbors before allocating the direct 1:3 transition atlas.
        sampling_spacing = self.settings.spacing
        graded = shortest < sampling_spacing/48.
        length = max(np.linalg.norm(np.diff(p, axis=0), axis=1).sum() for p in (bottom, top))
        base = max(2, 2*int(np.ceil(length/(2*self.settings.spacing))))
        field, physical = ExtendedCrown(self), ExtendedCrown(self, physical=True)
        failures = []
        biases = [(0,0),(2,0),(0,2),(4,0),(0,4),(2,2)]
        lengths = np.linalg.norm(np.diff(original_front, axis=0), axis=1)
        initial_counts = np.maximum(1, np.ceil(lengths/sampling_spacing).astype(np.int64))
        chosen = np.argsort(-lengths, kind='stable')[:2]
        remaining = len(back)-1-int(initial_counts.sum()-initial_counts[chosen].sum())
        if len(chosen) == 1:
            biases.append((len(back)-1-int(initial_counts[chosen[0]]), 0))
        elif remaining >= 2 and not graded:
            first = int(round(remaining*lengths[chosen[0]]/lengths[chosen].sum()))
            for delta in (0, -1, 1, -2, 2, -3, 3, -4, 4):
                counts = np.array((first+delta, remaining-first-delta))
                if counts.min() > 0:
                    biases.append(tuple(int(v) for v in counts-initial_counts[chosen]))
        options = ((rows, bias) for bias in dict.fromkeys(biases)
                   for rows in dict.fromkeys((base, base+2, 2*base, 3*base)))
        for rows, bias in options:
            try:
                front = terminal_samples(original_front, sampling_spacing, len(back), bias, graded=graded)
                lower = sample_chain(bottom, rows, self.settings.curve_tolerance*.2)
                upper = sample_chain(top, rows, self.settings.curve_tolerance*.2)
                rim = np.vstack((lower[:-1], front[:-1], upper[:-1], back[:-1, :2]))
                corners = (0, rows, rows+len(front)-1, 2*rows+len(front)-1)
                mesh, counts = FourSidedQuadChart(context.patterns, max_cells=10000).map(rim, corners)
                seed = dict(needed=False, passes=0, moved_vertices=0)
                if graded:
                    mesh, seed = seed_terminal_chart(mesh, len(rim))
                mesh, correction = QuadShapeCorrection(max_cells=10000, max_iterations=1000).improve(mesh)
                if not correction['qualified']:
                    raise ValueError('Terminal cap needs better fixed-boundary quad grading')
                mesh.xyz[:, 2] = field(mesh.xyz[:, :2])
                outer_ids = np.r_[np.arange(rows+1), np.arange(corners[2], corners[3]+1)]
                mesh.xyz[outer_ids, 2] = self.datum
                boundary_ids = np.arange(len(rim)-len(back)+1, len(rim))
                # The final back corner is the first rim vertex.
                boundary_ids = np.r_[boundary_ids, 0]
                if len(boundary_ids) != len(back):
                    raise ValueError('Terminal atlas lost its complete body cross-section')
                if np.max(np.linalg.norm(mesh.xyz[boundary_ids, :2]-back[:, :2], axis=1)) > 1e-8:
                    raise ValueError('Terminal atlas changed its body boundary')
                mesh.xyz[boundary_ids] = back
                velocities = np.broadcast_to(self.axis, (len(mesh.xyz), 2))
                heights = physical(mesh.xyz[:, :2])
                heights[outer_ids] = self.datum
                join = EntranceSlopeGuard().validate(mesh, boundary_ids, velocities, heights)
                physical_mesh = Mesh(np.column_stack((mesh.xyz[:, :2], heights)), mesh.quads)
                slopes = physical_mesh.slope_statistics()
                incoming = max(2*float(np.max(abs(context[Key('profile', sid)].derivative(np.array((-1., 1.)))))) /
                               context.network.section_index[sid].width for sid in self.source_ids)
                limit = max(.1, 1.1*incoming)
                if slopes['maximum'] > limit+1e-8:
                    raise ValueError('Terminal cap exceeds its physical slope limit')
                validate_overlaps({self.bundle.key: mesh}, cell_size=max(1., self.settings.spacing*2.))
                key = (-1, end)
                self.cap_meshes[key] = mesh
                self.cap_queries.append(QuadHeightField(mesh, datum=self.datum,
                    ceiling=self.datum+self.rise*self.settings.exaggeration))
                self.cap_keys.append(key)
                self.cap_entrance_bounds[key] = join
                self.union_profiles[end] = back
                self.union_cap_metadata[end] = dict(quads=len(mesh.quads), rows=rows,
                    terminal_count_bias=bias, atlas_intervals=counts, correction=correction, physical_slopes=slopes,
                    terminal_approximation=approximation, terminal_sampling_spacing=sampling_spacing,
                    locally_graded_terminal=graded,
                    terminal_interior_seed=seed,
                    maximum_body_join_longitudinal_slope=join, prior_failures=failures)
                return
            except ValueError as exc:
                failures.append(dict(rows=rows, terminal_count_bias=bias, reason=str(exc)))
        raise ValueError('No qualified direct terminal atlas: ' + '; '.join(dict.fromkeys(f['reason'] for f in failures)))

    def keys(self, mesh):
        keys = super().keys(mesh)
        for end, profile in self.union_profiles.items():
            plane = (self.first, self.last)[end]
            ids = np.flatnonzero(abs((mesh.xyz[:, :2]-self.origin)@self.axis-plane) < 1e-8)
            for i in ids:
                distance = np.linalg.norm(profile-mesh.xyz[i], axis=1)
                at = int(np.argmin(distance))
                if distance[at] > 1e-8:
                    raise ValueError('Terminal/body seam has an undeclared cross-section vertex')
                keys[i] = ('merged terminal/body', end, at)
        return keys

    def pavement(self, context):
        if not self.union_caps_prepared:
            body = AttachedAxialPlan.pavement(self, context).mesh
            for end in self.merged_ends:
                self.make_union_cap(context, body, end)
            self.union_caps_prepared = True
        part = super().pavement(context)
        part.metadata.update(merged_terminal_caps=self.union_cap_metadata,
            merged_terminal_ends=list(self.merged_ends),
            construction='common-fibre body with original nonparallel terminal unions')
        return part

    def separator(self, context):
        return super().separator(context) if self.islands else None
