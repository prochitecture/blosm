"""Finite, profile-preserving quad surfaces across dimensional Section seams.

This numerical primitive owns the source trims, both analytic curb curves and
the query data needed by downstream separators and source-relative anchors.
It does not rewrite authored Section dimensions or triangulate its domain.
"""
from dataclasses import dataclass
import math
import numpy as np
from .continuations import WidthContinuation
from .mesh import Mesh, cross
from .model import Key
from .providers import Part
from .surface_query import QuadHeightField


def product_key(owner):
    return Key('seam_surface' if owner.kind == 'connection' else 'profile_transition', owner.id)


def output_key(owner):
    return Key('seam_mesh' if owner.kind == 'connection' else 'junction_mesh', owner.id)


def composite_owns_transition(context, owner, bundle):
    if bundle is not None:
        plan = context[Key('bundle_plan', bundle.id)]
        if plan.mode == 'union' and owner in plan.transition_owners:
            return True
    return False


def pavement_signature(context, owner, bundle=None):
    return context.digest(product_key(owner)), composite_owns_transition(context,owner,bundle)


def published_pavement(context, owner, bundle=None):
    """One physical surface owner, while retaining source transition queries."""
    if composite_owns_transition(context,owner,bundle):
        return None
    product=context[product_key(owner)]
    if product is None:
        return None  # An explicit junction footprint owns this dimensional seam.
    if product.pavement.mesh is None:
        # A composite may stop owning an otherwise unchanged transition after
        # Blender has released its old numerical buffers. Reconstruct from the
        # certified plan without mutating a committed cached product.
        return Part(product.pavement.owner,product.plan.mesh(context.patterns,context.network.settings.exaggeration),
                    dict(product.pavement.metadata))
    return product.pavement


def mixed_bend_required(first,first_end,last,last_end):
    """Predict the mixed endpoint frames produced by PathProvider.

    Smooth paths use the common bisector, while an Exact endpoint retains its
    source segment direction. Equal dimensions therefore do not guarantee a
    compatible cross-section. Collinear and same-mode seams keep their existing
    dispatch, without constructing paths to decide the task graph.
    """
    if first.interpolation==last.interpolation:
        return False
    from .paths import unit
    def outward(source,end):
        points=np.asarray(source.points)[:,:2]
        return unit(points[1]-points[0]) if end==0 else unit(points[-2]-points[-1])
    a,b=outward(first,first_end),outward(last,last_end)
    common=unit(a-b)
    a=a if first.interpolation==0 else common
    b=b if last.interpolation==0 else -common
    return bool(np.linalg.norm(a+b)>1e-5)


def required(network, site):
    """Source-only dispatch; matching dimensions and frames need no area."""
    if len(site.ends) != 2:
        return False
    from .profiles import memberships, ProfileProvider
    members, sources = memberships(network), network.section_index
    ends = site.ends
    a, b = (sources[e.section_id] for e in ends)
    pa, pb = (ProfileProvider.build(sources[e.section_id], members.get(e.section_id), sources).oriented(e.end)
              for e in ends)
    oa, ob = (sources[e.section_id].offset * (1 if e.end == 0 else -1) for e in ends)
    return (abs(a.width - b.width) > 1e-6 or abs(oa + ob) > 1e-6
            or abs(pa.rise - pb.rise) > 1e-7 or abs(pa.datum - pb.datum) > 1e-8
            or pa.half + pb.half != 0
            or mixed_bend_required(a,ends[0].end,b,ends[1].end))


def neighboring_junctions(attachments, site):
    """Potential physical constraints through explicit compatible-route incidence.

    Candidate routes can stop before reaching this Section when dimensions
    change. Their actual footprint mouths decide whether a constraint exists.
    No proximity search or incidental cache contents establish ownership.
    """
    result = {}
    for own in site.ends:
        for junction, end in attachments[own.section_id]:
            if len(junction.ends) > 2 and end != own.end:
                result[junction.key] = junction
    return tuple(result[key] for key in sorted(result))


@dataclass(frozen=True)
class TransitionClearance:
    limits: np.ndarray
    constraints: tuple

    @classmethod
    def build(cls, context, approaches, neighbors=()):
        from .ownership import footprint_key
        limits = .4 * np.array([a.length for a in approaches])
        constraints = []
        for junction in neighbors:
            footprint = context[footprint_key(junction.key)]
            for i, approach in enumerate(approaches):
                mouth = footprint.mouths.get((approach.source.id, 1 - approach.end))
                if mouth is None:
                    continue
                remaining = mouth.station if approach.end == 0 else approach.length - mouth.station
                gap = min(.25, .05 * approach.width)
                available = max(0., remaining - gap)
                limits[i] = min(limits[i], available)
                constraints.append(dict(section=approach.source.id, end=approach.end,
                    junction=junction.key.text(), mouth_station=mouth.station,
                    remaining_source_length=float(remaining), clearance=gap, maximum_distance=float(available)))
        if np.any(limits <= 1e-7):
            raise ValueError('A neighboring junction consumes the Section transition span; shared region ownership is required')
        limits.flags.writeable = False
        return cls(limits, tuple(constraints))


def _certify_orientation(first, second, *, budget=8192):
    """Certify the complete ruled XY map, not just its sampled mesh corners.

    det(X_t, X_v) is affine across a ruled cross-section. At either curb it is
    a Bernstein polynomial. Subdivision bounds both complete edge
    polynomials and therefore every transverse coordinate between them.
    """
    if hasattr(first, 'segments'):
        if not hasattr(second, 'segments') or not np.array_equal(first.knots, second.knots):
            raise ValueError('Piecewise transition curbs need identical parameter intervals')
        checks, lower = 0, math.inf
        for i,(a,b) in enumerate(zip(first.segments,second.segments)):
            report = _certify_orientation(a,b,budget=budget-checks)
            checks += report['orientation_checks']
            lower = min(lower,report['minimum_jacobian_bound']/(first.knots[i+1]-first.knots[i]))
        return dict(minimum_jacobian_bound=lower,orientation_checks=checks,
                    orientation_pieces=len(first.segments))
    width = second.controls - first.controls
    pending = []
    for curb in (first, second):
        dv, dw = len(curb.first_derivative)-1, len(width)-1
        coefficients = np.zeros(dv+dw+1)
        for i, velocity in enumerate(curb.first_derivative):
            for j, delta in enumerate(width):
                coefficients[i + j] -= (math.comb(dv, i) * math.comb(dw, j)
                    / math.comb(dv+dw, i + j)) * cross(velocity, delta)
        pending.append((coefficients, 0))
    checks, lower = 0, math.inf
    while pending:
        coefficients, depth = pending.pop()
        checks += 1
        if checks > budget:
            raise ValueError('Section transition orientation exceeded its proof budget')
        if coefficients.min() > 1e-10:
            lower = min(lower, float(coefficients.min()))
            continue
        if coefficients.max() <= 1e-10 or depth >= 16:
            raise ValueError('Section transition has no certified positive longitudinal chart')
        work, left, right = coefficients, [coefficients[0]], [coefficients[-1]]
        while len(work) > 1:
            work = .5 * (work[:-1] + work[1:])
            left.append(work[0])
            right.append(work[-1])
        pending.extend(((np.asarray(left), depth + 1), (np.asarray(right[::-1]), depth + 1)))
    return dict(minimum_jacobian_bound=lower, orientation_checks=checks)


@dataclass
class SectionTransition:
    owner: Key
    approaches: tuple
    distances: tuple
    curbs: tuple
    rows: int
    across: int
    parameters: np.ndarray
    certification: dict
    height_corrections: np.ndarray = None

    @classmethod
    def solve(cls, context, site, approaches, clearance=None, *, _local_curvature=False):
        first, last = approaches
        if abs(first.profile.datum - last.profile.datum) > 1e-8:
            raise ValueError('A Section transition requires a common minimum height')
        from .exact_transitions import exact_path, ExactTransitionGeometry, paired_parameters
        exact = tuple(exact_path(a.path) for a in approaches)
        if any(exact) and not all(exact) and np.linalg.norm(first.direction + last.direction) > 1e-5:
            from .mixed_transitions import MixedTransitionFitter
            return MixedTransitionFitter.solve(context, site, approaches, clearance)
        # Keep a complete profile at both new mouths. Each trim remains inside
        # its own finite source, leaving room for the other end's features.
        differences = np.array([first.profile.value(-q) - last.profile.value(q) for q in (-1., 0., 1.)])
        a = .5 * (differences[0] + differences[2]) - differences[1]
        b = .5 * (differences[2] - differences[0])
        profile_delta = float(np.max(abs(differences)))
        if abs(a) > 1e-15:
            vertex = -b / (2. * a)
            if abs(vertex) < 1.:
                profile_delta = max(profile_delta, abs(float(first.profile.value(-vertex) - last.profile.value(vertex))))
        desired = max(2., 4. * abs(first.width - last.width) + 2. * abs(first.offset + last.offset),
                      1.875 * profile_delta / .03)
        clearance = clearance if clearance is not None else TransitionClearance.build(context, approaches)
        limits = clearance.limits
        from .intersections import half_intervals
        across = 2 * half_intervals(context.network.settings)
        failures = []
        previous = None
        for factor in (1., 1.5, 2., 3., 4.):
            distances = np.minimum(.5 * desired * factor, limits)
            if previous is not None and np.array_equal(distances, previous):
                continue
            previous = distances.copy()
            try:
                geometry = {}
                use_miter = all(exact) and ExactTransitionGeometry.required(approaches,distances)
                if use_miter:
                    curbs, geometry = ExactTransitionGeometry.build(approaches,distances)
                else:
                    curbs = []
                    for q in (-1., 1.):
                        a, va = first.jet(float(distances[0]), first.offset - q * .5 * first.width)
                        b, vb = last.jet(float(distances[1]), last.offset + q * .5 * last.width)
                        curbs.append(WidthContinuation(a, b, -va, vb, getattr(site, 'radius', 2.),
                                                       refine_curvature=_local_curvature))
                certification = _certify_orientation(*curbs)
                certification.update(geometry)
                if _local_curvature:
                    certification['curvature_refinement'] = [dict(intervals=c.curvature_intervals,
                        bound=c.curvature_bound) for c in curbs]
                settings = context.network.settings
                spacing = min(settings.spacing, 6. * min(first.width, last.width) / across)
                rows = (0 if use_miter else max(c.intervals(spacing, settings.curve_tolerance) for c in curbs))
                if rows * across > 50_000:
                    raise ValueError('Section transition exceeds its cell budget')
                join = float(distances[0] / distances.sum())
                # Include the ownership join without inserting an arbitrarily
                # thin row next to an existing uniform sample. Neither side's
                # parameter interval exceeds the certified 1/rows bound.
                count_a, count_b = max(1, math.ceil(rows * join)), max(1, math.ceil(rows * (1. - join)))
                # Each independently owned half-rim must have even boundary
                # parity, just like a road span or a certified junction curb.
                # Otherwise an attached median can have an odd closed boundary
                # and admit no conforming quad-only fill. Refinement preserves
                # every analytical curve and improves its sampling bound.
                count_a += count_a % 2
                count_b += count_b % 2
                parameters = np.r_[np.linspace(0., join, count_a + 1)[:-1],
                                    np.linspace(join, 1., count_b + 1)]
                if use_miter:
                    parameters = paired_parameters(curbs,spacing,settings.curve_tolerance,
                                                   max_intervals=50000//across)
                parameters.flags.writeable = False
                rows = len(parameters) - 1
                if rows * across > 50_000:
                    raise ValueError('Section transition exceeds its cell budget')
                result = cls(site.key, tuple(approaches), tuple(map(float, distances)), tuple(curbs),
                             rows, across, parameters, certification)
                return cls.certify(context, result, clearance)
            except ValueError as error:
                failures.append(str(error))
        if all(exact):
            from .transition_fitting import ExactTransitionFitter
            fitted = ExactTransitionFitter.fit(context, site, approaches, clearance, desired, across, failures)
            if fitted is not None:
                return fitted
        elif not any(exact) and not _local_curvature:
            # Preserve every previously accepted primary candidate. Only a
            # failed Smooth scan retries the same finite trims/polynomials
            # with tighter complete-curve proofs, retaining all mesh guards.
            return cls.solve(context, site, approaches, clearance, _local_curvature=True)
        raise ValueError('No bounded dimensional Section transition fits its source spans: ' + '; '.join(dict.fromkeys(failures)))

    @staticmethod
    def certify(context, result, clearance):
        """The same physical mesh guards for regular and fitted row partitions."""
        approaches, curbs = result.approaches, result.curbs
        rows, across = result.rows, result.across
        mesh = result.mesh(context.patterns, exaggeration=1.)
        quality = mesh.validate()
        if quality['maximum_edge_ratio'] > 8.:
            raise ValueError('Section transition quad edge ratio exceeds 8')
        incoming = max(float(np.max(abs(a.profile.derivative(np.array((-1., 1.))))))
                       * 2. / a.width for a in approaches)
        limit = max(.03, 1.1 * incoming)
        slopes = mesh.slope_statistics()
        if slopes['maximum'] > limit + 1e-8:
            raise ValueError('Section transition exceeds the incoming-profile slope limit')
        from .slope_projection import corner_matrices
        gradients = np.einsum('fcij,fj->fci', corner_matrices(mesh), mesh.xyz[mesh.quads, 2])
        join = 0.
        for end, cells, corners in ((0, slice(0, across), (0, 3)),
                                    (1, slice(-across, None), (1, 2))):
            t = float(end)
            q = 2. * context.patterns.grid(rows, across).uv[mesh.quads[cells]][:, corners, 1] - 1.
            velocity = ((1. - .5 * (q + 1.))[..., None] * curbs[0].derivative(t)
                        + (.5 * (q + 1.))[..., None] * curbs[1].derivative(t))
            velocity /= np.linalg.norm(velocity, axis=-1)[..., None]
            join = max(join, float(np.max(abs(np.sum(gradients[cells][:, corners] * velocity, axis=-1)))))
        if join > .03 + 1e-8:
            raise ValueError('Section transition entrance longitudinal slope exceeds 3 percent')
        result.certification.update(quality=quality, physical_slopes=slopes, maximum_slope=limit,
                             maximum_entrance_longitudinal_slope=join,
                             clearance_constraints=list(clearance.constraints),
                             maximum_trim_distances=list(map(float, clearance.limits)))
        return result

    @property
    def datum(self):
        return self.approaches[0].profile.datum

    @property
    def ceiling(self):
        return self.datum + max(a.profile.rise for a in self.approaches)

    def surface(self, t, q):
        t, q = np.broadcast_arrays(np.asarray(t), np.asarray(q))
        v = .5 * (q + 1.)
        xy = (1. - v[..., None]) * self.curbs[0](t) + v[..., None] * self.curbs[1](t)
        weight = t ** 3 * (10. + t * (-15. + 6. * t))
        a, b = self.approaches
        z = self.datum + (1. - weight) * a.profile.value(-q) + weight * b.profile.value(q)
        if self.height_corrections is not None:
            # Interpolate only the interior residual. Its two end rows, curbs
            # and crest are zero, preserving complete analytical mouth profiles
            # for arbitrary q as well as the copied published mouth samples.
            i = np.clip(np.searchsorted(self.parameters, t, side='right') - 1, 0, self.rows - 1)
            along = np.clip((t - self.parameters[i]) / (self.parameters[i + 1] - self.parameters[i]), 0., 1.)
            transverse = np.clip(v * self.across, 0., self.across)
            j = np.clip(transverse.astype(int), 0, self.across - 1)
            across = transverse - j
            h = self.height_corrections
            z = z + ((1. - along) * ((1. - across) * h[i, j] + across * h[i, j + 1])
                     + along * ((1. - across) * h[i + 1, j] + across * h[i + 1, j + 1]))
        return xy, z

    def mesh(self, patterns, exaggeration):
        pattern = patterns.grid(self.rows, self.across)
        t, v = pattern.uv.T
        t = np.interp(t, np.linspace(0., 1., self.rows + 1), self.parameters)
        xy, z = self.surface(t, 2. * v - 1.)
        result = Mesh(np.column_stack((xy, self.datum + exaggeration * (z - self.datum))), pattern.quads[:, ::-1])
        if self.certification.get('exact_miter') or self.certification.get('mixed_miter'):
            # Copy every mouth sample from the source frame, avoiding tiny
            # interpolation differences across separately published objects.
            q = np.linspace(-1.,1.,self.across+1)
            for i,approach in enumerate(self.approaches):
                side = -q if i==0 else q
                row = slice(0,self.across+1) if i==0 else slice(-self.across-1,None)
                result.xyz[row,:2] = approach.point(self.distances[i],approach.offset+.5*approach.width*side)
                result.xyz[row,2] = approach.profile.datum+exaggeration*approach.profile.value(side)
        result.seams = {0: np.arange(self.across + 1),
                        1: np.arange(self.rows * (self.across + 1), (self.rows + 1) * (self.across + 1))}
        return result

    def source_index(self, section_id):
        return next(i for i, a in enumerate(self.approaches) if a.source.id == section_id)

    def source_half_rim(self, section_id, source_side, exaggeration):
        i = self.source_index(section_id)
        a = self.approaches[i]
        outward_side = source_side * (1 if a.end == 0 else -1)
        q = -outward_side if i == 0 else outward_side
        join = self.distances[0] / sum(self.distances)
        t = self.parameters[self.parameters <= join + 1e-12] if i == 0 else self.parameters[self.parameters >= join - 1e-12]
        # First orient from the source seam into this Section; then orient in
        # the Section's authored direction. The common join is a real mesh row.
        if i == 0:
            t = t[::-1]
        if a.end == 1:
            t = t[::-1]
        xy, z = self.surface(t, q)
        return np.column_stack((xy, self.datum + exaggeration * (z - self.datum)))


@dataclass
class TransitionProduct:
    plan: SectionTransition
    pavement: Part
    field: QuadHeightField

    @classmethod
    def build(cls, context, plan):
        settings = context.network.settings
        mesh = plan.mesh(context.patterns, settings.exaggeration)
        field = QuadHeightField(mesh, datum=plan.datum,
            ceiling=plan.datum + settings.exaggeration * (plan.ceiling - plan.datum))
        construction = ('certified mixed-source polynomial Section transition' if plan.certification.get('mixed_miter')
                        else 'certified piecewise Exact miter Section transition' if plan.certification.get('exact_miter')
                        else 'certified quintic Section transition')
        metadata = dict(plan.certification, construction=construction,
            source_trims=[dict(section=a.source.id, end=a.end, distance=d) for a, d in zip(plan.approaches, plan.distances)],
            crown_ceiling=plan.ceiling, rows=plan.rows, transverse_intervals=plan.across)
        return cls(plan, Part(plan.owner, mesh, metadata), field)

    def height(self, xy):
        return self.field.evaluate(xy, default=self.plan.datum)

    def anchor(self, placement, exaggeration):
        from .providers import AnchorTransform
        plan = self.plan
        i = plan.source_index(placement.section_id)
        approach = plan.approaches[i]
        distance = placement.station if approach.end == 0 else approach.length - placement.station
        if distance < -1e-8 or distance > plan.distances[i] + 1e-8:
            raise ValueError('Furniture station lies outside its owned Section transition')
        if plan.certification.get('mixed_miter'):
            from .mixed_transitions import MixedTransitionGeometry
            t = MixedTransitionGeometry.parameter_at_station(plan,i,placement.station)
        else:
            join = plan.distances[0] / sum(plan.distances)
            t = join * (1. - distance / plan.distances[0]) if i == 0 else join + (1. - join) * distance / plan.distances[1]
        t = float(np.clip(t, 0., 1.))
        direction = (-1 if i == 0 else 1) * (1 if approach.end == 0 else -1)
        side = (1 if placement.side == 'RIGHT' else -1) * direction
        left, right = (curb(np.asarray(t)) for curb in plan.curbs)
        transverse = right - left
        transverse /= np.linalg.norm(transverse)
        position = (right if side == 1 else left) + side * placement.offset * transverse
        _, boundary_height = plan.surface(np.asarray(t), np.asarray(float(side)))
        default = plan.datum + exaggeration * (float(boundary_height) - plan.datum)
        height = float(self.field.evaluate(position[None], default=default)[0]) + placement.height
        tangent = direction * sum(curb.derivative(np.asarray(t)) for curb in plan.curbs)
        return AnchorTransform(tuple(map(float, position)) + (height,),
                               float(math.atan2(tangent[1], tangent[0])), placement.asset)

    def release_geometry(self):
        if self.pavement.mesh is not None:
            self.pavement.metadata['bounds'] = self.pavement.mesh.bounds
            self.pavement.mesh = None


def separator_range(path, metadata):
    """Keep centerline contact proofs through transition-owned source spans."""
    first, last = metadata['start'], metadata['end']
    for span in metadata.get('transition_spans', ()):
        if span['source_end'] == 0:
            first = 0.
        else:
            last = path.length
    return first, last


def extend_rim(context, section_id, side, rim):
    """Include both halves of a seam without gaps, duplicate faces or T joins."""
    metadata = context[Key('road', section_id)].metadata
    from .bundles import _join
    for span in metadata.get('transition_spans', ()):
        half = context[Key.parse(span['owner'])].plan.source_half_rim(
            section_id, side, context.network.settings.exaggeration)
        rim = _join(half, rim) if span['source_end'] == 0 else _join(rim, half)
    return rim
