"""Topology-owned Bundle regions sharing the actual pavement boundary samples.

The separated branch uses ruled members and four-sided direct quad patterns.
It never constructs a hand-drawn median outline or triangulates a polygon.
Intersection incidence selects the rounded caps; spatial proximity cannot
silently attach a separator to an unrelated road. More complex union branches
must explicitly replace this branch when pavement interiors overlap.
"""
from dataclasses import dataclass
import math
import numpy as np
from .intersections import IntersectionProvider
from .mesh import Mesh, cross
from .model import End, Key
from .providers import Part, RoadSampling
from .profiles import plateau_extension
from .updates import FeatureProvider, Task


def _join(*parts):
    result = []
    for part in parts:
        part = np.asarray(part)
        if not len(part):
            continue
        if result:
            if np.linalg.norm(result[-1][-1] - part[0]) > 1e-7:
                raise ValueError('Explicit Bundle boundary incidences do not meet')
            part = part[1:]
        if len(part):
            result.append(part)
    return np.concatenate(result)


def _beziers(path):
    if isinstance(path, CenterlineRange):
        return path.curves
    if path.exact:
        return [np.linspace(a, b, 4) for a, b in zip(path.points, path.points[1:])]
    return [np.array((p, p + a / 3., p + (2. * a + b) / 3., p + a + b + c))
            for p, a, b, c in path.coefficients]


@dataclass(frozen=True)
class CenterlineRange:
    """Exact polynomial restrictions to the visible, junction-trimmed span."""
    curves: tuple

    @classmethod
    def build(cls, path, first, last):
        if first < -1e-8 or last > path.length + 1e-8 or last <= first + 1e-8:
            raise ValueError('Centerline separation needs a positive visible station span')
        a, b = float(path.parameter(first)), float(path.parameter(last))
        curves = []
        for i, (lo, hi) in enumerate(zip(path.knots, path.knots[1:])):
            start, end = max(a, lo), min(b, hi)
            if end <= start + 1e-12:
                continue
            u, v = (start - lo) / (hi - lo), (end - lo) / (hi - lo)
            if path.exact:
                p, q = path.jet(np.array((start, end)))[0]
                curves.append(np.linspace(p, q, 4))
            else:
                p, d, e, f = path.coefficients[i]
                point = lambda t: p + t * (d + t * (e + t * f))
                derivative = lambda t: d + t * (2. * e + 3. * t * f)
                curves.append(np.array((point(u), point(u) + (v - u) / 3. * derivative(u),
                    point(v) - (v - u) / 3. * derivative(v), point(v))))
        return cls(tuple(curves))


def _split(curve):
    a, b, c = .5 * (curve[:-1] + curve[1:])
    d, e = .5 * (a + b), .5 * (b + c)
    f = .5 * (d + e)
    return np.array((curve[0], a, d, f)), np.array((f, e, c, curve[-1]))


def _point_distance(points, a, b):
    delta = b - a
    length2 = float(delta @ delta)
    t = np.clip((points - a) @ delta / max(length2, 1e-30), 0., 1.)
    return np.linalg.norm(points - a - t[:, None] * delta, axis=1)


def _chord_distance(a, b):
    p, r, q, s = a[0], a[-1] - a[0], b[0], b[-1] - b[0]
    denominator = float(cross(r, s))
    if abs(denominator) > 1e-20:
        t, u = cross(q - p, s) / denominator, cross(q - p, r) / denominator
        if 0. <= t <= 1. and 0. <= u <= 1.:
            return 0.
    return min(float(_point_distance(a[[0, -1]], b[0], b[-1]).min()),
               float(_point_distance(b[[0, -1]], a[0], a[-1]).min()))


def certify_disjoint(first, second, tolerance=1e-5):
    """Conservative cubic-hull proof; never accept an unresolved contact.

    The convex hull is within the maximal control-to-chord distance of its
    chord. Chord distance minus both bounds therefore certifies separation.
    A bounded de Casteljau search resolves close curved cases independently
    of the much coarser visual mesh samples.
    """
    seeds = ((a, b, 0) for p in first for q in second for a in _beziers(p) for b in _beziers(q))
    pairs = []
    checks = 0
    while True:
        if not pairs:
            seed = next(seeds, None)
            if seed is None:
                break
            pairs.append(seed)
        a, b, depth = pairs.pop()
        checks += 1
        if checks > 200_000:
            raise ValueError('Bundle centerline separation exceeded its proof budget')
        gap = np.maximum(a.min(axis=0) - b.max(axis=0), b.min(axis=0) - a.max(axis=0))
        if np.any(gap > tolerance):
            continue
        error_a = float(_point_distance(a[1:-1], a[0], a[-1]).max())
        error_b = float(_point_distance(b[1:-1], b[0], b[-1]).max())
        distance = _chord_distance(a, b)
        if distance - error_a - error_b > tolerance:
            continue
        if distance + error_a + error_b <= tolerance:
            raise ValueError('Bundle centerlines cross or touch')
        if depth >= 28:
            raise ValueError('Bundle centerline contact cannot be safely resolved')
        if error_a >= error_b:
            pairs.extend((part, b, depth + 1) for part in _split(a))
        else:
            pairs.extend((a, part, depth + 1) for part in _split(b))
    return checks


def _road_rim(context, section_id, side):
    section = context.network.section_index[section_id]
    path, road = context[Key('path', section_id)], context[Key('road', section_id)]
    profile = context[Key('profile', section_id)]
    samples = RoadSampling.from_metadata(road.metadata).samples(path, context.network.settings)
    xy = path.point(samples, section.offset + side * .5 * section.width)
    z = np.full(len(samples), profile.datum + context.network.settings.exaggeration * float(profile.value(side)))
    from .transitions import extend_rim
    return extend_rim(context, section_id, side, np.column_stack((xy, z)))


@dataclass
class Member:
    path: tuple
    rim: np.ndarray
    ends: tuple

    @classmethod
    def build(cls, context, member, role):
        paths, rims = [], []
        previous = None
        for traversal in member.path:
            path = context[Key('path', traversal.section_id)]
            road = context[Key('road', traversal.section_id)]
            from .transitions import separator_range
            paths.append(CenterlineRange.build(path, *separator_range(path, road.metadata)))
            if previous is not None:
                a = previous.points[-1 if last_forward else 0]
                b = path.points[0 if traversal.forward else -1]
                if np.linalg.norm(a - b) > 1e-7:
                    raise ValueError('Consecutive Bundle Sections have inconsistent endpoint copies')
            side = (1 if role == 0 else -1) * (1 if traversal.forward else -1)
            rim = _road_rim(context, traversal.section_id, side)
            rims.append(rim if traversal.forward else rim[::-1])
            previous, last_forward = path, traversal.forward
        a, b = member.path[0], member.path[-1]
        ends = (End(a.section_id, 0 if a.forward else 1), End(b.section_id, 1 if b.forward else 0))
        return cls(tuple(paths), _join(*rims), ends)


@dataclass
class Cap:
    first: np.ndarray
    middle: np.ndarray
    last: np.ndarray
    attached: bool

    @classmethod
    def free(cls, a, b, count):
        return cls(a[None, :], np.linspace(a, b, count + 1), b[None, :], False)

    @classmethod
    def rounded(cls, curve):
        if len(curve) < 3:
            raise ValueError('A rounded Bundle cap needs an interior curve sample')
        middle = len(curve) // 2
        return cls(curve[:middle + 1], curve[middle:middle + 1], curve[middle:], True)

    @property
    def rim(self):
        return _join(self.first, self.middle, self.last)


class BundleProvider(FeatureProvider):
    source_kinds = frozenset(('bundle',))

    def tasks(self, network):
        from .ownership import PLAN, SURFACES_READY, has_junctions
        from .transitions import required, product_key
        from .separator_ports import PORTS, SeparatorPortIndex
        if network.bundles:
            port_deps = tuple(sorted({*(j.layout_key for j in network.intersections),
                                      *((PLAN, SURFACES_READY) if has_junctions(network) else ())}))
            yield Task(PORTS, port_deps, SeparatorPortIndex.build)
        from .bundle_merged_entrances import MergedEntrancePreparation, preparation_key
        merged_candidates = {candidate.id for candidate, crossing, runs in MergedEntrancePreparation.candidates(network)}
        for bundle in network.bundles:
            from .bundle_regions import continuation_sections
            boundary_sections = set(bundle.sections) | continuation_sections(network, bundle)
            attached = [j for j in (*network.intersections, *network.connections)
                        if any(e.section_id in boundary_sections for e in j.ends)]
            sections = boundary_sections.copy()
            sections.update(e.section_id for j in attached for e in j.ends)
            from .bundle_axis import plan
            from .bundle_ends import mouth_dependencies
            from .bundle_boundaries import dimensional_sites
            members = {t.section_id for m in bundle.members for t in m.path}
            ends = {End(visit.section_id, end) for member in bundle.members for visit in member.path for end in (0, 1)}
            full_ends = {end for junction in attached if len(junction.ends) > 2 for end in junction.ends}
            angular_ports = (len(bundle.members) == 2 and not bundle.interior_section_ids
                             and all(member.complete and len(member.path) == 1 for member in bundle.members)
                             and ends <= full_ends)
            internal_ports = (bool(bundle.interior_section_ids) and len(bundle.members) == 2
                              and all(member.complete for member in bundle.members))
            plan_deps = tuple(sorted({bundle.key,
                *(Key(kind, sid) for sid in members for kind in ('path', 'profile')),
                *mouth_dependencies(network, bundle),
                *(j.layout_key for j in dimensional_sites(network, bundle)),
                *(j.key for j in attached),
                *((PORTS,) if angular_ports or internal_ports else ()),
                *((preparation_key(bundle),) if bundle.id in merged_candidates else ()),
                *((PLAN,) if has_junctions(network) else ())}))
            yield Task(Key('bundle_plan', bundle.id), plan_deps,
                lambda context, source=bundle: plan(context, source),
                signature=lambda context, source=bundle, dependencies=plan_deps:
                    self.plan_signature(context, source, dependencies))
            yield Task(Key('bundle_pavement', bundle.id), (Key('bundle_plan', bundle.id),),
                lambda context, source=bundle: context[Key('bundle_plan', source.id)].pavement(context)
                    if context[Key('bundle_plan', source.id)].mode == 'union' else None)
            deps = tuple(sorted({bundle.key, PORTS, Key('bundle_plan', bundle.id), *(Key(kind, sid) for sid in sections
                for kind in ('path', 'profile', 'road')), *(j.layout_key for j in attached),
                *(product_key(j.key) for j in attached if required(network, j)),
                *((PLAN, SURFACES_READY) if has_junctions(network) else ())}))
            yield Task(Key('separator', bundle.id), deps,
                lambda context, source=bundle: self.build(context, source),
                signature=(lambda context, source=bundle, dependencies=deps:
                    self.signature(context, source, dependencies)))

    @staticmethod
    def plan_signature(context, bundle, dependencies):
        from .ownership import PLAN
        from .separator_ports import PORTS
        from .bundle_compounds import compound_owners
        # Finite strips must yield when a compound consumes their mouth. The
        # graph records that ownership dependency, while a distant compound
        # edit must not invalidate an unrelated Bundle plan.
        return ([(key.text(), context.digest(key)) for key in dependencies if key not in (PLAN, PORTS)],
                [owner.text() for owner in compound_owners(context, bundle)],
                context[PORTS].signature(bundle.sections) if PORTS in dependencies else None)

    @staticmethod
    def signature(context, bundle, dependencies):
        from .ownership import PLAN, SURFACES_READY, surface_key
        from .bundle_compounds import compound_owners
        from .bundle_regions import continuation_sections
        from .separator_ports import PORTS
        sections = set(bundle.sections) | continuation_sections(context.network, bundle)
        return ([(key.text(), context.digest(key)) for key in dependencies if key not in (PLAN, SURFACES_READY, PORTS)],
                context[PORTS].signature(sections),
                [(owner.text(), context.digest(surface_key(owner))) for owner in compound_owners(context, bundle)])

    @staticmethod
    def cap(context, first, second, end, rim_coordinates=None, rim_source=None):
        from .separator_ports import PORTS
        ends = (first.ends[end], second.ends[end])
        ports = [context[PORTS].resolve(e) for e in ends]
        matches = [port.site if port is not None else None for port in ports]
        if not any(matches):
            if any(e in c.ends for c in (*context.network.intersections, *context.network.connections) for e in ends):
                raise ValueError('Attached Bundle endpoint has no visible pavement port; explicit continuation is required')
            return None
        if not all(matches):
            raise ValueError('A partly attached Bundle cap requires explicit junction continuation')
        if matches[0].key == matches[1].key:
            junction = matches[0]
            layout = context[junction.layout_key]
            indices = [port.approach_index for port in ports]
            from .ownership import CompoundProduct, surface_key
            product = context.products.get(surface_key(junction.key))
            if isinstance(product, CompoundProduct):
                actual = [list(layout.approaches[i].path.trims(layout.outward_stations[i]).items())[-1][0]
                          for i in indices]
                return Cap.rounded(product.cap(*actual))
            for side, corner in enumerate(layout.corners):
                a, b = layout.approaches[indices[0]], layout.approaches[indices[1]]
                forward = corner.first is a and corner.second is b and a.profile.peak == -1. and b.profile.peak == 1.
                reverse = corner.first is b and corner.second is a and b.profile.peak == -1. and a.profile.peak == 1.
                if forward or reverse:
                    t = np.linspace(0., 1., layout.rows + 1)
                    xy, z = IntersectionProvider.band_surface(layout, side, t, np.zeros_like(t))
                    xyz = np.column_stack((xy, a.profile.datum + context.network.settings.exaggeration * z))
                    return Cap.rounded(xyz if forward else xyz[::-1])
            raise ValueError('A merge junction has no adjacent inner-profile curb cap')
        corners, others, coordinates = [], [], []
        for member_index, (endpoint, junction) in enumerate(zip(ends, matches)):
            layout = context[junction.layout_key]
            index = ports[member_index].approach_index
            approach = layout.approaches[index]
            outward_q = approach.profile.peak if rim_coordinates is None else rim_coordinates[member_index] * (1 if endpoint.end == 0 else -1)
            if outward_q not in (-1., 1.):
                raise ValueError('A Bundle cap must select an explicit physical carriageway side')
            preceding = outward_q == 1.
            side = (index - 1) % len(layout.approaches) if preceding else index
            corner = layout.corners[side]
            other = corner.first if preceding else corner.second
            def actual_end(index):
                a = layout.approaches[index]
                return list(a.path.trims(layout.outward_stations[index]).items())[-1][0]
            own_end = actual_end(index)
            other_end = actual_end(next(i for i, a in enumerate(layout.approaches) if a is other))
            from .ownership import CompoundProduct, surface_key
            product = context.products.get(surface_key(junction.key))
            if isinstance(product, CompoundProduct):
                corners.append(product.edge_cap(own_end, other_end, 1. if preceding else 0., 0. if preceding else 1.))
            else:
                t = np.linspace(0., 1., layout.rows + 1)
                xy, z = IntersectionProvider.band_surface(layout, side, t, np.zeros_like(t))
                ceiling = max(a.profile.rise for a in layout.approaches)
                xyz = np.column_stack((xy, approach.profile.datum + context.network.settings.exaggeration * np.clip(z, 0., ceiling)))
                # Orient both rims from the member peak to its crossroad.
                corners.append(xyz[::-1] if preceding else xyz)
            others.append(other_end)
            coordinates.append((-1 if preceding else 1) * (1 if other_end[1] == 0 else -1))
        if others[0][0] != others[1][0] or others[0][1] == others[1][1] or coordinates[0] != coordinates[1]:
            raise ValueError('Bundle end cap requires one explicit connecting crossroad Section')
        middle = (rim_source or _road_rim)(context, others[0][0], coordinates[0])
        if others[0][1] == 1:
            middle = middle[::-1]
        a, b = corners[0], corners[1][::-1]
        return Cap(a, middle, b, True)

    @classmethod
    def build(cls, context, bundle):
        plan = context[Key('bundle_plan', bundle.id)]
        if plan.mode == 'union':
            return plan.separator(context)
        if len(bundle.members) != 2 or not all(m.complete for m in bundle.members):
            raise ValueError('Separated Bundle branch requires two complete ordered members')
        from .bundle_regions import interior_crossings
        interior_crossings(context.network, bundle)
        from .bundle_compounds import (CompoundBundleBuilder, compound_owners, has_internal_junction,
                                       has_consumed_sections)
        if has_consumed_sections(context, bundle) or compound_owners(context, bundle) or has_internal_junction(context.network, bundle):
            return CompoundBundleBuilder.build(context, bundle)
        members = tuple(Member.build(context, member, role) for role, member in enumerate(bundle.members))
        proof_checks = certify_disjoint(members[0].path, members[1].path)
        profiles = [context[Key('profile', t.section_id)] for m in bundle.members for t in m.path]
        if any(abs(p.datum - profiles[0].datum) > 1e-7 for p in profiles):
            raise ValueError('Bundle members require a common minimum height')
        caps = [cls.cap(context, *members, end) for end in (0, 1)]
        return cls.region(context, bundle, members, caps, proof_checks)

    @staticmethod
    def region(context, bundle, members, caps, proof_checks):
        profiles = [context[Key('profile', t.section_id)] for m in bundle.members for t in m.path]
        width = max(np.linalg.norm(members[1].rim[i, :2] - members[0].rim[i, :2]) for i in (0, -1))
        free_count = max(4, 2 * int(math.ceil(width / (2. * context.network.settings.spacing))))
        # Free rims have no neighboring mesh and can adopt the attached cap's
        # partition. Every pavement-adjacent rim keeps its original samples.
        for i in (0, 1):
            if caps[i] is None:
                count = len(caps[1 - i].rim) - 1 if caps[1 - i] is not None else free_count
                caps[i] = Cap.free(members[0].rim[0 if i == 0 else -1], members[1].rim[0 if i == 0 else -1], count)
        head, tail = caps

        def partition(cap, shift=0):
            if not cap.attached:
                return cap.first, cap.rim, cap.last
            a, b = len(cap.first) // 2 + shift, len(cap.last) // 2
            if not 0 <= a < len(cap.first) or not 0 <= b < len(cap.last):
                return None
            # Atlas corners lie on existing curved samples. The physical curb
            # is unchanged; a 180-degree mouth vertex is not a grid corner.
            return cap.first[:a + 1], _join(cap.first[a:], cap.middle, cap.last[:b + 1]), cap.last[b:]

        corner_shifts = None
        for head_shift, tail_shift in ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)):
            if (head_shift and not head.attached) or (tail_shift and not tail.attached):
                continue
            front_parts, back_parts = partition(head, head_shift), partition(tail, tail_shift)
            if front_parts is None or back_parts is None:
                continue
            ha, front, hb = front_parts
            ta, back, tb = back_parts
            left = _join(ha[::-1], members[0].rim, ta)
            right = _join(hb, members[1].rim, tb[::-1])
            rims = (left, back, right, front)
            counts = tuple(len(r) - 1 for r in rims)
            if (counts[0] - counts[2]) % 2 == (counts[1] - counts[3]) % 2 == 0:
                corner_shifts = head_shift, tail_shift
                break
        if corner_shifts is None:
            primitive_counts = (len(members[0].rim) - 1, len(members[1].rim) - 1,
                                len(head.rim) - 1, len(tail.rim) - 1)
            raise ValueError(f'Bundle separator boundary parity needs coordinated curb sampling: '
                             f'members/head/tail={primitive_counts}, total={sum(primitive_counts)}')
        pattern = context.patterns.quad(*counts)
        # u follows the Bundle; v goes left -> right in the common direction.
        u, v = pattern.uv.T
        samples = []
        for side, rim in enumerate(rims):
            ids = context.patterns.boundary(pattern, side)
            parameter = pattern.uv[ids, 0 if side in (0, 2) else 1]
            query = u if side in (0, 2) else v
            samples.append(np.column_stack([np.interp(query, parameter, rim[:, c]) for c in range(3)]))
        L, B, R, F = samples
        a, b, c, d = front[0], front[-1], back[-1], back[0]
        bilinear = ((1. - u) * (1. - v))[:, None] * a + ((1. - u) * v)[:, None] * b
        bilinear += (u * v)[:, None] * c + (u * (1. - v))[:, None] * d
        xyz = (1. - v[:, None]) * L + v[:, None] * R + (1. - u[:, None]) * F + u[:, None] * B - bilinear
        for side, rim in enumerate(rims):
            xyz[context.patterns.boundary(pattern, side)] = rim
        mesh = Mesh(xyz, pattern.quads[:, ::-1], 'grass')
        try:
            quality = mesh.validate()
            if quality['maximum_edge_ratio'] > 8.:
                raise ValueError(f"Bundle quad edge ratio {quality['maximum_edge_ratio']:.4g} exceeds 8")
        except ValueError as error:
            from .separator_surfaces import FixedRimSeparator
            loop = _join(left, back, right[::-1], front[::-1])[:-1][::-1]
            surface = FixedRimSeparator(context.patterns, spacing=context.network.settings.spacing).build(
                loop, datum=profiles[0].datum, height=profiles[0].rise * context.network.settings.exaggeration,
                exaggeration=context.network.settings.exaggeration)
            return Part(bundle.key, surface.mesh, dict(surface.metadata,
                initial_atlas_error=str(error), centerline_proof_checks=proof_checks, rim_intervals=counts,
                crown_ceiling=profiles[0].rise, attached_caps=[cap.attached for cap in caps]))
        boundary = np.concatenate(rims)
        xyz[:, 2], extension = plateau_extension(xyz[:, :2], boundary, profiles[0].datum,
            profiles[0].rise * context.network.settings.exaggeration)
        display_bound = extension.pop('continuous_slope_bound')
        extension['physical_continuous_slope_bound'] = (display_bound / context.network.settings.exaggeration
                                                       if context.network.settings.exaggeration else 0.)
        # Assign the saved boundary values exactly, including all corner keys.
        for side, rim in enumerate(rims):
            xyz[context.patterns.boundary(pattern, side)] = rim
        # Increasing longitudinal u then rightward v is clockwise in XY.
        mesh = Mesh(xyz, pattern.quads[:, ::-1], 'grass')
        mesh.seams = {side: context.patterns.boundary(pattern, side) for side in range(4)}
        quality = mesh.validate()
        return Part(bundle.key, mesh, dict(quality=quality, construction='separated boundary quad atlas',
            centerline_proof_checks=proof_checks, rim_intervals=counts,
            crown_ceiling=profiles[0].rise, attached_caps=[cap.attached for cap in caps],
            height_extension=extension,
            atlas_corner_shifts=corner_shifts,
            physical_slopes=mesh.slope_statistics(context.network.settings.exaggeration)))
