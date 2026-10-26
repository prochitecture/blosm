"""Mouth-aligned common-fibre unions retaining original rounded junction caps.

This branch restricts original polynomials, not their controls. Junction trims
may pass through explicit compatible Section seams. Skew mouths require another
chart; no radius or trim is changed to make them fit this one.
"""
from dataclasses import replace
import math
import numpy as np
from .bundle_axis import AxialPlan, AxialChartUnavailable, Boundary
from .mesh import Assembler
from .model import End, Key
from .providers import Part
from .updates import BuildContext


def member_ends(member):
    first, last = member.path[0], member.path[-1]
    return End(first.section_id, 0 if first.forward else 1), End(last.section_id, 1 if last.forward else 0)


def mouth_dependencies(network, bundle):
    if (len(bundle.members) != 2 or bundle.interior_section_ids or
            any(not m.complete for m in bundle.members)):
        return ()
    ends = {end for member in bundle.members for end in member_ends(member)}
    return tuple(sorted(j.layout_key for j in (*network.intersections, *network.connections)
                        if ends.intersection(j.ends)))


class PathRange:
    """A metric view preserving the source polynomial and persistent Section."""
    def __init__(self, path, first, last):
        if first < -1e-8 or last > path.length + 1e-8 or last <= first + 1e-7:
            raise AxialChartUnavailable('Attached Bundle has no positive visible member range')
        self.path, self.first = path, max(0., first)
        self.length = min(path.length, last) - self.first
        self.section, self.exact = path.section, path.exact
        self.parameters, self.stations = path.parameters, path.stations - self.first
        a, b = path.parameter(self.first), path.parameter(self.first + self.length)
        self.knots = path.knots[(path.knots > a + 1e-10) & (path.knots < b - 1e-10)]
        self.points = np.vstack((path.frame(self.first)[0], path.jet(self.knots)[0],
                                 path.frame(self.first + self.length)[0]))
        if not self.exact:
            # Certifying the complete source is conservative for this restriction.
            self.coefficients = path.coefficients

    def frame(self, station):
        return self.path.frame(np.asarray(station) + self.first)

    def point(self, station, distance):
        point, _, normal = self.frame(station)
        return point + np.asarray(distance)[..., None] * normal

    def samples(self, first, last, spacing, **kwargs):
        return self.path.samples(first + self.first, last + self.first, spacing, **kwargs) - self.first


class AttachedAxialPlan(AxialPlan):
    @classmethod
    def build(cls, context, bundle):
        if not mouth_dependencies(context.network, bundle):
            return None
        plan = cls(context, bundle)
        if plan.minimum_gap > 1e-8 and not plan.fill_micro_gap and not plan.physical_roles_swapped:
            # Keep the inexpensive established rectangular separator for
            # straight separated members. Curved members retain this chart
            # across the overlap event, independently of edit history.
            straight = []
            for points in plan.visible_points:
                delta = points[-1] - points[0]
                offsets = points - points[0]
                error = abs(offsets[:, 0] * delta[1] - offsets[:, 1] * delta[0]) / np.linalg.norm(delta)
                straight.append(bool(np.max(error) <= 1e-8))
            if all(straight):
                return None
        if plan.physical_roles_swapped or np.min(plan.gap(np.array((plan.first, plan.last)))) <= 1e-7:
            raise AxialChartUnavailable('Attached overlap needs separated, consistently ordered mouth profiles')
        plan._prefer_fibre_strips(context)
        return plan

    def __init__(self, context, bundle):
        from .bundles import CenterlineRange, certify_disjoint
        from .intersections import half_intervals
        sites = {e: j for j in (*context.network.intersections, *context.network.connections) for e in j.ends}
        products, self.ranges, intervals = dict(context.products), {}, set()
        visible, visible_members, points = [], [], []
        for member in bundle.members:
            for previous, following in zip(member.path, member.path[1:]):
                exit = End(previous.section_id, 1 if previous.forward else 0)
                entry = End(following.section_id, 0 if following.forward else 1)
                seam = sites.get(exit)
                if seam is None or entry not in seam.ends:
                    raise ValueError('Attached Bundle chains need explicit two-ended Section seams')
                if len(seam.ends) != 2:
                    raise AxialChartUnavailable('Attached Bundle member crosses an internal physical junction')
            ends = member_ends(member)
            layouts = []
            for endpoint in ends:
                site = sites.get(endpoint)
                if site is None:
                    layouts.append(None)
                    continue
                if len(site.ends) < 3:
                    raise AxialChartUnavailable('Attached union needs a free end or a multi-arm mouth')
                layouts.append(context[site.layout_key])
            active, ranges, controls, indices = [], [], [], []
            for index, traversal in enumerate(member.path):
                sid = traversal.section_id
                path = context[Key('path', sid)]
                first, last = 0., path.length
                for end, layout in enumerate(layouts):
                    if layout is None:
                        continue
                    local_end = (0 if traversal.forward else 1) if end == 0 else (1 if traversal.forward else 0)
                    mouth = layout.mouths.get((sid, local_end))
                    if mouth is None:
                        continue
                    intervals.add(mouth.transverse_intervals)
                    if local_end == 0:
                        first = max(first, mouth.station)
                    else:
                        last = min(last, mouth.station)
                self.ranges[sid] = first, last
                if last <= first + 1e-7:
                    continue  # This original Section remains junction-owned.
                view = PathRange(path, first, last) if first > 1e-8 or last < path.length - 1e-8 else path
                active.append(traversal)
                indices.append(index)
                ranges.append(CenterlineRange.build(path, first, last))
                controls.append(view.points if traversal.forward else view.points[::-1])
                products[Key('path', sid)] = view
            if not active or indices != list(range(indices[0], indices[-1] + 1)):
                raise AxialChartUnavailable('Attached Bundle needs one contiguous visible range per member')
            visible_members.append(replace(member, path=tuple(active)))
            visible.append(tuple(ranges))
            points.append(np.concatenate(controls))
        self.proof_checks = certify_disjoint(*visible)
        if len(intervals) > 1:
            raise AxialChartUnavailable('Attached mouths need one compatible transverse partition')
        self.port_intervals = next(iter(intervals), 2 * half_intervals(context.network.settings))
        self._pavement_part = None
        self.strip_selection = None
        self.visible_members, self.visible_points = tuple(visible_members), tuple(points)
        visible_bundle = replace(bundle, members=self.visible_members)
        super().__init__(BuildContext(context.network, products, context.patterns), visible_bundle, allow_separated=True)
        self.bundle = bundle  # Derived ranges never replace authored membership.

    def _prefer_fibre_strips(self, context):
        if not self.islands or self.separator_correction is not None:
            return
        from .fibre_strips import CappedFibreStrips
        previous = self.peak_a, self.peak_b, self._separator_cache, self.separator_correction
        try:
            candidate = CappedFibreStrips().build(context, self)
            if len(candidate.mesh.quads) >= .8 * len(self._separator_cache.quads):
                return
            for name, extras in (('peak_a', candidate.extra_first), ('peak_b', candidate.extra_second)):
                boundary = getattr(self, name)
                points = np.vstack((boundary.points, [p for p, key in extras])) if extras else boundary.points
                keys = boundary.keys + [key for p, key in extras]
                order = np.argsort((points - self.origin) @ self.axis, kind='stable')
                setattr(self, name, Boundary(points[order], [keys[i] for i in order]))
            self._separator_cache, self.separator_correction = candidate.mesh, candidate.metadata
            # A denser curb partition also changes pavement. Select it only if
            # both surfaces pass their unchanged quality and physical guards.
            self._pavement_part = super().pavement(context)
            self.strip_selection = dict(selected=True, previous_quads=len(previous[2].quads),
                                        selected_quads=len(candidate.mesh.quads))
        except ValueError as exc:
            self.peak_a, self.peak_b, self._separator_cache, self.separator_correction = previous
            self._pavement_part = None
            self.strip_selection = dict(selected=False, reason=str(exc))

    def pavement(self, context):
        if self._pavement_part is None:
            self._pavement_part = super().pavement(context)
        self._pavement_part.metadata.update(attached_mouth_ranges=self.ranges,
            mouth_transverse_intervals=self.port_intervals, separator_strip_selection=self.strip_selection,
            visible_member_sections=[[t.section_id for t in member.path] for member in self.visible_members])
        return self._pavement_part

    @staticmethod
    def seam_keys(mesh, joins):
        keys = [None] * len(mesh.xyz)
        for end, join in enumerate(joins):
            if join is None:
                continue
            a, b, count = join
            delta = b - a
            t = (mesh.xyz - a) @ delta / float(delta @ delta)
            distance = np.linalg.norm(mesh.xyz - (a + t[:, None] * delta), axis=1)
            found = set()
            for i in np.flatnonzero((distance < 1e-8) & (t >= -1e-8) & (t <= 1. + 1e-8)):
                k = int(round(float(t[i]) * count))
                if np.linalg.norm(mesh.xyz[i] - (a + k / count * delta)) > 1e-8:
                    raise ValueError('Attached cap does not retain the declared body cross-section')
                keys[i] = ('attached cap/body', end, k)
                found.add(k)
            if found and found != set(range(count + 1)):
                raise ValueError('Attached cap/body boundary has incompatible edge subdivisions')
        return keys

    def separator(self, context):
        from .bundles import BundleProvider, Member, Cap
        body = super().separator(context)
        if body is None:
            raise ValueError('Attached overlap has no positive end separator islands')
        members = tuple(Member.build(context, member, role) for role, member in enumerate(self.visible_members))
        caps, joins, attached = [], [], []
        for end in (0, 1):
            cap = BundleProvider.cap(context, *members, end)
            attached.append(cap is not None)
            island = self.islands[0 if end == 0 else -1]
            across = max(4, 4 * int(math.ceil(np.linalg.norm(island.first.points - island.second.points,
                                                           axis=1).max() / (4. * self.settings.spacing))))
            pieces = tuple(replace(m, rim=m.rim[[0 if end == 0 else -1]]) for m in members)
            flat = Cap.free(pieces[0].rim[0], pieces[1].rim[0], across)
            joins.append((flat.rim[0], flat.rim[-1], across) if cap is not None else None)
            if cap is not None:
                caps.append(BundleProvider.region(context, self.bundle, pieces,
                    [cap, flat] if end == 0 else [flat, cap], self.proof_checks))
        assembly = Assembler()
        assembly.add(body.mesh.xyz, body.mesh.quads, self.seam_keys(body.mesh, joins))
        for cap in caps:
            assembly.add(cap.mesh.xyz, cap.mesh.quads, self.seam_keys(cap.mesh, joins))
        mesh = assembly.mesh('grass')
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        if counts.max() > 2 or len(mesh.xyz) - len(unique) + len(mesh.quads) != len(self.islands):
            raise ValueError('Attached separator caps and body do not form the declared islands')
        return Part(self.bundle.key, mesh, dict(quality=mesh.validate(), islands=len(self.islands),
            construction='mouth-aligned Bundle union with original rounded caps', attached_mouth_ranges=self.ranges,
            visible_member_sections=[[t.section_id for t in member.path] for member in self.visible_members],
            attached_caps=attached,
            welded_cap_body_vertices=len(assembly.lookup), caps=[p.metadata for p in caps],
            separator_strip_selection=self.strip_selection, physical_slopes=mesh.slope_statistics(self.settings.exaggeration)))

    def release_geometry(self):
        super().release_geometry()
        self._pavement_part = None
