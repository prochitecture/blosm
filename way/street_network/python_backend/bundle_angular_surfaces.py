"""Finite physical road caps and crown-aligned direct quad surface patches."""
import math
import json
import copy
from types import SimpleNamespace
import numpy as np
from .mesh import Mesh, Assembler, cross
from .regions import signed_area
from .bundle_axis import _simplify
from .bundle_angular_curves import RayCurb, polyline_distance
from .four_side_charts import FourSidedQuadChart, quality as chart_quality
from .bundle_terminal_caps import seed_terminal_chart, graded_segment
from .quality_correction import QuadShapeCorrection

class AngularSurfaceBuilder:
    """Reuse bounded patterns with explicit curb, cap, crown and island interfaces."""

    def __init__(self, patterns, *, spacing, curve_tolerance, maximum_slope=0.1,
                 maximum_longitudinal_slope=0.03, max_cells=50000, cap_scale=None, mandatory_points=None):
        self.patterns = patterns
        (self.spacing, self.curve_tolerance) = (spacing, curve_tolerance)
        self.maximum_slope = maximum_slope
        self.maximum_longitudinal_slope = maximum_longitudinal_slope
        self.max_cells = max_cells
        self.cap_scale = spacing if cap_scale is None else cap_scale
        self.mandatory_points = mandatory_points
        self.boundary_sampling = []
        self.extended_grading = False

    @staticmethod
    def _longitudinal_bound(frame, mesh):
        return frame.longitudinal_bound(mesh)

    def split_value(self, field, xy):
        uv = field.frame.local(np.asarray(xy))
        return uv[..., 1] - self.crest(field, uv[..., 0])

    def split_outer(self, xy, labels, field):
        """Insert the two crown/end-cut intersections without changing a segment."""
        values = self.split_value(field, xy)
        (points, kinds, root_ids) = ([], [], [])
        for (i, (a, b)) in enumerate(zip(xy, np.roll(xy, -1, axis=0))):
            points.append(a)
            kinds.append(labels[i])
            (first, last) = (values[i], values[(i + 1) % len(xy)])
            if abs(first) < 1e-10:
                root_ids.append(len(points) - 1)
            if first * last < 0.0:
                (lo, hi) = (0.0, 1.0)
                for _ in range(48):
                    middle = 0.5 * (lo + hi)
                    value = float(self.split_value(field, a + middle * (b - a)))
                    if value * first > 0.0:
                        lo = middle
                    else:
                        hi = middle
                root_ids.append(len(points))
                points.append(a + 0.5 * (lo + hi) * (b - a))
                kinds.append(labels[i])
        if len(root_ids) != 2:
            raise ValueError(f'The physical outer cap chain has {len(root_ids)} crown crossings')
        points = np.asarray(points)
        root_ids = sorted(root_ids, key=lambda i: field.frame.local(points[i:i + 1])[0, 0])
        (first, last) = root_ids
        indices = np.array([(first + i) % len(points) for i in range((last - first) % len(points) + 1)])
        other = np.array([(last + i) % len(points) for i in range((first - last) % len(points) + 1)])[::-1]
        chains = [(points[indices], [kinds[i] for i in indices[:-1]]), (points[other], [kinds[i] for i in other[1:]])]
        chains.sort(key=lambda p: -float(np.mean(self.split_value(field, p[0][1:-1]))))
        return chains

    def simplify_sample(self, points, labels, tolerance, prefix):
        """Approximate only between protected semantic boundary changes/corners."""
        points = np.asarray(points)
        if self.mandatory_points is not None:
            from .polyline_sampling import ProtectedPolylineSampler, landmark_indices
            protected = landmark_indices(points,self.mandatory_points)
            protected.update(i for i in range(1,len(labels)) if labels[i]!=labels[i-1])
            sampled,owners,report = ProtectedPolylineSampler(spacing=self.spacing,
                tolerance=.6*self.curve_tolerance).sample(points,protected)
            self.boundary_sampling.append(dict(prefix=prefix,**report))
            return sampled,[labels[i] for i in owners],[(prefix,i) for i in range(len(sampled))]
        edge = np.diff(points, axis=0)
        lengths = np.linalg.norm(edge, axis=1)
        directions = edge / lengths[:, None]
        turns = abs(np.arctan2(cross(directions[:-1], directions[1:]), np.sum(directions[:-1] * directions[1:], axis=1)))
        protected = {0, len(points) - 1}
        protected.update((i + 1 for (i, t) in enumerate(turns) if t > 0.05 and (tolerance < 1e-08 or min(lengths[i], lengths[i + 1]) > self.curve_tolerance)))
        protected.update((i for i in range(1, len(labels)) if labels[i] != labels[i - 1]))
        simplified = _simplify(points, tolerance, protected)
        lookup = {tuple(p): i for (i, p) in enumerate(points)}
        original = [lookup[tuple(p)] for p in simplified]
        (result, kinds) = ([], [])
        sizes = np.linalg.norm(np.diff(simplified, axis=0), axis=1)
        for segment, (a, b, ia, ib) in enumerate(zip(simplified, simplified[1:], original, original[1:])):
            assert len(set(labels[ia:ib])) == 1
            length = sizes[segment]
            first = min(self.spacing, 16.*sizes[segment-1]) if segment else self.spacing
            last = min(self.spacing, 16.*sizes[segment+1]) if segment+1<len(sizes) else self.spacing
            if min(first, last) < self.spacing:
                fractions = graded_segment(length, self.spacing, first, last)
            else:
                fractions = np.linspace(0., 1., max(1, int(math.ceil(length/self.spacing)))+1)
            count = len(fractions)-1
            result.extend((a+fractions[:-1, None]*(b-a)))
            kinds.extend([labels[ia]] * count)
        result.append(simplified[-1])
        return (np.asarray(result), kinds, [(prefix, i) for i in range(len(result))])

    def ridge(self, frame, field, first, last, prefix):
        s = frame.local(np.asarray((first, last)))[:, 0]
        count = max(1, int(math.ceil((s[1] - s[0]) / self.spacing)))
        if count > 4096:
            raise ValueError('Ridge exceeds the bounded sample count')
        station = np.linspace(*s, count + 1)
        xy = frame.world(np.c_[station, self.crest(field, station)])
        xy[[0, -1]] = (first, last)
        return (xy, [(prefix, i) for i in range(len(xy))])

    def island_halves(self, loop, frame):
        xy = loop.xy
        station = frame.local(xy)[:, 0]
        (first, last) = (int(station.argmin()), int(station.argmax()))
        a = np.array([(first + i) % len(xy) for i in range((last - first) % len(xy) + 1)])
        b = np.array([(last + i) % len(xy) for i in range((first - last) % len(xy) + 1)])[::-1]
        chains = [xy[a], xy[b]]
        chains.sort(key=lambda p: -float(frame.local(p)[:, 1].mean()))
        return chains

    def join_chains(self, chains):
        (points, keys) = ([], [])
        for (xy, labels) in chains:
            if points:
                if np.linalg.norm(points[-1] - xy[0]) > 1e-08:
                    raise ValueError('Semantic patch chains disagree at an endpoint')
                labels = list(labels)
                labels[0] = keys[-1]
                points.extend(xy[1:])
                keys.extend(labels[1:])
            else:
                points.extend(xy)
                keys.extend(labels)
        return (np.asarray(points), keys)

    def _direct_single(self, frame, field, outer, seam, name, *, grass=False):
        try:
            return self._mapped_single(frame,field,outer,seam,name,grass=grass)
        except ValueError as original_error:
            if grass or not self.extended_grading or not name.endswith('/body'):
                raise
            from .bundle_angular_road_tiles import AngularRoadTiles
            try:
                return AngularRoadTiles(self,frame,field).build(outer,seam,name)
            except ValueError as error:
                raise ValueError(str(original_error)+'; local body: '+str(error)) from error

    def _mapped_single(self, frame, field, outer, seam, name, *, grass=False):
        (physical, kinds, boundary_keys) = outer
        (shared, shared_keys) = seam
        if np.linalg.norm(physical[0] - shared[0]) > 1e-08 or np.linalg.norm(physical[-1] - shared[-1]) > 1e-08:
            raise ValueError('The patch has unmatched crown endpoints')
        boundary_keys = list(boundary_keys)
        (boundary_keys[0], boundary_keys[-1]) = (shared_keys[0], shared_keys[-1])
        semantic = None
        if not grass:
            curb_edges = [i for (i, kind) in enumerate(kinds) if kind[0] == 'outer']
            (first, last) = (min(curb_edges), max(curb_edges) + 1)
            if any((kind[0] != 'outer' for kind in kinds[first:last])):
                raise ValueError('A semantic half patch has a disconnected outer curb')

            def insert(at, address):
                nonlocal physical, kinds
                physical = np.insert(physical, at + 1, 0.5 * (physical[at] + physical[at + 1]), axis=0)
                kinds = kinds[:at + 1] + [kinds[at]] + kinds[at + 1:]
                boundary_keys.insert(at + 1, (name, 'semantic parity', address))
            if (last - first - (len(shared) - 1)) % 2:
                at = first + int(np.argmax(np.linalg.norm(np.diff(physical[first:last + 1], axis=0), axis=1)))
                insert(at, 'curb')
                last += 1
            if (first - (len(physical) - 1 - last)) % 2:
                at = int(np.argmax(np.linalg.norm(np.diff(physical[:first + 1], axis=0), axis=1)))
                insert(at, 'cap')
                first += 1
                last += 1
            semantic = (0, first, last, len(physical) - 1)
        if (len(physical) + len(shared) - 2) % 2:
            at = int(np.argmax(np.linalg.norm(np.diff(physical, axis=0), axis=1)))
            physical = np.insert(physical, at + 1, 0.5 * (physical[at] + physical[at + 1]), axis=0)
            kinds = kinds[:at + 1] + [kinds[at]] + kinds[at + 1:]
            boundary_keys.insert(at + 1, (name, 'parity', at))
        rim = np.vstack((physical, shared[-2:0:-1]))
        keys = boundary_keys + shared_keys[-2:0:-1]
        values = field.values(rim)
        curb_ids = set()
        for (i, kind) in enumerate(kinds):
            if kind[0] == 'outer':
                curb_ids.update((i, i + 1))
        values[list(curb_ids)] = field.datum
        values[len(physical):] = field.ceiling
        values[0] = values[len(physical) - 1] = field.ceiling
        if grass:
            values[:] = field.ceiling
        if signed_area(rim) < 0.0:
            (rim, keys, values) = (rim[::-1].copy(), keys[::-1], values[::-1].copy())
            if semantic is not None:
                semantic = tuple(sorted((len(rim) - 1 - i for i in semantic)))
        uv = frame.local(rim)
        builder = FourSidedQuadChart(self.patterns, max_candidates=24, quality_iterations=100, maximum_ratio=8.0, max_cells=self.max_cells)
        semantic_failure = None
        chart = None
        if semantic is not None:
            try:
                (mapped, counts) = builder.map(uv, semantic)
                (mapped, seed) = seed_terminal_chart(mapped, len(uv))
                measured = chart_quality(mapped)
                correction = None
                if not builder.accepted(measured):
                    (mapped, correction) = QuadShapeCorrection(max_cells=self.max_cells, max_iterations=250).improve(mapped)
                    measured = chart_quality(mapped)
                if not builder.accepted(measured):
                    raise ValueError(f'Semantic initial quality {measured}')
                chart = SimpleNamespace(mesh=mapped, metadata=dict(quality=measured, atlas_corners=semantic, rim_intervals=counts, construction='semantic curb/cap/crown initial direct quad pattern', candidates=1, improvement=correction, terminal_seed=seed))
            except ValueError as error:
                semantic_failure = str(error)
        if chart is None:
            try:
                chart = builder.build(uv)
            except ValueError as error:
                raise ValueError(f'Semantic: {semantic_failure}; generic: {error}') from error
        chart.metadata['semantic_first_attempt_failure'] = semantic_failure
        xyz = np.c_[frame.world(chart.mesh.xyz[:, :2]), np.zeros(len(chart.mesh.xyz))]
        xyz[:len(rim), :2] = rim
        xyz[:, 2] = field.ceiling if grass else field.values(xyz[:, :2])
        xyz[:len(rim), 2] = values
        keys.extend(((name, 'inside', i) for i in range(len(rim), len(xyz))))
        mesh = Mesh(xyz, chart.mesh.quads, 'grass' if grass else 'pavement')
        quality = mesh.validate()
        world_correction = None
        if self.extended_grading and quality['maximum_edge_ratio']>8.:
            from .surface_quality import ProfiledQuadCorrection
            mesh,world_correction = ProfiledQuadCorrection(max_cells=self.max_cells,max_iterations=500).improve(
                mesh,(lambda xy:np.full(len(xy),field.ceiling)) if grass else field.values,
                datum=field.datum,ceiling=field.ceiling,exaggeration=1.,maximum_slope=self.maximum_slope)
            quality = mesh.validate()
        if quality['maximum_edge_ratio'] > 8.0:
            raise ValueError(f"World quad edge ratio {quality['maximum_edge_ratio']} exceeds 8")
        slope_correction = None
        if self.extended_grading and not grass:
            from .surface_quality import ProfiledSlopeCorrection
            mesh, slope_correction = ProfiledSlopeCorrection(max_cells=self.max_cells).improve(
                mesh, field.values, datum=field.datum, ceiling=field.ceiling,
                maximum_slope=self.maximum_slope, longitudinal_bound=frame.longitudinal_bound,
                maximum_longitudinal_slope=self.maximum_longitudinal_slope)
            quality = mesh.validate()
        details = dict(chart=chart.metadata, world_quality=quality, boundary_vertices=len(rim), boundary_accuracy_budget=.7*self.curve_tolerance, provenance_labels=True)
        if world_correction is not None:
            details['world_correction'] = world_correction
        if slope_correction is not None and slope_correction['applied']:
            details['profile_slope_correction'] = slope_correction
        return mesh,keys,details

    def bounded_subchain(self, ray, keys, first, last, first_key, last_key):
        middle = (ray.stations > first + 1e-08) & (ray.stations < last - 1e-08)
        points = np.vstack((ray(first), ray.points[middle], ray(last)))
        names = [first_key] + [key for (key, hit) in zip(keys, middle) if hit] + [last_key]
        for (slot, station) in ((0, first), (-1, last)):
            at = int(np.argmin(abs(ray.stations - station)))
            if abs(ray.stations[at] - station) < 1e-08:
                points[slot] = ray.points[at]
                names[slot] = keys[at]
        return (points, names)

    def balanced_physical(self, outer, seam, name):
        """Fix atlas parity only on unshared source edges, before any patch build."""
        (physical, kinds, keys) = outer
        (physical, kinds, keys) = (physical.copy(), list(kinds), list(keys))
        curb = [i for (i, kind) in enumerate(kinds) if kind[0] == 'outer']
        (first, last) = (min(curb), max(curb) + 1)

        def insert(at, role):
            nonlocal physical, kinds
            physical = np.insert(physical, at + 1, 0.5 * (physical[at] + physical[at + 1]), axis=0)
            kinds = kinds[:at + 1] + [kinds[at]] + kinds[at + 1:]
            keys.insert(at + 1, (name, 'source parity', role))
        if (last - first - (len(seam[0]) - 1)) % 2:
            at = first + int(np.argmax(np.linalg.norm(np.diff(physical[first:last + 1], axis=0), axis=1)))
            insert(at, 'curb')
            last += 1
        if (first - (len(physical) - 1 - last)) % 2:
            allowed = [i for i in range(len(kinds)) if (i < first or i >= last) and kinds[i][0] != 'cut']
            if not allowed:
                raise ValueError('A computational cut requires coordinated parity, not a unilateral extra point')
            at = max(allowed, key=lambda i: float(np.linalg.norm(physical[i + 1] - physical[i])))
            insert(at, 'terminal')
        return (physical, kinds, keys)

    def finite_half(self, frame, field, outer, seam, name, margin):
        (physical, kinds, keys) = outer
        (shared, shared_keys) = seam
        curb = [i for (i, kind) in enumerate(kinds) if kind[0] == 'outer']
        (f, l) = (min(curb), max(curb) + 1)
        kind = kinds[f]
        if any((k != kind for k in kinds[f:l])):
            raise ValueError('A finite road half has more than one outer-curb owner')
        curb_ray = RayCurb(frame, physical[f:l + 1])
        seam_ray = RayCurb(frame, shared)
        a = max(curb_ray.stations[0], seam_ray.stations[0]) + margin
        b = min(curb_ray.stations[-1], seam_ray.stations[-1]) - margin
        available = np.flatnonzero((seam_ray.stations >= a) & (seam_ray.stations <= b))
        if len(available) < 3:
            raise ValueError('A finite half has too few existing shared crown samples')
        (a, b) = seam_ray.stations[available[[0, -1]]]
        if b <= a + 4.0 * self.spacing:
            raise ValueError('Finite road caps leave no positive body')
        (curb_head, ck_head) = self.bounded_subchain(curb_ray, keys[f:l + 1], curb_ray.stations[0], a, keys[f], (name, 'head outer'))
        (curb_body, ck_body) = self.bounded_subchain(curb_ray, keys[f:l + 1], a, b, (name, 'head outer'), (name, 'tail outer'))
        (curb_tail, ck_tail) = self.bounded_subchain(curb_ray, keys[f:l + 1], b, curb_ray.stations[-1], (name, 'tail outer'), keys[l])
        (seam_head, sk_head) = self.bounded_subchain(seam_ray, shared_keys, seam_ray.stations[0], a, shared_keys[0], (name, 'head crown'))
        (seam_body, sk_body) = self.bounded_subchain(seam_ray, shared_keys, a, b, (name, 'head crown'), (name, 'tail crown'))
        (seam_tail, sk_tail) = self.bounded_subchain(seam_ray, shared_keys, b, seam_ray.stations[-1], (name, 'tail crown'), shared_keys[-1])
        (cuts, cut_keys) = ([], [])
        for (role, c, s, ckey, skey) in (('head', curb_body[0], seam_body[0], ck_body[0], sk_body[0]), ('tail', curb_body[-1], seam_body[-1], ck_body[-1], sk_body[-1])):
            count = max(2, int(math.ceil(np.linalg.norm(s - c) / self.spacing)))
            if self.extended_grading:
                if role=='head':
                    adjacent = (curb_head[-1]-curb_head[-2],curb_body[1]-curb_body[0],
                                seam_head[-1]-seam_head[-2],seam_body[1]-seam_body[0])
                else:
                    adjacent = (curb_body[-1]-curb_body[-2],curb_tail[1]-curb_tail[0],
                                seam_body[-1]-seam_body[-2],seam_tail[1]-seam_tail[0])
                minimum = min(np.linalg.norm(edge) for edge in adjacent)
                count = max(count,int(math.ceil(np.linalg.norm(s-c)/min(self.spacing,6.*minimum))))
                if count>64:
                    raise ValueError('Road cap cut exceeds its bounded profile family')
            count += count % 2
            points = np.linspace(c, s, count + 1)
            names = [(name, role, 'cut', i) for i in range(count + 1)]
            (names[0], names[-1]) = (ckey, skey)
            cuts.append(points)
            cut_keys.append(names)
        source_head = np.vstack((physical[:f], curb_head))
        source_head_keys = keys[:f] + ck_head
        source_head_kinds = kinds[:f] + [kind] * (len(curb_head) - 1)
        source_tail = np.vstack((curb_tail, physical[l + 1:]))
        source_tail_keys = ck_tail + keys[l + 1:]
        source_tail_kinds = [kind] * (len(curb_tail) - 1) + kinds[l:]
        head = (np.vstack((source_head, cuts[0][1:])), source_head_kinds + [('cut', 'head')] * (len(cuts[0]) - 1), source_head_keys + cut_keys[0][1:])
        body = (np.vstack((cuts[0][::-1], curb_body[1:], cuts[1][1:])), [('cut', 'head')] * (len(cuts[0]) - 1) + [kind] * (len(curb_body) - 1) + [('cut', 'tail')] * (len(cuts[1]) - 1), cut_keys[0][::-1] + ck_body[1:] + cut_keys[1][1:])
        tail = (np.vstack((cuts[1][::-1], source_tail[1:])), [('cut', 'tail')] * (len(cuts[1]) - 1) + source_tail_kinds, cut_keys[1][::-1] + source_tail_keys[1:])
        (assembly, details) = (Assembler(), [])
        for (role, rim, shared_part) in (('head', head, (seam_head, sk_head)), ('body', body, (seam_body, sk_body)), ('tail', tail, (seam_tail, sk_tail))):
            patch_name = name + '/' + role
            candidates = self.terminal_variants(rim, patch_name) if role != 'body' else (rim,)
            failed = []
            for candidate in candidates:
                candidate = self.balanced_physical(candidate, shared_part, patch_name)
                try:
                    mesh, addresses, report = self._direct_single(frame, field, candidate, shared_part, patch_name)
                    break
                except ValueError as error:
                    failed.append(str(error))
            else:
                raise ValueError(f'{role}: '+ '; '.join(dict.fromkeys(failed)))
            report['terminal_partition_failures'] = failed
            if any(('semantic parity' in key or 'parity' in key for key in addresses if isinstance(key, tuple))):
                raise ValueError('A finite half modified a shared computational cut partition')
            assembly.add(mesh.xyz, mesh.quads, addresses)
            details.append(report)
        mesh = assembly.mesh()
        mesh.validate()
        base_keys = [None] * len(mesh.xyz)
        for (key, index) in assembly.lookup.items():
            base_keys[index] = key
        return (mesh, base_keys, dict(construction='finite physical road caps and crown-aligned body', margin=margin, patches=details, world_quality=mesh.validate()))

    @staticmethod
    def terminal_variants(rim, name):
        """Bounded direct boundary partitions, without changing any source edge.

        Additional pairs of samples change local 1:3 pattern placement. Shared
        computational cuts and original curb/cap corners remain authoritative.
        """
        yield rim
        points, kinds, keys = rim
        front = [i for i, kind in enumerate(kinds) if kind[0] not in ('outer', 'cut')]
        chosen = sorted(front, key=lambda i:-float(np.linalg.norm(points[i+1]-points[i])))[:2]
        for additions in ((2, 0), (0, 2), (4, 0), (0, 4), (2, 2)):
            counts = dict(zip(chosen, additions))
            if not any(counts.values()):
                continue
            xy, labels, addresses = [], [], []
            for i, (a, b) in enumerate(zip(points, points[1:])):
                extra = counts.get(i, 0)
                samples = np.linspace(a, b, extra+2)[:-1]
                xy.extend(samples)
                labels.extend([kinds[i]]*(extra+1))
                addresses.append(keys[i])
                addresses.extend((name, 'terminal partition', i, j, extra) for j in range(extra))
            xy.append(points[-1])
            addresses.append(keys[-1])
            yield np.asarray(xy), labels, addresses

    def direct_patch(self, frame, field, outer, seam, name, *, grass=False):
        try:
            return self._patch(frame,field,outer,seam,name,grass=grass)
        except ValueError as original_error:
            if self.extended_grading:
                raise
            candidate = copy.copy(self)
            candidate.extended_grading = True
            mesh,keys,report = candidate._patch(frame,field,outer,seam,name,grass=grass)
            report['extended_physical_grading'] = True
            report['established_construction_failure'] = str(original_error)[:2000]
            return mesh,keys,report

    def _patch(self, frame, field, outer, seam, name, *, grass=False):
        if grass:
            return self._direct_single(frame, field, outer, seam, name, grass=True)
        failures = []
        for margin in (self.cap_scale, 2.*self.cap_scale, 3.*self.cap_scale,
                       4.*self.cap_scale, 6.*self.cap_scale):
            try:
                (mesh, keys, details) = self.finite_half(frame, field, outer, seam, name, margin)
                maximum = float(mesh.quad_slopes().max())
                if maximum > self.maximum_slope + 1e-08:
                    raise ValueError(f'Finite-half physical slope exceeds the physical slope limit: {maximum}')
                longitudinal = self._longitudinal_bound(frame, mesh)
                if longitudinal > self.maximum_longitudinal_slope + 1e-08:
                    raise ValueError(f'Finite-half angular longitudinal bound exceeds the longitudinal slope limit: {longitudinal}')
                details['angular_longitudinal_bound'] = longitudinal
                details['failures'] = failures
                return (mesh, keys, details)
            except ValueError as error:
                failures.append(dict(margin=margin, error=str(error)))
        try:
            (mesh, keys, details) = self._direct_single(frame, field, outer, seam, name)
            maximum = float(mesh.quad_slopes().max())
            if maximum > self.maximum_slope + 1e-08:
                raise ValueError(f'Whole-half physical slope exceeds the physical slope limit: {maximum}')
            longitudinal = self._longitudinal_bound(frame, mesh)
            if longitudinal > self.maximum_longitudinal_slope + 1e-08:
                raise ValueError(f'Whole-half angular longitudinal bound exceeds the longitudinal slope limit: {longitudinal}')
            details['angular_longitudinal_bound'] = longitudinal
            details['finite_failures'] = failures
            return (mesh, keys, details)
        except ValueError as error:
            failures.append(dict(branch='whole half', error=str(error)))
        raise ValueError('No qualified finite road half: ' + json.dumps(failures))

    def boundary_edges(self, mesh):
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        (edges, counts) = np.unique(edges, axis=0, return_counts=True)
        if counts.max() > 2:
            raise ValueError('The assembled surface contains a nonmanifold edge')
        return edges[counts == 1]

    def boundary_contract(self, pavement, high, low, grass=None):
        """Exclude internal cracks before treating physical boundaries as outlets."""
        edges = self.boundary_edges(pavement)
        ids = np.unique(edges)
        degrees = np.bincount(edges.ravel(), minlength=len(pavement.xyz))
        if np.any(degrees[ids] != 2):
            raise ValueError('The physical pavement boundary has an unresolved branch or T-junction')
        references = [high, low]
        conformity = None
        if grass is not None:
            ge = self.boundary_edges(grass)
            gi = np.unique(ge)
            mapping = {}
            for start in range(0, len(gi), 64):
                batch = gi[start:start + 64]
                distance = np.linalg.norm(grass.xyz[batch, None, :2] - pavement.xyz[ids, :2], axis=-1)
                nearest = np.argmin(distance, axis=1)
                if np.max(distance[np.arange(len(batch)), nearest]) > 1e-08:
                    raise ValueError('The pavement lost a shared grass boundary vertex')
                mapping.update(((int(g), int(p)) for (g, p) in zip(batch, ids[nearest])))
            available = {tuple(e) for e in edges}
            if any((tuple(sorted((mapping[int(a)], mapping[int(b)]))) not in available for (a, b) in ge)):
                raise ValueError('A grass boundary edge is split differently in the pavement')
            conformity = dict(shared_vertices=len(gi), shared_edges=len(ge), whole_edges_exact=True)
            references.extend((grass.xyz[e, :2] for e in ge))
        queries = np.vstack((pavement.xyz[ids, :2], pavement.xyz[edges, :2].mean(axis=1)))
        nearest = np.full(len(queries), np.inf)
        for chain in references:
            nearest = np.minimum(nearest, polyline_distance(queries, np.asarray(chain)))
        if nearest.max() > 1e-07:
            raise ValueError(f'An internal crack was exposed as a boundary: {nearest.max()} m')
        return (ids, dict(vertices=len(ids), edges=len(edges), degree_two=True, maximum_rim_distance=float(nearest.max()), grass_conformity=conformity))

    def crest(self, field, station):
        (first, last) = field.peaks(station)
        return 0.5 * (first + last)
