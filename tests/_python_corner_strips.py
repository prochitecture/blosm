"""Experimental normal corner strips with finite, profiled core mouths.

Not registered in the runtime. Physical curbs/complete incoming profiles remain
fixed while a finite central profile interval prevents crown-guide core cusps.
"""
from dataclasses import dataclass
import math
import numpy as np
from tested_python_streets.conditions import BoundaryData, CurbEdge, ProfileEdge, Port, _profile_gradient
from tested_python_streets.intersections import half_intervals, smooth
from tested_python_streets.mesh import Mesh, Assembler, cross
from tested_python_streets.paths import right
from tested_python_streets.quad_charts import LayeredQuadChart, FixedRimQuality, ChartResult
from tested_python_streets.four_side_charts import FourSidedQuadChart
from tested_python_streets.separator_surfaces import PartitionedSeparatorChart
from tested_python_streets.sampling import CertifiedCurbSampler
from tested_python_streets.surface_fields import BoundaryJetField
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.slope_projection import SlopeProjection
from tested_python_streets.regions import _classify
from tested_python_streets.quad_rotations import QuadEdgeRotations


class SourceCoreChart:
    """Direct quad blocks between shared source curves and a finite hub."""
    def __init__(self, patterns, spacing, max_cells):
        self.patterns, self.spacing, self.max_cells = patterns, spacing, max_cells

    def local(self, boundary, remaining, corners=None):
        errors = []
        if corners is not None:
            try:
                mesh, counts = FourSidedQuadChart(self.patterns, max_cells=remaining).map(boundary, corners)
                mesh.validate()
                mesh, improvement = FixedRimQuality(iterations=300).improve(mesh)
                quality = mesh.validate()
                if quality['maximum_edge_ratio'] <= 8.:
                    return ChartResult(mesh, dict(construction='source-anchored four-sided atlas',
                        side_intervals=counts, quality=quality, improvement=improvement))
                errors.append(f'Source-anchored atlas edge ratio {quality["maximum_edge_ratio"]:.3f}')
            except ValueError as exc:
                errors.append(str(exc))
        for builder in (LayeredQuadChart(self.patterns, spacing=self.spacing, max_cells=remaining, quality_iterations=100),
                        FourSidedQuadChart(self.patterns, max_cells=remaining),
                        PartitionedSeparatorChart(self.patterns, spacing=self.spacing, max_cells=remaining)):
            try:
                return builder.build(boundary)
            except ValueError as exc:
                errors.append(str(exc))
                if isinstance(builder, LayeredQuadChart) and 'needs better grading' in str(exc):
                    raw = builder.build(boundary, require_quality=False)
                    initial = raw.mesh.validate()['maximum_edge_ratio']
                    if initial <= 16.:
                        repaired, weighted = FixedRimQuality(iterations=300, power=2).improve(raw.mesh)
                        repaired, rotations = QuadEdgeRotations(max_passes=16, max_cells=remaining).improve(
                            repaired, FixedRimQuality(iterations=300, power=2))
                        quality = repaired.validate()
                        if quality['maximum_edge_ratio'] <= 8.:
                            return ChartResult(repaired, dict(construction='source block with bounded weighted repair',
                                quality=quality, weighted=weighted, rotations=rotations))
                        errors.append(f'Weighted source-block repair remains at {quality["maximum_edge_ratio"]:.3f}')
        raise ValueError('Source core block has no accepted chart: ' + '; '.join(errors))

    def split_sector(self, polygon, keys, cut_a, keys_a, outer, outer_keys, cut_b, keys_b,
                     arc, arc_keys, remaining):
        from tested_python_streets.decomposition import VisibilitySplitter
        middle = (len(outer) - 1) // 2
        hub_middle = (len(arc) - 1) // 2
        a = len(cut_a) - 1 + middle
        b = len(cut_a) - 1 + len(outer) - 1 + len(cut_b) - 1 + hub_middle
        if not VisibilitySplitter().visible(polygon, a, b):
            raise ValueError('Corner-midpoint core cut is not inside its source sector')
        first, last = arc[hub_middle], outer[middle]
        length = float(np.linalg.norm(last - first))
        step = min(.5 * self.spacing, np.linalg.norm(arc[hub_middle] - arc[hub_middle - 1]) * .5)
        distances, cursor = [0.], 0.
        while step < self.spacing and cursor + step < length:
            cursor += step;distances.append(cursor);step *= 1.7
        count = max(1, int(math.ceil((length - cursor) / self.spacing)))
        distances = np.r_[distances, np.linspace(cursor, length, count + 1)[1:]]
        if (len(distances) - 1) % 2:
            at = int(np.argmax(np.diff(distances)))
            distances = np.sort(np.r_[distances, .5 * (distances[at] + distances[at + 1])])
        cross_cut = first + (distances / length)[:, None] * (last - first)
        cross_keys = [arc_keys[hub_middle]] + [('corner-cut', k) for k in range(1, len(cross_cut) - 1)] + [outer_keys[middle]]
        parts = (
            (np.vstack((cut_a[:-1], outer[:middle], cross_cut[::-1][:-1], arc[:hub_middle + 1][::-1][:-1])),
             keys_a[:-1] + outer_keys[:middle] + cross_keys[::-1][:-1] + arc_keys[:hub_middle + 1][::-1][:-1],
             (0, len(cut_a) - 1, len(cut_a) - 1 + middle, len(cut_a) - 1 + middle + len(cross_cut) - 1)),
            (np.vstack((cross_cut[:-1], outer[middle:-1], cut_b[::-1][:-1], arc[hub_middle:][::-1][:-1])),
             cross_keys[:-1] + outer_keys[middle:-1] + keys_b[::-1][:-1] + arc_keys[hub_middle:][::-1][:-1],
             (0, len(cross_cut) - 1, len(cross_cut) - 1 + middle,
              len(cross_cut) - 1 + middle + len(cut_b) - 1)))
        assembly, reports = Assembler(), []
        assembly.add(np.column_stack((polygon, np.zeros(len(polygon)))), np.empty((0, 4), np.int32), keys)
        for i, (boundary, local_keys, corners) in enumerate(parts):
            try:
                result = self.local(boundary, remaining - len(assembly.quads), corners)
            except ValueError as exc:
                raise ValueError(f'Corner half {i}: {exc}') from exc
            local_keys.extend(('half-core', i, k) for k in range(len(boundary), len(result.mesh.xyz)))
            assembly.add(result.mesh.xyz, result.mesh.quads, local_keys)
            reports.append(result.metadata)
        mesh = assembly.mesh()
        return ChartResult(mesh, dict(construction='paired source/corner-midpoint core blocks',
                                      quality=mesh.validate(), blocks=reports))

    def build(self, context, layout, ring, peak_indices):
        from tested_python_streets.quad_charts import ChartResult
        approaches = layout.approaches
        if any(abs(a.offset) > 1e-9 for a in approaches):
            raise ValueError('Source core experiment currently needs zero source offsets')
        origin = approaches[0].origin
        inside, boundary = _classify(origin[None], ring, 1e-8)
        if not inside[0] or boundary[0] >= 0:
            raise ValueError('Common source root lies outside the retained core')
        edge = np.roll(ring, -1, axis=0) - ring
        t = np.clip(np.sum((origin - ring) * edge, axis=1) / np.sum(edge ** 2, axis=1), 0., 1.)
        clearance = float(np.linalg.norm(ring + t[:, None] * edge - origin, axis=1).min())
        angles = np.array([math.atan2(a.direction[1], a.direction[0]) % math.tau for a in approaches])
        gaps = np.diff(np.r_[angles, angles[0] + math.tau])
        radius = min(.5 * clearance, max(.25, self.spacing / (5. * gaps.min())))
        if radius <= 1e-3 or radius >= min(layout.outward_stations):
            raise ValueError('No finite common hub fits the source core')
        cuts, cut_keys = [], []
        for i, (a, stop) in enumerate(zip(approaches, layout.outward_stations)):
            curvature = max(getattr(p.path, 'curvature_bound', 0.) for p in a.path.pieces)
            step = min(self.spacing, math.sqrt(8. * context.network.settings.curve_tolerance / max(curvature, 1e-12)))
            near_step = min(.5 * step, radius * gaps.min() * .5)
            near = [radius]
            cursor = radius
            while near_step < step and cursor + near_step < stop:
                cursor += near_step;near.append(cursor);near_step *= 1.7
            count = max(1, int(math.ceil((stop - cursor) / step)))
            s = np.r_[near, np.linspace(cursor, stop, count + 1)[1:]]
            if (len(s) - 1) % 2:
                k = int(np.argmax(np.diff(s)))
                s = np.sort(np.r_[s, .5 * (s[k] + s[k + 1])])
            count = len(s) - 1
            xy = a.frame(s)[0]
            xy[-1] = ring[peak_indices[i]]
            samples = np.vstack((xy[:-1], .5 * (xy[:-1] + xy[1:])))
            hit, on = _classify(samples, ring, 1e-7)
            if not np.all(hit | (on >= 0)):
                raise ValueError(f'Source centerline cut {i} exits the core')
            cuts.append(xy)
            cut_keys.append([('hub', i)] + [('source-cut', i, j) for j in range(1, count)] + [('ring', peak_indices[i])])
        hub_points = [xy[0] for xy in cuts]
        hub_angles = [math.atan2(*(p - origin)[::-1]) % math.tau for p in hub_points]
        arcs, arc_keys = [], []
        for i in range(len(approaches)):
            j = (i + 1) % len(approaches)
            gap = (hub_angles[j] - hub_angles[i]) % math.tau
            if gap < 1e-4 or gap > gaps[i] + .1:
                raise ValueError('Hub ordering disagrees with source incidence')
            outer_count = (peak_indices[j] - peak_indices[i]) % len(ring)
            count = max(1, int(math.ceil(radius * gap / min(.5 * self.spacing, .2))))
            if (len(cuts[i]) - 1 + outer_count + len(cuts[j]) - 1 + count) % 2:
                count += 1
            if (count // 2 - outer_count // 2) % 2:
                count += 2
            u = np.linspace(0., 1., count + 1)
            r = ((1. - u) * np.linalg.norm(hub_points[i] - origin)
                 + u * np.linalg.norm(hub_points[j] - origin))
            angle = hub_angles[i] + u * gap
            xy = origin + r[:, None] * np.column_stack((np.cos(angle), np.sin(angle)))
            xy[0], xy[-1] = hub_points[i], hub_points[j]
            arcs.append(xy)
            arc_keys.append([('hub', i)] + [('hub-arc', i, k) for k in range(1, count)] + [('hub', j)])
        assembly = Assembler()
        assembly.add(np.column_stack((ring, np.zeros(len(ring)))), np.empty((0, 4), np.int32),
                     [('ring', i) for i in range(len(ring))])
        reports = []
        hub_xy, hub_keys = [], []
        for i in range(len(approaches)):
            j = (i + 1) % len(approaches)
            outer_ids = [(peak_indices[i] + k) % len(ring)
                         for k in range((peak_indices[j] - peak_indices[i]) % len(ring) + 1)]
            polygon = np.vstack((cuts[i][:-1], ring[outer_ids[:-1]], cuts[j][::-1][:-1], arcs[i][::-1][:-1]))
            keys = (cut_keys[i][:-1] + [('ring', k) for k in outer_ids[:-1]]
                    + cut_keys[j][::-1][:-1] + arc_keys[i][::-1][:-1])
            corners = (0, len(cuts[i]) - 1,
                       len(cuts[i]) - 1 + len(outer_ids) - 1,
                       len(cuts[i]) - 1 + len(outer_ids) - 1 + len(cuts[j]) - 1)
            try:
                result = self.local(polygon, self.max_cells - len(assembly.quads), corners)
            except ValueError as exc:
                try:
                    result = self.split_sector(polygon, keys, cuts[i], cut_keys[i], ring[outer_ids],
                        [('ring', k) for k in outer_ids], cuts[j], cut_keys[j], arcs[i], arc_keys[i],
                        self.max_cells - len(assembly.quads))
                except ValueError as alternate:
                    raise ValueError(f'Source sector {i}: {exc}; midpoint alternative: {alternate}') from alternate
            keys.extend(('source-core', i, k) for k in range(len(polygon), len(result.mesh.xyz)))
            assembly.add(result.mesh.xyz, result.mesh.quads, keys)
            reports.append(result.metadata)
            hub_xy.extend(arcs[i][:-1]);hub_keys.extend(arc_keys[i][:-1])
        hub = self.local(np.asarray(hub_xy), self.max_cells - len(assembly.quads))
        hub_keys.extend(('hub-core', k) for k in range(len(hub_xy), len(hub.mesh.xyz)))
        assembly.add(hub.mesh.xyz, hub.mesh.quads, hub_keys)
        mesh = assembly.mesh()
        mesh.seams = {0: np.arange(len(ring))}
        return ChartResult(mesh, dict(construction='source centerline core sectors and finite hub',
            hub_radius=radius, blocks=reports, hub=hub.metadata, quality=mesh.validate()))


class CornerStrips:
    def __init__(self, patterns, *, max_cells=50000):
        self.patterns, self.max_cells = patterns, max_cells

    @staticmethod
    def normal(boundary, t):
        t = np.asarray(t)
        step = 1e-6
        a, b = np.maximum(0., t - step), np.minimum(1., t + step)
        tangent = boundary(b) - boundary(a)
        norm = np.linalg.norm(tangent, axis=-1)
        if np.any(norm <= 1e-12):
            raise ValueError('Normal strip has an undefined curb tangent')
        return -right(tangent / norm[..., None])

    def guide(self, layout, side, t, fraction):
        a, b = layout.approaches[side], layout.approaches[(side + 1) % len(layout.approaches)]
        boundary = layout.boundaries[side]
        mix = smooth(np.asarray(t))
        width = .5 * (1. - fraction) * ((1. - mix) * a.width + mix * b.width)
        return boundary(t) + self.normal(boundary, t) * width[..., None]

    def build(self, context, site, layout, *, collar_intervals=2):
        m = half_intervals(context.network.settings)
        if not 1 <= collar_intervals < m:
            raise ValueError('Core collar must retain a finite part of each half profile')
        fraction = collar_intervals / m
        radial = m - collar_intervals
        approaches = layout.approaches
        datum = approaches[0].profile.datum
        rise = approaches[0].profile.rise
        if any(a.profile.half or abs(a.profile.rise - rise) > 1e-10 or abs(a.profile.datum - datum) > 1e-8
               for a in approaches):
            raise ValueError('This experiment needs full equal crowns and one datum')
        core_floor, ceiling = datum + rise * (1. - fraction ** 2), datum + rise
        sampler = CertifiedCurbSampler(context.network.settings.curve_tolerance, context.network.settings.spacing, 20000)
        assembly = Assembler()
        core_xy, core_keys, rules, port_specs, ring_z = [], [], [], [], []
        band_reports = []
        curvature_samples = []
        for side, boundary in enumerate(layout.boundaries):
            a, b = approaches[side], approaches[(side + 1) % len(approaches)]
            t = sampler.sample(boundary, 0., 1.)
            if (len(t) - 1) % 2:
                at = int(np.argmax(np.diff(t)))
                t = np.sort(np.r_[t, .5 * (t[at] + t[at + 1])])
            rows = len(t) - 1
            pattern = self.patterns.grid(rows, radial)
            tt = np.interp(pattern.uv[:, 0], np.linspace(0., 1., rows + 1), t)
            v = pattern.uv[:, 1]
            mix = smooth(tt)
            widths = .5 * (1. - fraction) * ((1. - mix) * a.width + mix * b.width)
            xy = boundary(tt) + self.normal(boundary, tt) * (widths * v)[:, None]
            q = -1. + (1. - fraction) * v
            z = datum + rise * (1. - q ** 2)
            def key(row, column):
                if row == 0:
                    return ('mouth', side, column)
                if row == rows:
                    return ('mouth', (side + 1) % len(approaches), 2 * m - column)
                return ('guide', side, row) if column == radial else ('sector', side, row, column)
            keys = []
            for i in range(rows + 1):
                for j in range(radial + 1):
                    index = i * (radial + 1) + j
                    if i in (0, rows):
                        approach, distance = (a, layout.outward_stations[side]) if i == 0 else (b, layout.outward_stations[(side + 1) % len(approaches)])
                        qq = -1. + j / m if i == 0 else 1. - j / m
                        xy[index] = approach.point(distance, approach.offset + .5 * approach.width * qq)
                        z[index] = datum + approach.profile.value(qq)
                    keys.append(key(i, j))
            band = Mesh(np.column_stack((xy, z)), pattern.quads)
            quality = band.validate()
            if quality['maximum_edge_ratio'] > 8.:
                raise ValueError(f'Normal strip {side} needs better grading: {quality["maximum_edge_ratio"]:.3f}')
            band_reports.append(dict(side=side, rows=rows, radial=radial, quality=quality))
            assembly.add(band.xyz, band.quads, keys)
            guide_ids = np.arange(rows + 1) * (radial + 1) + radial
            guide_xy = xy[guide_ids]
            # Gradients along the core guide are normal to its actual chords.
            # Full mouth jets below remain exact ProfileEdge rules.
            normal_b = self.normal(boundary, t)
            dt = 1e-5
            lo, hi = np.maximum(0., t - dt), np.minimum(1., t + dt)
            velocity = (self.guide(layout, side, hi, fraction) - self.guide(layout, side, lo, fraction)) / (hi - lo)[:, None]
            width = .5 * (1. - fraction) * ((1. - smooth(t)) * a.width + smooth(t) * b.width)
            denominator = cross(velocity, width[:, None] * normal_b)
            if denominator.min() <= 1e-10:
                raise ValueError('Normal strip inner guide folds')
            height_v = 2. * rise * fraction * (1. - fraction)
            gradients = height_v * -right(velocity) / denominator[:, None]
            for j in range(rows):
                vector = guide_xy[j + 1] - guide_xy[j]
                length = float(np.linalg.norm(vector))
                tangent = vector / length
                inward = -right(tangent)
                slopes = np.maximum(0., gradients[j:j + 2] @ inward)
                core_xy.append(guide_xy[j])
                ring_z.append(core_floor)
                core_keys.append(key(j, radial))
                rules.append(CurbEdge(core_floor, core_floor, float(slopes[0]), float(slopes[1]),
                                      length, 0., length, tuple(tangent)))
            next_side = (side + 1) % len(approaches)
            port_start = len(core_xy)
            count = 2 * collar_intervals
            q_values = np.linspace(fraction, -fraction, count + 1)
            for j, qq in enumerate(q_values[:-1]):
                core_xy.append(b.point(layout.outward_stations[next_side], b.offset + .5 * b.width * qq))
                ring_z.append(datum + b.profile.value(qq))
                core_keys.append(('mouth', next_side, m + collar_intervals - j))
                rules.append(ProfileEdge(b, layout.outward_stations[next_side], .5 * (1. + qq),
                                         .5 * (1. + q_values[j + 1])))
            port_specs.append((next_side, port_start, count, q_values))
        ring = np.asarray(core_xy)
        ring_heights = np.asarray(ring_z)
        ports = tuple(Port(site.key, i, approaches[i].source.id, approaches[i].end,
                           layout.outward_stations[i], tuple((start + j) % len(ring) for j in range(count + 1)),
                           tuple(.5 * (1. + q))) for i, start, count, q in port_specs)
        gradients = np.asarray([rule.sample(np.asarray(0.))[1] for rule in rules])
        conditions = BoundaryData(ring, ring_heights, gradients, tuple(rules), ports,
                                  core_floor, ceiling, False)
        field = BoundaryJetField((conditions,))
        remaining = self.max_cells - len(assembly.quads)
        failures, core = [], None
        class SourceBuilder:
            def build(self, boundary):
                peak_indices = {i: (start + collar_intervals) % len(ring) for i, start, count, q in port_specs}
                return SourceCoreChart(context.patterns, context.network.settings.spacing, remaining).build(
                    context, layout, boundary, peak_indices)
        for name, builder in (
            ('layered', LayeredQuadChart(self.patterns, spacing=context.network.settings.spacing, max_cells=remaining)),
            ('four-sided', FourSidedQuadChart(self.patterns, max_cells=remaining)),
            ('partitioned', PartitionedSeparatorChart(self.patterns, spacing=context.network.settings.spacing, max_cells=remaining)),
            ('source-sectors', SourceBuilder())):
            try:
                result = builder.build(ring)
                core = result.mesh
                core.seams = {0: np.arange(len(ring))}
                core.material = 'pavement'
                break
            except ValueError as exc:
                failures.append(f'{name}: {exc}')
        if core is None:
            error = ValueError('Core has no accepted quad chart: ' + '; '.join(failures))
            error.ring = ring
            error.bands = assembly.mesh()
            raise error
        core.xyz[:, 2] = field(core.xyz[:, :2])
        core.xyz[:len(ring), 2] = ring_heights
        core_fixed = np.arange(len(core.xyz)) < len(ring)
        limit = max(.10, 1.1 * max(4. * rise / a.width for a in approaches))
        projection = None
        if core.quad_slopes().max() > limit + 1e-8:
            core, projection = SlopeProjection(max_cells=remaining).solve(core, core_fixed,
                datum=core_floor, ceiling=ceiling, maximum_slope=limit)
        core_keys.extend(('core', i) for i in range(len(ring), len(core.xyz)))
        assembly.add(core.xyz, core.quads, core_keys)
        mesh = assembly.mesh()
        quality = mesh.validate()
        slopes = mesh.slope_statistics()
        if quality['maximum_edge_ratio'] > 8. or slopes['maximum'] > limit + 1e-8:
            raise ValueError(f'Corner-strip surface quality/slopes exceed bounds: {quality}, {slopes}')
        if mesh.xyz[:, 2].min() < datum - 1e-8 or mesh.xyz[:, 2].max() > ceiling + 1e-8:
            raise ValueError('Corner-strip surface exceeds the full crown envelope')
        return mesh, dict(quality=quality, physical_slopes=slopes, core_chart=name, core_failures=failures,
            core_projection=projection, strips=band_reports, collar_fraction=fraction,
            core_lower_height=core_floor, ceiling=ceiling), ring
