"""Alternative direct quad charts retaining a separator's complete XYZ rim."""
from dataclasses import dataclass, replace
import math
import numpy as np
from .quad_charts import LayeredQuadChart, ChartResult, QuadChartQualityError, FixedRimQuality
from .four_side_charts import FourSidedQuadChart, quality as quad_quality
from .profiles import plateau_extension
from .regions import signed_area, Region, RegionUnion
from .mesh import Assembler
from .model import Key


class SeparatorStripCuts:
    """Bound visible block length with noncrossing transverse computational cuts.

    Cut ends are existing physical vertices. Axial projection only proposes
    candidates; the visibility test certifies that each complete chord is inside
    its polygon. A final visibility pass qualifies every resulting local block.
    """
    def __init__(self, spacing):
        self.spacing = spacing

    def cuts(self, group, context):
        from .decomposition import VisibilitySplitter
        from .partition import BoundaryRef, JunctionCut
        boundary = group.loops[0].xy
        splitter = VisibilitySplitter()
        pieces, cuts = [tuple(range(len(boundary)))], []

        def apply(first, last):
            if tuple(sorted((first, last))) in {tuple(sorted(c)) for c in cuts}:
                return
            for pi, indices in enumerate(pieces):
                if first not in indices or last not in indices:
                    continue
                i, j = indices.index(first), indices.index(last)
                distance = (j - i) % len(indices)
                if min(distance, len(indices) - distance) < 2:
                    continue
                a = tuple(indices[(i + k) % len(indices)] for k in range(distance + 1))
                b = tuple(indices[(j + k) % len(indices)] for k in range(len(indices) - distance + 1))
                pieces[pi:pi + 1] = (a, b)
                cuts.append((first, last))
                if len(cuts) > 128:
                    raise ValueError('Separator strip decomposition exceeds its cut budget')
                return
            raise ValueError('Separator computational cut has no unique containing block')

        def visible_split(indices):
            for first, last in splitter.split(boundary[list(indices)]):
                apply(indices[first], indices[last])

        visible_split(pieces[0])
        for indices in tuple(pieces):
            polygon = boundary[list(indices)]
            centered = polygon - polygon.mean(axis=0)
            _, vectors = np.linalg.eigh(centered.T @ centered)
            axis = vectors[:, -1]
            axial = centered @ axis
            length = float(np.ptp(axial))
            step = max(12., 24. * self.spacing)
            count = int(math.ceil(length / step))
            if count > 65:
                raise ValueError('Separator strip decomposition exceeds its longitudinal budget')
            edge = np.roll(polygon, -1, axis=0) - polygon
            delta = np.roll(axial, -1) - axial
            used, previous = set(), -math.inf
            for station in np.linspace(axial.min(), axial.max(), count + 1)[1:-1]:
                t = np.divide(station - axial, delta, out=np.full(len(delta), np.inf), where=abs(delta) > 1e-12)
                hits = []
                for ei in np.flatnonzero((t >= 0.) & (t < 1.)):
                    point = polygon[ei] + t[ei] * edge[ei]
                    if not any(np.linalg.norm(point - q[1]) < 1e-7 for q in hits):
                        hits.append((int(ei), point))
                if len(hits) != 2:
                    continue
                options = []
                for a in (hits[0][0], (hits[0][0] + 1) % len(indices)):
                    for b in (hits[1][0], (hits[1][0] + 1) % len(indices)):
                        distance = abs(a - b)
                        if min(distance, len(indices) - distance) < 3 or a in used or b in used:
                            continue
                        lo, hi = sorted((float(axial[a]), float(axial[b])))
                        if lo <= previous + 1e-7 or not splitter.visible(polygon, a, b):
                            continue
                        score = abs(axial[a] - station) + abs(axial[b] - station)
                        options.append((float(score), a, b, hi))
                if options:
                    _, a, b, previous = min(options)
                    apply(indices[a], indices[b])
                    used.update((a, b))
        for indices in tuple(pieces):
            visible_split(indices)
        for first, last in cuts:
            yield JunctionCut(None, group.members[0], group.members[0],
                              BoundaryRef(0, first, 0.), BoundaryRef(0, last, 0.))


class PartitionedSeparatorChart:
    """Split a simple rim into visible blocks and weld declared internal seams."""
    def __init__(self, patterns, *, spacing=.7, max_cells=50000, interior_spacing=None):
        self.patterns, self.spacing, self.max_cells = patterns, spacing, max_cells
        self.interior_spacing = spacing if interior_spacing is None else interior_spacing
        if not math.isfinite(self.interior_spacing) or self.interior_spacing <= 0.:
            raise ValueError('Separator interior spacing must be finite and positive')

    def build(self, boundary):
        from .compounds import CompoundGroup, SampledFootprint, RimSegment
        from .partition import CompoundPartitioner
        owner = Key('separator_chart', 0)
        footprint = SampledFootprint(Region(owner, boundary),
            tuple(RimSegment('fixed', i, 0., 1.) for i in range(len(boundary))))
        group = CompoundGroup((owner,), (), (footprint,), RegionUnion().build((footprint.region,)))
        partition = CompoundPartitioner(policy=SeparatorStripCuts(self.spacing)).build(group, None)
        samples = partition.sample(self.spacing)
        original = {tuple(samples.xy[i]): i for i in samples.boundary_loops[0]}
        restoration = None
        if len(original) != len(boundary) or any(tuple(p) not in original for p in boundary):
            samples, original, restoration = self.source_boundary(group, samples, boundary)
        assembly = Assembler()
        assembly.add(np.column_stack((boundary, np.zeros(len(boundary)))), np.empty((0, 4), dtype=np.int32),
                     [('partition', int(original[tuple(p)])) for p in boundary])
        reports, failures, total = [], [], 0
        for bi, rim in enumerate(samples.rims):
            remaining = self.max_cells - total
            if remaining < 4:
                raise ValueError('Separator decomposition exceeds its cell budget')
            child, unfinished = None, []
            for name, builder in (
                ('layered', LayeredQuadChart(self.patterns, spacing=self.interior_spacing, max_cells=remaining, quality_iterations=100)),
                ('four-sided', FourSidedQuadChart(self.patterns, max_cells=remaining, quality_iterations=100))):
                try:
                    child = builder.build(samples.xy[list(rim)])
                    break
                except ValueError as error:
                    failures.append(f'block {bi} {name}: {error}')
                    if (isinstance(error, QuadChartQualityError)
                            and error.result.metadata['quality']['maximum_edge_ratio'] <= 16.):
                        unfinished.append(error.result)
            # Keep the existing successful chart choices byte-exact. Only if
            # both fail, let a positive block finish a short quality pass and
            # stop at the unchanged acceptance limit, before shape-energy
            # descent can worsen the maximum edge ratio again.
            if child is None:
                for previous in unfinished:
                    candidate, tail = FixedRimQuality(iterations=25, target_ratio=8.).improve(previous.mesh)
                    report = quad_quality(candidate)
                    if report['maximum_edge_ratio'] <= 8. and report['minimum_corner_sine'] >= .03:
                        if not np.array_equal(candidate.xyz[:len(rim)], previous.mesh.xyz[:len(rim)]):
                            raise ValueError('Separator quality continuation changed its physical rim')
                        child = ChartResult(candidate, dict(previous.metadata, quality=report,
                            acceptance_continuation=tail, acceptance_continuation_limit=25))
                        break
            if child is None:
                raise ValueError('A separator block has no accepted quad chart: ' + '; '.join(failures[-2:]))
            total += len(child.mesh.quads)
            keys = [('partition', int(node)) for node in rim]
            keys.extend(('separator-block', bi, i) for i in range(len(rim), len(child.mesh.xyz)))
            assembly.add(child.mesh.xyz, child.mesh.quads, keys)
            reports.append(child.metadata)
        mesh = assembly.mesh('grass')
        metadata = dict(construction='visibility-partitioned direct quad separator',
            blocks=reports, computational_cuts=len(partition.cuts), area_error=partition.area_error,
            quality=mesh.validate(), failed_block_charts=failures)
        if restoration is not None:
            metadata['source_boundary_restoration'] = restoration
        return ChartResult(mesh, metadata)

    @staticmethod
    def source_boundary(group, samples, boundary):
        """Restore arithmetic roundoff from exact source endpoint incidence.

        Coordinate proximity never chooses a source vertex. Each external node
        must refer to one original segment endpoint, with a complete bijection
        onto the prescribed rim. New samples and changed topology stay errors.
        """
        message = 'Separator decomposition changed the physical boundary samples'
        if len(samples.boundary_loops) != 1 or len(samples.boundary_loops[0]) != len(boundary):
            raise ValueError(message)
        nodes = {}
        for node in samples.boundary_loops[0]:
            ref = samples.boundary_refs.get(node)
            if ref is None or ref.loop != 0:
                raise ValueError(message)
            span = group.loops[ref.loop].spans[ref.edge]
            if span.owner != group.members[0] or not 0 <= span.segment < len(boundary):
                raise ValueError(message)
            fraction = span.first + ref.fraction * (span.last - span.first)
            if fraction == 0.:
                index = span.segment
            elif fraction == 1.:
                index = (span.segment + 1) % len(boundary)
            else:
                raise ValueError(message)
            if index in nodes:
                raise ValueError(message)
            nodes[index] = int(node)
        if set(nodes) != set(range(len(boundary))):
            raise ValueError(message)
        ordered = np.asarray([nodes[i] for i in range(len(boundary))])
        displacement = np.linalg.norm(samples.xy[ordered] - boundary, axis=1)
        limit = 256. * np.finfo(float).eps * max(1., float(abs(boundary).max()))
        if not np.all(np.isfinite(displacement)) or displacement.max() > limit:
            raise ValueError(message)
        xy = samples.xy.copy()
        xy[ordered] = boundary
        xy.flags.writeable = False
        restored = replace(samples, xy=xy)
        original = {tuple(boundary[i]): nodes[i] for i in range(len(boundary))}
        return restored, original, dict(source_identity_bijection=True,
            changed_vertices=int(np.any(samples.xy[ordered] != boundary, axis=1).sum()),
            maximum_roundoff=float(displacement.max()), roundoff_limit=limit,
            no_new_physical_samples=True)


@dataclass
class SeparatorSurface:
    mesh: object
    metadata: dict


class FixedRimSeparator:
    """Reuse accepted planar charts; choose only interior topology/coordinates.

    The existing road/junction boundary samples, heights and footprint remain
    authoritative. Both chart families construct quads directly. Their normal
    quality/cell budgets remain mandatory; an unavailable chart is diagnostic.
    """
    def __init__(self, patterns, *, spacing=.7, max_cells=50000):
        if not math.isfinite(spacing) or spacing <= 0. or type(max_cells) is not int or max_cells < 4:
            raise ValueError('Invalid fixed-rim separator chart limits')
        self.patterns, self.spacing, self.max_cells = patterns, spacing, max_cells

    def build(self, boundary, *, datum, height, exaggeration):
        boundary = np.asarray(boundary, dtype=float)
        if (boundary.ndim != 2 or boundary.shape[1] != 3 or len(boundary) < 4 or len(boundary) % 2
                or not np.all(np.isfinite(boundary)) or not all(map(math.isfinite, (datum, height, exaggeration)))
                or min(height, exaggeration) < 0. or signed_area(boundary[:, :2]) <= 0.):
            raise ValueError('Fixed-rim separator needs an even counterclockwise XYZ boundary and finite heights')
        if boundary[:, 2].min() < datum - 1e-8 or boundary[:, 2].max() > datum + height + 1e-8:
            raise ValueError('Separator boundary exceeds its plateau envelope')
        failures = []
        from .separator_angular import AngularSeparatorChart
        builders = (('layered', LayeredQuadChart(self.patterns, spacing=self.spacing,
                         max_cells=self.max_cells, quality_iterations=100)),
                    ('four-sided', FourSidedQuadChart(self.patterns, max_cells=self.max_cells,
                         quality_iterations=100)),
                    ('partitioned', PartitionedSeparatorChart(self.patterns, spacing=self.spacing,
                         interior_spacing=2. * self.spacing, max_cells=self.max_cells)),
                    ('angular', AngularSeparatorChart(self.patterns, spacing=self.spacing,
                         max_cells=self.max_cells)))
        for name, builder in builders:
            try:
                chart = builder.build(boundary[:, :2])
            except ValueError as error:
                failures.append(f'{name}: {error}')
                continue
            mesh = chart.mesh
            n = len(boundary)
            if not np.array_equal(mesh.xyz[:n, :2], boundary[:, :2]):
                raise ValueError('Alternative separator chart changed an authoritative boundary vertex')
            mesh.xyz[:, 2], extension = plateau_extension(mesh.xyz[:, :2], boundary, datum, height)
            mesh.xyz[:n] = boundary
            mesh.material = 'grass'
            mesh.seams = {0: np.arange(n)}
            bound = extension.pop('continuous_slope_bound')
            extension['physical_continuous_slope_bound'] = bound / exaggeration if exaggeration else 0.
            return SeparatorSurface(mesh, dict(chart=chart.metadata, height_extension=extension,
                quality=mesh.validate(), physical_slopes=mesh.slope_statistics(exaggeration) if exaggeration else
                    dict(maximum=0., percentile_99=0., samples_per_quad=9),
                boundary_vertices=n, failed_charts=failures, selected_chart=name,
                interior_spacing=2. * self.spacing if name == 'partitioned' else self.spacing,
                construction='fixed-rim direct quad separator chart'))
        raise ValueError('No accepted fixed-rim separator chart: ' + '; '.join(failures))
