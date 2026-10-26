"""Source-aware planar partitions of compound junctions, without triangulation.

Cuts join boundary ownership changes of neighboring source junctions. Derived
seams carry no prescribed carriageway profile: the compound height field must
be reconciled from its external conditions. Every original source ID remains.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass, replace
import math
import numpy as np
from .mesh import Assembler, cross
from .model import Key
from .quad_charts import LayeredQuadChart, visibility_center
from .regions import _Nodes, _classify, signed_area
from .compounds import rim_label


@dataclass(frozen=True)
class BoundaryRef:
    loop: int
    edge: int
    fraction: float


@dataclass(frozen=True)
class JunctionCut:
    section_id: int | None
    first_owner: object
    last_owner: object
    first: BoundaryRef
    last: BoundaryRef
    section_ids: tuple = ()

    @property
    def source_sections(self):
        """All original Sections represented by one computational cut."""
        return self.section_ids or (() if self.section_id is None else (self.section_id,))


@dataclass(frozen=True)
class BoundarySpan:
    loop: int
    edge: int
    first: float
    last: float


@dataclass(frozen=True)
class CutSpan:
    cut: int
    first: float
    last: float


@dataclass(frozen=True)
class PartitionEdge:
    first: int
    last: int
    boundary: BoundarySpan | None = None
    cut: CutSpan | None = None


@dataclass(frozen=True)
class PartitionBlock:
    vertices: tuple
    edges: tuple
    owners: tuple
    center: tuple
    area: float


class CutPolicy(ABC):
    @abstractmethod
    def cuts(self, group, context):
        """Return source-linked cuts with endpoints on the original boundary."""


class BoundaryOwnershipCuts(CutPolicy):
    def cuts(self, group, context):
        merged = {}
        for cut in self._cuts(group, context):
            refs = tuple(sorted(((cut.first.loop, cut.first.edge, cut.first.fraction),
                                 (cut.last.loop, cut.last.edge, cut.last.fraction))))
            token = tuple(sorted((cut.first_owner, cut.last_owner))), refs
            if token not in merged:
                merged[token] = cut
            else:
                previous = merged[token]
                aliases = tuple(sorted(set(previous.source_sections) | set(cut.source_sections)))
                merged[token] = replace(previous, section_ids=aliases)
        yield from merged.values()

    def _cuts(self, group, context):
        if group.physical_layout is not None or not group.claims:
            if (len(group.loops) != 1 or group.loops[0].hole
                    or (group.physical_layout is None and len(group.members) != 1)):
                raise ValueError('Non-source compound cuts need another region policy')
            from .decomposition import VisibilitySplitter
            loop = group.loops[0]
            sources = {p.region.owner: p for p in group.footprints}
            representative = {}
            for i, span in enumerate(loop.spans):
                label = rim_label(span, sources)
                if label.kind == 'curb':
                    token = span.owner, label.index
                    distance = abs(.5 * (label.first + label.last) - .5)
                    if token not in representative or distance < representative[token][0]:
                        representative[token] = distance, i
            preferred = [i for distance, i in representative.values()]
            for first, last in VisibilitySplitter().split(loop.xy, preferred=preferred):
                yield JunctionCut(None, group.members[0], group.members[0],
                                  BoundaryRef(0, first, 0.), BoundaryRef(0, last, 0.))
            return
        for claim in group.claims:
            owners = {claim.first_owner, claim.last_owner}
            switches = []
            for li, loop in enumerate(group.loops):
                for ei, span in enumerate(loop.spans):
                    previous = loop.spans[ei - 1]
                    if previous.owner != span.owner and {previous.owner, span.owner} == owners:
                        switches.append(BoundaryRef(li, ei, 0.))
            if len(switches) != 2:
                raise ValueError('Neighboring junctions need another cut policy: no unique pair of boundary ownership changes')
            path = context[Key('path', claim.section_id)]
            station = .5 * (claim.first_trim + claim.last_trim)
            normal = path.frame(station)[2]
            first, last = (group.loops[ref.loop].xy[ref.edge] for ref in switches)
            if (last - first) @ normal < 0.:
                switches.reverse()
            yield JunctionCut(claim.section_id, claim.first_owner, claim.last_owner, *switches)


@dataclass
class SampledPartition:
    xy: np.ndarray
    rims: tuple
    boundary_refs: dict
    boundary_loops: tuple
    cut_intervals: dict


@dataclass
class CompoundPartition:
    group: object
    xy: np.ndarray
    edges: tuple
    blocks: tuple
    cuts: tuple
    boundary_refs: dict
    boundary_loops: tuple
    area_error: float

    def sample(self, spacing):
        """Choose shared seam intervals with compatible quad-boundary parity."""
        if not math.isfinite(spacing) or spacing <= 0.:
            raise ValueError('Partition seam spacing must be finite and positive')
        incidence = {}
        parity = []
        for bi, block in enumerate(self.blocks):
            parity.append(sum(self.edges[ei].boundary is not None for ei, forward in block.edges) % 2)
            for ei, forward in block.edges:
                if self.edges[ei].cut is not None:
                    incidence.setdefault(ei, []).append(bi)
        adjacent = {i: [] for i in range(len(self.blocks))}
        for ei, owners in incidence.items():
            if len(owners) != 2 or owners[0] == owners[1]:
                raise ValueError('An internal compound seam needs two distinct pavement blocks')
            a, b = owners
            adjacent[a].append((b, ei))
            adjacent[b].append((a, ei))
        parent, order = {0: (None, None)}, [0]
        for node in order:
            for child, edge in sorted(adjacent[node]):
                if child not in parent:
                    parent[child] = (node, edge)
                    order.append(child)
        if len(parent) != len(self.blocks):
            raise ValueError('Compound block adjacency is disconnected')
        desired = {ei: 0 for ei in incidence}
        for child in order[:0:-1]:
            ancestor, ei = parent[child]
            desired[ei] = parity[child]
            parity[ancestor] ^= parity[child]
            parity[child] = 0
        if parity[0]:
            raise ValueError('Compound external boundaries have incompatible quad parity')
        points = list(self.xy.copy())
        sequences, intervals = {}, {}
        for ei, edge in enumerate(self.edges):
            count = 1
            if edge.cut is not None:
                length = np.linalg.norm(self.xy[edge.last] - self.xy[edge.first])
                count = max(1, int(math.ceil(length / spacing)))
                count += (count % 2) != desired[ei]
                if count > 16384 or len(points) + count > 100000:
                    raise ValueError('Compound seam sampling exceeds its point budget')
                intervals[ei] = count
            ids = [edge.first]
            for j in range(1, count):
                ids.append(len(points))
                points.append((1. - j / count) * self.xy[edge.first] + (j / count) * self.xy[edge.last])
            ids.append(edge.last)
            sequences[ei] = tuple(ids)
        rims = []
        for block in self.blocks:
            rim = []
            for ei, forward in block.edges:
                sequence = sequences[ei] if forward else sequences[ei][::-1]
                rim.extend(sequence[:-1])
            if len(rim) % 2:
                raise ValueError('Shared seam sampling did not produce an even block rim')
            rims.append(tuple(rim))
        xy = np.asarray(points)
        xy.flags.writeable = False
        return SampledPartition(xy, tuple(rims), dict(self.boundary_refs), self.boundary_loops, intervals)


class CompoundPartitioner:
    def __init__(self, policy=None, *, epsilon=1e-7, max_nodes=16384, max_cuts=256):
        if (not math.isfinite(epsilon) or epsilon <= 0. or type(max_nodes) is not int or max_nodes < 4
                or type(max_cuts) is not int or max_cuts < 1):
            raise ValueError('Invalid compound partition limits')
        self.policy = BoundaryOwnershipCuts() if policy is None else policy
        self.epsilon, self.max_nodes, self.max_cuts = epsilon, max_nodes, max_cuts

    def build(self, group, context):
        epsilon = self.epsilon
        cuts = tuple(self.policy.cuts(group, context))
        if len(cuts) > self.max_cuts:
            raise ValueError('Compound partition exceeds its cut budget')
        source_sections = [sid for cut in cuts for sid in cut.source_sections]
        if len(set(source_sections)) != len(source_sections):
            raise ValueError('Several compound cuts claim the same source Section')
        registry = _Nodes(epsilon)
        refs = {}
        loops = [loop.xy for loop in group.loops]
        if sum(map(len, loops)) > self.max_nodes:
            raise ValueError('Compound partition exceeds its graph-node budget')
        for li, loop in enumerate(loops):
            for ei, point in enumerate(loop):
                node = registry.add(point)
                refs[node] = BoundaryRef(li, ei, 0.)
        fractions = [[[0., 1.] for point in loop] for loop in loops]
        sources = {p.region.owner: p for p in group.footprints}
        # Each physical loop is made even on an owned curb, never by adding a
        # new sample to an incoming profile shared with an independent road.
        for li, loop in enumerate(group.loops):
            if len(loop.xy) % 2:
                candidates = [i for i, span in enumerate(loop.spans)
                    if rim_label(span, sources).kind == 'curb'
                    and np.linalg.norm(loop.xy[(i + 1) % len(loop.xy)] - loop.xy[i]) > 4. * epsilon]
                if not candidates:
                    raise ValueError('No owned curb can resolve compound boundary parity')
                at = max(candidates, key=lambda i: np.linalg.norm(loop.xy[(i + 1) % len(loop.xy)] - loop.xy[i]))
                fractions[li][at].append(.5)

        def filled(point):
            total = 0
            for loop in group.loops:
                inside, boundary = _classify(np.asarray(point)[None], loop.xy, epsilon * .1)
                if boundary[0] >= 0:
                    return None
                total += (-1 if loop.hole else 1) * int(inside[0])
            return total > 0

        def position(ref):
            loop = loops[ref.loop]
            return (1. - ref.fraction) * loop[ref.edge] + ref.fraction * loop[(ref.edge + 1) % len(loop)]

        coordinates, parameters = [], []
        for cut in cuts:
            a, b = position(cut.first), position(cut.last)
            d = b - a
            if np.linalg.norm(d) <= epsilon or not filled(.5 * (a + b)):
                raise ValueError('A compound source cut does not pass through pavement')
            for loop in loops:
                x, e = loop, np.roll(loop, -1, axis=0) - loop
                determinant = cross(d, e)
                valid = abs(determinant) > 1e-12 * np.linalg.norm(d) * np.linalg.norm(e, axis=1)
                t = np.divide(cross(x - a, e), determinant, out=np.zeros(len(loop)), where=valid)
                u = np.divide(cross(x - a, d), determinant, out=np.zeros(len(loop)), where=valid)
                interior = (t * np.linalg.norm(d) > epsilon) & ((1. - t) * np.linalg.norm(d) > epsilon)
                if np.any(valid & interior & (u >= -1e-9) & (u <= 1. + 1e-9)):
                    raise ValueError('A compound source cut crosses or touches a physical boundary internally')
                for at in np.flatnonzero(~valid & (abs(cross(x - a, d)) <= epsilon * np.linalg.norm(d))):
                    projected = np.sort(np.array(((x[at] - a) @ d, (x[at] + e[at] - a) @ d)) / (d @ d))
                    if min(1., projected[-1]) - max(0., projected[0]) > epsilon / np.linalg.norm(d):
                        raise ValueError('A compound source cut follows a physical boundary')
            coordinates.append(np.array((a, b)))
            parameters.append([0., 1.])
            fractions[cut.first.loop][cut.first.edge].append(cut.first.fraction)
            fractions[cut.last.loop][cut.last.edge].append(cut.last.fraction)
        for i, first in enumerate(coordinates):
            a, d = first[0], first[1] - first[0]
            for j in range(i + 1, len(coordinates)):
                b, e = coordinates[j][0], coordinates[j][1] - coordinates[j][0]
                determinant = cross(d, e)
                if abs(determinant) <= 1e-12 * np.linalg.norm(d) * np.linalg.norm(e):
                    if abs(cross(b - a, d)) / np.linalg.norm(d) <= epsilon:
                        projected = np.sort(np.array(((b - a) @ d, (b + e - a) @ d)) / (d @ d))
                        if min(1., projected[-1]) - max(0., projected[0]) > epsilon / np.linalg.norm(d):
                            raise ValueError('Overlapping source cuts need a combined seam policy')
                    continue
                t, u = cross(b - a, e) / determinant, cross(b - a, d) / determinant
                if -1e-9 <= t <= 1. + 1e-9 and -1e-9 <= u <= 1. + 1e-9:
                    parameters[i].append(float(np.clip(t, 0., 1.)))
                    parameters[j].append(float(np.clip(u, 0., 1.)))
        edges, lookup, boundary_loops = [], {}, []

        def add_edge(first, last, boundary=None, cut=None):
            key = tuple(sorted((first, last)))
            if key in lookup or first == last:
                raise ValueError('Compound partition has duplicate or collapsed graph edges')
            lookup[key] = len(edges)
            edges.append(PartitionEdge(first, last, boundary, cut))

        for li, loop in enumerate(loops):
            rim = []
            for ei, (a, b) in enumerate(zip(loop, np.roll(loop, -1, axis=0))):
                ts = sorted(set(fractions[li][ei]))
                for lo, hi in zip(ts, ts[1:]):
                    if (hi - lo) * np.linalg.norm(b - a) <= epsilon:
                        continue
                    first, last = registry.add(a + lo * (b - a)), registry.add(a + hi * (b - a))
                    refs.setdefault(first, BoundaryRef(li, ei, float(lo)))
                    refs.setdefault(last, BoundaryRef(li, ei, float(hi)))
                    add_edge(first, last, BoundarySpan(li, ei, float(lo), float(hi)))
                    rim.append(first)
            if len(rim) % 2:
                raise ValueError('Compound physical-loop parity refinement failed')
            boundary_loops.append(tuple(rim))
        for ci, (a, b) in enumerate(coordinates):
            ts = sorted(set(parameters[ci]))
            for lo, hi in zip(ts, ts[1:]):
                if (hi - lo) * np.linalg.norm(b - a) <= epsilon:
                    continue
                first, last = registry.add(a + lo * (b - a)), registry.add(a + hi * (b - a))
                add_edge(first, last, cut=CutSpan(ci, float(lo), float(hi)))
        if len(registry.points) > self.max_nodes:
            raise ValueError('Compound partition exceeds its graph-node budget')
        adjacency, directed = {}, {}
        for ei, edge in enumerate(edges):
            adjacency.setdefault(edge.first, set()).add(edge.last)
            adjacency.setdefault(edge.last, set()).add(edge.first)
            directed[(edge.first, edge.last)] = (ei, True)
            directed[(edge.last, edge.first)] = (ei, False)
        for node, neighbors in adjacency.items():
            adjacency[node] = sorted(neighbors, key=lambda i: math.atan2(*(registry.points[i] - registry.points[node])[::-1]))
        pending, blocks = set(directed), []
        while pending:
            first = min(pending)
            current, vertices, labels = first, [], []
            while True:
                if current not in pending:
                    raise ValueError('Compound partition has an unresolved half-edge incidence')
                pending.remove(current)
                u, v = current
                vertices.append(u)
                labels.append(directed[current])
                neighbors = adjacency[v]
                current = (v, neighbors[(neighbors.index(u) - 1) % len(neighbors)])
                if current == first:
                    break
            xy = np.asarray([registry.points[i] for i in vertices])
            area = signed_area(xy)
            if area <= epsilon ** 2:
                continue
            vector = np.roll(xy, -1, axis=0) - xy
            lengths = np.linalg.norm(vector, axis=1)
            at = int(np.argmax(lengths))
            tangent = vector[at] / lengths[at]
            probe = .5 * (xy[at] + xy[(at + 1) % len(xy)]) + min(1e-5, .01 * lengths[at]) * np.array((-tangent[1], tangent[0]))
            if not filled(probe):
                continue
            if len(set(vertices)) != len(vertices):
                raise ValueError('Compound partition still has an uncut hole or boundary bridge')
            center = visibility_center(xy)
            owners = set()
            for ei, forward in labels:
                if edges[ei].cut is not None:
                    cut = cuts[edges[ei].cut.cut]
                    owners.add(cut.last_owner if forward else cut.first_owner)
            blocks.append(PartitionBlock(tuple(vertices), tuple(labels), tuple(sorted(owners)), tuple(center), area))
        expected = sum(loop.area for loop in group.loops)
        error = abs(sum(block.area for block in blocks) - expected)
        if not blocks or error > max(1e-7, expected * 1e-10):
            raise ValueError('Compound partition does not preserve its pavement area and holes')
        xy = np.asarray(registry.points)
        xy.flags.writeable = False
        return CompoundPartition(group, xy, tuple(edges), tuple(blocks), cuts, refs, tuple(boundary_loops), error)


@dataclass
class CompoundChartResult:
    mesh: object
    holes: dict
    partition: CompoundPartition
    samples: SampledPartition
    boundary_refs: dict
    metadata: dict


class QuadChartGradingError(ValueError):
    """A positive assembled chart exhausts the bounded shape correction."""


class CompoundQuadChart:
    """Conforming planar pavement and hole charts over shared source cuts.

    Output is preparation geometry at Z=0. External profile conditions and the
    compound height field must be applied before a surface provider publishes it.
    """
    def __init__(self, patterns, *, spacing=.7, partitioner=None, max_cells=50000,
                 quality_iterations=300):
        self.patterns, self.spacing = patterns, spacing
        self.partitioner = CompoundPartitioner() if partitioner is None else partitioner
        self.max_cells, self.quality_iterations = max_cells, quality_iterations

    def build(self, group, context, *, require_quality=True):
        result = self._build(group, context, require_quality=require_quality)
        if not require_quality or not result.metadata.get('initial_chart_failures'):
            return result
        # A partition cut is computational: after assembly its vertices are
        # interior degrees of freedom. Keeping that cut fixed during each
        # leaf's repair can prevent an otherwise valid global quad chart.
        # Reuse the first positive charts; do not repeat their construction or
        # proxy relaxation merely to suppress a final leaf quality exception.
        # Actual physical boundaries remain exact.
        try:
            mesh, details = self._correct(result.mesh)
        except QuadChartGradingError as original:
            from .corner_seed import CornerBandSeed
            try:
                alternative = CornerBandSeed(self.patterns, max_cells=self.max_cells).build(group, context, self.partitioner)
            except ValueError as failure:
                raise QuadChartGradingError(f'{original}; corner seed: {failure}') from original
            alternative.metadata = dict(alternative.metadata, previous_chart_failure=str(original))
            return alternative
        result.mesh = mesh
        for hole in result.holes.values():
            hole.mesh, correction = self._correct(hole.mesh)
            hole.metadata = dict(hole.metadata, assembled_quality_correction=correction,
                                 quality=hole.mesh.validate())
        result.metadata = dict(result.metadata, quality=mesh.validate(),
            assembled_quality_correction=details,
            block_metadata_phase='BEFORE_ASSEMBLED_CORRECTION')
        return result

    def _correct(self, mesh):
        from .quad_rotations import QuadEdgeRotations
        from .quad_charts import FixedRimQuality
        from .quality_correction import QuadShapeCorrection
        correction = QuadShapeCorrection(max_cells=self.max_cells)
        _, initial, _ = correction.evaluate(mesh.xyz[:, :2], mesh.quads)
        if correction.qualified(initial):
            return mesh, dict(qualified=True, iterations=0, initial=initial, final=initial,
                              exact_fixed_rim=True, stop='already_qualified')
        seed, rotations = QuadEdgeRotations(max_passes=8, max_cells=self.max_cells).improve(
            mesh, FixedRimQuality(iterations=self.quality_iterations, power=2))
        candidate, details = correction.improve(seed)
        if not details['qualified']:
            raise QuadChartGradingError('Assembled direct quad chart needs better grading: correction exhausted its finite budget')
        return candidate, dict(details, seed_rotations=rotations)

    def _build(self, group, context, *, require_quality):
        from .sampling import AnalyticalBoundarySampler
        from .quad_charts import QuadChartQualityError
        group = AnalyticalBoundarySampler(tolerance=context.network.settings.curve_tolerance,
                                           spacing=self.spacing).build(group)
        partition = self.partitioner.build(group, context)
        samples = partition.sample(self.spacing)
        assembly, reports, failures, total = Assembler(), [], [], 0
        for index, (block, rim) in enumerate(zip(partition.blocks, samples.rims)):
            chart = LayeredQuadChart(self.patterns, spacing=self.spacing,
                max_cells=self.max_cells - total, quality_iterations=self.quality_iterations)
            try:
                result = chart.build(samples.xy[list(rim)], require_quality=require_quality)
            except QuadChartQualityError as exc:
                if not require_quality:
                    raise
                from .four_side_charts import FourSidedQuadChart
                try:
                    result = FourSidedQuadChart(self.patterns, max_cells=self.max_cells - total,
                        quality_iterations=min(100, self.quality_iterations)).build(samples.xy[list(rim)])
                except ValueError as alternative:
                    result = exc.result
                    failures.append(dict(block=index, reason=f'{exc}; {alternative}'))
            total += len(result.mesh.quads)
            keys = [('partition', int(node)) for node in rim]
            keys.extend(('block', index, i) for i in range(len(rim), len(result.mesh.xyz)))
            assembly.add(result.mesh.xyz, result.mesh.quads, keys)
            reports.append(dict(owners=[owner.text() for owner in block.owners], **result.metadata))
        mesh = assembly.mesh()
        mesh.seams = {li: np.array([assembly.lookup[('partition', int(node))] for node in loop], dtype=np.int32)
                      for li, loop in enumerate(samples.boundary_loops)}
        boundary_refs = {assembly.lookup[('partition', int(node))]: ref for node, ref in samples.boundary_refs.items()}
        quality = mesh.validate()
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        if (len(mesh.xyz) - len(unique) + len(mesh.quads) != 2 - len(group.loops)
                or not set(counts).issubset({1, 2})
                or int(np.sum(counts == 1)) != sum(map(len, samples.boundary_loops))):
            raise ValueError('Compound quad charts do not form the required conforming pavement region')
        holes = {}
        for li, loop in enumerate(group.loops):
            if not loop.hole:
                continue
            rim = samples.boundary_loops[li][::-1]
            chart = LayeredQuadChart(self.patterns, spacing=self.spacing,
                max_cells=self.max_cells - total, quality_iterations=self.quality_iterations)
            try:
                result = chart.build(samples.xy[list(rim)], require_quality=require_quality)
            except QuadChartQualityError as exc:
                if not require_quality:
                    raise
                result = exc.result
                failures.append(dict(hole=li, reason=str(exc)))
            result.mesh.material = 'grass'
            total += len(result.mesh.quads)
            holes[li] = result
        return CompoundChartResult(mesh, holes, partition, samples, boundary_refs,
            dict(planar_only=True, members=[key.text() for key in group.members], blocks=reports,
                 quality=quality, total_quads=total, holes=len(holes), area_error=partition.area_error,
                 initial_chart_failures=failures,
                 construction='source boundary ownership cuts with direct layered quad charts'))
