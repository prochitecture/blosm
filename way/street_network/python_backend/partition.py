"""Source-aware planar partitions of compound junctions, without triangulation.

Cuts join boundary ownership changes of neighboring source junctions. Derived
seams carry no prescribed carriageway profile: the compound height field must
be reconciled from its external conditions. Every original source ID remains.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass, replace
import math
import numpy as np
from .mesh import Assembler, Mesh, cross, publication_coordinates
from .model import Key
from .quad_charts import LayeredQuadChart, visibility_center
from .regions import _Nodes, _classify, signed_area
from .compounds import rim_label
from .source_exterior_sampling import chart_sampling_scope, chart_sampling_steps, try_primary_sampling
from .updates import drain_steps


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
        self.visibility_source_proof = None
        try:
            candidates = tuple(self._cuts(group, context))
        except ValueError as failure:
            if (str(failure) != 'Neighboring junctions need another cut policy: no unique pair of boundary ownership changes'
                    or context is None or group.physical_layout is not None or not group.claims
                    or len(group.loops) != 1 or group.loops[0].hole):
                raise
            from .compound_ports import CompoundPortPlanner
            if CompoundPortPlanner.incomplete_ports(group):
                raise
            from .source_strip_containment import certify
            proof = certify(context, group, {claim.section_id for claim in group.claims})
            candidates = tuple(self._physical_cuts(group))
            self.visibility_source_proof = proof
        merged = {}
        for cut in candidates:
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
            if group.physical_layout is not None and len(group.loops) > 1:
                from .physical_hole_cuts import PhysicalHoleCuts
                yield from PhysicalHoleCuts().cuts(group, context)
                return
            if (len(group.loops) != 1 or group.loops[0].hole
                    or (group.physical_layout is None and len(group.members) != 1)):
                raise ValueError('Non-source compound cuts need another region policy')
            yield from self._physical_cuts(group)
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

    @staticmethod
    def _physical_cuts(group):
        from .decomposition import VisibilitySplitter
        loop = group.loops[0]
        sources = {p.region.owner: p for p in group.footprints}
        representative = {}
        for i, span in enumerate(loop.spans):
            label = rim_label(span, sources)
            if label.kind == 'curb':
                token = span.owner, label.index
                distance = abs(.5*(label.first+label.last)-.5)
                if token not in representative or distance < representative[token][0]:
                    representative[token] = distance, i
        preferred = [i for distance, i in representative.values()]
        try:
            diagonals = VisibilitySplitter().split(loop.xy, preferred=preferred)
        except ValueError as failure:
            if str(failure) != 'Visibility decomposition needs another compatible cut family':
                raise
            diagonals = VisibilitySplitter(max_candidates=3072).split(loop.xy, preferred=preferred)
        for first, last in diagonals:
            yield JunctionCut(None, group.members[0], group.members[0],
                              BoundaryRef(0, first, 0.), BoundaryRef(0, last, 0.))


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
    def __init__(self, policy=None, *, epsilon=1e-7, max_nodes=16384, max_cuts=256, require_visibility=True):
        if (not math.isfinite(epsilon) or epsilon <= 0. or type(max_nodes) is not int or max_nodes < 4
                or type(max_cuts) is not int or max_cuts < 1):
            raise ValueError('Invalid compound partition limits')
        self.policy = BoundaryOwnershipCuts() if policy is None else policy
        self.require_visibility = bool(require_visibility)
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
        # Fractional cut endpoints can change the external loop parity.
        # Resolve this after every cut is known, on an owned curb only.
        for li, loop in enumerate(loops):
            counts = [sum((hi-lo)*np.linalg.norm(loop[(ei+1)%len(loop)]-loop[ei]) > epsilon
                          for lo,hi in zip(sorted(set(values)), sorted(set(values))[1:]))
                      for ei,values in enumerate(fractions[li])]
            if sum(counts) % 2:
                candidates = [(np.linalg.norm(loop[(ei+1)%len(loop)]-loop[ei])*(hi-lo),ei,lo,hi)
                    for ei,values in enumerate(fractions[li])
                    if rim_label(group.loops[li].spans[ei], sources).kind == 'curb'
                    for lo,hi in zip(sorted(set(values)),sorted(set(values))[1:])]
                if not candidates or max(candidates)[0] <= 4.*epsilon:
                    raise ValueError('No owned curb can resolve cut-endpoint parity')
                _,ei,lo,hi=max(candidates)
                fractions[li][ei].append(.5*(lo+hi))
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
            center = visibility_center(xy) if self.require_visibility else xy.mean(axis=0)
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
    def __init__(self, message, *, candidate=None, details=None):
        self.candidate, self.details = candidate, details
        super().__init__(message)


class GrassChartGradingError(QuadChartGradingError):
    """A qualified physical pavement family still has an ungraded grass disk."""
    def __init__(self, message, *, family=None, candidate=None, details=None):
        self.family = family
        super().__init__(message, candidate=candidate, details=details)


class CompoundQuadChart:
    """Conforming planar pavement and hole charts over shared source cuts.

    Output is preparation geometry at Z=0. External profile conditions and the
    compound height field must be applied before a surface provider publishes it.
    """
    def __init__(self, patterns, *, spacing=.7, partitioner=None, max_cells=50000,
                 quality_iterations=300, four_sided=False, piecewise_sampling=False, recover_grass=False,
                 source_height_admission=None):
        if source_height_admission is not None and not callable(source_height_admission):
            raise ValueError('Source-height chart admission requires an optional callable')
        self.patterns, self.spacing = patterns, spacing
        self.partitioner = CompoundPartitioner() if partitioner is None else partitioner
        self.max_cells, self.quality_iterations = max_cells, quality_iterations
        self.four_sided = bool(four_sided)
        self.piecewise_sampling = bool(piecewise_sampling)
        self.recover_grass = bool(recover_grass)
        self.source_height_admission = source_height_admission

    def _fibre_hole(self,boundary,max_cells):
        from .hole_fibre_chart import build_fibre_hole
        return build_fibre_hole(self.patterns,boundary,self.spacing,max_cells)

    def _ordinary_source_prefix_steps(self,group,context,previous):
        """Last grading retry for an independently admitted original root union."""
        from .ordinary_source_prefix_atlas import eligible,native_eligible,build_steps
        if native_eligible(group) or eligible(group):
            return (yield from build_steps(self,group,context,previous))
        return None

    def _piecewise_ordinary(self, group, context, previous):
        """Retry a rejected ordinary rim on the same lead/arc error budget."""
        if (self.piecewise_sampling or self.four_sided or len(group.members) != 1
                or len(group.footprints) != 1 or len(group.loops) != 1
                or group.loops[0].hole or group.claims or group.physical_layout is not None):
            raise ValueError('Piecewise ordinary retry needs one unchanged local junction footprint')
        result = CompoundQuadChart(self.patterns, spacing=self.spacing, partitioner=self.partitioner,
            max_cells=self.max_cells, quality_iterations=self.quality_iterations,
            piecewise_sampling=True).build(group, context)
        result.metadata['piecewise_curb_sampling'] = dict(previous_chart_failure=previous,
            original_analytical_curves=True, incoming_profiles_unchanged=True,
            curve_tolerance=context.network.settings.curve_tolerance,
            cell_budget=self.max_cells)
        return result

    @staticmethod
    def _composed_scope(group):
        layout = group.physical_layout
        if layout is None or not getattr(layout, 'finite_scope_source_proof', None):
            return False
        from .finite_scope_compounds import ExpandedScopePort
        return any(isinstance(port, ExpandedScopePort) for port in getattr(layout, 'finite_scope_ports', ()))

    def _repairable_leaf(self, group):
        if self._composed_scope(group):
            return True
        from .physical_hole_cuts import PhysicalHoleCuts
        policy = self.partitioner.policy
        return (isinstance(policy, PhysicalHoleCuts) and policy.allow_filled
                and policy.filled_eligible(group))

    def _filled_source(self, group, context, sampled_group, previous, *, source_corridors=True):
        return drain_steps(self._filled_source_steps(group,context,sampled_group,previous,source_corridors=source_corridors))

    def _filled_source_steps(self, group, context, sampled_group, previous, *, source_corridors=True):
        """Build a four-sided atlas over the certified, filled source region."""
        from .physical_hole_cuts import PhysicalHoleCuts
        from .source_strip_containment import certify
        if not PhysicalHoleCuts.filled_eligible(group) or not group.claims:
            raise ValueError('Filled-source cuts require a certified complete source union')
        proof = certify(context, sampled_group, {claim.section_id for claim in group.claims})
        result = yield from CompoundQuadChart(self.patterns, spacing=self.spacing,
            partitioner=CompoundPartitioner(policy=PhysicalHoleCuts(spacing=self.spacing,
                source_corridors=source_corridors, allow_filled=True), require_visibility=False),
            max_cells=self.max_cells, quality_iterations=min(100, self.quality_iterations),
            four_sided=True).build_steps(group, context, require_quality=True)
        from .four_side_charts import quality
        stored = quality(result.mesh.with_positions(publication_coordinates(result.mesh.xyz, 1., physical=True)))
        if stored['maximum_edge_ratio'] > 8. or stored['minimum_corner_sine'] < .03:
            raise ValueError('Filled-source pavement atlas fails stored-coordinate quality')
        result.metadata = dict(result.metadata, filled_source_partition=True,
            source_width_visibility=proof, stored_quality=stored, previous_chart_failure=previous)
        return result

    def _longitudinal(self, group, context, previous):
        return drain_steps(self._longitudinal_steps(group, context, previous))

    def _longitudinal_steps(self, group, context, previous):
        from .longitudinal_cuts import LongitudinalCuts
        failures = []
        attempts = tuple(dict.fromkeys((min(20, self.quality_iterations), self.quality_iterations)))
        for iterations in attempts:
            try:
                chart = CompoundQuadChart(self.patterns, spacing=self.spacing,
                    partitioner=CompoundPartitioner(policy=LongitudinalCuts(self.spacing), require_visibility=False),
                    max_cells=self.max_cells, quality_iterations=iterations,
                    four_sided=True, piecewise_sampling=self.piecewise_sampling)
                result = yield from chart._build_steps(group, context, require_quality=True)
            except ValueError as failure:
                if iterations == attempts[-1]:
                    raise
                failures.append(str(failure))
                continue
            result.metadata = dict(result.metadata, longitudinal_partition=True,
                previous_chart_failures=previous, previous_longitudinal_failures=failures,
                local_quality_iterations=iterations)
            return result

    def _reflex(self, group, context, previous):
        return drain_steps(self._reflex_steps(group, context, previous))

    def _reflex_steps(self, group, context, previous):
        from .reflex_cuts import ReflexCurbCuts
        chart = CompoundQuadChart(self.patterns, spacing=self.spacing,
            partitioner=CompoundPartitioner(policy=ReflexCurbCuts(), require_visibility=False),
            max_cells=self.max_cells, quality_iterations=min(100, self.quality_iterations),
            four_sided=True, piecewise_sampling=self.piecewise_sampling)
        result = yield from chart._build_steps(group, context, require_quality=True)
        result.metadata = dict(result.metadata, reflex_curb_partition=True, previous_chart_failures=previous)
        return result

    @chart_sampling_scope
    def build(self, group, context, *, require_quality=True):
        return drain_steps(self.build_steps(group,context,require_quality=require_quality))

    @chart_sampling_steps
    def build_steps(self, group, context, *, require_quality=True):
        """Return control between footprint candidates, with the same ordering and checks."""
        try:
            return (yield from self._build_with_recovery(group, context, require_quality=require_quality))
        except QuadChartGradingError as error:
            from .ordinary_half_charts import applicable, recover
            if self.four_sided or not require_quality or not applicable(group):
                raise
            yield
            return recover(self, group, context, error)

    def _build_with_recovery(self, group, context, *, require_quality):
        from .physical_hole_cuts import PhysicalHoleCuts
        if require_quality and not self.four_sided:
            if getattr(group.physical_layout, 'native_original_circle_fillets', ()):
                from .bank_circle_preparation import corridor_chart_steps
                try:
                    result = yield from corridor_chart_steps(self, group, context)
                except QuadChartGradingError:
                    result = None
                if result is not None:
                    return result
            from .source_pool_atlas import admission as pool_admission, build as build_large_pool
            certificate=pool_admission(group,context)
            if certificate is not None:
                return build_large_pool(self,group,context,certificate)
        if not self.four_sided and PhysicalHoleCuts.eligible(group):
            source_corridors = bool(getattr(group.physical_layout, 'source_cycle_proofs', None))
            previous = []
            grass_failure = None
            # A source corridor and the principal physical frame are bounded
            # computational alternatives over the same exact physical rims.
            alternatives = [(True, False, False), (False, False, False)] if source_corridors else [(False, False, False)]
            fragments = getattr(group.physical_layout, 'internal_profile_fragments', {})
            source_aggregate = (getattr(group.physical_layout, 'source_region_union', False)
                and isinstance(fragments, dict) and bool(fragments.get('complete_source_width')))
            if hasattr(group.physical_layout, 'source_cycle_proofs') or source_aggregate:
                alternatives.append((True, True, False))
            alternatives.append((False, False, True))
            if getattr(group.physical_layout, 'source_cycle_full_width_proof', None):
                alternatives.append((True, False, True))
            for attempt, (use_sources, hole_bridges, existing_vertices) in enumerate(alternatives):
                yield
                try:
                    trial = CompoundQuadChart(self.patterns, spacing=self.spacing,
                        partitioner=CompoundPartitioner(policy=PhysicalHoleCuts(spacing=self.spacing,
                            source_corridors=use_sources, hole_bridges=hole_bridges,
                            existing_vertices=existing_vertices), require_visibility=False),
                        max_cells=self.max_cells, quality_iterations=min(100, self.quality_iterations),
                        four_sided=True)
                    result = yield from trial.build_steps(group, context, require_quality=require_quality)
                    if require_quality:
                        from .four_side_charts import quality as stored_quality
                        stored = stored_quality(result.mesh.with_positions(
                            publication_coordinates(result.mesh.xyz, 1., physical=True)))
                        if stored['maximum_edge_ratio'] > 8. or stored['minimum_corner_sine'] < .03:
                            raise ValueError('Physical-hole pavement atlas fails stored-coordinate quality')
                except ValueError as failure:
                    if isinstance(failure, GrassChartGradingError) and failure.family is not None:
                        grass_failure = failure
                    if attempt==0 and require_quality:
                        yield
                        sampled=try_primary_sampling(self,group,context,failure,'first physical-hole chart')
                        if sampled is not None:return sampled
                    if attempt == len(alternatives)-1:
                        yield
                        from .physical_hole_recovery import build_steps as recover_source_normals
                        try:
                            return (yield from recover_source_normals(self,group,context,
                                [*previous,dict(source_corridors=use_sources,hole_bridges=hole_bridges,
                                    existing_vertices=existing_vertices,reason=str(failure))]))
                        except ValueError as recovery_failure:
                            last_recovery_failure = recovery_failure
                            if isinstance(recovery_failure, GrassChartGradingError) and recovery_failure.family is not None:
                                grass_failure = recovery_failure
                            if grass_failure is not None and require_quality:
                                yield
                                from .hole_banded_chart import retry_failed_family
                                try:
                                    return retry_failed_family(self, group, context, grass_failure,
                                        [*previous, dict(source_corridors=use_sources, hole_bridges=hole_bridges,
                                            existing_vertices=existing_vertices, reason=str(failure))], recovery_failure)
                                except ValueError as grass_recovery_failure:
                                    last_recovery_failure = grass_recovery_failure
                            from .paired_interval_atlas import eligible as paired_eligible, build_steps as recover_paired_intervals
                            if paired_eligible(group):
                                yield
                                try:
                                    return (yield from recover_paired_intervals(self,group,context,
                                        [*previous,dict(reason=str(failure)),dict(reason=str(last_recovery_failure))]))
                                except ValueError as paired_failure:
                                    last_recovery_failure = paired_failure
                            from .source_cycle_atlas import eligible as cycle_eligible, build as recover_original_corridors
                            if cycle_eligible(group,context):
                                yield
                                try:
                                    return recover_original_corridors(self,group,context,
                                        [*previous,dict(reason=str(failure)),dict(reason=str(last_recovery_failure))])
                                except ValueError as cycle_failure:
                                    last_recovery_failure = cycle_failure
                            from .owned_bundle_normal_atlas import eligible as owned_normals_eligible,build as recover_owned_normals
                            if owned_normals_eligible(group,context):
                                yield
                                try:
                                    return recover_owned_normals(self,group,context,
                                        [*previous,dict(reason=str(failure)),dict(reason=str(last_recovery_failure))])
                                except ValueError as owned_normal_failure:
                                    last_recovery_failure=owned_normal_failure
                            if require_quality:
                                from .source_exterior_sampling import eligible as sampling_eligible, recover_sampling
                                if sampling_eligible(group,context):
                                    return recover_sampling(self,group,context,failure)
                            raise failure from last_recovery_failure
                    previous.append(dict(source_corridors=use_sources, hole_bridges=hole_bridges,
                                         existing_vertices=existing_vertices, reason=str(failure)))
                    continue
                result.metadata = dict(result.metadata, physical_hole_partition=True,
                    source_corridors=use_sources, source_hole_bridges=hole_bridges,
                    existing_cut_vertices=existing_vertices,
                    previous_physical_atlas_failures=previous)
                if require_quality:
                    result.metadata['stored_quality'] = stored
                return result
        if not self.four_sided and self._composed_scope(group):
            from .finite_scope_compounds import ScopePlaneCuts
            result = CompoundQuadChart(self.patterns, spacing=self.spacing,
                partitioner=CompoundPartitioner(policy=ScopePlaneCuts(), require_visibility=False),
                max_cells=self.max_cells, quality_iterations=self.quality_iterations,
                four_sided=True).build(group, context, require_quality=require_quality)
            result.metadata = dict(result.metadata, composed_scope_partition=True)
            return result
        if require_quality and not self.four_sided:
            from .longitudinal_cuts import LongitudinalCuts
            from .reflex_cuts import ReflexCurbCuts
            longitudinal = LongitudinalCuts.eligible(group, self.spacing)
            reflex = ReflexCurbCuts.eligible(group)
            if longitudinal or reflex:
                # Retain a radial result already qualified by local repair.
                # Final edge rotations and global repair can be avoided when
                # a direct partition qualifies the same physical boundary.
                try:
                    probe = yield from self._build_steps(group, context, require_quality=False)
                    if probe.metadata['quality']['maximum_edge_ratio'] <= 8.:
                        return probe
                    previous = [dict(reason='Radial seed before final repair exceeds edge-ratio limit',
                                     quality=probe.metadata['quality'])]
                except ValueError as failure:
                    previous = [dict(reason=str(failure))]
                for eligible, construct in ((longitudinal, self._longitudinal), (reflex, self._reflex)):
                    if not eligible:
                        continue
                    try:
                        return construct(group, context, previous)
                    except ValueError:
                        pass  # The existing bounded construction remains available.
        from .quad_charts import QuadChartCellBudgetError,QuadChartLayerBudgetError,QuadChartCoreError,VisibilityKernelError
        from .decomposition import VisibilityCutFamilyError
        from .sampling import AnalyticalContactError
        try:
            result = yield from self._build_steps(group, context, require_quality=require_quality)
        except AnalyticalContactError as original:
            if self.four_sided or not require_quality:
                raise
            from .source_contact_recovery import recover
            try:
                return recover(self,group,context,original)
            except ValueError as recovery_failure:
                raise original from recovery_failure
        except VisibilityKernelError as original:
            if self.four_sided or not require_quality:
                raise
            from .visibility_recovery import recover_kernel_steps
            try:
                return (yield from recover_kernel_steps(self, group, context, original))
            except ValueError as recovery_failure:
                raise original from recovery_failure
        except VisibilityCutFamilyError as original:
            if self.four_sided or not require_quality:
                raise
            sampled=try_primary_sampling(self,group,context,original,'primary visibility partition')
            if sampled is not None:return sampled
            from .visibility_recovery import recover_splitter_steps
            try:
                return (yield from recover_splitter_steps(self,group,context,original))
            except ValueError as recovery_failure:
                if recovery_failure is original:raise
                raise original from recovery_failure
        except (QuadChartCellBudgetError,QuadChartLayerBudgetError,QuadChartCoreError) as original:
            if self.four_sided or not require_quality:
                raise
            if getattr(group.physical_layout,'original_lead_joins',()) and PhysicalHoleCuts.filled_eligible(group):
                # The short auxiliary lead can concentrate radial layers on
                # one side of a large compound. Try the principal-axis cuts
                # with the exact same source, seams and cell/quality limits.
                try:
                    return (yield from self._filled_source_steps(group,context,group,str(original),source_corridors=False))
                except ValueError:
                    pass
            sampled=try_primary_sampling(self,group,context,original,'primary chart budget')
            if sampled is not None:return sampled
            # A finite radial-layer budget is not a physical boundary defect.
            # The existing reflex cuts can also help ordinary multi-arm rims
            # below their cheap preselection threshold. Try that same bounded
            # construction after a typed layer failure, preserving every rim.
            try:
                result = self._reflex(group, context, [dict(reason=str(original))])
            except ValueError as alternative:
                from .compound_chart_recovery import recover
                try:
                    return recover(self,group,context,[dict(reason=str(original)),dict(reason=str(alternative))])
                except ValueError:
                    pass
                try:
                    return self._piecewise_ordinary(group, context, str(original)+'; reflex: '+str(alternative))
                except ValueError as sampling_failure:
                    from .source_exterior_sampling import eligible as sampling_eligible, recover_sampling
                    if sampling_eligible(group,context):
                        return recover_sampling(self,group,context,original)
                    from .ordinary_source_prefix_atlas import native_eligible
                    if native_eligible(group):
                        prefix = yield from self._ordinary_source_prefix_steps(group,context,str(original))
                        if prefix is not None:return prefix
                    raise original from sampling_failure
            name=('budget_reflex_fallback' if isinstance(original,(QuadChartCellBudgetError,QuadChartLayerBudgetError))
                  else 'radial_core_reflex_fallback')
            result.metadata[name] = dict(radial_failure=str(original),
                physical_boundary_changed=False, cell_budget=self.max_cells)
            return result
        if not require_quality or not result.metadata.get('initial_chart_failures'):
            return result
        # Elongated forks need longitudinal charts rather than radial fans.
        # This bounded alternative is tried only after the original chart has
        # failed its quality contract; successful existing charts are unchanged.
        if not self.four_sided:
            try:
                return self._longitudinal(group, context, result.metadata['initial_chart_failures'])
            except ValueError as failure:
                result.metadata['longitudinal_chart_failure'] = str(failure)
            try:
                return self._reflex(group, context, result.metadata['initial_chart_failures'])
            except ValueError as failure:
                result.metadata['reflex_chart_failure'] = str(failure)
        # A partition cut is computational: after assembly its vertices are
        # interior degrees of freedom. Keeping that cut fixed during each
        # leaf's repair can prevent an otherwise valid global quad chart.
        # Reuse the first positive charts; do not repeat their construction or
        # proxy relaxation merely to suppress a final leaf quality exception.
        # Actual physical boundaries remain exact.
        try:
            mesh, details = yield from self._correct_steps(result.mesh, retry_rotations=(self.four_sided
                and group.physical_layout is not None and len(group.loops) > 1))
        except QuadChartGradingError as original:
            if not self.four_sided:
                from .source_exterior_sampling import tiny_auxiliary_edges
                failed_cells = ((original.details or {}).get('final') or {}).get('cells_above_ratio_limit',0)
                if (failed_cells > 256 and original.candidate is not None
                        and PhysicalHoleCuts.filled_eligible(group) and group.claims
                        and tiny_auxiliary_edges(group,context.network.settings)):
                    # A short generated contact can stall a global line search.
                    # Repair the existing positive chart before rebuilding its
                    # source-proven perimeter with a different sampling phase.
                    from .local_quad_grading import grade_dense_steps
                    try:
                        mesh, repair = yield from grade_dense_steps(original.candidate,self.max_cells)
                    except ValueError as failure:
                        original.details = dict(original.details, independent_vertex_failure=str(failure))
                    else:
                        result.mesh = mesh
                        result.metadata = dict(result.metadata, quality=mesh.validate(),
                            total_quads=len(mesh.quads), assembled_quality_correction=repair,
                            previous_chart_failure=str(original), block_metadata_phase='BEFORE_ASSEMBLED_CORRECTION')
                        return result
                sampled=try_primary_sampling(self,group,context,original,'primary assembled grading')
                if sampled is not None:return sampled
            if not self.four_sided and PhysicalHoleCuts.filled_eligible(group):
                try:
                    return (yield from self._filled_source_steps(group, context, result.partition.group, str(original)))
                except ValueError as failure:
                    original.details = dict(original.details or {}, filled_source_failure=str(failure))
            from .corner_seed import CornerBandSeed
            try:
                alternative = CornerBandSeed(self.patterns, max_cells=self.max_cells).build(group, context, self.partitioner)
            except ValueError as failure:
                try:
                    return self._piecewise_ordinary(group, context, str(original)+'; corner seed: '+str(failure))
                except ValueError as sampling_failure:
                    if not self.four_sided:
                        from .compound_chart_recovery import recover_longitudinal_steps
                        try:
                            return (yield from recover_longitudinal_steps(self,group,context,
                                [dict(reason=str(original)),dict(reason=str(failure)),
                                 dict(reason=str(sampling_failure))]))
                        except ValueError:
                            pass
                        from .compound_chart_recovery import recover_source_reflex
                        try:
                            return recover_source_reflex(self,group,context,
                                [dict(reason=str(original)),dict(reason=str(failure)),
                                 dict(reason=str(sampling_failure))])
                        except ValueError:
                            pass
                        from .source_exterior_sampling import eligible as sampling_eligible, recover_sampling
                        if sampling_eligible(group,context):
                            return recover_sampling(self,group,context,original)
                    previous=f'{original}; corner seed: {failure}; piecewise sampling: {sampling_failure}'
                    try:
                        prefix = yield from self._ordinary_source_prefix_steps(group,context,previous)
                    except ValueError as prefix_failure:
                        raise QuadChartGradingError(previous+'; ordinary source prefix: '+str(prefix_failure)) from original
                    if prefix is not None:return prefix
                    raise QuadChartGradingError(previous) from original
            alternative.metadata = dict(alternative.metadata, previous_chart_failure=str(original))
            return alternative
        result.mesh = mesh
        for li, hole in tuple(result.holes.items()):
            if (self.four_sided and group.physical_layout is not None and len(group.loops) > 1
                    and hole.mesh.validate()['maximum_edge_ratio'] <= 8.):
                stored = hole.mesh.with_positions(publication_coordinates(hole.mesh.xyz, 1., physical=True)).validate()
                if stored['maximum_edge_ratio'] <= 8.:
                    hole.metadata = dict(hole.metadata, stored_quality=stored,
                        grass_quality_contract='positive corners and edge ratio <=8')
                    continue
            try:
                hole.mesh, correction = yield from self._correct_steps(hole.mesh, retry_rotations=(self.four_sided
                    and group.physical_layout is not None and len(group.loops) > 1))
            except QuadChartGradingError as original:
                if not (self.four_sided and group.physical_layout is not None and len(group.loops) > 1):
                    raise
                refinement = result.metadata.get('separator_rim_refinement') or {}
                if refinement.get('retained_qualified_rims') and original.candidate is not None:
                    from .local_quad_grading import grade
                    remaining = self.max_cells - len(mesh.quads) - sum(len(other.mesh.quads)
                        for index, other in result.holes.items() if index != li)
                    for narrow in (False, True):
                        try:
                            disk, correction = grade(original.candidate, remaining, narrow_grass_disk=narrow)
                        except ValueError:
                            continue
                        hole.mesh = disk
                        hole.metadata = dict(hole.metadata, assembled_quality_correction=correction,
                                             quality=disk.validate())
                        break
                    else:
                        disk = None
                    if disk is not None:
                        result.metadata['total_quads'] = len(mesh.quads) + sum(len(h.mesh.quads) for h in result.holes.values())
                        continue
                if not getattr(self, 'recover_grass', False):
                    from .physical_hole_cuts import PhysicalHoleCuts
                    policy = self.partitioner.policy
                    family = ((policy.source_corridors, policy.hole_bridges, policy.existing_vertices)
                              if isinstance(policy, PhysicalHoleCuts) else None)
                    raise GrassChartGradingError(str(original), family=family,
                        candidate=original.candidate, details=original.details) from original
                from .hole_banded_chart import recover_failed_grass
                try:
                    mesh = recover_failed_grass(self, result, li, original)
                except ValueError as failure:
                    original.details = dict(original.details or {}, banded_grass_failure=str(failure))
                    raise original from failure
                continue
            hole.metadata = dict(hole.metadata, assembled_quality_correction=correction,
                                 quality=hole.mesh.validate())
        result.metadata = dict(result.metadata, quality=mesh.validate(),
            assembled_quality_correction=details,
            block_metadata_phase='BEFORE_ASSEMBLED_CORRECTION')
        if details.get('positive_assembled_repair'):
            result.metadata['positive_assembled_repair'] = details['positive_assembled_repair']
        return result

    def _correct(self, mesh, *, retry_rotations=False):
        return drain_steps(self._correct_steps(mesh, retry_rotations=retry_rotations))

    def _correct_steps(self, mesh, *, retry_rotations=False):
        from .quad_rotations import QuadEdgeRotations
        from .quad_charts import FixedRimQuality
        from .quality_correction import QuadShapeCorrection
        correction = QuadShapeCorrection(max_cells=self.max_cells)
        _, initial, _ = correction.evaluate(mesh.xyz[:, :2], mesh.quads)
        if correction.qualified(initial):
            return mesh, dict(qualified=True, iterations=0, initial=initial, final=initial,
                              exact_fixed_rim=True, stop='already_qualified')
        seed, rotations = yield from QuadEdgeRotations(max_passes=8, max_cells=self.max_cells).improve_steps(
            mesh, FixedRimQuality(iterations=self.quality_iterations, power=2))
        candidate, details = yield from correction.improve_steps(seed)
        if (retry_rotations and not details['qualified'] and details['final'] is not None
                and details['final']['cells_above_ratio_limit'] <= 8
                and details['final']['maximum_edge_ratio'] <= 10.):
            rotated, final_rotations = yield from QuadEdgeRotations(max_passes=4, max_cells=self.max_cells).improve_steps(
                candidate, FixedRimQuality(iterations=min(100, self.quality_iterations), power=2))
            improved, final_correction = yield from correction.improve_steps(rotated)
            if final_correction['qualified']:
                return improved, dict(final_correction, seed_rotations=rotations,
                    previous_correction=details, final_rotations=final_rotations)
        if not details['qualified']:
            try:
                repaired, repair = yield from self._repair_positive_candidate_steps(candidate)
            except ValueError as failure:
                details = dict(details, bounded_repair_failure=str(failure))
            else:
                return repaired, dict(qualified=True, previous_correction=details,
                    seed_rotations=rotations, positive_assembled_repair=repair, exact_fixed_rim=True)
            raise QuadChartGradingError('Assembled direct quad chart needs better grading: correction exhausted its finite budget',
                candidate=candidate, details=dict(details, seed_rotations=rotations))
        return candidate, dict(details, seed_rotations=rotations)

    def _cell_collars_steps(self, mesh):
        """Repair a small number of positive cells with all old vertices fixed."""
        bad = np.zeros(len(mesh.quads), bool)
        for xy in (mesh.xyz[:, :2], publication_coordinates(mesh.xyz, 1., physical=True)[:, :2]):
            points = xy[mesh.quads]
            edge = np.roll(points, -1, axis=1)-points
            lengths = np.linalg.norm(edge, axis=-1)
            jac = cross(edge, np.roll(edge, -1, axis=1))
            denominator = lengths*np.roll(lengths, -1, axis=1)
            sine = np.divide(jac, denominator, out=np.zeros_like(jac), where=denominator > 0.)
            ratio = np.divide(lengths.max(axis=1), lengths.min(axis=1),
                out=np.full(len(lengths), np.inf), where=lengths.min(axis=1) > 0.)
            bad |= (ratio > 8.) | (jac.min(axis=1) <= 1e-10) | (sine.min(axis=1) < .03)
        if not 1 <= np.count_nonzero(bad) <= 64:
            raise ValueError('Fixed-cell collar repair exceeds its candidate budget')
        from .finite_scope_charts import grade_fixed_cells_steps
        return (yield from grade_fixed_cells_steps(mesh, self.patterns, max_cells=self.max_cells))

    def _repair_positive_candidate(self, mesh):
        return drain_steps(self._repair_positive_candidate_steps(mesh))

    def _repair_positive_candidate_steps(self, mesh):
        try:
            candidate, mapping, details = yield from self._cell_collars_steps(mesh)
        except ValueError:
            from .interior_quad_refinement import InteriorQuadRefinement
            candidate, mapping, details = yield from InteriorQuadRefinement(self.patterns,
                max_cells=self.max_cells).improve_steps(mesh)
        if len(mapping) != len(mesh.xyz) or len(np.unique(mapping)) != len(mapping):
            raise ValueError('Positive chart repair lost original vertex identity')
        # Keep every pre-existing vertex ID so external profile aliases and
        # boundary references need no reinterpretation after topology repair.
        extra = np.flatnonzero(~np.isin(np.arange(len(candidate.xyz)), mapping))
        order = np.r_[mapping, extra]
        inverse = np.empty(len(order), dtype=np.int32)
        inverse[order] = np.arange(len(order), dtype=np.int32)
        result = Mesh(candidate.xyz[order], inverse[candidate.quads], candidate.material,
            {key: inverse[ids] for key, ids in candidate.seams.items()}, dict(candidate.uv_layers))
        from .four_side_charts import quality
        report = quality(result)
        stored = quality(result.with_positions(publication_coordinates(result.xyz, 1., physical=True)))
        if (max(report['maximum_edge_ratio'], stored['maximum_edge_ratio']) > 8.
                or min(report['minimum_corner_sine'], stored['minimum_corner_sine']) < .03):
            raise ValueError('Positive chart repair fails physical or stored quality')
        return result, dict(details, quality=report, stored_quality=stored,
                            original_vertex_ids_preserved=True)

    def _build(self, group, context, *, require_quality):
        return drain_steps(self._build_steps(group,context,require_quality=require_quality))

    def _build_steps(self, group, context, *, require_quality):
        from .sampling import AnalyticalBoundarySampler
        from .quad_charts import QuadChartQualityError
        sampler = AnalyticalBoundarySampler(tolerance=context.network.settings.curve_tolerance,
            spacing=self.spacing, piecewise=self.piecewise_sampling)
        group = sampler.build(group)
        yield
        partition = self.partitioner.build(group, context)
        if self.four_sided:
            from .longitudinal_cuts import even_sided_partition
            partition=even_sided_partition(partition, max_nodes=self.partitioner.max_nodes)
        samples = partition.sample(self.spacing)
        assembly, reports, failures, total = Assembler(), [], [], 0
        for index, (block, rim) in enumerate(zip(partition.blocks, samples.rims)):
            yield
            chart = LayeredQuadChart(self.patterns, spacing=self.spacing,
                max_cells=self.max_cells - total, quality_iterations=self.quality_iterations)
            try:
                if self.four_sided:
                    from .four_side_charts import FourSidedQuadChart
                    try:
                        result=FourSidedQuadChart(self.patterns,max_cells=self.max_cells-total,
                            quality_iterations=self.quality_iterations).build(samples.xy[list(rim)])
                    except ValueError as four_failure:
                        result = None
                        if self._repairable_leaf(group):
                            try:
                                from .finite_scope_charts import fixed_rim_chart_steps
                                result = yield from fixed_rim_chart_steps(samples.xy[list(rim)], self.patterns,
                                    max_cells=self.max_cells-total)
                            except ValueError:
                                pass
                        if result is None:
                            try:
                                result = yield from chart.build_steps(samples.xy[list(rim)], require_quality=True)
                            except QuadChartQualityError as radial_failure:
                                if (not self._repairable_leaf(group)
                                        and (group.physical_layout is None or len(group.loops) == 1)):
                                    raise ValueError(f'No qualified direct leaf chart: {four_failure}; {radial_failure}') from radial_failure
                                # Computational cuts remain interior degrees
                                # of freedom after assembly. The full chart
                                # must still pass its final shape correction.
                                result = radial_failure.result
                                failures.append(dict(block=index, reason=f'{four_failure}; {radial_failure}'))
                else:
                    result = yield from chart.build_steps(samples.xy[list(rim)], require_quality=require_quality)
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
        separator_refinement = None
        if self.four_sided and group.physical_layout is not None and len(group.loops) > 1:
            from .hole_rim_refinement import grade_hole_mesh_steps
            if require_quality:
                try:
                    yield
                    mesh, pavement_seed_correction = yield from self._correct_steps(mesh, retry_rotations=True)
                except QuadChartGradingError as failure:
                    from .physical_hole_cuts import PhysicalHoleCuts
                    policy = self.partitioner.policy
                    if (failure.candidate is None or len(group.loops) < 3 or
                            not getattr(group.physical_layout, 'source_cycle_full_width_proof', None) or
                            not isinstance(policy, PhysicalHoleCuts) or policy.source_corridors or
                            policy.existing_vertices or policy.hole_bridges):
                        raise
                    from .local_quad_grading import grade
                    mesh, pavement_seed_correction = grade(failure.candidate, self.max_cells)
            if require_quality and 'local_adjustments' in pavement_seed_correction:
                # This pavement already meets both precision contracts on its
                # original island rims. Let the disk builders consume those
                # rims directly. Any disk requesting more points must still go
                # through the shared conforming attachment and exact seam check.
                separator_refinement = dict(rounds=[], exact_polygon=True, retained_qualified_rims=True)
                fixed_guides = np.empty(0, dtype=np.int32)
            else:
                mesh, boundary_refs, separator_refinement, fixed_guides = yield from grade_hole_mesh_steps(mesh, boundary_refs, group,
                    self.patterns, max_cells=self.max_cells)
            if require_quality:
                from .quality_correction import QuadShapeCorrection
                mesh, local_correction = yield from QuadShapeCorrection(max_cells=self.max_cells).improve_steps(mesh,
                    fixed_vertices=fixed_guides)
                local_refinements = []
                attachment_recovery=getattr(self,'recover_separator_attachment',None)
                attachment_seed=(mesh.with_positions(mesh.xyz.copy()),dict(boundary_refs),fixed_guides.copy()) if attachment_recovery is not None else None
                for _ in range(2):
                    if local_correction['qualified']:
                        break
                    points = mesh.xyz[mesh.quads, :2]
                    lengths = np.linalg.norm(np.roll(points, -1, axis=1)-points, axis=2)
                    selected = np.flatnonzero(lengths.max(axis=1)/lengths.min(axis=1) > 8.)
                    if not len(selected) or len(selected) > 64:
                        break
                    from .refinement import QuadRefiner
                    mesh, remap = yield from QuadRefiner(self.patterns, max_cells=self.max_cells).refine_steps(mesh, selected,
                        lambda xy: np.zeros(len(xy)))
                    boundary_refs = {int(remap[node]): ref for node, ref in boundary_refs.items()}
                    fixed_guides = remap[fixed_guides]
                    mesh, local_correction = yield from QuadShapeCorrection(max_cells=self.max_cells).improve_steps(mesh,
                        fixed_vertices=fixed_guides)
                    local_refinements.append(dict(selected=len(selected), quads=len(mesh.quads),
                        correction=local_correction))
                if not local_correction['qualified']:
                    if attachment_recovery is None:
                        raise ValueError('Separator-rim quad refinement needs better local grading')
                    mesh,boundary_refs,fixed_guides,local_correction=attachment_recovery(
                        *attachment_seed,local_correction,local_refinements)
                separator_refinement = dict(separator_refinement, pavement_seed_correction=pavement_seed_correction,
                    local_correction=local_correction, local_refinements=local_refinements)
            total = len(mesh.quads)
        quality = mesh.validate()
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        if (len(mesh.xyz) - len(unique) + len(mesh.quads) != 2 - len(group.loops)
                or not set(counts).issubset({1, 2})
                or int(np.sum(counts == 1)) != sum(map(len, mesh.seams.values()))):
            raise ValueError('Compound quad charts do not form the required conforming pavement region')
        holes, prescribed_hole_rims = {}, {}
        for li, loop in enumerate(group.loops):
            if not loop.hole:
                continue
            rim = samples.boundary_loops[li][::-1]
            yield
            hole_boundary = mesh.xyz[mesh.seams[li], :2][::-1]
            chart = LayeredQuadChart(self.patterns, spacing=self.spacing,
                max_cells=self.max_cells - total, quality_iterations=self.quality_iterations)
            try:
                if self.four_sided and group.physical_layout is not None:
                    from .four_side_charts import FourSidedQuadChart
                    try:
                        result = FourSidedQuadChart(self.patterns, max_cells=self.max_cells-total,
                            quality_iterations=self.quality_iterations).build(hole_boundary)
                    except ValueError:
                        from .hole_disk_chart import build_hole_disk
                        try:
                            result = build_hole_disk(self.patterns, hole_boundary, self.spacing,
                                self.max_cells-total, self.quality_iterations)
                        except ValueError:
                            try:
                                result = self._fibre_hole(hole_boundary,self.max_cells-total)
                            except ValueError:
                                result = yield from chart.build_steps(hole_boundary, require_quality=require_quality)
                            else:
                                if result.metadata['added_boundary_samples']:
                                    prescribed_hole_rims[li] = result.mesh.xyz[result.mesh.seams[0], :2]
                else:
                    result = yield from chart.build_steps(hole_boundary, require_quality=require_quality)
            except QuadChartQualityError as exc:
                if not require_quality:
                    raise
                result = exc.result
                failures.append(dict(hole=li, reason=str(exc)))
            result.mesh.material = 'grass'
            total += len(result.mesh.quads)
            holes[li] = result
        prescribed_refinement = None
        if prescribed_hole_rims:
            from .conforming_boundary_refinement import refine
            grass_cells = sum(len(hole.mesh.quads) for hole in holes.values())
            mesh, boundary_refs, _, prescribed_refinement = refine(mesh, boundary_refs, group,
                self.patterns, prescribed_hole_rims, max_cells=self.max_cells-grass_cells)
            if require_quality:
                yield
                mesh, correction = yield from self._correct_steps(mesh, retry_rotations=True)
                prescribed_refinement = dict(prescribed_refinement, assembled_correction=correction)
            quality = mesh.validate()
            total = len(mesh.quads)+grass_cells
            edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
            unique, counts = np.unique(edges, axis=0, return_counts=True)
            if (len(mesh.xyz)-len(unique)+len(mesh.quads) != 2-len(group.loops)
                    or not set(counts).issubset({1, 2})
                    or int(np.sum(counts == 1)) != sum(map(len, mesh.seams.values()))):
                raise ValueError('Prescribed separator refinement changed pavement topology')
            for li, hole in holes.items():
                expected = mesh.xyz[mesh.seams[li], :2][::-1]
                rim = hole.mesh.xyz[hole.mesh.seams[0], :2]
                starts = np.flatnonzero(np.all(rim == expected[0], axis=1))
                if len(starts) != 1:
                    raise ValueError('Prescribed separator chart lost its shared physical seam')
                ids = np.roll(hole.mesh.seams[0], -int(starts[0]))
                if not np.array_equal(hole.mesh.xyz[ids, :2], expected):
                    raise ValueError('Prescribed separator and pavement rims do not match exactly')
                hole.mesh.seams[0] = ids
        return CompoundChartResult(mesh, holes, partition, samples, boundary_refs,
            dict(planar_only=True, members=[key.text() for key in group.members], blocks=reports,
                 quality=quality, total_quads=total, holes=len(holes), area_error=partition.area_error,
                 initial_chart_failures=failures,
                 chart_policy=type(self.partitioner.policy).__name__,
                 source_width_visibility=getattr(self.partitioner.policy, 'visibility_source_proof', None),
                 separator_rim_refinement=separator_refinement,
                 prescribed_hole_refinement=prescribed_refinement,
                 **(dict(piecewise_curve_bounds=sampler.curve_reports) if self.piecewise_sampling else {}),
                 construction=('physical-rim cuts with direct four-sided and local radial quad charts' if self.four_sided
                               else 'source boundary ownership cuts with direct layered quad charts')))
