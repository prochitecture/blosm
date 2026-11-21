"""Computational bands for a prescribed simple grass disk.

Keep every original physical vertex and permit additional samples only on its
existing chords. Grass retains positive corner Jacobians and edge ratio <=8,
including storage. Small fixed tip angles are reported, not rejected by the
separate pavement sine criterion. No source topology or curb is reconstructed.
"""
from types import SimpleNamespace
import math
import numpy as np
from .compound_cycles import _footprint, PolylineCurb
from .compounds import CompoundGroup, RimSpan
from .four_side_charts import quality
from .mesh import publication_coordinates
from .model import Key
from .partition import BoundaryRef, CompoundPartitioner, CompoundQuadChart, JunctionCut
from .quad_charts import ChartResult
from .regions import BoundaryLoop, Region, RegionUnion, signed_area


def audit(mesh, boundary):
    """The existing grass contract, exact source lineage and a conforming disk."""
    double = quality(mesh)
    stored = quality(mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True)))
    if max(double['maximum_edge_ratio'], stored['maximum_edge_ratio']) > 8.:
        raise ValueError('Banded grass disk exceeds its physical or stored edge-ratio limit')
    rim = mesh.xyz[mesh.seams[0], :2]
    positions = []
    for point in boundary:
        hits = np.flatnonzero(np.all(rim == point, axis=1))
        if len(hits) != 1:
            raise ValueError('Banded grass disk changed an original physical rim vertex')
        positions.append(int(hits[0]))
    if positions != sorted(positions) or positions[0] != 0:
        raise ValueError('Banded grass disk changed the prescribed physical rim order')
    for i, start in enumerate(positions):
        end = positions[(i+1) % len(positions)]+(len(rim) if i+1 == len(positions) else 0)
        a, b = boundary[i], boundary[(i+1) % len(boundary)]
        d = b-a
        points = rim[np.arange(start, end+1) % len(rim)]
        t = (points-a)@d/(d@d)
        if np.max(np.linalg.norm(points-a-t[:, None]*d, axis=1)) > 1e-8 or np.any(np.diff(t) <= 0.):
            raise ValueError('Banded grass samples leave an original boundary chord')
    if abs(signed_area(rim)-signed_area(boundary)) > max(1e-8, signed_area(boundary)*1e-10):
        raise ValueError('Banded grass disk changed its physical area')
    edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    if (len(mesh.xyz)-len(unique)+len(mesh.quads) != 1 or not set(counts).issubset({1, 2})
            or np.sum(counts == 1) != len(rim)):
        raise ValueError('Banded grass atlas is not a conforming direct-quad disk')
    return dict(quality=double, stored_quality=stored, protected_vertices=len(boundary),
        added_boundary_samples=len(rim)-len(boundary), exact_physical_polygon=True,
        grass_quality_contract='positive corners and edge ratio <=8')


class _Bands:
    def __init__(self, count):
        self.count = count

    def cuts(self, group, context):
        boundary = group.loops[0].xy
        center = boundary.mean(axis=0)
        _, vectors = np.linalg.eigh((boundary-center).T@(boundary-center))
        along = (boundary-center)@vectors[:, -1]
        following = np.roll(along, -1)
        cuts, seen = [], set()
        for station in np.linspace(along.min(), along.max(), self.count+1)[1:-1]:
            selected = ((along <= station) & (station < following)) | ((following <= station) & (station < along))
            hits = []
            for edge in np.flatnonzero(selected):
                fraction = float((station-along[edge])/(following[edge]-along[edge]))
                hits.append(BoundaryRef(0, int((edge+int(fraction > .5)) % len(boundary)), 0.))
            if len(hits) != 2 or hits[0] == hits[1]:
                raise ValueError('Grass band needs two distinct original bank vertices')
            token = tuple(sorted(ref.edge for ref in hits))
            if token in seen:
                continue
            seen.add(token)
            cuts.append(JunctionCut(None, group.members[0], group.members[0], *hits))
        return tuple(cuts)


class _GrassChart(CompoundQuadChart):
    def _repairable_leaf(self, group):
        return True  # A positive leaf still requires assembled grass acceptance.


class _ExactPolyline(PolylineCurb):
    def __call__(self, parameter):
        parameter = np.asarray(parameter)
        value = super().__call__(parameter)
        knots = self.stations/self.length
        index = np.minimum(np.searchsorted(knots, parameter), len(knots)-1)
        return np.where((knots[index] == parameter)[..., None], self.points[index], value)


def build_banded_hole(patterns, boundary, spacing, max_cells, quality_iterations=100):
    boundary = np.asarray(boundary, dtype=float)
    if (boundary.ndim != 2 or boundary.shape[1] != 2 or not 4 <= len(boundary) <= 8192
            or not np.isfinite(boundary).all() or signed_area(boundary) <= 0.
            or not math.isfinite(spacing) or spacing <= 0. or type(max_cells) is not int
            or not 4 <= max_cells <= 50000 or type(quality_iterations) is not int or quality_iterations < 0):
        raise ValueError('Invalid bounded grass-disk chart inputs')
    owner = Key('grass_computational_chart', 0)
    RegionUnion().build((Region(owner, boundary),))
    source = _footprint(owner, boundary)
    if not np.array_equal(source.region.xy, boundary):
        raise ValueError('Grass band input contains unresolved physical boundary intervals')
    source.layout.boundaries = tuple(_ExactPolyline(curve.points) for curve in source.layout.boundaries)
    loop = BoundaryLoop(boundary, tuple(RimSpan(owner, label) for label in source.segments))
    group = CompoundGroup((owner,), (), (source,), (loop,))
    # These polylines are the exact chart domain; this context carries no
    # authored Street, Section, junction, source claim or physical ownership.
    context = SimpleNamespace(patterns=patterns,
        network=SimpleNamespace(settings=SimpleNamespace(curve_tolerance=1e-9)))
    center = boundary.mean(axis=0)
    _, vectors = np.linalg.eigh((boundary-center).T@(boundary-center))
    length = float(np.ptp((boundary-center)@vectors[:, -1]))
    base = max(2, int(math.ceil(length/max(signed_area(boundary)/length, 4*spacing))))
    if base > 64:
        raise ValueError('Grass disk exceeds its bounded longitudinal band budget')
    power = min(64, 2**int(math.ceil(math.log2(base))))
    failures = []
    for count in dict.fromkeys((power, base, min(64, 2*power))):
        chart = _GrassChart(patterns, spacing=spacing, max_cells=max_cells,
            quality_iterations=quality_iterations, four_sided=True,
            partitioner=CompoundPartitioner(policy=_Bands(count), require_visibility=False))
        try:
            output = chart._build(group, context, require_quality=True)
            try:
                report = audit(output.mesh, boundary)
                correction = None
            except ValueError:
                output.mesh, correction = chart._correct(output.mesh, retry_rotations=True)
                report = audit(output.mesh, boundary)
        except ValueError as failure:
            failures.append(str(failure))
            continue
        output.mesh.material = 'grass'
        return ChartResult(output.mesh, dict(report, bands=count, previous_failures=failures,
            assembled_correction=correction, construction='original-rim computational bands with direct quad leaves'))
    raise ValueError('No qualified bounded grass band atlas: '+'; '.join(failures))


def recover_failed_grass(chart, result, loop_index, failure):
    """Replace one failed grass chart, then reconcile its exact pavement seam."""
    previous = result.holes[loop_index]
    boundary = previous.mesh.xyz[previous.mesh.seams[0], :2]
    others = sum(len(hole.mesh.quads) for li, hole in result.holes.items() if li != loop_index)
    replacement = build_banded_hole(chart.patterns, boundary, chart.spacing,
        chart.max_cells-len(result.mesh.quads)-others, chart.quality_iterations)
    pavement, refs = result.mesh, result.boundary_refs
    refinement, correction = None, None
    if replacement.metadata['added_boundary_samples']:
        from .conforming_boundary_refinement import refine
        grass_cells = others+len(replacement.mesh.quads)
        pavement, refs, _, refinement = refine(pavement, refs, result.partition.group,
            chart.patterns, {loop_index: replacement.mesh.xyz[replacement.mesh.seams[0], :2]},
            max_cells=chart.max_cells-grass_cells)
        pavement, correction = chart._correct(pavement, retry_rotations=True)
    double = quality(pavement)
    stored = quality(pavement.with_positions(publication_coordinates(pavement.xyz, 1., physical=True)))
    if (max(double['maximum_edge_ratio'], stored['maximum_edge_ratio']) > 8.
            or min(double['minimum_corner_sine'], stored['minimum_corner_sine']) < .03):
        raise ValueError('Banded grass reconciliation fails unchanged pavement quality')
    expected = pavement.xyz[pavement.seams[loop_index], :2][::-1]
    rim = replacement.mesh.xyz[replacement.mesh.seams[0], :2]
    starts = np.flatnonzero(np.all(rim == expected[0], axis=1))
    if len(starts) != 1:
        raise ValueError('Banded grass reconciliation lost its shared physical seam')
    ids = np.roll(replacement.mesh.seams[0], -int(starts[0]))
    if not np.array_equal(replacement.mesh.xyz[ids, :2], expected):
        raise ValueError('Banded grass and pavement rims do not match exactly')
    replacement.mesh.seams[0] = ids
    total = len(pavement.quads)+others+len(replacement.mesh.quads)
    if total > chart.max_cells:
        raise ValueError('Banded grass reconciliation exceeds the compound cell budget')
    details = dict(previous_grass_failure=str(failure), grass=replacement.metadata,
        pavement_refinement=refinement, pavement_correction=correction,
        pavement_quality=double, stored_pavement_quality=stored)
    result.mesh, result.boundary_refs = pavement, refs
    result.holes[loop_index] = replacement
    result.metadata = dict(result.metadata, total_quads=total, quality=double,
        banded_grass_recovery={**result.metadata.get('banded_grass_recovery', {}), loop_index: details})
    return pavement


def retry_failed_family(parent, group, context, failure, previous, source_failure):
    """Only after all old physical families fail, retry one typed grass family."""
    from .physical_hole_cuts import PhysicalHoleCuts
    from .source_strip_containment import certify
    sources = {claim.section_id for claim in group.claims}
    original_proof = certify(context, group, sources) if sources else None
    use_sources, bridges, vertices = failure.family
    chart = CompoundQuadChart(parent.patterns, spacing=parent.spacing, max_cells=parent.max_cells,
        quality_iterations=min(100, parent.quality_iterations), four_sided=True, recover_grass=True,
        partitioner=CompoundPartitioner(policy=PhysicalHoleCuts(spacing=parent.spacing,
            source_corridors=use_sources, hole_bridges=bridges, existing_vertices=vertices), require_visibility=False))
    result = chart.build(group, context)
    final_proof = certify(context, result.partition.group, sources) if sources else None
    checks = []
    for mesh in (result.mesh, *(hole.mesh for hole in result.holes.values())):
        double = quality(mesh)
        stored = quality(mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True)))
        if (max(double['maximum_edge_ratio'], stored['maximum_edge_ratio']) > 8.
                or (mesh.material != 'grass' and min(double['minimum_corner_sine'], stored['minimum_corner_sine']) < .03)):
            raise ValueError('Recovered grass family fails unchanged physical or stored quality')
        checks.append(dict(material=mesh.material, physical=double, stored=stored))
    result.metadata = dict(result.metadata, physical_hole_partition=True,
        source_corridors=use_sources, source_hole_bridges=bridges, existing_cut_vertices=vertices,
        previous_physical_atlas_failures=previous,
        grass_family_recovery=dict(previous_source_normal_failure=str(source_failure),
            original_width_proof=original_proof, final_width_proof=final_proof,
            stored_checks=checks, physical_boundary_changed=False, cell_budget=parent.max_cells))
    if final_proof is not None:
        result.metadata['source_width_visibility'] = final_proof
    return result
