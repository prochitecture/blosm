"""Opt-in fixed-chart source-body seeds and existing affine height guides.

This module is not registered with a provider or builder. It changes heights
only, returns constraints for the existing slope/drainage/publication loop,
and never promotes body guide nodes to physical drainage outlets.
"""
from dataclasses import dataclass
import math
import numpy as np
from .updates import drain_steps
from .model import Key
from .mesh import cross
from .regions import _classify
from .height_slabs import HeightSlabs
from .source_body_references import (BodyGuideLimits, BodyWork, OriginalBodyReference,
    PairedBodyReference, SourceBodyUnavailable, NoPairedBodyDomain, declared_frame)

SOURCE_BODY_HEIGHT_POLICY = '1-opt-in-original-owned-body-height-guides'


def _smooth(value):
    t = np.clip(value, 0., 1.)
    return t**3*(10.+t*(-15.+6.*t))


def _weight(station, source, cores):
    """C2 soft transitions are confined to measured source support intervals."""
    result = np.ones(len(station))
    first, last = cores[0][0], cores[-1][1]
    if first > source.first:
        result *= _smooth((station-source.first)/(first-source.first))
    if last < source.last:
        result *= _smooth((source.last-station)/(source.last-last))
    for left, right in zip(cores, cores[1:]):
        a, b = left[1], right[0]
        selected = (station > a) & (station < b)
        # The midpoint is a computational seed knot inside a proven influence
        # interval. It is never a source End, mouth or physical ownership cut.
        result[selected] *= _smooth(2.*abs(station[selected]-.5*(a+b))/(b-a))
    return result


def _core_mask(station, cores):
    result = np.zeros(len(station), bool)
    for first, last in cores:
        result |= (station >= first) & (station <= last)
    return result


def _represented(context, source):
    owner, layout = source.region.owner, source.layout
    if owner.kind == 'cycle_strip':
        if owner.id not in context.network.section_index:
            raise ValueError('Source-body strip lost its original Section')
        return frozenset((owner.id,))
    if owner.kind == 'cycle_port' and len(getattr(layout, 'approaches', ())) == 1:
        a = layout.approaches[0]
        if hasattr(a.path, 'trims') and hasattr(a.path, 'pieces'):
            if any(piece.path is not context[Key('path', piece.section.id)] for piece in a.path.pieces):
                raise ValueError('Source-body port changed an original Route piece')
            return frozenset(sid for sid, end in a.path.trims(layout.outward_stations[0]))
    return frozenset()


def _classify_owned(points, polygon, work):
    work.charge('classification_work', len(points))
    result = np.zeros(len(points), bool)
    selected = np.flatnonzero(np.all(points >= polygon.min(0), axis=1)
                              & np.all(points <= polygon.max(0), axis=1))
    work.charge('classification_work', len(selected)*len(polygon))
    if len(selected):
        inside, boundary = _classify(points[selected], polygon, 1e-7)
        result[selected] = inside | (boundary >= 0)
    return result


def _owned_references(context, group, work):
    proof = getattr(group.physical_layout, 'source_cycle_full_width_proof', None)
    if proof is None or not isinstance(proof.get('sections'), (list, tuple)):
        raise SourceBodyUnavailable('Source-body guides need an explicit original owned-range proof')
    by_id = {}
    for row in proof['sections']:
        sid = row['section']
        first, last = map(float, row['interval'])
        by_id.setdefault(sid, []).append((first, last))
        if sum(map(len, by_id.values())) > work.limits.maximum_cells:
            raise SourceBodyUnavailable('Source-body original interval budget exhausted')
    if not by_id:
        raise SourceBodyUnavailable('Source-body guide has no certified original owned range')
    originals = []
    for sid, ranges in sorted(by_id.items()):
        ordered = sorted(ranges)
        if any(b[0] < a[1] for a, b in zip(ordered, ordered[1:])):
            raise ValueError('Source-body ownership must retain distinct original disjoint ranges')
        originals.extend(OriginalBodyReference.build(context, sid, a, b, work) for a, b in ordered)
    return tuple(originals)


def _contact_bounds(context, sources, originals):
    from .generated_junction_curbs import is_generated_junction_area
    result = {source: [source.first, source.last] for source in originals}
    seen = set()
    for area in sources:
        if is_generated_junction_area(area):continue
        for corner in getattr(area.layout, 'corners', ()):
            if id(corner) in seen:
                continue
            seen.add(id(corner))
            for approach, station in zip((corner.first, corner.second), corner.stations):
                if station <= 0.:
                    continue
                if hasattr(approach.path, 'trims'):
                    trims = approach.path.trims(station)
                elif approach.path is context[Key('path', approach.source.id)]:
                    trims = {(approach.source.id, approach.end): float(station if approach.end == 0
                                                               else approach.path.length-station)}
                else:
                    # Prepared multi-source area remains a foreign influence
                    # operand. Do not invent an original station from its XY.
                    continue
                for (sid, end), at in trims.items():
                    for source in originals:
                        if source.section_id != sid:
                            continue
                        if not math.isfinite(at) or not 0. <= at <= source.path.length:
                            raise ValueError('Source-body corner escaped its actual original source')
                        result[source][end] = (max(result[source][end], at) if end == 0
                                               else min(result[source][end], at))
    return result


def _influence(context, group, originals, representations, work):
    return drain_steps(_influence_steps(context, group, originals, representations, work))


def _influence_steps(context, group, originals, representations, work):
    from .source_mouth_clearance import _positive_overlap
    from .finite_containment import separated_segment
    contacts = _contact_bounds(context, group.area_footprints, originals)
    core, reports = {}, []
    for source in originals:
        same = {source.section_id}
        for bundle in context.network.bundles:
            roles = [role for role, member in enumerate(bundle.members)
                     if any(v.section_id == source.section_id for v in member.path)]
            if roles:
                if len(bundle.members) != 2 or len(roles) != 1:
                    raise SourceBodyUnavailable('Source-body influence needs one unique declared lateral role')
                # A same-role/nonlocal strip remains a foreign influence. Only
                # an actual opposite original body is paired, never a nearby
                # road or every path sharing this Bundle key.
                same |= {v.section_id for v in bundle.members[1-roles[0]].path}
        foreign = [s for s in group.area_footprints
                   if not representations[s.region.owner] or not representations[s.region.owner] <= same]
        first = np.array([s.region.xy.min(0) for s in foreign]).reshape(-1, 2)
        last = np.array([s.region.xy.max(0) for s in foreign]).reshape(-1, 2)
        at, path = source.stations, source.path
        work.charge('cells', len(at)-1)
        right = path.point(at, path.section.offset+.5*path.section.width)
        left = path.point(at, path.section.offset-.5*path.section.width)
        bad = np.zeros(len(at)-1, bool)
        witnesses = []
        clearance = 2.*work.limits.source_tolerance+1e-7
        for i in range(len(bad)):
            if i % 8 == 0:
                yield
            cell = np.array((right[i], right[i+1], left[i+1], left[i]))
            edge = np.roll(cell, -1, axis=0)-cell
            turns = cross(edge, np.roll(edge, -1, axis=0))
            if not (np.all(turns > 1e-12) or np.all(turns < -1e-12)):
                raise SourceBodyUnavailable('Source-body influence needs unfolded original finite source cells')
            nearby = np.flatnonzero(np.all(cell.max(0)+clearance >= first, axis=1)
                                     & np.all(cell.min(0)-clearance <= last, axis=1))
            for j in nearby:
                work.charge('pairs', 1)
                polygon = foreign[j].region.xy
                if (_positive_overlap(cell, polygon, tolerance=1e-7)
                        or not all(separated_segment(a, b, polygon, clearance=clearance)
                                   for a, b in zip(cell, np.roll(cell, -1, axis=0)))):
                    bad[i] = True
                    witnesses.append(dict(interval=(float(at[i]), float(at[i+1])),
                                          owner=foreign[j].region.owner.text()))
        ranges, begin = [], None
        for i in range(len(bad)+1):
            free = i < len(bad) and not bad[i]
            if free and begin is None:
                begin = i
            if not free and begin is not None:
                a = max(float(at[begin]), contacts[source][0])
                b = min(float(at[i]), contacts[source][1])
                if a < b:
                    ranges.append((a, b))
                begin = None
        core[source] = tuple(ranges)
        reports.append(dict(section=source.section_id, owned=(source.first, source.last),
            corner_contacts=contacts[source], core_ranges=ranges, cells=len(bad),
            source_chord_error=work.limits.source_tolerance, full_cell_clearance=clearance,
            influence_witnesses=witnesses, paired_body_sources_excluded=sorted(same-{source.section_id})))
    return core, reports


@dataclass(frozen=True)
class SourceBodyGuides:
    seed: object
    height_slabs: object
    core_nodes: object
    core_targets: object
    physical_fixed: object
    references: tuple
    paired_references: tuple
    dependencies: tuple
    ledger: object

    def physical_field(self, final_mesh, *, datum, ceiling):
        """Reconstruct only after the enclosing existing solver checks pass."""
        from .surface_query import QuadHeightField
        self.assert_unchanged(final_mesh)
        if self.height_slabs is not None and self.height_slabs.violation(final_mesh.xyz[:, 2]) > 1e-8:
            raise ValueError('Source-body final field lost an original hard body guide')
        return QuadHeightField(final_mesh, datum=datum, ceiling=ceiling,
            metadata=dict(source_body_height_policy=SOURCE_BODY_HEIGHT_POLICY),
            max_cells=self.ledger['maximum_mesh_cells'])

    def assert_unchanged(self, mesh):
        if (not np.array_equal(mesh.xyz[:, :2], self.seed.xyz[:, :2])
                or not np.array_equal(mesh.quads, self.seed.quads)
                or not np.array_equal(mesh.xyz[self.physical_fixed], self.seed.xyz[self.physical_fixed])
                or mesh.material != self.seed.material or set(mesh.seams) != set(self.seed.seams)
                or any(not np.array_equal(mesh.seams[k], self.seed.seams[k]) for k in self.seed.seams)
                or set(mesh.uv_layers) != set(self.seed.uv_layers)
                or any(not np.array_equal(mesh.uv_layers[k], self.seed.uv_layers[k]) for k in self.seed.uv_layers)):
            raise ValueError('Source-body guide changed physical boundary, XY, quads, seams, material or UVs')


def prepare(context, group, physical_mesh, boundary_refs, *, limits=None, require_complete=True):
    return drain_steps(prepare_steps(context, group, physical_mesh, boundary_refs,
        limits=limits, require_complete=require_complete))


def prepare_steps(context, group, physical_mesh, boundary_refs, *, limits=None, require_complete=True):
    """Return height-only seeds/rows; never run or register a global solver.

    A typed unsupported paired scope is fatal by default. Review callers may
    set require_complete=False to inspect supported rows plus explicit failures;
    that report never qualifies the complete source-body policy.
    """
    from .source_strip_containment import certify_intervals_steps
    from .bundle_axis import AxialChartUnavailable
    limits = BodyGuideLimits() if limits is None else limits
    if not isinstance(limits, BodyGuideLimits) or type(require_complete) is not bool:
        raise ValueError('Source-body preparation needs explicit bounded limits and Boolean admission')
    mesh = physical_mesh
    if len(mesh.quads) > limits.maximum_mesh_cells or len(mesh.xyz) > 4*limits.maximum_mesh_cells:
        raise SourceBodyUnavailable('Source-body guide exceeds its physical chart allocation budget')
    mesh.validate()
    fixed = np.zeros(len(mesh.xyz), bool)
    ids = np.asarray(tuple(boundary_refs))
    if (ids.ndim != 1 or ids.dtype.kind not in 'iu' or not len(ids)
            or np.any(ids < 0) or np.any(ids >= len(fixed))):
        raise ValueError('Source-body outlets must be the actual physical boundary references')
    edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    physical_boundary = np.unique(unique[counts == 1])
    if np.any(counts > 2) or not np.array_equal(np.unique(ids), physical_boundary):
        raise ValueError('Source-body outlets must cover exactly the physical mesh boundary')
    fixed[ids] = True
    work = BodyWork(limits)
    originals = _owned_references(context, group, work)
    ranges = {}
    for source in originals:
        ranges.setdefault(source.section_id, []).append((source.first, source.last))
    before = yield from certify_intervals_steps(context, group, ranges, maximum_error=limits.maximum_source_error,
                                maximum_cells=limits.maximum_cells)
    sources = group.area_footprints
    if sum(len(s.region.xy) for s in sources) > limits.maximum_polygon_vertices:
        raise SourceBodyUnavailable('Source-body area operand vertex budget exhausted')
    representations = {s.region.owner: _represented(context, s) for s in sources}
    if len(representations) != len(sources):
        raise ValueError('Source-body area operands require unique explicit owners')
    core, influence = yield from _influence_steps(context, group, originals, representations, work)
    masks = {}
    def available(source_ids):
        token = frozenset(source_ids)
        if token not in masks:
            selected = np.zeros(len(mesh.xyz), bool)
            for source in sources:
                represented = representations[source.region.owner]
                if represented and represented <= token:
                    selected |= _classify_owned(mesh.xyz[:, :2], source.region.xy, work)
            selected &= ~fixed
            for source in sources:
                represented = representations[source.region.owner]
                if represented and represented <= token:
                    continue
                at = np.flatnonzero(selected)
                if not len(at):
                    break
                selected[at[_classify_owned(mesh.xyz[at, :2], source.region.xy, work)]] = False
            masks[token] = selected
        return masks[token].copy()
    xyz, targets, seeds, reports = mesh.xyz.copy(), {}, {}, []
    def add(nodes, heights, weights, hard, scope, *, replace_soft=False):
        for node, height, weight, required in zip(nodes, heights, weights, hard):
            node = int(node)
            if replace_soft and required and node not in targets:
                seeds.pop(node, None)
            if node in seeds and weight > 1e-12 and abs(seeds[node][0]-height) > 1e-8:
                raise SourceBodyUnavailable('Overlapping source-body soft domains need an explicit junction blend', scope=scope)
            if weight > 1e-12:
                seeds[node] = (float(height), scope)
            xyz[node, 2] = (1.-weight)*mesh.xyz[node, 2]+weight*height
            if required:
                if node in targets and abs(targets[node]-height) > 1e-8:
                    raise ValueError('Original hard body height guides disagree')
                targets[node] = float(height)
    def original_nodes(source, *, remaining=False):
        selected = available((source.section_id,))
        if remaining:
            selected[np.asarray(tuple(targets), dtype=np.int64)] = False
        nodes = np.flatnonzero(selected)
        if not len(nodes):
            return None
        work.charge('classification_work', len(nodes)*(len(source.stations)-1))
        station, q, valid, proof = source.inverse(mesh.xyz[nodes, :2], max_cells=limits.maximum_cells)
        valid &= np.isfinite(q) & (abs(q) < 1.-1e-8)
        nodes, station, q = nodes[valid], station[valid], q[valid]
        hard = _core_mask(station, core[source])
        if remaining:
            keep = np.array([required or int(node) not in seeds for node, required in zip(nodes, hard)])
            nodes, station, q, hard = nodes[keep], station[keep], q[keep], hard[keep]
        return nodes, source.profile.datum+source.profile.value(q), _weight(station, source, core[source]), hard, proof
    for source in originals:
        if source.profile.half or not core[source]:
            continue
        yield
        values = original_nodes(source)
        if values is not None:
            nodes, heights, weights, hard, proof = values
            add(nodes, heights, weights, hard, dict(kind='original_full', section=source.section_id))
            reports.append(dict(kind='original_full', section=source.section_id, core_nodes=int(hard.sum()), inverse=proof))
    pairs, unsupported, dependencies = [], [], set(key for s in originals for key in s.dependencies)
    for bundle in context.network.bundles:
        owned = [[s for v in member.path for s in originals
                  if s.section_id == v.section_id and core[s]] for member in bundle.members]
        if len(owned) != 2 or not all(owned):
            continue
        try:
            frame, frame_proof = declared_frame(context, bundle)
        except AxialChartUnavailable as error:
            failure = SourceBodyUnavailable(str(error), scope=dict(bundle=bundle.id, kind='paired_coordinates'))
            if require_complete:
                raise failure from error
            unsupported.append(failure.scope)
            continue
        for a in owned[0]:
            for b in owned[1]:
                yield
                work.charge('field_pairs', 1)
                try:
                    pair = PairedBodyReference.build(context, bundle, (a, b), frame, frame_proof)
                    shared = pair.whole_body_intervals(core)
                except NoPairedBodyDomain:
                    reports.append(dict(kind='disjoint_finite_paired_ranges', bundle=bundle.id,
                                        sections=(a.section_id, b.section_id)))
                    continue
                except (SourceBodyUnavailable, AxialChartUnavailable) as error:
                    scope = dict(bundle=bundle.id, sections=(a.section_id, b.section_id), reason=str(error))
                    if require_complete:
                        raise SourceBodyUnavailable(str(error), scope=scope) from error
                    unsupported.append(scope)
                    continue
                dependencies.update(pair.dependencies)
                if not shared:
                    reports.append(dict(kind='no_complete_paired_core', bundle=bundle.id,
                                        sections=(a.section_id, b.section_id)))
                    continue
                pairs.append(pair)
                local = frame.local(mesh.xyz[:, :2])
                selected = available((a.section_id, b.section_id))
                selected &= (local[:, 0] > pair.first) & (local[:, 0] < pair.last)
                nodes = np.flatnonzero(selected)
                if not len(nodes):
                    continue
                values = local[nodes]
                edges = pair.field.edges(values[:, 0])
                valid = ((values[:, 1] > np.minimum(edges[1], edges[3])+1e-7)
                         & (values[:, 1] < np.maximum(edges[0], edges[2])-1e-7))
                nodes, values = nodes[valid], values[valid]
                weights = _weight(values[:, 0], pair, shared)
                hard = _core_mask(values[:, 0], shared)
                scope = dict(kind='declared_paired', bundle=bundle.id, sections=(a.section_id, b.section_id))
                add(nodes, pair.values(mesh.xyz[nodes, :2]), weights, hard, scope)
                fibre_report = (dict(whole_body_fibres=shared[0]) if len(shared) == 1
                                else dict(whole_body_fibre_ranges=shared))
                reports.append(dict(scope, core_nodes=int(hard.sum()), **fibre_report,
                                    finite_field=(pair.first, pair.last), frame=frame_proof))
    for source in originals:
        if not source.profile.half or not core[source]:
            continue
        yield
        values = original_nodes(source, remaining=True)
        if values is not None:
            nodes, heights, weights, hard, proof = values
            add(nodes, heights, weights, hard, dict(kind='unique_original_half', section=source.section_id),
                replace_soft=True)
            reports.append(dict(kind='unique_original_half', section=source.section_id,
                                core_nodes=int(hard.sum()), inverse=proof))
    nodes = np.array(sorted(targets), dtype=np.int32)
    heights = np.array([targets[node] for node in nodes])
    constraints = None
    if len(nodes):
        boundary = np.flatnonzero(fixed)
        low = int(boundary[np.argmin(mesh.xyz[boundary, 2])])
        high = int(boundary[np.argmax(mesh.xyz[boundary, 2])])
        extras = [int(i) for i in boundary if i not in (low, high)]
        if low == high or not extras:
            raise SourceBodyUnavailable('Exact body rows need distinct actual bounding physical profile anchors')
        alpha = (heights-mesh.xyz[low, 2])/(mesh.xyz[high, 2]-mesh.xyz[low, 2])
        if np.any(alpha < 0.) or np.any(alpha > 1.):
            raise SourceBodyUnavailable('Original body target exceeds the actual fixed profile height envelope')
        if len(nodes) > 8*len(mesh.quads):
            raise SourceBodyUnavailable('Original body rows exceed the existing HeightSlabs row budget')
        vertices = np.column_stack((nodes, np.full(len(nodes), low), np.full(len(nodes), high),
                                    np.full(len(nodes), extras[0]))).astype(np.int32)
        coefficients = np.column_stack((np.ones(len(nodes)), -(1.-alpha), -alpha, np.zeros(len(nodes))))
        constraints = HeightSlabs(vertices, coefficients, np.zeros(len(nodes)))
    seed = mesh.with_positions(xyz)
    after = yield from certify_intervals_steps(context, group, ranges, maximum_error=limits.maximum_source_error,
                               maximum_cells=limits.maximum_cells)
    if before != after:
        raise ValueError('Source-body preparation changed its complete original width proof')
    ledger = dict(policy=SOURCE_BODY_HEIGHT_POLICY, opt_in=True, runtime_registered=False,
        complete_supported_scope=not unsupported, unsupported=unsupported, source_ranges=ranges,
        source_influence=influence, fields=reports, work=dict(work.used),
        source_width_proof_before=before, source_width_proof_after=after,
        physical_only_outlets=True, physical_boundary_nodes=int(fixed.sum()), core_nodes=len(nodes),
        no_hard_guides_in_junction_support=True, exact_xy_quads_uv_material=True,
        maximum_mesh_cells=limits.maximum_mesh_cells, qualification='seed_and_constraints_only')
    for value in (nodes, heights, fixed):
        value.flags.writeable = False
    result = SourceBodyGuides(seed, constraints, nodes, heights, fixed, originals, tuple(pairs),
                              tuple(sorted(dependencies)), ledger)
    result.assert_unchanged(mesh)
    return result
