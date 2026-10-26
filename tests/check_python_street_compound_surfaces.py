"""Qualify compound surface primitives separately from default import dispatch."""
import argparse
from dataclasses import replace
from html import escape
import json
from pathlib import Path
import time
import numpy as np
import _python_street_package
from _python_street_fixtures import compound_cross_network, short_junction_pair
from check_python_street_footprints import ApproachInputs
from tested_python_streets import VERSION, model
from tested_python_streets.compounds import CompoundGroup, CompoundPlanner, rim_label, sample_footprint
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.mesh import cross
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.regions import RegionUnion
from tested_python_streets.updates import BuildContext, UpdateCoordinator
from tested_python_streets.validation import validate_overlaps

CASE_TITLES = {'compound_divided_cross_preparation': 'Divided crossing — four source junctions',
               'short_pair_preparation': 'Two closely spaced T-junctions',
               'unequal_crown_compound': 'T-junction pair — unequal incoming crowns',
               'nonzero_common_datum': 'T-junction pair — common datum +12.5 m'}


def curve_error(group, layouts):
    sources = {p.region.owner: p for p in group.footprints}
    maximum = 0.
    for loop in group.loops:
        for i, span in enumerate(loop.spans):
            label = rim_label(span, sources)
            if label.kind != 'curb':
                continue
            curve = layouts[span.owner].boundaries[label.index]
            actual = curve(np.linspace(label.first, label.last, 101))
            a, b = loop.xy[i], loop.xy[(i + 1) % len(loop.xy)]
            edge = b - a
            fraction = np.clip((actual - a) @ edge / (edge @ edge), 0., 1.)
            maximum = max(maximum, float(np.linalg.norm(actual - a - fraction[:, None] * edge, axis=1).max()))
    return maximum


def entrance_slopes(mesh, conditions, exaggeration):
    owners = {}
    for fi, q in enumerate(mesh.quads):
        for side in range(4):
            key = tuple(sorted((int(q[side]), int(q[(side + 1) % 4]))))
            owners.setdefault(key, []).append((fi, side))
    records = []
    for boundary in conditions:
        for port in boundary.ports:
            maximum = 0.
            for first, last in zip(port.vertices, port.vertices[1:]):
                xy = boundary.xy[[first, last]]
                ids = [int(np.argmin(np.linalg.norm(mesh.xyz[:, :2] - point, axis=1))) for point in xy]
                assert np.max(np.linalg.norm(mesh.xyz[ids, :2] - xy, axis=1)) < 1e-8
                hits = owners[tuple(sorted(ids))]
                assert len(hits) == 1
                fi, side = hits[0]
                p = mesh.xyz[mesh.quads[fi]].copy()
                p[:, 2] /= exaggeration
                u, v = ((.5, 0.), (1., .5), (.5, 1.), (0., .5))[side]
                a = (1. - v) * (p[1] - p[0]) + v * (p[2] - p[3])
                b = (1. - u) * (p[3] - p[0]) + u * (p[2] - p[1])
                determinant = cross(a[:2], b[:2])
                gradient = np.array(((a[2] * b[1] - b[2] * a[1]) / determinant,
                                     (a[0] * b[2] - b[0] * a[2]) / determinant))
                tangent = (xy[1] - xy[0]) / np.linalg.norm(xy[1] - xy[0])
                normal = np.array((tangent[1], -tangent[0]))
                maximum = max(maximum, abs(float(gradient @ normal)))
            records.append(dict(section=port.section_id, end=port.end, maximum=maximum))
    return records


def draw(cases, destination):
    height = 120 + 340 * len(cases)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="1160" height="{height}" viewBox="0 0 1160 {height}">',
           '<rect width="100%" height="100%" fill="#f3f6fa"/>',
           '<g font-family="Segoe UI, sans-serif" fill="#17324d">',
           '<text x="28" y="36" font-size="25" font-weight="600">Compound surfaces — numerical construction</text>',
           '<text x="28" y="65" font-size="16">Full entrance profiles · zero low curbs · highest-crown ceiling · direct quads</text>',
           '<text x="28" y="90" font-size="16">Left: height map. Right: exaggerated relief. Construction checks; normal import coverage is reported separately.</text>']
    for index, (name, surface, record) in enumerate(cases):
        y = 115 + 340 * index
        svg.append(f'<rect x="18" y="{y}" width="1124" height="325" rx="12" fill="white" stroke="#cbd5e1"/>')
        svg.append(f'<text x="34" y="{y + 26}" font-size="17" font-weight="600">{escape(name)}</text>')
        meshes = (surface.mesh, *surface.holes.values())
        for view in (0, 1):
            projected = []
            for mesh in meshes:
                xy, z = mesh.xyz[:, :2], mesh.xyz[:, 2]
                p = xy if view == 0 else np.column_stack((.819152 * xy[:, 0] - .573576 * xy[:, 1],
                                                         .286788 * xy[:, 0] + .409576 * xy[:, 1] + .866025 * z))
                projected.append(p)
            all_points = np.row_stack(projected)
            lo, hi = all_points.min(axis=0), all_points.max(axis=0)
            center = .5 * (lo + hi)
            scale = min(495. / (hi[0] - lo[0]), 225. / (hi[1] - lo[1]))
            for mesh, p in zip(meshes, projected):
                order = np.argsort(p[mesh.quads].mean(axis=1)[:, 1]) if view else range(len(mesh.quads))
                for fi in order:
                    face = mesh.quads[fi]
                    z = float(mesh.xyz[face, 2].mean())
                    fraction = np.clip((z - surface.metadata['datum']) / max(.001,
                        (surface.metadata['ceiling'] - surface.metadata['datum']) * record['exaggeration']), 0., 1.)
                    low, high = np.array((61, 98, 134)), np.array((197, 224, 239))
                    rgb = (low + fraction * (high - low)).astype(int)
                    color = '#b9d4a7' if mesh.material == 'grass' else '#' + ''.join(f'{v:02x}' for v in rgb)
                    path = 'M' + ' L'.join(f'{294 + view * 560 + scale * (q[0] - center[0]):.3f},{y + 163 - scale * (q[1] - center[1]):.3f}' for q in p[face]) + ' Z'
                    svg.append(f'<path d="{path}" fill="{color}" stroke="#5f7e94" stroke-width=".35"/>')
        svg.append(f'<text x="34" y="{y + 300}" font-size="14">{record["total_quads"]} quads; max slope {100. * record["physical_slopes"]["maximum"]:.2f}%; incoming boundary {100. * record["field"]["maximum_boundary_slope"]:.2f}%; height ceiling {record["ceiling"]:.3f} m.</text>')
    svg.append('</g></svg>')
    destination.write_text('\n'.join(svg), encoding='utf8')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sites', type=int, nargs='*', default=[],
                        help='Additional standalone junction IDs from the optional numerical input')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    mixed = short_junction_pair(model)
    mixed = replace(mixed, uuid='unequal_crown_compound', sections=tuple(replace(s, rise=.04 if s.id % 2 else .08) for s in mixed.sections))
    raised = short_junction_pair(model)
    raised = replace(raised, uuid='nonzero_common_datum', sections=tuple(replace(s, datum=12.5) for s in raised.sections))
    networks = [compound_cross_network(model), short_junction_pair(model), mixed, raised]
    if args.input:
        networks.extend(model.Network.from_dict(source) for source in json.loads(args.input.read_text(encoding='utf8')))
    records, pictures = [], []
    for network in networks:
        patterns = PatternLibrary()
        prepared = UpdateCoordinator((PathProvider(), ProfileProvider(), ApproachInputs()), patterns).prepare(network)
        context = BuildContext(network, prepared.products, patterns)
        layouts = {}
        for junction in network.intersections:
            try:
                layouts[junction.key] = IntersectionProvider.footprint(context, junction)
            except ValueError:
                pass  # Complete footprint support is audited independently.
        groups = list(CompoundPlanner().build(layouts))
        covered = {key for group in groups for key in group.members}
        for junction in network.intersections:
            if junction.id in args.sites and junction.key in layouts and junction.key not in covered:
                footprint = sample_footprint(junction.key, layouts[junction.key])
                groups.append(CompoundGroup((junction.key,), (), (footprint,), RegionUnion().build((footprint.region,))))
        for group in groups:
            record = dict(network=network.uuid, members=[key.text() for key in group.members])
            start = time.perf_counter()
            try:
                surface = CompoundSurfaceBuilder(patterns).build(group, context)
                sampled = surface.chart.partition.group
                conditions = [BoundaryConditions().build(sampled, loop, layouts) for loop in sampled.loops]
                ex = network.settings.exaggeration
                boundary_error, position_error = 0., 0.
                for vertex, ref in surface.boundary_refs.items():
                    xy, height, gradient = conditions[ref.loop].edge_samples(ref.edge, ref.fraction)
                    actual_height = surface.field.datum + (surface.mesh.xyz[vertex, 2] - surface.field.datum) / ex
                    boundary_error = max(boundary_error, float(abs(actual_height - height)))
                    position_error = max(position_error, float(np.linalg.norm(surface.mesh.xyz[vertex, :2] - xy)))
                assert boundary_error < 1e-8
                assert position_error < 1e-8
                approximation = curve_error(sampled, layouts)
                assert approximation <= network.settings.curve_tolerance + 1e-7
                meshes = {model.Key('compound', 0): surface.mesh}
                meshes.update({model.Key('hole', li): mesh for li, mesh in surface.holes.items()})
                overlap = validate_overlaps(meshes)
                joins = entrance_slopes(surface.mesh, conditions, ex)
                # This check covers entrance midpoint longitudinal slopes, not
                # every possible trajectory through the intersection interior.
                assert max(p['maximum'] for p in joins) < .03
                physical = (surface.mesh.xyz[:, 2] - surface.field.datum) / ex
                edges = np.stack((surface.mesh.quads, np.roll(surface.mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2)
                a, b = np.r_[edges[:, 0], edges[:, 1]], np.r_[edges[:, 1], edges[:, 0]]
                neighbor_minimum = np.full(len(physical), np.inf)
                np.minimum.at(neighbor_minimum, a, physical[b])
                interior = np.ones(len(physical), bool)
                interior[list(surface.boundary_refs)] = False
                minima = np.flatnonzero(interior & (physical < neighbor_minimum - 1e-6))
                assert not len(minima), 'Isolated interior vertex minima require a different height construction'
                records.append(dict(**record, status='SURFACE_PRIMITIVE_QUALIFIED',
                    **surface.metadata, boundary_height_error=boundary_error,
                    boundary_position_error=position_error, maximum_curve_error=approximation,
                    entrance_longitudinal_slopes=joins, overlap_check=overlap,
                    isolated_interior_vertex_minima=len(minima), minimum_depth_tolerance=1e-6,
                    exaggeration=ex, seconds=time.perf_counter() - start))
                title = CASE_TITLES.get(network.uuid, 'Real OSM surface — sites ' + ', '.join(str(key.id) for key in group.members))
                pictures.append((title, surface, records[-1]))
            except (ValueError, AssertionError) as exc:
                records.append(dict(**record, status='UNQUALIFIED', reason=str(exc), seconds=time.perf_counter() - start))
            print(json.dumps(records[-1]), flush=True)
    report = dict(implementation=VERSION, default_provider_integration=False, complete_network_coverage=False, cases=records)
    (args.output / 'compound_surface_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    draw(pictures, args.output / 'compound_surfaces.svg')


if __name__ == '__main__':
    main()
