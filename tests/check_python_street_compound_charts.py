"""Conforming planar compound pavement/hole charts; no surface coverage claim."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import _python_street_package
from _python_street_fixtures import compound_cross_network, short_junction_pair
from check_python_street_footprints import ApproachInputs
from tested_python_streets import VERSION, model
from tested_python_streets.compounds import CompoundPlanner
from tested_python_streets.compounds import rim_label
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.updates import BuildContext, UpdateCoordinator
from tested_python_streets.validation import validate_overlaps


def draw(cases, destination):
    from html import escape
    width, height = 1120, 115 + 405 * ((len(cases) + 1) // 2)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f3f6fa"/>', '<g font-family="Segoe UI, sans-serif" fill="#17324d">',
        '<text x="25" y="36" font-size="25" font-weight="600">Compound quad charts — planar preparation</text>',
        '<text x="25" y="65" font-size="16">Source-linked seams join the original boundaries; the separator hole remains.</text>',
        '<text x="25" y="90" font-size="16">Shared seam vertices match. Profiles and the height field remain a separate construction stage.</text>']
    for index, (name, result, report) in enumerate(cases):
        x, y = 15 + 550 * (index % 2), 110 + 405 * (index // 2)
        lo, hi = result.mesh.xyz[:, :2].min(axis=0), result.mesh.xyz[:, :2].max(axis=0)
        center = .5 * (lo + hi)
        scale = min(475. / (hi[0] - lo[0]), 285. / (hi[1] - lo[1]))
        def path(points):
            return 'M' + ' L'.join(f'{x + 265 + scale * (p[0] - center[0]):.3f},{y + 205 - scale * (p[1] - center[1]):.3f}' for p in points) + ' Z'
        svg.extend([f'<rect x="{x}" y="{y}" width="530" height="390" rx="12" fill="white" stroke="#cbd5e1"/>',
            f'<text x="{x + 15}" y="{y + 27}" font-size="17" font-weight="600">{escape(name)}</text>'])
        offset = 0
        for bi, block in enumerate(result.metadata['blocks']):
            count = block['quality']['quads']
            color = ('#edf4fc', '#e7f1f8', '#f1f5fb', '#eaf3f0')[bi % 4]
            for face in result.mesh.quads[offset:offset + count]:
                svg.append(f'<path d="{path(result.mesh.xyz[face, :2])}" fill="{color}" stroke="#7392a9" stroke-width=".45"/>')
            offset += count
        for hole in result.holes.values():
            for face in hole.mesh.quads:
                svg.append(f'<path d="{path(hole.mesh.xyz[face, :2])}" fill="#c3dcba" stroke="#7c9e71" stroke-width=".45"/>')
        for loop in result.mesh.seams.values():
            svg.append(f'<path d="{path(result.mesh.xyz[loop, :2])}" fill="none" stroke="#234c70" stroke-width="1.5"/>')
        for edge_index, count in result.samples.cut_intervals.items():
            edge = result.partition.edges[edge_index]
            svg.append(f'<path d="{path(result.partition.xy[[edge.first, edge.last]])}" fill="none" stroke="#2d6372" stroke-width="1.4"/>')
        svg.append(f'<text x="{x + 15}" y="{y + 370}" font-size="14">{report["total_quads"]} quads; {report["blocks"]} blocks; max edge ratio {report["maximum_edge_ratio"]:.2f}</text>')
    svg.append('</g></svg>')
    destination.write_text('\n'.join(svg), encoding='utf8')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    networks = [compound_cross_network(model), short_junction_pair(model)]
    if args.input:
        networks.extend(model.Network.from_dict(data) for data in json.loads(args.input.read_text(encoding='utf8')))
    records, pictures = [], []
    for network in networks:
        patterns = PatternLibrary(network.settings.pattern_cache_bytes)
        prepared = UpdateCoordinator((PathProvider(), ProfileProvider(), ApproachInputs()), patterns).prepare(network)
        context = BuildContext(network, prepared.products, patterns)
        layouts = {}
        for junction in network.intersections:
            try:
                layouts[junction.key] = IntersectionProvider.footprint(context, junction)
            except ValueError:
                continue  # Reported by the full footprint/import qualifier.
        for group in CompoundPlanner().build(layouts):
            record = dict(network=network.uuid, members=[key.text() for key in group.members])
            start = time.perf_counter()
            try:
                result = CompoundQuadChart(patterns, spacing=network.settings.spacing).build(group, context, require_quality=False)
                sampled_group = result.partition.group
                conditions = [BoundaryConditions().build(sampled_group, loop, layouts) for loop in sampled_group.loops]
                sources = {p.region.owner: p for p in sampled_group.footprints}
                chord_error = 0.
                for loop in sampled_group.loops:
                    for at, span in enumerate(loop.spans):
                        label = rim_label(span, sources)
                        if label.kind != 'curb':
                            continue
                        curve = layouts[span.owner].boundaries[label.index]
                        fractions = np.linspace(0., 1., 33)
                        points = curve(label.first + fractions * (label.last - label.first))
                        chords = loop.xy[at] + fractions[:, None] * (loop.xy[(at + 1) % len(loop.xy)] - loop.xy[at])
                        chord_error = max(chord_error, float(np.linalg.norm(points - chords, axis=1).max()))
                assert chord_error <= network.settings.curve_tolerance + 1e-7
                original_conditions = [BoundaryConditions().build(group, loop, layouts) for loop in group.loops]
                ports = {(p.section_id, p.end): (b, p) for b in conditions for p in b.ports}
                for boundary in original_conditions:
                    for port in boundary.ports:
                        new_boundary, new_port = ports[port.section_id, port.end]
                        np.testing.assert_array_equal(boundary.xy[list(port.vertices)], new_boundary.xy[list(new_port.vertices)])
                        np.testing.assert_array_equal(boundary.heights[list(port.vertices)], new_boundary.heights[list(new_port.vertices)])
                for node, ref in result.boundary_refs.items():
                    xy, height, gradient = conditions[ref.loop].edge_samples(ref.edge, ref.fraction)
                    assert np.linalg.norm(result.mesh.xyz[node, :2] - xy) < 1e-7
                meshes = {model.Key('compound', 0): result.mesh}
                meshes.update({model.Key('hole', li): hole.mesh for li, hole in result.holes.items()})
                overlap = validate_overlaps(meshes, cell_size=max(1., network.settings.spacing * 2.))
                maximum = max(mesh.validate()['maximum_edge_ratio'] for mesh in meshes.values())
                record.update(status='PLANAR_QUALIFIED' if maximum <= 8. else 'NEEDS_GRADING',
                    total_quads=result.metadata['total_quads'], pavement_quads=len(result.mesh.quads),
                    grass_quads={str(li): len(hole.mesh.quads) for li, hole in result.holes.items()},
                    blocks=len(result.partition.blocks), maximum_edge_ratio=maximum,
                    area_error=result.partition.area_error, boundary_condition_references=True,
                    maximum_sampled_curve_error=chord_error, full_external_profiles_preserved=True,
                    sampled_boundary_vertices=sum(len(loop.xy) for loop in sampled_group.loops),
                    boundary_area_change=sum(loop.area for loop in sampled_group.loops) - sum(loop.area for loop in group.loops),
                    overlap_check=overlap, block_reports=result.metadata['blocks'],
                    cut_intervals=result.samples.cut_intervals)
                pictures.append(('Sites ' + ', '.join(str(key.id) for key in group.members), result, record))
            except ValueError as exc:
                record.update(status='UNSUPPORTED', reason=str(exc))
            record['seconds'] = time.perf_counter() - start
            records.append(record)
            print(json.dumps({key: value for key, value in record.items() if key not in ('block_reports', 'cut_intervals')}), flush=True)
    report = dict(implementation=VERSION, planar_only=True, pavement_coverage_asserted=False, cases=records)
    (args.output / 'compound_chart_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    draw(pictures, args.output / 'compound_quad_charts.svg')


if __name__ == '__main__':
    main()
