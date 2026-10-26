"""Bounded planar-chart qualification, separate from pavement surface coverage.

The optional input is numerical_inputs.json from the normal-import check.
No heights, separators, compound decomposition or Blender objects are generated.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import _python_street_package
from _python_street_fixtures import fixture_networks
from check_python_street_footprints import ApproachInputs
from tested_python_streets import VERSION, model
from tested_python_streets.compounds import sample_footprint
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.mesh import cross
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.quad_charts import LayeredQuadChart
from tested_python_streets.updates import BuildContext, UpdateCoordinator
from tested_python_streets.validation import validate_overlaps


def draw(cases, destination):
    from html import escape
    width, height = 1140, 108 + 355 * ((len(cases) + 2) // 3)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f3f6fa"/>', '<g font-family="Segoe UI, sans-serif" fill="#17324d">',
        '<text x="24" y="35" font-size="25" font-weight="600">Direct quad charts — planar preparation</text>',
        '<text x="24" y="64" font-size="16">Physical rim fixed. Cached 3:1 blocks terminate fine sampling inside the mesh.</text>',
        '<text x="24" y="88" font-size="16">No triangle-first stage. Height reconciliation and compound surfaces are still separate work.</text>']
    for index, (name, mesh, report) in enumerate(cases):
        x, y = 15 + 375 * (index % 3), 105 + 355 * (index // 3)
        lo, hi = mesh.xyz[:, :2].min(axis=0), mesh.xyz[:, :2].max(axis=0)
        center = .5 * (lo + hi)
        scale = min(315. / (hi[0] - lo[0]), 245. / (hi[1] - lo[1]))
        def path(points):
            return 'M' + ' L'.join(f'{x + 180 + scale * (p[0] - center[0]):.3f},{y + 176 - scale * (p[1] - center[1]):.3f}' for p in points) + ' Z'
        svg.extend([f'<rect x="{x}" y="{y}" width="360" height="340" rx="12" fill="white" stroke="#cad6e1"/>',
            f'<text x="{x + 14}" y="{y + 25}" font-size="17" font-weight="600">{escape(name)}</text>'])
        for face in mesh.quads:
            svg.append(f'<path d="{path(mesh.xyz[face, :2])}" fill="none" stroke="#688eaa" stroke-width=".55"/>')
        svg.append(f'<path d="{path(mesh.xyz[mesh.seams[0], :2])}" fill="none" stroke="#244c70" stroke-width="1.7"/>')
        svg.append(f'<text x="{x + 14}" y="{y + 321}" font-size="14">{len(mesh.quads)} quads; max edge ratio {report["quality"]["maximum_edge_ratio"]:.2f}</text>')
    svg.append('</g></svg>')
    destination.write_text('\n'.join(svg), encoding='utf8')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    networks = list(fixture_networks(model))
    if args.input:
        networks.extend(model.Network.from_dict(data) for data in json.loads(args.input.read_text(encoding='utf8')))
    args.output.mkdir(parents=True, exist_ok=True)
    records, pictures = [], []
    for network in networks:
        patterns = PatternLibrary(network.settings.pattern_cache_bytes)
        prepared = UpdateCoordinator((PathProvider(), ProfileProvider(), ApproachInputs()), patterns).prepare(network)
        context = BuildContext(network, prepared.products, patterns)
        chart = LayeredQuadChart(patterns, spacing=network.settings.spacing)
        for junction in network.intersections:
            record = dict(network=network.uuid, site=junction.id)
            start = time.perf_counter()
            try:
                layout = IntersectionProvider.footprint(context, junction)
                rim = sample_footprint(junction.key, layout).region.xy
                result = chart.build(rim, require_quality=False)
                mesh = result.mesh
                assert np.array_equal(mesh.xyz[mesh.seams[0], :2], rim)
                edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
                unique, counts = np.unique(edges, axis=0, return_counts=True)
                assert len(mesh.xyz) - len(unique) + len(mesh.quads) == 1
                assert set(counts).issubset({1, 2}) and np.sum(counts == 1) == len(rim)
                overlap = validate_overlaps({junction.key: mesh}, cell_size=max(1., network.settings.spacing * 2.))
                points = mesh.xyz[mesh.quads, :2]
                edge_lengths = np.linalg.norm(np.roll(points, -1, axis=1) - points, axis=-1)
                areas = .5 * np.sum(cross(points, np.roll(points, -1, axis=1)), axis=1)
                qualified = (result.metadata['quality']['maximum_edge_ratio'] <= 8.
                             and result.metadata['quality']['minimum_corner_sine'] >= .03)
                record.update(status='PLANAR_QUALIFIED' if qualified else 'NEEDS_GRADING', **result.metadata,
                    fixed_rim=True, quad_disk=True, overlap_check=overlap,
                    edge_length_percentiles=np.percentile(edge_lengths, (1, 50, 99)).tolist(),
                    quad_area_percentiles=np.percentile(areas, (1, 50, 99)).tolist())
                if network.uuid in ('T', 'three_way', 'oblique', 'five_way', 'unequal_widths', 'divided_T'):
                    name = network.uuid + (f' / {junction.id}' if network.uuid == 'divided_T' else '')
                    pictures.append((name, mesh, record))
            except ValueError as exc:
                record.update(status='UNSUPPORTED', reason=str(exc))
            record['seconds'] = time.perf_counter() - start
            records.append(record)
    report = dict(implementation=VERSION, geometry_only=True, pavement_coverage_asserted=False,
        cases=records, summary={status: sum(r['status'] == status for r in records)
            for status in sorted({r['status'] for r in records})})
    (args.output / 'quad_chart_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    draw(pictures, args.output / 'quad_charts.svg')
    print(json.dumps(dict(geometry_only=True, summary=report['summary'])), flush=True)


if __name__ == '__main__':
    main()
