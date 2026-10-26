"""Reproduce physical junction extents separately from surface guide choices.

Input is numerical_inputs.json exported by check_python_street_import.py.
This is a diagnostic, never a mesh acceptance or full-network coverage test.
"""
import argparse
import json
from pathlib import Path
import _python_street_package
from tested_python_streets import VERSION, model
from tested_python_streets.compounds import CompoundPlanner
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.updates import BuildContext, FeatureProvider, UpdateCoordinator


class ApproachInputs(FeatureProvider):
    def tasks(self, network):
        return (task for task in IntersectionProvider().tasks(network) if task.key.kind == 'approach')


def inspect(network):
    patterns = PatternLibrary(network.settings.pattern_cache_bytes)
    coordinator = UpdateCoordinator((PathProvider(), ProfileProvider(), ApproachInputs()), patterns)
    prepared = coordinator.prepare(network)
    context = BuildContext(network, prepared.products, patterns)
    layouts, records = {}, []
    for junction in network.intersections:
        record = dict(id=junction.id, sections=[end.section_id for end in junction.ends])
        try:
            layout = IntersectionProvider.footprint(context, junction)
            layouts[junction.key] = layout
            record.update(footprint='SOLVED', approaches=[dict(section=a.source.id, end=a.end,
                route=[p.section.id for p in a.path.pieces], mouth=distance, available=a.length)
                for a, distance in zip(layout.approaches, layout.outward_stations)],
                trims=[dict(section=sid, end=end, station=mouth.station)
                    for (sid, end), mouth in sorted(layout.mouths.items())],
                corners=[dict(first=c.first.source.id, second=c.second.source.id,
                    radius=c.radius, contact_stations=c.stations) for c in layout.corners])
        except ValueError as exc:
            record.update(footprint='FAILED', reason=str(exc))
        if junction.key in layouts:
            try:
                IntersectionProvider.layout(context, junction)
                record['guides'] = 'SOLVED'
            except ValueError as exc:
                record.update(guides='FAILED', guide_reason=str(exc))
        records.append(record)
    # Flag source intervals claimed by two real sites. These require ownership
    # coordination, not a tangent extension through an actual intersection.
    claims = {}
    for key, layout in layouts.items():
        for (sid, end), mouth in layout.mouths.items():
            claims.setdefault(sid, {0: [], 1: []})[end].append((key, mouth.station))
    conflicts = []
    for sid, ends in claims.items():
        for first, lo in ends[0]:
            for second, hi in ends[1]:
                if first != second and lo > hi + 1e-7:
                    conflicts.append(dict(section=sid, junctions=[first.text(), second.text()],
                        first_trim=lo, last_trim=hi, overlap=lo - hi))
    report = dict(network=network.uuid, sections=len(network.sections), junctions=records,
                  overlapping_source_claims=conflicts)
    try:
        groups = CompoundPlanner().build(layouts)
        report['compound_boundary_status'] = 'SOLVED'
        report['compound_groups'] = []
        for group in groups:
            record = dict(members=[key.text() for key in group.members],
                internal_sections=sorted({c.section_id for c in group.claims}),
                footprints=[p.region.xy.tolist() for p in group.footprints],
                loops=[dict(hole=loop.hole, area=loop.area, xy=loop.xy.tolist(),
                    source_owners=sorted({span.owner.text() for span in loop.spans})) for loop in group.loops])
            try:
                data = [BoundaryConditions().build(group, loop, layouts) for loop in group.loops]
                record.update(boundary_conditions='SOLVED', datum=data[0].datum, ceiling=data[0].ceiling,
                    ports=[dict(section=port.section_id, end=port.end, station=port.station, samples=len(port.vertices),
                        xy=boundary.xy[list(port.vertices)].tolist()) for boundary in data for port in boundary.ports],
                    boundary_ranges=[dict(hole=boundary.hole, minimum=float(boundary.heights.min()),
                                          maximum=float(boundary.heights.max())) for boundary in data])
            except ValueError as exc:
                record.update(boundary_conditions='FAILED', boundary_condition_reason=str(exc))
            report['compound_groups'].append(record)
    except ValueError as exc:
        report.update(compound_boundary_status='FAILED', compound_boundary_reason=str(exc))
    return report


def draw(reports, destination):
    from html import escape
    groups = [group for report in reports for group in report.get('compound_groups', ())]
    if not groups:
        return
    width, height = 1120, 120 + 390 * ((len(groups) + 1) // 2)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f3f6fa"/>',
        '<g font-family="Segoe UI, sans-serif" fill="#17324d">',
        '<text x="28" y="36" font-size="25" font-weight="600">Compound boundary planning</text>',
        '<text x="28" y="64" font-size="16">Blue: combined footprint. Green: preserved hole. Grey: original junction boundaries.</text>',
        '<text x="28" y="88" font-size="16">Preparation only — no compound pavement mesh has been accepted.</text>']
    for index, group in enumerate(groups):
        x, y = 20 + 550 * (index % 2), 110 + 390 * (index // 2)
        points = [point for loop in group['loops'] for point in loop['xy']]
        lo = [min(p[axis] for p in points) for axis in (0, 1)]
        hi = [max(p[axis] for p in points) for axis in (0, 1)]
        scale = min(470. / (hi[0] - lo[0]), 285. / (hi[1] - lo[1]))
        center = [.5 * (lo[k] + hi[k]) for k in (0, 1)]
        def path(points):
            mapped = [(x + 265 + scale * (p[0] - center[0]), y + 202 - scale * (p[1] - center[1])) for p in points]
            return 'M' + ' L'.join(f'{a:.3f},{b:.3f}' for a, b in mapped) + ' Z'
        svg.extend([f'<rect x="{x}" y="{y}" width="530" height="370" rx="14" fill="white" stroke="#cbd5e1"/>',
            f'<text x="{x + 18}" y="{y + 27}" font-size="17" font-weight="600">Sites {escape(", ".join(k.split(":")[1] for k in group["members"]))}</text>'])
        for loop in group['loops']:
            fill = '#b9d8ad' if loop['hole'] else '#dceafa'
            stroke = '#53884b' if loop['hole'] else '#285a89'
            svg.append(f'<path d="{path(loop["xy"])}" fill="{fill}" stroke="{stroke}" stroke-width="2"/>')
        for polygon in group['footprints']:
            svg.append(f'<path d="{path(polygon)}" fill="none" stroke="#8b97a3" stroke-width=".8" stroke-dasharray="4 3"/>')
        for port in group.get('ports', ()):
            svg.append(f'<path d="{path(port["xy"])}" fill="none" stroke="#172b40" stroke-width="2.8"/>')
        svg.append(f'<text x="{x + 18}" y="{y + 351}" font-size="14">{len(group["members"])} sites; {len(group["internal_sections"])} overlapping spans; {sum(loop["hole"] for loop in group["loops"])} holes</text>')
    svg.append('</g></svg>')
    destination.write_text('\n'.join(svg), encoding='utf8')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    inputs = json.loads(args.input.read_text(encoding='utf8'))
    reports = [inspect(model.Network.from_dict(data)) for data in inputs]
    result = dict(implementation=VERSION, diagnostic_only=True, input=str(args.input.resolve()), networks=reports)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'footprint_checks.json').write_text(json.dumps(result, indent=2), encoding='utf8')
    draw(reports, args.output / 'compound_boundaries.svg')
    for report in reports:
        records = report['junctions']
        print(json.dumps(dict(network=report['network'], sections=report['sections'],
            footprints_solved=sum(r['footprint'] == 'SOLVED' for r in records), junctions=len(records),
            guide_failures=[r['id'] for r in records if r.get('guides') == 'FAILED'],
            overlapping_claims=len(report['overlapping_source_claims']),
            compound_boundary_status=report['compound_boundary_status'],
            compound_groups=len(report.get('compound_groups', ())))), flush=True)


if __name__ == '__main__':
    main()
