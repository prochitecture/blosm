"""Capture a real planar chart gate and compare existing computational cuts.

Diagnostic only: unqualified cells are drawn explicitly, never published.
Optional strip cuts reuse an existing policy with the same physical boundary.
"""
import argparse
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
import _python_street_package
from check_python_street_footprints import ApproachInputs
from tested_python_streets import model, VERSION
from tested_python_streets.compounds import CompoundGroup, sample_footprint
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.ownership import footprint_key, OwnershipPlan
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.partition import CompoundPartitioner, CompoundQuadChart
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.quad_charts import LayeredQuadChart
from tested_python_streets.regions import RegionUnion
from tested_python_streets.sampling import AnalyticalBoundarySampler
from tested_python_streets.separator_surfaces import SeparatorStripCuts
from tested_python_streets.updates import UpdateCoordinator, BuildContext
from tested_python_streets.decomposition import VisibilitySplitter


class AngleSplitter(VisibilitySplitter):
    def __init__(self, minimum_angle, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.minimum_angle = np.deg2rad(minimum_angle)

    def cut_angle(self, polygon, first, last):
        angle = super().cut_angle(polygon, first, last)
        return angle if angle >= self.minimum_angle else 0.


def context_for(network):
    patterns = PatternLibrary(network.settings.pattern_cache_bytes)
    update = UpdateCoordinator((PathProvider(), ProfileProvider(), ApproachInputs()), patterns).prepare(network)
    context = BuildContext(network, update.products, patterns)
    for site in network.intersections:
        context.products[footprint_key(site.key)] = IntersectionProvider.footprint(context, site)
    return context


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--junction', type=int, required=True)
    parser.add_argument('--try-strip-cuts', action='store_true')
    parser.add_argument('--minimum-cut-angle', type=float, action='append', default=[])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    networks = [model.Network.from_dict(data) for data in json.loads(args.input.read_text(encoding='utf8'))]
    network, = [n for n in networks if any(j.id == args.junction for j in n.intersections)]
    site = next(j for j in network.intersections if j.id == args.junction)
    context = context_for(network)
    plan = OwnershipPlan.build(context, network.intersections)
    if site.key in plan.by_site:
        group = plan.by_site[site.key]
    else:
        footprint = sample_footprint(site.key, context[footprint_key(site.key)])
        group = CompoundGroup((site.key,), (), (footprint,), RegionUnion().build((footprint.region,)))
    original_group = group
    group = AnalyticalBoundarySampler(tolerance=network.settings.curve_tolerance,
                                     spacing=network.settings.spacing).build(group)
    partition = CompoundPartitioner().build(group, context)
    samples = partition.sample(network.settings.spacing)
    args.output.mkdir(parents=True, exist_ok=True)
    report = dict(implementation=VERSION, junction=site.id, planar_diagnostic_only=True,
                  source_input=str(args.input.resolve()), blocks=[], strip_candidate=None,
                  full_pavement_qualified=False)
    report['angle_candidates'] = []
    report['initial_cuts'] = [dict(first=c.first.edge, last=c.last.edge,
        minimum_angle_degrees=float(np.rad2deg(VisibilitySplitter.cut_angle(group.loops[0].xy,
            c.first.edge, c.last.edge)))) for c in partition.cuts if c.first.loop == c.last.loop == 0]
    pictures = []
    for bi, rim in enumerate(samples.rims):
        boundary = samples.xy[list(rim)]
        result = LayeredQuadChart(context.patterns, spacing=network.settings.spacing).build(boundary, require_quality=False)
        xyz, faces = result.mesh.xyz, result.mesh.quads
        lengths = np.linalg.norm(np.roll(xyz[faces, :2], -1, axis=1) - xyz[faces, :2], axis=-1)
        ratios = lengths.max(axis=1) / lengths.min(axis=1)
        np.savez_compressed(args.output / f'block_{bi}.npz', boundary=boundary, xyz=xyz, quads=faces)
        record = dict(block=bi, boundary_samples=len(boundary), cuts=len(partition.cuts),
                      cells_above_ratio_8=int(np.sum(ratios > 8.)), maximum_ratio=float(ratios.max()),
                      construction=result.metadata)
        report['blocks'].append(record)
        pictures.append((f'Existing block {bi}', result.mesh, ratios))
    if args.try_strip_cuts:
        try:
            chart = CompoundQuadChart(context.patterns, spacing=network.settings.spacing,
                partitioner=CompoundPartitioner(policy=SeparatorStripCuts(network.settings.spacing))).build(original_group, context)
            report['strip_candidate'] = dict(planar_qualified=True, quads=len(chart.mesh.quads),
                cuts=len(chart.partition.cuts), blocks=len(chart.partition.blocks), quality=chart.mesh.validate(),
                area_error=chart.partition.area_error)
            np.savez_compressed(args.output / 'strip_candidate.npz', xyz=chart.mesh.xyz, quads=chart.mesh.quads)
            lengths = np.linalg.norm(np.roll(chart.mesh.xyz[chart.mesh.quads, :2], -1, axis=1)
                                     - chart.mesh.xyz[chart.mesh.quads, :2], axis=-1)
            pictures.append(('Existing strip-cut policy candidate', chart.mesh, lengths.max(axis=1) / lengths.min(axis=1)))
        except ValueError as exc:
            report['strip_candidate'] = dict(planar_qualified=False, reason=str(exc))
    for angle in args.minimum_cut_angle:
        try:
            with patch('tested_python_streets.decomposition.VisibilitySplitter',
                       lambda *a, **kw: AngleSplitter(angle, *a, **kw)):
                chart = CompoundQuadChart(context.patterns, spacing=network.settings.spacing).build(original_group, context)
            record = dict(minimum_angle_degrees=angle, planar_qualified=True, quads=len(chart.mesh.quads),
                cuts=len(chart.partition.cuts), blocks=len(chart.partition.blocks), quality=chart.mesh.validate(),
                area_error=chart.partition.area_error)
            np.savez_compressed(args.output / f'angle_{angle:g}_candidate.npz', xyz=chart.mesh.xyz, quads=chart.mesh.quads)
            lengths = np.linalg.norm(np.roll(chart.mesh.xyz[chart.mesh.quads, :2], -1, axis=1)
                                     - chart.mesh.xyz[chart.mesh.quads, :2], axis=-1)
            pictures.append((f'Minimum cut angle {angle:g}° candidate', chart.mesh, lengths.max(axis=1) / lengths.min(axis=1)))
        except ValueError as exc:
            record = dict(minimum_angle_degrees=angle, planar_qualified=False, reason=str(exc))
        report['angle_candidates'].append(record)
    (args.output / 'chart_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    fig, axes = plt.subplots(1, len(pictures), figsize=(8. * len(pictures), 8.), squeeze=False)
    fig.suptitle(f'Junction {site.id}: planar chart diagnostic — no pavement accepted', fontsize=16, fontweight='bold')
    for ax, (title, mesh, ratio) in zip(axes[0], pictures):
        color = np.where((ratio > 8.)[:, None], np.array((.92, .32, .27, .65)), np.array((.76, .86, .94, .75)))
        ax.add_collection(PolyCollection(mesh.xyz[mesh.quads, :2], facecolors=color, edgecolors='#426c87', linewidths=.35))
        for s in network.sections:
            points = np.asarray(s.points)
            ax.plot(points[:, 0], points[:, 1], color='#25333d', lw=.7, alpha=.6)
        ax.autoscale_view()
        lo, hi = mesh.xyz[:, :2].min(axis=0), mesh.xyz[:, :2].max(axis=0)
        margin = .08 * max(hi - lo)
        ax.set_xlim(lo[0] - margin, hi[0] + margin)
        ax.set_ylim(lo[1] - margin, hi[1] + margin)
        ax.set_aspect('equal')
        ax.set_title(f'{title}\n{len(mesh.quads)} quads; {int(np.sum(ratio > 8.))} red cells above edge ratio 8')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.grid(alpha=.2)
    fig.tight_layout(rect=(0., 0., 1., .95))
    fig.savefig(args.output / 'chart_failure.svg')
    fig.savefig(args.output / 'chart_failure.png', dpi=150)
    plt.close(fig)
    print(json.dumps(dict(junction=site.id, blocks=len(report['blocks']),
                         maximum_ratio=max(r['maximum_ratio'] for r in report['blocks']),
                         initial_cuts=report['initial_cuts'], angle_candidates=report['angle_candidates'],
                         strip_candidate=report['strip_candidate'], full_pavement_qualified=False)), flush=True)


if __name__ == '__main__':
    main()
