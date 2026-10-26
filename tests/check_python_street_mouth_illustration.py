"""Illustrate a compatible-seam remnant; reference disables only absorption."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _python_street_package
from _python_street_fixtures import terminal_fixture_networks
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.intersections import IntersectionProvider


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference-road', type=Path)
    parser.add_argument('--result-junction', type=Path)
    parser.add_argument('--profile-intervals', type=int, default=16)
    args = parser.parse_args()
    if bool(args.reference_road) != bool(args.result_junction) or args.profile_intervals < 1:
        parser.error('Provide both real-mesh captures and a positive profile interval count')
    args.output.mkdir(parents=True, exist_ok=True)
    network, = terminal_fixture_networks(model)
    end = network.sections[0].points[-1]
    network = replace(network, sections=(*network.sections,
        model.Section(4, 1, (end, (30., 0., 0.)), width=5.4)),
        connections=(model.Connection(99, (model.End(1, 1), model.End(4, 0))),))
    saved = IntersectionProvider.mouth_station
    try:
        IntersectionProvider.mouth_station = staticmethod(lambda approach, station, settings: station)
        before = PythonMeshBackend(network.settings).prepare(network)
    finally:
        IntersectionProvider.mouth_station = staticmethod(saved)
    after = PythonMeshBackend(network.settings).prepare(network)
    strip = before.products[model.Key('road', 1)].mesh
    for sid in (2, 3, 4):
        a, b = before.products[model.Key('road', sid)].mesh, after.products[model.Key('road', sid)].mesh
        np.testing.assert_array_equal(a.xyz, b.xyz)
        np.testing.assert_array_equal(a.quads, b.quads)
    assert getattr(after.products[model.Key('road', 1)], 'mesh', None) is None
    rows = strip.xyz.reshape(-1, 17, 3)
    cap = rows[-1]
    junction = after.products[model.Key('junction_mesh', 10)].mesh
    distance = np.linalg.norm(junction.xyz[:, None] - cap[None], axis=-1)
    maximum_error = float(distance.min(axis=0).max())
    assert maximum_error < 1e-8
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.titleweight': 'bold', 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
    for ax, update, title, color in zip(axes[0], (before, after),
            ('Before: a separate 18 mm road remnant', 'After: the full profile meets the Section seam'),
            ('#e87926', '#16836b')):
        for key, part in update.products.items():
            mesh = getattr(part, 'mesh', None)
            if mesh is None:
                continue
            highlight = key == model.Key('road', 1) if update is before else key == model.Key('junction_mesh', 10)
            ax.add_collection(PolyCollection(mesh.xyz[mesh.quads, :2], linewidths=.4,
                facecolors=color if highlight else '#dce7ee', edgecolors='#6b8291', alpha=.85))
        ax.plot((end[0], end[0]), (-2.7, 2.7), color='#7939a8', lw=1.5, ls='--', label='Authored Section seam')
        ax.set(xlim=(-7., 6.), ylim=(-3.5, 7.), aspect='equal', title=title, xlabel='Metres', ylabel='Metres')
        ax.annotate('Same circular curbs', xy=(3.29, 3.29), xytext=(-1.8, 5.8),
                    arrowprops={'arrowstyle': '->', 'color': '#314b61'}, color='#314b61')
        ax.legend(loc='lower left', fontsize=9)
    ax = axes[1, 0]
    xy = strip.xyz[:, :2].copy()
    xy[:, 0] = (xy[:, 0] - rows[0, 0, 0]) * 1000.
    ax.add_collection(PolyCollection(xy[strip.quads], facecolors='#fff0dc', edgecolors='#d76a1c', linewidths=.8))
    ax.autoscale()
    ax.set(title='Why this strip needs ownership adjustment', xlabel='Along the road (millimetres)',
           ylabel='Across the road (metres)')
    ax.text(.5, .02, 'Axes use different units to expose the 32 thin cells.',
            transform=ax.transAxes, ha='center', fontsize=9, color='#6f491d')
    ax = axes[1, 1]
    closest = junction.xyz[distance.argmin(axis=0)]
    ax.plot(cap[:, 1], cap[:, 2] / network.settings.exaggeration, 'o', color='#e87926',
            ms=5, label='Complete road profile before')
    ax.plot(closest[:, 1], closest[:, 2] / network.settings.exaggeration, '-', color='#16836b',
            lw=2, label='Junction mouth after')
    ax.set(title='The crown is retained at the same source seam', xlabel='Transverse distance (metres)',
           ylabel='Physical elevation (metres)', ylim=(-.005, .09))
    ax.grid(alpha=.2)
    ax.legend(fontsize=9)
    fig.suptitle('A tiny compatible-Section remnant becomes junction-owned pavement', fontsize=17, color='#173a50')
    fig.savefig(args.output / 'mouth_alignment.svg')
    fig.savefig(args.output / 'mouth_alignment.png', dpi=150)
    report = dict(passed=True, source_sections=[s.id for s in network.sections],
        reference_policy='Absorption disabled; all other construction unchanged',
        remnant_length=float(rows[-1, 0, 0] - rows[0, 0, 0]), remnant_quads=len(strip.quads),
        remnant_maximum_edge_ratio=strip.validate()['maximum_edge_ratio'],
        maximum_profile_error=maximum_error, continuation_exact=True,
        after_junction_quality=junction.validate())
    if args.reference_road:
        with np.load(args.reference_road) as saved_road, np.load(args.result_junction) as saved_junction:
            original = saved_road['xyz'][-(args.profile_intervals + 1):]
            result = saved_junction['xyz']
            error = float(np.linalg.norm(result[:, None] - original[None], axis=-1).min(axis=0).max())
            assert error <= 2e-5, ('Real terminal profile mismatch', error)
            report.update(real_reference_road=str(args.reference_road.resolve()),
                          real_result_junction=str(args.result_junction.resolve()),
                          real_profile_error=error, real_profile_vertices=len(original))
    (args.output / 'mouth_checks.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report))


if __name__ == '__main__':
    main()
