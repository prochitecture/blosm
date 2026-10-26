"""Scoped corner-strip/core prototype on effective numerical source inputs."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import _python_street_package
from _python_corner_strips import CornerStrips
from check_python_street_chart_failure import context_for
from tested_python_streets import model
from tested_python_streets.ownership import footprint_key
from tested_python_streets.validation import validate_overlaps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--junction', type=int, required=True)
    parser.add_argument('--collar', type=int, action='append', default=[])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    network, = [n for n in map(model.Network.from_dict, json.loads(args.input.read_text()))
                if any(j.id == args.junction for j in n.intersections)]
    site = next(j for j in network.intersections if j.id == args.junction)
    context = context_for(network)
    layout = context[footprint_key(site.key)]
    results, drawings = [], []
    for collar in args.collar or (2, 3, 4):
        start = time.perf_counter()
        try:
            mesh, report, ring = CornerStrips(context.patterns).build(context, site, layout, collar_intervals=collar)
            report['overlaps'] = validate_overlaps({site.key: mesh}, cell_size=1.4)
            report.update(collar_intervals=collar, status='SCOPED_SURFACE', seconds=time.perf_counter() - start)
            np.savez_compressed(args.output / f'collar_{collar}.npz', xyz=mesh.xyz, quads=mesh.quads, ring=ring)
            drawings.append((collar, mesh, ring))
        except ValueError as exc:
            report = dict(collar_intervals=collar, status='FAILED', reason=str(exc), seconds=time.perf_counter() - start)
            if hasattr(exc, 'ring'):
                np.savez_compressed(args.output / f'collar_{collar}_failed_core.npz', boundary=exc.ring,
                                    xyz=exc.bands.xyz, quads=exc.bands.quads)
                drawings.append((collar, exc.bands, exc.ring))
        results.append(report)
        print(json.dumps(report), flush=True)
    (args.output / 'strip_checks.json').write_text(json.dumps(dict(junction=site.id,
        experimental_only=True, full_network_qualified=False, cases=results), indent=2))
    if drawings:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        fig, axes = plt.subplots(1, len(drawings), figsize=(8. * len(drawings), 8.), squeeze=False)
        fig.suptitle('Corner strips and finite core mouths — scoped experiment', fontweight='bold')
        for ax, (collar, mesh, ring) in zip(axes[0], drawings):
            ax.add_collection(PolyCollection(mesh.xyz[mesh.quads, :2], facecolors='#dceafa', edgecolors='#426c87', linewidths=.3))
            ax.plot(*np.vstack((ring, ring[:1])).T, color='#c43a39', lw=1.3)
            ax.autoscale_view();ax.set_aspect('equal');ax.grid(alpha=.2)
            ax.set_title(f'Collar {collar}: {len(mesh.quads)} quads; red = core boundary')
            ax.set_xlabel('X (m)');ax.set_ylabel('Y (m)')
        fig.tight_layout()
        fig.savefig(args.output / 'corner_strips.svg')
        fig.savefig(args.output / 'corner_strips.png', dpi=150)


if __name__ == '__main__':
    main()
