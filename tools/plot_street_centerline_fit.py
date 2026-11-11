"""Draw selected saved pavement meshes from inspect_street_centerlines reports."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, required=True)
    parser.add_argument('--after', type=Path, required=True)
    parser.add_argument('--section', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = [json.loads(path.read_text()) for path in (args.before, args.after)]
    section = next(s for s in data[0]['sections'] if s['section'] == args.section)
    source = np.asarray(section['source'])[:, :2]
    margin = 1.5 * section['width']
    low, high = source.min(axis=0)-margin, source.max(axis=0)+margin
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), layout='constrained')
    fig.set_facecolor('#f5f7fa')
    for ax, report, title in zip(axes, data, ('Before', 'After')):
        ax.set_facecolor('#eef2f5')
        mesh = next(m for m in report['pavements'] if m['key'] == f'road:{args.section}')
        points = np.asarray(mesh['xy'])
        faces = [points[face] for face in mesh['faces']]
        ax.add_collection(PolyCollection(faces, facecolor='#7a99aa', edgecolor='none', antialiaseds=False, rasterized=True))
        ax.plot(*source.T, color='#172b3a', linewidth=1.7, label='Unchanged input centerline')
        ax.scatter(*source.T, color='#e19b35', edgecolor='#ffffff', s=30, zorder=4)
        ax.set_title(title, loc='left', fontweight='bold', color='#172b3a', fontsize=13)
        ax.set(xlim=(low[0], high[0]), ylim=(low[1], high[1]), aspect='equal', ylabel='Metres')
        ax.grid(color='#ffffff', linewidth=.8, zorder=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    axes[0].legend(loc='upper left', frameon=False)
    axes[1].set_xlabel('Metres')
    fig.suptitle('Saved carriageway mesh and source alignment', color='#172b3a', fontsize=17, fontweight='bold')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == '__main__':
    main()
