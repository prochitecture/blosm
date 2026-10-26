"""Inspect two captured, unpublished surfaces with independent polygon clipping."""
import argparse
import json
from pathlib import Path
import numpy as np


def area(polygon):
    if len(polygon) < 3:
        return 0.
    return .5 * abs(float(np.sum(np.cross(polygon, np.roll(polygon, -1, axis=0)))))


def clip(subject, boundary):
    result = np.asarray(subject)
    for a, b in zip(boundary, np.roll(boundary, -1, axis=0)):
        if not len(result):
            break
        values = np.cross(b - a, result - a)
        points = []
        for i, p in enumerate(result):
            q = result[(i + 1) % len(result)]
            first, last = values[i], values[(i + 1) % len(result)]
            if first >= 0.:
                points.append(p)
            if (first >= 0.) != (last >= 0.):
                points.append(p + first / (first - last) * (q - p))
        result = np.asarray(points)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--first', required=True)
    parser.add_argument('--second', required=True)
    parser.add_argument('--input', type=Path, required=True)
    args = parser.parse_args()
    sources = [np.load(args.directory / (name + '.npz')) for name in (args.first, args.second)]
    faces = [source['xyz'][source['quads'], :2] for source in sources]
    bounds = [(p.min(axis=1), p.max(axis=1)) for p in faces]
    intersections, pairs = [], []
    for i, p in enumerate(faces[0]):
        ids = np.flatnonzero(np.all(bounds[1][1] > bounds[0][0][i], axis=1) &
                            np.all(bounds[1][0] < bounds[0][1][i], axis=1))
        for j in ids:
            polygon = clip(p, faces[1][j])
            value = area(polygon)
            if value > 1e-10:
                intersections.append(polygon)
                pairs.append(dict(first=i, second=int(j), area=value))
    report = dict(first=args.first, second=args.second, published=False,
                  overlap_area=sum(p['area'] for p in pairs), pairs=pairs)
    (args.directory / 'overlap_detail.json').write_text(json.dumps(report, indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    fig, ax = plt.subplots(figsize=(11, 9))
    for polys, color, label in zip(faces, ('#416f98', '#db8949'), (args.first, args.second)):
        ax.add_collection(PolyCollection(polys, facecolors=color, edgecolors=color, alpha=.28, linewidths=.5))
        ax.plot([], [], color=color, linewidth=5, label=label)
    ax.add_collection(PolyCollection(intersections, facecolors='#e7373d', edgecolors='#ab1724', linewidths=.4))
    lo, hi = faces[1].reshape(-1, 2).min(axis=0), faces[1].reshape(-1, 2).max(axis=0)
    margin = .15 * max(hi - lo)
    for network in json.loads(args.input.read_text()):
        for section in network['sections']:
            xy = np.asarray(section['points'])[:, :2]
            if np.any(np.all(xy >= lo - margin, axis=1) & np.all(xy <= hi + margin, axis=1)):
                ax.plot(xy[:, 0], xy[:, 1], 'k.-', linewidth=.8, markersize=3)
                ax.text(*xy[len(xy) // 2], str(section['id']), fontsize=9)
    ax.set_xlim(lo[0] - margin, hi[0] + margin)
    ax.set_ylim(lo[1] - margin, hi[1] + margin)
    ax.set_aspect('equal'); ax.grid(alpha=.2); ax.legend()
    ax.set_title(f'Unpublished overlap diagnosis: {report["overlap_area"]:.6f} m² in red')
    ax.set_xlabel('X (m)'); ax.set_ylabel('Y (m)')
    fig.tight_layout()
    fig.savefig(args.directory / 'overlap_detail.svg')
    fig.savefig(args.directory / 'overlap_detail.png', dpi=130)
    plt.close(fig)
    print(json.dumps({key: value for key, value in report.items() if key != 'pairs'}), flush=True)


if __name__ == '__main__':
    main()
