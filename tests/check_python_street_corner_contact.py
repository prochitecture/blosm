"""Independent bounded contact oracle and original outside-support experiment.

SciPy and Matplotlib are diagnostic dependencies only, never addon imports.
Expanded cases use explicit original OSM segments and are not full mesh builds.
"""
import argparse
from dataclasses import replace
import json
import math
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _python_street_package
from tested_python_streets.model import Network, Key
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import UpdateCoordinator, BuildContext
from tested_python_streets.routes import RouteTopology, Route
from tested_python_streets.intersections import Approach, _corner


def approaches(network, site):
    patterns = PatternLibrary()
    update = UpdateCoordinator((PathProvider(), ProfileProvider()), patterns).prepare(network)
    context = BuildContext(network, update.products, patterns)
    topology = RouteTopology(network)
    result = []
    for end in site.ends:
        route = Route.build(context, topology.trace(end))
        result.append(Approach(route, end.end, context[Key('profile', end.section_id)]))
    return sorted(result, key=lambda a: math.atan2(a.direction[1], a.direction[0]) % math.tau)


def expand(network, pair, context, hops):
    """Preserve cut controls; add only matching original outside segment data."""
    if context.get('version') != 1 or context.get('issues'):
        raise ValueError('Outside context must have unambiguous original provenance')
    sections = network.section_index.copy()
    changes = []
    for approach in pair:
        last = approach.path.pieces[-1]
        identity, end = last.section.id, 1 if last.forward else 0
        matches = [(node, opened) for node in context['nodes'] for opened in node.get('open_ends', ())
                   if opened['section_id'] == identity and opened['end'] == end and opened['owner_scope'] == 'ACTIVE']
        if len(matches) != 1:
            raise ValueError('Approach does not have one explicit original rectangle cut')
        node, opened = matches[0]
        section = sections[identity]
        cut = np.asarray(section.points[0 if end == 0 else -1])[:2]
        incoming = [i for i in node['incidences'] if i['source_way_id'] == opened['source_way_id']
                    and i['source_nodes'] == opened['source_nodes'] and abs(i['parameter'] - opened['parameter']) < 1e-8]
        if len(incoming) != 1 or np.linalg.norm(cut - incoming[0]['clipped_point']) > .001:
            raise ValueError('Cut provenance does not match the authored endpoint')
        extra = [tuple(node['position']) + (0.,)]
        segment_ids = [opened['source_nodes']]
        if hops > 1:
            legs = [leg for leg in node['legs'] if leg['source_way_id'] == opened['source_way_id']
                    and leg['source_nodes'] != opened['source_nodes'] and leg['category'] == section.category]
            if node['topology_kind'] != 'CONTINUATION' or len(legs) != 1:
                raise ValueError('Further support lacks a unique same-way continuation')
            extra.append(tuple(legs[0]['other_position']) + (0.,))
            segment_ids.append(legs[0]['source_nodes'])
        points = tuple(reversed(extra)) + section.points if end == 0 else section.points + tuple(extra)
        sections[identity] = replace(section, points=points)
        changes.append(dict(section=identity, end=end, source_way=opened['source_way_id'],
                            source_segments=segment_ids, extra_controls=extra))
    return replace(network, sections=tuple(sections[s.id] for s in network.sections)), changes


def oracle(first, second, radius):
    gap = math.atan2(np.cross(first.direction, second.direction), first.direction @ second.direction) % math.tau
    side = 1. if gap < math.pi else -1.
    distances = (first.offset - first.width * .5 - side * radius,
                 second.offset + second.width * .5 + side * radius)
    lower = -4. * max(first.width + radius, second.width + radius)
    upper = np.array((first.length, second.length))
    stations = [np.linspace(lower, a.length, 257) for a in (first, second)]
    curves = [a.point(s, d) for a, s, d in zip((first, second), stations, distances)]
    squared = np.sum((curves[0][:, None] - curves[1][None]) ** 2, axis=-1)
    closest = np.unravel_index(np.argpartition(squared.ravel(), 16)[:16], squared.shape)
    seeds = list(zip(stations[0][closest[0]], stations[1][closest[1]]))
    seeds += [(a, b) for a in np.linspace(lower, upper[0], 7) for b in np.linspace(lower, upper[1], 7)]
    def residual(x):
        return first.point(x[0], distances[0]) - second.point(x[1], distances[1])
    best, best_x, roots = float('inf'), None, []
    for seed in seeds:
        fit = least_squares(residual, np.clip(seed, lower + 1e-8, upper - 1e-8),
                            bounds=(lower, upper), jac='3-point', max_nfev=150,
                            ftol=1e-12, xtol=1e-12, gtol=1e-12)
        error = float(np.linalg.norm(fit.fun))
        if error < best:
            best, best_x = error, fit.x.copy()
        if error <= 1e-7 and all(np.linalg.norm(fit.x - r['stations']) > 1e-4 for r in roots):
            center = first.point(fit.x[0], distances[0])
            contacts = np.array((first.point(fit.x[0], first.offset - first.width * .5),
                                 second.point(fit.x[1], second.offset + second.width * .5)))
            vectors = contacts - center
            sweep = math.atan2(np.cross(vectors[0], vectors[1]), vectors[0] @ vectors[1])
            roots.append(dict(stations=fit.x.tolist(), residual=error, center=center.tolist(),
                radius_error=float(np.max(abs(np.linalg.norm(vectors, axis=1) - radius))),
                start_angle=math.atan2(vectors[0, 1], vectors[0, 0]), sweep=sweep,
                correct_sector=bool(sweep * side <= 1e-7 and abs(sweep) <= math.pi + 1e-7)))
    try:
        actual = _corner(first, second, radius)
        production = dict(stations=actual.stations, radius=actual.radius)
    except ValueError as exc:
        production = dict(error=str(exc))
    return dict(minimum_found_residual=best, best_stations=best_x.tolist(), roots=roots,
                production=production), curves


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--source-report', type=Path, required=True)
    parser.add_argument('--junction', type=int, required=True)
    parser.add_argument('--first', type=int, required=True)
    parser.add_argument('--second', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--trial-radius', type=float, action='append', default=[])
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    network = next(n for n in map(Network.from_dict, json.loads(args.input.read_text()))
                   if any(j.id == args.junction for j in n.intersections))
    site = next(j for j in network.intersections if j.id == args.junction)
    original = approaches(network, site)
    first = next(a for a in original if a.source.id == args.first)
    second = next(a for a in original if a.source.id == args.second)
    assert original[(original.index(first) + 1) % len(original)] is second
    radius = dict(zip((e.section_id for e in site.ends), site.radii)).get(args.first, site.radius)
    context = json.loads(args.source_report.read_text())['source']['clipped_source_context']
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none', 'font.family': 'DejaVu Sans'})
    fig, axes = plt.subplots(1, 3, figsize=(15, 7), layout='constrained')
    records = []
    for hops, ax in enumerate(axes):
        try:
            candidate, additions = (network, []) if hops == 0 else expand(network, (first, second), context, hops)
        except ValueError as exc:
            records.append(dict(extra_source_hops=hops, support_unavailable=str(exc)))
            ax.text(.5, .5, str(exc), ha='center', va='center', wrap=True, transform=ax.transAxes)
            ax.set_axis_off()
            continue
        current = approaches(candidate, site)
        pair = [next(a for a in current if a.source.id == identity) for identity in (args.first, args.second)]
        result, curves = oracle(*pair, radius)
        records.append(dict(extra_source_hops=hops, lengths=[a.length for a in pair], additions=additions, **result))
        print(json.dumps(dict(hops=hops, lengths=records[-1]['lengths'],
                              residual=result['minimum_found_residual'], roots=len(result['roots']))), flush=True)
        for a, color in zip(pair, ('#c26324', '#168699')):
            s = np.linspace(0., a.length, 500)
            xy = a.point(s, a.offset)
            ax.plot(*xy.T, color=color, label=f'Section {a.source.id} approach')
            for sign in (-1., 1.):
                curb = a.point(s, a.offset + sign * a.width * .5)
                ax.plot(*curb.T, color=color, alpha=.35, lw=.8)
        for curve in curves:
            ax.plot(*curve.T, ':', color='#557085', lw=.9)
        for root in result['roots']:
            angle = root['start_angle'] + np.linspace(0., root['sweep'], 100)
            xy = np.array(root['center']) + radius * np.column_stack((np.cos(angle), np.sin(angle)))
            ax.plot(*xy.T, color='#1e8f57', lw=2)
        y = context['bounds_m'][2]
        ax.axhline(y, color='#757575', ls='--', lw=1, label='Import rectangle boundary')
        ax.set(aspect='equal', xlabel='Metres', ylabel='Metres', title=('Retained source', 'First original outside node',
                                                                    'Next original same-way segment')[hops])
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(alpha=.15)
    fig.suptitle(f'Junction {args.junction}: independently searching the requested {radius:g} m corner', fontsize=16)
    fig.savefig(args.output / 'contact_support.svg')
    fig.savefig(args.output / 'contact_support.png', dpi=150)
    report = dict(junction=args.junction, sections=[args.first, args.second], radius=radius,
        scope='Bounded multistart numerical oracle, not a nonexistence proof or a complete mesh qualification', cases=records)
    report['trial_radii'] = []
    for trial in args.trial_radius:
        if not math.isfinite(trial) or trial <= 0.:
            raise ValueError('Trial radii must be positive and finite')
        result, _ = oracle(first, second, trial)
        report['trial_radii'].append(dict(radius=trial, **result))
        print(json.dumps(dict(trial_radius=trial, residual=result['minimum_found_residual'],
                              roots=len(result['roots']))), flush=True)
    (args.output / 'contact_checks.json').write_text(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
