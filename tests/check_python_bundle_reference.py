"""Compare Python meshes to the independent GN-era numerical geometry oracle.

The oracle directory is an explicit qualification input, not an addon runtime
dependency. It must supply bundle_reference.py and robust_bend_validation.py.
Run with ordinary Python/NumPy, or the bounded worker's --python option.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets import VERSION


def combined(parts):
    vertices, faces, kinds, lookup = [], [], [], {}
    for part in parts:
        mesh = part.mesh
        ids = []
        for point in mesh.xyz:
            # Oracle-only weld within ONE declared Bundle. Production keeps
            # separately owned pavement and grass objects with identical rims.
            key = tuple(np.round(point, 7))
            if key not in lookup:
                lookup[key] = len(vertices)
                vertices.append(point)
            ids.append(lookup[key])
        local = np.asarray(ids)[mesh.quads]
        faces.extend(local)
        kinds.extend([int(mesh.material == 'grass')] * len(local))
    return np.asarray(vertices), np.asarray(faces), np.asarray(kinds)


def sampled_boundary_distance(xyz, faces, reference):
    edges, counts = np.unique(np.sort(np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1),
                              return_counts=True, axis=0)
    output = xyz[edges[counts == 1], :2]
    cells = reference['corridor_cells']
    target = np.asarray([edge for cell in cells for edge in (cell[[0, 1]], cell[[3, 2]])]
                        + [cells[0][[0, 3]], cells[-1][[1, 2]]])

    def distance(queries, segments):
        maximum = 0.
        start, delta = segments[:, 0], segments[:, 1] - segments[:, 0]
        length2 = np.sum(delta * delta, axis=1)
        for at in range(0, len(queries), 128):
            difference = queries[at:at + 128, None, :] - start
            parameter = np.clip(np.sum(difference * delta, axis=-1) / length2, 0., 1.)
            distances = np.linalg.norm(difference - parameter[:, :, None] * delta, axis=-1)
            maximum = max(maximum, float(np.max(distances.min(axis=1))))
        return maximum

    a = np.concatenate((output[:, 0], output[:, 1], output.mean(axis=1)))
    b = np.concatenate((target[:, 0], target[:, 1], target.mean(axis=1)))
    forward, backward = distance(a, target), distance(b, output)
    return dict(maximum=max(forward, backward), output_to_reference=forward, reference_to_output=backward,
                sampling='boundary vertices and edge midpoints; independent reference at 128 subdivisions')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source = args.reference_root / 'bundle_reference.py'
    sys.path.insert(0, str(args.reference_root.resolve()))
    spec = importlib.util.spec_from_file_location('independent_bundle_reference', source)
    oracle = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(oracle)
    geometry = Path(__file__).resolve().parents[1] / 'way/street_network/python_backend'
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in geometry.glob('*.py')}
    report = dict(evaluation='pure numerical double precision', implementation=VERSION, geometry_sha256=hashes,
        cases=[], passed=False,
        oracle=str(source.resolve()), oracle_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    for fixture in oracle.fixture_cases():
        started = time.perf_counter()
        sections = (model.Section(1, 1, fixture['points_a'], width=fixture['width_a'], rise=fixture['rise'], interpolation=fixture['interpolation']),
                    model.Section(2, 2, fixture['points_b'], width=fixture['width_b'], rise=fixture['rise'], interpolation=fixture['interpolation']))
        bundle = model.Bundle(7, (model.BundleMember((model.Traversal(1),)),
                                 model.BundleMember((model.Traversal(2, not fixture['reverse_b']),))))
        network = model.Network(fixture['name'], sections, bundles=(bundle,))
        item = dict(name=fixture['name'], passed=False)
        try:
            result = PythonMeshBackend(network.settings).prepare(network)
        except ValueError as exc:
            item.update(passed=not fixture['supported'] and 'centerlines cross or touch' in str(exc), reason=str(exc))
            report['cases'].append(item)
            print(json.dumps(item), flush=True)
            continue
        parts = [p for p in result.products.values() if hasattr(p, 'mesh')]
        xyz, faces, kinds = combined(parts)
        topology = oracle.audit_mesh(xyz, faces, kinds, rise=fixture['rise'],
            expected_islands=fixture['expected_islands'], check_overlap=True)
        polygons = [xyz[f, :2] for f in faces]
        pavement = [p for p, kind in zip(polygons, kinds) if kind == 0]
        reference = oracle.smooth_pair_reference(fixture) if fixture['interpolation'] else oracle.exact_pair_reference(fixture)
        # Smooth boundaries are approximated to a physical 5 mm target. Their
        # coverage area budget is perimeter * tolerance, not exact equality to
        # a separately sampled polygon. Exact offsets retain the strict budget.
        perimeter = sum(np.sum(np.linalg.norm(np.diff(reference[m][side], axis=0), axis=1))
                        for m in ('member_a', 'member_b') for side in ('left', 'right'))
        area_budget = max(1e-5, float(perimeter) * network.settings.curve_tolerance) if fixture['interpolation'] else 1e-5
        coverage = oracle.coverage_audit(polygons, reference['corridor_cells'], area_tolerance=area_budget)
        preservation = oracle.source_coverage(pavement, reference, area_tolerance=area_budget)
        boundary = sampled_boundary_distance(xyz, faces, reference)
        boundary['passed'] = boundary['maximum'] <= (network.settings.curve_tolerance if fixture['interpolation'] else 1e-7)
        item.update(passed=fixture['supported'] and topology['passed'] and coverage['passed'] and preservation['passed'] and boundary['passed'],
            topology=topology, corridor=coverage, source_coverage=preservation,
            sampled_boundary_distance=boundary,
            quad_count=len(faces), seconds=time.perf_counter() - started,
            maximum_edge_ratio=max(p.metadata['quality']['maximum_edge_ratio'] for p in parts))
        report['cases'].append(item)
        print(json.dumps(dict(name=item['name'], passed=item['passed'], quads=len(faces),
            islands=topology['separator_components'], corridor_uncovered=coverage['uncovered_area'],
            corridor_outside=coverage['outside_area'], seconds=item['seconds'])), flush=True)
    report['passed'] = all(c['passed'] for c in report['cases'])
    (args.output / 'bundle_reference_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(dict(passed=report['passed'], cases=len(report['cases']))), flush=True)
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
