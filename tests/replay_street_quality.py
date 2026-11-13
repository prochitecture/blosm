"""Compare complete fixed-rim descents using captured real chart seeds.

Only the gradient workspace is substituted. Boundary pins, energy, iteration
limits, line search and final mesh checks are the same on both sides.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import statistics
import time
from unittest.mock import patch
import numpy as np
import _python_street_package
from _python_surface_performance_reference import ReferenceQuadEnergy, ReferenceQuadGradient
from tested_python_streets import quad_charts
from tested_python_streets.mesh import Mesh
from tested_python_streets.quad_workspace import QuadEnergy, QuadGradient


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--capture', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = [r for r in json.loads((args.capture / 'construction.json').read_text())['records']
            if 'quality_replay' in r]
    results = []
    with np.load(args.capture / 'quality_inputs.npz', allow_pickle=False) as inputs:
        for row in rows:
            index = row['quality_replay']
            mesh = Mesh(inputs[f'q{index}_xyz'], inputs[f'q{index}_quads'])
            expected, expected_report = None, None
            timings = {'reference': [], 'candidate': []}
            for label in ('reference', 'candidate', 'candidate', 'reference'):
                gradient = ReferenceQuadGradient if label == 'reference' else QuadGradient
                energy = ReferenceQuadEnergy if label == 'reference' else QuadEnergy
                with patch.object(quad_charts, 'QuadGradient', gradient), \
                        patch.object(quad_charts, 'QuadEnergy', energy):
                    start = time.perf_counter()
                    actual, report = quad_charts.FixedRimQuality(**row['quality_settings']).improve(mesh)
                    timings[label].append(time.perf_counter() - start)
                if expected is None:
                    expected, expected_report = actual, report
                else:
                    np.testing.assert_array_equal(actual.xyz.view(np.uint64), expected.xyz.view(np.uint64))
                    np.testing.assert_array_equal(actual.quads, expected.quads)
                    assert report == expected_report, (index, report, expected_report)
            results.append(dict(case=index, quads=len(mesh.quads), settings=row['quality_settings'],
                steps=expected_report['steps'], seconds=timings, exact_coordinates_and_decisions=True))
    medians = {name: sum(statistics.median(r['seconds'][name]) for r in results)
               for name in ('reference', 'candidate')}
    report = dict(passed=True, cases=results, total_median_seconds=medians,
        speedup=medians['reference'] / medians['candidate'],
        qualification='Paired same-process complete chart repairs; not full Import')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps({key: value for key, value in report.items() if key != 'cases'}), flush=True)


if __name__ == '__main__':
    path = Path(__file__).resolve().parents[1] / 'app/import_performance.py'
    spec = importlib.util.spec_from_file_location('street_quality_performance', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with module.ImportPerformance():
        main()
