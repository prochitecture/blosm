"""Compare optimized primitives with dev78 on captured, real geometry inputs.

Run with Blender's standalone Python (no bpy needed). Does not generate street
meshes or change the capture. Timings are paired microbenchmarks, not Import.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import statistics
import time

import numpy as np
import _python_street_package
from _python_surface_performance_reference import ReferencePatternLibrary, visibility_center as original_center
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.quad_charts import visibility_center


def kernels(function, cases):
    result = []
    for polygon, epsilon in cases:
        try:
            value = function(polygon, epsilon)
        except ValueError as error:
            value = str(error)
        result.append(value)
    return result


def patterns(factory, counts):
    cache = factory()
    return [cache.quad(*item) for item in counts]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--capture', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    events = json.loads((args.capture / 'construction.json').read_text())['records']
    with np.load(args.capture / 'kernel_inputs.npz', allow_pickle=False) as data:
        cases = [(data['p' + str(row['replay'])], row['epsilon']) for row in events if 'replay' in row]
    counts = [row['counts'] for row in events if row['operation'] == 'PatternLibrary.quad']
    expected, times, reports = {}, {}, {}
    for name, runner, baseline, candidate, inputs in (
            ('kernel', kernels, original_center, visibility_center, cases),
            ('pattern', patterns, ReferencePatternLibrary, PatternLibrary, counts)):
        times[name] = {'reference': [], 'candidate': []}
        for label, function in (('reference', baseline), ('candidate', candidate),
                                ('candidate', candidate), ('reference', baseline)):
            start = time.perf_counter()
            result = runner(function, inputs)
            times[name][label].append(time.perf_counter() - start)
            if name not in expected:
                expected[name] = result
            else:
                for i, (a, b) in enumerate(zip(expected[name], result)):
                    if name == 'pattern':
                        assert np.array_equal(a.uv, b.uv) and np.array_equal(a.quads, b.quads), (name, i)
                    elif isinstance(a, str) or isinstance(b, str):
                        assert a == b, (name, i, a, b)
                    else:
                        assert np.array_equal(a, b), (name, i, a, b)
        medians = {key: statistics.median(value) for key, value in times[name].items()}
        reports[name] = dict(cases=len(inputs), exact_equal=True, seconds=times[name], median=medians,
                             speedup=medians['reference'] / medians['candidate'])
        print(name, json.dumps(reports[name]), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(passed=True, comparisons=reports,
        qualification='Paired ABBA primitive replay; not full-import performance'), indent=2), encoding='utf8')


if __name__ == '__main__':
    # Match the production benchmark's scoped Windows preference without
    # importing bpy or changing process priority/global power settings.
    path = Path(__file__).resolve().parents[1] / 'app/import_performance.py'
    spec = importlib.util.spec_from_file_location('street_replay_performance', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with module.ImportPerformance():
        main()
