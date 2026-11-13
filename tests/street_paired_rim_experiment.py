"""Unselected separator construction experiment; not imported by the addon.

Pair visible opposite rim vertices before seeking local visibility kernels.
The production generator remains unchanged until geometry and visual gates pass.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import time

import numpy as np
import _python_street_package
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.separator_surfaces import (
    FixedRimSeparator, PartitionedSeparatorChart, SeparatorStripCuts)
from tested_python_streets.profiles import plateau_extension


class PairedRimChart(PartitionedSeparatorChart):
    four_sided_first = True
    def cut_policy(self):
        return SeparatorStripCuts(self.spacing, initial_visibility=False)

    def block_builders(self, remaining):
        builders = super().block_builders(remaining)
        return tuple(reversed(builders)) if self.four_sided_first else builders


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--capture', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--layered-first', action='store_true')
    args = parser.parse_args()
    PairedRimChart.four_sided_first = not args.layered_first
    args.output.mkdir(parents=True, exist_ok=True)
    rows = [r for r in json.loads((args.capture / 'construction.json').read_text())['records']
            if 'separator_replay' in r]
    results = []
    with np.load(args.capture / 'separator_inputs.npz', allow_pickle=False) as inputs:
        for row in rows:
            index = row['separator_replay']
            boundary = inputs['p' + str(index)]
            settings = row['settings']
            result = dict(case=index)
            saved = {}
            for name in ('reference', 'paired-rim'):
                start = time.perf_counter()
                try:
                    if name == 'reference':
                        surface = FixedRimSeparator(PatternLibrary(), spacing=settings['spacing'],
                            max_cells=settings['max_cells'], paired_rim=False).build(boundary, datum=settings['datum'],
                            height=settings['height'], exaggeration=settings['exaggeration'])
                    else:
                        surface = PairedRimChart(PatternLibrary(), spacing=settings['spacing'],
                            interior_spacing=2. * settings['spacing'], max_cells=settings['max_cells']).build(boundary[:, :2])
                        surface.mesh.xyz[:, 2], _ = plateau_extension(surface.mesh.xyz[:, :2], boundary,
                            settings['datum'], settings['height'])
                        surface.mesh.xyz[:len(boundary)] = boundary
                    seconds = time.perf_counter() - start
                    mesh = surface.mesh
                    np.testing.assert_array_equal(mesh.xyz[:len(boundary)], boundary)
                    result[name] = dict(seconds=seconds, accepted=True, quality=mesh.validate(),
                                        metadata=surface.metadata)
                    saved[name + '_xyz'], saved[name + '_quads'] = mesh.xyz, mesh.quads
                except ValueError as error:
                    result[name] = dict(seconds=time.perf_counter() - start, accepted=False, error=str(error))
            np.savez_compressed(args.output / f'case_{index}.npz', boundary=boundary, **saved)
            results.append(result)
            print(json.dumps({k: ({a:b for a,b in v.items() if a != 'metadata'} if isinstance(v,dict) else v)
                              for k,v in result.items()}), flush=True)
    (args.output / 'report.json').write_text(json.dumps(dict(cases=results,
        layered_first=args.layered_first,
        qualification='Construction experiment only; not selected in runtime'), indent=2), encoding='utf8')


if __name__ == '__main__':
    path = Path(__file__).resolve().parents[1] / 'app/import_performance.py'
    spec = importlib.util.spec_from_file_location('street_experiment_performance', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with module.ImportPerformance():
        main()
