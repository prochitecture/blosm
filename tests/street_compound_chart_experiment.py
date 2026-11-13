"""Unselected four-sided-first compound experiment using real source contracts.

The normal backend supplies each actual union and its prepared dependencies.
Alternative construction runs separately; its result is never published by the
reference backend. Times are diagnostic (warm dependencies, cold chart patterns),
not a whole-Import speed claim. Output comparisons can reject a candidate; they
do not replace the full numerical and rendered acceptance gates.
"""
import argparse
from dataclasses import replace
import copy
import importlib.util
import json
from pathlib import Path
import time
from unittest.mock import patch

import numpy as np
import _python_street_package
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.four_side_charts import FourSidedQuadChart
from tested_python_streets import partition
from tested_python_streets.model import Network
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.surface_query import QuadHeightField


Layered = partition.LayeredQuadChart


class FourSidedFirst(Layered):
    def build(self, boundary, *, require_quality=True):
        if require_quality:
            try:
                return FourSidedQuadChart(self.patterns, max_cells=self.max_cells,
                    quality_iterations=min(100, self.quality.iterations),
                    maximum_ratio=self.maximum_ratio).build(boundary)
            except ValueError:
                pass
        return super().build(boundary, require_quality=require_quality)


def compare(reference, candidate):
    first, second = reference.mesh, candidate.mesh
    exact = (np.array_equal(first.xyz, second.xyz)
             and np.array_equal(first.quads, second.quads))
    result = dict(exact_mesh=exact, reference_quads=len(first.quads),
                  candidate_quads=len(second.quads))
    if exact:
        result['maximum_height_difference'] = 0.
        return result
    # Symmetric samples include both meshes' cell centres and vertices. The
    # bound is diagnostic; a different chart still needs independent visual QA.
    samples = []
    for mesh in (first, second):
        for points in (mesh.xyz[:, :2], mesh.xyz[mesh.quads, :2].mean(axis=1)):
            samples.append(points[np.linspace(0, len(points) - 1,
                              min(1024, len(points)), dtype=int)])
    samples = np.unique(np.concatenate(samples), axis=0)
    fields = [QuadHeightField(surface.mesh, datum=surface.field.datum,
                ceiling=surface.field.ceiling) for surface in (reference, candidate)]
    try:
        difference = np.abs(fields[0](samples) - fields[1](samples))
        result.update(samples=len(samples), maximum_height_difference=float(difference.max()),
                      p95_height_difference=float(np.quantile(difference, .95)))
    except ValueError as error:
        result['query_error'] = str(error)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--minimum-seconds', type=float, default=.25)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    original = CompoundSurfaceBuilder.build

    def build(owner, group, context):
        start = time.perf_counter()
        reference = original(owner, group, context)
        elapsed = time.perf_counter() - start
        if elapsed < args.minimum_seconds:
            return reference
        record = dict(index=len(records), members=[key.text() for key in group.members],
                      reference_seconds=elapsed, reference_quads=len(reference.mesh.quads))
        trial = copy.copy(owner)
        trial.patterns = PatternLibrary()
        candidate_context = replace(context, patterns=trial.patterns)
        start = time.perf_counter()
        try:
            with patch.object(partition, 'LayeredQuadChart', FourSidedFirst):
                candidate = original(trial, group, candidate_context)
            record.update(candidate_seconds=time.perf_counter() - start, accepted=True)
            record['comparison'] = compare(reference, candidate)
            np.savez_compressed(args.output / f'case_{len(records)}.npz',
                reference_xyz=reference.mesh.xyz, reference_quads=reference.mesh.quads,
                candidate_xyz=candidate.mesh.xyz, candidate_quads=candidate.mesh.quads)
            record['candidate_metadata'] = candidate.metadata
        except ValueError as error:
            record.update(candidate_seconds=time.perf_counter() - start,
                          accepted=False, error=str(error))
        records.append(record)
        print(json.dumps({k:v for k,v in record.items() if k != 'candidate_metadata'}), flush=True)
        # A failure in the experiment must not replace accepted reference output.
        return reference

    with patch.object(CompoundSurfaceBuilder, 'build', build):
        for source in args.source:
            network = Network.from_dict(json.loads(source.read_text(encoding='utf8')))
            prepared = PythonMeshBackend(network.settings).prepare(network)
            print('REFERENCE_COMPLETE', source.name, len(prepared.products), flush=True)
    report = dict(cases=records, sources=[str(p.resolve()) for p in args.source],
        selected_in_runtime=False, qualification='Unselected construction/shape experiment')
    (args.output / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf8')


if __name__ == '__main__':
    path = Path(__file__).resolve().parents[1] / 'app/import_performance.py'
    spec = importlib.util.spec_from_file_location('compound_experiment_performance', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with module.ImportPerformance():
        main()
