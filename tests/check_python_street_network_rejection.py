"""Retain rejected numerical buffers for an independent overlap diagnosis.

This diagnostic never publishes or commits a candidate. The normal overlap
validator still decides whether it passes; saved meshes may be invalid.
"""
import argparse
import json
from pathlib import Path
import time
from unittest.mock import patch
import numpy as np
import _python_street_package
from check_python_street_targeted_quality import TargetedChart
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.providers import Part
from tested_python_streets.validation import validate_overlaps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--junction', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    network, = [n for n in map(model.Network.from_dict, json.loads(args.input.read_text()))
                if any(j.id == args.junction for j in n.intersections)]
    start = time.perf_counter()
    backend = PythonMeshBackend(network.settings)
    with patch('tested_python_streets.compound_surfaces.CompoundQuadChart', TargetedChart):
        update = backend.updates.prepare(network)
    parts = {key: value for key, value in update.products.items()
             if isinstance(value, Part) and value.mesh is not None}
    for key, part in parts.items():
        np.savez_compressed(args.output / (key.text().replace(':', '_') + '.npz'),
                            xyz=part.mesh.xyz, quads=part.mesh.quads)
    report = dict(published=False, committed=False, input=str(args.input.resolve()),
                  geometry_seconds=time.perf_counter() - start,
                  parts={key.text(): dict(quality=p.mesh.validate(), metadata=p.metadata)
                         for key, p in parts.items()}, timings=update.timings)
    try:
        report['overlaps'] = validate_overlaps({key: p.mesh for key, p in parts.items()},
            changed=update.rebuilt, cell_size=max(1., network.settings.spacing * 2.))
        report['passed'] = True
    except ValueError as exc:
        report.update(passed=False, reason=str(exc))
    report['total_seconds'] = time.perf_counter() - start
    (args.output / 'rejection_checks.json').write_text(json.dumps(report, indent=2))
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('parts', 'timings')}), flush=True)


if __name__ == '__main__':
    main()
