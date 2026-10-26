"""Compare bounded existing shape/quad-rotation repairs on captured chart blocks."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.quad_charts import FixedRimQuality
from tested_python_streets.quad_rotations import QuadEdgeRotations
from tested_python_streets.validation import validate_overlaps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--capture', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    for path in sorted(args.capture.glob('block_*.npz')):
        data = np.load(path)
        boundary = data['boundary']
        mesh = Mesh(data['xyz'], data['quads'], seams={0: np.arange(len(boundary))})
        record = dict(block=path.stem, initial=mesh.validate(), trials=[])
        best = mesh
        for name, solve in (
            ('weighted-energy-300', lambda source: FixedRimQuality(iterations=300, power=2).improve(source)),
            ('weighted-energy-900', lambda source: FixedRimQuality(iterations=900, power=2).improve(source)),
            ('quad-rotations-weighted', lambda source: QuadEdgeRotations(max_passes=8).improve(source,
                FixedRimQuality(iterations=300, power=2)))):
            start = time.perf_counter()
            result, details = solve(best)
            quality = result.validate()
            assert np.array_equal(result.xyz[result.seams[0], :2], boundary)
            if quality['maximum_edge_ratio'] < best.validate()['maximum_edge_ratio']:
                best = result
            record['trials'].append(dict(name=name, seconds=time.perf_counter() - start, quality=quality,
                                        details=details, exact_fixed_rim=True))
            print(json.dumps(dict(block=path.stem, trial=name, quality=quality)), flush=True)
            if best.validate()['maximum_edge_ratio'] <= 8.:
                break
        record['final'] = best.validate()
        record['overlaps'] = validate_overlaps({path.stem: best}, cell_size=1.4)
        record['planar_qualified'] = record['final']['maximum_edge_ratio'] <= 8.
        np.savez_compressed(args.output / (path.stem + '_candidate.npz'), xyz=best.xyz, quads=best.quads, boundary=boundary)
        records.append(record)
    report = dict(planar_only=True, pavement_surface_qualified=False, blocks=records,
                  all_blocks_within_ratio_8=bool(records) and all(r['planar_qualified'] for r in records))
    (args.output / 'repair_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(dict(planar_only=True, all_blocks_within_ratio_8=report['all_blocks_within_ratio_8'])), flush=True)


if __name__ == '__main__':
    main()
