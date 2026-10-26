"""Compare NumPy fixed-chart projection with an external test-only oracle."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import _python_street_package
from tested_python_streets import VERSION
from tested_python_streets.mesh import Mesh
from tested_python_streets.slope_projection import SlopeProjection


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', type=Path, required=True, help='Directory containing oracle_checks.json and oracle_mesh.npz')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    oracle = json.loads((args.reference / 'oracle_checks.json').read_text(encoding='utf8'))
    path = args.reference / 'oracle_mesh.npz'
    with np.load(path, allow_pickle=False) as data:
        xyz, faces = data['xyz'].copy(), data['quads'].copy()
        expected = xyz[:, 2].copy()
        initial = data['initial_heights'].copy()
        fixed = data['boundary'].copy()
    xyz[:, 2] = initial
    source = Mesh(xyz.copy(), faces)
    start = time.perf_counter()
    result, projection = SlopeProjection().solve(source, fixed,
        datum=oracle['height_range'][0], ceiling=oracle['height_range'][1], maximum_slope=oracle['slope_limit'])
    seconds = time.perf_counter() - start
    slopes = result.quad_slopes()
    edges = np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2)
    a, b = np.r_[edges[:, 0], edges[:, 1]], np.r_[edges[:, 1], edges[:, 0]]
    minimum = np.full(len(xyz), np.inf)
    np.minimum.at(minimum, a, result.xyz[b, 2])
    minima = np.flatnonzero((~fixed) & (result.xyz[:, 2] < minimum - 1e-6))
    difference = float(np.max(abs(result.xyz[:, 2] - expected)))
    assert np.array_equal(source.xyz, xyz)
    assert np.array_equal(result.quads, faces)
    assert np.array_equal(result.xyz[:, :2], xyz[:, :2])
    assert np.array_equal(result.xyz[fixed, 2], initial[fixed])
    assert slopes.max() <= oracle['slope_limit'] + 1e-8
    assert result.xyz[:, 2].min() >= oracle['height_range'][0]
    assert result.xyz[:, 2].max() <= oracle['height_range'][1]
    assert difference < 1e-5
    assert not len(minima)
    report = dict(passed=True, implementation=VERSION, seconds=seconds,
        reference_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), projection=projection,
        maximum_physical_slope=float(slopes.max()), oracle_maximum_height_difference=difference,
        isolated_interior_minima=len(minima), quads=len(faces), vertices=len(xyz))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'projection_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    np.savez(args.output / 'projected_mesh.npz', xyz=result.xyz, quads=faces)
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
