"""Native saved-network overlap replay against the frozen dev78 validator."""
import argparse
import importlib.util
import json
from pathlib import Path
import statistics
import sys
import time

import bpy
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.model import Key
from tested_python_streets.validation import validate_overlaps
from _python_validation_performance_reference import validate_overlaps as original_validate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()), load_ui=False, use_scripts=False)
    networks = []
    for controller in bpy.data.objects:
        if not controller.get('st_py_source_data'):
            continue
        meshes = {}
        for obj in controller['st_py_outputs'].objects:
            if obj.type != 'MESH':
                continue
            xyz = np.empty((len(obj.data.vertices), 3), np.float32)
            ids = np.empty(len(obj.data.loops), np.int32)
            obj.data.vertices.foreach_get('co', xyz.ravel())
            obj.data.loops.foreach_get('vertex_index', ids)
            meshes[Key.parse(obj['st_py_key'])] = Mesh(xyz, ids.reshape(-1, 4))
        networks.append(meshes)
    times = dict(reference=[], candidate=[])
    expected = None
    for label, function in (('reference', original_validate), ('candidate', validate_overlaps),
                            ('candidate', validate_overlaps), ('reference', original_validate)):
        start = time.perf_counter()
        result = [function(meshes) for meshes in networks]
        times[label].append(time.perf_counter() - start)
        if expected is None:
            expected = result
        else:
            assert result == expected, (result, expected)
    report = dict(passed=True, exact_validation_reports=True, reports=expected, seconds=times,
        speedup=statistics.median(times['reference']) / statistics.median(times['candidate']),
        qualification='Paired ABBA validation replay on saved coordinates; not full Import')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    spec = importlib.util.spec_from_file_location('street_validation_performance', ROOT / 'app/import_performance.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with module.ImportPerformance():
        main()
