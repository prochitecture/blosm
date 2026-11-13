"""Native BVH equivalence and timings for deduplicated terrain ray queries.

Read a benchmark scene in isolated factory Blender; do not regenerate streets.
Compare every existing quad-safety sample and output vertex with original rays.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import time

import bpy
import numpy as np
from mathutils import Vector

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
import _python_street_package
from tested_python_streets.terrain.blender import TerrainSampler
from tested_python_streets.terrain.alignment import quad_samples


def reference(sampler, xy):
    result = np.empty(len(xy))
    for i, (x, y) in enumerate(xy):
        hit = sampler.bvh.ray_cast(Vector((float(x), float(y), sampler.top)),
                                  Vector((0, 0, -1)), sampler.reach)[0]
        if hit is None:
            raise ValueError(f'Original terrain has no surface below corridor point ({x:.3f}, {y:.3f}); extend the terrain or revise the street footprint')
        result[i] = hit.z
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()), load_ui=False, use_scripts=False)
    controllers = [o for o in bpy.data.objects if o.get('st_py_source_data')]
    scale = bpy.context.scene.unit_settings.scale_length
    seconds, queries, rays, batches = dict(reference=0., candidate=0.), 0, 0, 0
    for controller in controllers:
        sampler = TerrainSampler(bpy.data.objects['Terrain'], controller, scale)
        assert sampler._bvh is None
        assert sampler(np.empty((0, 2))).shape == (0,)
        assert sampler._bvh is None, 'An empty query should not build the terrain index'

        def check(xy):
            nonlocal queries, batches
            results = {}
            # Alternate the order to avoid consistently favoring warm BVH data.
            order = ('reference', 'candidate') if batches % 2 else ('candidate', 'reference')
            for label in order:
                start = time.perf_counter()
                results[label] = sampler(xy) if label == 'candidate' else reference(sampler, xy)
                seconds[label] += time.perf_counter() - start
            assert np.array_equal(results['reference'], results['candidate']), batches
            queries += len(xy)
            batches += 1

        for obj in controller['st_py_outputs'].objects:
            if obj.type != 'MESH':
                continue
            xyz = np.empty((len(obj.data.vertices), 3), np.float32)
            obj.data.vertices.foreach_get('co', xyz.ravel())
            xyz = xyz.astype(float) * scale
            ids = np.empty(len(obj.data.loops), np.int32)
            obj.data.loops.foreach_get('vertex_index', ids)
            quads = ids.reshape(-1, 4)
            for _, points in quad_samples(xyz, quads, 1.):
                check(points[:, :, :2].reshape(-1, 2))
            check(xyz[:, :2])
        rays += sampler.ray_count
        # Empty/small calls, identical float32 origins from distinct doubles,
        # cache-free batch boundaries, and first-miss diagnostics.
        points = np.tile(xyz[:min(100, len(xyz)), :2], (300, 1))
        for count in (0, 1, 31, 32, 16385):
            data = points[:count].copy()
            if count:
                data[::2] = np.nextafter(data[::2], np.inf)
            assert np.array_equal(reference(sampler, data), sampler(data))
        outside = np.vstack((points[:5], np.tile(((1e7, 1e7), (-1e7, -1e7)), (32, 1))))
        messages = []
        for call in (lambda: reference(sampler, outside), lambda: sampler(outside)):
            try:
                call()
            except ValueError as error:
                messages.append(str(error))
        assert len(messages) == 2 and messages[0] == messages[1], messages
    report = dict(passed=True, exact_heights=True, requested_queries=queries,
                  unique_bvh_rays=rays, batches=batches, seconds=seconds,
                  lazy_index=True, empty_query_without_index=True,
                  speedup=seconds['reference'] / seconds['candidate'],
                  qualification='Same-process terrain query comparison; not full Import')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding='utf8')
    print('TERRAIN_QUERY_EQUIVALENCE', json.dumps(report), flush=True)


if __name__ == '__main__':
    spec = importlib.util.spec_from_file_location('terrain_query_performance', ROOT / 'app/import_performance.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with module.ImportPerformance():
        main()
