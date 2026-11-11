"""Compare benchmark scenes independently of the timed generation process.

Launch in isolated factory Blender with autoexec disabled. Read only the input
scenes; put the report in artifacts. UUIDs and terrain-reference metadata are
excluded, but published float32 positions, prepared relief, face indices, UVs,
material names and construction metadata are checked. No generator is executed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import bpy
import numpy as np


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def read_scene(path):
    bpy.ops.wm.open_mainfile(filepath=str(path.resolve()), load_ui=False, use_scripts=False)
    sources, meshes = {}, {}
    for obj in bpy.data.objects:
        if obj.get('st_py_source_data'):
            source = json.loads(obj['st_py_source_data'])
            source.pop('uuid', None)
            identity = str(min(s['id'] for s in source['sections']))
            assert identity not in sources
            sources[identity] = digest(source)
        key = obj.get('st_py_key')
        if not key or obj.type != 'MESH':
            continue
        assert key not in meshes, key
        xyz = np.empty((len(obj.data.vertices), 3), np.float32)
        obj.data.vertices.foreach_get('co', xyz.ravel())
        final_positions = hashlib.sha256(xyz.tobytes()).hexdigest()
        relief = obj.data.attributes.get('st_terrain_prepared_z')
        if relief:
            z = np.empty(len(xyz), np.float32)
            relief.data.foreach_get('value', z)
            xyz[:, 2] = z
        ids = np.empty(len(obj.data.loops), np.int32)
        obj.data.loops.foreach_get('vertex_index', ids)
        sizes = np.empty(len(obj.data.polygons), np.int32)
        obj.data.polygons.foreach_get('loop_total', sizes)
        assert np.all(sizes == 4), key
        uv = {}
        for layer in obj.data.uv_layers:
            values = np.empty((len(obj.data.loops), 2), np.float32)
            layer.uv.foreach_get('vector', values.ravel())
            uv[layer.name] = hashlib.sha256(values.tobytes()).hexdigest()
        metadata = json.loads(obj.get('st_py_metadata', '{}'))
        for name in ('terrain_datum', 'terrain_reference', 'bounds'):
            metadata.pop(name, None)
        meshes[key] = dict(prepared_positions=hashlib.sha256(xyz.tobytes()).hexdigest(),
            final_positions=final_positions, faces=hashlib.sha256(ids.tobytes()).hexdigest(),
            quads=len(sizes), uv=uv, materials=[m.name if m else None for m in obj.data.materials],
            construction_metadata=digest(metadata))
    assert sources and meshes
    return dict(sources=sources, meshes=meshes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--allow-terrain-differences', action='store_true')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    if len(args.input) < 2:
        parser.error('At least two input scenes are required')
    evidence = {str(p.resolve()): read_scene(p) for p in args.input}
    baseline = evidence[str(args.input[0].resolve())]
    differences = {}
    for name, candidate in evidence.items():
        keys = sorted(baseline['meshes'].keys() | candidate['meshes'].keys())
        changed = []
        for key in keys:
            a, b = dict(baseline['meshes'].get(key, {})), dict(candidate['meshes'].get(key, {}))
            if args.allow_terrain_differences:
                a.pop('final_positions', None)
                b.pop('final_positions', None)
            if a != b:
                changed.append(key)
        differences[name] = dict(source_equal=baseline['sources'] == candidate['sources'], changed_meshes=changed)
    report = dict(passed=all(d['source_equal'] and not d['changed_meshes'] for d in differences.values()),
        baseline=str(args.input[0].resolve()), final_terrain_positions_checked=not args.allow_terrain_differences,
        differences=differences, evidence=evidence)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding='utf8')
    print('STREET_TIMING_EQUIVALENCE', json.dumps(dict(passed=report['passed'], differences=differences)), flush=True)
    assert report['passed'], differences


if __name__ == '__main__':
    main()
