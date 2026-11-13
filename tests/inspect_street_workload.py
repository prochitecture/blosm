"""Read saved street output sizes and construction choices without regeneration.

Launch in an isolated factory Blender with autoexec disabled. This is a workload
inspection, not a timing benchmark. Reports belong in the artifacts directory.
"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

import bpy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--export-sources', type=Path,
                        help='Export normalized source contracts for numerical experiments')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()), load_ui=False, use_scripts=False)
    groups = defaultdict(lambda: dict(objects=0, vertices=0, quads=0, uv_bytes_float32=0))
    charts, meshes, sources = Counter(), [], []

    def choices(value):
        if isinstance(value, dict):
            if 'selected_chart' in value:
                charts[str(value['selected_chart'])] += 1
            for child in value.values():
                choices(child)
        elif isinstance(value, list):
            for child in value:
                choices(child)

    for obj in bpy.data.objects:
        if obj.get('st_py_source_data'):
            source = json.loads(obj['st_py_source_data'])
            if args.export_sources:
                args.export_sources.mkdir(parents=True, exist_ok=True)
                first = min(section['id'] for section in source['sections'])
                (args.export_sources / f'network_{first}.json').write_text(
                    json.dumps(source, indent=2), encoding='utf8')
            counts = {key: len(source.get(key, [])) for key in
                      ('sections', 'intersections', 'bundles', 'furniture')}
            counts['streets'] = len({section['street_id'] for section in source['sections']})
            sources.append(counts)
        key = obj.get('st_py_key')
        if not key or obj.type != 'MESH':
            continue
        assert all(p.loop_total == 4 for p in obj.data.polygons), key
        record = dict(key=key, vertices=len(obj.data.vertices), quads=len(obj.data.polygons),
                      uv_bytes_float32=len(obj.data.loops) * 2 * 4 * len(obj.data.uv_layers))
        meshes.append(record)
        group = groups[key.split(':', 1)[0]]
        group['objects'] += 1
        for field in ('vertices', 'quads', 'uv_bytes_float32'):
            group[field] += record[field]
        choices(json.loads(obj.get('st_py_metadata', '{}')))
    assert meshes, 'No saved Python street outputs found'
    report = dict(scene=str(args.scene.resolve()), sources=sources, groups=dict(groups),
                  selected_charts_in_saved_metadata=dict(charts),
                  largest_meshes=sorted(meshes, key=lambda m: m['quads'], reverse=True)[:20],
                  note='Saved accepted outputs only; discarded chart attempts require runtime instrumentation.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding='utf8')
    print('STREET_WORKLOAD', json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
