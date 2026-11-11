"""Inspect numerical centerline fitting from a saved street scene, without rebuilding meshes."""
import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import bpy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _python_street_package
from tested_python_streets.model import Network
from tested_python_streets.providers import PathProvider


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--section', type=int, action='append', default=[], help='Also export the saved pavement mesh for these Sections')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()), load_ui=False, use_scripts=False)
    results, networks, pavements = [], [], []
    for obj in bpy.data.objects:
        if not obj.get('st_py_source_data'):
            continue
        network = Network.from_dict(json.loads(obj['st_py_source_data']))
        networks.append(network.to_dict())
        outputs = obj.get('st_py_outputs')
        if outputs:
            for mesh_obj in outputs.objects:
                if mesh_obj.type == 'MESH' and mesh_obj.get('st_py_key') in {f'road:{s}' for s in args.section}:
                    xyz = np.empty((len(mesh_obj.data.vertices), 3), np.float32)
                    mesh_obj.data.vertices.foreach_get('co', xyz.ravel())
                    pavements.append(dict(network=network.uuid, key=mesh_obj['st_py_key'],
                        xy=xyz[:, :2].tolist(), faces=[list(p.vertices) for p in mesh_obj.data.polygons]))
        for section in network.sections:
            seams = [j for j in (*network.intersections, *network.connections) if len(j.ends) == 2
                     and any(e.section_id == section.id for e in j.ends)]
            record = dict(network=network.uuid, section=section.id, street=section.street_id,
                          width=section.width, source=section.points)
            try:
                path = PathProvider.build(SimpleNamespace(network=network), section, seams)
                query = path.point(np.linspace(0., path.length, 1001))
                points = np.asarray(section.points)[:, :2]
                delta = np.diff(points, axis=0)
                t = np.clip(np.sum((query[:, None] - points[None, :-1]) * delta[None], axis=2)
                            / np.sum(delta * delta, axis=1)[None], 0., 1.)
                distance = np.linalg.norm(query[:, None] - points[None, :-1]
                                          - t[..., None] * delta[None], axis=2).min(axis=1)
                record.update(localized=path.localized, alignment_fallback=path.alignment_fallback, deviation=float(distance.max()),
                              path=query[::4].tolist(), length=path.length)
            except ValueError as exc:
                record['error'] = str(exc)
            results.append(record)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(networks=networks, sections=results, pavements=pavements), indent=2), encoding='utf8')
    largest = sorted((r for r in results if 'deviation' in r), key=lambda r: r['deviation'], reverse=True)[:20]
    print(json.dumps([dict(section=r['section'], street=r['street'], width=r['width'],
                           deviation=r['deviation'], localized=r['localized'], source=r['source'])
                      for r in largest]), flush=True)


if __name__ == '__main__':
    main()
