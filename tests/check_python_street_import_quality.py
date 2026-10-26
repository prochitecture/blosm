"""Audit published Blender buffers, including float32 coordinate precision."""
import argparse
import json
from pathlib import Path
import sys
import bpy
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT / 'tests'))
from check_python_street_scene import mesh_buffers
from blosm.way.street_network.python_backend.model import Network
from blosm.way.street_network.python_backend.model import Key
from blosm.way.street_network.python_backend.mesh import Mesh
from blosm.way.street_network.python_backend.validation import validate_overlaps


def saved_components(key_mode='section-set'):
    """Use Section sets for fresh imports, or persistent UUIDs for galleries."""
    result = {}
    for controller in (o for o in bpy.data.objects if o.get('st_py_source_data')):
        network = Network.from_dict(json.loads(controller['st_py_source_data']))
        sections = tuple(sorted(s.id for s in network.sections))
        key = (network.uuid, sections) if key_mode == 'network-uuid' else sections
        assert key not in result, ('Ambiguous reference component identity', key_mode, key)
        result[key] = mesh_buffers(controller)
    return result


def boundary_segments(buffers):
    """Canonical XYZ segments, independent of interior vertex/face numbering."""
    xyz, quads = buffers[0].reshape(-1, 3), buffers[1].reshape(-1, 4)
    edges = np.sort(np.stack((quads, np.roll(quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    return tuple(sorted(tuple(sorted((tuple(xyz[a]), tuple(xyz[b])))) for a, b in unique[counts == 1]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference-scene', type=Path)
    parser.add_argument('--reference-key', choices=('section-set', 'network-uuid'), default='section-set')
    parser.add_argument('--allow-new-network', action='append', default=[])
    parser.add_argument('--capture-elongated', action='store_true')
    parser.add_argument('--capture-part', action='append', default=[])
    parser.add_argument('--allow-changed-part', action='append', default=[])
    parser.add_argument('--preserve-boundary-part', action='append', default=[])
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    reference = None
    if args.reference_scene:
        bpy.ops.wm.open_mainfile(filepath=str(args.reference_scene.resolve()))
        reference = saved_components(args.reference_key)
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()))
    changed_parts = set()
    preserved_boundaries = set()
    unchanged_parts = 0
    added_networks = []
    if reference is not None:
        current = saved_components(args.reference_key)
        assert reference.keys() <= current.keys(), 'A reference component disappeared'
        additional = current.keys() - reference.keys()
        if args.allow_new_network:
            assert args.reference_key == 'network-uuid', 'Named additions require network UUID comparison'
            added_networks = sorted(key[0] for key in additional)
            assert added_networks == sorted(args.allow_new_network), 'Unexpected added reference components'
        else:
            assert not additional, 'Component membership changed'
        for component in reference:
            meshes = current[component]
            previous = reference[component]
            for key in meshes.keys() | previous.keys():
                if key in args.preserve_boundary_part:
                    assert key in meshes and key in previous
                    assert boundary_segments(meshes[key]) == boundary_segments(previous[key]), (
                        'Prescribed physical boundary changed', key)
                    preserved_boundaries.add(key)
                if (key not in meshes or key not in previous or
                        not np.array_equal(meshes[key][0], previous[key][0]) or
                        not np.array_equal(meshes[key][1], previous[key][1])):
                    changed_parts.add(key)
                else:
                    unchanged_parts += 1
            assert changed_parts <= set(args.allow_changed_part), ('Unexpected geometry/ownership changes', changed_parts)
        del current, reference
        assert preserved_boundaries == set(args.preserve_boundary_part), 'Requested boundary part was not compared'
    elif args.preserve_boundary_part:
        raise ValueError('Boundary preservation comparison needs a reference scene')
    results = []
    for component, controller in enumerate(o for o in bpy.data.objects if o.get('st_py_source_data')):
        network = Network.from_dict(json.loads(controller['st_py_source_data']))
        datum = min(s.datum for s in network.sections)
        rise = max([s.rise for s in network.sections] + [b.profile_rise or 0. for b in network.bundles])
        ceiling = max(s.datum for s in network.sections) + controller['st_py_exaggeration'] * rise
        parts, published_meshes = [], {}
        objects = {o['st_py_key']: o for o in controller['st_py_outputs'].objects if o.type == 'MESH'}
        for key, (xyz, quads) in mesh_buffers(controller).items():
            counts = np.empty(len(objects[key].data.polygons), dtype=np.int32)
            objects[key].data.polygons.foreach_get('loop_total', counts)
            assert np.all(counts == 4), key
            xyz, quads = xyz.reshape(-1, 3), quads.reshape(-1, 4)
            assert np.all(np.isfinite(xyz)) and quads.shape[1] == 4
            minimum, maximum = np.inf, 0.
            over = 0
            for start in range(0, len(quads), 4096):
                p = xyz[quads[start:start + 4096], :2].astype(float)
                edges = np.roll(p, -1, axis=1) - p
                following = np.roll(edges, -1, axis=1)
                determinants = edges[..., 0] * following[..., 1] - edges[..., 1] * following[..., 0]
                minimum = min(minimum, float(determinants.min()))
                lengths = np.sqrt(np.sum(edges * edges, axis=-1))
                ratios = lengths.max(axis=1) / lengths.min(axis=1)
                maximum = max(maximum, float(ratios.max()))
                over += int(np.count_nonzero(ratios > 8.))
            assert minimum > 1e-10, (key, minimum)
            assert float(xyz[:, 2].min()) >= datum - 1e-6, key
            assert float(xyz[:, 2].max()) <= ceiling + 1e-6, key
            published_meshes[Key.parse(key)] = Mesh(xyz, quads)
            parts.append(dict(key=key, quads=len(quads), minimum_corner=minimum,
                maximum_edge_ratio=maximum, cells_above_ratio_8=over,
                minimum_height=float(xyz[:, 2].min()), maximum_height=float(xyz[:, 2].max())))
            if (args.capture_elongated and over) or key in args.capture_part:
                capture = args.output / f'component_{component}_{key.replace(":", "_")}.npz'
                np.savez_compressed(capture, xyz=xyz, quads=quads)
                parts[-1].update(capture=str(capture.resolve()),
                                metadata=json.loads(objects[key]['st_py_metadata']))
        # Audit the actual float32 publication, including inter-part seams and
        # self-overlaps, rather than relying only on the pre-publication check.
        overlap = validate_overlaps(published_meshes, cell_size=max(1., network.settings.spacing * 2.))
        results.append(dict(network=network.uuid, sections=len(network.sections), datum=datum, ceiling=ceiling,
            quads=sum(p['quads'] for p in parts), parts=parts, overlap_check=overlap))
    passed = all(p['cells_above_ratio_8'] == 0 for r in results for p in r['parts'])
    report = dict(passed=passed, scene=str(args.scene.resolve()), all_convex=True, all_quads=True,
        all_within_quad_ratio_limit=passed, no_published_overlaps=True,
        within_global_crown_envelope=True, networks=results)
    if args.reference_scene:
        report.update(reference_scene=str(args.reference_scene.resolve()),
                      reference_key=args.reference_key,
                      added_networks=added_networks,
                      exact_reference_positions=not changed_parts, exact_reference_topology=not changed_parts,
                      reference_changed_parts=sorted(changed_parts), unchanged_reference_parts=unchanged_parts,
                      exact_preserved_boundaries=sorted(preserved_boundaries),
                      unchanged_reference_parts_exact=True)
    (args.output / 'import_quality_checks.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(dict(passed=passed, networks=len(results), quads=sum(r['quads'] for r in results),
        maximum_edge_ratio=max(p['maximum_edge_ratio'] for r in results for p in r['parts']),
        cells_above_ratio_8=sum(p['cells_above_ratio_8'] for r in results for p in r['parts']))), flush=True)
    return 0 if passed else 1


if __name__ == '__main__':
    sys.exit(main())
