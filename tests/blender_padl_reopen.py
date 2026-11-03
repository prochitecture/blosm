"""Reopen both Berlin outputs, verify assets/edits, and compare both backends."""
import json
from pathlib import Path
import sys
import argparse


def top(owner, fid):
    import bpy
    bpy.context.view_layer.update()
    evaluated = owner.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh = evaluated.to_mesh()
    try:
        ids = mesh.attributes['footprint_id']
        return max((evaluated.matrix_world @ mesh.vertices[i].co).z for p in mesh.polygons
                   if ids.data[p.index].value == fid for i in p.vertices)
    finally:
        evaluated.to_mesh_clear()


def main():
    import bpy
    import addon_utils
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root.parent))
    addon_utils.enable(root.name, default_set=True)
    from blosm.building.padl_integration.backends import BuildingGenerationService
    from blosm.building.padl_integration.gn_runtime import read_status
    from blosm.padl.api import edit_parameter
    from blosm.padl.persistence import read_footprints
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',type=Path,default=root.parent / 'artifacts/padl_parts/packed')
    args = parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    folder = args.folder
    reports = {name: json.loads((folder / ('berlin_' + name + '.json')).read_text())
               for name in ('geometry_nodes', 'python')}
    left, right = (reports[name] for name in ('geometry_nodes', 'python'))
    assert left['source_sha256'] == right['source_sha256']
    assert left['bounds'].keys() <= right['bounds'].keys()
    common = left['bounds'].keys() & right['bounds'].keys()
    maximum = max(abs(a-b) for key in common for a, b in zip(left['bounds'][key], right['bounds'][key]))
    assert maximum < 0.0001, maximum
    for backend in reports:
        bpy.ops.wm.open_mainfile(filepath=str(folder / ('berlin_' + backend + '.blend')))
        owner = next(obj for obj in bpy.data.objects if obj.get('padl_source_keys'))
        output = owner if backend == 'geometry_nodes' else bpy.data.objects[owner.name + ' Volumes']
        count = len(bpy.data.node_groups)
        service = BuildingGenerationService(backend=backend)
        if backend == 'geometry_nodes':
            assert service.library.tree('simple_3d') == owner.modifiers['PADL Buildings'].node_group
            assert len(bpy.data.node_groups) == count
        # Both a constrained edit and a new explicit level count after reopening.
        footprints = read_footprints(owner)
        if backend == 'python':
            main_part = next((f for f in footprints if json.loads(f.identity)[1] == 'osm:way:840186547'),None)
            if main_part is not None:
                assert abs(top(output,main_part.footprint_id)-right['bounds'][str(main_part.footprint_id)][1]) < .0001
        successes = {int(k) for k in reports[backend]['bounds']}
        used = set()
        for supplied in (True, False):
            footprint = next(f for f in footprints if f.footprint_id in successes and
                f.footprint_id not in used and f.attributes.get('roof_type', 'flat') == 'flat' and
                f.attributes['building_use'] in ('apartments', 'house', 'office') and
                'max_elevation' not in f.attributes)
            fid = footprint.footprint_id
            used.add(fid)
            before = top(output, fid)
            if supplied:
                # Berlin's supplied heights are mostly on generic/part records.
                # Capture a current level-model building's H as an explicit edit.
                baseline = reports[backend]['bounds'][str(fid)][0]
                edit_parameter(owner, fid, 'max_elevation', before - baseline)
            old = footprint.attributes.get('levels_wall_count', 2)
            edit_parameter(owner, fid, 'levels_wall_count', old + 2)
            if backend == 'python':
                result = service.regenerate(owner, collection=output.users_collection[0],
                    terrain=bpy.data.objects['Terrain'], output_object=output)
                assert result.output == output
                assert not next(s for s in result.statuses if s['footprint_id'] == fid)['error']
            else:
                record = next(s for s in read_status(owner, owner.modifiers['PADL Buildings']) if s['footprint_id'] == fid)
                assert not record['error'], record
                assert record['levels_wall_count'] == old + 2
            current = next(f for f in read_footprints(owner) if f.footprint_id == fid)
            assert current.attributes['levels_wall_count'] == old + 2
            assert ('max_elevation' in current.attributes) == supplied
            after = top(output, fid)
            assert abs(after - before) < .0001 if supplied else abs(after - before) > .1
        print('PADL_REOPEN_EDIT_OK', backend)
    (folder / 'reopen_comparison.json').write_text(json.dumps({'common_building_count': len(common),
        'python_only_footprints':sorted(right['bounds'].keys()-left['bounds'].keys(),key=int),
        'maximum_vertical_bound_difference_metres': maximum, 'reopen_and_edit': 'passed'}, indent=2) + '\n')
    print('PADL_BERLIN_PARITY_OK', maximum)


if __name__ == '__main__':
    main()
