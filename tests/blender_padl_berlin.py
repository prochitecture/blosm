"""End-to-end simple 3D import against the user's Berlin terrain scene.

Run with Blender --background --factory-startup --python this_file --
--backend geometry_nodes|python. Writes only into workspace artifacts.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import time


def main():
    import bpy
    import addon_utils
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--backend', choices=('geometry_nodes', 'python'), required=True)
    parser.add_argument('--source', default=r'C:\prochitecture\setup\bldg2_berlin_karl_marx_allee.blend')
    parser.add_argument('--osm', default=r'C:\dev\workspace\osm_extracts\facade_visibility\berlin_karl_marx_allee.osm')
    parser.add_argument('--output', default=r'C:\dev\workspace\artifacts\padl_berlin')
    parser.add_argument('--terrain-offset-z', type=float, default=0)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    source, destination = Path(args.source), Path(args.output)
    destination.mkdir(parents=True, exist_ok=True)
    before = hashlib.sha256(source.read_bytes()).hexdigest()
    started = time.monotonic()
    bpy.ops.wm.open_mainfile(filepath=str(source))
    # The source contains a legacy service credential; do not copy it to tests.
    if 'api_key' in bpy.context.scene:
        del bpy.context.scene['api_key']
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root.parent))
    module = addon_utils.enable(root.name, default_set=True)
    assert module is not None, 'Addon registration failed'
    preferences = bpy.context.preferences.addons[root.name].preferences
    preferences.dataDir = str(destination)
    settings = bpy.context.scene.blosm
    settings.dataType, settings.mode, settings.osmSource = 'osm', '3Dsimple', 'file'
    settings.osmFilepath, settings.terrainObject = str(Path(args.osm)), 'Terrain'
    settings.buildingBackend, settings.buildingSeed = args.backend, 17
    settings.buildings, settings.singleObject = True, True
    settings.setupScript = ''
    settings.relativeToInitialImport = False
    settings.assetLibraryMode = 'APPEND'
    for name in ('water', 'forests', 'vegetation', 'highways', 'railways', 'aviation'):
        if hasattr(settings, name):
            setattr(settings, name, False)
    terrain = bpy.data.objects['Terrain']
    # The source hides Terrain globally; enable it in test copies so transform
    # edits are evaluated before recording the matrix and testing preservation.
    terrain.hide_viewport = False
    terrain.hide_set(False)
    terrain.location.z += args.terrain_offset_z
    bpy.context.view_layer.update()
    terrain_matrix = [list(row) for row in terrain.matrix_world]
    assert bpy.ops.blosm.import_data() == {'FINISHED'}
    bpy.context.view_layer.update()
    from blosm.app.blender import app
    from blosm.padl.persistence import read_footprints
    from blosm.padl.terrain import BlenderTerrainSampler
    result, report = app.padlBuildingResult, app.padlBuildingReport
    assert result is not None and report['input_count'] > 0
    assert terrain_matrix == [list(row) for row in terrain.matrix_world], ('Import moved the terrain', terrain_matrix, [list(row) for row in terrain.matrix_world])
    footprints = read_footprints(result.owner)
    assert len(footprints) == report['input_count']
    # Initial imports must leave policy/default parameters absent.
    for footprint in footprints:
        assert not {'base_elevation', 'top_height', 'level_height', 'ground_floor_height', 'ground_floor_factor'} & footprint.attributes.keys()
    sampler = BlenderTerrainSampler(terrain)
    samples = {}
    for footprint in footprints:
        try:
            sample = sampler.sample(footprint)
            samples[footprint.footprint_id] = sample.minimum
        except Exception:
            pass  # Coverage failures are required to appear in backend statuses.
    bpy.context.view_layer.update()
    evaluated = result.output.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh = evaluated.to_mesh()
    try:
        assert len(mesh.polygons) > 0, 'No building volumes generated'
        attributes = {a.name for a in mesh.attributes if not a.name.startswith('.')}
        assert {'footprint_id', 'face_type', 'component_type'} <= attributes
        assert not {'building_use', 'roof_type', 'max_elevation', 'has_max_elevation', 'padl_random_key'} & attributes
        ids = mesh.attributes['footprint_id']
        bounds = {}
        for polygon in mesh.polygons:
            fid = ids.data[polygon.index].value
            assert fid > 0
            assert mesh.attributes['face_type'].data[polygon.index].value == 4
            assert mesh.attributes['component_type'].data[polygon.index].value in (1, 2)
            zs = [(evaluated.matrix_world @ mesh.vertices[i].co).z for i in polygon.vertices]
            low, high = bounds.get(fid, (float('inf'), -float('inf')))
            bounds[fid] = min(low, *zs), max(high, *zs)
        success = {s['footprint_id'] for s in report['statuses'] if not s.get('error')}
        assert set(bounds) == success, (set(bounds) - success, success - set(bounds))
        assert all(fid in samples for fid in bounds), 'Generated a building without terrain coverage'
        for fid, (low, high) in bounds.items():
            footprint = next(f for f in footprints if f.footprint_id == fid)
            assert abs(low - samples[fid] - footprint.attributes.get('min_elevation', 0)) < 0.002, (fid, low, samples[fid])
            if 'max_elevation' in footprint.attributes:
                assert abs(high - samples[fid] - footprint.attributes['max_elevation']) < 0.01, (fid, high, samples[fid])
        messages = [message.message for mod in result.output.modifiers
                    for message in getattr(mod, 'node_warnings', [])]
        # Native solver warnings for rejected real-world footprints are valid
        # only when the separate Status output also reports those failures.
        assert not messages or any(row.get('error') for row in report['statuses']), messages
        assert all(message == 'ERR_NON_CONVEX' for message in messages), messages
        details = {'backend': args.backend, 'source_sha256': before,
            'elapsed_seconds': time.monotonic() - started,
            'input_count': len(footprints), 'generated_count': len(bounds),
            'vertices': len(mesh.vertices), 'faces': len(mesh.polygons),
            'native_solver_messages': messages,
            'mapping_issues': dict(Counter(row['code'] for row in report['issues'])),
            'generation_errors': dict(Counter(str(row['error']) for row in report['statuses'] if row.get('error'))),
            'building_uses': dict(Counter(f.attributes['building_use'] for f in footprints)),
            'roof_types': dict(Counter(f.attributes.get('roof_type', 'default') for f in footprints)),
            'bounds': bounds, 'report': report}
    finally:
        evaluated.to_mesh_clear()
    # Check reuse through a completely new service instance.
    if args.backend == 'geometry_nodes':
        from blosm.building.padl_integration.assets import BuildingProgramLibrary
        count = len(bpy.data.node_groups)
        assert BuildingProgramLibrary().tree('simple_3d') == result.output.modifiers['PADL Buildings'].node_group
        assert count == len(bpy.data.node_groups), 'Repeated load duplicated GN groups'
    else:
        assert not any(g.get('padl_compiler_version') for g in bpy.data.node_groups), 'Python imported the GN compiler/library'
    output = destination / ('berlin_' + args.backend + '.blend')
    bpy.ops.wm.save_as_mainfile(filepath=str(output))
    assert hashlib.sha256(source.read_bytes()).hexdigest() == before, 'Original scene was modified'
    (destination / ('berlin_' + args.backend + '.json')).write_text(json.dumps(details, indent=2) + '\n')
    print('PADL_BERLIN_OK', args.backend, details['generated_count'], 'buildings;', details['faces'], 'faces;', output)


if __name__ == '__main__':
    main()
