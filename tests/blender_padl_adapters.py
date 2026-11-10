"""Installed-addon strategies: all roofs, holes, packing, reuse and strict terrain."""
from pathlib import Path
import importlib
import os
import sys


def bounds(output):
    import bpy
    bpy.context.view_layer.update()
    obj = output.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh = obj.to_mesh()
    try:
        result = {}
        for face in mesh.polygons:
            fid = mesh.attributes['footprint_id'].data[face.index].value
            points = [obj.matrix_world @ mesh.vertices[i].co for i in face.vertices]
            low, high = result.get(fid, (float('inf'), -float('inf')))
            result[fid] = min(low, *(v.z for v in points)), max(high, *(v.z for v in points))
        return result
    finally:
        obj.to_mesh_clear()


def main():
    import bpy
    import addon_utils
    root = Path(__file__).resolve().parents[1]
    package = os.environ.get('BLOSM_TEST_PACKAGE', root.name)
    if not package.startswith('bl_ext.'):
        sys.path.insert(0, str(root.parent))
    assert addon_utils.enable(package, default_set=True) is not None
    BuildingGenerationService = importlib.import_module(package + '.building.padl_integration.backends').BuildingGenerationService
    Footprint = importlib.import_module(package + '.padl.geometry').Footprint
    read_footprints = importlib.import_module(package + '.padl.persistence').read_footprints
    OsmFaceAttributeMapper = importlib.import_module(package + '.building.padl_integration.mapping').OsmFaceAttributeMapper
    terrain_mesh = bpy.data.meshes.new('Terrain')
    terrain_mesh.from_pydata([(-10,-10,0), (150,-10,0), (150,30,0), (-10,30,0)], [], [(0,1,2,3)])
    terrain = bpy.data.objects.new('Terrain', terrain_mesh)
    bpy.context.collection.objects.link(terrain)
    terrain.location.z = 36
    bpy.context.view_layer.update()
    footprints = []
    for i, roof in enumerate(('flat', 'skillion', 'gabled', 'hipped'), 1):
        x = (i - 1) * 20
        hole = (((x+3,2), (x+3,4), (x+5,4), (x+5,2)),) if roof in ('flat', 'skillion') else ()
        footprints.append(Footprint('osm:way:' + str(i), i,
            ((x,0), (x+12,0), (x+12,8), (x,8)), hole,
            {'building_use': 'house', 'roof_type': roof, 'max_elevation': 14.0,
             'roof_height': 4.0, 'levels_wall_count': 3}))
    outputs = {}
    for i,shape in enumerate(('gambrel','round','saltbox'),5):
        x=(i-1)*20
        attrs=OsmFaceAttributeMapper().map({'building':'house','roof:shape':shape,
            'roof:direction':'NE','height':'14','roof:height':'4','building:levels':'3'}).attributes
        footprints.append(Footprint('osm:way:'+str(i),i,((x,0),(x+12,0),(x+12,8),(x,8)),attributes=attrs))
    for backend in ('geometry_nodes', 'python'):
        service = BuildingGenerationService(backend=backend)
        packed = service.generate(footprints, collection=bpy.context.collection, terrain=terrain, seed=71)
        assert all(not s['error'] for s in packed.statuses), packed.statuses
        assert len(packed.statuses) == len(footprints)
        outputs[backend] = bounds(packed.output)
        assert all(abs(low-36) < .0001 and abs(high-50) < .0001 for low, high in outputs[backend].values())
        group_count = len(bpy.data.node_groups)
        random_keys = {f.identity: f.attributes for f in read_footprints(packed.owner)}
        for footprint in footprints:
            one = service.generate([footprint], collection=bpy.context.collection, terrain=terrain, seed=71)
            assert all(not s['error'] for s in one.statuses), one.statuses
            pair = bounds(one.output)[footprint.footprint_id]
            assert max(abs(a-b) for a,b in zip(pair, outputs[backend][footprint.footprint_id])) < .0001
            f = read_footprints(one.owner)[0]
            assert random_keys[f.identity] == f.attributes
        assert len(bpy.data.node_groups) == group_count, 'Dedicated objects duplicated GN trees'
        # A missing hit must omit only the uncovered building, without fallback.
        outside = Footprint('osm:way:outside', 20, ((160,0), (172,0), (172,8), (160,8)), (),
                            {'building_use': 'house', 'roof_type': 'flat'})
        mixed = service.generate([footprints[0], outside], collection=bpy.context.collection, terrain=terrain, seed=71)
        assert set(bounds(mixed.output)) == {1}
        assert next(s for s in mixed.statuses if s['footprint_id'] == 20)['error']
        if backend == 'geometry_nodes':
            assert not list(packed.output.modifiers['PADL Buildings'].node_warnings)
        print('PADL_ADAPTERS_OK', backend)
    assert outputs['geometry_nodes'] == outputs['python']


if __name__ == '__main__':
    main()
