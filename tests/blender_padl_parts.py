"""Building parts: both backends, common datum, courtyards, packing and reopen."""
from pathlib import Path
from dataclasses import replace
import importlib
import os
import sys

def rect(x, y, w, h):
    return ((x,y),(x+w,y),(x+w,y+h),(x,y+h))

def assert_bounds(a,b):
    assert set(a) == set(b), (set(a),set(b))
    assert max((abs(x-y) for key in a for x,y in zip(a[key],b[key])),default=0) < .0001, (a,b)

def bounds(output):
    import bpy
    bpy.context.view_layer.update()
    evaluated = output.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh = evaluated.to_mesh()
    try:
        result = {}
        for face in mesh.polygons:
            fid = mesh.attributes['footprint_id'].data[face.index].value
            low, high = result.get(fid, (float('inf'), -float('inf')))
            zs = [(evaluated.matrix_world @ mesh.vertices[i].co).z for i in face.vertices]
            result[fid] = min(low, *zs), max(high, *zs)
        assert not {'parent_id', 'has_holes', 'building_use', 'padl_random_key'} & {a.name for a in mesh.attributes}
        return result
    finally:
        evaluated.to_mesh_clear()

def main():
    import bpy, addon_utils
    root = Path(__file__).resolve().parents[1]
    package = os.environ.get('BLOSM_TEST_PACKAGE', root.name)
    if not package.startswith('bl_ext.'):
        sys.path.insert(0, str(root.parent))
    assert addon_utils.enable(package, default_set=True) is not None
    Footprint = importlib.import_module(package + '.padl.geometry').Footprint
    read_footprints = importlib.import_module(package + '.padl.persistence').read_footprints
    api = importlib.import_module(package + '.padl.api')
    edit_parameter, create_input_object = api.edit_parameter, api.create_input_object
    PadlError = importlib.import_module(package + '.padl.errors').PadlError
    BuildingGenerationService = importlib.import_module(package + '.building.padl_integration.backends').BuildingGenerationService
    group_buildings = importlib.import_module(package + '.building.padl_integration.source').group_buildings
    collection = bpy.context.collection
    terrain_mesh = bpy.data.meshes.new('Parts terrain')
    terrain_mesh.from_pydata([(-10,-10,8),(80,-10,26),(80,30,26),(-10,30,8)], [], [(0,1,2,3)])
    terrain = bpy.data.objects.new('Parts terrain', terrain_mesh)
    collection.objects.link(terrain)
    def f(key, ring, role=1, parent=0, holes=(), **attrs):
        return Footprint('test:' + str(key), key, ring, holes,
            dict(building_use='apartments', roof_type='flat', max_elevation=20., levels_wall_count=4, **attrs), role, parent)
    footprints = [
        f(100, rect(0,0,20,20), 5),
        f(101, rect(10,0,10,20), 3, 100, (rect(12,3,3,4),), min_elevation=5.),
        f(102, rect(0,0,10,20), 3, 100),
        f(200, rect(30,0,20,20), 6, holes=(rect(42,4,3,3),)),
        f(201, rect(40,0,6,20), 3, 200, min_elevation=10.),
        f(300, rect(60,0,5,8), 6),
    ]
    # Invalid ownership must be refused before either backend sees the mesh.
    for rows in ([f(901, rect(0,0,10,10),3,999)],
                 [f(901,rect(0,0,10,10)),f(902,rect(1,1,2,2),3,901)],
                 [f(901,rect(0,0,10,10),5)]):
        try:
            create_input_object(rows)
        except PadlError:
            pass
        else:
            raise AssertionError('Invalid part ownership was admitted')
    outputs = {}
    for backend in ('geometry_nodes', 'python'):
        service = BuildingGenerationService(backend=backend)
        packed = service.generate(footprints, collection=collection, terrain=terrain, seed=17, name='Parts ' + backend)
        assert len(packed.statuses) == 5 and all(not s['error'] for s in packed.statuses), packed.statuses
        before = bounds(packed.output)
        assert set(before) == {101,102,200,201,300}, before
        for fid, expected in {101:(15,30),102:(10,30),200:(16,36),201:(26,36),300:(22,42)}.items():
            assert max(abs(a-b) for a,b in zip(before[fid], expected)) < .0001, (fid, before[fid], expected)
        read = read_footprints(packed.owner)
        assert [(f.face_type,f.parent_id) for f in read] == [(f.face_type,f.parent_id) for f in footprints]
        assert next(f for f in read if f.footprint_id == 101).elevation_reference.footprint_id == 100
        group_count = len(bpy.data.node_groups)
        with service.batch():
            for group in group_buildings(footprints, False):
                one = service.generate(group, collection=collection, terrain=terrain, seed=17)
                for fid, pair in bounds(one.output).items():
                    assert max(abs(a-b) for a,b in zip(pair,before[fid])) < .0001
        assert len(bpy.data.node_groups) == group_count, 'Batch leaked or duplicated GN groups'
        edit_parameter(packed.owner, 101, 'levels_wall_count', 7)
        if backend == 'python':
            service.regenerate(packed.owner, collection=collection, terrain=terrain, output_object=packed.output)
        assert_bounds(bounds(packed.output),before)
        # A courtyard is open: no roof face may contain a point inside its hole.
        obj = packed.output.evaluated_get(bpy.context.evaluated_depsgraph_get())
        mesh = obj.to_mesh()
        try:
            from mathutils import Vector
            from mathutils.geometry import intersect_point_tri_2d
            point = Vector((13.,4.))
            mesh.calc_loop_triangles()
            for triangle in mesh.loop_triangles:
                face = mesh.polygons[triangle.polygon_index]
                if mesh.attributes['footprint_id'].data[face.index].value != 101 or mesh.attributes['component_type'].data[face.index].value != 2:
                    continue
                a,b,c = [obj.matrix_world @ mesh.vertices[i].co for i in triangle.vertices]
                assert not intersect_point_tri_2d(point,a.xy,b.xy,c.xy), 'Roof covers part courtyard'
        finally:
            obj.to_mesh_clear()
        outputs[backend] = before
        if backend == 'geometry_nodes':
            assert not list(packed.output.modifiers['PADL Buildings'].node_warnings)
        invalid = replace(footprints[1], identity='test:103', footprint_id=103,
                          attributes=dict(footprints[1].attributes, roof_type='gabled', roof_height=3.))
        mixed = service.generate([*footprints,invalid], collection=collection, terrain=terrain, seed=17)
        assert next(s for s in mixed.statuses if s['footprint_id'] == 103)['error']
        assert_bounds(bounds(mixed.output),before)
        # A covered part cannot use an incompletely projected building datum.
        limited_mesh = bpy.data.meshes.new('Limited terrain')
        limited_mesh.from_pydata([(1,-1,10),(21,-1,10),(21,21,10),(1,21,10)],[],[(0,1,2,3)])
        limited = bpy.data.objects.new('Limited terrain',limited_mesh)
        collection.objects.link(limited)
        failed = service.generate(footprints[:2], collection=collection, terrain=limited, seed=17)
        assert not bounds(failed.output)
        assert all(s['error'] for s in failed.statuses)
        outside = f(104,rect(75,0,10,10),3,100)
        partial = service.generate([footprints[0],footprints[2],outside], collection=collection, terrain=terrain, seed=17)
        assert set(bounds(partial.output)) == {102}
        assert next(s for s in partial.statuses if s['footprint_id'] == 104)['error']
        print('PARTS_FIXTURE_OK', backend, flush=True)
    assert_bounds(outputs['geometry_nodes'],outputs['python'])
    destination = Path(os.environ.get('BLOSM_TEST_ARTIFACTS', root.parent / 'artifacts/padl_parts')) / 'parts_fixtures.blend'
    destination.parent.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(destination))
    bpy.ops.wm.open_mainfile(filepath=str(destination))
    for backend in ('geometry_nodes','python'):
        owner = bpy.data.objects['Parts ' + backend]
        assert len(read_footprints(owner)) == 6
        output = owner if backend == 'geometry_nodes' else bpy.data.objects['Parts python Volumes']
        assert_bounds(bounds(output),outputs[backend])
    print('PARTS_REOPEN_OK', flush=True)

if __name__ == '__main__':
    main()
