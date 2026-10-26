"""Background-Blender publication/edit/save checks and a combined inspection.

Run with --python-exit-code 1 --python this_file -- --output <directory>.
No GN solver is evaluated. The optional text bootstrap registers only the
Python editing panel when inspecting the result without installing this branch.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import bpy
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT / 'tests'))
import blosm
from blosm.way.street_network.python_backend import model
from blosm.way.street_network.python_backend.scene import NetworkScene, SourceObjects, clear_runtime
from blosm.way.street_network.python_backend import ui
from _python_angular_fixtures import angular_fixture_networks, angular_chain_fixture_networks
from _python_street_fixtures import (fixture_networks, ownership_fixture_networks, height_fixture_networks,
                                    transition_fixture_networks, clearance_fixture_networks, side_junction_fixture_networks,
                                    interior_fixture_networks, consumed_cap_fixture_networks, locality_fixture_networks,
                                    terminal_fixture_networks, radius_fixture_networks, source_support_fixture_networks,
                                    attached_overlap_fixture_networks, multisection_attached_fixture_networks,
                                    skew_bundle_fixture_networks, dimensional_bundle_fixture_networks,
                                    one_sided_free_cap_fixture_networks, absorbed_transition_fixture_networks,
                                    nonparallel_terminal_fixture, higher_valence_fixture_networks, absorbed_profile_fixture_networks,
                                    two_merged_terminal_fixture, exact_bent_seam_fixture, sharp_exact_seam_fixture,
                                    compound_seam_fixture_networks, compound_internal_crown_fixture,
                                    compound_route_fixture_networks, sharp_smooth_seam_fixture)


def mesh_buffers(controller):
    result = {}
    for obj in controller['st_py_outputs'].objects:
        if obj.type != 'MESH':
            continue
        xyz = np.empty(len(obj.data.vertices) * 3, dtype=np.float32)
        ids = np.empty(len(obj.data.loops), dtype=np.int32)
        obj.data.vertices.foreach_get('co', xyz)
        obj.data.loops.foreach_get('vertex_index', ids)
        result[obj['st_py_key']] = (xyz, ids)
    return result


def section_control(controller, section_id):
    """Find an interior control through edge ownership, independent of grouping."""
    for source in controller['st_py_sources'].objects:
        degree = {}
        ids = source.data.attributes['st_section_id']
        for edge in source.data.edges:
            if ids.data[edge.index].value == section_id:
                for vertex in edge.vertices:
                    degree[vertex] = degree.get(vertex, 0) + 1
        controls = [i for i, count in degree.items() if count == 2]
        if controls:
            return source, source.data.vertices[controls[len(controls) // 2]]
    raise AssertionError(f'No interior source control for Section {section_id}')


def verify_seamed_approach(controller):
    before = mesh_buffers(controller)
    assert 'road:1' not in before, 'Fully consumed Section retained duplicate pavement'
    assert {'road:4', 'junction_mesh:10'} <= before.keys()
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:5'].as_pointer()
    anchor = tuple(outputs['placement:81'].location)
    source, vertex = section_control(controller, 4)
    original = vertex.co.copy()
    vertex.co.x += .2
    source.data.update()
    NetworkScene.update(controller)
    after = mesh_buffers(controller)
    assert not np.array_equal(before['junction_mesh:10'][0], after['junction_mesh:10'][0])
    assert not np.array_equal(before['road:4'][0], after['road:4'][0])
    edited_outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert tuple(edited_outputs['placement:81'].location) != anchor
    assert remote == edited_outputs['road:5'].as_pointer()
    rebuilt = json.loads(controller['st_py_last_update'])['rebuilt']
    assert 'road:5' not in rebuilt
    vertex.co = original
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0])
               and np.array_equal(before[k][1], restored[k][1]) for k in before)


def verify_compound(controller):
    before = mesh_buffers(controller)
    assert {'junction_mesh:4', 'junction_grass:4', 'separator:1', 'separator:2'} <= before.keys()
    assert all('road:' + str(sid) not in before for sid in (2, 5, 8, 11))
    poses = {o['st_py_key']: tuple(o.location) for o in controller['st_py_outputs'].objects if o.get('st_py_key') == 'placement:81'}
    source, vertex = section_control(controller, 8)
    original = vertex.co.copy()
    vertex.co.x += .1
    source.data.update()
    NetworkScene.update(controller)
    after = mesh_buffers(controller)
    assert not np.array_equal(before['junction_mesh:4'][0], after['junction_mesh:4'][0])
    rebuilt = json.loads(controller['st_py_last_update'])['rebuilt']
    assert {'junction_mesh:4', 'junction_grass:4', 'separator:1', 'separator:2', 'placement:81'} <= set(rebuilt)
    changed_poses = {o['st_py_key']: tuple(o.location) for o in controller['st_py_outputs'].objects if o.get('st_py_key') == 'placement:81'}
    assert poses != changed_poses
    vertex.co = original
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[key][0], restored[key][0]) and np.array_equal(before[key][1], restored[key][1]) for key in before)


def verify_projected(controller):
    from blosm.way.street_network.python_backend.mesh import Mesh
    from blosm.way.street_network.python_backend.surface_query import QuadHeightField
    before = mesh_buffers(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert json.loads(outputs['junction_mesh:1']['st_py_metadata'])['height_projection'] is not None
    field = QuadHeightField(Mesh(before['junction_mesh:1'][0].reshape(-1, 3),
        before['junction_mesh:1'][1].reshape(-1, 4)), datum=0., ceiling=.81)
    pose = tuple(outputs['placement:81'].location)
    assert abs(field(np.array((pose[:2],)))[0] - pose[2]) < 2e-6
    pointers = {key: obj.as_pointer() for key, obj in outputs.items() if obj.type == 'MESH'}
    source_data = json.loads(controller['st_py_source_data'])
    station = source_data['furniture'][0]['station']
    source_data['furniture'][0]['station'] += .15
    controller['st_py_source_data'] = json.dumps(source_data)
    NetworkScene.update(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert pointers == {key: obj.as_pointer() for key, obj in outputs.items() if obj.type == 'MESH'}
    assert json.loads(controller['st_py_last_update'])['rebuilt'] == ['placement:81']
    changed = tuple(outputs['placement:81'].location)
    assert changed != pose and abs(field(np.array((changed[:2],)))[0] - changed[2]) < 2e-6
    source_data['furniture'][0]['station'] = station
    controller['st_py_source_data'] = json.dumps(source_data)
    NetworkScene.update(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert tuple(outputs['placement:81'].location) == pose
    source, vertex = section_control(controller, 3)
    original = vertex.co.copy()
    vertex.co.x += .05
    source.data.update()
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    assert not np.array_equal(before['junction_mesh:1'][0], edited['junction_mesh:1'][0])
    NetworkScene.update(controller, force=True)
    fresh = mesh_buffers(controller)
    assert edited.keys() == fresh.keys()
    assert all(np.array_equal(edited[k][0], fresh[k][0]) and np.array_equal(edited[k][1], fresh[k][1]) for k in edited)
    vertex.co = original
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0]) and np.array_equal(before[k][1], restored[k][1]) for k in before)


def verify_transition(controller):
    """Transition ownership, publication-buffer release and width dependencies."""
    before = mesh_buffers(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert {'seam_mesh:6', 'seam_mesh:7', 'separator:7'} <= outputs.keys()
    pose = tuple(outputs['placement:81'].location)
    pointers = {key: obj.as_pointer() for key, obj in outputs.items() if obj.type == 'MESH'}
    source_data = json.loads(controller['st_py_source_data'])
    station = source_data['furniture'][0]['station']
    source_data['furniture'][0]['station'] += .2
    controller['st_py_source_data'] = json.dumps(source_data)
    NetworkScene.update(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert pointers == {key: obj.as_pointer() for key, obj in outputs.items() if obj.type == 'MESH'}
    assert json.loads(controller['st_py_last_update'])['rebuilt'] == ['placement:81']
    assert tuple(outputs['placement:81'].location) != pose
    source_data['furniture'][0]['station'] = station
    controller['st_py_source_data'] = json.dumps(source_data)
    NetworkScene.update(controller)
    source = next(o for o in controller['st_py_sources'].objects
                  if 2 in {v.value for v in o.data.attributes['st_section_id'].data})
    ids, widths = source.data.attributes['st_section_id'], source.data.attributes['st_width']
    changed = [(i, widths.data[i].value) for i in range(len(ids.data)) if ids.data[i].value == 2]
    for i, value in changed:
        widths.data[i].value = value + .4
    source.data.update()
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    assert not np.array_equal(before['seam_mesh:6'][0], edited['seam_mesh:6'][0])
    assert not np.array_equal(before['separator:7'][0], edited['separator:7'][0])
    NetworkScene.update(controller, force=True)
    full = mesh_buffers(controller)
    assert edited.keys() == full.keys()
    assert all(np.array_equal(edited[k][0], full[k][0]) and np.array_equal(edited[k][1], full[k][1]) for k in full)
    for i, value in changed:
        widths.data[i].value = value
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0]) and np.array_equal(before[k][1], restored[k][1]) for k in before)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert tuple(outputs['placement:81'].location) == pose


def verify_clearance(controller):
    """A remote arm's width controls the adjacent transition's physical budget."""
    before = mesh_buffers(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:5'].as_pointer()
    pose = tuple(outputs['placement:81'].location)
    source = next(o for o in controller['st_py_sources'].objects
                  if 3 in {v.value for v in o.data.attributes['st_section_id'].data})
    ids, widths = source.data.attributes['st_section_id'], source.data.attributes['st_width']
    changed = [(i, widths.data[i].value) for i in range(len(ids.data)) if ids.data[i].value == 3]
    for i, value in changed:
        widths.data[i].value = value + .4
    source.data.update()
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    rebuilt = json.loads(controller['st_py_last_update'])['rebuilt']
    assert 'seam_surface:6' in rebuilt and 'path:1' not in rebuilt and 'road:5' not in rebuilt
    assert not np.array_equal(before['seam_mesh:6'][0], edited['seam_mesh:6'][0])
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert outputs['road:5'].as_pointer() == remote
    assert tuple(outputs['placement:81'].location) != pose
    metadata = json.loads(outputs['road:1']['st_py_metadata'])
    assert abs(metadata['end'] - metadata['start'] - .25) < 1e-8
    NetworkScene.update(controller, force=True)
    full = mesh_buffers(controller)
    assert edited.keys() == full.keys()
    assert all(np.array_equal(edited[k][0], full[k][0]) and np.array_equal(edited[k][1], full[k][1]) for k in full)
    for i, value in changed:
        widths.data[i].value = value
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0]) and np.array_equal(before[k][1], restored[k][1]) for k in before)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert tuple(outputs['placement:81'].location) == pose


def verify_curb_continuation(controller):
    before = mesh_buffers(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:3'].as_pointer()
    metadata = json.loads(outputs['separator:7']['st_py_metadata'])
    assert len(metadata['runs']) == 1
    source = next(o for o in controller['st_py_sources'].objects
                  if 4 in {v.value for v in o.data.attributes['st_section_id'].data})
    ids, widths = source.data.attributes['st_section_id'], source.data.attributes['st_width']
    changed = [(i, widths.data[i].value) for i in range(len(ids.data)) if ids.data[i].value == 4]
    for i, value in changed:
        widths.data[i].value = value + .4
    source.data.update()
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    rebuilt = json.loads(controller['st_py_last_update'])['rebuilt']
    assert 'separator:7' in rebuilt and 'road:3' not in rebuilt
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert outputs['road:3'].as_pointer() == remote
    assert len(json.loads(outputs['separator:7']['st_py_metadata'])['runs']) == 1
    NetworkScene.update(controller, force=True)
    fresh = mesh_buffers(controller)
    assert edited.keys() == fresh.keys()
    assert all(np.array_equal(edited[k][0], fresh[k][0]) and np.array_equal(edited[k][1], fresh[k][1]) for k in fresh)
    for i, value in changed:
        widths.data[i].value = value
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0]) and np.array_equal(before[k][1], restored[k][1]) for k in before)


def verify_interior_continuation(controller):
    before = mesh_buffers(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote, pose = outputs['road:99'].as_pointer(), tuple(outputs['placement:82'].location)
    metadata = json.loads(outputs['separator:7']['st_py_metadata'])
    assert len(metadata['runs']) == 2 and len(metadata['derived_continuations']) == 1
    assert metadata['derived_continuations'][0]['boundary_street_id'] == 5
    source = next(o for o in controller['st_py_sources'].objects
                  if 5 in {v.value for v in o.data.attributes['st_section_id'].data})
    ids, widths = source.data.attributes['st_section_id'], source.data.attributes['st_width']
    changed = [(i, widths.data[i].value) for i in range(len(ids.data)) if ids.data[i].value == 5]
    for i, value in changed:
        widths.data[i].value = value + .8
    source.data.update()
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    rebuilt = json.loads(controller['st_py_last_update'])['rebuilt']
    assert 'separator:7' in rebuilt and 'road:99' not in rebuilt
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert outputs['road:99'].as_pointer() == remote
    assert tuple(outputs['placement:82'].location) != pose
    NetworkScene.update(controller, force=True)
    fresh = mesh_buffers(controller)
    assert edited.keys() == fresh.keys()
    assert all(np.array_equal(edited[k][0], fresh[k][0]) and np.array_equal(edited[k][1], fresh[k][1]) for k in fresh)
    for i, value in changed:
        widths.data[i].value = value
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0]) and np.array_equal(before[k][1], restored[k][1]) for k in before)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert tuple(outputs['placement:82'].location) == pose


def verify_consumed_separator_cap(controller):
    before = mesh_buffers(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote, pose = outputs['road:99'].as_pointer(), tuple(outputs['placement:83'].location)
    metadata = json.loads(outputs['separator:7']['st_py_metadata'])
    assert len(metadata['runs']) == 1 and metadata['runs'][0]['attached_caps'] == [False, True]
    changed = []
    for source in controller['st_py_sources'].objects:
        ids, widths = source.data.attributes['st_section_id'], source.data.attributes['st_width']
        for i, value in enumerate(ids.data):
            if value.value in (2, 3):
                changed.append((source, i, widths.data[i].value))
                widths.data[i].value += .8
        source.data.update()
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    rebuilt = json.loads(controller['st_py_last_update'])['rebuilt']
    assert 'separator:7' in rebuilt and 'road:99' not in rebuilt
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert outputs['road:99'].as_pointer() == remote
    assert tuple(outputs['placement:83'].location) != pose
    NetworkScene.update(controller, force=True)
    fresh = mesh_buffers(controller)
    assert edited.keys() == fresh.keys()
    assert all(np.array_equal(edited[k][0], fresh[k][0]) and np.array_equal(edited[k][1], fresh[k][1]) for k in fresh)
    for source, i, value in changed:
        source.data.attributes['st_width'].data[i].value = value
        source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0]) and np.array_equal(before[k][1], restored[k][1]) for k in before)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert tuple(outputs['placement:83'].location) == pose


def verify_terminal_mouth(controller):
    before = mesh_buffers(controller)
    assert 'road:1' not in before
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:3'].as_pointer()
    source = next(o for o in controller['st_py_sources'].objects
                  if 1 in {v.value for v in o.data.attributes['st_section_id'].data})
    vertex = max(source.data.vertices, key=lambda v: v.co.x)
    original = vertex.co.copy()
    vertex.co.x += .2
    source.data.update()
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    assert 'road:1' in edited
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert outputs['road:3'].as_pointer() == remote
    assert outputs['placement:85'].location.z > .2
    assert json.loads(outputs['road:1']['st_py_metadata'])['quality']['maximum_edge_ratio'] < 8.
    NetworkScene.update(controller, force=True)
    fresh = mesh_buffers(controller)
    assert edited.keys() == fresh.keys()
    assert all(np.array_equal(edited[k][0], fresh[k][0]) and np.array_equal(edited[k][1], fresh[k][1]) for k in fresh)
    vertex.co = original
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0]) and np.array_equal(before[k][1], restored[k][1]) for k in before)


def verify_preferred_radius(controller):
    before = mesh_buffers(controller)
    saved_data = controller['st_py_source_data']
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    identities = {key: o.as_pointer() for key, o in outputs.items()}
    fits = json.loads(outputs['junction_mesh:1']['st_py_metadata'])['fitted_corner_radii']
    assert len(fits) == 2 and all(f['actual_radius'] < f['requested_radius'] for f in fits)
    network = model.Network.from_dict(json.loads(saved_data))
    override = network.intersections[0].with_corner_radius(model.End(1, 0), 10.)
    from dataclasses import replace
    controller['st_py_source_data'] = json.dumps(replace(network, intersections=(override,)).to_dict())
    try:
        NetworkScene.update(controller)
    except ValueError as exc:
        assert 'Requested-radius' in str(exc)
    else:
        raise AssertionError('An infeasible explicit override was silently fitted')
    assert identities == {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    assert all(np.array_equal(before[k][0], mesh_buffers(controller)[k][0]) for k in before)
    controller['st_py_source_data'] = saved_data
    NetworkScene.update(controller)
    source = next(o for o in controller['st_py_sources'].objects
                  if 2 in {v.value for v in o.data.attributes['st_section_id'].data})
    ids = source.data.attributes['st_section_id']
    vertices = {vi for edge in source.data.edges if ids.data[edge.index].value == 2 for vi in edge.vertices}
    vertex = max((source.data.vertices[i] for i in vertices), key=lambda v: v.co.y)
    original = vertex.co.copy()
    vertex.co.y += 4.
    source.data.update()
    NetworkScene.update(controller)
    after = mesh_buffers(controller)
    assert not np.array_equal(after['junction_mesh:1'][0], before['junction_mesh:1'][0])
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert outputs['road:4'].as_pointer() == identities['road:4']
    assert all(f['actual_radius'] > fits[0]['actual_radius'] for f in
               json.loads(outputs['junction_mesh:1']['st_py_metadata'])['fitted_corner_radii'])
    NetworkScene.update(controller, force=True)
    fresh = mesh_buffers(controller)
    assert after.keys() == fresh.keys()
    assert all(np.array_equal(after[k][0], fresh[k][0]) and np.array_equal(after[k][1], fresh[k][1]) for k in after)
    vertex.co = original
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0]) and np.array_equal(before[k][1], restored[k][1]) for k in before)


def verify_source_support(controller):
    before = mesh_buffers(controller)
    network = SourceObjects.read(controller)
    assert network.section_index[2].selection_ends == (0,)
    assert (0., 2.5, 0.) in network.section_index[2].points, 'Original cut control was removed'
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:4'].as_pointer()
    pose = tuple(outputs['placement:85'].location)
    source = next(o for o in controller['st_py_sources'].objects
                  if 2 in {v.value for v in o.data.attributes['st_section_id'].data})
    ids = source.data.attributes['st_section_id']
    own = {vi for edge in source.data.edges if ids.data[edge.index].value == 2 for vi in edge.vertices}
    vertex = max((source.data.vertices[i] for i in own), key=lambda v: v.co.y)
    original = vertex.co.copy()
    vertex.co.y += 1.
    source.data.update()
    assert all(np.array_equal(before[k][0], mesh_buffers(controller)[k][0]) for k in before)
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    assert not np.array_equal(edited['road:2'][0], before['road:2'][0])
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert outputs['road:4'].as_pointer() == remote
    assert tuple(outputs['placement:85'].location) != pose
    assert SourceObjects.read(controller).section_index[2].selection_ends == (0,)
    NetworkScene.update(controller, force=True)
    fresh = mesh_buffers(controller)
    assert edited.keys() == fresh.keys()
    assert all(np.array_equal(edited[k][0], fresh[k][0]) and np.array_equal(edited[k][1], fresh[k][1]) for k in fresh)
    vertex.co = original
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0]) and np.array_equal(before[k][1], restored[k][1]) for k in before)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert tuple(outputs['placement:85'].location) == pose


def verify_remote_mouth_locality(controller):
    before = mesh_buffers(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    retained = {key: outputs[key].as_pointer() for key in ('road:4', 'separator:7', 'placement:84')}
    source = next(o for o in controller['st_py_sources'].objects
                  if 8 in {v.value for v in o.data.attributes['st_section_id'].data})
    ids, widths = source.data.attributes['st_section_id'], source.data.attributes['st_width']
    changed = [(i, widths.data[i].value) for i in range(len(ids.data)) if ids.data[i].value == 8]
    for i, value in changed:
        widths.data[i].value = value + .4
    source.data.update()
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    rebuilt = json.loads(controller['st_py_last_update'])['rebuilt']
    assert 'junction_mesh:2' in rebuilt
    assert not set(retained).intersection(rebuilt)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert all(outputs[key].as_pointer() == pointer for key, pointer in retained.items())
    NetworkScene.update(controller, force=True)
    fresh = mesh_buffers(controller)
    assert edited.keys() == fresh.keys()
    assert all(np.array_equal(edited[k][0], fresh[k][0]) and np.array_equal(edited[k][1], fresh[k][1]) for k in fresh)
    for i, value in changed:
        widths.data[i].value = value
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0]) and np.array_equal(before[k][1], restored[k][1]) for k in before)


def verify_attached_union(controller):
    before = mesh_buffers(controller)
    assert 'bundle_pavement:7' in before and 'road:1' not in before and 'road:2' not in before
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:99'].as_pointer(), tuple(outputs['placement:99'].location)
    initial_anchor = tuple(outputs['placement:83'].location)
    source, vertex = section_control(controller, 1)
    original = vertex.co.copy()
    vertex.co.y += .2
    source.data.update()
    NetworkScene.update(controller)
    after = mesh_buffers(controller)
    assert not np.array_equal(after['bundle_pavement:7'][0], before['bundle_pavement:7'][0])
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert outputs['road:99'].as_pointer() == remote[0]
    assert tuple(outputs['placement:99'].location) == remote[1]
    assert tuple(outputs['placement:83'].location) != initial_anchor
    NetworkScene.update(controller, force=True)
    full = mesh_buffers(controller)
    assert after.keys() == full.keys()
    assert all(np.array_equal(after[k][0], full[k][0]) and np.array_equal(after[k][1], full[k][1]) for k in after)
    vertex.co = original
    source.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert restored.keys() == before.keys()
    assert all(np.array_equal(restored[k][0], before[k][0]) and np.array_equal(restored[k][1], before[k][1]) for k in before)
    # Change both middle controls to separate the carriageways physically.
    restore = []
    for sid, value in ((1, 5.), (2, -5.)):
        obj, control = section_control(controller, sid)
        restore.append((obj, control, control.co.copy()))
        control.co.y = value
        obj.data.update()
    NetworkScene.update(controller)
    grass = next(o for o in controller['st_py_outputs'].objects if o.get('st_py_key') == 'separator:7')
    assert json.loads(grass['st_py_metadata'])['islands'] == 1
    separated = mesh_buffers(controller)
    NetworkScene.update(controller, force=True)
    full = mesh_buffers(controller)
    assert all(np.array_equal(separated[k][0], full[k][0]) and np.array_equal(separated[k][1], full[k][1]) for k in separated)
    for obj, control, value in restore:
        control.co = value
        obj.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert all(np.array_equal(restored[k][0], before[k][0]) and np.array_equal(restored[k][1], before[k][1]) for k in before)
    identities = {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(controller)
    assert identities == {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    cold = mesh_buffers(controller)
    assert all(np.array_equal(cold[k][0], before[k][0]) and np.array_equal(cold[k][1], before[k][1]) for k in before)


def verify_multisection_union(controller):
    before = mesh_buffers(controller)
    assert 'bundle_pavement:7' in before
    assert not any(f'road:{sid}' in before for sid in (1, 11, 12, 2, 21, 22))
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:99'].as_pointer(), tuple(outputs['placement:99'].location)
    anchor = tuple(outputs['placement:83'].location)

    def equal(first, second):
        assert first.keys() == second.keys()
        assert all(np.array_equal(first[k][0], second[k][0]) and np.array_equal(first[k][1], second[k][1]) for k in first)

    source, control = section_control(controller, 11)
    original = control.co.copy()
    control.co.y += .2
    source.data.update()
    equal(before, mesh_buffers(controller))  # No continuous edit observer.
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    assert not np.array_equal(before['bundle_pavement:7'][0], edited['bundle_pavement:7'][0])
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert tuple(outputs['placement:83'].location) != anchor
    assert outputs['road:99'].as_pointer() == remote[0]
    assert tuple(outputs['placement:99'].location) == remote[1]
    NetworkScene.update(controller, force=True)
    equal(edited, mesh_buffers(controller))
    control.co = original
    source.data.update()
    NetworkScene.update(controller)
    equal(before, mesh_buffers(controller))

    saved = []
    for obj in controller['st_py_sources'].objects:
        ids, widths = obj.data.attributes['st_section_id'], obj.data.attributes['st_width']
        for index, sid in enumerate(ids.data):
            if sid.value in (1, 11, 12, 2, 21, 22):
                saved.append((obj, index, widths.data[index].value))
                widths.data[index].value = 6.2
        obj.data.update()
    NetworkScene.update(controller)
    wider = mesh_buffers(controller)
    assert not np.array_equal(before['bundle_pavement:7'][0], wider['bundle_pavement:7'][0])
    NetworkScene.update(controller, force=True)
    equal(wider, mesh_buffers(controller))
    for obj, index, value in saved:
        obj.data.attributes['st_width'].data[index].value = value
        obj.data.update()
    NetworkScene.update(controller)
    equal(before, mesh_buffers(controller))
    identities = {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(controller)
    equal(before, mesh_buffers(controller))
    assert identities == {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}


def verify_skew_bundle(controller):
    before = mesh_buffers(controller)
    pavement = next(o for o in controller['st_py_outputs'].objects if o.get('st_py_key') == 'bundle_pavement:7')
    metadata = json.loads(pavement['st_py_metadata'])
    assert metadata['finite_caps'] == 4
    assert metadata['entrance_longitudinal_slope_limit'] == .03
    assert metadata['maximum_entrance_longitudinal_slope_bound'] <= .03+1e-8
    assert len(metadata['cap_entrance_longitudinal_slope_bounds']) == 4
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:99'].as_pointer(), tuple(outputs['placement:99'].location)
    anchor = tuple(outputs['placement:84'].location)
    selected = None
    for obj in controller['st_py_sources'].objects:
        vertices = {v for edge in obj.data.edges if obj.data.attributes['st_section_id'].data[edge.index].value == 1
                    for v in edge.vertices}
        if vertices:
            index = min(vertices, key=lambda i: abs(obj.data.vertices[i].co.x+28.))
            selected = obj, obj.data.vertices[index]
            break
    assert selected is not None
    source, control = selected
    original = control.co.copy()

    def equal(first, second):
        assert first.keys() == second.keys()
        assert all(np.array_equal(first[k][0],second[k][0]) and np.array_equal(first[k][1],second[k][1]) for k in first)

    control.co.y += .05
    source.data.update()
    equal(before,mesh_buffers(controller))
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    assert not np.array_equal(edited['bundle_pavement:7'][0],before['bundle_pavement:7'][0])
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert tuple(outputs['placement:84'].location) != anchor
    assert outputs['road:99'].as_pointer() == remote[0]
    assert tuple(outputs['placement:99'].location) == remote[1]
    NetworkScene.update(controller,force=True)
    equal(edited,mesh_buffers(controller))
    control.co = original
    source.data.update()
    NetworkScene.update(controller)
    equal(before,mesh_buffers(controller))
    identities = {o['st_py_key']:o.as_pointer() for o in controller['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(controller)
    equal(before,mesh_buffers(controller))
    assert identities == {o['st_py_key']:o.as_pointer() for o in controller['st_py_outputs'].objects}


def verify_dimensional_bundle(controller):
    """Native attribute/geometry edits, composite ownership and cold rebuild."""
    before = mesh_buffers(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    metadata = json.loads(outputs['bundle_pavement:7']['st_py_metadata'])
    attached = controller['st_py_network_uuid'] == 'attached_dimensional_bundle'
    assert len(metadata['dimensional_transition_owners']) == (4 if attached else 2)
    assert not any(key.startswith('seam_mesh:') for key in before)
    remote = outputs['road:99'].as_pointer(), tuple(outputs['placement:99'].location)
    anchor_ids = (85,86) if attached else (81,82)

    def poses():
        return {o['st_py_key']: tuple(o.location) for o in controller['st_py_outputs'].objects
                if o.get('st_py_key') in {f'placement:{i}' for i in anchor_ids}}

    anchors = poses()

    def equal(a,b):
        assert a.keys() == b.keys()
        assert all(np.array_equal(a[k][0],b[k][0]) and np.array_equal(a[k][1],b[k][1]) for k in a)

    values = []
    selected = (11,21) if attached else (2,4)
    for source in controller['st_py_sources'].objects:
        ids, widths = source.data.attributes['st_section_id'], source.data.attributes['st_width']
        for i, value in enumerate(ids.data):
            if value.value in selected:
                values.append((source,i,widths.data[i].value))
                widths.data[i].value += .1
        source.data.update()
    equal(before,mesh_buffers(controller))
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    assert not np.array_equal(before['bundle_pavement:7'][0],edited['bundle_pavement:7'][0])
    assert poses() != anchors
    assert not any(key.startswith('seam_mesh:') for key in edited)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert outputs['road:99'].as_pointer() == remote[0]
    assert tuple(outputs['placement:99'].location) == remote[1]
    NetworkScene.update(controller,force=True)
    equal(edited,mesh_buffers(controller))
    for source,i,value in values:
        source.data.attributes['st_width'].data[i].value = value
        source.data.update()
    NetworkScene.update(controller)
    equal(before,mesh_buffers(controller))
    assert poses() == anchors

    source,control = section_control(controller,11 if attached else 1)
    original = control.co.copy()
    control.co.y += .05
    source.data.update()
    equal(before,mesh_buffers(controller))
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    assert not np.array_equal(before['bundle_pavement:7'][0],edited['bundle_pavement:7'][0])
    NetworkScene.update(controller,force=True)
    equal(edited,mesh_buffers(controller))
    control.co = original
    source.data.update()
    NetworkScene.update(controller)
    equal(before,mesh_buffers(controller))
    identities = {o['st_py_key']:o.as_pointer() for o in controller['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(controller)
    equal(before,mesh_buffers(controller))
    assert identities == {o['st_py_key']:o.as_pointer() for o in controller['st_py_outputs'].objects}


def verify_junction_owned_transition(controller):
    """New junction ownership survives width edits, reopening and release."""
    from check_python_street_reopen import equal
    before = mesh_buffers(controller)
    assert 'road:1' not in before and 'seam_mesh:6' not in before
    outputs = {o['st_py_key']:o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:99'].as_pointer()
    positions = {k:tuple(o.location) for k,o in outputs.items() if k.startswith('placement:')}
    widths = []
    for source in controller['st_py_sources'].objects:
        ids, values = source.data.attributes['st_section_id'], source.data.attributes['st_width']
        for edge in source.data.edges:
            if ids.data[edge.index].value == 2:
                widths.append((source,edge.index,values.data[edge.index].value))
                values.data[edge.index].value += .05
        source.data.update()
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    assert not equal(before,edited)
    now = {o['st_py_key']:o for o in controller['st_py_outputs'].objects}
    assert now['road:99'].as_pointer() == remote
    assert tuple(now['placement:82'].location) != positions['placement:82']
    NetworkScene.update(controller,force=True)
    assert equal(edited,mesh_buffers(controller))
    for source,index,value in widths:
        source.data.attributes['st_width'].data[index].value = value
        source.data.update()
    NetworkScene.update(controller)
    assert equal(before,mesh_buffers(controller))
    restore = []
    for source in controller['st_py_sources'].objects:
        ids=source.data.attributes['st_section_id']
        vertices={v for e in source.data.edges if ids.data[e.index].value in (1,2) for v in e.vertices}
        for vi in vertices:
            vertex=source.data.vertices[vi]
            if abs(vertex.co.y)<1e-6 and (abs(vertex.co.x-1.)<1e-6 or abs(vertex.co.x-2.)<1e-6):
                restore.append((source,vi,vertex.co.copy()))
                vertex.co.x=6. if abs(vertex.co.x-1.)<1e-6 else 12.
        source.data.update()
    NetworkScene.update(controller)
    assert 'seam_mesh:6' in mesh_buffers(controller)
    separate=mesh_buffers(controller)
    NetworkScene.update(controller,force=True)
    assert equal(separate,mesh_buffers(controller))
    for source,vi,point in restore:
        source.data.vertices[vi].co=point
        source.data.update()
    NetworkScene.update(controller)
    assert equal(before,mesh_buffers(controller))
    clear_runtime()
    ids={o['st_py_key']:o.as_pointer() for o in controller['st_py_outputs'].objects}
    NetworkScene.update(controller)
    assert equal(before,mesh_buffers(controller))
    assert ids == {o['st_py_key']:o.as_pointer() for o in controller['st_py_outputs'].objects}


def verify_junction_owned_profile(controller):
    """Actual rise edits retain ownership, unrelated output and anchored furniture."""
    from check_python_street_reopen import equal
    verify_junction_owned_transition(controller)
    before = mesh_buffers(controller)
    objects = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = objects['road:99'].as_pointer()
    anchor = tuple(objects['placement:82'].location)
    values = []
    selected = None
    for source in controller['st_py_sources'].objects:
        ids, rises = source.data.attributes['st_section_id'], source.data.attributes['st_rise']
        for edge in source.data.edges:
            if ids.data[edge.index].value == 2:
                values.append((source, edge.index, rises.data[edge.index].value))
                rises.data[edge.index].value += .01
                selected = source
        source.data.update()
    assert values and selected is not None
    assert equal(before, mesh_buffers(controller)), 'Rise edit bypassed explicit Update'
    bpy.context.view_layer.objects.active = selected
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    changed = mesh_buffers(controller)
    assert not equal(before, changed)
    objects = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert objects['road:99'].as_pointer() == remote
    assert tuple(objects['placement:82'].location) != anchor
    assert 'seam_mesh:6' not in changed and 'road:1' not in changed
    NetworkScene.update(controller, force=True)
    assert equal(changed, mesh_buffers(controller))
    for source, index, value in values:
        source.data.attributes['st_rise'].data[index].value = value
        source.data.update()
    NetworkScene.update(controller)
    assert equal(before, mesh_buffers(controller))


def verify_compound_owned_seam(controller, *, consumed_sections=(1, 6), physical_exterior=True,
                              claimed_sections=None):
    """Complete shared entrances survive button edits and cold reconstruction."""
    from check_python_street_reopen import equal
    before = mesh_buffers(controller)
    objects = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    poses = {key: tuple(objects[key].location) for key in ('placement:81', 'placement:82')}
    identities = {key: o.as_pointer() for key, o in objects.items()}
    metadata = json.loads(objects['junction_mesh:1']['st_py_metadata'])
    assert metadata['physical_exterior_ports'] == physical_exterior
    assert metadata['member_junctions'] == ['intersection:1', 'intersection:2']
    assert metadata['consumed_sections'] == list(consumed_sections if claimed_sections is None else claimed_sections)
    assert not {*(f'road:{sid}' for sid in consumed_sections), 'seam_mesh:6', 'junction_mesh:2'} & before.keys()
    clear_runtime()
    NetworkScene.update(controller)
    assert equal(before, mesh_buffers(controller))
    assert identities == {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    attributes = []
    source = control = None
    for obj in controller['st_py_sources'].objects:
        ids = obj.data.attributes['st_section_id']
        vertices = set()
        for edge in obj.data.edges:
            if ids.data[edge.index].value == 4:
                vertices.update(edge.vertices)
                for name, delta in (('st_width', .1), ('st_rise', .01)):
                    item = obj.data.attributes[name].data[edge.index]
                    attributes.append((obj, name, edge.index, item.value))
                    item.value += delta
        if vertices:
            source = obj
            control = max((obj.data.vertices[i] for i in vertices), key=lambda v: v.co.x)
        obj.data.update()
    assert source is not None and attributes and control is not None
    assert equal(before, mesh_buffers(controller)), 'A parameter edit bypassed explicit Update'
    remote = next(o for o in controller['st_py_outputs'].objects if o['st_py_key'] == 'road:99').as_pointer()
    bpy.context.view_layer.objects.active = source
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    changed = mesh_buffers(controller)
    assert not equal(before, changed)
    objects = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert objects['road:99'].as_pointer() == remote
    assert any(tuple(objects[key].location) != pose for key, pose in poses.items())
    NetworkScene.update(controller, force=True)
    assert equal(changed, mesh_buffers(controller))
    original_control = control.co.copy()
    control.co.y += .05
    source.data.update()
    assert equal(changed, mesh_buffers(controller)), 'A source edit bypassed explicit Update'
    remote = next(o for o in controller['st_py_outputs'].objects if o['st_py_key'] == 'road:99').as_pointer()
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    moved = mesh_buffers(controller)
    assert not equal(changed, moved)
    assert next(o for o in controller['st_py_outputs'].objects if o['st_py_key'] == 'road:99').as_pointer() == remote
    NetworkScene.update(controller, force=True)
    assert equal(moved, mesh_buffers(controller))
    control.co = original_control
    for obj, name, index, value in attributes:
        obj.data.attributes[name].data[index].value = value
        obj.data.update()
    NetworkScene.update(controller)
    assert equal(before, mesh_buffers(controller))
    objects = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert all(tuple(objects[key].location) == pose for key, pose in poses.items())


def verify_internal_compound_crown(controller, *, section_id=1, consumed_sections=(1,),
                                   physical_exterior=False, claimed_sections=None):
    """Editing an internal crown retains geometry; exposed crowns still update."""
    from check_python_street_reopen import equal
    verify_compound_owned_seam(controller, consumed_sections=consumed_sections,
                              physical_exterior=physical_exterior, claimed_sections=claimed_sections)
    before = mesh_buffers(controller)
    objects = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    identities = {key: obj.as_pointer() for key, obj in objects.items()}
    poses = {key: tuple(objects[key].location) for key in ('placement:81', 'placement:82')}
    values = []
    selected = None
    for source in controller['st_py_sources'].objects:
        ids, rises = source.data.attributes['st_section_id'], source.data.attributes['st_rise']
        for edge in source.data.edges:
            if ids.data[edge.index].value == section_id:
                values.append((source, edge.index, rises.data[edge.index].value))
                rises.data[edge.index].value += .1
                selected = source
        source.data.update()
    assert values and selected is not None
    bpy.context.view_layer.objects.active = selected
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    assert equal(before, mesh_buffers(controller))
    updated = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert identities.keys() == updated.keys()
    # The affected compound/furniture tasks may be conservatively rebuilt.
    # Exterior roads are unchanged dependencies and must retain their objects.
    assert all(updated[key].as_pointer() == pointer for key, pointer in identities.items() if key.startswith('road:'))
    assert all(tuple(updated[key].location) == pose for key, pose in poses.items())
    for source, index, value in values:
        assert source.data.attributes['st_rise'].data[index].value > value + .09
    NetworkScene.update(controller, force=True)
    assert equal(before, mesh_buffers(controller))
    for source, index, value in values:
        source.data.attributes['st_rise'].data[index].value = value
        source.data.update()
    NetworkScene.update(controller)
    assert equal(before, mesh_buffers(controller))


def verify_one_sided_free_cap(controller):
    verify_dimensional_bundle(controller)
    before=mesh_buffers(controller)
    outputs={o['st_py_key']:o for o in controller['st_py_outputs'].objects}
    metadata=json.loads(outputs['bundle_pavement:7']['st_py_metadata'])
    assert metadata['finite_caps']==2 and metadata['merged_free_end']==1
    assert metadata['finite_strip_ends']==[0]
    assert metadata['maximum_entrance_longitudinal_slope_bound']<=.03+1e-8
    anchor=tuple(outputs['placement:84'].location)
    remote=outputs['road:99'].as_pointer()
    source=control=None
    for obj in controller['st_py_sources'].objects:
        ids=obj.data.attributes['st_section_id']
        vertices={v for edge in obj.data.edges if ids.data[edge.index].value==1 for v in edge.vertices}
        if vertices:
            source=obj
            control=min((obj.data.vertices[i] for i in vertices),key=lambda v:abs(v.co.x+30.))
            break
    assert control is not None
    original=control.co.copy()
    control.co.y+=.05
    source.data.update()
    assert all(np.array_equal(before[k][0],mesh_buffers(controller)[k][0]) for k in before)
    NetworkScene.update(controller)
    edited=mesh_buffers(controller)
    outputs={o['st_py_key']:o for o in controller['st_py_outputs'].objects}
    assert tuple(outputs['placement:84'].location)!=anchor
    assert outputs['road:99'].as_pointer()==remote
    NetworkScene.update(controller,force=True)
    full=mesh_buffers(controller)
    assert all(np.array_equal(edited[k][0],full[k][0]) and np.array_equal(edited[k][1],full[k][1]) for k in edited)
    control.co=original
    source.data.update()
    NetworkScene.update(controller)
    restored=mesh_buffers(controller)
    assert all(np.array_equal(before[k][0],restored[k][0]) and np.array_equal(before[k][1],restored[k][1]) for k in before)


def verify_nonparallel_terminal(controller):
    from check_python_street_reopen import equal
    before = mesh_buffers(controller)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    metadata = json.loads(outputs['bundle_pavement:7']['st_py_metadata'])
    assert metadata['merged_terminal_ends'] == [1] and metadata['finite_caps'] == 3
    assert metadata['maximum_entrance_longitudinal_slope_bound'] <= .03+1e-8
    remote = outputs['road:99'].as_pointer()
    anchor = tuple(outputs['placement:85'].location)
    source = control = None
    for obj in controller['st_py_sources'].objects:
        ids = obj.data.attributes['st_section_id']
        vertices = {v for edge in obj.data.edges if ids.data[edge.index].value == 2 for v in edge.vertices}
        if vertices:
            source = obj
            control = min((obj.data.vertices[i] for i in vertices), key=lambda v: abs(v.co.x-30.))
            break
    assert control is not None
    original = control.co.copy()
    for y in (original.y+.05, 4.0001, 4.001, 3.99, 4.02, 4.):
        # A forced rebuild below intentionally replaces all generated objects.
        # Compare identity across this incremental Update, not across that
        # separate full-rebuild diagnostic.
        remote = next(o for o in controller['st_py_outputs'].objects if o['st_py_key']=='road:99').as_pointer()
        control.co.y = y
        source.data.update()
        previous = mesh_buffers(controller)
        NetworkScene.update(controller)
        edited = mesh_buffers(controller)
        assert not equal(previous, edited)
        outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
        assert outputs['road:99'].as_pointer() == remote
        if y == 4.:
            metadata = json.loads(outputs['bundle_pavement:7']['st_py_metadata'])
            assert 'merged_terminal_caps' not in metadata and metadata['merged_free_end'] == 1
        else:
            assert tuple(outputs['placement:85'].location) != anchor
        NetworkScene.update(controller, force=True)
        assert equal(edited, mesh_buffers(controller))
    control.co = original
    source.data.update()
    NetworkScene.update(controller)
    assert equal(before, mesh_buffers(controller))
    identities = {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(controller)
    assert equal(before, mesh_buffers(controller))
    assert identities == {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}


def verify_extended_fixture(controller):
    """Position/width edit, local/full equality, restore and cache identity."""
    from check_python_street_reopen import equal
    before = mesh_buffers(controller)
    network = SourceObjects.read(controller)
    target = np.asarray(network.section_index[1].points[1])
    source = control = None
    widths = []
    for obj in controller['st_py_sources'].objects:
        ids = obj.data.attributes['st_section_id']
        edges = [edge for edge in obj.data.edges if ids.data[edge.index].value == 1]
        if edges:
            source = obj
            vertices = {v for edge in edges for v in edge.vertices}
            control = min((obj.data.vertices[i] for i in vertices), key=lambda v: np.linalg.norm(np.asarray(v.co)-target))
            widths = [(edge.index, obj.data.attributes['st_width'].data[edge.index].value) for edge in edges]
            break
    assert source is not None and control is not None
    original = control.co.copy()
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:99'].as_pointer() if 'road:99' in outputs else None
    pose = tuple(outputs['placement:85'].location) if 'placement:85' in outputs else None
    control.co.y += .02
    for edge, value in widths:
        source.data.attributes['st_width'].data[edge].value = value+.05
    source.data.update()
    assert equal(before, mesh_buffers(controller)), 'Mesh changed before explicit Update'
    NetworkScene.update(controller)
    edited = mesh_buffers(controller)
    assert not equal(before, edited)
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    if remote is not None:
        assert outputs['road:99'].as_pointer() == remote
        assert tuple(outputs['placement:85'].location) != pose
    NetworkScene.update(controller, force=True)
    assert equal(edited, mesh_buffers(controller))
    control.co = original
    for edge, value in widths:
        source.data.attributes['st_width'].data[edge].value = value
    source.data.update()
    NetworkScene.update(controller)
    assert equal(before, mesh_buffers(controller))
    identities = {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(controller)
    assert equal(before, mesh_buffers(controller))
    assert identities == {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}


def verify_exact_bent_seam(controller):
    """Switch the original bent Bundle to Exact, edit, restore and return to Smooth."""
    from check_python_street_reopen import equal
    before = mesh_buffers(controller)
    network = SourceObjects.read(controller)
    connections = (2,6) if network.bundles else (8,)
    modes = []
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:99'].as_pointer()
    for source in controller['st_py_sources'].objects:
        ids = source.data.attributes['st_section_id']
        values = source.data.attributes['st_interpolation']
        for edge in source.data.edges:
            if ids.data[edge.index].value in (1,2,3,4,5,6):
                modes.append((source,edge.index,values.data[edge.index].value))
                values.data[edge.index].value = 0
        source.data.update()
    update = NetworkScene.update(controller)
    for connection in connections:
        assert update.products[model.Key('seam_layout',connection)].transition.certification['exact_miter']
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    assert outputs['road:99'].as_pointer() == remote
    assert all(np.isfinite(tuple(outputs['placement:87' if network.bundles else 'placement:85'].location)))
    if network.bundles:
        meta = json.loads(outputs['bundle_pavement:7']['st_py_metadata'])
        assert meta['merged_terminal_ends'] == [0,1]
    verify_extended_fixture(controller)
    for source,edge,value in modes:
        source.data.attributes['st_interpolation'].data[edge].value = value
    for source in controller['st_py_sources'].objects:
        source.data.update()
    NetworkScene.update(controller)
    assert equal(before,mesh_buffers(controller))


def verify_fitted_transition(controller, *, smooth=False):
    """Edit a finite sharp seam, compare local/full output and restore its cache."""
    from check_python_street_reopen import equal, furniture_poses
    before=mesh_buffers(controller);poses=furniture_poses(controller)
    outputs={o['st_py_key']:o for o in controller['st_py_outputs'].objects}
    remote=outputs['road:99'].as_pointer();widths=[];moved=None
    for source in controller['st_py_sources'].objects:
        ids=source.data.attributes['st_section_id'];values=source.data.attributes['st_width']
        endpoints=set()
        for edge in source.data.edges:
            sid=ids.data[edge.index].value
            if sid==1:
                widths.append((source,edge.index,values.data[edge.index].value))
                values.data[edge.index].value+=.02
            if sid==2:endpoints.update(edge.vertices)
        if endpoints:
            vertex=max((source.data.vertices[i] for i in endpoints),key=lambda v:v.co.length)
            moved=(source,vertex.index,vertex.co.copy());vertex.co.x+=.02
        source.data.update()
    assert moved is not None and widths
    assert equal(before,mesh_buffers(controller)), 'Source edit updated mesh before Update'
    result=NetworkScene.update(controller)
    plan=result.products[model.Key('seam_layout',8)].transition
    if smooth:
        assert 'exact_miter' not in plan.certification
        assert any(c.curvature_intervals>1 for c in plan.curbs)
        assert all(c.curvature_bound<=.5+1e-8 for c in plan.curbs)
    else:
        assert plan.certification['sampling_policy']=='certified_curve_error'
        assert plan.certification['certified_curve_error']<=.005
        assert plan.certification['extended_free_ends']
    outputs={o['st_py_key']:o for o in controller['st_py_outputs'].objects}
    assert outputs['road:99'].as_pointer()==remote
    local=mesh_buffers(controller);changed_poses=furniture_poses(controller)
    assert not equal(before,local) and changed_poses['placement:85']!=poses['placement:85']
    NetworkScene.update(controller,force=True)
    assert equal(local,mesh_buffers(controller)) and changed_poses==furniture_poses(controller)
    for source,index,value in widths:source.data.attributes['st_width'].data[index].value=value
    source,index,value=moved;source.data.vertices[index].co=value
    for source in controller['st_py_sources'].objects:source.data.update()
    NetworkScene.update(controller)
    assert equal(before,mesh_buffers(controller)) and poses==furniture_poses(controller)
    identities={o['st_py_key']:o.as_pointer() for o in controller['st_py_outputs'].objects}
    clear_runtime();NetworkScene.update(controller)
    assert identities=={o['st_py_key']:o.as_pointer() for o in controller['st_py_outputs'].objects}
    assert equal(before,mesh_buffers(controller))


def verify_angular_bundle(controller, *, vertex_section=29, width_sections=(17,)):
    """Actual button edits, local/full buffers and poses, then exact recovery."""
    def state():
        meshes = mesh_buffers(controller)
        poses = {o['st_py_key']: (tuple(o.location), tuple(o.rotation_euler))
                 for o in controller['st_py_outputs'].objects if o.type == 'EMPTY'}
        return meshes, poses

    def equal(first, second):
        a, poses_a = first
        b, poses_b = second
        assert a.keys() == b.keys() and poses_a == poses_b
        assert all(np.array_equal(a[k][0], b[k][0]) and np.array_equal(a[k][1], b[k][1]) for k in a)

    def update_button(source):
        bpy.context.view_layer.objects.active = source
        source.select_set(True)
        assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}

    before = state()
    outputs = {o['st_py_key']: o for o in controller['st_py_outputs'].objects}
    remote = outputs['road:99'].as_pointer()
    source, vertex = section_control(controller, vertex_section)
    position = vertex.co.copy()
    vertex.co.x += .1
    source.data.update()
    equal(before, state())  # Input changes do not publish until Update.
    update_button(source)
    edited = state()
    assert not np.array_equal(before[0]['bundle_pavement:13'][0], edited[0]['bundle_pavement:13'][0])
    rebuilt = json.loads(controller['st_py_last_update'])['rebuilt']
    assert {'bundle_pavement:13', 'placement:901'} <= set(rebuilt)
    assert 'road:99' not in rebuilt and 'placement:999' not in rebuilt
    assert remote == next(o for o in controller['st_py_outputs'].objects if o['st_py_key'] == 'road:99').as_pointer()
    NetworkScene.update(controller, force=True)
    equal(edited, state())
    vertex.co = position
    source.data.update()
    update_button(source)
    equal(before, state())

    widths = []
    for obj in controller['st_py_sources'].objects:
        ids, values = obj.data.attributes['st_section_id'], obj.data.attributes['st_width']
        for edge in obj.data.edges:
            if ids.data[edge.index].value in width_sections:
                widths.append((obj, edge.index, values.data[edge.index].value))
                values.data[edge.index].value += .1
        obj.data.update()
    assert widths
    equal(before, state())
    update_button(widths[0][0])
    widened = state()
    assert not np.array_equal(before[0]['bundle_pavement:13'][0], widened[0]['bundle_pavement:13'][0])
    assert before[1]['placement:901'] != widened[1]['placement:901']
    NetworkScene.update(controller, force=True)
    equal(widened, state())
    for obj, index, value in widths:
        obj.data.attributes['st_width'].data[index].value = value
        obj.data.update()
    update_button(source)
    equal(before, state())
    identities = {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    clear_runtime()
    update_button(source)
    equal(before, state())
    assert identities == {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    return dict(actual_update_operator=True, explicit_publication=True, vertex_update=True,
                width_update=True, exact_local_full_mesh_and_furniture=True,
                exact_restoration=True, unrelated_outputs_reused=True, fresh_reconstruction=True)


def verify_angular_chain(controller):
    """Middle-Section edits retain each original source and complete ownership."""
    source = SourceObjects.read(controller)
    members = source.bundles[0].members
    assert tuple(len(m.path) for m in members) == (3,4)
    assert any(not visit.forward for member in members for visit in member.path)
    first_ids = tuple(t.section_id for t in members[0].path)
    result = verify_angular_bundle(controller, vertex_section=229, width_sections=first_ids)
    update = NetworkScene.update(controller)
    plan = update.products[model.Key('bundle_plan',13)]
    assert tuple(plan.source_ids) == tuple(t.section_id for m in members for t in m.path)
    outputs = {o['st_py_key'] for o in controller['st_py_outputs'].objects}
    assert all('road:'+str(sid) not in outputs for sid in plan.source_ids)
    result.update(explicit_section_chains=True, reversed_section_storage=True,
                  member_sections=list(map(lambda m:len(m.path), members)),
                  middle_section_edit=True, original_source_ids=list(plan.source_ids))
    return result


def verify_angular_resolution(controller):
    """Audit real float32 R2/R4 publication, zero exaggeration and restoration."""
    from blosm.way.street_network.python_backend.mesh import Mesh
    before = mesh_buffers(controller)
    settings = (controller['st_py_resolution'], controller['st_py_exaggeration'])
    records = []
    for resolution, exaggeration in ((2, 0.), (4, 10.)):
        clear_runtime()
        controller['st_py_resolution'] = resolution
        controller['st_py_exaggeration'] = exaggeration
        NetworkScene.update(controller)
        current = mesh_buffers(controller)
        quality = {}
        for key, (xyz, ids) in current.items():
            mesh = Mesh(xyz.reshape(-1, 3), ids.reshape(-1, 4))
            quality[key] = mesh.validate()
            assert quality[key]['maximum_edge_ratio'] <= 8., (key, resolution)
            assert float(mesh.xyz[:, 2].min()) >= -1e-6
            assert float(mesh.xyz[:, 2].max()) <= .08 * exaggeration + 1e-6
        poses = {o['st_py_key']: (tuple(o.location), tuple(o.rotation_euler))
                 for o in controller['st_py_outputs'].objects if o.type == 'EMPTY'}
        NetworkScene.update(controller, force=True)
        rebuilt = mesh_buffers(controller)
        assert current.keys() == rebuilt.keys()
        assert all(np.array_equal(current[k][0], rebuilt[k][0]) and
                   np.array_equal(current[k][1], rebuilt[k][1]) for k in current)
        assert poses == {o['st_py_key']: (tuple(o.location), tuple(o.rotation_euler))
                         for o in controller['st_py_outputs'].objects if o.type == 'EMPTY'}
        records.append(dict(resolution=resolution, exaggeration=exaggeration, quality=quality,
                            exact_local_full_mesh_and_furniture=True))
    controller['st_py_resolution'], controller['st_py_exaggeration'] = settings
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert before.keys() == restored.keys()
    assert all(np.array_equal(before[k][0], restored[k][0]) and
               np.array_equal(before[k][1], restored[k][1]) for k in before)
    clear_runtime()
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    assert Path(blosm.__file__).resolve().parent == ROOT
    for obj in tuple(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    pavement = bpy.data.materials.new('Python streets / pavement')
    pavement.diffuse_color = (.18, .25, .29, 1.)
    grass = bpy.data.materials.new('Python streets / grass')
    grass.diffuse_color = (.22, .4, .17, 1.)
    materials = {'pavement': pavement, 'grass': grass}
    ui.register()
    controllers, report = [], {'backend': 'python', 'completed': False, 'cases': [], 'bundle_qualified': False}
    base_networks = (*fixture_networks(model), *ownership_fixture_networks(model),
                                *height_fixture_networks(model), *transition_fixture_networks(model),
                                *clearance_fixture_networks(model), *side_junction_fixture_networks(model),
                                *interior_fixture_networks(model), *consumed_cap_fixture_networks(model),
                                *locality_fixture_networks(model), *terminal_fixture_networks(model),
                                *radius_fixture_networks(model), *source_support_fixture_networks(model),
                                *attached_overlap_fixture_networks(model), *multisection_attached_fixture_networks(model),
                                *skew_bundle_fixture_networks(model), *dimensional_bundle_fixture_networks(model),
                                *one_sided_free_cap_fixture_networks(model), *absorbed_transition_fixture_networks(model),
                                nonparallel_terminal_fixture(model), *higher_valence_fixture_networks(model),
                                two_merged_terminal_fixture(model), two_merged_terminal_fixture(model,bent_seams=True),
                                exact_bent_seam_fixture(model), sharp_exact_seam_fixture(model,length=30.),
                                sharp_exact_seam_fixture(model,length=60.), *absorbed_profile_fixture_networks(model),
                                *compound_seam_fixture_networks(model), compound_internal_crown_fixture(model),
                                *compound_route_fixture_networks(model), sharp_smooth_seam_fixture(model))
    angular_networks = angular_fixture_networks(model)
    angular_chains = angular_chain_fixture_networks(model)
    def append_fixture(network, i):
        start = time.perf_counter()
        controller = NetworkScene.create(network, name=network.uuid, grouping='grouped' if i % 2 else 'street', materials=materials)
        controller.location = ((i % 3) * 180., -(i // 3) * 180., 0.)
        if i >= len(base_networks):
            angular_index = i - len(base_networks)
            bottom = -((len(base_networks) - 1) // 3) * 180. - 240.
            controller.location = ((angular_index % 3) * 230., bottom - (angular_index // 3) * 230., 0.)
        controllers.append(controller)
        summary = json.loads(controller['st_py_last_update'])
        report['cases'].append(dict(name=network.uuid, generation_seconds=time.perf_counter() - start,
            mesh_objects=sum(o.type == 'MESH' for o in controller['st_py_outputs'].objects),
            quads=sum(len(o.data.polygons) for o in controller['st_py_outputs'].objects if o.type == 'MESH'),
            update=summary, surface_metadata={o['st_py_key']: json.loads(o['st_py_metadata'])
                for o in controller['st_py_outputs'].objects if o.type == 'MESH'}))
        print(json.dumps(dict(fixture=network.uuid, quads=report['cases'][-1]['quads'])), flush=True)
        return controller

    # Existing dependency tests intentionally start with populated caches.
    # Finish them before the angular cold-reconstruction checks clear caches.
    for i, network in enumerate(base_networks):
        append_fixture(network, i)
    compound = next(c for c in controllers if c['st_py_network_uuid'] == 'compound_divided_cross')
    verify_compound(compound)
    report['compound_update_verified'] = True
    projected = next(c for c in controllers if c['st_py_network_uuid'] == 'profiled_corridor_crossing')
    verify_projected(projected)
    report['projected_surface_update_verified'] = True
    transition = next(c for c in controllers if c['st_py_network_uuid'] == 'divided_dimensional_seams')
    verify_transition(transition)
    report['transition_update_verified'] = True
    clearance = next(c for c in controllers if c['st_py_network_uuid'] == 'junction_dimensional_clearance')
    verify_clearance(clearance)
    report['junction_clearance_update_verified'] = True
    side = next(c for c in controllers if c['st_py_network_uuid'] == 'bundle_outward_side_junction')
    verify_curb_continuation(side)
    report['curb_continuation_update_verified'] = True
    interior = next(c for c in controllers if c['st_py_network_uuid'] == 'bundle_crossing_free_continuation')
    verify_interior_continuation(interior)
    report['interior_continuation_update_verified'] = True
    consumed = next(c for c in controllers if c['st_py_network_uuid'] == 'bundle_consumed_terminal_cap')
    verify_consumed_separator_cap(consumed)
    report['consumed_separator_cap_update_verified'] = True
    locality = next(c for c in controllers if c['st_py_network_uuid'] == 'bundle_remote_mouth_locality')
    verify_remote_mouth_locality(locality)
    report['remote_mouth_locality_verified'] = True
    terminal = next(c for c in controllers if c['st_py_network_uuid'] == 'terminal_profiled_stub')
    verify_terminal_mouth(terminal)
    report['terminal_mouth_update_verified'] = True
    preferred = next(c for c in controllers if c['st_py_network_uuid'] == 'preferred_radius_short_T')
    verify_preferred_radius(preferred)
    report['preferred_radius_update_verified'] = True
    recovered = next(c for c in controllers if c['st_py_network_uuid'] == 'clipped_original_support_T')
    verify_source_support(recovered)
    report['source_support_update_verified'] = True
    attached = next(c for c in controllers if c['st_py_network_uuid'] == 'attached_curved_overlap')
    verify_attached_union(attached)
    report['attached_overlap_update_verified'] = True
    for identity in ('attached_section_chain', 'attached_consumed_section_chain'):
        verify_multisection_union(next(c for c in controllers if c['st_py_network_uuid'] == identity))
    report['multisection_bundle_updates_verified'] = True
    for identity in ('skew_attached_overlap','skew_free_overlap'):
        verify_skew_bundle(next(c for c in controllers if c['st_py_network_uuid'] == identity))
    report['skew_bundle_updates_verified'] = True
    for identity in ('dimensional_bundle','bent_dimensional_bundle','attached_dimensional_bundle'):
        verify_dimensional_bundle(next(c for c in controllers if c['st_py_network_uuid'] == identity))
    report['dimensional_bundle_updates_verified'] = True
    verify_one_sided_free_cap(next(c for c in controllers if c['st_py_network_uuid']=='one_sided_free_cap'))
    report['one_sided_free_cap_updates_verified']=True
    verify_nonparallel_terminal(next(c for c in controllers if c['st_py_network_uuid']=='nonparallel_merged_terminal'))
    report['nonparallel_terminal_updates_verified']=True
    for identity in ('six_way', 'seven_way', 'two_merged_terminals'):
        verify_extended_fixture(next(c for c in controllers if c['st_py_network_uuid']==identity))
    report['extended_fixture_updates_verified'] = ['six_way', 'seven_way', 'two_merged_terminals']
    verify_exact_bent_seam(next(c for c in controllers if c['st_py_network_uuid']=='two_merged_terminals_bent_seams'))
    verify_exact_bent_seam(next(c for c in controllers if c['st_py_network_uuid']=='exact_bent_dimensional_seam'))
    report['exact_bent_seam_updates_verified'] = True
    for name in ('sharp_exact_seam_30','sharp_exact_seam_60'):
        verify_fitted_transition(next(c for c in controllers if c['st_py_network_uuid']==name))
    report['fitted_exact_seam_updates_verified'] = True
    verify_fitted_transition(next(c for c in controllers if c['st_py_network_uuid']=='sharp_smooth_seam_60'), smooth=True)
    report['sharp_smooth_seam_updates_verified'] = True
    verify_junction_owned_transition(next(c for c in controllers if c['st_py_network_uuid']=='junction_owned_width_seam'))
    for identity in ('junction_owned_lower_crown', 'junction_owned_higher_crown'):
        verify_junction_owned_profile(next(c for c in controllers if c['st_py_network_uuid'] == identity))
    report['junction_owned_profile_updates_verified'] = True
    report['junction_owned_transition_updates_verified']=True
    compound_seam_names = ('compound_owned_width_seam', 'compound_owned_profile_seam')
    for identity in compound_seam_names:
        verify_compound_owned_seam(next(c for c in controllers if c['st_py_network_uuid'] == identity))
    report['compound_owned_seam_updates_verified'] = list(compound_seam_names)
    verify_internal_compound_crown(next(c for c in controllers if c['st_py_network_uuid'] == 'compound_internal_crown'))
    report['internal_compound_crown_updates_verified'] = True
    verify_compound_owned_seam(next(c for c in controllers if c['st_py_network_uuid'] == 'compound_split_same_profile'),
                              physical_exterior=False)
    for identity, claims in (('compound_split_internal_crown', (1, 6)), ('compound_short_prefix', (6,))):
        verify_internal_compound_crown(next(c for c in controllers if c['st_py_network_uuid'] == identity),
            section_id=6, consumed_sections=(1, 6), physical_exterior=True, claimed_sections=claims)
    report['compound_route_updates_verified'] = ['compound_split_same_profile', 'compound_split_internal_crown',
                                                 'compound_short_prefix']
    seamed = next(c for c in controllers if c['st_py_network_uuid'] == 'seamed_approach')
    verify_seamed_approach(seamed)
    report['seamed_approach_update_verified'] = True
    divided = next(c for c in controllers if c['st_py_network_uuid'] == 'divided_T')
    saved = mesh_buffers(divided)
    member = next(o for o in divided['st_py_sources'].objects
                  if 1 in {v.value for v in o.data.attributes['st_section_id'].data})
    original_member = member.data.vertices[1].co.copy()
    member.data.vertices[1].co.x -= .3
    member.data.update()
    NetworkScene.update(divided)
    rebuilt = json.loads(divided['st_py_last_update'])['rebuilt']
    assert 'separator:7' in rebuilt and 'road:1' in rebuilt
    assert 'road:5' not in rebuilt, 'A local member bend rebuilt the remote external road'
    assert not np.array_equal(saved['separator:7'][0], mesh_buffers(divided)['separator:7'][0])
    member.data.vertices[1].co = original_member
    member.data.update()
    NetworkScene.update(divided)
    assert all(np.array_equal(saved[k][0], mesh_buffers(divided)[k][0]) for k in saved)
    report['divided_fixture_update_verified'] = True
    free = next(c for c in controllers if c['st_py_network_uuid'] == 'Bundle / curved_separated')
    original_free = mesh_buffers(free)
    restore = []
    for source in free['st_py_sources'].objects:
        owners = {}
        ids = source.data.attributes['st_section_id']
        for edge in source.data.edges:
            for vi in edge.vertices:
                owners[vi] = ids.data[edge.index].value
        for vertex in source.data.vertices:
            restore.append((source, vertex.index, vertex.co.copy()))
            vertex.co.y += -2. if owners[vertex.index] == 1 else 2.
        source.data.update()
    NetworkScene.update(free)
    union_keys = {o['st_py_key'] for o in free['st_py_outputs'].objects}
    assert union_keys == {'bundle_pavement:7'}, 'Overlap branch retained duplicate roads or grass'
    assert json.loads(free['st_py_last_update'])['timings']['removed_objects'] == 3
    saved_union = mesh_buffers(free)
    identities = {o['st_py_key']: o.as_pointer() for o in free['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(free)
    assert identities == {o['st_py_key']: o.as_pointer() for o in free['st_py_outputs'].objects}
    assert all(np.array_equal(saved_union[k][0], mesh_buffers(free)[k][0]) for k in saved_union)
    for source, vi, value in restore:
        source.data.vertices[vi].co = value
    for source in free['st_py_sources'].objects:
        source.data.update()
    NetworkScene.update(free)
    restored_free = mesh_buffers(free)
    assert original_free.keys() == restored_free.keys()
    assert all(np.array_equal(original_free[k][0], restored_free[k][0])
               and np.array_equal(original_free[k][1], restored_free[k][1]) for k in original_free)
    report['free_bundle_branch_switch_verified'] = True
    controller = next(c for c in controllers if c['st_py_network_uuid'] == 'merging_bundle')
    sources = controller['st_py_sources']
    obj = next(o for o in sources.objects if 1 in {v.value for v in o.data.attributes['st_section_id'].data})
    # The grouped source orders Sections by their input record order.
    original = obj.data.vertices[1].co.copy()
    before = mesh_buffers(controller)
    obj.data.vertices[1].co.y += .2
    obj.data.update()
    unchanged = mesh_buffers(controller)
    assert all(np.array_equal(before[k][0], unchanged[k][0]) for k in before), 'Outputs changed before Update'
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    after = mesh_buffers(controller)
    assert any(not np.array_equal(before[k][0], after[k][0]) for k in before), 'Update did not regenerate geometry'
    report['manual_update_verified'] = True
    report['edit_rebuilt'] = json.loads(controller['st_py_last_update'])['rebuilt']
    obj.data.vertices[1].co = original
    obj.data.update()
    NetworkScene.update(controller)
    restored = mesh_buffers(controller)
    assert all(np.array_equal(before[k][0], restored[k][0]) for k in before), 'Restore was not exact'
    report['exact_restoration'] = True
    # A disconnected end copy must produce a diagnostic and preserve the last
    # successful generated buffers and manifest.
    manifest = controller['st_py_manifest']
    original = obj.data.vertices[0].co.copy()
    obj.data.vertices[0].co.y += 1.
    obj.data.update()
    try:
        NetworkScene.update(controller)
    except ValueError:
        pass
    else:
        raise AssertionError('Inconsistent endpoint copies were accepted')
    assert controller['st_py_manifest'] == manifest
    preserved = mesh_buffers(controller)
    assert all(np.array_equal(restored[k][0], preserved[k][0]) for k in restored)
    obj.data.vertices[0].co = original
    obj.data.update()
    NetworkScene.update(controller)
    report['failed_update_preserves_last_good'] = True
    # Reconstruct all runtime data from saved scene input/manifest properties.
    clear_runtime()
    NetworkScene.update(controller)
    hydrated = mesh_buffers(controller)
    assert all(np.array_equal(restored[k][0], hydrated[k][0]) for k in restored)
    report['cache_reconstruction'] = True
    angular_controllers = []
    for i, network in enumerate(angular_networks, start=len(base_networks)):
        angular_controllers.append(append_fixture(network, i))
        clear_runtime()
    assert len(angular_controllers) == len(angular_networks) == 12
    report['angular_bundle_updates'] = []
    for controller in angular_controllers:
        record = dict(network=controller['st_py_network_uuid'], **verify_angular_bundle(controller))
        report['angular_bundle_updates'].append(record)
        print(json.dumps(dict(angular_update=record)), flush=True)
        clear_runtime()
    report['angular_resolution_checks'] = verify_angular_resolution(
        next(c for c in angular_controllers if c['st_py_network_uuid'] == 'angular_180_1_3'))
    chain_controllers = []
    for i, network in enumerate(angular_chains, start=len(base_networks)+len(angular_networks)):
        chain_controllers.append(append_fixture(network,i))
        clear_runtime()
    assert len(chain_controllers) == 12
    report['angular_chain_updates'] = []
    for controller in chain_controllers:
        record = dict(network=controller['st_py_network_uuid'], **verify_angular_chain(controller))
        report['angular_chain_updates'].append(record)
        print(json.dumps(dict(angular_chain_update=record)),flush=True)
        clear_runtime()
    report['angular_chain_resolution_checks'] = verify_angular_resolution(
        next(c for c in chain_controllers if c['st_py_network_uuid']=='angular_chain_180_1_3'))
    report['gn_modifiers'] = sum(m.type == 'NODES' for o in bpy.data.objects for m in o.modifiers)
    assert report['gn_modifiers'] == 0
    bootstrap = bpy.data.texts.new('OPTIONAL - Python street editing.py')
    bootstrap.write('''# Register the Update panel from the dedicated development worktree.\nimport importlib.util, sys, bpy\nfrom pathlib import Path\np = Path(r"''' + str(ROOT / 'way/street_network/python_backend') + '''")\nname = "blosm_python_street_preview"\nif name not in sys.modules:\n    spec = importlib.util.spec_from_file_location(name, p / "__init__.py", submodule_search_locations=[str(p)])\n    module = importlib.util.module_from_spec(spec)\n    sys.modules[name] = module\n    spec.loader.exec_module(module)\nif not hasattr(bpy.types, "BLOSM_PT_python_streets"):\n    from blosm_python_street_preview import ui\n    ui.register()\n'''.replace('import importlib.util, sys', 'import importlib.util, sys, bpy'))
    demo = bpy.data.texts.new('DEMO - Exact bent dimensional seams.py')
    demo.write('''# Enable the current editing panel first. Change MODE to 3 to restore Smooth.\nimport bpy\nMODE = 0\nTARGET = 'two_merged_terminals_bent_seams'  # Or 'exact_bent_dimensional_seam'.\nif bpy.context.object and bpy.context.object.mode != 'OBJECT':\n    bpy.ops.object.mode_set(mode='OBJECT')\nc = next(o for o in bpy.data.objects if o.get('st_py_network_uuid') == TARGET)\nfor source in c['st_py_sources'].objects:\n    ids = source.data.attributes['st_section_id']\n    modes = source.data.attributes['st_interpolation']\n    for edge in source.data.edges:\n        if ids.data[edge.index].value != 99:\n            modes.data[edge.index].value = MODE\n    source.data.update()\nbpy.context.view_layer.objects.active = c\nbpy.ops.blosm.update_python_streets()\n''')
    instructions = bpy.data.texts.new('READ ME - Python streets')
    instructions.write('Angular Bundle cases appear below the earlier fixtures: angular_180/225/270, gap 0 (merged throughout) or 1 (grass island), interpolation 0 (Exact) or 3 (Smooth). Edit Section 17 or 29, then press Update. The shared crown, zero outer curbs, source-relative arrow and unrelated Street 99 are retained.\n\n')
    instructions.write('angular_chain_180/225/270 adds the same six geometry/interpolation combinations with three Sections in one member and four in the other. Middle Sections 117 and 229 are stored backwards. Edit an interior control on Section 229, then Update; the shared pavement, separator and arrow placements regenerate together. Set width consistently on every Section of a member (17/117/217 or 29/129/229/329). Smooth seams use the prepared common tangents. Exact seams are collinear source connections, with curved bends retained inside the Sections; bent Exact joins and dimensional angular seams still require their separate construction. External junction attachments remain further work.\n\n')
    instructions.write('The sharp_exact_seam_30 and sharp_exact_seam_60 regressions deliberately open in Exact mode. They demonstrate a 150-degree dimensional seam with fitted finite trims and certified curve sampling. The sharp_smooth_seam_60 fixture uses Smooth interpolation and certifies the same curb polynomial using bounded local curvature proofs. Edit a free endpoint or Section width and press Update. The other fixtures retain their saved interpolation. The 30 m Smooth and mixed-mode cases still need further construction work.\n\n')
    instructions.write('junction_owned_lower_crown and junction_owned_higher_crown retain the short Section and rise change inside the rounded junction. Their external mouth keeps Section 2\'s complete 4 cm or 12 cm crown. Edit st_rise on all Section 2 edges and press Update; the junction and near anchors follow. Both have equal 6 m source widths, so this is a rise-only transition.\n\n')
    instructions.write('compound_owned_width_seam and compound_owned_profile_seam combine two T-junctions 6 m apart into one rounded physical intersection. The original connecting Sections 1/6, their Connection and both logical junctions remain source records. Section 6 changes width, and the profile fixture also changes its rise. Edit the free end, st_width or st_rise of exterior Section 4 and press Update. The complete outer profiles, intersection and both near arrow anchors reconcile together; road 99 is an unrelated control.\n\n')
    instructions.write('compound_internal_crown keeps the connecting Section 1\'s authored 18 cm crown, but that Section lies wholly inside the compound. Its four exposed entrances have 8 cm crowns, so the intersection stays below 8 cm. Changing st_rise on the internal Section and pressing Update retains the same mesh and furniture. Changing an exterior Section\'s rise or geometry still updates the surface.\n\n')
    instructions.write('Higher-valence fixtures: six_way and seven_way have unequal widths and crown heights. Edit a middle source control or width, then Update.\nThe two_merged_terminals fixture has two nonparallel merged ends around one closed separator island. Its arrows follow the two end caps; Street 99 is an unrelated control road. Internal width-changing seams are straight in this fixture. The two_merged_terminals_bent_seams and exact_bent_dimensional_seam fixtures exercise bent dimensional seams. All sources are saved in Smooth mode. Run DEMO - Exact bent dimensional seams.py to choose Exact or restore Smooth; its Update is checked against a full build. Sharp-angle and mixed-mode cases remain unfinished.\n\n')
    instructions.write('Python numerical construction checkpoint.\nSelect a wire centerline and edit its vertices.\nEnable the dedicated worktree addon or run OPTIONAL - Python street editing.py once.\nOpen the Streets sidebar and press Update.\nSection properties change width/offset/interpolation; then press Update.\nArrow empties demonstrate source-relative furniture anchors.\nThe scene covers six junction families, eleven free Bundle regression cases, an approach crossing a Section seam and a two-junction bent street.\nThe seamed_approach junction consumes Section 1 and meets the continuing Section 4; edit its middle control and press Update.\nThe divided T has a generated separator sharing rounded pavement boundaries.\nFree common-axis Bundles merge overlaps and generate separator islands.\nThe offset_reordered case intentionally places a carriageway away from its original centerline.\nThe compound_divided_cross fixture combines four junctions and retains its central island and four outer separator runs. Edit its Section 8 middle control and press Update. The compound_short_pair fixture includes an anchor on the consumed connecting Section. The merging_bundle fixture terminates its median at one shared junction. General interior Street islands, remaining local chart variants and full import qualification remain unfinished.\n')
    instructions.write('\nprofiled_corridor_crossing demonstrates bounded interior height adjustment with source-relative anchors.\nThe dimensional_seam, divided_dimensional_seams and curved_dimensional_seam fixtures join different authored Section widths through certified quintic curb transitions.\nIn divided_dimensional_seams, change Section 2 width and press Update: the transition, adjacent road, separator and furniture reconcile together.\n')
    instructions.write('\njunction_dimensional_clearance demonstrates clearance-aware trims. Change the incoming Section 3 width and press Update: the nearby width transition follows the new junction mouth while retaining a complete profiled strip and the unrelated road.\n')
    instructions.write('\nbundle_outward_side_junction continues its median along the accepted junction curb. Change Section 4 width and press Update.\nbundle_crossing_side_junctions leaves the crossing carriageway clear and produces two separate grass areas.\n')
    instructions.write('\nbundle_crossing_free_continuation retains interior Street 4. One Bundle member ends there; the continuing member and the full-profile, nonmember Street 5 bound a second separator. Change Street 5 width and press Update: its pavement, the derived separator and its furniture update together. Street 5 is not added to the Bundle.\n')
    instructions.write('\nbundle_consumed_terminal_cap keeps the original junction on Section 3 while its rounded mouth extends across the seam into Section 2. Change widths on both Sections of Street 2 together and press Update: the actual cap, separator and furniture update while the unrelated road stays unchanged.\n')
    instructions.write('\nbundle_remote_mouth_locality: change Section 8 width at the far T-junction. Its potential route through Sections 6 and 4 never reaches the Bundle. Update retains road 4, the separator and the unmoved furniture while rebuilding the far junction.\n')
    instructions.write('\nterminal_profiled_stub retains the complete short Section and full terminal profile in its junction, avoiding a separate 18 mm road strip. Extend its free endpoint by 0.2 m and press Update: a regular road span reappears and its inset furniture follows the crown. Restore the endpoint to return to the absorbed terminal cap.\n')
    instructions.write('\npreferred_radius_short_T has a preferred 10 m style radius and an 8.5 m short arm. Its two circular corners fit smaller radii while retaining the full incoming profiles and zero curbs. Extend the short arm and press Update: the fitted radii increase. Explicit radius overrides remain strict.\n')
    instructions.write('\nclipped_original_support_T keeps the old 2.5 m selection control and adds bounded authentic same-way support. Its new outer endpoint remains a selection cut. Extend that endpoint and press Update: the road and source-relative furniture update. The original Section ID, full crown and rounded zero curbs are preserved.\n')
    instructions.write('\nattached_curved_overlap joins two curved Bundle members to four rounded T-junctions. Their middle pavement overlaps and the separator forms two islands. Move middle controls from Y = 2 / -2 to 5 / -5 and press Update to open a continuous separator; restore them to return to the overlap.\n')
    instructions.write('\nattached_section_chain and attached_consumed_section_chain split the same members into Sections 1/11/12 and 2/21/22. The latter keeps its short end Sections as junction-owned source records. Edit a middle control on Section 11 and press Update.\n')
    instructions.write('\nskew_attached_overlap and skew_free_overlap move the upper control at X = -28 by 0.2 m. Finite quad strips preserve the actual mouth profiles and curbs while joining the shared Bundle body. Move that control slightly in Y and press Update. The physical entrance-longitudinal slope bound is limited to 3 percent, independently of inspection exaggeration.\n')
    instructions.write('\ndimensional_bundle and bent_dimensional_bundle widen both carriageways from 6 m to 10 m through certified Section transitions: their median closes and overlapping pavement is merged. attached_dimensional_bundle keeps rounded junction caps while its middle Sections widen from 6 m to 8 m. Edit st_width on the appropriate Section edges, or move an interior control, then press Update. Transition-owned furniture follows the merged height field.\n')
    instructions.write('\none_sided_free_cap reproduces a 5 cm move of the upper control at X = -30. Its separated end becomes skew while the far end remains merged. Two finite profile strips join that separated end to the shared Bundle body. Move the same control slightly in Y and press Update; the nearby arrow anchor follows the cap.\n')
    instructions.write('\njunction_owned_width_seam keeps the original 2 m inner Section while its rounded junction enters the slightly wider Section 2. Change Section 2 width and press Update. Both near arrows query the junction. Move the two seam copies from X = 2 to 12 (and the inner control from 1 to 6), then Update: the isolated dimensional transition reappears. Restore the controls to return ownership to the junction.\n')
    instructions.write('\nnonparallel_merged_terminal moves the Section 2 control at X = 30 to Y = 4.05. A finite union cap retains both actual terminal cuts and the exposed short curb. Move that control slightly and press Update; its nearby arrow follows the surface. Set Y = 4 to return to the aligned merged-end branch, then restore 4.05.\n')
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type == 'VIEW_3D':
                space = area.spaces.active
                space.shading.color_type = 'MATERIAL'
                rows = (len(controllers) + 2) // 3
                space.region_3d.view_distance = max(350., rows * 180.)
                space.region_3d.view_location = (180., -.5 * 180. * (rows - 1), 0.)
    bpy.context.preferences.filepaths.file_preview_type = 'NONE'
    bpy.ops.wm.save_as_mainfile(filepath=str((args.output / 'python_streets_checkpoint.blend').resolve()))
    report['scene'] = str((args.output / 'python_streets_checkpoint.blend').resolve())
    report['completed'] = True
    (args.output / 'scene_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(dict(completed=True, scene=report['scene'], cases=len(report['cases']),
        total_quads=sum(r['quads'] for r in report['cases']), gn_modifiers=report['gn_modifiers'])), flush=True)


if __name__ == '__main__':
    main()
