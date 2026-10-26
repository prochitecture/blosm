"""Fresh-process save/reopen, Edit Mode, width and Update verification."""
import argparse
import json
from pathlib import Path
import sys
import bpy
import bmesh
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT / 'tests'))
from check_python_street_scene import (mesh_buffers, verify_seamed_approach, verify_compound,
                                       verify_projected, verify_transition, verify_clearance, verify_curb_continuation,
                                       verify_interior_continuation, verify_consumed_separator_cap,
                                       verify_remote_mouth_locality, verify_terminal_mouth, verify_preferred_radius,
                                       verify_source_support, verify_attached_union, verify_multisection_union, verify_skew_bundle,
                                       verify_dimensional_bundle, verify_one_sided_free_cap, verify_junction_owned_transition,
                                       verify_nonparallel_terminal, verify_extended_fixture, verify_exact_bent_seam,
                                       verify_fitted_transition, verify_junction_owned_profile, verify_compound_owned_seam,
                                       verify_internal_compound_crown, verify_angular_bundle,
                                       verify_angular_resolution, verify_angular_chain)
from blosm.way.street_network.python_backend.scene import NetworkScene, clear_runtime
from blosm.way.street_network.python_backend import ui


def equal(a, b):
    return a.keys() == b.keys() and all(np.array_equal(a[k][0], b[k][0])
                                      and np.array_equal(a[k][1], b[k][1]) for k in a)


def furniture_poses(controller):
    return {o['st_py_key']: (tuple(o.location), tuple(o.rotation_euler), o.get('st_py_asset'))
            for o in controller['st_py_outputs'].objects if o.type == 'EMPTY'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()))
    ui.register()
    controllers = [o for o in bpy.data.objects if o.get('st_py_source_data')]
    angular = [c for c in controllers if c['st_py_network_uuid'].startswith('angular_')
               and not c['st_py_network_uuid'].startswith('angular_chain_')]
    chains = [c for c in controllers if c['st_py_network_uuid'].startswith('angular_chain_')]
    angular_records = []
    if angular:
        expected = {f'angular_{sweep}_{gap}_{mode}' for sweep in (180, 225, 270)
                    for gap in (0, 1) for mode in (0, 3)}
        assert {c['st_py_network_uuid'] for c in angular} == expected and len(angular) == 12
        for controller in angular:
            clear_runtime()
            angular_records.append(dict(network=controller['st_py_network_uuid'], **verify_angular_bundle(controller)))
            print(json.dumps(dict(angular_reopen=angular_records[-1])), flush=True)
        angular_resolutions = verify_angular_resolution(
            next(c for c in angular if c['st_py_network_uuid'] == 'angular_180_1_3'))
    chain_records = []
    if chains:
        expected = {f'angular_chain_{sweep}_{gap}_{mode}' for sweep in (180,225,270)
                    for gap in (0,1) for mode in (0,3)}
        assert {c['st_py_network_uuid'] for c in chains} == expected and len(chains)==12
        for controller in chains:
            clear_runtime()
            chain_records.append(dict(network=controller['st_py_network_uuid'], **verify_angular_chain(controller)))
            print(json.dumps(dict(angular_chain_reopen=chain_records[-1])),flush=True)
        chain_resolutions = verify_angular_resolution(
            next(c for c in chains if c['st_py_network_uuid']=='angular_chain_180_1_3'))
    for identity in ('dimensional_bundle','bent_dimensional_bundle','attached_dimensional_bundle'):
        verify_dimensional_bundle(next(c for c in controllers if c['st_py_network_uuid'] == identity))
    verify_one_sided_free_cap(next(c for c in controllers if c['st_py_network_uuid']=='one_sided_free_cap'))
    verify_junction_owned_transition(next(c for c in controllers if c['st_py_network_uuid']=='junction_owned_width_seam'))
    profile_names = ('junction_owned_lower_crown', 'junction_owned_higher_crown')
    profile_controllers = [c for c in controllers if c['st_py_network_uuid'] in profile_names]
    for controller in profile_controllers:
        verify_junction_owned_profile(controller)
    compound_seam_controllers = [c for c in controllers if c['st_py_network_uuid'] in
                                ('compound_owned_width_seam', 'compound_owned_profile_seam')]
    for controller in compound_seam_controllers:
        verify_compound_owned_seam(controller)
    internal_crown_controllers = [c for c in controllers if c['st_py_network_uuid'] == 'compound_internal_crown']
    for controller in internal_crown_controllers:
        verify_internal_compound_crown(controller)
    for controller in (c for c in controllers if c['st_py_network_uuid'] == 'compound_split_same_profile'):
        verify_compound_owned_seam(controller, physical_exterior=False)
    for identity, claims in (('compound_split_internal_crown', (1, 6)), ('compound_short_prefix', (6,))):
        for controller in (c for c in controllers if c['st_py_network_uuid'] == identity):
            verify_internal_compound_crown(controller, section_id=6, consumed_sections=(1, 6),
                                          physical_exterior=True, claimed_sections=claims)
    verify_nonparallel_terminal(next(c for c in controllers if c['st_py_network_uuid']=='nonparallel_merged_terminal'))
    for identity in ('six_way', 'seven_way', 'two_merged_terminals'):
        verify_extended_fixture(next(c for c in controllers if c['st_py_network_uuid']==identity))
    verify_exact_bent_seam(next(c for c in controllers if c['st_py_network_uuid']=='two_merged_terminals_bent_seams'))
    verify_exact_bent_seam(next(c for c in controllers if c['st_py_network_uuid']=='exact_bent_dimensional_seam'))
    for name in ('sharp_exact_seam_30','sharp_exact_seam_60'):
        verify_fitted_transition(next(c for c in controllers if c['st_py_network_uuid']==name))
    for controller in (c for c in controllers if c['st_py_network_uuid']=='sharp_smooth_seam_60'):
        verify_fitted_transition(controller, smooth=True)
    attached = next(c for c in controllers if c['st_py_network_uuid'] == 'attached_curved_overlap')
    verify_attached_union(attached)
    for identity in ('attached_section_chain', 'attached_consumed_section_chain'):
        verify_multisection_union(next(c for c in controllers if c['st_py_network_uuid'] == identity))
    for identity in ('skew_attached_overlap','skew_free_overlap'):
        verify_skew_bundle(next(c for c in controllers if c['st_py_network_uuid'] == identity))
    recovered = next(c for c in controllers if c['st_py_network_uuid'] == 'clipped_original_support_T')
    saved_recovered = mesh_buffers(recovered)
    recovered_ids = {o['st_py_key']: o.as_pointer() for o in recovered['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(recovered)
    assert equal(saved_recovered, mesh_buffers(recovered))
    assert recovered_ids == {o['st_py_key']: o.as_pointer() for o in recovered['st_py_outputs'].objects}
    verify_source_support(recovered)
    preferred = next(c for c in controllers if c['st_py_network_uuid'] == 'preferred_radius_short_T')
    saved_preferred = mesh_buffers(preferred)
    preferred_ids = {o['st_py_key']: o.as_pointer() for o in preferred['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(preferred)
    assert equal(saved_preferred, mesh_buffers(preferred))
    assert preferred_ids == {o['st_py_key']: o.as_pointer() for o in preferred['st_py_outputs'].objects}
    verify_preferred_radius(preferred)
    terminal = next(c for c in controllers if c['st_py_network_uuid'] == 'terminal_profiled_stub')
    saved_terminal = mesh_buffers(terminal)
    terminal_ids = {o['st_py_key']: o.as_pointer() for o in terminal['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(terminal)
    assert equal(saved_terminal, mesh_buffers(terminal))
    assert terminal_ids == {o['st_py_key']: o.as_pointer() for o in terminal['st_py_outputs'].objects}
    verify_terminal_mouth(terminal)
    locality = next(c for c in controllers if c['st_py_network_uuid'] == 'bundle_remote_mouth_locality')
    saved_locality = mesh_buffers(locality)
    locality_ids = {o['st_py_key']: o.as_pointer() for o in locality['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(locality)
    assert equal(saved_locality, mesh_buffers(locality))
    assert locality_ids == {o['st_py_key']: o.as_pointer() for o in locality['st_py_outputs'].objects}
    verify_remote_mouth_locality(locality)
    consumed = next(c for c in controllers if c['st_py_network_uuid'] == 'bundle_consumed_terminal_cap')
    saved_consumed = mesh_buffers(consumed)
    consumed_ids = {o['st_py_key']: o.as_pointer() for o in consumed['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(consumed)
    assert equal(saved_consumed, mesh_buffers(consumed))
    assert consumed_ids == {o['st_py_key']: o.as_pointer() for o in consumed['st_py_outputs'].objects}
    verify_consumed_separator_cap(consumed)
    interior = next(c for c in controllers if c['st_py_network_uuid'] == 'bundle_crossing_free_continuation')
    saved_interior = mesh_buffers(interior)
    interior_ids = {o['st_py_key']: o.as_pointer() for o in interior['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(interior)
    assert equal(saved_interior, mesh_buffers(interior))
    assert interior_ids == {o['st_py_key']: o.as_pointer() for o in interior['st_py_outputs'].objects}
    verify_interior_continuation(interior)
    side = next(c for c in controllers if c['st_py_network_uuid'] == 'bundle_outward_side_junction')
    saved_side = mesh_buffers(side)
    side_ids = {o['st_py_key']: o.as_pointer() for o in side['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(side)
    assert equal(saved_side, mesh_buffers(side))
    assert side_ids == {o['st_py_key']: o.as_pointer() for o in side['st_py_outputs'].objects}
    verify_curb_continuation(side)
    clearance = next(c for c in controllers if c['st_py_network_uuid'] == 'junction_dimensional_clearance')
    saved_clearance = mesh_buffers(clearance)
    clearance_ids = {o['st_py_key']: o.as_pointer() for o in clearance['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(clearance)
    assert equal(saved_clearance, mesh_buffers(clearance))
    assert clearance_ids == {o['st_py_key']: o.as_pointer() for o in clearance['st_py_outputs'].objects}
    verify_clearance(clearance)
    transition = next(c for c in controllers if c['st_py_network_uuid'] == 'divided_dimensional_seams')
    saved_transition = mesh_buffers(transition)
    transition_ids = {o['st_py_key']: o.as_pointer() for o in transition['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(transition)
    assert equal(saved_transition, mesh_buffers(transition))
    assert transition_ids == {o['st_py_key']: o.as_pointer() for o in transition['st_py_outputs'].objects}
    verify_transition(transition)
    projected = next(c for c in controllers if c['st_py_network_uuid'] == 'profiled_corridor_crossing')
    saved_projected = mesh_buffers(projected)
    projected_ids = {o['st_py_key']: o.as_pointer() for o in projected['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(projected)
    assert equal(saved_projected, mesh_buffers(projected))
    assert projected_ids == {o['st_py_key']: o.as_pointer() for o in projected['st_py_outputs'].objects}
    verify_projected(projected)
    compound = next(c for c in controllers if c['st_py_network_uuid'] == 'compound_divided_cross')
    saved_compound = mesh_buffers(compound)
    compound_ids = {o['st_py_key']: o.as_pointer() for o in compound['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(compound)
    assert equal(saved_compound, mesh_buffers(compound))
    assert compound_ids == {o['st_py_key']: o.as_pointer() for o in compound['st_py_outputs'].objects}
    verify_compound(compound)
    seamed = next(c for c in controllers if c['st_py_network_uuid'] == 'seamed_approach')
    saved_seamed = mesh_buffers(seamed)
    seamed_ids = {o['st_py_key']: o.as_pointer() for o in seamed['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(seamed)
    assert equal(saved_seamed, mesh_buffers(seamed))
    assert seamed_ids == {o['st_py_key']: o.as_pointer() for o in seamed['st_py_outputs'].objects}
    verify_seamed_approach(seamed)
    controller = next(c for c in controllers if c['st_py_network_uuid'] == 'bent_two_intersections')
    before = mesh_buffers(controller)
    anchors_before = furniture_poses(controller)
    assert len(anchors_before) == 2
    identities = {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(controller)
    assert equal(before, mesh_buffers(controller)), 'Fresh rehydration changed the saved output'
    assert identities == {o['st_py_key']: o.as_pointer() for o in controller['st_py_outputs'].objects}, 'Unchanged outputs were replaced after reopening'
    source = next(o for o in controller['st_py_sources'].objects
                  if 1 in {v.value for v in o.data.attributes['st_section_id'].data})
    bpy.ops.object.select_all(action='DESELECT')
    source.select_set(True)
    bpy.context.view_layer.objects.active = source
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(source.data)
    bm.verts.ensure_lookup_table()
    original = bm.verts[1].co.copy()
    bm.verts[1].co.y += .15
    bmesh.update_edit_mesh(source.data, loop_triangles=False, destructive=False)
    assert equal(before, mesh_buffers(controller)), 'Edit Mode moved generated output before Update'
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    assert not equal(before, mesh_buffers(controller)), 'Edit Mode input was not read'
    assert anchors_before != furniture_poses(controller), 'Attached furniture did not follow the source edit'
    bm.verts[1].co = original
    bmesh.update_edit_mesh(source.data, loop_triangles=False, destructive=False)
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    assert equal(before, mesh_buffers(controller)), 'Edit Mode restoration changed the original buffers'
    assert anchors_before == furniture_poses(controller)
    bpy.ops.object.mode_set(mode='OBJECT')
    # Exercise an authored attribute edit, then explicit regeneration. The
    # neighboring intersection meshes must be included in the affected closure.
    attribute = source.data.attributes['st_width']
    old = np.array([value.value for value in attribute.data])
    section_ids = [value.value for value in source.data.attributes['st_section_id'].data]
    for value, identity in zip(attribute.data, section_ids):
        if identity == 1:
            value.value = 6.2
    source.data.update()
    assert equal(before, mesh_buffers(controller))
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    rebuilt = json.loads(controller['st_py_last_update'])['rebuilt']
    assert 'junction_mesh:1' in rebuilt and 'junction_mesh:2' in rebuilt
    attribute.data.foreach_set('value', old)
    source.data.update()
    assert bpy.ops.blosm.update_python_streets() == {'FINISHED'}
    assert equal(before, mesh_buffers(controller))
    divided = next(c for c in controllers if c['st_py_network_uuid'] == 'divided_T')
    before_divided = mesh_buffers(divided)
    unchanged_ids = {o['st_py_key']: o.as_pointer() for o in divided['st_py_outputs'].objects}
    NetworkScene.update(divided)
    assert equal(before_divided, mesh_buffers(divided))
    assert unchanged_ids == {o['st_py_key']: o.as_pointer() for o in divided['st_py_outputs'].objects}
    member = next(o for o in divided['st_py_sources'].objects
                  if 1 in {v.value for v in o.data.attributes['st_section_id'].data})
    original_member = member.data.vertices[1].co.copy()
    member.data.vertices[1].co.x -= .25
    member.data.update()
    NetworkScene.update(divided)
    assert not equal(before_divided, mesh_buffers(divided))
    assert 'separator:7' in json.loads(divided['st_py_last_update'])['rebuilt']
    member.data.vertices[1].co = original_member
    member.data.update()
    NetworkScene.update(divided)
    assert equal(before_divided, mesh_buffers(divided))
    union = next(c for c in controllers if c['st_py_network_uuid'] == 'Bundle / curved_overlap')
    original_union = mesh_buffers(union)
    union_ids = {o['st_py_key']: o.as_pointer() for o in union['st_py_outputs'].objects}
    clear_runtime()
    NetworkScene.update(union)
    assert original_union.keys() == {'bundle_pavement:7'}
    assert equal(original_union, mesh_buffers(union))
    assert union_ids == {o['st_py_key']: o.as_pointer() for o in union['st_py_outputs'].objects}
    member = next(o for o in union['st_py_sources'].objects
                  if 1 in {v.value for v in o.data.attributes['st_section_id'].data})
    attribute = member.data.attributes['st_width']
    old_widths = [v.value for v in attribute.data]
    for value, sid in zip(attribute.data, member.data.attributes['st_section_id'].data):
        if sid.value == 1:
            value.value += .2
    member.data.update()
    NetworkScene.update(union)
    assert not equal(original_union, mesh_buffers(union))
    attribute.data.foreach_set('value', old_widths)
    member.data.update()
    NetworkScene.update(union)
    assert equal(original_union, mesh_buffers(union))
    report = dict(passed=True, scene=str(args.scene.resolve()), controllers=len(controllers),
        fresh_reconstruction=True, edit_mode_update=True, attributes_update=True,
        exact_restoration=True, two_junction_dependency_update=True, divided_separator_reopen_update=True,
        union_reopen_width_update=True, furniture_anchor_reopen_update=True, seamed_approach_reopen_update=True,
        compound_reopen_update=True,
        projected_surface_reopen_update=True,
        gn_modifiers=sum(m.type == 'NODES' for o in bpy.data.objects for m in o.modifiers))
    report['dimensional_transition_reopen_update'] = True
    report['junction_clearance_reopen_update'] = True
    report['curb_continuation_reopen_update'] = True
    report['interior_continuation_reopen_update'] = True
    report['consumed_separator_cap_reopen_update'] = True
    report['remote_mouth_locality_reopen_update'] = True
    report['terminal_mouth_reopen_update'] = True
    report['preferred_radius_reopen_update'] = True
    report['source_support_reopen_update'] = True
    report['attached_overlap_reopen_update'] = True
    report['multisection_union_reopen_update'] = True
    report['skew_bundle_reopen_update'] = True
    report['dimensional_bundle_reopen_update'] = True
    report['one_sided_free_cap_reopen_update'] = True
    report['junction_owned_transition_reopen_update'] = True
    report['junction_owned_profile_reopen_update'] = [c['st_py_network_uuid'] for c in profile_controllers]
    report['compound_owned_seam_reopen_update'] = [c['st_py_network_uuid'] for c in compound_seam_controllers]
    report['internal_compound_crown_reopen_update'] = [c['st_py_network_uuid'] for c in internal_crown_controllers]
    report['nonparallel_terminal_reopen_update'] = True
    report['extended_fixture_reopen_update'] = ['six_way', 'seven_way', 'two_merged_terminals']
    report['exact_bent_seam_reopen_update'] = True
    report['fitted_exact_seam_reopen_update'] = True
    report['angular_bundle_reopen_updates'] = angular_records
    if angular:
        report['angular_resolution_checks'] = angular_resolutions
    report['angular_chain_reopen_updates'] = chain_records
    if chains:
        report['angular_chain_resolution_checks'] = chain_resolutions
    (args.output / 'reopen_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
