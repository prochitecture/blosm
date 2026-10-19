"""Evaluate newly authored independent providers before their final assembly.

This schedules normal native dependency-graph evaluations. It creates no mesh
replacement or cache owned by Python and changes no solver inputs. A source
edit after import still invalidates Blender's native dependencies normally.
"""
import time
import bpy


PROVIDER_ORDER = {
    'original_patch_provider.1': 0,
    'original_strip_provider.1': 1,
    'free_bundle_provider.1': 2,
    'separator_ownership_provider.1': 3,
    'separator_fill_provider.1': 4,
    'original_intake_provider.1': 5,
    'aggregate_provider.1': 5,
}


def prepared_providers(controller):
    """Enabled native producers, ordered before their downstream consumers."""
    providers = []
    for obj in bpy.data.objects:
        if obj.parent != controller: continue
        for modifier in obj.modifiers:
            if modifier.type != 'NODES' or modifier.node_group is None: continue
            schema = modifier.node_group.get('blosm_street_input_schema')
            if schema in PROVIDER_ORDER and modifier.show_viewport:
                providers.append((PROVIDER_ORDER[schema], obj.name, obj, modifier))
    providers.sort(key=lambda row: row[:2])
    return providers


def evaluate_prepared_network(controller, *, observe=None):
    providers = prepared_providers(controller)
    # Failure labels read the final controller. Keep them from requesting an
    # intermediate assembly while the independent providers are being enabled.
    displays = [m for obj in bpy.data.objects if obj.parent == controller
                for m in obj.modifiers if m.type == 'NODES' and m.node_group is not None
                and m.node_group.get('blosm_street_input_schema') == 'original_failure_display.1'
                and m.show_viewport]
    modifier = next(m for m in controller.modifiers if m.type == 'NODES')
    original = modifier.show_viewport
    rows = []
    try:
        if modifier.show_viewport:
            modifier.show_viewport = False
        for display in displays: display.show_viewport = False
        for _, _, obj, child in providers: child.show_viewport = False
        for kind, name, obj, child in providers:
            if observe: observe(dict(phase='starting_provider', name=name, kind=kind))
            child.show_viewport = True
            started = time.perf_counter(); bpy.context.view_layer.update()
            mesh = obj.evaluated_get(bpy.context.evaluated_depsgraph_get()).data
            row = dict(name=name, kind=kind, seconds=time.perf_counter()-started,
                       points=len(mesh.vertices), faces=len(mesh.polygons))
            rows.append(row)
            if observe: observe(dict(phase='evaluated_provider', **row))
        if observe: observe(dict(phase='starting_assembly', name=controller.name))
        if modifier.show_viewport != original:
            modifier.show_viewport = original
        started = time.perf_counter(); bpy.context.view_layer.update()
        if observe: observe(dict(phase='evaluated_assembly', name=controller.name, seconds=time.perf_counter()-started))
        if displays:
            for display in displays: display.show_viewport = True
            bpy.context.view_layer.update()
    finally:
        # RNA assignments tag native geometry dirty even when the value is
        # unchanged. Do not invalidate the assembly we have just evaluated.
        # On exceptions, providers not yet reached still need restoration.
        for _, _, obj, child in providers:
            if not child.show_viewport:
                child.show_viewport = True
        if modifier.show_viewport != original:
            modifier.show_viewport = original
        for display in displays:
            if not display.show_viewport:
                display.show_viewport = True
    return rows
