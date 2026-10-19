"""Bind optional native ownership validation without a second public solver."""
import bpy
from ..gn.kernel.pipeline import set_input


class SeparatorOwnershipProvider:
    result_input = 'Separator Ownership Result'
    key_input = 'Separator Ownership Key'
    marker = 'blosm_separator_ownership_provider'

    def __init__(self, solver):
        from .assets import validate_original_separator_ownership_provider
        validate_original_separator_ownership_provider(solver)
        self.solver = solver

    def install(self, controller, collection):
        from ..gn.editable_input.curved_patch_provider_bindings import input_values, set_inputs
        from .free_bundle_provider import FreeBundleProvider
        network = controller.modifiers[0].node_group
        if network.get('blosm_separator_ownership_provider_input') != 1:
            raise ValueError('The shared solver has no separator ownership provider input')
        if self.solver.get('blosm_street_asset_version') != network.get('blosm_street_asset_version'):
            raise ValueError('Ownership provider and shared network revisions must match')
        inputs = input_values(controller)
        rows = inputs['Separators']
        if (inputs['Separator Count'] <= 0 or rows is None or rows.type != 'MESH'
                or len(rows.data.vertices) != inputs['Separator Count']):
            raise ValueError('Only a consistent nonempty separator declaration needs a provider')
        if any(o.get(self.marker) and o.parent == controller for o in bpy.data.objects):
            raise ValueError('Rebind or remove the existing separator ownership provider first')
        used = {o.get('blosm_separator_ownership_key') for o in bpy.data.objects}
        key = 1
        while key in used: key += 1
        if key >= 2**24:
            raise ValueError('Separator ownership provider identity capacity exceeded')
        child = bpy.data.objects.new(controller.name+' / separator ownership', controller.data)
        collection.objects.link(child); child.parent = controller
        child.hide_render = True; child.hide_select = True; child.hide_set(True)
        child[self.marker] = 1
        mod = child.modifiers.new('Shared native separator ownership', 'NODES')
        mod.node_group = self.solver; mod.show_viewport = False
        controller['blosm_separator_ownership_key'] = key
        set_input(controller, self.key_input, key)
        set_input(controller, self.result_input, child)
        set_inputs(child, input_values(controller))
        # Reuse the native numeric-control binding. ID references are the same
        # authored tables/collections; coordinate edits require no Python timer.
        FreeBundleProvider._drivers(controller, child)
        mod.show_viewport = True
        return child
