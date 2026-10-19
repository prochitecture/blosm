"""Author one shared native provider for a network's declared free Bundles.

Only source references and numeric controls are bound here. Surface creation,
interactive source edits, validation and diagnostics are Geometry Nodes work.
"""
import bpy
from ..gn.kernel.pipeline import set_input


class FreeBundleProvider:
    def __init__(self, solver):
        from .assets import validate_original_free_bundle_provider
        validate_original_free_bundle_provider(solver)
        self.solver = solver

    def install(self, controller, collection):
        from ..gn.editable_input.curved_patch_provider_bindings import input_values, set_inputs
        network = controller.modifiers[0].node_group
        if network.get('blosm_free_bundle_provider_input') != 1:
            raise ValueError('The common solver has no free-Bundle provider input')
        if self.solver.get('blosm_street_asset_version') != network.get('blosm_street_asset_version'):
            raise ValueError('The free-Bundle provider and network revisions must match')
        if input_values(controller)['Bundle Count'] <= 0:
            raise ValueError('Only a declared free Bundle needs this provider')
        if any(o.get('blosm_free_bundle_provider') and o.parent == controller for o in bpy.data.objects):
            raise ValueError('Rebind or remove the existing free-Bundle provider before installing another')
        child = bpy.data.objects.new(controller.name+' / free Bundle provider', controller.data)
        collection.objects.link(child); child.parent = controller
        child.hide_render = True; child.hide_select = True; child.hide_set(True)
        child['blosm_free_bundle_provider'] = 1
        mod = child.modifiers.new('Shared native free Bundle construction', 'NODES')
        mod.node_group = self.solver; mod.show_viewport = False
        used = {o.get('blosm_free_bundle_provider_key') for o in bpy.data.objects}
        key = 1
        while key in used: key += 1
        if key >= 2**24:
            raise ValueError('Free-Bundle provider identity capacity exceeded')
        controller['blosm_free_bundle_provider_key'] = key
        set_input(controller, 'Free Bundle Provider Key', key)
        set_input(controller, 'Free Bundle Result', child)
        set_inputs(child, input_values(controller))
        self._drivers(controller, child)
        mod.show_viewport = True
        return child

    @staticmethod
    def _drivers(controller, child, *, exclude=()):
        source_mod, target_mod = controller.modifiers[0], child.modifiers[0]
        source = {s.name: getattr(source_mod.properties.inputs, s.identifier)
                  for s in source_mod.node_group.interface.items_tree
                  if s.item_type == 'SOCKET' and s.in_out == 'INPUT'
                  and s.socket_type in ('NodeSocketFloat', 'NodeSocketInt', 'NodeSocketBool')}
        for socket in target_mod.node_group.interface.items_tree:
            if (socket.item_type != 'SOCKET' or socket.in_out != 'INPUT' or
                    socket.name not in source or socket.name in exclude):
                continue
            prop = getattr(target_mod.properties.inputs, socket.identifier)
            curve = prop.driver_add('value'); curve.driver.type = 'AVERAGE'
            curve.extrapolation = 'LINEAR'
            for point in curve.keyframe_points: point.interpolation = 'LINEAR'
            variable = curve.driver.variables.new(); variable.name = 'control'; variable.type = 'SINGLE_PROP'
            variable.targets[0].id = controller
            variable.targets[0].data_path = source[socket.name].path_from_id('value')
