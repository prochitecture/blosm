"""Bind a shared original-span producer to editable native input dependencies.

Coordinate and numeric edits use Blender dependencies and scalar drivers.
Replacing source IDs is an explicit rebind operation after topology preparation;
this class neither generates surfaces nor installs a polling handler.
"""
import bpy
from mathutils import Matrix
from ..gn.kernel.pipeline import set_input


class OriginalIntakeProvider:
    result_input = 'Original Span Result'
    key_input = 'Original Span Provider Key'
    marker = 'blosm_original_intake_provider'
    key_property = 'blosm_original_intake_key'
    scalar_types = ('NodeSocketFloat', 'NodeSocketInt', 'NodeSocketBool')

    def __init__(self, solver):
        from .assets import validate_original_intake_provider
        validate_original_intake_provider(solver)
        self.solver = solver

    @staticmethod
    def _inputs(obj):
        mod = obj.modifiers[0]
        return {s.name:(s, getattr(mod.properties.inputs, s.identifier))
                for s in mod.node_group.interface.items_tree
                if s.item_type == 'SOCKET' and s.in_out == 'INPUT'
                and s.socket_type != 'NodeSocketGeometry'}

    def _validate_controller(self, controller):
        from .assets import validate_original_intake_provider
        validate_original_intake_provider(self.solver, controller.modifiers[0].node_group)
        inputs = self._inputs(controller)
        if inputs['Prepared Network Kind'][1].value != 0:
            raise ValueError('Only original_spans.1 needs an original-span provider')
        if controller.type != 'MESH':
            raise ValueError('The original-span input carrier must be a mesh')
        return inputs

    def install(self, controller, collection):
        inputs = self._validate_controller(controller)
        if inputs[self.result_input][1].value is not None or any(
                obj.get(self.marker) and obj.parent == controller for obj in bpy.data.objects):
            raise ValueError('Rebind or remove the existing original-span provider before installing another')
        used = {obj.get(self.key_property) for obj in bpy.data.objects}
        key = next((i for i in range(1, 2**24) if i not in used), None)
        if key is None:
            raise ValueError('Original-span provider identity capacity exceeded')
        child = bpy.data.objects.new(controller.name+' / original span provider', controller.data)
        collection.objects.link(child)
        child.parent = controller
        child.matrix_parent_inverse = Matrix.Identity(4)
        child.matrix_basis = Matrix.Identity(4)
        child.hide_render = True; child.hide_select = True; child.hide_set(True)
        child[self.marker] = 1
        mod = child.modifiers.new('Shared native original-span construction', 'NODES')
        mod.node_group = self.solver; mod.show_viewport = False
        controller[self.key_property] = key
        set_input(controller, self.key_input, key)
        set_input(controller, self.result_input, child)
        try:
            self.rebind(controller)
        except Exception:
            set_input(controller, self.result_input, None)
            set_input(controller, self.key_input, 0)
            del controller[self.key_property]
            bpy.data.objects.remove(child, do_unlink=True)
            raise
        mod.show_viewport = True
        return child

    def rebind(self, controller):
        """Refresh IDs and scalar drivers after supported source regrouping.

        The caller must have already rebuilt topology tables and the dependent
        patch/strip/separator providers. This is not a topology-editing command.
        Materials are assigned live by the common entry, independently of this
        snapshot of their IDs on the producer.
        """
        inputs = self._validate_controller(controller)
        child = inputs[self.result_input][1].value
        key = inputs[self.key_input][1].value
        if (child is None or child.parent != controller or child.get(self.marker) != 1
                or child.type != 'MESH' or len(child.modifiers) != 1
                or child.modifiers[0].node_group != self.solver
                or type(key) is not int or not 0 < key < 2**24
                or controller.get(self.key_property) != key):
            raise ValueError('The original-span provider binding does not belong to this controller')
        child.data = controller.data
        destinations = self._inputs(child)
        for name, (socket, prop) in destinations.items():
            value = inputs[name][1].value
            if value == child or (isinstance(value, bpy.types.Collection) and child.name in value.all_objects):
                raise ValueError('Original-span inputs cannot include their own provider')
        for name, (socket, prop) in destinations.items():
            origin = inputs[name][1]
            prop.value = origin.value
            if socket.socket_type not in self.scalar_types:
                continue
            prop.driver_remove('value')
            curve = prop.driver_add('value'); curve.driver.type = 'AVERAGE'
            curve.extrapolation = 'LINEAR'
            for point in curve.keyframe_points: point.interpolation = 'LINEAR'
            variable = curve.driver.variables.new()
            variable.name = 'control'; variable.type = 'SINGLE_PROP'
            variable.targets[0].id = controller
            variable.targets[0].data_path = origin.path_from_id('value')
        child.update_tag()
        return child
