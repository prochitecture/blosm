"""Attach and inspect precompiled groups; never import the GN compiler."""
from ...util.geometry_nodes import getGnInput, iterGnInputs, setGnInput
from ...padl.errors import PadlError


def socket_id(tree, name):
    return next(item.identifier for item in tree.interface.items_tree
                if item.item_type == 'SOCKET' and item.in_out == 'INPUT' and item.name == name)


def terrain_limits(tree, terrain, collection):
    import bpy
    if terrain is None:
        return None
    pending, visited, helper = [tree], set(), None
    while pending:
        group = pending.pop()
        if group.as_pointer() in visited:
            continue
        visited.add(group.as_pointer())
        if group.get('prochitecture_terrain_z_limits_v1'):
            helper = group
            break
        pending.extend(n.node_tree for n in group.nodes if getattr(n, 'node_tree', None))
    if helper is None:
        raise PadlError('Building asset lacks its terrain bounds helper', code='incompatible_asset')
    for owner in bpy.data.objects:
        if owner.type == 'MESH' and len(owner.modifiers) == 1:
            mod = owner.modifiers[0]
            if mod.type == 'NODES' and mod.node_group == helper and getGnInput(mod, socket_id(helper, 'terrain')) == terrain:
                return owner
    owner = bpy.data.objects.new('PADL Terrain Limits', bpy.data.meshes.new('PADL terrain bounds'))
    collection.objects.link(owner)
    owner.hide_render = True
    owner.hide_set(True)
    mod = owner.modifiers.new('Shared terrain bounds', 'NODES')
    mod.node_group = helper
    setGnInput(mod, socket_id(helper, 'terrain'), terrain)
    return owner


def attach(owner, tree, terrain, collection):
    axes = owner.matrix_world.to_3x3()
    if (any(abs(axes.col[i].length - 1) > 1e-6 for i in range(3)) or
            any(abs(axes.col[i].dot(axes.col[j])) > 1e-6 for i, j in ((0, 1), (0, 2), (1, 2))) or
            abs(axes.col[2].x) > 1e-6 or abs(axes.col[2].y) > 1e-6 or
            axes.col[2].z < 1 - 1e-6 or axes.determinant() < 0):
        raise PadlError('PADL GN input owners must be upright with unit scale and no mirroring/shear', code='unsupported_transform')
    mod = owner.modifiers.get('PADL Buildings')
    if mod is None:
        mod = owner.modifiers.new('PADL Buildings', 'NODES')
    elif mod.type != 'NODES':
        raise PadlError('PADL Buildings modifier name is already in use', code='asset_conflict')
    mod.node_group = tree
    for name, value in (('Seed', owner['padl_seed']), ('Terrain', terrain),
                        ('Terrain limits', terrain_limits(tree, terrain, collection))):
        setGnInput(mod, socket_id(tree, name), value)
    return mod


class StatusReader:
    """Read diagnostics, optionally in a reusable graph outside the import scene."""
    def __init__(self, tree, *, private_scene=False):
        import bpy
        self.private_scene = private_scene
        self.source_scene = bpy.context.scene if private_scene else None
        self.scene = self.owner = self.empty_mesh = self.modifier = None
        self.wrapper = bpy.data.node_groups.new('PADL temporary status reader', 'GeometryNodeTree')
        try:
            for item in tree.interface.items_tree:
                if item.item_type == 'SOCKET' and item.in_out == 'INPUT':
                    self.wrapper.interface.new_socket(name=item.name, in_out='INPUT', socket_type=item.socket_type)
            self.wrapper.interface.new_socket(name='Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
            incoming = self.wrapper.nodes.new('NodeGroupInput')
            outgoing = self.wrapper.nodes.new('NodeGroupOutput')
            group = self.wrapper.nodes.new('GeometryNodeGroup')
            group.node_tree = tree
            for item in tree.interface.items_tree:
                if item.item_type == 'SOCKET' and item.in_out == 'INPUT':
                    self.wrapper.links.new(incoming.outputs[item.name],group.inputs[item.name])
            self.wrapper.links.new(group.outputs['Status'],outgoing.inputs['Geometry'])
            if private_scene:
                from ...progressive.ownership import TRANSIENT
                # Rebuilding relations in the growing import scene for every
                # temporary owner is expensive. Keep this small graph alive
                # for the batch and discard its IDs when undo restores a turn.
                self.wrapper[TRANSIENT] = True
                self.scene = bpy.data.scenes.new('PADL private status evaluation')
                self.scene[TRANSIENT] = True
                self.scene.unit_settings.scale_length = self.source_scene.unit_settings.scale_length
                self.empty_mesh = bpy.data.meshes.new('PADL temporary status empty mesh')
                self.empty_mesh[TRANSIENT] = True
                self.owner = bpy.data.objects.new('PADL temporary status owner', self.empty_mesh)
                self.owner[TRANSIENT] = True
                self.owner.hide_render = True
                self.scene.collection.objects.link(self.owner)
                self.modifier = self.owner.modifiers.new('PADL private Status', 'NODES')
                self.modifier.node_group = self.wrapper
        except Exception:
            self.close()
            raise

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        import bpy
        from ...progressive.ownership import alive
        # The input mesh and terrain objects are borrowed. Remove only the
        # reader's own IDs, including after user deletion or history cleanup.
        for kind, item in (('objects', self.owner), ('scenes', self.scene),
                           ('meshes', self.empty_mesh), ('node_groups', self.wrapper)):
            if item is not None and alive(item):
                getattr(bpy.data, kind).remove(item, do_unlink=True)

    @staticmethod
    def _statuses(mesh):
        return [{attr.name: attr.data[i].value for attr in mesh.attributes
                 if attr.domain == 'POINT' and not attr.name.startswith('.') and hasattr(attr.data[i], 'value')}
                for i in range(len(mesh.vertices))]

    def _read_private(self, owner, modifier):
        import bpy
        from ...progressive.ownership import alive
        if self.source_scene is None or not alive(self.source_scene):
            raise PadlError('The import scene used for building status evaluation was deleted', code='missing_asset')
        if (any(item is None or not alive(item)
                for item in (self.scene, self.owner, self.empty_mesh, self.wrapper, self.modifier)) or
                self.owner.name not in self.scene.objects or self.modifier.node_group != self.wrapper):
            raise PadlError('The temporary building status evaluator was deleted or changed during the import', code='missing_asset')
        evaluated = None
        try:
            self.owner.data = owner.data
            self.owner.matrix_world = owner.matrix_world
            for value in iterGnInputs(modifier):
                bound = getGnInput(modifier, value.identifier)
                setGnInput(self.modifier, socket_id(self.wrapper, value.name), bound)
                if isinstance(bound, bpy.types.Object) and bound.name not in self.scene.objects:
                    self.scene.collection.objects.link(bound)
            self.owner.update_tag()
            # Terrain modifiers can be animated or use Scene Time. Evaluate
            # them at the import scene's time, including its seconds-per-frame.
            source = self.source_scene
            changed = (self.scene.frame_current != source.frame_current or
                       self.scene.frame_subframe != source.frame_subframe)
            for name in ('fps', 'fps_base', 'frame_map_old', 'frame_map_new'):
                value = getattr(source.render, name)
                if getattr(self.scene.render, name) != value:
                    setattr(self.scene.render, name, value)
                    changed = True
            # Subdivision terrain must obey the same viewport simplification
            # cap as its visible counterpart in the import scene.
            for name in ('use_simplify', 'simplify_subdivision'):
                value = getattr(source.render, name)
                if getattr(self.scene.render, name) != value:
                    setattr(self.scene.render, name, value)
            if changed:
                self.scene.frame_set(source.frame_current, subframe=source.frame_subframe)
            layer = self.scene.view_layers[0]
            layer.update()
            # context.evaluated_depsgraph_get() can select the main scene's
            # graph even with a context override. Use the private layer itself.
            graph = layer.depsgraph
            if graph.scene.original != self.scene:
                raise PadlError('Unable to evaluate building status in its temporary scene', code='invalid_input')
            evaluated = self.owner.evaluated_get(graph)
            return self._statuses(evaluated.to_mesh())
        finally:
            if evaluated is not None:
                evaluated.to_mesh_clear()
            # Never retain an imported input mesh through the private owner
            # between UI turns (including if the user then deletes that input).
            if alive(self.owner) and alive(self.empty_mesh):
                self.owner.data = self.empty_mesh

    def read(self, owner, modifier):
        if self.private_scene:
            return self._read_private(owner, modifier)
        import bpy
        temporary, evaluated = None, None
        try:
            temporary = owner.copy()
            temporary.name = 'PADL temporary status owner'
            temporary.hide_render = True
            bpy.context.collection.objects.link(temporary)
            found = False
            for item in list(temporary.modifiers):
                if item.name == modifier.name:
                    item.node_group = self.wrapper
                    for value in iterGnInputs(modifier):
                        setGnInput(item, socket_id(self.wrapper,value.name),getGnInput(modifier,value.identifier))
                    found = True
                elif found:
                    temporary.modifiers.remove(item)
            bpy.context.view_layer.update()
            evaluated = temporary.evaluated_get(bpy.context.evaluated_depsgraph_get())
            mesh = evaluated.to_mesh()
            return self._statuses(mesh)
        finally:
            if evaluated:
                evaluated.to_mesh_clear()
            if temporary:
                bpy.data.objects.remove(temporary,do_unlink=True)


def read_status(owner, modifier):
    """Evaluate separate Status geometry without touching the visible group."""
    with StatusReader(modifier.node_group) as reader:
        return reader.read(owner,modifier)
