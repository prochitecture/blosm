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
    """One temporary graph for an import batch; inputs remain per modifier."""
    def __init__(self, tree):
        import bpy
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
        except Exception:
            bpy.data.node_groups.remove(self.wrapper)
            raise

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        import bpy
        bpy.data.node_groups.remove(self.wrapper)

    def read(self, owner, modifier):
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
            return [{attr.name: attr.data[i].value for attr in mesh.attributes
                     if attr.domain == 'POINT' and not attr.name.startswith('.') and hasattr(attr.data[i], 'value')}
                    for i in range(len(mesh.vertices))]
        finally:
            if evaluated:
                evaluated.to_mesh_clear()
            if temporary:
                bpy.data.objects.remove(temporary,do_unlink=True)


def read_status(owner, modifier):
    """Evaluate separate Status geometry without touching the visible group."""
    with StatusReader(modifier.node_group) as reader:
        return reader.read(owner,modifier)
