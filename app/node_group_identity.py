"""Runtime graph identity checks independent of the generator source checkout."""
import hashlib
import json


def scalar(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if hasattr(value, 'name_full'):
        return {'id_type': value.bl_rna.identifier, 'name': value.name_full}
    try:
        return [scalar(item) for item in value]
    except TypeError:
        return str(value)


def properties(obj):
    ignored = {'rna_type', 'name', 'label', 'select', 'location', 'width', 'height',
               'color', 'use_custom_color', 'hide', 'show_options', 'show_preview', 'show_texture'}
    result = {}
    for prop in obj.bl_rna.properties:
        if prop.identifier in ignored or prop.is_readonly:
            continue
        if prop.type in {'BOOLEAN', 'INT', 'FLOAT', 'STRING', 'ENUM'}:
            result[prop.identifier] = scalar(getattr(obj, prop.identifier))
    return result


def graph_fingerprint(group, memo=None, active=None):
    """Include dependency graphs, constants, interface and dynamic item collections."""
    memo = {} if memo is None else memo
    active = set() if active is None else active
    key = group.as_pointer()
    if key in memo:
        return memo[key]
    if key in active:
        raise ValueError('Recursive node-group dependency')
    active.add(key)
    nodes = []
    for node in group.nodes:
        dynamic = {}
        for prop in node.bl_rna.properties:
            if prop.type == 'COLLECTION' and prop.identifier.endswith('_items'):
                dynamic[prop.identifier] = [
                    (getattr(item, 'name', ''), properties(item)) for item in getattr(node, prop.identifier)]
        child = getattr(node, 'node_tree', None)
        sockets = lambda seq: [(s.identifier, s.bl_idname, scalar(getattr(s, 'default_value', None))) for s in seq]
        nodes.append((node.name, node.bl_idname, properties(node), sockets(node.inputs), sockets(node.outputs), dynamic,
                      graph_fingerprint(child, memo, active) if child else None))
    record = {'nodes': nodes,
              'links': sorted((link.from_node.name, link.from_socket.identifier,
                               link.to_node.name, link.to_socket.identifier, link.multi_input_sort_id) for link in group.links),
              'interface': [(s.name, s.identifier, s.in_out, s.socket_type, properties(s), scalar(getattr(s, 'default_value', None)))
                            for s in group.interface.items_tree if s.item_type == 'SOCKET']}
    digest = hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    active.remove(key)
    memo[key] = digest
    return digest
