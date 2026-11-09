"""Material resolution and publication, independent of geometry construction.

Roles are semantic identifiers, never material-slot indices. A future PADL
adapter can supply the same role mapping. Current generated meshes use slot 0.
"""


def resolve(controller, role, overrides=None):
    if overrides and overrides.get(role) is not None:
        return overrides[role]
    if role == 'pavement':
        selection = getattr(controller, 'blosm_street_pavement', None)
        if selection is not None:
            return selection
    return controller.get('st_py_' + role + '_material')


def assign(mesh, material):
    if material is not None:
        if len(mesh.materials):
            if mesh.materials[0] != material:
                mesh.materials[0] = material
        else:
            mesh.materials.append(material)


def refresh(controller, products, overrides=None):
    """Apply material-only edits to existing outputs without rebuilding meshes."""
    from .model import Key
    count = 0
    for obj in controller['st_py_outputs'].objects:
        if obj.type != 'MESH' or obj.get('st_py_controller') != controller:
            continue
        product = products.get(Key.parse(obj['st_py_key']))
        if product is None:
            continue
        role = product.mesh.material if product.mesh is not None else product.metadata.get('terrain_material', 'pavement')
        material = resolve(controller, role, overrides)
        if material is not None and (not obj.data.materials or obj.data.materials[0] != material):
            assign(obj.data, material)
            count += 1
    return count
