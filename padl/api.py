"""Shared footprint operations and the Python geometry API; no GN dependency."""
from .bindings import FACE_TYPES, FootprintBindings
from .errors import PadlError
from .model import Program
from .parser import PadlParser
from .persistence import load_program, prepare_owner, read_footprints


def parse(document):
    if isinstance(document, Program):
        return document
    if isinstance(document, dict):
        return PadlParser().parse(document)
    return PadlParser().load(document)


def create_input_object(footprints, *, name='PADL Footprints', collection=None):
    import bpy
    footprints = tuple(footprints)
    vertices, faces, rows = [], [], []
    seen = set()
    bindings = FootprintBindings()
    for footprint in footprints:
        footprint.validate()
        if footprint.footprint_id in seen:
            raise PadlError('Duplicate footprint ID', code='invalid_identity')
        seen.add(footprint.footprint_id)
        for i, ring in enumerate((footprint.outer, *footprint.holes)):
            start = len(vertices)
            vertices.extend((v[0], v[1], v[2] if len(v) > 2 else 0) for v in ring)
            faces.append(tuple(range(start, len(vertices))))
            if i == 0:
                rows.append(bindings.encode(footprint.attributes, footprint_id=footprint.footprint_id, has_holes=bool(footprint.holes)))
            else:
                rows.append({'face_type': 2, 'footprint_id': 0, 'parent_id': footprint.footprint_id, 'has_holes': False})
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(vertices, [], faces)
    for name in sorted({key for row in rows for key in row}):
        attr = mesh.attributes.new(name, FACE_TYPES[name], 'FACE')
        attr.data.foreach_set('value', [row.get(name, False if FACE_TYPES[name] == 'BOOLEAN' else 0) for row in rows])
    owner = bpy.data.objects.new(mesh.name, mesh)
    (collection or bpy.context.collection).objects.link(owner)
    return owner


def generate_python(owner, document=None, *, source_keys=None, namespace=None, seed=None, terrain=None,
                    terrain_policy='error', output_object=None, available=None):
    """Generate a result object; the input owner and current inputs remain intact."""
    import bpy
    from .python_backend import PythonBackend
    from .terrain import BlenderTerrainSampler
    program = parse(document) if document is not None else load_program(owner)
    if source_keys is not None:
        prepare_owner(owner, program, source_keys, namespace=namespace, seed=seed)
    elif seed is not None:
        set_seed(owner, seed)
    result = PythonBackend(available=available).generate(program, read_footprints(owner),
                    seed=owner.get('padl_seed', seed), terrain=BlenderTerrainSampler(terrain, policy=terrain_policy))
    mesh = result.mesh.to_blender_mesh(owner.name + ' Buildings')
    if output_object is None:
        output_object = bpy.data.objects.new(mesh.name, mesh)
        bpy.context.collection.objects.link(output_object)
    elif output_object is owner:
        bpy.data.meshes.remove(mesh)
        raise PadlError('Python output must be separate from the editable input owner', code='invalid_input')
    else:
        output_object.data = mesh
    # MeshData is in document/world coordinates, including terrain translation.
    from mathutils import Matrix
    output_object.matrix_world = Matrix.Identity(4)
    return output_object, result


def edit_parameter(owner, footprint_id, name, value=None, *, clear=False):
    """Set current value/presence together; never mark a derived value as supplied."""
    bindings = FootprintBindings()
    rows, outer, _ = bindings.validate_mesh(owner.data)
    if footprint_id not in outer or name not in bindings.parameters:
        raise PadlError('Unknown footprint or editable parameter', code='invalid_input')
    binding = bindings.parameters[name]
    index = outer[footprint_id]
    if name.startswith('has_'):
        raise PadlError('Use clear/set on the associated value so presence changes atomically', code='invalid_input')
    if clear:
        if binding.presence:
            owner.data.attributes[binding.presence].data[index].value = False
        elif name == 'roof_type':
            owner.data.attributes[name].data[index].value = 0
        elif name == 'min_elevation':
            owner.data.attributes[name].data[index].value = 0
        else:
            raise PadlError(f'{name} cannot be cleared', code='invalid_input')
    else:
        candidate = dict(rows[index])
        if binding.categories:
            reverse = {v: k for k, v in binding.categories.items()}
            if value not in reverse:
                raise PadlError(f'Unsupported {name} {value!r}', code='invalid_category')
            value = reverse[value]
        bindings.check(name, value, binding.data_type)
        if name == 'levels_wall_count' and value <= 0:
            raise PadlError('Wall level count must be positive', code='invalid_count')
        candidate[name] = value
        if binding.presence:
            candidate[binding.presence] = True
        bindings.decode(candidate)
        attr = owner.data.attributes.get(name) or owner.data.attributes.new(name, binding.data_type, 'FACE')
        attr.data[index].value = value
        if binding.presence:
            owner.data.attributes[binding.presence].data[index].value = True
    owner.data.update()
    owner.update_tag()


def set_seed(owner, seed):
    """Set the Python execution seed; regenerate to apply the new samples."""
    FootprintBindings.check('padl_seed', seed, 'INT')
    owner['padl_seed'] = seed
