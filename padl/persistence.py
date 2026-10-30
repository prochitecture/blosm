"""User-approved minimal Blender persistence; source inputs remain editable."""
import json

from .bindings import FACE_TYPES, OPTIONAL, FootprintBindings
from .errors import PadlError
from .geometry import Footprint
from .parser import PadlParser
from .randomness import stable_key

def source_identity(namespace, key):
    if not isinstance(namespace, str) or not namespace or not isinstance(key, str) or not key:
        raise PadlError('Dataset namespace and source keys must be nonempty strings', code='invalid_identity')
    return json.dumps([namespace, key], ensure_ascii=False, separators=(',', ':'))

def _text(name, data, prior=None):
    import bpy
    content = json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n'
    if isinstance(prior, bpy.types.Text) and prior.as_string() == content:
        return prior
    text = bpy.data.texts.new(name)
    text.write(content)
    return text

def prepare_owner(owner, program, source_keys, *, namespace, seed=None):
    """Validate first, then initialize application metadata and absent flags.

    source_keys maps integer association IDs to full source keys (e.g.
    ``osm:way:123``). Dataset namespace is part of the evaluation identity.
    """
    if owner.type != 'MESH':
        raise PadlError('PADL input owner must be a mesh object', code='invalid_input')
    bindings = FootprintBindings()
    rows, outer, holes = bindings.validate_mesh(owner.data)
    if not isinstance(source_keys, dict):
        raise PadlError('Source keys must be a map from integer footprint IDs to strings', code='invalid_identity')
    try:
        keys = {int(k): v for k, v in source_keys.items()
                if (type(k) is int or isinstance(k, str) and str(int(k)) == k)}
    except (TypeError, ValueError, OverflowError) as exc:
        raise PadlError('Invalid source-key footprint ID', code='invalid_identity') from exc
    if any(not isinstance(v, str) or not v for v in keys.values()):
        raise PadlError('Source keys must be nonempty strings', code='invalid_identity')
    if len(keys) != len(source_keys) or set(keys) != set(outer) or len(set(keys.values())) != len(keys):
        raise PadlError('Source-key map must identify each outer footprint exactly once', code='invalid_identity')
    identities = {fid: source_identity(namespace, key) for fid, key in keys.items()}
    seed = program.context.get('seed', 0) if seed is None else seed
    bindings.check('padl_seed', seed, 'INT')
    mesh = owner.data
    for name in ('min_elevation', 'building_use', 'roof_type', *('has_' + key for key in OPTIONAL), 'padl_random_key'):
        if mesh.attributes.get(name) is None:
            mesh.attributes.new(name, FACE_TYPES[name], 'FACE')
    random_keys = [0] * len(mesh.polygons)
    for fid, index in outer.items():
        random_keys[index] = stable_key(('footprint', identities[fid]))
    mesh.attributes['padl_random_key'].data.foreach_set('value', random_keys)
    owner['padl_seed'] = seed
    owner['padl_dataset_namespace'] = namespace
    owner['padl_source_keys'] = _text('PADL source identities', {str(k): v for k, v in keys.items()}, owner.get('padl_source_keys'))
    owner['padl_document'] = _text('PADL input document', program.document, owner.get('padl_document'))
    mesh.update()
    owner.update_tag()
    return identities


def load_program(owner):
    """Read the authored document saved on a Python input owner."""
    import bpy
    text = owner.get('padl_document')
    if not isinstance(text, bpy.types.Text):
        raise PadlError('PADL document Text reference is missing', code='missing_metadata')
    return PadlParser().loads(text.as_string(), filename=text.name)


def read_footprints(owner):
    """Read current inputs without retaining any original imported-value history."""
    import bpy
    bindings = FootprintBindings()
    rows, outer, holes = bindings.validate_mesh(owner.data)
    record = owner.get('padl_source_keys')
    if not isinstance(record, bpy.types.Text):
        raise PadlError('Source-identity Text reference is missing', code='missing_metadata')
    try:
        keys = json.loads(record.as_string())
    except (TypeError, ValueError) as exc:
        raise PadlError('Invalid source-identity map', code='invalid_identity') from exc
    if not isinstance(keys, dict) or set(keys) != {str(fid) for fid in outer}:
        raise PadlError('Source-identity map does not match current footprints', code='invalid_identity')
    if any(not isinstance(v, str) or not v for v in keys.values()) or len(set(keys.values())) != len(keys):
        raise PadlError('Duplicate source identities', code='invalid_identity')
    mesh, namespace = owner.data, owner.get('padl_dataset_namespace')
    # Build geometry in the document/world frame; placement in an owner's local
    # frame is performed by the backend adapter, not by changing PADL elevations.
    def ring(index):
        return tuple(tuple(owner.matrix_world @ mesh.vertices[i].co) for i in mesh.polygons[index].vertices)
    result = []
    for fid, index in outer.items():
        identity = source_identity(namespace, keys[str(fid)])
        expected = stable_key(('footprint', identity))
        if rows[index].get('padl_random_key') != expected:
            raise PadlError('Random identity binding does not match source identity; prepare the owner again', code='invalid_identity')
        result.append(Footprint(identity, fid, ring(index), tuple(ring(i) for i in holes.get(fid, [])), bindings.decode(rows[index])))
    return result
