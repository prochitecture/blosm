"""Blender input transport and explicit publication for the numerical backend.

Only authored edge meshes and a saved controller are read. Output ownership is
by network UUID and typed key, independent of names and collection iteration.
"""
from dataclasses import replace
import json
import math
import bpy
import bmesh
import numpy as np
from bpy.app.handlers import persistent
from .backend import PythonMeshBackend, DEFAULT_FEATURES
from .model import Network, Key, fingerprint
from .providers import Part, AnchorTransform
from .mesh import Mesh, publication_coordinates
from .validation import bounds_overlap, validate_overlaps
from . import VERSION

_ENGINES = {}
_ATTRIBUTES = (('st_section_id', 'INT'), ('st_width', 'FLOAT'),
               ('st_offset', 'FLOAT'), ('st_rise', 'FLOAT'), ('st_interpolation', 'INT'))


def _scale():
    value = bpy.context.scene.unit_settings.scale_length
    if not math.isfinite(value) or value <= 0:
        raise ValueError('Invalid scene unit scale')
    return value


def controller_for(obj):
    if obj is None:
        return None
    if obj.get('st_py_network_uuid') and obj.get('st_py_source_data'):
        return obj
    return obj.get('st_py_controller')


def clear_runtime():
    _ENGINES.clear()


@persistent
def _clear_on_file_load(_unused):
    clear_runtime()


def register_runtime():
    if _clear_on_file_load not in bpy.app.handlers.load_pre:
        bpy.app.handlers.load_pre.append(_clear_on_file_load)


def unregister_runtime():
    if _clear_on_file_load in bpy.app.handlers.load_pre:
        bpy.app.handlers.load_pre.remove(_clear_on_file_load)
    clear_runtime()


def _prune_runtime():
    """Release engines whose persistent controller has been deleted.

    This runs on explicit Update, never from drawing or geometry-edit events.
    Other scenes in the same blend file can still own live controllers.
    """
    if _ENGINES:
        live = {o.get('st_py_network_uuid') for o in bpy.data.objects if o.get('st_py_source_data')}
        for identity in tuple(_ENGINES):
            if identity not in live:
                del _ENGINES[identity]


class SourceObjects:
    @staticmethod
    def create(network, collection, controller, grouping='street'):
        groups = {}
        for section in network.sections:
            owner = section.street_id if grouping == 'street' else 0
            groups.setdefault(owner, []).append(section)
        for identity, sections in sorted(groups.items()):
            name = f'Street {identity}' if grouping == 'street' else 'All Streets'
            points, edges, roles = [], [], []
            attributes = {name: [] for name, kind in _ATTRIBUTES}
            for section in sections:
                first = len(points)
                points.extend(np.asarray(section.points) / _scale())
                edges.extend((first + i, first + i + 1) for i in range(len(section.points) - 1))
                roles.extend([1] + [0] * (len(section.points) - 2) + [2])
                values = (section.id, section.width, section.offset, section.rise, section.interpolation)
                for (name_attr, kind), value in zip(_ATTRIBUTES, values):
                    attributes[name_attr].extend([value] * (len(section.points) - 1))
            mesh = bpy.data.meshes.new(name + ' input')
            mesh.from_pydata(points, edges, [])
            for name_attr, kind in _ATTRIBUTES:
                mesh.attributes.new(name_attr, kind, 'EDGE').data.foreach_set('value', attributes[name_attr])
            mesh.attributes.new('st_endpoint_role', 'INT', 'POINT').data.foreach_set('value', roles)
            obj = bpy.data.objects.new(name, mesh)
            obj.parent = controller
            collection.objects.link(obj)
            obj['st_py_controller'] = controller
            obj.display_type = 'WIRE'
            obj.show_in_front = True

    @staticmethod
    def read(controller):
        base = Network.from_dict(json.loads(controller['st_py_source_data']))
        known = base.section_index
        sections = {}
        collection = controller['st_py_sources']
        scale = _scale()
        for obj in collection.all_objects:
            if obj.type != 'MESH' or obj.get('st_py_controller') != controller:
                continue
            if any(abs(v - 1.) > 1e-6 for v in obj.scale):
                raise ValueError('Source objects require unit scale; edit vertices or dimensions')
            temporary = obj.mode != 'EDIT'
            bm = bmesh.new() if temporary else bmesh.from_edit_mesh(obj.data)
            if temporary:
                bm.from_mesh(obj.data)
            try:
                layers = {}
                for name, kind in _ATTRIBUTES:
                    layer = (bm.edges.layers.int if kind == 'INT' else bm.edges.layers.float).get(name)
                    if layer is None:
                        raise ValueError('Missing authored attribute: ' + name)
                    layers[name] = layer
                role = bm.verts.layers.int.get('st_endpoint_role')
                if role is None or bm.faces or any(not v.link_edges for v in bm.verts):
                    raise ValueError('Inputs require marked open edge chains without faces or loose metadata vertices')
                transform = controller.matrix_world.inverted() @ obj.matrix_world
                groups = {}
                for edge in bm.edges:
                    groups.setdefault(edge[layers['st_section_id']], []).append(edge)
                for identity, edges in groups.items():
                    if identity not in known or identity in sections:
                        raise ValueError('Unknown or duplicated Section; use a managed topology operation')
                    vertices = {v for e in edges for v in e.verts}
                    starts = [v for v in vertices if v[role] == 1]
                    ends = [v for v in vertices if v[role] == 2]
                    if len(starts) != 1 or len(ends) != 1:
                        raise ValueError('Section endpoint roles were changed')
                    allowed = set(edges)
                    current, previous = starts[0], None
                    ordered = [current]
                    while current != ends[0]:
                        options = [e for e in current.link_edges if e in allowed and e is not previous]
                        if len(options) != 1:
                            raise ValueError('Section is branched, disconnected or cyclic')
                        edge = options[0]
                        current, previous = edge.other_vert(current), edge
                        if current in ordered:
                            raise ValueError('Section is cyclic')
                        ordered.append(current)
                    if len(ordered) != len(vertices) or len(edges) != len(ordered) - 1:
                        raise ValueError('Section topology is inconsistent')
                    if len(ordered) != len(known[identity].points):
                        raise ValueError('Section topology changed; use a managed topology operation')
                    values = {}
                    for name, field in (('st_width', 'width'), ('st_offset', 'offset'),
                                        ('st_rise', 'rise'), ('st_interpolation', 'interpolation')):
                        data = [e[layers[name]] for e in edges]
                        if max(data) - min(data) > 1e-6:
                            raise ValueError('Section attributes must be uniform along its edges')
                        values[field] = data[0]
                    points = tuple(tuple(float(v) * scale for v in transform @ vertex.co) for vertex in ordered)
                    sections[identity] = replace(known[identity], points=points, **values)
            finally:
                if temporary:
                    bm.free()
        if sections.keys() != known.keys():
            raise ValueError('An authored Section is missing from the source collection')
        settings = replace(base.settings, cell_size=controller['st_py_cell_size'],
                           resolution=controller['st_py_resolution'], exaggeration=controller['st_py_exaggeration'])
        return replace(base, sections=tuple(sections[s.id] for s in base.sections), settings=settings)


class SceneWriter:
    def __init__(self, materials=None):
        self.materials = materials or {}

    @staticmethod
    def _mesh(name, mesh, scale):
        result = bpy.data.meshes.new(name)
        try:
            result.vertices.add(len(mesh.xyz))
            result.vertices.foreach_set('co', publication_coordinates(mesh.xyz, scale).ravel())
            result.loops.add(mesh.quads.size)
            result.loops.foreach_set('vertex_index', mesh.quads.ravel())
            result.polygons.add(len(mesh.quads))
            result.polygons.foreach_set('loop_start', np.arange(len(mesh.quads), dtype=np.int32) * 4)
            result.polygons.foreach_set('loop_total', np.full(len(mesh.quads), 4, dtype=np.int32))
            result.polygons.foreach_set('use_smooth', np.ones(len(mesh.quads), dtype=bool))
            result.update(calc_edges=True)
            return result
        except Exception:
            bpy.data.meshes.remove(result)
            raise

    def publish(self, controller, network, engine, update):
        collection = controller['st_py_outputs']
        owned = {Key.parse(obj['st_py_key']): obj for obj in collection.objects
                 if obj.get('st_py_controller') == controller and obj.get('st_py_key')}
        if len(owned) != sum(obj.get('st_py_controller') == controller for obj in collection.objects):
            raise ValueError('Duplicate or invalid generated ownership')
        changed = {key: p.mesh for key, p in update.products.items()
                   if key in update.rebuilt and isinstance(p, Part) and p.mesh is not None}
        scale = _scale()
        # Validate the representation that will actually be stored. Comparing
        # a rebuilt double-precision seam with an unchanged float32 seam can
        # manufacture an overlap when both round to identical stored points.
        nearby = {key: Mesh(publication_coordinates(mesh.xyz, scale, physical=True), mesh.quads)
                  for key, mesh in changed.items()}
        for key, mesh in nearby.items():
            if mesh.validate()['maximum_edge_ratio'] > 8.:
                raise ValueError(f'Published quad edge ratio exceeds 8: {key.text()}')
        changed_bounds = tuple(mesh.bounds for mesh in nearby.values())
        for key, obj in owned.items():
            if key in update.rebuilt or key in update.removed or obj.type != 'MESH':
                continue
            product = update.products.get(key)
            if not isinstance(product, Part) or not any(bounds_overlap(
                    publication_coordinates(product.metadata['bounds'], scale, physical=True), bounds)
                                          for bounds in changed_bounds):
                continue
            xyz = np.empty(len(obj.data.vertices) * 3, dtype=np.float32)
            ids = np.empty(len(obj.data.loops), dtype=np.int32)
            obj.data.vertices.foreach_get('co', xyz)
            obj.data.loops.foreach_get('vertex_index', ids)
            nearby[key] = Mesh(xyz.reshape(-1, 3).astype(float) * scale, ids.reshape(-1, 4))
        # Retain the existing 20 um publication tolerance. Every comparison
        # now uses actual stored precision, independently of coordinate size.
        update.timings['publication_overlap_check'] = validate_overlaps(nearby, changed=changed,
            cell_size=max(1., network.settings.spacing * 2.), tolerance=2e-5)
        staged = []
        snapshots = {name: controller.get(name) for name in
                     ('st_py_source_data', 'st_py_manifest', 'st_py_last_update', 'st_py_error', 'st_py_features')}
        try:
            for key in update.rebuilt:
                product = update.products[key]
                if isinstance(product, Part):
                    mesh = self._mesh(key.text(), product.mesh, _scale())
                    try:
                        obj = bpy.data.objects.new(key.text(), mesh)
                    except Exception:
                        bpy.data.meshes.remove(mesh)
                        raise
                    # Register immediately: metadata/material/property errors
                    # below must not leak a detached object or its mesh.
                    staged.append(obj)
                    material = self.materials.get(product.mesh.material)
                    if material is None:
                        material = controller.get('st_py_' + product.mesh.material + '_material')
                    if material is not None:
                        mesh.materials.append(material)
                    obj['st_py_owner'] = product.owner.text()
                    product.metadata['bounds'] = product.mesh.bounds
                    obj['st_py_metadata'] = json.dumps(product.metadata, separators=(',', ':'), allow_nan=False)
                elif isinstance(product, AnchorTransform):
                    obj = bpy.data.objects.new(key.text(), None)
                    staged.append(obj)
                    obj.empty_display_type = 'SINGLE_ARROW'
                    obj.empty_display_size = 3. / _scale()
                    obj.location = np.asarray(product.position) / _scale()
                    obj.rotation_euler.z = product.heading
                    obj['st_py_asset'] = product.asset
                else:
                    continue
                obj.parent = controller
                obj['st_py_controller'] = controller
                obj['st_py_key'] = key.text()
            for obj in staged:
                collection.objects.link(obj)
            update.timings['published_objects'] = len(staged)
            retired = set(update.removed) | {key for key in update.rebuilt
                if not isinstance(update.products[key], (Part, AnchorTransform))}
            update.timings['removed_objects'] = sum(key in owned for key in retired)
            controller['st_py_source_data'] = json.dumps(network.to_dict(), separators=(',', ':'), allow_nan=False)
            controller['st_py_features'] = json.dumps(engine.feature_names)
            # Compute a prospective manifest without advancing the engine. A
            # publication failure can still restore the complete previous state.
            manifest = dict(version=1, implementation=VERSION, network_uuid=update.network_id, revision=engine.updates.revision + 1,
                source_fingerprints={k.text(): v for k, v in update.source_fingerprints.items()},
                task_fingerprints={k.text(): v for k, v in update.task_fingerprints.items()},
                dependencies={k.text(): [d.text() for d in ds] for k, ds in update.dependencies.items()})
            controller['st_py_manifest'] = json.dumps(manifest, separators=(',', ':'))
            controller['st_py_last_update'] = json.dumps(dict(rebuilt=[k.text() for k in update.rebuilt],
                removed=[k.text() for k in update.removed], timings=update.timings), separators=(',', ':'))
            controller['st_py_error'] = ''
            engine.commit(update)
        except Exception:
            for obj in staged:
                data = obj.data
                bpy.data.objects.remove(obj, do_unlink=True)
                if isinstance(data, bpy.types.Mesh) and data.users == 0:
                    bpy.data.meshes.remove(data)
            for name, value in snapshots.items():
                if value is None:
                    if name in controller:
                        del controller[name]
                else:
                    controller[name] = value
            raise
        for key in set(update.rebuilt) | set(update.removed):
            if key not in owned:
                continue
            old = owned[key]
            data = old.data
            bpy.data.objects.remove(old, do_unlink=True)
            if isinstance(data, bpy.types.Mesh) and data.users == 0:
                bpy.data.meshes.remove(data)
        engine.release_published_meshes()
        return update


class NetworkScene:
    @classmethod
    def create(cls, network, parent=None, *, name='Python streets', grouping='street', materials=None,
               generate=True, features=None):
        if any(o.get('st_py_network_uuid') == network.uuid for o in bpy.data.objects):
            raise ValueError('A controller already owns this network UUID')
        parent = parent or bpy.context.scene.collection
        collection = bpy.data.collections.new(name)
        parent.children.link(collection)
        controller = bpy.data.objects.new(name, None)
        collection.objects.link(controller)
        controller['st_py_network_uuid'] = network.uuid
        controller['st_py_source_data'] = json.dumps(network.to_dict(), separators=(',', ':'))
        controller['st_py_backend'] = 'python'
        controller['st_py_features'] = json.dumps(tuple(DEFAULT_FEATURES if features is None else features))
        controller['st_py_cell_size'] = network.settings.cell_size
        controller['st_py_resolution'] = network.settings.resolution
        controller['st_py_exaggeration'] = network.settings.exaggeration
        controller.id_properties_ui('st_py_cell_size').update(min=.05, soft_min=.1, description='Physical cell size in metres')
        controller.id_properties_ui('st_py_resolution').update(min=1, max=4)
        controller.id_properties_ui('st_py_exaggeration').update(min=0., soft_max=10.)
        for role in ('sources', 'outputs'):
            child = bpy.data.collections.new(name + ' / ' + role)
            collection.children.link(child)
            controller['st_py_' + role] = child
        for key, material in (materials or {}).items():
            controller['st_py_' + key + '_material'] = material
        SourceObjects.create(network, controller['st_py_sources'], controller, grouping)
        if generate:
            cls.update(controller, materials=materials)
        return controller

    @staticmethod
    def update(controller, *, force=False, materials=None):
        try:
            register_runtime()
            _prune_runtime()
            network = SourceObjects.read(controller)
            engine = _ENGINES.get(network.uuid)
            saved = json.loads(controller.get('st_py_manifest', '{}'))
            features = tuple(json.loads(controller.get('st_py_features', json.dumps(DEFAULT_FEATURES))))
            if (features == ('paths', 'intersections', 'carriageways', 'furniture')
                    and saved.get('implementation') == '0.1.0.dev1'):
                features = DEFAULT_FEATURES
            if (engine is None or engine.feature_names != features
                    or (saved and fingerprint(saved) != fingerprint(engine.manifest()))):
                engine = PythonMeshBackend(network.settings, feature_names=features)
                if saved:
                    engine.updates.restore_manifest(saved)
                    # Recover persistent output identities and small metadata,
                    # not whole geometry buffers. Derived paths/layouts are
                    # reconstructed as needed; unchanged meshes retain their
                    # existing Blender objects even on the first Update after
                    # reopening or Undo/Redo.
                    for obj in controller['st_py_outputs'].objects:
                        if obj.get('st_py_controller') != controller or not obj.get('st_py_key'):
                            continue
                        key = Key.parse(obj['st_py_key'])
                        if obj.type == 'MESH' and obj.get('st_py_metadata') and obj.get('st_py_owner'):
                            engine.updates.products[key] = Part(Key.parse(obj['st_py_owner']), None,
                                json.loads(obj['st_py_metadata']))
                        elif key.kind == 'placement':
                            engine.updates.products[key] = AnchorTransform(
                                tuple(float(v) * _scale() for v in obj.location),
                                float(obj.rotation_euler.z), obj.get('st_py_asset', 'placeholder'))
                _ENGINES[network.uuid] = engine
            update = engine.prepare(network, force=force)
            return SceneWriter(materials).publish(controller, network, engine, update)
        except Exception as exc:
            controller['st_py_error'] = str(exc)
            raise
