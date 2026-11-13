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
from .validation import bounds_overlap, validate_overlaps_steps
from .updates import drain_steps
from ....progressive.ownership import Ownership, alive
from . import VERSION

_ENGINES = {}
_UPDATING = set()
_ATTRIBUTES = (('st_section_id', 'INT'), ('st_width', 'FLOAT'),
               ('st_offset', 'FLOAT'), ('st_rise', 'FLOAT'), ('st_interpolation', 'INT'))


class StreetInputChanged(RuntimeError):
    """A pending publication lost its inputs; do not accept diagnostic output."""


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


class _SceneRevision:
    """Validate authored inputs across checkpoints without observing selection.

    Source geometry is checked once more before publication. Cheap controller
    checks run on each turn; the OSM session also supplies its terrain/scene
    validation. All producers retain an evaluated terrain signature.
    """
    def __init__(self, controller, network, validate=None):
        self.controller, self.network, self.validate_external = controller, network, validate
        self.scene = bpy.context.scene
        self.scale = _scale()
        self.source = controller['st_py_sources']
        self.output = controller['st_py_outputs']
        self.parents = tuple(controller.users_collection)
        self.signature = self._signature()
        self.network_signature = fingerprint(network)
        self.terrain = getattr(getattr(controller, 'blosm_street_terrain', None), 'target', None)
        # The shared OSM guard detects ordinary source edits cheaply, but
        # modifier dependency geometry/topology needs an evaluated final check.
        self.terrain_signature = self._terrain_signature() if self.terrain else None

    def _signature(self):
        controller = self.controller
        props = getattr(controller, 'blosm_street_terrain', None)
        terrain = tuple(getattr(props, field) for field in (
            'target', 'mode', 'clearance', 'maximum_grade', 'control_spacing',
            'sample_spacing', 'maximum_lift')) if props else ()
        return (tuple(controller.get(key) for key in (
            'st_py_network_uuid', 'st_py_source_data', 'st_py_manifest', 'st_py_features',
            'st_py_cell_size', 'st_py_resolution', 'st_py_exaggeration',
            'st_py_pavement_material', 'st_py_grass_material')),
            getattr(controller, 'blosm_street_pavement', None), terrain,
            tuple(v for row in controller.matrix_world for v in row))

    def _terrain_signature(self):
        from .terrain.blender import TerrainSampler
        return TerrainSampler(self.terrain, self.controller, self.scale).signature

    def _linked_to_scene(self):
        # Alive datablocks can still have been unlinked through the Outliner.
        # Do not accept a component whose controls or outputs became invisible
        # merely because its collection pointers are unchanged.
        required = {self.source, self.output, *self.parents}
        pending, visited = [self.scene.collection], set()
        while pending and required:
            collection = pending.pop()
            if collection in visited:
                continue
            visited.add(collection)
            required.discard(collection)
            pending.extend(collection.children)
        return not required

    def __call__(self, *, final=False):
        if self.validate_external:
            self.validate_external()
        controller = self.controller
        if not all(alive(item) for item in (self.scene, controller, self.source, self.output)):
            raise StreetInputChanged('Street update targets were deleted')
        if bpy.context.scene != self.scene or _scale() != self.scale:
            raise StreetInputChanged('Street update scene or units changed')
        if (controller.get('st_py_sources') != self.source or controller.get('st_py_outputs') != self.output
                or tuple(controller.users_collection) != self.parents or self._signature() != self.signature):
            raise StreetInputChanged('Street update inputs changed; prepare the streets again')
        if not self._linked_to_scene():
            raise StreetInputChanged('Street source or output collection was unlinked from the import scene')
        if self.terrain is not None and not alive(self.terrain):
            raise StreetInputChanged('Street terrain was deleted')
        if final:
            try:
                source_signature = fingerprint(SourceObjects.read(controller))
            except (ValueError, RuntimeError) as exc:
                raise StreetInputChanged('Authored street sources changed during generation') from exc
            if source_signature != self.network_signature:
                raise StreetInputChanged('Authored street sources changed during generation')
            if self.terrain_signature is not None and self._terrain_signature() != self.terrain_signature:
                raise StreetInputChanged('Street terrain changed during generation')

    def final(self):
        self(final=True)


def _guard_steps(steps, validate):
    """Close nested work on cancellation and validate each resumed checkpoint."""
    try:
        while True:
            validate()
            try:
                command = next(steps)
            except StopIteration as finished:
                return finished.value
            yield command
    finally:
        steps.close()


class _SceneOwnership(Ownership):
    """Share explicit allocation records with the containing import session."""
    def __init__(self, parent=None):
        super().__init__()
        self.parent = parent

    def own(self, kind, item):
        super().own(kind, item)
        if self.parent is not None:
            self.parent.own(kind, item)
        return item


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
                # Blender matrices are float32. Inverting/multiplying them in
                # float32 can introduce a spurious source Z after translating
                # a controller, especially in centimetre scenes. Cancel the
                # shared frame in double precision before reading controls.
                transform = np.linalg.solve(np.asarray(controller.matrix_world, dtype=float),
                                            np.asarray(obj.matrix_world, dtype=float))
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
                    coordinates = np.asarray([vertex.co for vertex in ordered], dtype=float)
                    points = tuple(map(tuple, (coordinates @ transform[:3,:3].T + transform[:3,3]) * scale))
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
    def _mesh(name, mesh, scale, *, ownership=None):
        from .uv import publication_uv_coordinates
        mesh.validate()
        result = bpy.data.meshes.new(name)
        if ownership is not None:
            ownership.own('meshes', result)
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
            for name, values in mesh.uv_layers.items():
                layer = result.uv_layers.new(name=name, do_init=False)
                if layer.name != name:
                    raise ValueError('UV map name cannot be represented in Blender: ' + name)
                layer.uv.foreach_set('vector', publication_uv_coordinates(values).reshape(-1))
            if result.uv_layers.get('UVMap') is not None:
                result.uv_layers.active = result.uv_layers['UVMap']
                result.uv_layers['UVMap'].active_render = True
            return result
        except Exception:
            bpy.data.meshes.remove(result)
            raise

    def publish(self, controller, network, engine, update, *, validate=None, on_progress=None,
                parent_ownership=None):
        return drain_steps(self.publish_steps(controller, network, engine, update,
            validate=validate or getattr(self, '_validate', None),
            on_progress=on_progress or getattr(self, '_on_progress', None), parent_ownership=parent_ownership))

    def publish_steps(self, controller, network, engine, update, *, validate=None, on_progress=None,
                      parent_ownership=None):
        """Stage detached owned outputs; perform the complete swap in one turn.

        Closing this producer removes only its pending IDs. No accepted output
        or controller metadata changes before the final non-yielding commit.
        """
        collection = controller['st_py_outputs']
        owned = {Key.parse(obj['st_py_key']): obj for obj in collection.objects
                 if obj.get('st_py_controller') == controller and obj.get('st_py_key')}
        if len(owned) != sum(obj.get('st_py_controller') == controller for obj in collection.objects):
            raise ValueError('Duplicate or invalid generated ownership')
        original_ids = {key: (obj, obj.data) for key, obj in owned.items()}
        ownership = _SceneOwnership(parent_ownership)
        staged = []
        committed = False

        def check():
            if validate:
                validate()
            if not alive(controller) or not alive(collection) or controller.get('st_py_outputs') != collection:
                raise StreetInputChanged('Street publication target was deleted or replaced')
            current = {Key.parse(obj['st_py_key']): (obj, obj.data) for obj in collection.objects
                       if obj.get('st_py_controller') == controller and obj.get('st_py_key')}
            if (current != original_ids or len(current) !=
                    sum(obj.get('st_py_controller') == controller for obj in collection.objects)):
                raise StreetInputChanged('Generated street ownership changed during publication')
            if any(not alive(item) for _, item in ownership.pending.values()):
                raise StreetInputChanged('Pending street output was deleted')

        def phase(message):
            check()
            if on_progress:
                on_progress(message)

        try:
            phase('Fitting street terrain')
            yield None
            check()
            from .terrain.blender import ScenePlacement
            placement = yield from _guard_steps(ScenePlacement.create_steps(
                controller, network, update, owned, _scale(), fit_cache=engine.terrain_constraints,
                ownership=ownership, on_progress=on_progress), check)
            update = placement.update
            yield None
            phase('Validating stored street geometry')
            changed = {key: p.mesh for key, p in update.products.items()
                       if key in update.rebuilt and isinstance(p, Part) and p.mesh is not None}
            scale = _scale()
            # Compare actual stored precision for both rebuilt and unchanged
            # seams, retaining the existing 20 um overlap tolerance.
            nearby = {}
            for key, mesh in changed.items():
                check()
                stored = Mesh(publication_coordinates(mesh.xyz, scale, physical=True), mesh.quads)
                if stored.validate()['maximum_edge_ratio'] > 8.:
                    raise ValueError(f'Published quad edge ratio exceeds 8: {key.text()}')
                nearby[key] = stored
                yield None
            changed_bounds = tuple(mesh.bounds for mesh in nearby.values())
            for key, obj in owned.items():
                check()
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
                yield None
            phase('Checking stored street overlaps (complete component)')
            yield None
            check()
            update.timings['publication_overlap_check'] = yield from _guard_steps(validate_overlaps_steps(
                nearby, changed=changed, cell_size=max(1., network.settings.spacing * 2.), tolerance=2e-5,
                on_progress=on_progress), check)
            del nearby
            yield None
            for key in update.rebuilt:
                phase('Staging street object ' + key.text())
                product = update.products[key]
                if isinstance(product, Part):
                    mesh = self._mesh(key.text(), product.mesh, scale, ownership=ownership)
                    obj = ownership.own('objects', bpy.data.objects.new(key.text(), mesh))
                    staged.append(obj)
                    from .materials import resolve, assign
                    assign(mesh, resolve(controller, product.mesh.material, self.materials))
                    obj['st_py_owner'] = product.owner.text()
                    product.metadata['bounds'] = product.mesh.bounds
                    product.metadata['terrain_material'] = product.mesh.material
                    obj['st_py_metadata'] = json.dumps(product.metadata, separators=(',', ':'), allow_nan=False)
                elif isinstance(product, AnchorTransform):
                    obj = ownership.own('objects', bpy.data.objects.new(key.text(), None))
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
                placement.decorate(obj, key)
                yield None
            phase('Committing street component')
            # Final input validation precedes serialization and scene mutation.
            # The callback supplied by NetworkScene includes source geometry.
            if validate:
                getattr(validate, 'final', validate)()
            update.timings['published_objects'] = len(staged)
            retired = set(update.removed) | {key for key in update.rebuilt
                if not isinstance(update.products[key], (Part, AnchorTransform))}
            update.timings['removed_objects'] = sum(key in owned for key in retired)
            manifest = dict(version=1, implementation=VERSION, network_uuid=update.network_id, revision=engine.updates.revision + 1,
                source_fingerprints={k.text(): v for k, v in update.source_fingerprints.items()},
                task_fingerprints={k.text(): v for k, v in update.task_fingerprints.items()},
                dependencies={k.text(): [d.text() for d in ds] for k, ds in update.dependencies.items()})
            properties = dict(st_py_source_data=json.dumps(network.to_dict(), separators=(',', ':'), allow_nan=False),
                st_py_features=json.dumps(engine.feature_names), st_py_manifest=json.dumps(manifest, separators=(',', ':')),
                st_py_error='', st_py_terrain_signature=placement.signature, st_py_terrain_mode=placement.mode)
            if placement.report is not None:
                properties['st_py_terrain_report'] = json.dumps(placement.report, separators=(',', ':'))
            retiring = [owned[key] for key in set(update.rebuilt) | set(update.removed) if key in owned]
            old_links = [(obj, tuple(obj.users_collection)) for obj in retiring]
            material_slots = [(obj.data, tuple(obj.data.materials)) for key, obj in owned.items()
                              if key not in set(update.rebuilt) | set(update.removed) and obj.type == 'MESH']
            snapshots = {name: controller.get(name) for name in (*properties, 'st_py_last_update')}
            coordinator = engine.updates
            engine_state = {name: getattr(coordinator, name) for name in (
                'network_id', 'source_fingerprints', 'task_fingerprints', 'dependencies', 'products', 'revision')}
            try:
                # No suspension from the first link through retirement. Old
                # objects remain allocated until every fallible commit step
                # succeeds, allowing complete restoration of the prior state.
                for obj in staged:
                    collection.objects.link(obj)
                for obj, links in old_links:
                    for owner in links:
                        owner.objects.unlink(obj)
                for name, value in properties.items():
                    controller[name] = value
                from .materials import refresh
                update.timings['material_objects_updated'] = refresh(controller, update.products, self.materials)
                controller['st_py_last_update'] = json.dumps(dict(rebuilt=[k.text() for k in update.rebuilt],
                    removed=[k.text() for k in update.removed], timings=update.timings), separators=(',', ':'))
                engine.commit(placement.raw)
                ownership.accept()
                committed = True
            except BaseException:
                for name, value in engine_state.items():
                    setattr(coordinator, name, value)
                for mesh, slots in material_slots:
                    mesh.materials.clear()
                    for material in slots:
                        mesh.materials.append(material)
                for obj, links in old_links:
                    for owner in links:
                        if obj.name not in owner.objects:
                            owner.objects.link(obj)
                for name, value in snapshots.items():
                    if value is None:
                        if name in controller:
                            del controller[name]
                    else:
                        controller[name] = value
                raise
            # Postcommit disposal errors retain the coherent accepted result.
            # Old IDs are detached, so a failed removal cannot show duplicates.
            cleanup_errors = []
            for old in retiring:
                try:
                    data = old.data
                    bpy.data.objects.remove(old, do_unlink=True)
                    if isinstance(data, bpy.types.Mesh) and data.users == 0:
                        bpy.data.meshes.remove(data)
                except Exception as exc:
                    cleanup_errors.append(str(exc))
            try:
                engine.release_published_meshes()
            except Exception as exc:
                cleanup_errors.append(str(exc))
            if cleanup_errors:
                update.timings['publication_cleanup_errors'] = cleanup_errors
            return update
        finally:
            if not committed:
                ownership.rollback()


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
        return drain_steps(NetworkScene._update_steps(controller, force=force, materials=materials,
                                                     blocking=True))

    @staticmethod
    def update_steps(controller, *, force=False, materials=None, validate=None, on_progress=None,
                     parent_ownership=None):
        """Cooperative Update with the same numerical and publication algorithms."""
        return (yield from NetworkScene._update_steps(controller, force=force, materials=materials,
                    validate=validate, on_progress=on_progress, parent_ownership=parent_ownership))

    @staticmethod
    def _update_steps(controller, *, force=False, materials=None, validate=None, on_progress=None,
                      blocking=False, parent_ownership=None):
        identity = None
        original_engine = engine = None
        completed = False
        try:
            register_runtime()
            _prune_runtime()
            if on_progress:
                on_progress('Reading street sources')
            identity = controller['st_py_network_uuid']
            if identity in _UPDATING:
                identity = None
                raise RuntimeError('This street network already has an update in progress')
            _UPDATING.add(identity)
            network = SourceObjects.read(controller)
            revision = _SceneRevision(controller, network, validate)
            yield None
            revision()
            original_engine = engine = _ENGINES.get(network.uuid)
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
                                tuple(obj.get('st_terrain_prepared_position',
                                    tuple(float(v) * _scale() for v in obj.location))),
                                float(obj.rotation_euler.z), obj.get('st_py_asset', 'placeholder'))
                _ENGINES[network.uuid] = engine
            update = yield from _guard_steps(engine.prepare_steps(network, force=force,
                                                                 on_progress=on_progress), revision)
            writer = SceneWriter(materials)
            if blocking:
                # Preserve the public synchronous publication hook used by
                # extensions and fault-injection tests. It drains publish_steps.
                writer._validate, writer._on_progress = revision, on_progress
                result = writer.publish(controller, network, engine, update)
            else:
                result = yield from writer.publish_steps(controller, network, engine, update,
                    validate=revision, on_progress=on_progress, parent_ownership=parent_ownership)
            completed = True
            return result
        except Exception as exc:
            if alive(controller):
                controller['st_py_error'] = str(exc)
            raise
        finally:
            if identity is not None:
                _UPDATING.discard(identity)
                if not completed and engine is not None and _ENGINES.get(identity) is engine:
                    if original_engine is None:
                        _ENGINES.pop(identity, None)
                    else:
                        _ENGINES[identity] = original_engine
