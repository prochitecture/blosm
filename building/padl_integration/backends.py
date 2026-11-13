"""Backend strategies share inputs, identity, PADL evaluation and terrain rules."""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from contextlib import contextmanager, nullcontext

from ...padl import services
from ...padl.api import create_input_object
from ...padl.errors import PadlError
from ...padl.bindings import FaceType
from ...padl.geometry import MeshData
from ...padl.persistence import prepare_owner, read_footprints
from ...padl.python_backend import PythonBackend
from ...padl.terrain import BlenderTerrainSampler
from .assets import BuildingProgramLibrary


@dataclass
class GenerationResult:
    owner: object
    output: object
    statuses: list = field(default_factory=list)


class BuildingBackend(ABC):
    @contextmanager
    def batch(self, library, role, *, progressive=False):
        yield

    @abstractmethod
    def generate(self, owner, program, *, library, role, terrain, collection, output_object=None):
        """Return an output and per-footprint diagnostics, retaining the input."""


class _ProgressiveGraph:
    """Revalidate shared assets after external edits, without hashing each turn."""
    def __init__(self, library, role, tree):
        self.library, self.role = library, role
        self.dirty = self.busy = False
        self._remember(tree)
        self.handler = self._updated

    def _remember(self, tree):
        import bpy
        groups, pending = {}, [tree]
        while pending:
            group = pending.pop()
            pointer = group.as_pointer()
            if pointer in groups:
                continue
            groups[pointer] = group
            pending.extend(node.node_tree for node in group.nodes
                           if getattr(node, 'node_tree', None) is not None)
        self.groups = tuple(groups.values())
        self.group_ids = frozenset(groups)
        self.documents = tuple((group, key, value, value.as_string())
                               for group in self.groups for key, value in group.items()
                               if isinstance(value, bpy.types.Text))

    def _updated(self, scene, graph):
        session = getattr(getattr(self.library, 'app', None), 'osm_session', None)
        if self.busy or (session is not None and session.in_step):
            return
        # Publishing a building also produces deferred Object/Mesh/Collection
        # updates after the timer returns. Those do not change the generator.
        # Blender reports edited root/nested NodeTrees by their original ID,
        # even when all is_updated_* flags are false (e.g. socket/node edits).
        if any(update.id.original.as_pointer() in self.group_ids for update in graph.updates):
            self.dirty = True

    def checked_tree(self, tree):
        from ...progressive.ownership import alive
        if any(not alive(group) for group in self.groups):
            raise PadlError('A building Geometry Nodes group was deleted during the import', code='missing_asset')
        # Text edits do not necessarily trigger a dependency-graph notification.
        for group, key, text, content in self.documents:
            if (not alive(text) or group.get(key) != text or text.as_string() != content):
                raise PadlError('A building Geometry Nodes document changed during the import', code='asset_conflict')
        if self.dirty:
            tree = self.library.tree(self.role)
            self._remember(tree)
            self.dirty = False
        return tree


class GeometryNodesBackend(BuildingBackend):
    @contextmanager
    def batch(self, library, role, *, progressive=False):
        import bpy
        from .gn_runtime import StatusReader
        tree = library.tree(role)
        # Reuse verification until an external edit invalidates it. Progressive
        # turns retain no temporary reader in the scene or undo snapshots.
        with nullcontext() if progressive else StatusReader(tree) as reader:
            self._batch_tree, self._batch_reader = tree, reader
            self._batch_guard = guard = _ProgressiveGraph(library, role, tree) if progressive else None
            if guard is not None:
                bpy.app.handlers.depsgraph_update_post.append(guard.handler)
            try:
                yield
            finally:
                if guard is not None and guard.handler in bpy.app.handlers.depsgraph_update_post:
                    bpy.app.handlers.depsgraph_update_post.remove(guard.handler)
                self._batch_tree = self._batch_reader = None
                self._batch_guard = None

    def generate(self, owner, program, *, library, role, terrain, collection, output_object=None):
        from .gn_runtime import attach, read_status
        tree = getattr(self, '_batch_tree', None)
        if tree is None:
            tree = library.tree(role)
        else:
            from ...progressive.ownership import alive
            if not alive(tree):
                raise PadlError('The building Geometry Nodes group was deleted during the import', code='missing_asset')
        guard = getattr(self, '_batch_guard', None)
        if guard is not None:
            tree = self._batch_tree = guard.checked_tree(tree)
        if tree.get('padl_program_sha256') != program.fingerprint:
            raise PadlError('This PADL document needs a matching precompiled GN program', code='program_mismatch')
        if guard is not None:
            guard.busy = True
        try:
            modifier = attach(owner, tree, terrain, collection)
            reader = getattr(self, '_batch_reader', None)
            statuses = reader.read(owner, modifier) if reader else read_status(owner, modifier)
        finally:
            if guard is not None:
                guard.busy = False
        return GenerationResult(owner, owner, statuses)


class PythonBuildingBackend(BuildingBackend):
    def generate(self, owner, program, *, library, role, terrain, collection, output_object=None):
        import bpy
        mesh, statuses = MeshData(), []
        sampler = BlenderTerrainSampler(terrain)
        backend = PythonBackend(available=library.available_roofs(role))
        for footprint in read_footprints(owner):
            if not FaceType.generates_volume(footprint.face_type):
                continue
            try:
                result = backend.generate(program, [footprint], seed=owner['padl_seed'], terrain=sampler)
                mesh.extend(result.mesh)
                statuses.append({'footprint_id': footprint.footprint_id, 'error': 0})
            except PadlError as exc:
                statuses.append({'footprint_id': footprint.footprint_id, 'error': exc.code, 'message': str(exc)})
        if output_object is owner:
            raise PadlError('Python output must be separate from the editable input', code='invalid_input')
        generated = mesh.to_blender_mesh(owner.name + ' Volumes')
        if output_object is None:
            output = bpy.data.objects.new(owner.name + ' Volumes', generated)
            collection.objects.link(output)
        else:
            output = output_object
            previous = output.data
            output.data = generated
            if previous.users == 0:
                bpy.data.meshes.remove(previous)
        from mathutils import Matrix
        output.matrix_world = Matrix.Identity(4)
        owner.hide_render = True
        owner.hide_set(True)
        return GenerationResult(owner, output, statuses)


class BuildingGenerationService:
    """Add future modes as program roles and future backends as strategies."""
    backends = {'geometry_nodes': GeometryNodesBackend, 'python': PythonBuildingBackend}

    def __init__(self, *, library=None, backend='geometry_nodes', role='simple_3d'):
        services.configure_blosm(package=__package__.removesuffix('.building.padl_integration'))
        self.library = library or BuildingProgramLibrary()
        self.backend = self.backends[backend]() if isinstance(backend, str) else backend
        self.role = role

    def batch(self, *, progressive=False):
        """Scope verified assets and temporary inspection helpers to one import."""
        return self.backend.batch(self.library, self.role, progressive=progressive)

    def generate(self, footprints, *, collection, terrain=None, seed=0, namespace='osm', name='Buildings'):
        import bpy
        footprints = tuple(footprints)
        program = self.library.program(self.role)
        owner = create_input_object(footprints, collection=collection, name=name)
        try:
            prepare_owner(owner, program, {f.footprint_id: f.identity for f in footprints}, namespace=namespace, seed=seed)
            return self.backend.generate(owner, program, library=self.library, role=self.role,
                                         terrain=terrain, collection=collection)
        except Exception:
            mesh = owner.data
            bpy.data.objects.remove(owner, do_unlink=True)
            if mesh.users == 0:
                bpy.data.meshes.remove(mesh)
            raise

    def regenerate(self, owner, *, collection, terrain=None, output_object=None):
        """Generate from current editable inputs; never reimport or release H."""
        from ...padl.persistence import load_program
        # Call on a Python owner. GN owners update interactively on mesh edits.
        return self.backend.generate(owner, load_program(owner), library=self.library,
                                     role=self.role, terrain=terrain, collection=collection, output_object=output_object)
