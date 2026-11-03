"""Backend strategies share inputs, identity, PADL evaluation and terrain rules."""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from contextlib import contextmanager

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
    def batch(self, library, role):
        yield

    @abstractmethod
    def generate(self, owner, program, *, library, role, terrain, collection, output_object=None):
        """Return an output and per-footprint diagnostics, retaining the input."""


class GeometryNodesBackend(BuildingBackend):
    @contextmanager
    def batch(self, library, role):
        from .gn_runtime import StatusReader
        tree = library.tree(role)
        with StatusReader(tree) as reader:
            self._batch_tree, self._batch_reader = tree, reader
            try:
                yield
            finally:
                self._batch_tree = self._batch_reader = None

    def generate(self, owner, program, *, library, role, terrain, collection, output_object=None):
        from .gn_runtime import attach, read_status
        tree = getattr(self, '_batch_tree', None) or library.tree(role)
        if tree.get('padl_program_sha256') != program.fingerprint:
            raise PadlError('This PADL document needs a matching precompiled GN program', code='program_mismatch')
        modifier = attach(owner, tree, terrain, collection)
        reader = getattr(self, '_batch_reader', None)
        statuses = reader.read(owner, modifier) if reader else read_status(owner, modifier)
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

    def batch(self):
        """Scope verified assets and temporary inspection helpers to one import."""
        return self.backend.batch(self.library, self.role)

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
