"""Python geometry backend using PADL evaluation and independent roof adapters.

The gabled adapter reuses Blosm's profile geometry; the hipped adapter reuses Blosm's existing bpypolyskel kernel directly, without
PML, legacy building orchestration, GN compilation, or silent flat-roof fallback.
Triangulation and skeleton math use Blender's mathutils at geometry-execution time.
"""
from dataclasses import dataclass
import math

from .errors import PadlError
from .bindings import FaceType
from .evaluation import Evaluator, FootprintContext
from .geometry import MeshData, contains, signed_area


def direction(footprint, options):
    ring = footprint.outer
    reference = options.get('reference_source', 'longest_edge')
    alignment = options.get('alignment', 'along')
    orientation = options.get('orientation', 'ccw').lower()
    if reference not in ('longest_edge', 'edge_index') or alignment not in ('along', 'perpendicular') or orientation not in ('ccw', 'cw'):
        raise PadlError('Unsupported roof direction configuration', code='unsupported_feature')
    lengths = [math.dist(a[:2], b[:2]) for a, b in zip(ring, ring[1:] + ring[:1])]
    if reference == 'longest_edge':
        maximum = max(lengths)
        index = next(i for i, length in enumerate(lengths) if length >= maximum - 0.0001)
    else:
        index = options.get('edge_index', 0)
        if type(index) not in (int, float) or int(index) != index:
            raise PadlError('roof.edge_index must be integral', code='invalid_input')
        index = int(index) % len(ring)
    a, b = ring[index], ring[(index + 1) % len(ring)]
    dx, dy = b[0] - a[0], b[1] - a[1]
    length = math.hypot(dx, dy)
    dx, dy = dx / length, dy / length
    if orientation == 'cw':
        dx, dy = -dx, -dy
    if alignment == 'perpendicular':
        dx, dy = -dy, dx
    return dy, -dx  # uphill is ridge direction cross Z, matching public GN groups


def wall_faces(mesh, ring, height, bottom, footprint_id):
    for a, b in zip(ring, ring[1:] + ring[:1]):
        mesh.polygon(((a[0], a[1], bottom), (b[0], b[1], bottom),
                      (b[0], b[1], height(b)), (a[0], a[1], height(a))), 1, footprint_id)


class PlanarRoofAdapter:
    def generate(self, footprint, volume, datum=0):
        if volume.roof_type == 'flat' and not footprint.holes:
            mesh = MeshData()
            ring = footprint.oriented_rings()[0]
            height = datum + volume.roof_min_elevation
            mesh.polygon([(p[0], p[1], height) for p in ring], 2, footprint.footprint_id)
            wall_faces(mesh, ring, lambda p: height, datum + volume.min_elevation, footprint.footprint_id)
            return mesh
        from mathutils import Vector
        from mathutils.geometry import delaunay_2d_cdt
        mesh = MeshData()
        rings = footprint.oriented_rings()
        points = [Vector(v[:2]) for ring in rings for v in ring]
        edges, start = [], 0
        for ring in rings:
            edges.extend((start + i, start + (i + 1) % len(ring)) for i in range(len(ring)))
            start += len(ring)
        verts, _, faces, *_ = delaunay_2d_cdt(points, edges, [], 0, 1e-7)
        base = datum + volume.roof_min_elevation
        if volume.roof_type == 'flat':
            height = lambda p: base
        else:
            axis = direction(footprint, volume.roof_options)
            project = lambda p: p[0] * axis[0] + p[1] * axis[1]
            low, high = min(map(project, footprint.outer)), max(map(project, footprint.outer))
            if high - low <= 0.0001 or volume.roof_height <= 0.0001:
                raise PadlError('Skillion direction span and rise must exceed 0.0001 m', code='invalid_geometry')
            height = lambda p: base + volume.roof_height * (project(p) - low) / (high - low)
        roof_area = 0
        for face in faces:
            polygon = [verts[i] for i in face]
            center = tuple(sum(p[i] for p in polygon) / len(polygon) for i in (0, 1))
            if not contains(center, rings[0]) or any(contains(center, ring) for ring in rings[1:]):
                continue
            if signed_area(polygon) < 0:
                polygon.reverse()
            roof_area += abs(signed_area(polygon))
            mesh.polygon([(p[0], p[1], height(p)) for p in polygon], 2, footprint.footprint_id)
        if abs(roof_area - footprint.area) > max(1e-4, footprint.area * 1e-5):
            raise PadlError('Triangulation does not cover the footprint minus its holes', code='solver_failure')
        for ring in rings:
            wall_faces(mesh, ring, height, datum + volume.min_elevation, footprint.footprint_id)
        return mesh


class GabledRoofAdapter:
    """Adapt Blosm's shared profile geometry to already resolved PADL dimensions."""
    profile = (((0., 0.), (.5, 1.), (1., 0.)), {'numSamples': 10, 'angleToHeight': .5})

    def generate(self, footprint, volume, datum=0):
        if footprint.holes:
            raise PadlError('Gabled roofs do not support holes in the initial profile', code='unsupported_geometry')
        from .profile_geometry import ProfiledRoofGeometry
        kernel = ProfiledRoofGeometry(self.profile)
        ring = footprint.oriented_rings()[0]
        axis = direction(footprint, volume.roof_options)
        # Keep kernel float32 coordinates near zero, independent of packing.
        ox = (min(p[0] for p in ring) + max(p[0] for p in ring)) / 2
        oy = (min(p[1] for p in ring) + max(p[1] for p in ring)) / 2
        local = tuple((p[0] - ox, p[1] - oy) for p in ring)
        try:
            kernel.prepare(local, bottom=volume.min_elevation,
                           roof_base=volume.roof_min_elevation,
                           roof_height=volume.roof_height, direction=axis)
            if not kernel.make(None):
                raise ValueError('Profile geometry generation did not complete')
        except Exception as exc:
            raise PadlError(f'Python gabled profile solver failed: {exc}', code='solver_failure') from exc
        mesh = MeshData()
        roof_area = 0.
        for kind, faces in ((1, kernel.wallIndices), (2, kernel.roofIndices)):
            for face in faces:
                points = [tuple(kernel.verts[i]) for i in face]
                if kind == 2:
                    area = signed_area(points)
                    roof_area += abs(area)
                    if area < 0:
                        points.reverse()
                mesh.polygon([(p[0] + ox, p[1] + oy, p[2] + datum) for p in points],
                             kind, footprint.footprint_id)
        if abs(roof_area - footprint.area) > max(1e-4, footprint.area * 1e-5):
            raise PadlError('Python gabled output failed footprint coverage validation', code='solver_failure')
        if not mesh.faces or abs(max(v[2] for v in mesh.vertices) - (datum + volume.max_elevation)) > max(1e-4, abs(datum + volume.max_elevation)*1e-6):
            raise PadlError('Python gabled output failed height validation', code='solver_failure')
        return mesh


class HippedRoofAdapter:
    def generate(self, footprint, volume, datum=0):
        from mathutils import Vector
        from .services import addon_service
        polygonize = addon_service('lib.bpypolyskel.bpypolyskel').polygonize
        if volume.roof_options:
            raise PadlError('Hipped direction/profile options are unsupported', code='unsupported_feature')
        rings = footprint.oriented_rings()
        base = datum + volume.roof_min_elevation
        verts = [Vector((v[0], v[1], base)) for ring in rings for v in ring]
        holes, start = [], len(rings[0])
        for ring in rings[1:]:
            holes.append((start, len(ring)))
            start += len(ring)
        try:
            faces = polygonize(verts, 0, len(rings[0]), holes or None, height=volume.roof_height)
        except Exception as exc:
            raise PadlError(f'Python hipped solver failed: {exc}', code='solver_failure') from exc
        if not faces or any(len(face) < 3 or len(set(face)) != len(face) or any(i < 0 or i >= len(verts) for i in face) for face in faces):
            raise PadlError('Python hipped solver produced invalid faces', code='solver_failure')
        roof_area = sum(abs(signed_area([verts[i] for i in face])) for face in faces)
        highest = max(v.z for v in verts)
        if abs(roof_area - footprint.area) > max(1e-4, footprint.area * 1e-5) or abs(highest - (base + volume.roof_height)) > 1e-4:
            raise PadlError('Python hipped output failed coverage/height validation', code='solver_failure')
        mesh = MeshData()
        for face in faces:
            polygon = [verts[i] for i in face]
            if signed_area(polygon) < 0:
                polygon.reverse()
            mesh.polygon(polygon, 2, footprint.footprint_id)
        for ring in rings:
            wall_faces(mesh, ring, lambda p: base, datum + volume.min_elevation, footprint.footprint_id)
        return mesh


@dataclass(frozen=True)
class PythonResult:
    mesh: MeshData
    evaluations: tuple
    terrain_samples: tuple = ()


class PythonBackend:
    def __init__(self, adapters=None, available=None):
        planar = PlanarRoofAdapter()
        self.adapters = {'flat': planar, 'skillion': planar, 'gabled': GabledRoofAdapter(), 'hipped': HippedRoofAdapter()}
        self.adapters.update(adapters or {})
        self.available = frozenset(self.adapters if available is None else available)

    def generate(self, program, footprints, *, seed=None, global_inputs=None, terrain=None):
        footprints = tuple(footprints)
        by_id = {f.footprint_id: f for f in footprints}
        evaluator = Evaluator(program, seed=seed, global_inputs=global_inputs)
        from .terrain import TerrainSample
        output, results, samples = MeshData(), [], []
        seen_ids, seen_keys = set(), set()
        for footprint in footprints:
            footprint.validate()
            if footprint.footprint_id in seen_ids or footprint.identity in seen_keys:
                raise PadlError('Repeated footprint association ID or source identity', code='invalid_identity')
            seen_ids.add(footprint.footprint_id)
            seen_keys.add(footprint.identity)
            if not FaceType.generates_volume(footprint.face_type):
                continue
            result = evaluator.evaluate(FootprintContext(footprint.identity, footprint.attributes))
            reference = footprint
            if footprint.face_type == FaceType.PART:
                reference = by_id.get(footprint.parent_id) or footprint.elevation_reference
                if reference is None or reference.footprint_id != footprint.parent_id or not FaceType.is_building(reference.face_type):
                    raise PadlError('Part elevation requires its owning building footprint', code='invalid_identity')
            # Coverage is required for both the part and the common datum rings.
            if terrain is not None and reference is not footprint:
                terrain.sample(footprint)
            sample = terrain.sample(reference) if terrain is not None else TerrainSample()
            datum = sample.minimum
            for volume in result.volumes:
                if volume.roof_type not in self.available:
                    raise PadlError(f'{volume.roof_type} is unavailable in this profile', code='unavailable_roof')
                output.extend(self.adapters[volume.roof_type].generate(footprint, volume, datum))
            results.append(result)
            samples.append(sample)
        return PythonResult(output, tuple(results), tuple(samples))
