"""Bounded positive seeds from exact one-vertex corner half-planes.

Physical rim vertices and quad connectivity are immutable. This is only an
initial computational atlas repair; shape, height and publication acceptance
remain the responsibility of the consuming chart.
"""
import numpy as np
from .mesh import Mesh, cross


def _jacobian(xy, faces):
    edges = np.roll(xy[faces], -1, axis=1)-xy[faces]
    return cross(edges, np.roll(edges, -1, axis=1))


def positive_seeds(mesh, *, max_bad_faces=4, max_candidates=4, max_cells=8000):
    if (type(max_bad_faces) is not int or not 1 <= max_bad_faces <= 4
            or type(max_candidates) is not int or not 1 <= max_candidates <= 4
            or type(max_cells) is not int or not 1 <= max_cells <= 50000):
        raise ValueError('Seed untangling requires bounded face and candidate budgets')
    if (len(mesh.quads) > max_cells or not np.isfinite(mesh.xyz).all()
            or np.ptp(mesh.xyz[:, 2]) > 1e-10):
        raise ValueError('Seed untangling requires a bounded planar quad chart')
    xy, faces = mesh.xyz[:, :2], mesh.quads
    initial = _jacobian(xy, faces)
    bad = np.flatnonzero(initial.min(axis=1) <= 1e-10)
    if not len(bad):
        mesh.validate()
        return ((mesh, dict(already_positive=True, original_topology=True)),)
    if len(bad) > max_bad_faces:
        raise ValueError('Seed untangling exceeds its nonpositive-face budget')
    edges = np.sort(np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    if counts.max() > 2:
        raise ValueError('Seed untangling requires manifold quad edge incidence')
    fixed = set(map(int, unique[counts == 1].ravel()))
    free = sorted(set(map(int, faces[bad].ravel()))-fixed)
    groups = []
    for vertex in free:
        incident = np.flatnonzero(np.any(faces == vertex, axis=1))
        if len(incident) > 64:
            continue
        local_faces = faces[incident]
        origin = xy[vertex]
        points = xy[local_faces]-origin
        selected = local_faces == vertex

        def jacobian(point):
            candidate = points.copy()
            candidate[selected] = point
            edge = np.roll(candidate, -1, axis=1)-candidate
            return cross(edge, np.roll(edge, -1, axis=1)).ravel()

        constant = jacobian((0., 0.))
        coefficients = np.c_[jacobian((1., 0.))-constant, jacobian((0., 1.))-constant]
        extent = float(np.linalg.norm(points, axis=-1).max())
        if not np.isfinite(extent) or extent <= 1e-8:
            continue
        polygon = [np.array((-extent, -extent)), np.array((extent, -extent)),
                   np.array((extent, extent)), np.array((-extent, extent))]
        margin = 1e-7
        for normal, value in zip(coefficients, constant-margin):
            if not polygon:
                break
            clipped = []
            for a, b in zip(polygon, polygon[1:]+polygon[:1]):
                first, last = float(normal @ a+value), float(normal @ b+value)
                if first >= 0.:
                    clipped.append(a)
                if (first >= 0.) != (last >= 0.):
                    clipped.append(a+first/(first-last)*(b-a))
            polygon = clipped
        if len(polygon) < 3:
            continue
        polygon = np.asarray(polygon)
        center = polygon.mean(axis=0)
        closest = []
        for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
            delta = b-a
            if delta @ delta > 0.:
                closest.append(a+np.clip(-a @ delta/(delta @ delta), 0., 1.)*delta)
        if not closest:
            continue
        nearest = min(closest, key=lambda p: float(p @ p))
        candidates = []
        for blend in (.05, 1.):
            point = (1.-blend)*nearest+blend*center
            if np.linalg.norm(point) > .5*extent:
                continue
            xyz = mesh.xyz.copy()
            xyz[vertex, :2] = origin+point
            candidate = Mesh(xyz, faces, mesh.material, dict(mesh.seams), dict(mesh.uv_layers))
            try:
                candidate.validate()
                stored = Mesh(xyz.astype(np.float32), faces, mesh.material)
                stored.validate()
            except ValueError:
                continue
            ids = np.asarray(sorted(fixed), dtype=np.int32)
            if not np.array_equal(candidate.xyz[ids], mesh.xyz[ids]):
                raise ValueError('Seed untangling changed a physical rim vertex')
            witness = dict(construction='exact affine one-vertex corner half-planes',
                vertex=vertex, original_xy=origin.tolist(), prepared_xy=xyz[vertex, :2].tolist(),
                displacement=float(np.linalg.norm(point)), displacement_bound=.5*extent,
                incident_faces=incident.tolist(), original_nonpositive_faces=bad.tolist(),
                halfplane_margin=margin, interior_blend=blend,
                minimum_incident_jacobian=float(_jacobian(xyz[:, :2], local_faces).min()),
                minimum_stored_incident_jacobian=float(_jacobian(stored.xyz[:, :2], local_faces).min()),
                complete_mesh_positive=True, complete_stored_mesh_positive=True,
                exact_physical_rim=True, original_topology=True, shape_qualified=False)
            candidates.append((candidate, witness))
        if candidates:
            groups.append(candidates)
    result = []
    for choice in range(2):
        for candidates in groups:
            if choice < len(candidates):
                result.append(candidates[choice])
                if len(result) == max_candidates:
                    return tuple(result)
    if not result:
        raise ValueError('No bounded one-vertex positive computational seed exists')
    return tuple(result)
