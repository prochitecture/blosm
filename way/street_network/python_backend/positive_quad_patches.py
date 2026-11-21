"""Bounded direct-quad seeds for a singular two-cell computational patch.

Some four-sided parity terminations give a collinear physical boundary vertex
no interior edge. Moving its fixed neighbours cannot fix that cell. Retile the
same six-node patch with two rotated quads or three quads and one interior
point. No triangles or new physical boundary samples are constructed.
"""
import numpy as np
from .mesh import Mesh, cross, publication_coordinates
from .quad_rotations import _PAIRS, _metrics


def _edges(faces):
    edges = np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2)
    return np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True)


def _centers(points, parity):
    origin = points.mean(axis=0)
    local = points-origin
    extent = float(np.linalg.norm(local, axis=1).max())
    triples = tuple(tuple((parity+2*k+j) % 6 for j in range(3)) for k in range(3))

    def jacobians(center):
        cells = np.asarray([[*local[list(ids)], center] for ids in triples])
        edges = np.roll(cells, -1, axis=1)-cells
        return cross(edges, np.roll(edges, -1, axis=1)).ravel()

    constant = jacobians(np.zeros(2))
    coefficients = np.c_[jacobians((1., 0.))-constant, jacobians((0., 1.))-constant]
    polygon = [np.array((-extent, -extent)), np.array((extent, -extent)),
               np.array((extent, extent)), np.array((-extent, extent))]
    # Each corner Jacobian is affine in the only new point. The clipped
    # polygon certifies all twelve corners together, rather than guessing a
    # center and hoping that a subsequent relaxation will unfold the mesh.
    for normal, value in zip(coefficients, constant-1e-7):
        clipped = []
        for a, b in zip(polygon, polygon[1:]+polygon[:1]):
            first, last = float(normal @ a+value), float(normal @ b+value)
            if first >= 0.:
                clipped.append(a)
            if (first >= 0.) != (last >= 0.):
                clipped.append(a+first/(first-last)*(b-a))
        polygon = clipped
        if not polygon:
            break
    return triples, None if len(polygon) < 3 else np.mean(polygon, axis=0)+origin


def positive_patches(mesh, *, max_cells=8000, max_candidates=4, include_stored_failures=False,
                     constrained_faces=()):
    """Return positive seeds; the consumer still owns final shape/height checks.

    Optional constrained faces name a separately proved fixed-height problem
    on a valid planar patch. No height limit is changed by this topology API.
    """
    if (type(max_cells) is not int or not 1 <= max_cells <= 50000
            or type(max_candidates) is not int or not 1 <= max_candidates <= 4
            or type(include_stored_failures) is not bool
            or not 1 <= len(mesh.quads) <= max_cells or mesh.uv_layers
            or not np.isfinite(mesh.xyz).all() or np.ptp(mesh.xyz[:, 2]) > 1e-10):
        raise ValueError('Positive patch recovery needs a bounded uncharted planar mesh')
    if (not isinstance(constrained_faces,tuple) or len(constrained_faces)>4
            or any(isinstance(face,(bool,np.bool_)) or not isinstance(face,(int,np.integer))
                   or not 0<=face<len(mesh.quads) for face in constrained_faces)
            or len(set(constrained_faces))!=len(constrained_faces)):
        raise ValueError('Positive patch constraints need at most four distinct existing faces')
    faces, xy = mesh.quads, mesh.xyz[:, :2]
    jacobians, _, _ = _metrics(xy, faces)
    bad = set(map(int, np.flatnonzero(jacobians <= 1e-10)))
    if include_stored_failures:
        stored = publication_coordinates(mesh.xyz, 1., physical=True)
        stored_jacobians, _, _ = _metrics(stored[:, :2], faces)
        bad.update(map(int, np.flatnonzero(stored_jacobians <= 1e-10)))
    bad.update(map(int,constrained_faces))
    if not 1 <= len(bad) <= 4:
        raise ValueError('Positive patch recovery needs one bounded singular patch')
    unique, inverse, counts = _edges(faces)
    if counts.max() > 2:
        raise ValueError('Positive patch recovery requires manifold shared edges')
    order = np.argsort(inverse, kind='stable')
    offsets = np.r_[0, np.cumsum(counts)]
    candidates = []

    def accept(xyz, quads, report):
        if len(quads) > max_cells:
            return
        candidate = Mesh(xyz, quads, mesh.material, dict(mesh.seams))
        try:
            quality = candidate.validate()
            candidate.with_positions(publication_coordinates(xyz, 1., physical=True)).validate()
        except ValueError:
            return
        actual, _, incidence = _edges(quads)
        added = len(xyz)-len(mesh.xyz)
        if (incidence.max() > 2 or not np.array_equal(actual[incidence == 1], unique[counts == 1])
                or len(actual) != len(unique)+2*added
                or not np.array_equal(xyz[:len(mesh.xyz)], mesh.xyz)):
            return
        details=dict(report, added_interior_vertices=added, exact_original_vertices=True,
                 exact_boundary_edges=True, complete_stored_mesh_positive=True,
                 triangles_created=False, shape_qualified=False,
                 admitted_stored_failures=include_stored_failures)
        if constrained_faces:details['admitted_constraint_faces']=constrained_faces
        candidates.append((quality['maximum_edge_ratio'],candidate,details))

    for edge in np.flatnonzero(counts == 2):
        first, last = map(int, order[offsets[edge]:offsets[edge+1]]//4)
        if not bad.intersection((first, last)) or bad.difference((first, last)):
            continue
        a, b = map(int, unique[edge])
        f, g = tuple(map(int, faces[first])), tuple(map(int, faces[last]))
        start = next(i for i in range(4) if {f[i], f[(i+1) % 4]} == {a, b})
        f = f[start:]+f[:start]
        a, b = f[:2]
        start = g.index(b)
        g = g[start:]+g[:start]
        if g[1] != a:
            raise ValueError('Positive patch shared edges need opposite orientations')
        ring = np.asarray((a, g[2], g[3], b, f[2], f[3]), dtype=np.int32)
        if len(set(ring)) != 6:
            continue
        for variant, indices in enumerate(_PAIRS):
            quads = faces.copy()
            quads[[first, last]] = ring[indices]
            accept(mesh.xyz, quads, dict(construction='two direct quad cells',
                                        original_faces=(first, last), variant=variant))
        for parity in (0, 1):
            triples, center = _centers(xy[ring], parity)
            if center is None:
                continue
            xyz = np.vstack((mesh.xyz, (*center, mesh.xyz[0, 2])))
            patch = np.asarray([(*ring[list(ids)], len(mesh.xyz)) for ids in triples], dtype=np.int32)
            quads = np.vstack((np.delete(faces, (first, last), axis=0), patch))
            accept(xyz, quads, dict(construction='three direct quad cells with certified interior point',
                                   original_faces=(first, last), parity=parity))
    if not candidates:
        raise ValueError('No bounded positive direct-quad patch exists')
    return tuple((candidate, report) for _, candidate, report in
                 sorted(candidates, key=lambda item: item[0])[:max_candidates])
