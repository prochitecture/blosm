"""Late fixed-boundary recovery for newly composed gapped Bundle grass.

The established positive stored mesh and height field return unchanged. A
stored-only singular two-cell patch may rotate its computational interior;
physical pavement edges, original vertices and source profiles stay fixed.
"""
import numpy as np
from .drainage import BoundaryDrainage
from .fixed_grass_heights import boundary_mask, reconcile
from .four_side_charts import quality
from .mesh import Mesh
from .positive_quad_patches import positive_patches
from .surface_query import QuadHeightField
from .publication_slopes import PublicationSlopeGuard


def recover(mesh, keys, field, *, namespace, max_cells=50000):
    mesh.validate()
    if len(keys) != len(mesh.xyz) or len(mesh.quads) > max_cells:
        raise ValueError('Gapped grass recovery lost its finite vertex identities')
    stored = Mesh(mesh.xyz.astype(np.float32).astype(float), mesh.quads, 'grass')
    topology = False
    try:
        stored.validate()
    except ValueError:
        topology = True
    if not topology:
        if max(quality(mesh)['maximum_edge_ratio'], quality(stored)['maximum_edge_ratio']) > 8.:
            return mesh, keys, None  # Preserve the existing shape failure.
        fixed = boundary_mask(mesh)
        drainage = BoundaryDrainage(max_cells=max_cells).fill(stored, fixed)[1]
        if (max(mesh.quad_slopes().max(), stored.quad_slopes().max()) <= .1+1e-7
                and drainage['maximum_fill'] <= 1e-6):
            return mesh, keys, None
        floors=PublicationSlopeGuard.fixed_gradient_floors(mesh,fixed,field.datum)
        constrained=tuple(map(int,np.flatnonzero(np.max(floors,axis=1)>.1+PublicationSlopeGuard.tolerance)))
        if constrained:
            # A three-boundary-vertex ear can force an excessive normal
            # gradient even when every physical edge slope is admissible.
            # Retile only its bounded two-cell computational patch; the
            # actual boundary heights and the ten-percent guard stay fixed.
            flat=mesh.with_positions(np.c_[mesh.xyz[:,:2],np.zeros(len(mesh.xyz))])
            candidates=positive_patches(flat,max_cells=max_cells,
                include_stored_failures=True,constrained_faces=constrained)
        else:
            candidates = ((mesh, None),)
    else:
        flat = mesh.with_positions(np.c_[mesh.xyz[:, :2], np.zeros(len(mesh.xyz))])
        candidates = positive_patches(flat, max_cells=max_cells, include_stored_failures=True)
    original_field = None
    failures = []
    for seed, patch in candidates:
        try:
            if patch is not None:
                seed.xyz[:len(mesh.xyz)] = mesh.xyz
                if len(seed.xyz) > len(mesh.xyz):
                    if original_field is None:
                        original_field = QuadHeightField(mesh, datum=field.datum,
                            ceiling=field.ceiling, max_cells=max_cells)
                    seed.xyz[len(mesh.xyz):, 2] = original_field.evaluate(seed.xyz[len(mesh.xyz):, :2])
            names = list(keys)+[(namespace, 'grass patch', i)
                               for i in range(len(mesh.xyz), len(seed.xyz))]
            output, query, physical, heights = reconcile(seed, datum=field.datum,
                ceiling=field.ceiling, max_cells=max_cells, policy='fixed_boundary_gapped_source_grass')
        except ValueError as error:
            failures.append(str(error))
            continue
        return output, names, dict(topology_patch=patch, height_reconciliation=heights,
            computational_internal_edges_may_change=patch is not None,
            physical_boundary_changed=False, original_profile_points_changed=False,
            prior_failures=failures)
    raise ValueError('No fixed-boundary gapped grass recovery qualifies: '+str(failures))
