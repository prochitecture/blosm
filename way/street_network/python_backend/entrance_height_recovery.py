"""Failed-only coupled height recovery with explicit entrance slope slabs."""
import numpy as np
from .drainage import BoundaryDrainage, ProtectedSpillUnavailable
from .height_slabs import HeightSlabs
from .slope_projection import ProjectionIterationBudgetError, SlopeProjection,corner_matrices
from .updates import drain_steps


def entrance_slabs(mesh, fixed, conditions, limit, *, at_corners=True):
    """Constrain both corner normal gradients along each real profile edge."""
    if type(at_corners) is not bool:
        raise ValueError('Entrance directional sampling requires an explicit Boolean')
    lookup={tuple(mesh.xyz[node,:2]):int(node) for node in np.flatnonzero(fixed)}
    edges={}
    for fi,face in enumerate(mesh.quads):
        for side in range(4):
            edges.setdefault(tuple(sorted((int(face[side]),int(face[(side+1)%4])))),[]).append((fi,side))
    matrices=corner_matrices(mesh)
    ids,coefficients=[],[]
    for boundary in conditions:
        for port in boundary.ports:
            nodes=[]
            for vertex in port.vertices:
                node=lookup.get(tuple(boundary.xy[vertex]))
                if node is None:
                    raise ValueError('Directional entrance recovery lost an original profile vertex')
                nodes.append(node)
            for first,last in zip(nodes,nodes[1:]):
                hits=edges.get(tuple(sorted((first,last))),())
                if len(hits) != 1:
                    raise ValueError('Directional entrance recovery lost a whole physical profile edge')
                fi,side=hits[0]
                tangent=mesh.xyz[last,:2]-mesh.xyz[first,:2]
                tangent/=np.linalg.norm(tangent)
                normal=np.array((tangent[1],-tangent[0]))
                if at_corners:
                    for corner in (side,(side+1)%4):
                        ids.append(mesh.quads[fi])
                        coefficients.append(normal@matrices[fi,corner])
                else:
                    # Match the authoritative incoming-edge midpoint check.
                    # Its bilinear XY Jacobian is not the average of the
                    # two corner inverse Jacobians on a skew quad.
                    p=mesh.xyz[mesh.quads[fi],:2]
                    u,v=((.5,0.),(1.,.5),(.5,1.),(0.,.5))[side]
                    a=(1.-v)*(p[1]-p[0])+v*(p[2]-p[3])
                    b=(1.-u)*(p[3]-p[0])+u*(p[2]-p[1])
                    determinant=a[0]*b[1]-a[1]*b[0]
                    if determinant<=1e-12:
                        raise ValueError('Entrance midpoint constraints need positive source quads')
                    du=np.array((-(1.-v),1.-v,v,-v))
                    dv=np.array((-(1.-u),-u,u,1.-u))
                    gradient=np.stack(((b[1]*du-a[1]*dv)/determinant,
                                       (a[0]*dv-b[0]*du)/determinant))
                    ids.append(mesh.quads[fi]);coefficients.append(normal@gradient)
    if not ids:
        raise ValueError('Directional entrance recovery has no original physical port')
    return HeightSlabs(np.asarray(ids),np.asarray(coefficients),np.full(len(ids),limit))


def combine_slabs(primary, additional):
    """Retain caller constraints alongside freshly measured entrance rows."""
    if not isinstance(primary,HeightSlabs) or (additional is not None and not isinstance(additional,HeightSlabs)):
        raise ValueError('Additional directional constraints require HeightSlabs')
    if additional is None:return primary
    return HeightSlabs(np.vstack((primary.vertices,additional.vertices)),
        np.vstack((primary.coefficients,additional.coefficients)),np.r_[primary.limits,additional.limits])


def _protected_drainage(candidate, body, fixed, drainage, slabs, field, physical_limit, tolerance):
    xyz = candidate.xyz.copy()
    xyz[body.core_nodes, 2] = body.core_targets
    protected = fixed.copy()
    protected[body.core_nodes] = True
    try:
        heights, details = drainage.fill_preserving(candidate.with_positions(xyz), fixed, protected)
    except ProtectedSpillUnavailable:
        return None
    if (heights.min() < field.datum or heights.max() > field.ceiling or
            slabs.violation(heights) > tolerance):
        return None
    gradients = np.einsum('fcij,fj->fci', corner_matrices(candidate), heights[candidate.quads])
    maximum = float(np.linalg.norm(gradients, axis=2).max())
    if maximum > physical_limit + tolerance:
        return None
    xyz[:, 2] = heights
    return candidate.with_positions(xyz), dict(details, maximum_slope=maximum)


def recover(builder, mesh, field, fixed, conditions, physical_limit, max_cells, *, at_corners=True,
            additional_slabs=None, body_guides=None):
    return drain_steps(recover_steps(builder, mesh, field, fixed, conditions, physical_limit, max_cells,
        at_corners=at_corners, additional_slabs=additional_slabs, body_guides=body_guides))


def recover_steps(builder, mesh, field, fixed, conditions, physical_limit, max_cells, *, at_corners=True,
                  additional_slabs=None, body_guides=None):
    """Return only after normal envelope, spill, profile and slope checks pass."""
    target=.95*builder.maximum_join_slope
    slabs=combine_slabs(entrance_slabs(mesh,fixed,conditions,target,at_corners=at_corners),additional_slabs)
    floors=np.full(len(mesh.xyz),field.datum)
    projector=SlopeProjection(max_cells=max_cells)
    drainage=BoundaryDrainage(max_cells=max_cells)
    candidate=mesh
    rounds=[]
    for _ in range(8):
        try:
            candidate,projection=yield from projector.solve_steps(candidate,fixed,datum=field.datum,ceiling=field.ceiling,
                maximum_slope=physical_limit,minimum_heights=floors,height_slabs=slabs)
        except ProjectionIterationBudgetError as failure:
            if body_guides is None or not np.any(floors[body_guides.core_nodes] > body_guides.core_targets):
                raise
            proposal = _protected_drainage(failure.candidate, body_guides, fixed, drainage,
                slabs, field, physical_limit, projector.slope_tolerance)
            if proposal is None:
                raise
            candidate, outlet = proposal
            rounds.append(dict(previous_projection=str(failure), protected_drainage=outlet))
            break
        filled,spill=drainage.fill(candidate,fixed)
        rounds.append(dict(projection=projection,drainage=spill))
        if spill['maximum_fill'] <= 1e-6:break
        if body_guides is not None and np.any(filled[body_guides.core_nodes] > body_guides.core_targets + 1e-8):
            proposal = _protected_drainage(candidate, body_guides, fixed, drainage,
                slabs, field, physical_limit, projector.slope_tolerance)
            if proposal is not None:
                candidate, outlet = proposal
                rounds[-1]['protected_drainage'] = outlet
                break
        raised=filled > candidate.xyz[:,2]+1e-10
        floors[raised]=np.maximum(floors[raised],filled[raised])
    else:
        raise ValueError('Directional entrance recovery exceeded its coupled spill budget')
    checks=builder.projection_checks(candidate,fixed,conditions)
    if (checks['isolated_interior_minima'] or checks['maximum_entrance_longitudinal_slope'] > builder.maximum_join_slope
            or slabs.violation(candidate.xyz[:,2]) > projector.slope_tolerance
            or not np.array_equal(candidate.xyz[fixed],mesh.xyz[fixed])
            or not np.array_equal(candidate.xyz[:,:2],mesh.xyz[:,:2])
            or not np.array_equal(candidate.quads,mesh.quads)):
        raise ValueError('Directional entrance recovery failed its unchanged physical contract')
    if additional_slabs is not None:
        candidate=mesh.with_positions(candidate.xyz)
    return candidate,dict(construction='fixed-chart entrance normal slabs with slope balls and spill floors',
        rounds=rounds,joins=checks,directional_target=target,physical_slope_limit=physical_limit,
        exact_fixed_xyz=True,exact_xy_and_topology=True,
        **(dict(directional_samples='profile_edge_midpoints') if not at_corners else {}),
        **(dict(additional_height_rows=len(additional_slabs.vertices),
            additional_height_violation=additional_slabs.violation(candidate.xyz[:,2]))
            if additional_slabs is not None else {}))
