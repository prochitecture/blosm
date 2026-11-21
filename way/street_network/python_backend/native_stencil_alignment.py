"""Align new quad vertices on one analytically verified carriageway profile."""
import numpy as np

from .four_side_charts import quality
from .intersections import Approach
from .local_quad_grading import _minimize, _objective
from .mesh import cross, publication_coordinates


def _grade_penalty(xy, heights, tangent, faces, limit):
    """Differentiate the four corner plane gradients as their XY and prescribed heights move."""
    points, z = xy[faces], heights[faces]
    ahead = np.roll(points,-1,axis=1)-points
    behind = np.roll(points,1,axis=1)-points
    det = cross(ahead,behind)
    if np.any(det<=1e-12):
        return 1e6,np.zeros_like(xy)
    za = np.roll(z,-1,axis=1)-z
    zb = np.roll(z,1,axis=1)-z
    g = np.stack((za*behind[:,:,1]-zb*ahead[:,:,1],zb*ahead[:,:,0]-za*behind[:,:,0]),axis=-1)/det[:,:,None]
    length = np.maximum(np.linalg.norm(g,axis=2),1e-14)
    excess = np.maximum(length/limit-.94,0.)
    # If A g = dz then A dg = d(dz) - dA g. The two weights solve
    # the transposed system once per corner for its scalar slope penalty.
    derivative = (2*excess/(limit*length))[:,:,None]*g
    wa = cross(derivative,behind)/det
    wb = cross(ahead,derivative)/det
    t = tangent[faces]
    contributions = np.roll(wa[:,:,None]*(np.roll(t,-1,axis=1)-g),1,axis=1)
    contributions += np.roll(wb[:,:,None]*(np.roll(t,1,axis=1)-g),-1,axis=1)
    contributions -= (wa+wb)[:,:,None]*(t-g)
    gradient = np.column_stack([np.bincount(faces.ravel(),weights=part.ravel(),minlength=len(xy))
                                for part in (contributions[:,:,0],contributions[:,:,1])])
    return float(np.sum(excess**2)),gradient


def align(mesh, original_nodes, fixed, body, reference, limit):
    """Keep every existing guide fixed. Move only added vertices belonging to this one exact body."""
    added = np.ones(len(mesh.xyz),bool)
    added[original_nodes] = False
    if not 0 < np.count_nonzero(added) <= 512:
        raise ValueError('Native alignment needs a bounded local conforming refinement')
    affected = mesh.quads[added[mesh.quads].any(axis=1)]
    known = np.zeros(len(mesh.xyz),bool)
    known[body.core_nodes] = True
    if len(affected)>2048 or not known[affected].all() or fixed[affected].any():
        raise ValueError('Native alignment requires original body heights on the complete local patch')
    nodes,indices = np.unique(affected,return_inverse=True)
    faces = indices.reshape(-1,4)
    mobile = np.flatnonzero(added[nodes])
    xy = mesh.xyz[nodes,:2].copy()
    base_z = mesh.xyz[nodes,2].copy()
    arm = Approach(reference.path,0,reference.profile)
    station,q,valid,_ = reference.inverse(xy)
    if not valid.all() or np.max(abs(reference.profile.datum+reference.profile.value(q)-base_z))>1e-8:
        raise ValueError('Native alignment cannot span multiple source profiles')
    scale = np.full(len(mobile),np.median(np.linalg.norm(np.roll(xy[faces],-1,axis=1)-xy[faces],axis=2)))
    tangent = np.zeros_like(xy)
    for i in mobile:
        tangent[i] = arm.physical_gradient(station[i],q[i])

    def evaluate(offsets):
        value,shape_gradient = _objective(xy,faces,mobile,scale,offsets,.06)
        trial = xy.copy()
        trial[mobile] += offsets.reshape(-1,2)*scale[:,None]
        local_z = base_z+np.sum(tangent*(trial-xy),axis=1)
        cost,gradient = _grade_penalty(trial,local_z,tangent,faces,limit)
        return value+cost,shape_gradient+(gradient[mobile]*scale[:,None]).ravel()

    offsets,steps = _minimize(evaluate,2*len(mobile))
    xy[mobile] += offsets.reshape(-1,2)*scale[:,None]
    xyz = mesh.xyz.copy()
    xyz[nodes[mobile],:2] = xy[mobile]
    # Use the original nonlinear field for every final height. The local
    # tangent plane above is only a search direction.
    xyz[nodes[mobile],2] = reference.values(xy[mobile])
    candidate = mesh.with_positions(xyz)
    physical = quality(candidate)
    stored = quality(candidate.with_positions(publication_coordinates(xyz,1.,physical=True)))
    if (max(physical['maximum_edge_ratio'],stored['maximum_edge_ratio'])>8.
            or min(physical['minimum_corner_sine'],stored['minimum_corner_sine'])<.03
            or not np.array_equal(xyz[original_nodes],mesh.xyz[original_nodes])):
        raise ValueError('Native alignment did not qualify original-guide and stored quad geometry')
    return candidate,dict(iterations=steps,movable_vertices=len(mobile),local_faces=len(faces),
                          physical_quality=physical,stored_quality=stored)
