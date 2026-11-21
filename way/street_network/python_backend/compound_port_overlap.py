"""Resolve intersecting entrance profiles through an existing source junction."""
import numpy as np
from .mesh import cross
from .model import End


def shared_port_neighbors(planner, members, layout):
    """An actual common far junction can be required by crossed full mouths.

    Two approaches may converge at a second original junction outside the
    current compound. If their complete transverse profiles intersect, both
    cannot be independent exposed entrances. Include their recorded common
    neighbor; equal coordinates alone never establish that ownership.
    """
    inventory=[]
    for approach,station in zip(layout.approaches,layout.outward_stations):
        route=planner.topology.trace(End(approach.source.id,approach.end))
        if (route.terminal_site not in planner.sites or route.terminal_site in members
                or route.terminal_scope or route.terminal_selection
                or getattr(approach.path,'absorbed_sites',())
                or tuple(p.section.id for p in approach.path.pieces)
                   !=tuple(v.section_id for v in route.traversals)):
            inventory.append(None);continue
        parameters=getattr(approach,'mouth_parameters',None)
        u=(parameters(layout.radial,approach.width/(2*layout.radial)) if parameters
           else np.linspace(1.,0.,2*layout.radial+1))
        xy=approach.point(station,approach.offset+approach.width*(u-.5))
        inventory.append((route.terminal_site,xy))
    neighbors=set()
    for i,first in enumerate(inventory):
        if first is None:continue
        for second in inventory[i+1:]:
            if second is None or first[0]!=second[0]:continue
            if _intersects(first[1],second[1]):neighbors.add(first[0])
    return neighbors


def _intersects(first, second, tolerance=1e-7):
    a,b=first[:-1],second[:-1]
    u,v=np.diff(first,axis=0),np.diff(second,axis=0)
    lu,lv=np.linalg.norm(u,axis=1),np.linalg.norm(v,axis=1)
    determinant=cross(u[:,None],v[None])
    delta=b[None]-a[:,None]
    valid=abs(determinant)>1e-12*lu[:,None]*lv[None]
    t=np.divide(cross(delta,v[None]),determinant,out=np.zeros_like(determinant),where=valid)
    s=np.divide(cross(delta,u[:,None]),determinant,out=np.zeros_like(determinant),where=valid)
    if np.any(valid & (t>=-tolerance/lu[:,None]) & (t<=1.+tolerance/lu[:,None])
            & (s>=-tolerance/lv[None]) & (s<=1.+tolerance/lv[None])):return True
    # Parallel or endpoint contacts also prevent two separate exposed mouths.
    for points,start,edge,length in ((first,b,v,lv),(second,a,u,lu)):
        offset=points[:,None]-start[None]
        q=np.clip(np.sum(offset*edge[None],axis=-1)/(length[None]**2),0.,1.)
        if np.any(np.linalg.norm(offset-q[...,None]*edge[None],axis=-1)<=tolerance):return True
    return False
