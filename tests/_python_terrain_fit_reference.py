"""Frozen dev78 terrain fit for exact cache/invalidation comparisons."""
import heapq
import numpy as np
from tested_python_streets.terrain.alignment import PlacementSettings, HeightCurve, PlacementResult, quad_samples

def fit(references, components, terrain, settings=PlacementSettings()):
    """terrain(xy[N,2]) returns original-terrain Z[N], or raises on any miss."""
    settings.validate()
    refs = {r.identity: r for r in references}
    if len(refs) != len(references):
        raise ValueError('Duplicate elevation reference')
    if len({c.identity for c in components}) != len(components):
        raise ValueError('Duplicate corridor component')
    if any(c.reference not in refs for c in components):
        raise ValueError('Corridor component refers to an unknown elevation reference')
    node_ids, parent, floors = {}, [], []
    for key, ref in refs.items():
        node_ids[key] = np.arange(len(parent), len(parent)+len(ref.stations))
        parent.extend(range(len(parent), len(parent)+len(ref.stations)))
        floors.extend([-np.inf]*len(ref.stations))
    parent, floors = np.asarray(parent), np.asarray(floors)
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    def join(a, b):
        a, b = root(a), root(b)
        parent[max(a,b)] = min(a,b)
    contacts, mappings, sample_count = {}, {}, 0
    for c in components:
        ref = refs[c.reference]
        if np.min(c.stations) < ref.stations[0]-1e-6 or np.max(c.stations) > ref.stations[-1]+1e-6:
            raise ValueError('Component stations lie outside its elevation reference: '+c.identity)
        for s, contact in c.contacts:
            j = int(np.argmin(np.abs(ref.stations-s)))
            if abs(ref.stations[j]-s) > 1e-6:
                raise ValueError('A contact must lie on an elevation control')
            node = node_ids[c.reference][j]
            if contact in contacts:
                join(node, contacts[contact])
            contacts[contact] = node
        lo = np.clip(np.searchsorted(ref.stations, c.stations, side='right')-1, 0, max(0,len(ref.stations)-2))
        hi = np.minimum(lo+1, len(ref.stations)-1)
        mappings[c.identity] = (node_ids[c.reference][lo], node_ids[c.reference][hi])
        for ids, samples in quad_samples(c.xyz, c.quads, settings.sample_spacing):
            height = np.asarray(terrain(samples[:,:,:2].reshape(-1,2))).reshape(samples.shape[:2])
            if not np.isfinite(height).all():
                raise ValueError('Terrain does not cover the complete corridor: '+c.identity)
            required = np.max(height - (samples[:,:,2]-c.datum), axis=1) + settings.clearance+c.clearance
            # Raise all controls influencing every vertex of this face. The
            # final face, not only the analytical reference, meets these bounds.
            for nodes in mappings[c.identity]:
                np.maximum.at(floors, nodes[c.quads[ids]].ravel(), np.repeat(required, 4))
            sample_count += height.size
    roots = np.array([root(i) for i in range(len(parent))])
    merged = np.full(len(parent), -np.inf)
    np.maximum.at(merged, roots, floors)
    edges = [[] for _ in parent]
    for key, ref in refs.items():
        ns = roots[node_ids[key]]
        for a,b,distance in zip(ns[:-1], ns[1:], np.diff(ref.stations)):
            # Capped Hermite derivatives are bounded by twice the largest
            # adjacent secant. This reserve guarantees the requested grade.
            cost = .5*settings.maximum_grade*distance
            if a != b:
                edges[a].append((b,cost)); edges[b].append((a,cost))
    queue = [(-h,i) for i,h in enumerate(merged) if np.isfinite(h)]
    heapq.heapify(queue)
    while queue:
        negative, a = heapq.heappop(queue)
        if -negative < merged[a]-1e-12:
            continue
        for b,cost in edges[a]:
            value = merged[a]-cost
            if value > merged[b]+1e-12:
                merged[b] = value
                heapq.heappush(queue,(-value,b))
    curves = {key: HeightCurve(ref.stations, merged[roots[node_ids[key]]]) for key,ref in refs.items()}
    if any(not np.isfinite(c.heights).all() for c in curves.values()):
        raise ValueError('An elevation reference has no terrain constraints')
    grade = max((c.maximum_grade() for c in curves.values()), default=0.)
    if grade > settings.maximum_grade+1e-9:
        raise ValueError('Terrain reference exceeds the requested longitudinal grade')
    coordinates, maximum_lift = {}, 0.
    for c in components:
        xyz = c.xyz.copy()
        xyz[:,2] += curves[c.reference](c.stations)-c.datum
        ground = np.asarray(terrain(xyz[:,:2]))
        lift = float(np.max(xyz[:,2]-ground))
        maximum_lift = max(maximum_lift, lift)
        if lift > settings.maximum_lift + settings.clearance+c.clearance:
            raise ValueError(f'{c.identity}: terrain fit lifts the surface {lift:.2f} m above ground; limit is {settings.maximum_lift:.2f} m. Reduce the protected region or revise the alignment/terrain settings')
        coordinates[c.identity] = xyz
    return PlacementResult(curves, coordinates, dict(references=len(refs), samples=sample_count,
        maximum_reference_grade=grade, maximum_surface_lift=maximum_lift,
        sample_spacing=settings.sample_spacing, clearance=settings.clearance,
        policy='sampled no-cut; unchanged original terrain'))
