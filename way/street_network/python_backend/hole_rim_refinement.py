"""Conforming 1:3 refinement on owned separator-boundary chords only."""
import numpy as np
from .compounds import rim_label
from .partition import BoundaryRef
from .mesh import Assembler
from .updates import drain_steps


def grade_hole_mesh(mesh, refs, group, patterns, *, maximum_ratio=4., max_cells=50000, max_rounds=8):
    return drain_steps(grade_hole_mesh_steps(mesh, refs, group, patterns, maximum_ratio=maximum_ratio,
        max_cells=max_cells, max_rounds=max_rounds))


def grade_hole_mesh_steps(mesh, refs, group, patterns, *, maximum_ratio=4., max_cells=50000, max_rounds=8):
    sources = {p.region.owner: p for p in group.footprints}
    history = []
    guides = set(range(len(mesh.xyz)))
    for _ in range(max_rounds):
        split = set()
        for li, ids in mesh.seams.items():
            if not group.loops[li].hole:
                continue
            points = mesh.xyz[ids, :2]
            lengths = np.linalg.norm(np.roll(points, -1, axis=0)-points, axis=1)
            for i in np.flatnonzero(lengths > maximum_ratio*np.minimum(np.roll(lengths, 1), np.roll(lengths, -1))):
                first, last = int(ids[i]), int(ids[(i+1) % len(ids)])
                ref = refs[first]
                if rim_label(group.loops[ref.loop].spans[ref.edge], sources).kind != 'curb':
                    raise ValueError('Separator refinement cannot insert an incoming profile sample')
                split.add(tuple(sorted((first, last))))
        if not split:
            return mesh, refs, dict(rounds=history, exact_polygon=True, prescribed_profiles_unchanged=True), np.array(sorted(guides), dtype=np.int32)
        affected = set()
        atlases, total = [], 0
        for index, quad in enumerate(mesh.quads):
            if index % 128 == 0:
                yield
            sides = tuple(3 if tuple(sorted((int(quad[i]), int(quad[(i+1)%4])))) in split else 1
                for i in range(4))
            if 3 in sides:
                affected.update(map(int, quad))
            pattern = patterns.quad(*sides)
            total += len(pattern.quads)
            if total > max_cells:
                raise ValueError('Separator-rim refinement exceeds its cell budget')
            atlases.append(pattern)
        assembly = Assembler()
        for fi, (quad, pattern) in enumerate(zip(mesh.quads, atlases)):
            if fi % 128 == 0:
                yield
            if affected.intersection(map(int, quad)):
                guides.difference_update(map(int, quad))
            u, v = pattern.uv.T
            p = mesh.xyz[quad]
            xyz = (((1.-u)*(1.-v))[:, None]*p[0]+(u*(1.-v))[:, None]*p[1]
                   +(u*v)[:, None]*p[2]+((1.-u)*v)[:, None]*p[3])
            keys = []
            for i, (x, y) in enumerate(pattern.uv):
                if x in (0., 1.) and y in (0., 1.):
                    corner = {(0.,0.):0,(1.,0.):1,(1.,1.):2,(0.,1.):3}[x,y]
                    keys.append(('vertex', int(quad[corner])))
                    xyz[i] = mesh.xyz[quad[corner]]
                else:
                    side, fraction = ((0,x) if y==0. else (1,y) if x==1. else
                                      (2,1.-x) if y==1. else (3,1.-y) if x==0. else (None,None))
                    if side is None:
                        keys.append(('face', fi, i))
                    else:
                        a, b = int(quad[side]), int(quad[(side+1) % 4])
                        t = fraction if a < b else 1.-fraction
                        keys.append(('edge', min(a,b), max(a,b), int(round(t*3.))))
            assembly.add(xyz, pattern.quads, keys)
        output = assembly.mesh(mesh.material)
        mapping = {i:assembly.lookup['vertex', i] for i in range(len(mesh.xyz))}
        guides = {mapping[i] for i in guides}
        updated = {mapping[i]:ref for i,ref in refs.items()}
        seams = {}
        for li, ids in mesh.seams.items():
            values = []
            for a, b in zip(ids, np.roll(ids, -1)):
                a, b = int(a), int(b)
                values.append(mapping[a])
                if tuple(sorted((a,b))) not in split:
                    continue
                first, last = refs[a], refs[b]
                end = last.fraction if first.edge == last.edge else 1.
                for part in (1,2):
                    index = part if a < b else 3-part
                    node = assembly.lookup['edge', min(a,b), max(a,b), index]
                    fraction = first.fraction+part/3.*(end-first.fraction)
                    updated[node] = BoundaryRef(first.loop, first.edge, fraction)
                    values.append(node)
            seams[li] = np.asarray(values, dtype=np.int32)
        output.seams = seams
        output.validate()
        history.append(dict(split_chords=len(split), quads=len(output.quads)))
        mesh, refs = output, updated
    raise ValueError('Separator-rim refinement exceeds its propagation budget')
