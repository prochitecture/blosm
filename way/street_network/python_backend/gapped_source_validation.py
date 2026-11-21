"""Whole original-source coverage of a new body and its actual query owners."""
from types import SimpleNamespace
import numpy as np
from .bundle_merged_body import MergedRunBody
from .regions import signed_area
from .source_strip_containment import certify


def mesh_union_loops(meshes):
    """Cancel only identical whole interfaces; construct no geometric repair."""
    edges = {}
    for mesh in meshes:
        boundary = {tuple(edge) for edge in MergedRunBody._boundary_vertices(mesh)}
        for face in mesh.quads:
            for a, b in zip(face, np.roll(face, -1)):
                if tuple(sorted((int(a), int(b)))) not in boundary: continue
                first, last = tuple(mesh.xyz[a, :2]), tuple(mesh.xyz[b, :2])
                token = tuple(sorted((first, last)))
                edges.setdefault(token, []).append((first, last))
    directed = {}
    for pairs in edges.values():
        if len(pairs) == 2:
            if pairs[0] != pairs[1][::-1]:
                raise ValueError('Source coverage has two equally directed overlapping physical edges')
            continue
        if len(pairs) != 1: raise ValueError('Source coverage has a nonmanifold whole physical interface')
        first, last = pairs[0]
        if first in directed: raise ValueError('Source coverage has an unresolved physical boundary branch')
        directed[first] = last
    loops = []
    while directed:
        first = min(directed); at = first; points = []
        while True:
            points.append(at)
            if at not in directed: raise ValueError('Source coverage has an open whole-edge boundary')
            at = directed.pop(at)
            if at == first: break
            if len(points) > 20000: raise ValueError('Source coverage exceeds its finite boundary allocation')
        xy = np.asarray(points)
        if len(xy) < 3 or abs(signed_area(xy)) < 1e-10:
            raise ValueError('Source coverage has a collapsed physical loop')
        loops.append(SimpleNamespace(xy=xy, hole=signed_area(xy) < 0.))
    return tuple(loops)


def complete_width(context, inventory, pavement, owner_meshes, *, maximum_cells=20000):
    """Every original member Section, including whole query-owned intervals."""
    loops = mesh_union_loops((pavement, *owner_meshes))
    sections = tuple(v.section_id for member in inventory.run.members for v in member.visits)
    if len(sections) != len(set(sections)):
        raise ValueError('Gapped source coverage repeats an original member Section')
    intervals = {row['section']: (row['first'], row['last'])
                 for own in inventory.spans_by_role for span in own for row in span.source_intervals}
    for span in inventory.owned_spans:
        intervals.update((row['section'], (row['first'], row['last'])) for row in span.intervals)
    proof = certify(context, SimpleNamespace(loops=loops), sections,
        intervals=intervals, maximum_cells=maximum_cells)
    return dict(proof, complete_original_source_ids=list(sections),
        interval_scope='visible original intervals plus separately certified wholly owned Sections',
        physical_interfaces='identical whole mesh edges only', no_geometric_union_repair=True)
