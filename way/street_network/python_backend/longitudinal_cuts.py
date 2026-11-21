"""Direct quad blocks along elongated physical pavement regions.

The cuts are computational. The sampled analytical boundary, incoming profile
rows and source ownership remain authoritative. A short cap is split locally
behind a reflex nose; long sides become four-sided strip charts.
"""
from dataclasses import replace
import math
import numpy as np
from .compounds import rim_label
from .mesh import cross
from .partition import (BoundaryOwnershipCuts, BoundaryRef, BoundarySpan,
                        JunctionCut, PartitionEdge)


class LongitudinalCuts(BoundaryOwnershipCuts):
    def __init__(self, spacing, *, maximum_cuts=64, minimum_aspect=6., all_noses=False):
        if (not math.isfinite(spacing) or spacing <= 0.
                or type(maximum_cuts) is not int or maximum_cuts < 1
                or not math.isfinite(minimum_aspect) or minimum_aspect<2.
                or type(all_noses) is not bool):
            raise ValueError('Invalid longitudinal chart budgets')
        self.spacing, self.maximum_cuts = float(spacing), maximum_cuts
        self.minimum_aspect=float(minimum_aspect)
        self.all_noses=all_noses

    @staticmethod
    def eligible(group, spacing):
        if len(group.loops) != 1 or group.loops[0].hole:
            return False
        points = np.asarray(group.loops[0].xy, dtype=float)
        local = points-points.mean(axis=0)
        _, vectors = np.linalg.eigh(local.T @ local)
        extents = np.ptp(local @ vectors, axis=0)
        width, length = map(float, extents)
        return width > 1e-7 and length >= max(6.*width, 8.*spacing)

    @staticmethod
    def reference(edge, fraction, count):
        # Splitting a physical chord at an arbitrary near-end fraction can
        # create a tiny boundary edge. Existing samples carry the same curve
        # accuracy and give a better-sized endpoint for a computational cut.
        return BoundaryRef(0, (int(edge) + int(fraction > .5)) % count, 0.)

    def cuts(self, group, context):
        if len(group.loops) != 1 or group.loops[0].hole:
            raise ValueError('Longitudinal partition requires one exterior pavement loop')
        loop = group.loops[0]
        points = np.asarray(loop.xy, dtype=float)
        center = points.mean(axis=0)
        _, vectors = np.linalg.eigh((points-center).T @ (points-center))
        axis = vectors[:, -1]
        if axis[int(np.argmax(abs(axis)))] < 0.:
            axis = -axis
        normal = np.array((-axis[1], axis[0]))
        along, across = (points-center) @ axis, (points-center) @ normal
        length, width = float(np.ptp(along)), float(np.ptp(across))
        if width <= 1e-7 or length < max(self.minimum_aspect*width, 8.*self.spacing):
            raise ValueError('Physical footprint is not an elongated corridor')
        count = max(2, int(math.ceil(length/max(width, 4.*self.spacing))))
        if count > self.maximum_cuts:
            raise ValueError('Longitudinal partition exceeds its cut budget')
        sources = {p.region.owner: p for p in group.footprints}
        labels = [rim_label(span, sources) for span in loop.spans]
        delta = np.roll(points, -1, axis=0)-points
        projected = np.roll(along, -1)-along
        # Keep semantic ownership cuts when they still define the boundary.
        # A certified replacement exterior has no independent star-kernel
        # requirement: each mapped quad is certified instead.
        cuts = list(super().cuts(group, context)) if group.claims and group.physical_layout is None else []
        additional = []
        for station in np.linspace(along.min(), along.max(), count+1)[1:-1]:
            # A half-open interval makes a boundary vertex one incidence.
            selected = ((along <= station) & (station < np.roll(along, -1))) | (
                        (np.roll(along, -1) <= station) & (station < along))
            hits = []
            for edge in np.flatnonzero(selected):
                fraction = float((station-along[edge])/projected[edge])
                position = points[edge]+fraction*delta[edge]
                hits.append((float((position-center) @ normal), int(edge), fraction))
            hits.sort(key=lambda item: item[0])
            if len(hits) % 2:
                raise ValueError('Longitudinal partition has an unresolved boundary crossing')
            for first, last in zip(hits[::2], hits[1::2]):
                if any(labels[hit[1]].kind != 'curb' for hit in (first, last)):
                    continue
                if last[0]-first[0] <= 4.*self.spacing:
                    continue
                if any(abs(projected[hit[1]]) < .2*np.linalg.norm(delta[hit[1]]) for hit in (first, last)):
                    continue  # Avoid a cut tangent to a physical curb.
                additional.append(JunctionCut(None, group.members[0], group.members[0],
                    self.reference(first[1], first[2], len(points)),
                    self.reference(last[1], last[2], len(points))))
        additional.extend(self.nose_cuts(points, delta, along, axis, normal, labels, group.members[0]))
        if not additional:
            raise ValueError('Elongated pavement has no transverse cuts on owned curbs')
        if len(cuts)+len(additional) > self.maximum_cuts:
            raise ValueError('Longitudinal partition exceeds its cut budget')
        return tuple(cuts+additional)

    def nose_cuts(self, points, delta, along, axis, normal, labels, owner):
        unit = delta/np.linalg.norm(delta, axis=1)[:, None]
        turns = np.arctan2(cross(np.roll(unit, 1, axis=0), unit),
                           np.sum(np.roll(unit, 1, axis=0)*unit, axis=1))
        clusters = []
        for vertex in np.flatnonzero(turns < -1e-6):
            if not clusters or vertex != clusters[-1][-1]+1:
                clusters.append([])
            clusters[-1].append(int(vertex))
        if not clusters:
            return
        if len(clusters) > 1 and clusters[0][0] == 0 and clusters[-1][-1] == len(points)-1:
            clusters[0] = clusters.pop()+clusters[0]
        selected_clusters=(clusters if self.all_noses else
            (max(clusters,key=lambda indices:-float(turns[indices].sum())),))
        for selected in selected_clusters:
            yield from self._nose_cluster(points,delta,along,axis,normal,labels,owner,turns,selected)

    def _nose_cluster(self,points,delta,along,axis,normal,labels,owner,turns,selected):
        if -float(turns[selected].sum()) < math.pi/3.:
            return  # A gentle bend does not define a fork nose.
        orientation = 1. if np.mean(along[selected]) > .5*(along.min()+along.max()) else -1.
        outward = orientation*axis
        index = min(selected, key=lambda i: orientation*along[i])
        origin = points[index]
        for direction in (-outward+normal, -outward-normal):
            broadcast = np.broadcast_to(direction, delta.shape)
            denom = cross(broadcast, delta)
            offset = points-origin
            valid = abs(denom) > 1e-12
            t = np.divide(cross(offset, delta), denom,
                          out=np.full(len(points), np.inf), where=valid)
            u = np.divide(cross(offset, broadcast), denom,
                          out=np.full(len(points), np.inf), where=valid)
            hits = np.flatnonzero((t > 1e-6) & (u >= 0.) & (u <= 1.))
            if not len(hits):
                continue
            edge = int(hits[np.argmin(t[hits])])
            if labels[edge].kind != 'curb':
                continue
            yield JunctionCut(None, owner, owner, BoundaryRef(0, index, 0.),
                              self.reference(edge, u[edge], len(points)))


def even_sided_partition(partition, *, max_nodes=16384):
    """Refine owned curb chords so structural atlas sides have even parity.

    Opposite sides of the cached four-sided patterns require equal parity.
    Every cut interval can then use an even count. Existing physical vertices
    and incoming profiles are untouched; a split adds a midpoint to an owned
    curb chord, leaving the complete physical polygon exactly unchanged.
    """
    points = list(partition.xy.copy())
    references = dict(partition.boundary_refs)
    sources = {p.region.owner: p for p in partition.group.footprints}
    marked = set()
    for edge in partition.edges:
        if edge.cut is not None:
            marked.update((edge.first, edge.last))
    for loop in partition.boundary_loops:
        p = partition.xy[list(loop)]
        e = np.roll(p, -1, axis=0)-p
        u = e/np.linalg.norm(e, axis=1)[:, None]
        turn = np.arctan2(cross(np.roll(u, 1, axis=0), u), np.sum(np.roll(u, 1, axis=0)*u, axis=1))
        marked.update(loop[i] for i in np.flatnonzero(abs(turn) > .35))
    split = set()
    for block in partition.blocks:
        starts = [i for i, node in enumerate(block.vertices) if node in marked]
        if len(starts) < 2:
            continue
        n = len(block.edges)
        for a, b in zip(starts, starts[1:]+[starts[0]+n]):
            side = [block.edges[i % n][0] for i in range(a, b)]
            if len(side) % 2 == 0 or any(partition.edges[i].cut is not None for i in side):
                continue
            candidates = []
            for ei in side:
                edge = partition.edges[ei]
                span = edge.boundary
                if rim_label(partition.group.loops[span.loop].spans[span.edge], sources).kind == 'curb':
                    candidates.append((np.linalg.norm(partition.xy[edge.last]-partition.xy[edge.first]), ei))
            if not candidates:
                raise ValueError('Physical profile has incompatible four-sided atlas parity')
            split.add(max(candidates)[1])
    if len(points)+len(split) > max_nodes:
        raise ValueError('Longitudinal parity refinement exceeds its graph-node budget')
    edges, mapping, middle = [], {}, {}
    for ei, edge in enumerate(partition.edges):
        if ei not in split:
            mapping[ei] = (len(edges),)
            edges.append(edge)
            continue
        node = len(points)
        points.append(.5*(partition.xy[edge.first]+partition.xy[edge.last]))
        span = edge.boundary
        t = .5*(span.first+span.last)
        references[node] = BoundaryRef(span.loop, span.edge, t)
        mapping[ei] = (len(edges), len(edges)+1)
        middle[(edge.first, edge.last)] = middle[(edge.last, edge.first)] = node
        edges.extend((PartitionEdge(edge.first, node, BoundarySpan(span.loop, span.edge, span.first, t)),
                      PartitionEdge(node, edge.last, BoundarySpan(span.loop, span.edge, t, span.last))))
    blocks = []
    for block in partition.blocks:
        labels, vertices = [], []
        for ei, forward in block.edges:
            ids = mapping[ei] if forward else mapping[ei][::-1]
            for index in ids:
                edge = edges[index]
                labels.append((index, forward))
                vertices.append(edge.first if forward else edge.last)
        blocks.append(replace(block, vertices=tuple(vertices), edges=tuple(labels)))
    loops = []
    for loop in partition.boundary_loops:
        values = []
        for a, b in zip(loop, loop[1:]+loop[:1]):
            values.append(a)
            if (a, b) in middle:
                values.append(middle[a, b])
        loops.append(tuple(values))
    return replace(partition, xy=np.asarray(points), edges=tuple(edges), blocks=tuple(blocks),
                   boundary_refs=references, boundary_loops=tuple(loops))
