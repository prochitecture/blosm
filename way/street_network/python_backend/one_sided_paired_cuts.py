"""Opt-in computational seams for one-sided paired physical junctions.

Existing merged-mouth peak vertices separate the two source half-profiles.
They remain one authoritative profile after shared quad assembly. A genuine
hole receives two connections to its existing longitudinal tips; no physical
boundary or incoming profile sample is added or moved.
"""
import numpy as np
from .compounds import rim_label
from .partition import BoundaryRef, CutPolicy, JunctionCut


class OneSidedPairedCuts(CutPolicy):
    def __init__(self, source_intervals=(), *, existing_vertices=False):
        self.source_intervals = tuple(source_intervals)
        self.existing_vertices = bool(existing_vertices)

    def cuts(self, group, context):
        if not getattr(group.physical_layout, 'one_sided_source_preparation', False):
            raise ValueError('Paired physical source cuts require explicit one-sided preparation')
        holes = [(i, loop) for i, loop in enumerate(group.loops) if loop.hole]
        if len(holes) > 1:
            raise ValueError('Paired physical source cuts exceed the one-hole policy')
        from .one_sided_paired_junctions import certify
        certify(group)
        inventory = {source.region.owner: source for source in group.footprints}
        selected = {}
        for li, loop in enumerate(group.loops):
            if loop.hole:
                continue
            for at, span in enumerate(loop.spans):
                label = rim_label(span, inventory)
                if label.kind != 'mouth':
                    continue
                layout = inventory[span.owner].layout
                approach = layout.approaches[label.index]
                if not hasattr(approach, 'near_role'):
                    continue
                peak = .5*(approach.profile.peak+1.)
                if abs(label.first-peak) <= 1e-10:
                    selected[approach.end_index] = BoundaryRef(li, at, 0.)
        if set(selected) != {0, 1}:
            raise ValueError('Paired physical source cuts need both existing complete mouth peaks')
        # end_index=1 is the preceding visit's tail, end_index=0 the following
        # visit's head, in the original Bundle member traversal direction.
        first, last = selected[1], selected[0]
        owner = group.members[0]
        if not holes:
            cuts = [JunctionCut(None, owner, owner, first, last)]
        else:
            li, hole = holes[0]
            a = group.loops[first.loop].xy[first.edge]
            b = group.loops[last.loop].xy[last.edge]
            direction = b-a
            projected = hole.xy@direction
            low, high = int(np.argmin(projected)), int(np.argmax(projected))
            if low == high:
                raise ValueError('Paired physical source hole has no resolved original tips')
            cuts = [JunctionCut(None, owner, owner, first, BoundaryRef(li, low, 0.)),
                    JunctionCut(None, owner, owner, last, BoundaryRef(li, high, 0.))]
        # One explicitly owned lateral source chain supplies the transverse
        # bands. Opposed nearly parallel chains are never superimposed.
        from .model import Key
        near_ids = {approach.source.id for approach in group.physical_layout.approaches
                    if hasattr(approach, 'near_role')}
        if self.source_intervals and (len(self.source_intervals) != len(near_ids)
                or {interval.section_id for interval in self.source_intervals} != near_ids):
            raise ValueError('Paired physical cuts require exactly the owned original near-source intervals')
        for interval in self.source_intervals:
            path = context[Key('path', interval.section_id)]
            if not 0. <= interval.first < interval.last <= path.length:
                raise ValueError('Paired physical source cuts cannot extrapolate an original source interval')
            count = max(2, int(np.ceil((interval.last-interval.first)/(2.*path.section.width))))
            if count > 64:
                raise ValueError('Paired physical source cuts exceed the finite band budget')
            for station in np.linspace(interval.first, interval.last, count+1)[1:-1]:
                point, tangent, normal = path.frame(station)
                hits = []
                for li, boundary in enumerate(group.loops):
                    along = (boundary.xy-point)@tangent
                    following = np.roll(along, -1)
                    crossing = ((along <= 0.) & (0. < following)) | ((following <= 0.) & (0. < along))
                    for at in np.flatnonzero(crossing):
                        fraction = float(-along[at]/(following[at]-along[at]))
                        xy = boundary.xy[at]+fraction*(boundary.xy[(at+1) % len(boundary.xy)]-boundary.xy[at])
                        hits.append((float((xy-point)@normal), li, int(at), fraction))
                hits.sort()
                if len(hits) % 2:
                    raise ValueError('Paired source normal has unresolved physical incidence')
                for a, b in zip(hits[::2], hits[1::2]):
                    if any(rim_label(group.loops[hit[1]].spans[hit[2]], inventory).kind != 'curb' for hit in (a, b)):
                        continue
                    def reference(hit):
                        _, li, edge, fraction = hit
                        if self.existing_vertices:
                            return BoundaryRef(li, (edge+int(fraction > .5)) % len(group.loops[li].xy), 0.)
                        return BoundaryRef(li, edge, fraction)
                    cuts.append(JunctionCut(None, owner, owner, reference(a), reference(b)))
                if len(cuts) > 128:
                    raise ValueError('Paired physical source cuts exceed their finite cut budget')
        return tuple(cuts)
