"""Physical ownership planning for overlapping junction footprints.

This preparatory stage is independent of internal crest guides. It is not yet
a compound pavement provider: holes, profile reconciliation, quad charts and
surface qualification must succeed before a plan can replace local surfaces.
"""
from dataclasses import dataclass
import numpy as np
from .regions import Region, RegionUnion


@dataclass(frozen=True)
class RimSegment:
    kind: str
    index: int
    first: float
    last: float


@dataclass(frozen=True)
class SampledFootprint:
    region: Region
    segments: tuple
    layout: object = None


@dataclass(frozen=True)
class RimSpan:
    """An analytical boundary interval after redistributing its samples."""
    owner: object
    label: RimSegment


def rim_label(span, sources):
    if isinstance(span, RimSpan):
        return span.label
    label = sources[span.owner].segments[span.segment]
    return RimSegment(label.kind, label.index,
        label.first + span.first * (label.last - label.first),
        label.first + span.last * (label.last - label.first))


def sample_footprint(owner, layout):
    """Retain the existing curb accuracy and every incoming profile vertex.

    mouth parameters are public transverse u in [0,1]; curb parameters belong
    to the corresponding BoundaryCurve. SourceSpan from RegionUnion addresses
    these sampled segments, so clipping does not erase source provenance.
    """
    rim, labels = [], []
    for side, approach in enumerate(layout.approaches):
        u = np.linspace(1., 0., 2 * layout.radial + 1)
        points = approach.point(layout.outward_stations[side], approach.offset + approach.width * (u - .5))
        rim.extend(points[:-1])
        labels.extend(RimSegment('mouth', side, float(a), float(b)) for a, b in zip(u, u[1:]))
        t = np.linspace(0., 1., layout.rows + 1)
        rim.extend(layout.boundaries[side](t)[:-1])
        labels.extend(RimSegment('curb', side, float(a), float(b)) for a, b in zip(t, t[1:]))
    return SampledFootprint(Region(owner, np.asarray(rim)), tuple(labels), layout)


@dataclass(frozen=True)
class SharedClaim:
    section_id: int
    first_owner: object
    last_owner: object
    first_trim: float
    last_trim: float


@dataclass(frozen=True)
class CompoundGroup:
    members: tuple
    claims: tuple
    footprints: tuple
    loops: tuple
    # Logical members/claims survive when a shared physical exterior replaces
    # invalid clipped internal mouths. Only derived geometry belongs here.
    physical_layout: object = None


class CompoundPlanner:
    def __init__(self, union=None):
        self.union = RegionUnion() if union is None else union

    def build(self, layouts):
        """Group actual sites by overlapping intervals on original Sections.

        Compatible two-ended seams are already accounted for by Route.trims.
        Proximity alone cannot connect unrelated or grade-separated streets.
        Singleton sites retain their existing local construction.
        """
        layouts = {key: value for key, value in layouts.items() if value.corners}
        by_section = {}
        parent = {key: key for key in layouts}

        def root(key):
            while parent[key] != key:
                parent[key] = parent[parent[key]]
                key = parent[key]
            return key

        for key, layout in layouts.items():
            for (sid, end), mouth in layout.mouths.items():
                by_section.setdefault(sid, {0: [], 1: []})[end].append((key, mouth.station))
        claims = []
        for sid, ends in sorted(by_section.items()):
            for first, lo in sorted(ends[0]):
                for last, hi in sorted(ends[1]):
                    if first == last or lo <= hi + self.union.epsilon:
                        continue
                    a, b = sorted((root(first), root(last)))
                    parent[b] = a
                    claims.append(SharedClaim(sid, first, last, float(lo), float(hi)))
        groups = {}
        for key in sorted(parent):
            groups.setdefault(root(key), []).append(key)
        result = []
        for members in groups.values():
            if len(members) < 2:
                continue
            footprints = tuple(sample_footprint(key, layouts[key]) for key in members)
            loops = self.union.build(footprint.region for footprint in footprints)
            if sum(not loop.hole for loop in loops) != 1:
                raise ValueError('Overlapping source claims do not form one connected compound footprint')
            result.append(CompoundGroup(tuple(members), tuple(c for c in claims if c.first_owner in members),
                                        footprints, loops))
        return tuple(result)
