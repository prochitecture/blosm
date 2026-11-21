"""Original bank/station ancestry through finite terminal polygon clipping."""
from dataclasses import dataclass, replace
import math
import numpy as np
from .model import Key
from .regions import signed_area


@dataclass(frozen=True)
class OriginalTerminalBank:
    role: int
    section: int
    side: int
    first: float
    last: float

    def __post_init__(self):
        Key('section', self.section)
        if (type(self.role) is not int or self.role not in (0, 1) or self.side not in (-1, 1)
                or not math.isfinite(self.first+self.last) or min(self.first, self.last) < 0.
                or self.first == self.last):
            raise ValueError('Terminal curb ancestry needs one finite original role/side/station interval')

    def subspan(self, first, last):
        return replace(self, first=self.first+first*(self.last-self.first),
                       last=self.first+last*(self.last-self.first))


@dataclass(frozen=True)
class ClipEdge:
    kind: str
    bank: object = None

    def __post_init__(self):
        if self.kind not in ('curb', 'internal', 'mouth') or (self.kind == 'curb') != isinstance(self.bank, OriginalTerminalBank):
            raise ValueError('Terminal edge is missing or has ambiguous original curb ancestry')

    def subspan(self, first, last):
        return replace(self, bank=self.bank.subspan(first, last) if self.bank is not None else None)


def original_strip(approach, role):
    """Use the declared ray/station arrays; no nearest-bank assignment."""
    path, visit = approach.paths[role], approach.visits[role]
    first, last = approach.field.curves[2*role:2*role+2]
    stations = np.asarray(approach.source_curve_stations[role])
    if (len(stations) != len(first.points) or len(stations) != len(last.points)
            or len(stations) < 2 or not np.isfinite(stations).all()
            or stations.min() < 0. or stations.max() > path.length
            or np.any((1. if visit.forward else -1.)*np.diff(stations) <= 0.)):
        raise ValueError('Terminal curb ancestry lost its declared original source station arrays')
    direction = 1 if visit.forward else -1
    count = len(stations)
    tags = [ClipEdge('curb', OriginalTerminalBank(role, path.section.id, -direction, float(a), float(b)))
            for a, b in zip(stations, stations[1:])]
    tags.append(ClipEdge('internal'))
    tags.extend(ClipEdge('curb', OriginalTerminalBank(role, path.section.id, direction, float(a), float(b)))
                for a, b in zip(stations[:0:-1], stations[-2::-1]))
    tags.append(ClipEdge('internal'))
    points = np.vstack((first.points, last.points[::-1]))
    if len(points) != 2*count or len(tags) != len(points):
        raise ValueError('Terminal curb ancestry changed its original bank inventory')
    return points, tuple(tags)


def validate(approach, bank, points):
    """Validate the declared lineage; this never searches for another owner."""
    if not isinstance(bank, OriginalTerminalBank):
        raise ValueError('Terminal curb is missing one unambiguous original source bank')
    path = approach.paths[bank.role]
    if path.section.id != bank.section:
        raise ValueError('Terminal curb ancestry changed its original Section/role identity')
    first, last = approach.source_ranges[bank.role]
    if not first-1e-9 <= min(bank.first, bank.last) < max(bank.first, bank.last) <= last+1e-9:
        raise ValueError('Terminal curb ancestry extrapolates its original finite source interval')
    offset = path.section.offset+.5*path.section.width*bank.side
    expected = path.point(np.asarray((bank.first, bank.last)), offset)
    error = float(np.linalg.norm(expected-np.asarray(points), axis=1).max())
    if error > approach.source_curve_tolerance+1e-8:
        raise ValueError('Terminal curb ancestry disagrees with its declared original offset side')
    return error


def clip(points, tags, center, normal):
    """The existing half-plane arithmetic, with outgoing-edge ancestry."""
    points = np.asarray(points)
    if len(tags) != len(points):
        raise ValueError('Terminal clipping needs one explicit ancestry tag per original edge')
    output, labels = [], []
    for a, b, tag in zip(points, np.roll(points, -1, axis=0), tags):
        da, db = float((a-center)@normal), float((b-center)@normal)
        crossing = (da < 0. < db) or (db < 0. < da)
        fraction = da/(da-db) if crossing else None
        if da <= 0.:
            output.append(a)
            labels.append(ClipEdge('mouth') if da == 0. and db >= 0. else
                          tag.subspan(0., fraction if crossing else 1.))
        if crossing:
            output.append(a+(b-a)*fraction)
            labels.append(ClipEdge('mouth') if da < 0. else tag.subspan(fraction, 1.))
    values = np.asarray(output, dtype=float)
    if len(values) < 3:
        raise ValueError('A terminal fibre removes an entire finite source strip')
    keep = np.r_[True, np.linalg.norm(np.diff(values, axis=0), axis=1) > 1e-8]
    selected = np.flatnonzero(keep)
    # Coalescing follows the established geometry exactly. The last outgoing
    # edge before the next retained vertex carries the surviving source tag.
    retained = [labels[(int(selected[(i+1) % len(selected)])-1) % len(labels)]
                for i in range(len(selected))]
    values = values[keep]
    if np.linalg.norm(values[0]-values[-1]) <= 1e-8:
        values = values[:-1]
        retained = retained[:-1]
    if signed_area(values) < 0.:
        values = values[::-1]
        retained = [retained[i].subspan(1., 0.) for i in np.roll(np.arange(len(retained))[::-1], -1)]
    if len(values) != len(retained):
        raise ValueError('Terminal clipping lost its finite source edge ancestry')
    return values, tuple(retained)
