"""Adapt Blosm's parsed OSM rings without using legacy building generators."""
from dataclasses import dataclass, field

from ...padl.errors import PadlError
from ...padl.geometry import Footprint, contains
from ...parse.osm import Osm


@dataclass
class ImportReport:
    issues: list = field(default_factory=list)
    statuses: list = field(default_factory=list)
    input_count: int = 0

    def issue(self, source, code, message):
        self.issues.append({'source': source, 'code': code, 'message': message})

    def as_dict(self):
        return {'input_count': self.input_count, 'issues': self.issues, 'statuses': self.statuses}


def clean_ring(points):
    ring = [tuple(float(x) for x in point) for point in points]
    # OSM occasionally repeats a closing coordinate through distinct node IDs.
    if len(ring) > 1 and ring[0][:2] == ring[-1][:2]:
        ring.pop()
    # Retain intermediate/collinear corners and winding for future edge-based
    # façade rules. Kernel-specific simplification belongs in a roof adapter.
    return tuple(ring)


class OsmFootprintSource:
    def __init__(self, osm, mapper, report):
        self.osm, self.mapper, self.report = osm, mapper, report

    def collect(self, ways, relations):
        candidates, consumed = [], set()
        for source, relation in sorted(relations.items()):
            lines = relation.ls if isinstance(relation.ls, list) else [relation.ls]
            for line in lines:
                if line.parts:
                    consumed.update(p[0] for p in line.parts.values())
                else:
                    consumed.add(line.wayId)
            outers = [line for line in lines if line.role == Osm.outer]
            inners = [line for line in lines if line.role == Osm.inner]
            if not outers:
                self.report.issue(source, 'invalid_geometry', 'Multipolygon has no closed outer boundary')
                continue
            rings = {id(line): clean_ring(relation.getLinestringData(line, self.osm)) for line in lines}
            assigned = {}
            for inner in inners:
                ring = rings[id(inner)]
                parents = [outer for outer in outers if ring and contains(ring[0], rings[id(outer)])]
                if len(parents) != 1:
                    self.report.issue(source, 'invalid_geometry', 'Inner boundary must belong to exactly one outer boundary')
                    break
                assigned.setdefault(id(parents[0]), []).append(ring)
            else:
                for outer in outers:
                    way_ids = [p[0] for p in outer.parts.values()] if outer.parts else [outer.wayId]
                    key = source if len(outers) == 1 else source + ':outer:' + min(way_ids, key=int)
                    candidates.append((key, rings[id(outer)], tuple(assigned.get(id(outer), ())), relation.tags))
        for source, way in sorted(ways.items()):
            if source.rsplit(':', 1)[-1] not in consumed:
                candidates.append((source, clean_ring(way.getData(self.osm)), (), way.tags))
        footprints = []
        for key, outer, holes, tags in sorted(candidates):
            mapped = self.mapper.map(tags)
            for issue in mapped.issues:
                self.report.issue(key, issue.code, issue.message)
            if not mapped.accepted:
                continue
            footprint = Footprint(key, len(footprints) + 1, outer, holes, mapped.attributes)
            try:
                footprint.validate()
            except PadlError as exc:
                self.report.issue(key, exc.code, str(exc))
                continue
            footprints.append(footprint)
        self.report.input_count = len(footprints)
        return footprints
