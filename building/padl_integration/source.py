"""Adapt Blosm's parsed OSM rings without using legacy building generators."""
from dataclasses import dataclass, field, replace

from ...padl.errors import PadlError
from ...padl.geometry import Footprint, contains, on_segment, cross
from ...padl.bindings import FaceType
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
                    candidates.append((key, rings[id(outer)], tuple(assigned.get(id(outer), ())), relation.tags, relation, source))
        for source, way in sorted(ways.items()):
            if source.rsplit(':', 1)[-1] not in consumed:
                candidates.append((source, clean_ring(way.getData(self.osm)), (), way.tags, way, source))
        footprints = BuildingFootprintAssembler(self.mapper, self.report).assemble(candidates)
        self.report.input_count = len(footprints)
        return footprints


def tag_enabled(tags, name):
    return tags.get(name) not in (None, '', 'no', 'demolished', 'razed', 'construction', 'proposed')


def enclosed_by(part, building):
    """Contain the complete polygon, allowing shared outline edges/corners."""
    def boundary(point, ring):
        return any(on_segment(a, b, point) for a, b in zip(ring, ring[1:] + ring[:1]))

    def inside(point, ring):
        return contains(point, ring) or boundary(point, ring)
    if not part.outer or not building.outer:
        return False
    edges = tuple(zip(part.outer, part.outer[1:] + part.outer[:1]))
    points = (*part.outer, *(tuple((a[i] + b[i]) / 2 for i in (0, 1)) for a, b in edges))
    if not all(inside(p, building.outer) and not any(contains(p, h) and not boundary(p, h)
                                                   for h in building.holes) for p in points):
        return False
    return not any(cross(a, b, c) * cross(a, b, d) < -1e-12 and
                   cross(c, d, a) * cross(c, d, b) < -1e-12
                   for a, b in edges for ring in (building.outer, *building.holes)
                   for c, d in zip(ring, ring[1:] + ring[:1]))


class BuildingFootprintAssembler:
    """Resolve source ownership before mapping or selecting volume footprints."""
    def __init__(self, mapper, report):
        self.mapper, self.report = mapper, report

    def assemble(self, candidates):
        records, rejected_geometry = {}, set()
        for key, outer, holes, tags, element, source in sorted(candidates):
            footprint = Footprint(key, len(records) + 1, outer, holes)
            try:
                footprint.validate()
            except PadlError as exc:
                self.report.issue(key, exc.code, str(exc))
                rejected_geometry.add(key)
            records[key] = (footprint, tags, element, source)
        buildings = {key: row for key, row in records.items() if tag_enabled(row[1], 'building')}
        parents = {}
        for key, (part, tags, element, _) in records.items():
            if key in buildings or not tag_enabled(tags, 'building:part'):
                continue
            explicit = getattr(element, 'o', None)
            choices = list(buildings)
            if explicit:
                source = ('osm:way:' if explicit[1] == Osm.way else 'osm:relation:') + str(explicit[0])
                choices = [k for k, row in buildings.items() if row[3] == source]
            if not explicit or len(choices) != 1:
                choices = [k for k in choices if enclosed_by(part, buildings[k][0])]
            if len(choices) != 1:
                self.report.issue(key, 'part_ownership', 'Building part requires exactly one building outline; explicit relation membership takes precedence')
                continue
            parents[key] = choices[0]
        children = set(parents.values())
        result = {}
        for key, (footprint, tags, _, _) in records.items():
            if key in rejected_geometry:
                continue
            if key not in buildings and key not in parents:
                continue
            role = (FaceType.BUILDING_ALSO_PART if tag_enabled(tags, 'building:part') else
                    FaceType.BUILDING_WITH_PARTS if key in children else FaceType.BUILDING) if key in buildings else FaceType.PART
            mapped_tags = dict(tags)
            if role == FaceType.PART:
                parent_tags = buildings[parents[key]][1]
                # Inherit classification, never an overall building height/count.
                if 'building:use' not in tags and not tag_enabled(tags, 'building'):
                    form = tags.get('building:part', 'yes')
                    mapped_tags['building'] = parent_tags.get('building', 'yes') if form == 'yes' else form
                    if form == 'yes' and parent_tags.get('building:use'):
                        mapped_tags['building:use'] = parent_tags['building:use']
            mapped = self.mapper.map(mapped_tags, footprint=footprint)
            if role == FaceType.PART and tags.get('building:part') == 'yes' and not any(
                    tags.get(name) for name in ('building:use', 'office', 'amenity')):
                mapped.attributes['building_use'] = self.mapper.map(parent_tags).attributes['building_use']
            for issue in mapped.issues:
                self.report.issue(key, issue.code, issue.message)
            if not mapped.accepted and role != FaceType.BUILDING_WITH_PARTS:
                continue
            result[key] = replace(footprint, attributes=mapped.attributes, face_type=int(role),
                parent_id=buildings[parents[key]][0].footprint_id if key in parents else 0)
        # Never substitute an outline volume for rejected child parts.
        valid_buildings = {k for k, f in result.items() if FaceType.is_building(f.face_type)}
        result = {k: f for k, f in result.items() if k not in parents or parents[k] in valid_buildings}
        valid_parents = {f.parent_id for f in result.values() if f.face_type == FaceType.PART}
        return [f for f in result.values() if f.face_type != FaceType.BUILDING_WITH_PARTS or f.footprint_id in valid_parents]


def group_buildings(footprints, single_object):
    """Packing is independent of part ownership and backend implementation."""
    footprints = tuple(footprints)
    if single_object:
        return (footprints,) if footprints else ()
    groups = {f.footprint_id: [] for f in footprints if FaceType.is_building(f.face_type)}
    for footprint in footprints:
        key = footprint.parent_id or footprint.footprint_id
        if key in groups:
            groups[key].append(footprint)
    return tuple(tuple(group) for group in groups.values())
