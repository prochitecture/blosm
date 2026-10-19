"""Retain original OSM incidence outside the centerline clipping rectangle.

This is sparse preparation data, not additional network geometry. It records
only clipped source ends and the original segments incident to their outside
nodes. Neither proximity nor matching XY positions create a connection.
"""
from collections import defaultdict
from copy import deepcopy
import math


def _point(value):
    point = tuple(float(v) for v in value[:2])
    if len(point) != 2 or not all(math.isfinite(v) for v in point):
        raise ValueError('Clipped source context requires finite XY coordinates')
    return point


def _identity(value):
    return None if value is None or isinstance(value, bool) else str(value)


class ClippedSourceContext:
    """Collect cuts once, then read only the needed outside node incidences."""
    def __init__(self, bounds, *, tolerance=.001, source_selection='UNSPECIFIED'):
        self.bounds = tuple(float(v) for v in bounds)
        if (len(self.bounds) != 4 or not all(math.isfinite(v) for v in self.bounds)
                or self.bounds[0] >= self.bounds[1] or self.bounds[2] >= self.bounds[3]
                or not math.isfinite(tolerance) or tolerance < 0):
            raise ValueError('Require finite clipping bounds and nonnegative tolerance')
        self.tolerance = tolerance
        if source_selection not in ('UNSPECIFIED','ALL_WAYS','SURFACE_VEHICLES'):
            raise ValueError('Unknown clipped source preparation selection')
        self.source_selection = source_selection
        self._nodes = {}
        self._cuts = defaultdict(dict)
        self._issues = []
        self._report = None

    def _distance(self, point):
        x,y = point; low_x,high_x,low_y,high_y = self.bounds
        return math.hypot(max(low_x-x,0.,x-high_x),max(low_y-y,0.,y-high_y))

    def _on_boundary(self, point):
        x,y = point; low_x,high_x,low_y,high_y = self.bounds; t = self.tolerance
        return (low_x-t <= x <= high_x+t and low_y-t <= y <= high_y+t
                and min(abs(x-low_x),abs(x-high_x),abs(y-low_y),abs(y-high_y)) <= t)

    def observe(self, segment, source_record, clipped_points):
        """Read the exact source interval already assigned to an accepted NetSection."""
        if self._report is not None: raise ValueError('Clipped context is already finalized')
        way,first,last,start,stop = source_record
        way,first,last = _identity(way),_identity(first),_identity(last)
        if not way or not first or not last:
            self._issues.append(dict(code='MISSING_OSM_SOURCE_ID', message='Outside context cannot be resolved'))
            return
        if (not all(type(v) in (int,float) and math.isfinite(v) and -1e-6 <= v <= 1.+1e-6
                    for v in (start,stop)) or start == stop):
            raise ValueError('Invalid clipped source interval')
        if (first,last) != (_identity(segment.id1),_identity(segment.id2)):
            raise ValueError('Clipped source record differs from its original segment')
        originals = (_point(segment.v1),_point(segment.v2))
        for end,parameter,inward in ((0,start,stop),(1,stop,start)):
            point = _point(clipped_points[end])
            if not 0. < parameter < 1. or not self._on_boundary(point): continue
            source_end = int(parameter > inward)
            node = (first,last)[source_end]; outside = originals[source_end]
            if self._distance(outside) <= self.tolerance: continue
            if node in self._nodes and self._nodes[node] != outside:
                raise ValueError('One OSM node has conflicting projected coordinates')
            self._nodes[node] = outside
            key = way,first,last,parameter
            row = dict(source_way_id=way, source_nodes=(first,last), source_interval=(start,stop),
                       parameter=parameter, source_end=source_end, clipped_point=point)
            existing = self._cuts[node].get(key)
            if existing and existing['clipped_point'] != point:
                raise ValueError('One clipped source parameter has conflicting positions')
            self._cuts[node][key] = row

    def finalize(self, ways):
        """Read adjacent original segments, including legs wholly outside the rectangle."""
        if self._report is not None: raise ValueError('Clipped context is already finalized')
        legs = defaultdict(dict)
        # The caller supplies the same preparation selection. No whole Way,
        # manager or parser object is retained by the resulting context.
        for way in ways:
            way_id = _identity(getattr(way,'source_id',None))
            if not way_id: continue
            for segment in way.segments:
                ids = (_identity(segment.id1),_identity(segment.id2))
                for end,node in enumerate(ids):
                    if node not in self._nodes: continue
                    position = _point((segment.v1,segment.v2)[end])
                    if position != self._nodes[node]:
                        raise ValueError('Outside OSM incidence has conflicting coordinates')
                    other = ids[1-end]
                    if not other: raise ValueError('Outside OSM incidence has no opposite node')
                    row = dict(source_way_id=way_id, other_node_id=other,
                        other_position=_point((segment.v1,segment.v2)[1-end]), category=way.category,
                        tags=deepcopy(way.element.tags), source_nodes=ids, source_end=end)
                    key = way_id,other
                    if key in legs[node] and legs[node][key] != row:
                        raise ValueError('Repeated outside OSM incidence conflicts')
                    legs[node][key] = row
        rows = []
        for node,position in sorted(self._nodes.items()):
            incidences = [self._cuts[node][key] for key in sorted(self._cuts[node])]
            adjacent = [legs[node][key] for key in sorted(legs[node])]
            degree = len({leg['other_node_id'] for leg in adjacent})
            rows.append(dict(osm_node_id=node, position=position, boundary_distance_m=self._distance(position),
                topology_kind='JUNCTION' if degree>2 else 'CONTINUATION' if degree==2 else 'ENDPOINT',
                incidences=incidences, legs=adjacent))
        self._report = dict(version=1, policy='ORIGINAL_OUTSIDE_SOURCE_CONTEXT_ONLY',
            source_selection=self.source_selection,bounds_m=self.bounds,
            tolerance_m=self.tolerance, nodes=rows, issues=self._issues)
        return self

    def report(self):
        if self._report is None: raise ValueError('Clipped context has not been finalized')
        return deepcopy(self._report)


def resolve_open_ends(context, streets, retained_streets=()):
    """Address final Street/Section ends by exact source incidence and parameter."""
    result = deepcopy(context)
    if not result: return result
    if result.get('version') != 1: raise ValueError('Unknown clipped source context version')
    index = {}
    for node in result['nodes']:
        node['open_ends'] = []
        for cut in node['incidences']:
            key = cut['source_way_id'],*cut['source_nodes'],cut['parameter']
            if key in index: raise ValueError('Ambiguous clipped source parameter')
            index[key] = node,cut
    for scope,inventory in (('ACTIVE',streets),('RETAINED_PHYSICAL',retained_streets)):
        for street in inventory:
            sections = [s for s in street['items'] if s['kind']=='Section']
            if not sections: continue
            for end,junction in enumerate(street['endpoints']):
                if junction is not None: continue
                section = sections[0 if end==0 else -1]
                source = section.get('source_segments',())
                if not source or not section['points']: continue
                way,first,last,start,stop = source[0 if end==0 else -1]
                key = str(way),str(first),str(last),start if end==0 else stop
                if key not in index: continue
                node,cut = index[key]
                point = _point(section['points'][0 if end==0 else -1])
                if math.dist(point,cut['clipped_point']) > result['tolerance_m']:
                    result['issues'].append(dict(code='CLIPPED_SOURCE_POSITION_CHANGED',street_id=street['id'],
                        section_id=section['id'],end=end,osm_node_id=node['osm_node_id']))
                    continue
                node['open_ends'].append(dict(street_id=street['id'],section_id=section['id'],end=end,
                    owner_scope=scope,source_way_id=cut['source_way_id'],source_nodes=cut['source_nodes'],
                    parameter=cut['parameter']))
    return result
