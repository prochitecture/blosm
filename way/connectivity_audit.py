"""Inspect the legacy XY graph's connection assumptions using OSM provenance.

This records evidence; it does not infer heights from layer tags or silently
rewrite the preparation graph. Grade-aware graph construction remains separate.
"""
from collections import Counter, defaultdict
from itertools import combinations


class ConnectivityAudit:
    @staticmethod
    def grade(tags):
        return tuple(str(tags.get(key, default)).strip().lower()
            for key, default in (('layer', '0'), ('bridge', 'no'), ('tunnel', 'no')))

    def read(self, ways, splits):
        events = defaultdict(list)
        for way in ways:
            for segment in way.segments:
                first, last = tuple(segment.v1), tuple(segment.v2)
                points = splits.get((first, last), splits.get((last, first), (first, last)))
                for point in points:
                    point = tuple(point)
                    node_ids = set()
                    if point == first: node_ids.add(str(segment.id1))
                    if point == last: node_ids.add(str(segment.id2))
                    events[point].append((way, node_ids))
        rows = []
        seen = set()
        for point, participants in events.items():
            for (a, an), (b, bn) in combinations(participants, 2):
                if an & bn:
                    continue  # Explicit OSM connectivity, including grade transitions.
                ids = tuple(sorted((str(a.source_id), str(b.source_id))))
                key = point, ids
                if key in seen: continue
                seen.add(key)
                grade_a, grade_b = self.grade(a.element.tags), self.grade(b.element.tags)
                rows.append(dict(location=point, osm_way_ids=ids,
                    code='GRADE_SEPARATED_XY' if grade_a != grade_b else 'INFERRED_XY_CONNECTION',
                    categories=(a.category, b.category), grades=(grade_a, grade_b)))
        return dict(policy='Audit only: legacy graph uses XY intersections; these events lack a shared OSM endpoint node',
                    counts=dict(Counter(row['code'] for row in rows)), events=rows)
