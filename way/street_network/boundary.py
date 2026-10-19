"""Explicit gaps around unsupported junctions keep neighboring roads inspectable."""
import math
from ..way_properties import wayCategoryProps


class UnresolvedBoundaryPolicy:
    def apply(self, snapshot, candidates, active_junctions):
        radii = {}
        def extent(section):
            return (wayCategoryProps.get(section['category'], {}).get('radius', 2.)
                + float(section['width'] or 0.)*.5 + abs(float(section['offset'] or 0.)))
        for street in snapshot.streets.values():
            sections = [i for i in street['items'] if i['kind']=='Section']
            if not sections: continue
            for end, jid in enumerate(street['endpoints']):
                radii[jid] = max(radii.get(jid, 0.), extent(sections[-1 if end else 0]))
            for item in street['items']:
                if item.get('surface_treatment') == 'NATIVE_JUNCTION':
                    radii[item['id']] = max(radii.get(item['id'],0.), *(extent(s) for s in sections))
        for sid in candidates:
            street = snapshot.streets[sid]
            sections = [i for i in street['items'] if i['kind']=='Section']
            cuts = []
            for end, jid in enumerate(street['endpoints']):
                row = snapshot.intersections.get(jid)
                if jid is None or jid in active_junctions or row and row['order'] < 3:
                    continue
                cuts.append((sections[-1 if end else 0], end, jid))
            by_id = {s['id']:s for s in sections}
            for item in street['items']:
                if item.get('surface_treatment') == 'NATIVE_JUNCTION' and item['id'] not in active_junctions:
                    cuts.extend((by_id[section], end, item['id']) for section,end in item['through'])
            for section, end, jid in cuts:
                # A visible exclusion around the unknown mouth, not an attempt
                # to replace its missing intersection with independent strips.
                row = snapshot.intersections.get(jid, {})
                radius = radii.get(jid, extent(section))
                point = section['points'][-1 if end else 0][:2]
                center = row.get('location', point)
                clearance = radius + math.dist(point, center[:2])
                key = 'clip_end' if end else 'clip_start'
                section[key] = max(section.get(key, 0.), clearance)
                street.setdefault('unresolved_boundaries', []).append(dict(intersection_id=jid,
                    section_id=section['id'], end=end, clearance_m=clearance,
                    reason=row.get('reason', 'A required mouth is outside the supported component domain')))
