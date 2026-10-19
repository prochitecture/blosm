"""Keep unverified XY graph connections out of generated carriageways.

The legacy preparation graph uses positions as node identities. An XY crossing
without a shared OSM node is not enough evidence for an at-grade connection.
Until that graph can represent separate coincident nodes, retain the affected
source Streets as diagnosed context rather than publishing a false junction.
This policy never converts OSM layer numbers into physical elevations.
"""
from collections import defaultdict


class SurfaceConnectivityReview:
    def apply(self, snapshot, categories):
        audit = snapshot.preparation_audit
        by_way = defaultdict(set)
        for sid, street in snapshot.streets.items():
            for item in street['items']:
                if item['kind'] == 'Section':
                    for wid in item.get('source_way_ids', ()):
                        by_way[str(wid)].add(sid)
        events = []
        for source in audit.get('events', ()):
            event = dict(source)
            affected = sorted({sid for wid in source['osm_way_ids'] for sid in by_way[str(wid)]
                if any(i.get('category') in categories for i in snapshot.streets[sid]['items'])})
            event['surface_street_ids'] = affected
            in_scope = any(c in categories for c in source.get('categories', ()))
            event['surface_treatment'] = ('UNVERIFIED_CONNECTIVITY' if affected else
                'UNMAPPED_SOURCE' if in_scope else 'DEFERRED_COMPONENT')
            if not affected and in_scope:
                snapshot.issue('ConnectivityEvent', len(events), 'UNMAPPED_SOURCE',
                    'Audited carriageway connection has no retained Street provenance: '
                    + ', '.join(source['osm_way_ids']))
            for sid in affected:
                street = snapshot.streets[sid]
                reason = ('Prepared XY connection has no shared OSM node; affected source ways: '
                          + ', '.join(source['osm_way_ids']))
                snapshot.issue('Street', sid, 'UNVERIFIED_CONNECTIVITY', reason)
                # Preserve any earlier scope reason and still record this issue.
                if street['status'] == 'SUPPORTED':
                    street.update(status='UNSUPPORTED', reason=reason)
            events.append(event)
        snapshot.preparation_audit = dict(audit, events=events,
            surface_policy='Require shared OSM nodes for graph connections; retain affected source Streets as unsupported context',
            surface_review_counts={kind.lower(): sum(e['surface_treatment'] == code for e in events)
                for kind, code in (('unverified', 'UNVERIFIED_CONNECTIVITY'),
                                   ('deferred', 'DEFERRED_COMPONENT'), ('unmapped', 'UNMAPPED_SOURCE'))})
