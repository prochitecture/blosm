"""Retain sources absorbed by logical Bundle aggregation for physical providers.

The ordinary logical inventory remains unchanged. A physical provider requests
an explicit merged view and must still prove surface ownership and coverage.
Neither retained source availability nor this view claims generated geometry.
"""
from types import SimpleNamespace


class PhysicalStreetInventory:
    @classmethod
    def from_report(cls, report):
        """Restore the complete physical view, including absorbed Bundle roles."""
        def index(key):
            rows=report.get(key,())
            result={row['id']:row for row in rows}
            if len(result)!=len(rows):
                raise ValueError('Duplicate original IDs in '+key)
            return result
        return cls.view(SimpleNamespace(streets=index('streets'),bundles=index('bundles'),
            retained_physical_streets=index('retained_physical_streets'),
            retained_physical_bundles=index('retained_physical_bundles'),
            junction_origins=report['junction_origins']))

    def review(self, snapshot):
        owners = {}
        for site in snapshot.junction_origins.get('physical', ()):
            for connector in site['ring']:
                for endpoint in connector['endpoints']:
                    owners.setdefault((endpoint['street_id'], endpoint['end']), []).append(
                        (site['junction_id'], endpoint))
        sections = {}
        all_sources = list(snapshot.streets.values()) + list(snapshot.retained_physical_streets.values())
        for street in all_sources:
            for item in street['items']:
                if item['kind'] == 'Section':
                    sections.setdefault(item['id'], []).append(street['id'])
        for sid, row in snapshot.retained_physical_streets.items():
            row['source_role'] = 'RETAINED_PHYSICAL'
            row['coverage'] = 'PENDING_PHYSICAL_PROVIDER'
            reason = row.get('reason')
            ends = [owners.get((sid, end), ()) for end in (0, 1)]
            original = [item for item in row['items'] if item['kind'] == 'Section']
            if any(len(matches) != 1 for matches in ends):
                reason = 'Retained Street needs two unique original physical endpoint owners'
            elif not original or any(ends[end][0][1]['terminal_kind'] != 'Section' or
                    ends[end][0][1]['terminal_id'] != original[0 if end == 0 else -1]['id']
                    for end in (0, 1)):
                reason = 'Retained Street no longer matches its original terminal Sections'
            elif any(len(sections[item['id']]) != 1 for item in original):
                reason = 'Retained physical Street Sections already have another source owner'
            else:
                # Obsolete linked connectors are not the ownership authority.
                row['endpoints'] = tuple(matches[0][0] for matches in ends)
                row['endpoint_space'] = 'PHYSICAL_JUNCTION'
            if reason:
                row.update(status='UNSUPPORTED', reason=reason)
            snapshot.issue('Street', sid, 'RETAINED_PHYSICAL_SOURCE', reason or
                'Original source absorbed by logical aggregation; physical surface coverage remains pending')

    @staticmethod
    def view(snapshot):
        """Read active and retained original Streets through an explicit view."""
        retained = getattr(snapshot, 'retained_physical_streets', {})
        if set(snapshot.streets) & set(retained):
            raise ValueError('A Street cannot have active and retained physical ownership')
        retained_bundles = getattr(snapshot, 'retained_physical_bundles', {})
        if set(snapshot.bundles) & set(retained_bundles):
            raise ValueError('A Bundle cannot have active and retained physical ownership')
        return SimpleNamespace(streets={**snapshot.streets, **retained},
                               bundles={**snapshot.bundles, **retained_bundles},
                               junction_origins=snapshot.junction_origins)


class PhysicalBundleInventory:
    """Keep original Bundle roles consumed by logical intersection aggregation.

    Input objects are still referenced by their original retained Streets.
    Their head/tail ordering and recorded locations are authoritative; no
    spatial pairing or replacement IDs are inferred. This exports ownership,
    not an extra rendered Bundle beside a future compound surface.
    """
    def read(self, snapshot, bundles, topology_reader):
        sources = PhysicalStreetInventory.view(snapshot).streets
        target = snapshot.retained_physical_bundles
        for bundle in bundles:
            bid = bundle.id
            if bid in snapshot.bundles or bid in target:
                raise ValueError('An original Bundle has duplicate inventory ownership')
            row = dict(id=bid, head=tuple(s.id for s in bundle.streetsHead),
                tail=tuple(s.id for s in bundle.streetsTail),
                head_locations=tuple(tuple(float(v) for v in p[:2]) for p in bundle.headLocs),
                tail_locations=tuple(tuple(float(v) for v in p[:2]) for p in bundle.tailLocs),
                owned_streets=tuple(sorted(sid for sid,s in sources.items() if s.get('bundle_id')==bid)),
                source_role='RETAINED_PHYSICAL',coverage='PENDING_PHYSICAL_PROVIDER')
            try:
                topology = topology_reader.read(row, sources, origins=snapshot.junction_origins).report()
                row['topology'] = topology
                members = topology['members']
                if (len(members)!=2 or any(not m['complete'] for m in members)
                        or any(m['head_order']!=m['tail_order'] for m in members)):
                    raise ValueError('Retained Bundle needs two complete members with consistent lateral order')
                row.update(status='SUPPORTED',reason=None,members=[dict(street_id=m['street_id'],
                    forward=m['forward'],head_order=m['head_order'],tail_order=m['tail_order']) for m in members])
            except ValueError as exc:
                row.update(status='UNSUPPORTED',reason=str(exc))
                row.setdefault('topology',dict(status='UNRESOLVED',reason=str(exc)))
            target[bid] = row
            snapshot.issue('Bundle',bid,'RETAINED_PHYSICAL_BUNDLE',row.get('reason') or
                'Original Bundle absorbed by logical aggregation; physical surface ownership remains pending')
