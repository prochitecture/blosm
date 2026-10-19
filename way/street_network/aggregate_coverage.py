"""Semantic surface accounting for an explicitly prepared aggregate network."""
from dataclasses import dataclass
import json


@dataclass(frozen=True)
class AggregateSurfaceCoverage:
    owners: tuple

    @classmethod
    def prepare(cls, plan, tables):
        rows = [dict(kind=0, owner=i, junctions=[site.junction_id], streets=[], bundles=[])
                for i, site in enumerate(plan.ports.sites)]
        for kind, table, prefix in ((1, 'Corridors', 'st_cs_'), (3, 'External Approaches', 'st_ca_')):
            for i, sid in enumerate(tables[table].attributes[prefix+'street_id'].data):
                rows.append(dict(kind=kind, owner=i, streets=[sid.value], junctions=[], bundles=[]))
        for bid in sorted({d.value for d in tables['Separators'].attributes['st_sb_bundle_id'].data}):
            rows.append(dict(kind=4, owner=bid, streets=[], junctions=[], bundles=[bid]))
        return cls(tuple(rows))

    @classmethod
    def from_object(cls, obj):
        return cls(tuple(json.loads(obj['blosm_aggregate_surface_owners'])))

    def encode(self):
        return json.dumps(self.owners)

    def evaluate(self, accepted):
        actual = {tuple(int(v) for v in key) for key in accepted}
        expected = {(row['kind'], row['owner']) for row in self.owners}
        missing = expected-actual; unexpected = actual-expected
        result = dict(complete=not missing and not unexpected,
            expected_owners=sorted(expected), accepted_owners=sorted(actual),
            missing_owners=sorted(missing), unexpected_owners=sorted(unexpected))
        for key in ('streets', 'junctions', 'bundles'):
            result['missing_'+key] = sorted({v for row in self.owners
                if (row['kind'], row['owner']) in missing for v in row[key]})
        result['failed_streets'] = result['missing_streets']
        result['affected_streets'] = result['missing_streets']
        return result
