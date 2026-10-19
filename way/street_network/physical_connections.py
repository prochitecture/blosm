"""Extend a physical surface selection using exact original endpoint owners.

No distance/angle matching, geometry cuts, profile averaging or fake free ends.
Every added junction brings its complete ring and original Street sources.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class PhysicalConnections:
    junctions: tuple
    added_sites: tuple
    internal_links: tuple  # Street ID, authored-head owner, authored-tail owner.
    cycle_rank: int


class PreparedPhysicalConnections:
    def __init__(self, physical_reader):
        self.reader=physical_reader

    @staticmethod
    def _incidences(junctions):
        rows={}
        for junction in junctions:
            for connector in junction.connectors:
                if len(connector.members)!=1:
                    raise ValueError('Physical connection needs individual original Street incidences')
                member,=connector.members;street=junction.streets[member.street_index]
                key=(street.source_id,0 if member.leaving else 1)
                if key in rows:raise ValueError('Original Street end has multiple physical owners')
                rows[key]=(junction,street,junction.profiles[connector.profile_index])
        return rows

    def extend(self, snapshot, junctions, street_ids):
        """Follow selected external Streets to their original far physical sites.

        The selection is explicit so qualification can grow by connected groups.
        Remaining external incidences are still pending; cycle_rank describes
        this selected connected ribbon graph, not full OSM coverage.
        """
        junctions=tuple(junctions);ids=tuple(street_ids)
        if not ids or len(set(ids))!=len(ids):raise ValueError('Connection selection must be nonempty and unique')
        old=self._incidences(junctions);existing={j.source_id for j in junctions}
        if len(existing)!=len(junctions):raise ValueError('Physical junction selection has duplicate owners')
        owners={}
        for site in snapshot.junction_origins['physical']:
            for connector in site['ring']:
                for endpoint in connector['endpoints']:
                    key=(endpoint['street_id'],endpoint['end'])
                    owners.setdefault(key,[]).append(site['junction_id'])
        added=set()
        for sid in ids:
            known=[key for key in old if key[0]==sid]
            if len(known)!=1:raise ValueError('Selected connection must have exactly one currently owned end')
            far=(sid,1-known[0][1]);matches=owners.get(far,())
            if len(matches)!=1:raise ValueError('Far Street end needs one exact recorded physical owner')
            if matches[0] in existing:raise ValueError('Far owner is present but its original incidence is missing')
            added.add(matches[0])
        extension=self.reader.prepare_sites(snapshot,tuple(sorted(added)))
        complete=junctions+extension.junctions;all_rows=self._incidences(complete);links=[]
        for sid in sorted({key[0] for key in all_rows}):
            if (sid,0) not in all_rows or (sid,1) not in all_rows:continue
            a,b=(all_rows[(sid,end)] for end in (0,1))
            if a[1]!=b[1]:raise ValueError('Original Street differs between its endpoint owners')
            pa,pb=a[2],b[2]
            if (pa.rise,pa.preset,pa.separator,pa.half_profile,a[0].datum)!=(pb.rise,pb.preset,pb.separator,-pb.half_profile,b[0].datum):
                raise ValueError('Connected Street needs explicit profile/height transition treatment')
            links.append((sid,a[0].source_id,b[0].source_id))
        neighbors={j.source_id:set() for j in complete}
        for _sid,a,b in links:neighbors[a].add(b);neighbors[b].add(a)
        reached=set();pending=[complete[0].source_id]
        while pending:
            at=pending.pop()
            if at in reached:continue
            reached.add(at);pending.extend(neighbors[at]-reached)
        if len(reached)!=len(complete):raise ValueError('Physical extension must form one connected source graph')
        return PhysicalConnections(complete,tuple(sorted(added)),tuple(links),len(links)-len(complete)+1)
