"""Topological contraction of two original junction rings across one span.

The plan carries only original Section ends and ordered incidences. Geometry,
corner overlap and transition feasibility remain live solver responsibilities.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class PortCompound:
    site_ids: tuple
    span_identity: tuple
    internal_ports: tuple
    external_ports: tuple
    # External side ordinal, original corner at its start, original corner
    # at its end. Each original corner is addressed by its first port.
    bridges: tuple


class PreparedPortCompound:
    def prepare(self, plan, span_identity):
        spans = [s for s in plan.spans if s.identity == tuple(span_identity)]
        if len(spans) != 1:
            raise ValueError('A compound needs one existing original span')
        span = spans[0]
        owners = (span.start.junction_id, span.end.junction_id)
        if None in owners or owners[0] == owners[1]:
            raise ValueError('The internal span needs two distinct junction owners')
        sites = []
        for jid in owners:
            matches = [s for s in plan.sites if s.junction_id == jid]
            if len(matches) != 1:
                raise ValueError('Both original junction rings must be present')
            sites.append(matches[0])
        internal, halves = [], []
        for site, end in zip(sites, (span.start, span.end)):
            indices = [i for i, p in enumerate(site.ports) if p.endpoint == end.endpoint]
            if len(indices) != 1 or len(site.ports) < 3:
                raise ValueError('Each ring needs exactly one end of the internal span')
            at = indices[0]; internal.append(site.ports[at])
            halves.append(tuple(site.ports[(at+i) % len(site.ports)] for i in range(1, len(site.ports))))
        external = halves[0]+halves[1]
        if len({p.endpoint for p in external}) != len(external):
            raise ValueError('A compound repeats an external original Section end')
        # A second connection forms a cycle/hole and is a different topology;
        # never absorb it silently or reject two unrelated spans of one Street.
        identities = [plan.span_for(p).identity for p in external]
        if len(set(identities)) != len(identities):
            raise ValueError('Two-ring contraction cannot absorb a second shared span')
        bridges = ((len(halves[0])-1, halves[0][-1], internal[1]),
                   (len(external)-1, halves[1][-1], internal[0]))
        return PortCompound(owners, span.identity, tuple(internal), external, bridges)
