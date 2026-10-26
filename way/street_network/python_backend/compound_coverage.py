"""Original source intervals wholly owned inside a physical compound.

An overlapping claim proves grouping, not complete coverage. A short prefix
may be wholly inside one member while only its continuation is shared. Keep
these two facts separate so authored Section splits do not create entrances.
"""
from .model import Key
from .routes import RouteTopology


class CompoundSourceCoverage:
    epsilon = 1e-7

    def __init__(self, context, group, layouts):
        self.context, self.group = context, group
        intervals = {}
        for key in group.members:
            for (sid, end), mouth in layouts[key].mouths.items():
                length = context[Key('path', sid)].length
                span = (0., mouth.station) if end == 0 else (mouth.station, length)
                intervals.setdefault(sid, []).append(span)
        complete = set()
        for sid, spans in intervals.items():
            covered = 0.
            for first, last in sorted(spans):
                if first > covered + self.epsilon:
                    break
                covered = max(covered, last)
            if covered >= context[Key('path', sid)].length - self.epsilon:
                complete.add(sid)
        self.fully_owned = frozenset(complete)

    def internal_sections(self):
        """Return complete routes ending at another actual member junction.

        A free terminal retains its exposed profile even if the source span is
        short. Conversely, connectivity alone cannot consume a long route that
        leaves and later re-enters the compound: every interval must be owned.
        Geometric containment is independently checked before replacing a rim.
        """
        topology = RouteTopology(self.context.network)
        sites = {site.key: site for site in self.context.network.intersections}
        result = set()
        for key in self.group.members:
            for end in sites[key].ends:
                route = topology.trace(end)
                ids = {visit.section_id for visit in route.traversals}
                if (route.terminal_site in self.group.members and route.terminal_site != key
                        and ids <= self.fully_owned):
                    result.update(ids)
        return frozenset(result)
