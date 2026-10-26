"""Visible pavement ports indexed by their original Section endpoint identity.

A mouth can move through compatible Section seams while its original physical
junction stays put. Record that derived incidence once, without a spatial
nearest-junction search or a new source endpoint.
"""
from dataclasses import dataclass
from .model import End, Key

PORTS = Key('separator_ports', 0)


@dataclass(frozen=True)
class VisiblePort:
    site: object
    approach_index: int


@dataclass
class SeparatorPortIndex:
    ports: dict
    signatures: dict

    @classmethod
    def build(cls, context):
        from .ownership import PLAN, CompoundProduct, surface_key
        ports, signatures = {}, {}
        plan = context.products.get(PLAN)
        source_sites = {end: site for site in context.network.intersections for end in site.ends}
        sites_by_key = {site.key: site for site in context.network.intersections}
        visited = set()
        for site in context.network.intersections:
            layout = context[site.layout_key]
            owner = (plan.owner(site.key) if plan is not None else None) or site.key
            group = plan.by_site.get(site.key) if plan is not None else None
            if group is not None and group.physical_layout is not None:
                if owner in visited:
                    continue
                visited.add(owner)
            surface = surface_key(owner)
            product = context.products.get(surface)
            merged = any(hasattr(approach, 'source_mouths') for approach in layout.approaches)
            representative = sites_by_key[owner] if merged else site
            stamp = (representative.key.text(), context.digest(representative.layout_key),
                     context.digest(surface) if surface in context.products else None)
            for index, (approach, station) in enumerate(zip(layout.approaches, layout.outward_stations)):
                from .source_mouths import source_mouths
                aliases = source_mouths(approach, station)
                for alias in aliases:
                    endpoint = alias.endpoint
                    if isinstance(product, CompoundProduct) and (endpoint.section_id, endpoint.end) not in product.peaks:
                        continue  # This original mouth is internal to the union.
                    # Multiple logical members share one physical layout. Each
                    # alias retains its original incidence, not the first site.
                    logical = source_sites.get(endpoint, site) if merged else site
                    ports.setdefault(endpoint, []).append(VisiblePort(logical, index))
                    signatures.setdefault(endpoint, []).append(stamp)
        return cls(ports, signatures)

    def resolve(self, endpoint):
        candidates = self.ports.get(endpoint, ())
        if len(candidates) > 1:
            raise ValueError('A visible separator endpoint has ambiguous pavement port ownership')
        return candidates[0] if candidates else None

    def signature(self, section_ids):
        result = []
        for sid in sorted(set(section_ids)):
            for end in (0, 1):
                endpoint = End(sid, end)
                if endpoint in self.signatures:
                    result.append((sid, end, sorted(self.signatures[endpoint])))
        return result
