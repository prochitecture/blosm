"""Visible pavement ports indexed by their original Section endpoint identity.

A mouth can move through compatible Section seams while its original physical
junction stays put. Record that derived incidence once, without a spatial
nearest-junction search or a new source endpoint.
"""
from dataclasses import dataclass, field
import math
from .model import End, Key

PORTS = Key('separator_ports', 0)
PHYSICAL_PORT_BINDING_POLICY = '1-published-source-mouth-identities'


@dataclass(frozen=True)
class VisiblePort:
    site: object
    approach_index: int


@dataclass(frozen=True)
class DerivedVisiblePort:
    owner: object
    approach_index: int
    alias: object


@dataclass(frozen=True)
class PhysicalPortBinding:
    owner: object
    product: object
    layout: object
    approach_index: int
    alias: object
    aliases: tuple


def bind_physical_port(context, visible, endpoint):
    """Bind a real published mouth by source identity in its current layout.

    The visible index addresses the projected source layout. A completed
    compound can retain a different, source-reconciled approach inventory;
    that inventory's positional indices are not persistent mouth identities.
    """
    from .ownership import PLAN, CompoundProduct, surface_key
    from .source_mouths import source_mouths
    if not isinstance(endpoint, End) or not isinstance(visible, (VisiblePort, DerivedVisiblePort)):
        raise ValueError('Physical mouth binding needs an explicit visible source incidence')
    derived = isinstance(visible, DerivedVisiblePort)
    plan = context.products.get(PLAN)
    owner = (visible.owner if derived else
             (plan.owner(visible.site.key) if plan is not None else None) or visible.site.key)
    product = context[surface_key(owner)]
    layout = (product.group.physical_layout if isinstance(product, CompoundProduct)
              and product.group.physical_layout is not None else None)
    if layout is None:
        if derived:
            raise ValueError('A derived physical mouth needs its accepted physical owner layout')
        layout = context[visible.site.layout_key]
    expected = None
    if isinstance(product, CompoundProduct):
        identity = ((owner, endpoint.section_id, endpoint.end) if derived else
                    (endpoint.section_id, endpoint.end))
        peaks = (product.derived_peaks if derived else product.peaks) or {}
        if identity not in peaks:
            raise ValueError('Physical mouth is absent from its published source profiles')
        expected = peaks[identity][2]
    matches = []
    for index, (approach, station) in enumerate(zip(layout.approaches, layout.outward_stations)):
        aliases = source_mouths(approach, station)
        for alias in aliases:
            if alias.original_incidence != derived and alias.endpoint == endpoint:
                matches.append(PhysicalPortBinding(owner, product, layout, index, alias, aliases))
    if len(matches) != 1:
        raise ValueError('Physical mouth needs one unique published source identity')
    binding = matches[0]
    if (not math.isfinite(binding.alias.station)
            or expected is not None and (not math.isfinite(expected) or abs(expected-binding.alias.station) > 1e-8)
            or derived and binding.alias != visible.alias):
        raise ValueError('Physical mouth changed its published original source station or alias')
    return binding


@dataclass
class SeparatorPortIndex:
    ports: dict
    signatures: dict
    derived_ports: dict = field(default_factory=dict)
    derived_signatures: dict = field(default_factory=dict)

    @classmethod
    def build(cls, context):
        from .ownership import PLAN, CompoundProduct, surface_key
        ports, signatures, derived_ports, derived_signatures = {}, {}, {}, {}
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
                    if not alias.original_incidence:
                        token = owner, endpoint.section_id, endpoint.end
                        if (isinstance(product, CompoundProduct)
                                and token not in (getattr(product, 'derived_peaks', None) or {})):
                            continue
                        derived_ports.setdefault(token, []).append(DerivedVisiblePort(owner, index, alias))
                        derived_signatures.setdefault(token, []).append(stamp)
                        continue
                    if isinstance(product, CompoundProduct) and (endpoint.section_id, endpoint.end) not in product.peaks:
                        continue  # This original mouth is internal to the union.
                    # Multiple logical members share one physical layout. Each
                    # alias retains its original incidence, not the first site.
                    logical = source_sites.get(endpoint, site) if merged else site
                    ports.setdefault(endpoint, []).append(VisiblePort(logical, index))
                    signatures.setdefault(endpoint, []).append(stamp)
        return cls(ports, signatures, derived_ports, derived_signatures)

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

    def resolve_derived(self, owner, section_id, end):
        """Resolve an explicitly owned physical cut, never a source End."""
        candidates = self.derived_ports.get((owner, section_id, end), ())
        if len(candidates) > 1:
            raise ValueError('A derived physical source cut has ambiguous pavement port ownership')
        return candidates[0] if candidates else None
