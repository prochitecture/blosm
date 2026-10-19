"""Whole, unattached Bundle ownership over the original prepared Streets.

This adapter proves topology only. Native Geometry Nodes still checks curves,
contact, pavement unions and separator islands from the live source controls.
"""
from dataclasses import dataclass, replace
from .bundle_topology import PreparedBundleTopology
from .physical_dimensions import PreparedPhysicalDimensions
from .port_network import PreparedPortNetwork
from .port_compounds import PreparedPortCompound


@dataclass(frozen=True)
class FreeBundleSurface:
    bundle_id: int
    members: tuple

    @property
    def streets(self):
        return tuple(member.street_id for member in self.members)


class PreparedFreeBundles:
    def __init__(self, reader):
        self.reader = reader

    def prepare(self, inventory, bundle_id, *, occupied=()):
        bundle = inventory.bundles[bundle_id]
        topology = PreparedBundleTopology(self.reader.categories).read(bundle, inventory.streets)
        members = topology.members
        if len(members) != 2 or any(not m.complete for m in members):
            raise ValueError('Free Bundle surfaces require two complete original member Streets')
        if topology.junctions or topology.interior_streets or any(len(m.spans) != 1 for m in members):
            raise ValueError('Internal Bundle junctions need the attached surface bridge')
        if any(m.head_order != m.tail_order for m in members):
            raise ValueError('Bundle centerlines cannot exchange lateral order')
        streets = tuple(m.street_id for m in members)
        if set(streets).intersection(occupied):
            raise ValueError('Bundle members already belong to an original junction component')
        for sid in streets:
            street = inventory.streets[sid]
            if street.get('bundle_id') != bundle_id:
                raise ValueError('Prepared Bundle and original Street ownership disagree')
            if tuple(street.get('endpoints', ())) != (None, None):
                raise ValueError('Attached Bundle members need the attached surface bridge')
            for item in street['items']:
                if item['kind'] != 'Section' and (item['kind'] not in ('SymLane', 'SideLane', 'AuthoredTransition')
                        or item.get('surface_treatment') not in ('COMPATIBLE_SEAM', 'WIDTH_TRANSITION')
                        or item.get('points') or item.get('transition_reason')):
                    raise ValueError('Free Bundle member item has no explicit continuous sweep treatment')
        ports = self.reader.prepare(inventory, street_ids=streets)
        network = PreparedPortNetwork(PreparedPortCompound()).prepare(ports)
        if ports.sites or ports.section_scopes or len(ports.spans) != 2 or any(
                end.owner_kind not in ('FREE', 'SELECTION_BOUNDARY')
                for span in ports.spans for end in (span.start, span.end)):
            raise ValueError('Free Bundle requires complete members without junction or category incidences')
        PreparedPhysicalDimensions().read_spans(inventory, ports.spans, source_plan=ports)
        network = replace(network, routes=tuple(replace(route, kind='FREE_BUNDLE') for route in network.routes))
        return ports, network, FreeBundleSurface(bundle_id, members)
