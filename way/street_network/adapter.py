"""Read prepared Blosm topology without rewriting source centerlines.

Feature providers consume the same snapshot. Unsupported input is retained in
the inventory and inspection geometry, never coerced into a compatible fixture.
"""
from dataclasses import dataclass, field, replace
import math

from ..gn.kernel import adapter as junction_adapter, inputs
from ..way_properties import wayCategoryProps
from .design import StreetDesign
from .physical_ports import GroundSurfacePolicy


@dataclass
class NetworkSnapshot:
    streets: dict = field(default_factory=dict)
    intersections: dict = field(default_factory=dict)
    bundles: dict = field(default_factory=dict)
    junctions: dict = field(default_factory=dict)
    issues: list = field(default_factory=list)
    preparation_audit: dict = field(default_factory=dict)
    ownership_audit: dict = field(default_factory=dict)
    coupled_bundle_networks: list = field(default_factory=list)
    junction_interior_links: dict = field(default_factory=dict)
    junction_origins: dict = field(default_factory=dict)
    retained_physical_streets: dict = field(default_factory=dict)
    retained_physical_bundles: dict = field(default_factory=dict)
    clipped_source_context: dict = field(default_factory=dict)

    def issue(self, kind, source_id, code, message):
        self.issues.append(dict(kind=kind, id=source_id, code=code, message=message))

    def report(self):
        return dict(streets=list(self.streets.values()),
                    intersections=list(self.intersections.values()),
                    bundles=list(self.bundles.values()), issues=self.issues,
                    qualified_junctions=sorted(self.junctions), preparation_audit=self.preparation_audit,
                    ownership_audit=self.ownership_audit,
                    coupled_bundle_networks=self.coupled_bundle_networks,
                    junction_interior_links=self.junction_interior_links,
                    junction_origins=self.junction_origins,
                    retained_physical_streets=list(self.retained_physical_streets.values()),
                    retained_physical_bundles=list(self.retained_physical_bundles.values()),
                    clipped_source_context=self.clipped_source_context)


class PreparedNetworkAdapter:
    def __init__(self, design=None, *, ground_passages=False, recover_boundary_junctions=True):
        self.design = design or StreetDesign()
        self.ground_passages = ground_passages
        self.recover_boundary_junctions = recover_boundary_junctions

    @staticmethod
    def items(street):
        item, seen = street.head, set()
        while item is not None:
            if id(item) in seen:
                raise ValueError('Cyclic Street item chain')
            seen.add(id(item))
            yield item
            if item is street.tail:
                return
            item = item.succ
        raise ValueError('Street chain does not reach its tail')

    def read(self, manager):
        snapshot = NetworkSnapshot()
        snapshot.preparation_audit = dict(getattr(manager, 'streetPreparationAudit', {}))
        origins = getattr(manager, 'junctionOrigins', None)
        if origins is not None:
            snapshot.junction_origins = origins.report()
        snapshot.preparation_audit['minor_cycles'] = getattr(manager, 'minorCycleAudit', [])
        active = tuple(manager.iterStreets())
        active_ids = {street.id for street in active}
        retained = tuple(street for sid, street in getattr(manager, 'physicalStreets', {}).items()
                         if sid not in active_ids)
        for street in active + retained:
            entries = []
            for item in self.items(street):
                points = tuple(tuple(float(v) for v in p[:2])+(0.,)
                               for p in getattr(item, 'centerline', ()))
                entries.append(dict(kind=type(item).__name__, id=getattr(item,'id',-1),
                    points=points, category=getattr(item,'category',None),
                    width=getattr(item,'width',None), offset=getattr(item,'offset',None),
                    tags=dict(getattr(item,'tags',{})),
                    source_way_ids=sorted(getattr(item,'sourceWayIds',())),
                    source_segments=getattr(item,'sourceSegments',()),
                    width_origin=getattr(item,'widthOrigin',None)))
                if type(item).__name__ == 'Intersection' and item.isMinor:
                    attachments = []
                    for side, head in (('LEFT', item.leftHead), ('RIGHT', item.rightHead)):
                        connector, seen = head, set()
                        while connector is not None and id(connector) not in seen:
                            seen.add(id(connector))
                            terminal = connector.item.head if connector.leaving else connector.item.tail
                            attachments.append(dict(street_id=connector.item.id, side=side,
                                category=getattr(terminal, 'category', None), connector_id=connector.id,
                                section_id=getattr(terminal, 'id', -1), end=0 if connector.leaving else 1))
                            connector = connector.succ
                    entries[-1]['attachments'] = attachments
                    entries[-1]['location'] = tuple(item.location)
                    entries[-1]['through'] = tuple((s.id, end) for s, end in ((item.pred, 1), (item.succ, 0))
                        if type(s).__name__ == 'Section')
                    if attachments and all(a['category'] not in self.design.categories for a in attachments):
                        entries[-1]['surface_treatment'] = 'DEFERRED_ATTACHMENT'
                        snapshot.issue('MinorIntersection', item.id, 'DEFERRED_COMPONENT',
                            'Through carriageway retained; attached non-carriageway components await their feature providers')
                    elif len(entries[-1]['through']) == 2:
                        entries[-1]['surface_treatment'] = 'NATIVE_JUNCTION'
                elif type(item).__name__ == 'SymLane':
                    # A preserved surface/tag boundary does not necessarily
                    # change the physical carriageway dimensions.
                    a, b = item.pred, item.succ
                    entries[-1]['dimensions'] = ((a.width, a.offset), (b.width, b.offset))
                    if abs(a.width-b.width) < 1e-7 and abs(a.offset-b.offset) < 1e-7:
                        entries[-1]['surface_treatment'] = 'COMPATIBLE_SEAM'
            from .transitions import PreparedTransitionAdapter
            PreparedTransitionAdapter().read(entries)
            reason = self.street_support(entries)
            inventory = snapshot.streets if street.id in active_ids else snapshot.retained_physical_streets
            inventory[street.id] = dict(id=street.id, items=entries,
                status='SUPPORTED' if reason is None else 'UNSUPPORTED', reason=reason,
                endpoints=tuple(getattr(getattr(c, 'intersection', None), 'id', None)
                                for c in (street.pred, street.succ)),
                bundle_id=getattr(getattr(street, 'bundle', None), 'id', None))
            if reason:
                snapshot.issue('Street', street.id, 'STREET_SCOPE', reason)
        from .physical_sources import PhysicalStreetInventory
        PhysicalStreetInventory().review(snapshot)
        from .ownership import PreparedOwnershipReview
        PreparedOwnershipReview().apply(snapshot)
        from .connectivity import SurfaceConnectivityReview
        SurfaceConnectivityReview().apply(snapshot, self.design.categories)
        for bundle in manager.iterBundles():
            members = tuple(s.id for s in bundle.streetsHead)
            tail = tuple(s.id for s in bundle.streetsTail)
            records, reason = [], None
            if len(members) != 2 or set(members) != set(tail):
                reason = 'Current Bundle backend requires two complete members at both ends'
            else:
                for index, (street, location) in enumerate(zip(bundle.streetsHead, bundle.headLocs)):
                    source = snapshot.streets.get(street.id)
                    if source is None:
                        reason = 'Bundle references a Street absent from the prepared ownership inventory'
                        break
                    sections = [i for i in source['items'] if i['kind']=='Section']
                    if not sections:
                        reason = 'Bundle member contains no Section'
                        break
                    start, end = sections[0]['points'][0][:2], sections[-1]['points'][-1][:2]
                    forward, reverse = math.dist(start, location[:2]) < 1e-5, math.dist(end, location[:2]) < 1e-5
                    if forward == reverse:
                        reason = 'Bundle membership does not cover the complete Street endpoints'
                        break
                    tail_index = tail.index(street.id)
                    if tail_index >= len(bundle.tailLocs) or math.dist(
                            end if forward else start, bundle.tailLocs[tail_index][:2]) >= 1e-5:
                        reason = 'Partial Bundle membership requires an explicit extent adapter'
                        break
                    records.append(dict(street_id=street.id, forward=forward,
                        head_order=index, tail_order=tail_index))
                if len(records) != 2 and reason is None:
                    reason = 'Bundle endpoint metadata is incomplete'
            snapshot.bundles[bundle.id] = dict(id=bundle.id, head=members, tail=tail,
                members=records, status='UNSUPPORTED' if reason else 'SUPPORTED', reason=reason,
                head_locations=tuple(tuple(float(v) for v in p[:2]) for p in bundle.headLocs),
                tail_locations=tuple(tuple(float(v) for v in p[:2]) for p in bundle.tailLocs),
                owned_streets=tuple(sorted(sid for sid,s in snapshot.streets.items()
                                          if s['bundle_id']==bundle.id)))
            if reason:
                snapshot.issue('Bundle', bundle.id, 'BUNDLE_SCOPE', reason)
        from .physical_sources import PhysicalBundleInventory
        from .bundle_topology import PreparedBundleTopology
        retained_bundles = {street.bundle.id:street.bundle for street in active + retained
            if getattr(street,'bundle',None) is not None and street.bundle.id not in snapshot.bundles}
        PhysicalBundleInventory().read(snapshot,retained_bundles.values(),PreparedBundleTopology(self.design.categories))
        from .junction_interiors import PreparedJunctionInteriorReview
        interior_review=PreparedJunctionInteriorReview()
        for intersection in manager.majorIntersections.values():
            row = dict(id=intersection.id, location=tuple(intersection.location),
                       order=intersection.order, status='ENDPOINT' if intersection.order < 3 else 'UNSUPPORTED')
            snapshot.intersections[intersection.id] = row
            if intersection.order < 3:
                if intersection.order == 2 and any(hasattr(c.item,'streetsHead') for c in intersection):
                    row.update(status='UNSUPPORTED',reason='Two-connector Bundle terminal needs an explicit transition contract')
                    snapshot.issue('Intersection',intersection.id,'BUNDLE_TERMINAL_TRANSITION',row['reason'])
                continue
            try:
                connectors = junction_adapter._connectors(intersection)
                links,retained=interior_review.read(snapshot,intersection,connectors)
                if links or retained:
                    snapshot.junction_interior_links[intersection.id]=dict(
                        links=[link.report() for link in links],retained=list(retained))
                profiles = {c: inputs.ProfileInput(self.design.profile_rise,
                    inputs.SMOOTH_CROWN if self.design.profile_rise else inputs.FLAT,
                    hasattr(c.item, 'streetsHead')) for c in connectors}
                junction = junction_adapter.from_blosm_intersection(intersection,
                    profiles_by_connector=profiles, category_properties=wayCategoryProps,
                    right_hand_traffic=manager.rightHandTraffic,
                    cell_size=self.design.cell_size, check_scope=False,
                    section_reader=lambda s: self.source_sections(snapshot, s.id))
                for street in junction.streets:
                    row_street = snapshot.streets.get(street.source_id)
                    if row_street is None or row_street['status'] != 'SUPPORTED':
                        raise ValueError(f'Street {street.source_id} is outside the current surface domain')
                # The editable solver sees full paths; only the mouth frame is
                # passed to the local intersection contract for this check.
                local_streets = []
                for index, street in enumerate(junction.streets):
                    member = next(m for c in junction.connectors for m in c.members if m.street_index == index)
                    section = street.sections[0 if member.leaving else -1]
                    points = section.centerline[:2] if member.leaving else section.centerline[-2:]
                    local_streets.append(replace(street, sections=(replace(section, centerline=points, chain_index=0),)))
                inputs.validate_initial_scope(replace(junction, streets=tuple(local_streets)))
                snapshot.junctions[intersection.id] = junction
                row['status'] = 'SUPPORTED'
                row['streets'] = tuple(s.source_id for s in junction.streets)
                row['bundle_ids'] = tuple(c.bundle_source_id for c in junction.connectors if c.item_kind==inputs.BUNDLE_ITEM)
            except (ValueError, KeyError, StopIteration) as exc:
                row['reason'] = str(exc)
                snapshot.issue('Intersection', intersection.id, 'INTERSECTION_SCOPE', str(exc))
        from .minor import MinorJunctionAdapter
        MinorJunctionAdapter(self.design).read(snapshot)
        from .bundle_topology import PreparedBundleTopology
        PreparedBundleTopology(self.design.categories).apply(snapshot)
        from .bundle_junctions import BundleJunctionPlanner, BundleJunctionError, BoundedBundlePlanner, CoupledBundlePlanner
        patch_planner = BundleJunctionPlanner()
        for bundle in snapshot.bundles.values():
            bundle['junction_patches'] = patch_planner.review(bundle['topology'])
            # Plans preserve original Section ownership and both endpoint
            # incidences; they do not enable a partially supported Bundle.
            try:
                bundle['incident_networks'] = [network.report()
                    for network in patch_planner.incident_networks(bundle['topology'])]
                bundle['member_networks'] = [network.report()
                    for network in patch_planner.incident_networks(bundle['topology'],include_member_spans=True)]
            except BundleJunctionError as exc:
                bundle['incident_networks'] = []
                bundle['member_networks'] = []
                bundle['incident_network_error'] = dict(code=exc.code, reason=str(exc))
                snapshot.issue('Bundle', bundle['id'], exc.code, str(exc))
            # Semantic eligibility only. The candidate bounded native provider
            # is not enabled in the released asset by exporting this plan.
            bundle['bounded_network'] = BoundedBundlePlanner(self.design.categories).review(snapshot, bundle['id'])
        snapshot.coupled_bundle_networks = CoupledBundlePlanner(self.design.categories).review(snapshot)
        from .bundles import PreparedBundleReview
        PreparedBundleReview().apply(snapshot)
        context = getattr(manager,'clippedSourceContext',None)
        if context is not None:
            from ..clipped_context import resolve_open_ends
            snapshot.clipped_source_context = resolve_open_ends(context.report(),snapshot.streets.values(),
                                                               snapshot.retained_physical_streets.values())
            if self.recover_boundary_junctions:
                from .boundary_junctions import BoundaryJunctionAdapter
                BoundaryJunctionAdapter(self.design,wayCategoryProps,
                    ground_passages=self.ground_passages).apply(snapshot)
        return snapshot

    def street_support(self, entries):
        for item in entries:
            if item.get('transition_reason'):
                return item['transition_reason']
            if item.get('surface_treatment') in ('DEFERRED_ATTACHMENT', 'NATIVE_JUNCTION', 'COMPATIBLE_SEAM', 'WIDTH_TRANSITION'):
                continue
            if item['kind'] != 'Section':
                return item.get('transition_reason') or 'Unimplemented Street item: '+item['kind']
            if item['category'] not in self.design.categories:
                return 'Reserved component/category: '+str(item['category'])
            reason = GroundSurfacePolicy.reason(item['tags'], building_passages=self.ground_passages)
            if reason:
                return reason
            points = item['points']
            if len(points) < 2 or points[0] == points[-1]:
                return 'Closed or empty Section requires a dedicated path treatment'
            if not item['width'] or item['width'] <= 0:
                return 'Missing or invalid physical width'
        sections = [item for item in entries if item['kind'] == 'Section']
        if len(sections) > 1 and sections[0]['points'][0] == sections[-1]['points'][-1]:
            return 'Closed Street requires a dedicated path treatment'
        return None

    @staticmethod
    def source_sections(snapshot, street_id):
        if street_id not in snapshot.streets:
            raise ValueError(f'Connector references Street {street_id} absent from the prepared ownership inventory')
        return PreparedNetworkAdapter.section_inputs(snapshot.streets[street_id])

    @staticmethod
    def section_inputs(street):
        if street['status'] != 'SUPPORTED':
            raise ValueError(f"Street {street['id']}: {street['reason']}")
        sections = tuple(inputs.SectionInput(item['id'], item['points'],
            float(item['width']), float(item['offset'] or 0.), str(item['category']),
            index, tuple(sorted(item['tags'].items())))
            for index, item in enumerate(i for i in street['items'] if i['kind']=='Section'))
        for a, b in zip(sections, sections[1:]):
            if math.dist(a.centerline[-1], b.centerline[0]) > 1e-6:
                raise ValueError(f"Street {street['id']}: disconnected Sections around attachment")
        return sections
