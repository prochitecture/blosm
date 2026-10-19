"""Feature coverage measured against authored input, not rendered face counts."""
from collections import Counter, defaultdict
import math


class CoverageReport:
    """Retain partial/invalid output as such; never count it as generated length."""
    @staticmethod
    def initialize_junctions(snapshot, selected_ids=None):
        """Reset evidence, distinguishing explicit selection from rejection.

        A prepared selection can include physical sites and their logical
        owners. None means the complete declared input; an empty set means an
        explicit selection with no eligible sites. Neither implies generation.
        """
        for jid,row in snapshot.intersections.items():
            excluded = selected_ids is not None and jid not in selected_ids
            row['surface_status'] = 'DELIBERATELY_EXCLUDED' if excluded else 'UNSUPPORTED'

    @staticmethod
    def bundle_status(bundle, streets, pending):
        """Use evaluated members for active and retained Bundle inventories alike.

        Original-span rendering marks a member GENERATED only after its required
        pavement and complete separator ownership are accepted. The retained
        Bundle's preparation-time `coverage` flag is not surface evidence.
        Interior deferred paths are not carriageway members, so do not substitute
        the broader `owned_streets` inventory for these member declarations.
        """
        members = set(bundle['head']) | set(bundle.get('tail', bundle['head']))
        if not members:
            return pending
        states = {streets.get(sid, {}).get('surface_status', pending) for sid in members}
        if 'surface_status' in bundle:
            # An explicit Bundle result may narrow accepted member evidence; it
            # cannot promote absent or partial member surfaces to GENERATED.
            states.add(bundle['surface_status'])
        return next(iter(states)) if len(states) == 1 else 'PARTIAL'

    def build(self, snapshot):
        lengths = defaultdict(float)
        counts = Counter()
        rows = []
        retained = getattr(snapshot, 'retained_physical_streets', {})
        if set(retained) & set(snapshot.streets):
            raise ValueError('Coverage cannot count one source as both active and retained')
        sources = {**snapshot.streets, **retained}
        for street in (*snapshot.streets.values(), *retained.values()):
            pending = 'PENDING_PHYSICAL_PROVIDER' if street['id'] in retained else 'NOT_EVALUATED'
            status = street.get('surface_status', pending)
            length = sum(math.dist(a, b) for item in street['items']
                for a, b in zip(item.get('points', ()), item.get('points', ())[1:]))
            counts[status] += 1
            lengths[status] += length
            rows.append(dict(street_id=street['id'], status=status,
                input_centerline_length_m=length, reason=street.get('surface_reason') or street['reason'],
                osm_way_ids=sorted({w for item in street['items'] for w in item.get('source_way_ids', ())}),
                endpoints=street['endpoints'], bundle_id=street['bundle_id']))
        junctions = Counter(row.get('surface_status', 'NOT_EVALUATED') for row in snapshot.intersections.values())
        bundles = []
        for bundle in snapshot.bundles.values():
            status = self.bundle_status(bundle, sources, 'NOT_EVALUATED')
            bundles.append(dict(bundle_id=bundle['id'], status=status,
                                streets=list(bundle['head'])))
        retained_bundles = getattr(snapshot,'retained_physical_bundles',{})
        if set(retained_bundles) & set(snapshot.bundles):
            raise ValueError('Coverage cannot count one Bundle as both active and retained')
        for bundle in retained_bundles.values():
            bundles.append(dict(bundle_id=bundle['id'],
                status=self.bundle_status(bundle, sources, 'PENDING_PHYSICAL_PROVIDER'),
                streets=list(bundle['head']),source_role='RETAINED_PHYSICAL'))
        return dict(length_basis='Prepared input centerlines, including retained transition spans; not pavement area or trimmed surface length',
                    street_counts=dict(counts), centerline_length_m=dict(lengths), streets=rows,
                    junction_counts=dict(junctions), bundles=bundles)
