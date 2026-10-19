"""Original-incidence candidates for an internal multi-Bundle separator.

A graph cycle is not evidence of an island. This planner only declares an
unbranched cycle of complete original Bundle members inside one surface owner.
The native construction must still check the live footprint, profiles, available
space and pavement partition. No coordinates, generated faces or coverage claims
are returned here.
"""
from collections import defaultdict
from dataclasses import dataclass

from .original_span_separators import OriginalSpanSeparatorPlanner
from .physical_ports import SectionEnd


@dataclass(frozen=True)
class InternalCycleMember:
    bundle_id: int
    span_id: tuple[int, int, int]
    from_end: SectionEnd
    to_end: SectionEnd
    from_junction: int
    to_junction: int


@dataclass(frozen=True)
class InternalCycleRegion:
    logical_junction_ids: tuple[int, ...]
    physical_junction_ids: tuple[int, ...]
    bundle_ids: tuple[int, ...]
    members: tuple[InternalCycleMember, ...]
    deferred_interiors: tuple[dict, ...] = ()


@dataclass(frozen=True)
class InternalCyclePreparation:
    regions: tuple[InternalCycleRegion, ...]
    pending: tuple[dict, ...]
    not_applicable: tuple[dict, ...] = ()


class InternalCycleRegionPlanner(OriginalSpanSeparatorPlanner):
    """Prepare the initial straight, two-member Bundle cycle branch.

    Authored Street direction and Bundle lateral order remain independent of
    cycle traversal. The traversal starts at the lowest original junction ID,
    then follows the lowest unused Street ID; source/table order cannot change
    it. Moving controls does not change these original incidence identities.
    """

    def prepare(self, snapshot, component):
        owners = self.owners(component)
        regions, pending, not_applicable, seen = [], [], [], set()
        # A cycle of ordinary intersection Streets is not, by itself, an
        # internal Bundle separator. Check both directions of the original
        # membership relation so a stale missing Street.bundle_id cannot hide
        # a declaration on the Bundle side. Index once for the whole component.
        bundle_references = self.bundle_references(snapshot)
        for cluster in sorted(component.network.logical_clusters,
                              key=lambda row: row.logical_junction_id):
            if cluster.source_cycle_rank == 0:
                continue
            key = tuple(sorted(cluster.physical_junction_ids))
            if key in seen:
                raise ValueError('An internal cycle cannot have duplicate surface owners')
            seen.add(key)
            try:
                source_ids = self.non_bundle_cycle(snapshot, component, cluster, bundle_references)
                if source_ids is not None:
                    not_applicable.append(dict(logical_junction_ids=tuple(cluster.owner_ids),
                        physical_junction_ids=key, street_ids=source_ids,
                        code='NO_INTERNAL_BUNDLE_MEMBERS',
                        reason='No original Bundle members or interior Streets in this cycle; '
                               'a centerline cycle alone does not declare a Bundle separator'))
                    continue
                regions.append(self.candidate(snapshot, component, cluster, owners))
            except ValueError as exc:
                pending.append(dict(logical_junction_ids=tuple(cluster.owner_ids),
                                    physical_junction_ids=key, reason=str(exc)))
        return InternalCyclePreparation(tuple(regions), tuple(pending), tuple(not_applicable))

    @staticmethod
    def bundle_references(snapshot):
        """Known original member/interior roles, including retained Bundles."""
        references = defaultdict(set)
        for bid, bundle in snapshot.bundles.items():
            topology = bundle.get('topology', {})
            ids = set(bundle.get('owned_streets', ()))
            for key in ('head', 'tail'):
                ids.update(bundle.get(key, ()))
            for records in (bundle.get('members', ()), topology.get('members', ())):
                ids.update(row['street_id'] for row in records)
            ids.update(topology.get('interior_streets', ()))
            for chain in topology.get('chains', ()):
                ids.update(chain)
            for sid in ids:
                references[sid].add(bid)
        return references

    @staticmethod
    def non_bundle_cycle(snapshot, component, cluster, references):
        """Report inapplicability without granting any surface coverage.

        Missing/unsupported sources and conflicting ownership remain pending.
        This makes no assertion about whether another feature (for example a
        roundabout) needs an island; it only scopes the Bundle-separator branch.
        """
        ids = tuple(sorted(span.street_id for span in cluster.internal_spans))
        if not ids or len(ids) != len(set(ids)):
            raise ValueError('An internal cycle needs distinct original source Streets')
        if any(sid not in component.ports.source_streets for sid in ids):
            raise ValueError('An internal cycle source is absent from its component')
        for sid in ids:
            source = snapshot.streets.get(sid)
            if (type(sid) is not int or sid < 0 or source is None or
                    source.get('id') != sid or source.get('status') != 'SUPPORTED'):
                raise ValueError('Every internal cycle source needs supported original ownership')
            if source.get('bundle_id') is not None or references.get(sid):
                return None
        return ids

    def candidate(self, snapshot, component, cluster, owners=None):
        if self.categories is None:
            raise ValueError('Internal cycles require an explicit carriageway category scope')
        owner_ids = tuple(sorted(cluster.owner_ids))
        site_ids = tuple(sorted(cluster.physical_junction_ids))
        if (not owner_ids or len(set(owner_ids)) != len(owner_ids) or
                any(type(v) is not int or v < 0 for v in (*owner_ids, *site_ids)) or
                len(site_ids) != len(set(site_ids)) or not 3 <= len(site_ids) <= 64):
            raise ValueError('Internal cycle needs unique original logical and physical owners')
        owners = self.owners(component) if owners is None else owners
        if any(owners.get(jid) != site_ids for jid in site_ids):
            raise ValueError('The whole internal cycle must belong to one generated surface owner')
        if cluster.source_cycle_rank != 1 or cluster.footprint_spans:
            raise ValueError('This branch needs one simple cycle without additional footprint spans')
        sites = {site.junction_id: site for site in component.ports.sites}
        if len(sites) != len(component.ports.sites) or not set(site_ids) <= sites.keys():
            raise ValueError('Internal cycle lost its complete original physical sites')
        spans = {}
        for span in component.ports.spans:
            identity = (span.street_id, span.section_ids[0], span.section_ids[-1])
            if identity in spans:
                raise ValueError('Duplicate original span identity')
            spans[identity] = span
        adjacency, records = defaultdict(list), {}
        for span in cluster.internal_spans:
            sid = span.street_id
            source = snapshot.streets.get(sid)
            if (type(sid) is not int or sid < 0 or sid in records or source is None or
                    source.get('id') != sid or source.get('status') != 'SUPPORTED' or
                    sid not in component.ports.source_streets):
                raise ValueError('Every cycle member needs one supported original source Street')
            sections = [s for s in source.get('items', ()) if s.get('kind') == 'Section']
            if (len(source['items']) != 1 or len(sections) != 1 or
                    len(sections[0].get('points', ())) != 2 or sections[0].get('offset') != 0 or
                    sections[0].get('category') not in self.categories or
                    tuple(span.section_ids) != (sections[0]['id'],) or not span.in_scope or
                    (span.first_control, span.last_control, span.source_control_count) != (0, 1, 2)):
                raise ValueError('Internal cycle needs complete two-control, zero-offset carriageway members')
            section = sections[0]['id']
            identity = (sid, section, section)
            if spans.get(identity) != span:
                raise ValueError('Cluster span differs from the original physical span inventory')
            bid = source.get('bundle_id')
            if type(bid) is not int or bid not in snapshot.bundles:
                raise ValueError('Every cycle member must retain its original Bundle identity')
            ends = (span.start, span.end)
            for at, end in enumerate(ends):
                if (end.owner_kind != 'PHYSICAL' or end.junction_id not in site_ids or
                        end.endpoint != SectionEnd(sid, section, at) or
                        source.get('endpoint_space') != 'PHYSICAL_JUNCTION' or
                        len(source.get('endpoints', ())) != 2 or
                        source['endpoints'][at] != end.junction_id):
                    raise ValueError('Cycle member lost its exact original Section endpoint')
                ports = sites[end.junction_id].ports
                if (len({p.endpoint for p in ports}) != len(ports) or
                        any(p.junction_id != end.junction_id for p in ports) or
                        sum(p.endpoint == end.endpoint for p in ports) != 1):
                    raise ValueError('Original physical incidence must contain each cycle end exactly once')
                adjacency[end.junction_id].append(sid)
            if ends[0].junction_id == ends[1].junction_id:
                raise ValueError('An internal cycle member cannot return to its starting junction')
            records[sid] = (bid, identity, ends)
        if set(adjacency) != set(site_ids) or any(len(v) != 2 for v in adjacency.values()):
            raise ValueError('Need one unbranched cycle covering every original owner site')

        bundle_ids = tuple(sorted({row[0] for row in records.values()}))
        deferred = []
        for bid in bundle_ids:
            bundle = snapshot.bundles[bid]
            if (bundle.get('id') != bid or bundle.get('status') != 'SUPPORTED' or
                    bundle.get('topology', {}).get('bundle_id') != bid):
                raise ValueError('Internal cycle lost its supported original Bundle mapping')
            members, ids = self.members(bundle, component, snapshot)
            expected = {sid for sid, row in records.items() if row[0] == bid}
            if set(ids) != expected or set(bundle.get('owned_streets', ())) != expected:
                raise ValueError('Each participating Bundle must keep both complete original members')
            for member in members:
                sid = member['street_id']; section = records[sid][1][1]
                extent, = member['spans']
                start, end = ((0., 1.) if member['forward'] else (1., 0.))
                expected_sections = [dict(start=dict(section_id=section, segment_index=0, factor=start),
                                          end=dict(section_id=section, segment_index=0, factor=end))]
                if (extent.get('ordinal') != 0 or tuple(extent.get('sections', ())) != tuple(expected_sections) or
                        extent.get('start') != dict(bundle_end='HEAD', junction_id=None) or
                        extent.get('end') != dict(bundle_end='TAIL', junction_id=None)):
                    raise ValueError('Bundle extent must cover the exact complete authored Section')
            extra = self.deferred_interiors(bundle, snapshot, component)
            if extra:
                deferred.append(dict(bundle_id=bid, streets=extra))
        for bid, bundle in snapshot.bundles.items():
            declared = {m['street_id'] for m in bundle.get('topology', {}).get('members', ())}
            declared.update(bundle.get('owned_streets', ()))
            if any(records[sid][0] != bid for sid in declared.intersection(records)):
                raise ValueError('An internal cycle member has competing Bundle ownership')

        at = start = min(site_ids)
        used, ordered = set(), []
        while len(used) < len(records):
            choices = sorted(sid for sid in adjacency[at] if sid not in used)
            if not choices:
                raise ValueError('Disconnected cycles cannot be treated as one separator region')
            sid = choices[0]; used.add(sid)
            bid, identity, pair = records[sid]
            first, last = pair if pair[0].junction_id == at else pair[::-1]
            ordered.append(InternalCycleMember(bid, identity, first.endpoint, last.endpoint,
                                               first.junction_id, last.junction_id))
            at = last.junction_id
            if at == start and len(used) != len(records):
                raise ValueError('Disconnected cycles cannot be treated as one separator region')
        if at != start:
            raise ValueError('Original internal cycle did not close')
        return InternalCycleRegion(owner_ids, site_ids, bundle_ids, tuple(ordered), tuple(deferred))
