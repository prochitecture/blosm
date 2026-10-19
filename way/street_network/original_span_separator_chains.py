"""Original identities for a separator bounded by several generated corridors.

This planner is independent of Blender and does not accept geometry. Native
asset capabilities select the qualified construction branches explicitly.
No nearest-point matching or input-specific Street IDs are used.
"""
from dataclasses import dataclass, field

from .physical_ports import SectionEnd
from .original_span_separators import OriginalSpanSeparatorPlanner, SeparatorPreparation


@dataclass(frozen=True)
class ChainedBundleSeparator:
    bundle_id: int
    member_street_ids: tuple[int, int]
    head_sites: tuple[int, ...]
    tail_sites: tuple[int, ...]
    member_span_ids: tuple[tuple[tuple[int, int, int], ...], ...]
    head_ends: tuple[SectionEnd, SectionEnd]
    tail_ends: tuple[SectionEnd, SectionEnd]
    intermediate_sites: tuple[tuple[tuple[int, ...], ...], ...]
    # Original Street IDs in each lateral side's head-to-tail traversal.
    # Empty retains the established one-Street-per-side contract.
    member_chains: tuple[tuple[int, ...], ...] = field(default=(), kw_only=True)


@dataclass(frozen=True)
class ChainedSeparatorPreparation(SeparatorPreparation):
    chained: tuple[ChainedBundleSeparator, ...] = ()


@dataclass(frozen=True)
class FreeChainedBundleSeparator(ChainedBundleSeparator):
    """Complete member chains with one attached end and two genuine free ends.

    An empty head_sites or tail_sites tuple denotes the free Bundle end. The
    SectionEnd identities still identify those original endpoints exactly.
    This declaration alone does not qualify a generated boundary or cap.
    """
    attached_at_head: bool

    @property
    def attached_sites(self):
        return self.head_sites if self.attached_at_head else self.tail_sites

    @property
    def free_ends(self):
        return self.tail_ends if self.attached_at_head else self.head_ends


@dataclass(frozen=True)
class FreeChainedSeparatorPreparation(ChainedSeparatorPreparation):
    free_chained: tuple[FreeChainedBundleSeparator, ...] = ()


class ChainedSeparatorPlanner(OriginalSpanSeparatorPlanner):
    """Qualify complete chains; the native graph must qualify their opening."""

    def __init__(self, categories=None, *, longitudinal_chains=False):
        super().__init__(categories)
        self.longitudinal_chains = longitudinal_chains

    def prepare(self, snapshot, component, *, free_ends=False, free_chains=False):
        """Retain established branches, then classify their pending chains.

        Selecting this planner is an explicit capability choice. Free-ended
        chains require their own opt-in; the candidate asset's closed-chain
        capability does not imply support for constructing a new free cap.
        """
        if free_chains and not free_ends:
            raise ValueError('Free-ended chains require the free-end capability')
        basic = OriginalSpanSeparatorPlanner(self.categories).prepare(
            snapshot, component, free_ends=free_ends)
        owners = self.owners(component)
        chained, free_chained, pending = [], [], []
        deferred = list(basic.deferred_interiors)
        for row in basic.pending:
            bundle = snapshot.bundles[row['bundle_id']]
            try:
                declaration = self.candidate(bundle, snapshot, component, owners)
            except ValueError as exc:
                reason = row['reason']+'; chained: '+str(exc)
                if not free_chains:
                    pending.append(dict(row, reason=reason))
                    continue
                try:
                    declaration = self.free_chain_candidate(bundle, snapshot, component, owners)
                except ValueError as free_exc:
                    pending.append(dict(row, reason=reason+'; free chained: '+str(free_exc)))
                    continue
                free_chained.append(declaration)
            else:
                chained.append(declaration)
            interiors = self.deferred_interiors(bundle, snapshot, component)
            if interiors:
                deferred.append(dict(bundle_id=bundle['id'], streets=interiors))
        # Base preparation fields can grow independently of these branches.
        # Positional construction would put chains into selection_ended and
        # free chains into chained after the selection-boundary extension.
        fields = dict(closed=basic.closed, pending=tuple(pending), free_ended=basic.free_ended,
            deferred_interiors=tuple(deferred), selection_ended=basic.selection_ended, chained=tuple(chained))
        return (FreeChainedSeparatorPreparation(**fields, free_chained=tuple(free_chained)) if free_chains else
                ChainedSeparatorPreparation(**fields))

    @staticmethod
    def route_end(route, at, *, allow_free=False):
        span = route.span
        end = (span.start, span.end)[at]
        port = (route.start_port, route.end_port)[at]
        expected = SectionEnd(span.street_id, span.section_ids[-1 if at else 0], at)
        if allow_free and end.owner_kind == 'FREE':
            if end.endpoint != expected or end.junction_id is not None or port is not None:
                raise ValueError('Free chain end must keep its exact original endpoint without a generated port')
            return end
        if (end.endpoint != expected or port is None or port.endpoint != expected or
                port.junction_id != end.junction_id or end.junction_id is None or end.owner_kind == 'FREE'):
            raise ValueError('Chained separator lost an original endpoint or generated port')
        return end

    def candidate(self, bundle, snapshot, component, owners=None):
        return self._candidate(bundle, snapshot, component, owners, free_ends=False)

    def free_chain_candidate(self, bundle, snapshot, component, owners=None):
        return self._candidate(bundle, snapshot, component, owners, free_ends=True)

    def _candidate(self, bundle, snapshot, component, owners, *, free_ends):
        if bundle.get('topology', {}).get('chains') and self.longitudinal_chains:
            from .original_span_longitudinal_separators import LongitudinalSeparatorPreparation
            return LongitudinalSeparatorPreparation(self).candidate(
                bundle, snapshot, component, owners, free_ends=free_ends)
        owners = self.owners(component) if owners is None else owners
        members, ids = self.members(bundle, component, snapshot, split_extents=True)
        junctions = bundle['topology'].get('junctions', ())
        declared = {(j['street_id'], j['junction_id']): j for j in junctions}
        if len(declared) != len(junctions):
            raise ValueError('Duplicate intermediate Bundle junction declaration')
        used = set(); chains = []; middle = []; terminals = []; ends = []
        for member in members:
            sid = member['street_id']
            routes = sorted((r for r in component.network.routes if r.span.street_id == sid),
                            key=lambda r: r.span.first_control)
            kinds = {'CORRIDOR', 'APPROACH'} if free_ends else {'CORRIDOR'}
            if not routes or any(r.kind not in kinds or not r.span.in_scope for r in routes):
                raise ValueError('Chained members require accepted in-scope routes throughout')
            sections = tuple(s['id'] for s in snapshot.streets[sid]['items'] if s['kind'] == 'Section')
            actual = tuple(s for r in routes for s in r.span.section_ids)
            spans = [r.span for r in routes]
            count = spans[0].source_control_count
            if (actual != sections or len(set(actual)) != len(actual) or
                    spans[0].first_control != 0 or spans[-1].last_control != count-1 or
                    any(s.source_control_count != count or not 0 <= s.first_control < s.last_control < count
                        for s in spans) or any(a.last_control != b.first_control for a,b in zip(spans,spans[1:]))):
                raise ValueError('Corridor chain must cover the complete original Street without gaps or repeats')
            authored_ends = [(self.route_end(r, 0, allow_free=free_ends),
                              self.route_end(r, 1, allow_free=free_ends)) for r in routes]
            if free_ends:
                free_positions = [(i, at) for i, pair in enumerate(authored_ends)
                                  for at, end in enumerate(pair) if end.owner_kind == 'FREE']
                if len(free_positions) != 1 or free_positions[0] not in ((0, 0), (len(routes)-1, 1)):
                    raise ValueError('Each free chain needs exactly one free endpoint at a complete Street end')
                free_route = free_positions[0][0]
                if any(r.kind != ('APPROACH' if i == free_route else 'CORRIDOR')
                       for i, r in enumerate(routes)):
                    raise ValueError('A free chain must end with one approach and otherwise contain corridors')
            sites = []
            for pair, following in zip(authored_ends, authored_ends[1:]):
                previous, after = pair[1], following[0]
                key = (sid, previous.junction_id)
                declaration = declared.get(key)
                if previous.junction_id != after.junction_id or declaration is None:
                    raise ValueError('Every chain join needs its exact original Bundle junction')
                expected = {(previous.endpoint.section_id, 1), (after.endpoint.section_id, 0)}
                through = [tuple(x) for x in declaration.get('through', ())]
                if len(through) != 2 or set(through) != expected:
                    raise ValueError('Intermediate through endpoints differ from the original route chain')
                if key in used or previous.junction_id not in owners:
                    raise ValueError('Intermediate junction is repeated or lacks an accepted surface owner')
                used.add(key); sites.append(owners[previous.junction_id])
            first, last = authored_ends[0][0], authored_ends[-1][1]
            ordered = routes if member['forward'] else routes[::-1]
            if not member['forward']:
                first, last = last, first
                sites.reverse()
            terminal_owners = []
            for endpoint in (first, last):
                if free_ends and endpoint.owner_kind == 'FREE':
                    terminal_owners.append(())
                elif endpoint.junction_id in owners:
                    terminal_owners.append(owners[endpoint.junction_id])
                else:
                    raise ValueError('Attached chained terminals need accepted surface owners')
            extents = member['spans']
            if len(extents) != len(ordered):
                raise ValueError('Original Bundle extents and generated corridor chain disagree')
            for ordinal, (extent, route) in enumerate(zip(extents, ordered)):
                # Topology mapping has already classified completeness. Verify
                # its ordered Section inventory and end labels again here.
                section_rows = extent.get('sections', ())
                expected_sections = route.span.section_ids if member['forward'] else route.span.section_ids[::-1]
                if (extent.get('ordinal') != ordinal or
                        tuple(s['start']['section_id'] for s in section_rows) != expected_sections or
                        any(s['start']['section_id'] != s['end']['section_id'] for s in section_rows)):
                    raise ValueError('Mapped Bundle extent does not match the directed original Sections')
                route_ends = authored_ends[routes.index(route)]
                if not member['forward']: route_ends = route_ends[::-1]
                for at, label in enumerate(('start', 'end')):
                    edge = ordinal == 0 and at == 0 or ordinal == len(extents)-1 and at == 1
                    endpoint = extent.get(label, {})
                    expected_label = ('HEAD' if at == 0 else 'TAIL') if edge else None
                    expected_site = None if edge else route_ends[at].junction_id
                    if endpoint.get('bundle_end') != expected_label or endpoint.get('junction_id') != expected_site:
                        raise ValueError('Bundle extent ends disagree with the original chain joins')
            chains.append(tuple((sid, r.span.section_ids[0], r.span.section_ids[-1]) for r in ordered))
            middle.append(tuple(sites))
            terminals.append(tuple(terminal_owners))
            ends.append((first.endpoint, last.endpoint))
        if not used or used != set(declared):
            raise ValueError('Chained branch needs complete intermediate junction coverage')
        if terminals[0] != terminals[1] or terminals[0][0] == terminals[0][1]:
            raise ValueError('Chained members must agree on attached owners and free-end placement')
        interior = [site for sites in middle for site in sites]
        if len(interior) != len(set(interior)) or set(interior).intersection(terminals[0]):
            raise ValueError('Shared or repeated intermediate surface owners need multiple separator regions')
        args = (bundle['id'], ids, *terminals[0], tuple(chains),
                tuple(pair[0] for pair in ends), tuple(pair[1] for pair in ends), tuple(middle))
        return (FreeChainedBundleSeparator(*args, bool(terminals[0][0])) if free_ends else
                ChainedBundleSeparator(*args))
