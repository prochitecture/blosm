"""Semantic dispatch for separators in the original-span provider.

This prepares identities only. Eligibility is not a claim that the live
boundary or resulting surface passes the native geometry checks. Other Bundle
configurations remain explicit pending work, never silently treated as filled.
"""
from dataclasses import dataclass

from .physical_ports import SectionEnd


@dataclass(frozen=True)
class ClosedBundleSeparator:
    bundle_id: int
    member_street_ids: tuple[int, int]
    head_sites: tuple[int, ...]
    tail_sites: tuple[int, ...]


@dataclass(frozen=True)
class FreeEndBundleSeparator:
    """Two complete members attached to one owner, with two real free ends.

    Member order follows the original Bundle head order. ``attached_at_head``
    says which Bundle end has the generated patch; it does not reverse either
    authored Street. A native boundary builder still has to qualify the space
    between the two generated approaches before creating a separator.
    """
    bundle_id: int
    member_street_ids: tuple[int, int]
    attached_sites: tuple[int, ...]
    attached_at_head: bool
    free_ends: tuple[SectionEnd, SectionEnd]


@dataclass(frozen=True)
class SelectionEndBundleSeparator:
    """One attached owner and explicitly bounded member extents.

    At least one terminal is a declared selection boundary. A terminal may
    instead be genuinely free; both original endpoint roles are retained.
    """
    bundle_id: int
    member_street_ids: tuple[int, int]
    attached_sites: tuple[int, ...]
    attached_at_head: bool
    terminal_ends: tuple


@dataclass(frozen=True)
class SeparatorPreparation:
    closed: tuple[ClosedBundleSeparator, ...]
    pending: tuple[dict, ...]
    free_ended: tuple[FreeEndBundleSeparator, ...] = ()
    deferred_interiors: tuple[dict, ...] = ()
    selection_ended: tuple[SelectionEndBundleSeparator, ...] = ()


class OriginalSpanSeparatorPlanner:
    """Use declared Bundle order and generated owners, independent of file IDs."""

    def __init__(self, categories=None):
        # A missing scope keeps the earlier conservative contract. An explicit
        # category scope can defer modal paths, but never an active carriageway.
        self.categories = None if categories is None else frozenset(categories)
        if self.categories is not None and (not self.categories or
                any(not isinstance(category, str) or not category for category in self.categories)):
            raise ValueError('Separator category scope must contain explicit category names')

    def deferred_interiors(self, bundle, snapshot, component):
        ids = tuple(bundle.get('topology', {}).get('interior_streets', ()))
        if not ids:
            return ()
        if self.categories is None:
            raise ValueError('Bundle interior Streets require an explicit feature category scope')
        if len(ids) != len(set(ids)) or any(type(sid) is not int or sid < 0 for sid in ids):
            raise ValueError('Bundle interior Street identities must be unique nonnegative integers')
        declared = {route.span.street_id for route in component.network.routes if route.span.in_scope}
        rows = []
        for sid in ids:
            street = snapshot.streets.get(sid)
            if street is None:
                raise ValueError('Bundle interior Street is missing from the original inventory')
            sections = [item for item in street.get('items', ()) if item.get('kind') == 'Section']
            categories = {section.get('category') for section in sections}
            if not categories or any(not isinstance(category, str) or not category for category in categories):
                raise ValueError('Bundle interior Street has incomplete category evidence')
            if categories.intersection(self.categories) or sid in declared:
                raise ValueError('Active Bundle interior carriageway needs a separator-region construction')
            rows.append(dict(street_id=sid, section_ids=tuple(item['id'] for item in sections),
                             categories=tuple(sorted(categories))))
        return tuple(rows)

    @staticmethod
    def owners(component):
        groups = [(site,) for site in component.network.ordinary_site_ids]
        groups += [tuple(c.site_ids) for c in component.network.compounds]
        groups += [tuple(c.physical_junction_ids) for c in component.network.logical_clusters]
        result = {}
        for group in groups:
            owner = tuple(sorted(group))
            if not owner or len(set(owner)) != len(owner):
                raise ValueError('Separator terminal owner has missing or duplicate sites')
            for site in owner:
                if site in result:
                    raise ValueError('Separator terminal site has multiple surface owners')
                result[site] = owner
        return result

    def members(self, bundle, component, snapshot, *, split_extents=False):
        bid = bundle['id']
        if type(bid) is not int or bid < 0:
            raise ValueError('Bundle identity must be a nonnegative integer')
        topology = bundle.get('topology', {})
        if topology.get('status') != 'MAPPED':
            raise ValueError('Bundle extents do not have an accepted original topology mapping')
        if topology.get('chains'):
            raise ValueError('Separator construction across multiple original Streets per lateral chain remains unsupported')
        members = sorted(topology.get('members', ()), key=lambda m: m['head_order'])
        if len(members) != 2:
            raise ValueError('Two-member branch needs exactly two ordered members')
        self.deferred_interiors(bundle, snapshot, component)
        features = {'INTERIOR_STREETS'} | ({'INTERNAL_JUNCTIONS'} if split_extents else set())
        if (topology.get('junctions') and not split_extents) or set(topology.get('required_features', ())) - features:
            raise ValueError('Bundle attachments or split extents need another construction branch')
        if ([m['head_order'] for m in members] != [0, 1] or
                [m['tail_order'] for m in members] != [0, 1]):
            raise ValueError('Member order must agree at both Bundle ends')
        ids = tuple(m['street_id'] for m in members)
        if len(set(ids)) != 2 or any(type(s) is not int or s < 0 for s in ids):
            raise ValueError('A separator needs two distinct original Street IDs')
        if not set(ids) <= set(component.ports.source_streets):
            raise ValueError('Both complete Bundle members must belong to this network component')
        for member in members:
            if (not member['complete'] or not member['spans'] or
                    (not split_extents and len(member['spans']) != 1) or type(member['forward']) is not bool):
                raise ValueError('Branch needs complete, directed member extents')
        return members, ids

    @staticmethod
    def complete_route(member, snapshot, component, kind):
        sid = member['street_id']
        routes = [r for r in component.network.routes if r.span.street_id == sid]
        if len(routes) != 1 or routes[0].kind != kind:
            label = 'closed member needs one complete generated corridor' if kind == 'CORRIDOR' else \
                    'free-ended member needs one complete generated approach'
            raise ValueError('Each '+label)
        span = routes[0].span
        sections = tuple(s['id'] for s in snapshot.streets[sid]['items'] if s['kind'] == 'Section')
        if (span.section_ids != sections or span.first_control != 0 or
                span.last_control != span.source_control_count-1):
            raise ValueError('Generated route does not cover the complete original Street')
        return routes[0]

    def candidate(self, bundle, snapshot, component, owners):
        members, ids = self.members(bundle, component, snapshot)
        terminals = []
        for member in members:
            span = self.complete_route(member, snapshot, component, 'CORRIDOR').span
            ends = (span.start.junction_id, span.end.junction_id)
            if not member['forward']:
                ends = ends[::-1]
            if any(end not in owners for end in ends):
                raise ValueError('Both separator terminals need an accepted generated surface owner')
            terminals.append(tuple(owners[end] for end in ends))
        if terminals[0] != terminals[1] or terminals[0][0] == terminals[0][1]:
            raise ValueError('Members must connect the same two distinct terminal surface owners')
        return ClosedBundleSeparator(bundle['id'], ids, *terminals[0])

    def free_candidate(self, bundle, snapshot, component, owners):
        members, ids = self.members(bundle, component, snapshot)
        terminals, free_ends = [], []
        for member in members:
            route = self.complete_route(member, snapshot, component, 'APPROACH')
            span = route.span
            ends = (span.start, span.end)
            ports = (route.start_port, route.end_port)
            if not member['forward']:
                ends, ports = ends[::-1], ports[::-1]
            attached = [i for i, end in enumerate(ends) if end.junction_id in owners]
            if len(attached) != 1:
                raise ValueError('Free-ended member needs exactly one generated terminal owner')
            at = attached[0]
            fixed, free = ends[at], ends[1-at]
            port = ports[at]
            if (fixed.owner_kind == 'FREE' or port is None or
                    port.endpoint != fixed.endpoint or port.junction_id != fixed.junction_id):
                raise ValueError('Attached separator end lost its original generated port')
            if free.owner_kind != 'FREE' or free.junction_id is not None or ports[1-at] is not None:
                raise ValueError('Separator opening needs explicit free ends, not unresolved or category connections')
            expected = SectionEnd(span.street_id, span.section_ids[-1] if free is span.end else
                                  span.section_ids[0], int(free is span.end))
            if free.endpoint != expected:
                raise ValueError('Free end must reference the original complete Street endpoint')
            terminals.append((at, owners[fixed.junction_id]))
            free_ends.append(free.endpoint)
        if terminals[0] != terminals[1]:
            raise ValueError('Both members need the same attached Bundle end and generated owner')
        return FreeEndBundleSeparator(bundle['id'], ids, terminals[0][1],
                                      terminals[0][0] == 0, tuple(free_ends))

    def prepare(self, snapshot, component, *, free_ends=False):
        """Free-ended dispatch is opt-in until its native provider is installed."""
        owners = self.owners(component)
        closed, pending, free_ended, deferred_interiors = [], [], [], []
        selected = set(component.ports.source_streets)
        whole={row.bundle_id for row in getattr(component,'free_bundles',())}
        for bid, bundle in sorted(snapshot.bundles.items()):
            # Kind 8 supplies both pavement and all automatic grass islands.
            # It must not also request a second separator-only owner.
            if bid in whole:continue
            members = bundle.get('topology', {}).get('members', ())
            related = {m['street_id'] for m in members} | set(bundle.get('owned_streets', ()))
            if not selected.intersection(related):
                continue
            if bundle['id'] != bid:
                raise ValueError('Bundle inventory key differs from its original identity')
            try:
                closed.append(self.candidate(bundle, snapshot, component, owners))
                deferred = self.deferred_interiors(bundle, snapshot, component)
                if deferred: deferred_interiors.append(dict(bundle_id=bid, streets=deferred))
            except ValueError as exc:
                if free_ends:
                    try:
                        free_ended.append(self.free_candidate(bundle, snapshot, component, owners))
                        deferred = self.deferred_interiors(bundle, snapshot, component)
                        if deferred: deferred_interiors.append(dict(bundle_id=bid, streets=deferred))
                        continue
                    except ValueError as free_exc:
                        reason = f'Closed: {exc}; free-ended: {free_exc}'
                else:
                    reason = str(exc)
                pending.append(dict(bundle_id=bid, status='PENDING', reason=reason))
        return SeparatorPreparation(tuple(closed), tuple(pending), tuple(free_ended), tuple(deferred_interiors))
