"""Semantic boundary roles for side junctions on prepared Bundle members.

This planner assigns profiles and span incidences, never surface coordinates.
The native kernel and span/median assembler have separate capability checks.
In particular, a successful plan must not mark the whole Bundle GENERATED.
"""
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class HalfProfileEnd:
    section_id: int
    end: int
    half_profile: int  # Crown at outward-frame left (-1) or right (+1).


@dataclass(frozen=True)
class BundleSideJunction:
    bundle_id: int
    junction_id: int
    street_id: int
    partner_street_id: int
    inner_u: float  # Left=0/right=1 in the authored through-Street direction.
    before_span: int  # Bundle traversal order, not authored Street direction.
    after_span: int
    through: tuple[HalfProfileEnd, ...]
    attachments: tuple[tuple[int, int, int], ...]  # Street, Section, end.

    def report(self):
        return asdict(self)


class BundleJunctionError(ValueError):
    def __init__(self, code, message):
        self.code = code
        super().__init__(message)


@dataclass(frozen=True)
class IncidentSection:
    """One authored Section, with independent owners of its two ends.

    None means unbound in this component, not a claim that the source end is
    free in the complete street network. Profile side is in authored direction.
    """
    section_id: int
    street_id: int
    half_profile: bool
    inner_u: float
    start_junction: int | None
    end_junction: int | None


@dataclass(frozen=True)
class IncidentSeam:
    street_id: int
    before_section: int
    after_section: int


@dataclass(frozen=True)
class BundleIncidentNetwork:
    bundle_id: int
    junctions: tuple[BundleSideJunction, ...]
    sections: tuple[IncidentSection, ...]
    seams: tuple[IncidentSeam, ...] = ()
    member_spans: tuple[tuple[int, int], ...] = ()  # Street ID, derived span ordinal.

    def report(self):
        return dict(status='PLANNED_MEMBER_NETWORK' if self.member_spans else 'PLANNED_INCIDENT_NETWORK', **asdict(self))


@dataclass(frozen=True)
class FreeBundleEnd:
    bundle_id: int
    side: int  # Bundle HEAD=0 / TAIL=1, independent of Street authoring direction.
    member_ends: tuple[tuple[int, int], ...]  # Complete original Section ID, authored end.


@dataclass(frozen=True)
class BoundedBundleNetwork:
    """Complete source ownership between paired junctions or explicit free ends.

    A plan is topological eligibility, not geometric validation or generated
    coverage. The native provider must still validate its terminal frames,
    curvature, offsets, separator gap and mesh before accepting the component.
    """
    bundle_id: int
    terminal_junctions: tuple[int | None, int | None]  # HEAD, TAIL; None requires an explicit free record.
    junctions: tuple[BundleSideJunction, ...]
    sections: tuple[IncidentSection, ...]
    seams: tuple[IncidentSeam, ...]
    member_spans: tuple[tuple[int, int], ...]
    deferred_streets: tuple[int, ...] = ()  # Other feature domains; not generated.
    free_ends: tuple[FreeBundleEnd, ...] = ()
    ordinary_junctions: tuple[int, ...] = ()  # Intermediate junctions on complete ordinary branch Streets.

    def report(self):
        return dict(status='PLANNED_BOUNDED_NETWORK', **asdict(self))


class BundleJunctionPlanner:
    """Use declared lateral order and attachment side, never proximity.

    Current patch domain: one outer-side vehicle junction on either member of
    a two-member Bundle. Inner attachments/openings require a different patch
    provider. Coordinates remain with the original Section sources.
    """
    def plan(self, topology, junction):
        if topology.get('status') != 'MAPPED':
            raise BundleJunctionError('UNMAPPED_EXTENTS', 'Bundle extents must be mapped before planning a patch')
        members = topology['members']
        if len(members) != 2:
            raise BundleJunctionError('MEMBER_COUNT', 'This patch provider requires two ordered outer members')
        if {m['head_order'] for m in members} != {0, 1} or any(m['head_order'] != m['tail_order'] for m in members):
            raise BundleJunctionError('MEMBER_ORDER', 'Bundle members must retain their declared lateral order')
        matched = [m for m in members if m['street_id'] == junction['street_id']]
        if len(matched) != 1:
            raise BundleJunctionError('MEMBER_REFERENCE', 'Junction must belong to exactly one outer member')
        member = matched[0]
        partner = next(m for m in members if m is not member)
        # The existing Bundle convention is head_order 0 toward +right in
        # Bundle direction. Reversing a Street reverses its local right vector.
        inner_u = float((member['head_order'] == 0) == member['forward'])
        outer_side = 'LEFT' if inner_u else 'RIGHT'
        vehicles = [a for a in junction['attachments'] if a['carriageway']]
        if not vehicles:
            raise BundleJunctionError('NO_CARRIAGEWAY_ATTACHMENT', 'No vehicle surface joins at this junction')
        if any(a.get('side') not in ('LEFT', 'RIGHT') for a in vehicles):
            raise BundleJunctionError('MISSING_SIDE', 'Every carriageway attachment needs its declared side')
        if any(a['side'] != outer_side for a in vehicles):
            raise BundleJunctionError('INNER_ATTACHMENT', 'A median-side attachment requires a separator opening patch')
        jid = junction['junction_id']
        before = [s for s in member['spans'] if s['end']['junction_id'] == jid]
        after = [s for s in member['spans'] if s['start']['junction_id'] == jid]
        if len(before) != 1 or len(after) != 1 or before[0]['ordinal']+1 != after[0]['ordinal']:
            raise BundleJunctionError('SPAN_INCIDENCE', 'The side patch needs consecutive member spans on both sides')
        refs = junction['through']
        if len(refs) != 2 or {end for _sid, end in refs} != {0, 1} or refs[0][0] == refs[1][0]:
            raise BundleJunctionError('THROUGH_INCIDENCE', 'Two distinct authored through Section ends are required')
        expected = {(before[0]['sections'][-1]['end']['section_id'], 1 if member['forward'] else 0),
                    (after[0]['sections'][0]['start']['section_id'], 0 if member['forward'] else 1)}
        if {tuple(ref) for ref in refs} != expected:
            raise BundleJunctionError('THROUGH_INCIDENCE', 'Through Sections do not match the adjacent member spans')
        sign = 1 if inner_u else -1
        return BundleSideJunction(topology['bundle_id'], jid, member['street_id'], partner['street_id'],
            inner_u, before[0]['ordinal'], after[0]['ordinal'],
            tuple(HalfProfileEnd(sid, end, sign if end == 0 else -sign) for sid, end in refs),
            tuple((a['street_id'], a['section_id'], a['end']) for a in vehicles))

    def review(self, topology):
        """Report each local capability without enabling native Bundle output."""
        result = []
        for junction in topology.get('junctions', ()):
            try:
                patch = self.plan(topology, junction)
                result.append(dict(status='PLANNED', **patch.report()))
            except BundleJunctionError as exc:
                result.append(dict(status='UNSUPPORTED', junction_id=junction['junction_id'],
                                   street_id=junction['street_id'], code=exc.code, reason=str(exc)))
        return result

    def incident_networks(self, topology, *, include_member_spans=False):
        """Group supported patches through shared *semantic* Section IDs.

        Emit each incident Section once, including both endpoint bindings.
        With include_member_spans, declared ranges between supported patches
        and complete member endpoints also contribute their Section seams.
        Partial terminal ranges, unsupported junctions and the opposite member
        are not inferred or silently declared generated by this local planner.
        Call review() for the explicit reasons for omitted unsupported patches.
        """
        patches = []
        for junction in topology.get('junctions', ()):
            try:
                patches.append(self.plan(topology, junction))
            except BundleJunctionError:
                continue
        by_junction = {}
        records = {}
        incident = {}
        for patch in patches:
            jid = patch.junction_id
            if jid in by_junction:
                raise BundleJunctionError('DUPLICATE_JUNCTION', 'An incident junction has more than one patch owner')
            by_junction[jid] = patch
            refs = [(r.section_id, r.end, patch.street_id, True, patch.inner_u) for r in patch.through]
            refs.extend((sid, end, street, False, .5) for street, sid, end in patch.attachments)
            for sid, end, street, half, inner in refs:
                if end not in (0, 1):
                    raise BundleJunctionError('SECTION_END', 'An incident Section end must be 0 or 1')
                role = (street, half, inner)
                if sid not in records:
                    records[sid] = (role, [None, None])
                previous, owners = records[sid]
                if previous != role:
                    raise BundleJunctionError('SECTION_PROFILE_CONFLICT', 'A shared Section has incompatible owner or profile roles')
                if owners[end] is not None:
                    raise BundleJunctionError('SECTION_END_OWNER', 'A Section end has more than one junction owner')
                owners[end] = jid
                incident.setdefault(sid, set()).add(jid)
        seams = {}; spans = []
        if include_member_spans:
            for member in topology.get('members', ()):
                street = member['street_id']
                inner = float((member['head_order'] == 0) == member['forward'])
                role = (street, True, inner)
                for span in member['spans']:
                    boundaries = (span['start'], span['end'])
                    ends = {b['junction_id'] for b in boundaries if b.get('junction_id') in by_junction}
                    if not ends or not all(b.get('junction_id') in by_junction or
                        (b.get('junction_id') is None and b.get('bundle_end') in ('HEAD','TAIL')
                         and member.get('complete',False)) for b in boundaries):
                        continue
                    ids = [r['start']['section_id'] for r in span['sections']]
                    if not ids or len(ids)!=len(set(ids)) or any(
                        r['start']['section_id']!=r['end']['section_id'] for r in span['sections']):
                        raise BundleJunctionError('SPAN_SECTION_RANGE', 'Each span range must address a unique original Section')
                    for sid in ids:
                        if sid in records and records[sid][0]!=role:
                            raise BundleJunctionError('SECTION_PROFILE_CONFLICT', 'A span Section has incompatible owner or profile roles')
                        records.setdefault(sid,(role,[None,None]));incident.setdefault(sid,set()).update(ends)
                    authored = ids if member['forward'] else ids[::-1]
                    for before,after in zip(authored,authored[1:]):
                        if records[before][1][1] is not None or records[after][1][0] is not None:
                            raise BundleJunctionError('SEAM_AT_JUNCTION', 'A Section seam cannot replace an incident junction')
                        for sid,end,partner in ((before,1,after),(after,0,before)):
                            key=(sid,end)
                            if key in seams and seams[key]!=partner:
                                raise BundleJunctionError('SECTION_SEAM_OWNER', 'A Section end has more than one seam partner')
                            seams[key]=partner
                    spans.append((street,span['ordinal'],tuple(ids)))
        adjacency = {jid: set() for jid in by_junction}
        for owners in incident.values():
            for jid in owners:
                adjacency[jid].update(owners - {jid})
        remaining = set(by_junction)
        networks = []
        while remaining:
            pending = [min(remaining)]; connected = set()
            while pending:
                jid = pending.pop()
                if jid in connected:
                    continue
                connected.add(jid); pending.extend(adjacency[jid] - connected)
            remaining.difference_update(connected)
            sections = []
            for sid in sorted(records):
                if not incident[sid] & connected:
                    continue
                (street, half, inner), owners = records[sid]
                sections.append(IncidentSection(sid, street, half, inner, *owners))
            included={s.section_id for s in sections}
            seam_records=tuple(IncidentSeam(records[sid][0][0],sid,partner)
                for (sid,end),partner in sorted(seams.items()) if end==1 and sid in included)
            span_records=tuple((street,ordinal) for street,ordinal,ids in spans if set(ids)<=included)
            networks.append(BundleIncidentNetwork(topology['bundle_id'],
                tuple(by_junction[jid] for jid in sorted(connected)), tuple(sections),seam_records,span_records))
        return tuple(networks)


class BoundedBundlePlanner:
    """Plan a complete two-member component from existing semantic incidences.

    No coordinates, nearest-neighbour grouping or invented endpoint IDs are
    used. With include_branch_junctions, complete ordinary branch Streets may
    contain prepared intermediate junctions, including recursively reached
    branches. A branch endpoint entering an unowned major junction still needs
    broader component closure. Partial extents and separator openings remain
    explicitly unsupported here.
    """
    def __init__(self, categories=None, *, allow_free_ends=False, include_branch_junctions=False):
        # Without an explicit feature selection, do not guess that an interior
        # Street may be omitted. Normal import supplies StreetDesign.categories.
        self.categories = None if categories is None else frozenset(categories)
        self.allow_free_ends = allow_free_ends
        self.include_branch_junctions = include_branch_junctions

    def _ordinary_junction_refs(self, snapshot, street_id, item):
        """Verify the complete declared minor-junction incidence, without geometry inference."""
        jid=item['id'];junction=snapshot.junctions.get(jid)
        if junction is None:
            reason=snapshot.intersections.get(jid,{}).get('reason','No prepared junction contract')
            raise BundleJunctionError('BRANCH_JUNCTION_INPUT',f'Street {street_id}, junction {jid}: {reason}')
        source=[s for s in snapshot.streets[street_id]['items'] if s['kind']=='Section']
        through={tuple(ref) for ref in item.get('through',())}
        pairs=[((a['id'],1),(b['id'],0)) for a,b in zip(source,source[1:])]
        if len(item.get('through',()))!=2 or not any(through==set(pair) for pair in pairs):
            raise BundleJunctionError('BRANCH_THROUGH_INCIDENCE',
                f'Junction {jid} needs the two consecutive original through Section ends')
        expected=[(street_id,sid,end) for sid,end in through]
        expected.extend((a['street_id'],a['section_id'],a['end']) for a in item.get('attachments',())
            if self.categories is None or a.get('category') in self.categories)
        actual=[]
        for c in junction.connectors:
            profile=junction.profiles[c.profile_index]
            if len(c.members)!=1 or c.bundle_source_id>=0 or profile.separator or getattr(profile,'half_profile',0):
                raise BundleJunctionError('BRANCH_JUNCTION_PROFILE',
                    f'Ordinary branch junction {jid} needs single, full-profile approaches')
            m=c.members[0];street=junction.streets[m.street_index]
            s=street.sections[0 if m.leaving else -1]
            actual.append((street.source_id,s.source_id,int(not m.leaving)))
        if len(expected)<3 or len(set(expected))!=len(expected) or len(set(actual))!=len(actual) or set(actual)!=set(expected):
            raise BundleJunctionError('BRANCH_JUNCTION_INCIDENCE',
                f'Junction {jid} must retain both through ends and every selected original attachment')
        return tuple(actual)

    def plan(self, snapshot, bundle_id):
        return self._plan(snapshot, bundle_id, frozenset((bundle_id,)), frozenset())

    def _plan(self, snapshot, bundle_id, component_bundles, component_junctions):
        topology = snapshot.bundles[bundle_id].get('topology', {})
        if topology.get('status') != 'MAPPED':
            raise BundleJunctionError('UNMAPPED_EXTENTS', 'Bundle extents must be mapped first')
        members = topology['members']
        if len(members) != 2:
            raise BundleJunctionError('MEMBER_COUNT', 'Bounded assembly requires two outer members')
        if {m['head_order'] for m in members} != {0, 1} or any(m['tail_order'] != m['head_order'] for m in members):
            raise BundleJunctionError('MEMBER_ORDER', 'The two members must retain their lateral order')
        if not all(m['complete'] for m in members):
            raise BundleJunctionError('PARTIAL_EXTENTS', 'Partial member extents need a separate range provider')
        deferred = []
        for sid in topology.get('interior_streets', ()):
            street = snapshot.streets.get(sid)
            sections = [s for s in street['items'] if s['kind']=='Section'] if street else []
            if (self.categories is not None and sections and all(s.get('category') and
                    s['category'] not in self.categories for s in sections)):
                deferred.append(sid)
            else:
                raise BundleJunctionError('INTERIOR_STREETS', 'Interior carriageway Streets need a separate opening provider')
        member_by_id = {m['street_id']: m for m in members}
        if len(member_by_id) != 2:
            raise BundleJunctionError('DUPLICATE_MEMBER', 'The member Street identities must be distinct')
        side_planner = BundleJunctionPlanner()
        patches = tuple(side_planner.plan(topology, j) for j in topology.get('junctions', ()))
        if len({p.junction_id for p in patches}) != len(patches):
            raise BundleJunctionError('DUPLICATE_JUNCTION', 'Each side junction needs one owner')

        sections_by_street = {}
        branch_items = {}  # Junction ID -> (original through Street ID, source item).
        rows = {}  # Section ID -> [Street ID, half-profile, inner_u, end owners].
        def include_street(sid):
            if sid in sections_by_street:
                return sections_by_street[sid]
            street = snapshot.streets.get(sid)
            if street is None or street['status'] != 'SUPPORTED':
                raise BundleJunctionError('STREET_SCOPE', f'Street {sid} is unavailable or unsupported')
            is_member = sid in member_by_id
            if not is_member and street.get('bundle_id') is not None:
                raise BundleJunctionError('OTHER_BUNDLE', 'An attachment belongs to another Bundle')
            for item in street['items']:
                if item['kind'] == 'Section':
                    continue
                allowed = ('COMPATIBLE_SEAM', 'WIDTH_TRANSITION', 'DEFERRED_ATTACHMENT')
                if is_member or self.include_branch_junctions:
                    allowed += ('NATIVE_JUNCTION',)
                if item.get('surface_treatment') not in allowed:
                    raise BundleJunctionError('BRANCH_TOPOLOGY',
                        f'Street {sid}, {item["kind"]} {item["id"]} '
                        f'({item.get("surface_treatment", "unclassified")}) requires another junction or feature provider')
                if not is_member and item.get('surface_treatment')=='NATIVE_JUNCTION':
                    jid=item['id']
                    if jid in branch_items and branch_items[jid][0]!=sid:
                        raise BundleJunctionError('BRANCH_JUNCTION_OWNER','An intermediate junction has two through-Street owners')
                    branch_items[jid]=(sid,item)
            sections = [s for s in street['items'] if s['kind'] == 'Section']
            if not sections:
                raise BundleJunctionError('EMPTY_STREET', 'A planned Street needs an original Section')
            member = member_by_id.get(sid)
            inner = float((member['head_order'] == 0) == member['forward']) if member else .5
            for section in sections:
                if section['id'] in rows:
                    raise BundleJunctionError('DUPLICATE_SECTION', 'An original Section has multiple Street owners')
                rows[section['id']] = [sid, is_member, inner, [None, None]]
            sections_by_street[sid] = sections
            return sections

        for member in members:
            source = include_street(member['street_id'])
            actual = [r['start']['section_id'] for s in member['spans'] for r in s['sections']]
            expected = [s['id'] for s in source]
            if not member['forward']:
                expected.reverse()
            if (actual != expected or any(r['start']['section_id'] != r['end']['section_id']
                    for s in member['spans'] for r in s['sections'])):
                raise BundleJunctionError('SPAN_SECTION_RANGE', 'Member spans must cover each complete original Section exactly once')

        def own(street_id, section_id, end, junction_id):
            include_street(street_id)
            if section_id not in rows or rows[section_id][0] != street_id or end not in (0, 1):
                raise BundleJunctionError('SECTION_REFERENCE', 'A junction references the wrong Street or Section end')
            if rows[section_id][3][end] is not None:
                raise BundleJunctionError('SECTION_END_OWNER', 'A Section end has more than one junction owner')
            rows[section_id][3][end] = junction_id

        terminals = {}
        for junction in snapshot.junctions.values():
            paired = [c for c in junction.connectors if c.bundle_source_id == bundle_id]
            if not paired:
                continue
            if len(paired) != 1 or len(paired[0].members) != 2:
                raise BundleJunctionError('PAIRED_TERMINAL', 'A terminal requires exactly one paired Bundle connector')
            connector = paired[0]
            sides, ids = set(), []
            for ref in connector.members:
                street = junction.streets[ref.street_index]
                member = member_by_id.get(street.source_id)
                if member is None:
                    raise BundleJunctionError('MEMBER_REFERENCE', 'A terminal references a different Bundle member')
                ids.append(street.source_id)
                sides.add(0 if ref.leaving == member['forward'] else 1)
            if len(sides) != 1 or set(ids) != set(member_by_id):
                raise BundleJunctionError('TERMINAL_INCIDENCE', 'Both members must enter the same declared Bundle end')
            side = sides.pop()
            if side in terminals:
                raise BundleJunctionError('TERMINAL_OWNER', 'A Bundle end belongs to multiple junctions')
            terminals[side] = junction.source_id
            for c in junction.connectors:
                if c is not connector and (len(c.members) != 1 or c.bundle_source_id >= 0):
                    if c.bundle_source_id in component_bundles and len(c.members)==2:
                        # Its own member planner supplies these original
                        # Sections. The complete component owns this patch once.
                        continue
                    raise BundleJunctionError('OTHER_BUNDLE', 'Multiple Bundles at one terminal need a broader network provider')
                for ref in c.members:
                    street = junction.streets[ref.street_index]
                    source = include_street(street.source_id)
                    section = street.sections[0 if ref.leaving else -1]
                    if section.source_id != source[0 if ref.leaving else -1]['id']:
                        raise BundleJunctionError('TERMINAL_SECTION', 'A terminal must address the complete Street endpoint')
                    own(street.source_id, section.source_id, int(not ref.leaving), junction.source_id)
        free_ends=[]
        if self.allow_free_ends:
            for side in sorted({0,1}-terminals.keys()):
                refs=[]
                for member in sorted(members,key=lambda m:m['head_order']):
                    sid=member['street_id'];end=side if member['forward'] else 1-side
                    declared=snapshot.streets[sid].get('endpoints',(None,None))[end]
                    event=snapshot.intersections.get(declared,{})
                    if declared is not None and not (event.get('status')=='ENDPOINT' and event.get('order')==1):
                        raise BundleJunctionError('MISSING_TERMINAL',
                            f'Bundle {bundle_id}, Street {sid}, authored end {end} still references '
                            f'intersection {declared}; an unresolved endpoint cannot be treated as free')
                    section=sections_by_street[sid][0 if end==0 else -1]
                    refs.append((section['id'],end))
                free_ends.append(FreeBundleEnd(bundle_id,side,tuple(refs)))
        if set(terminals)|{e.side for e in free_ends} != {0,1} or (len(terminals)==2 and terminals[0]==terminals[1]):
            raise BundleJunctionError('MISSING_TERMINAL',
                'Both Bundle ends need distinct paired junctions or explicitly permitted free-end records')
        if set(terminals.values()) & {p.junction_id for p in patches}:
            raise BundleJunctionError('DUPLICATE_JUNCTION', 'A junction cannot be both terminal and internal')
        for patch in patches:
            for ref in patch.through:
                own(patch.street_id, ref.section_id, ref.end, patch.junction_id)
            for street, section, end in patch.attachments:
                source = include_street(street)
                if section != source[0 if end == 0 else -1]['id']:
                    raise BundleJunctionError('ATTACHMENT_SECTION', 'An attachment must address a complete Street endpoint')
                own(street, section, end, patch.junction_id)

        ordinary=set()
        reserved=set(terminals.values())|{p.junction_id for p in patches}
        # include_street may discover another intermediate junction on a newly
        # included complete attachment. Process the closure, never truncate it.
        while set(branch_items)-ordinary:
            jid=min(set(branch_items)-ordinary)
            if jid in reserved:
                raise BundleJunctionError('JUNCTION_OWNER','A junction cannot be both an ordinary branch patch and a Bundle patch')
            street,item=branch_items[jid]
            for owner,sid,end in self._ordinary_junction_refs(snapshot,street,item):
                own(owner,sid,end,jid)
            ordinary.add(jid)

        seams = []
        for sid, source in sorted(sections_by_street.items()):
            for before, after in zip(source, source[1:]):
                a, b = before['id'], after['id']
                owners = rows[a][3][1], rows[b][3][0]
                if owners == (None, None):
                    seams.append(IncidentSeam(sid, a, b))
                elif owners[0] is None or owners[0] != owners[1]:
                    raise BundleJunctionError('CHAIN_INCIDENCE', 'Consecutive Sections need either a shared junction or a seam')
            endpoints = snapshot.streets[sid].get('endpoints', (None, None))
            for end, section in ((0, source[0]), (1, source[-1])):
                owner = rows[section['id']][3][end]
                external = endpoints[end]
                if external is not None and external != owner:
                    event = snapshot.intersections.get(external, {})
                    if event.get('status') != 'ENDPOINT' and external not in component_junctions:
                        raise BundleJunctionError('EXTERNAL_JUNCTION',
                            f'Street {sid}, authored end {end} enters unowned junction {external}; include it in a broader network')
                if sid in member_by_id and owner is None and not any(
                        (section['id'],end) in free.member_ends for free in free_ends):
                    raise BundleJunctionError('UNBOUND_MEMBER', 'Each complete member end needs a paired terminal or explicit free-end record')
        return BoundedBundleNetwork(bundle_id, (terminals.get(0), terminals.get(1)),
            tuple(sorted(patches, key=lambda p: p.junction_id)),
            tuple(IncidentSection(section, street, half, inner, *owners)
                for section, (street, half, inner, owners) in sorted(rows.items())),
            tuple(seams), tuple((m['street_id'], s['ordinal']) for m in
                sorted(members, key=lambda m: m['street_id']) for s in m['spans']),tuple(sorted(deferred)),tuple(free_ends),
            tuple(sorted(ordinary)))

    def review(self, snapshot, bundle_id):
        try:
            return self.plan(snapshot, bundle_id).report()
        except BundleJunctionError as exc:
            return dict(status='UNSUPPORTED', code=exc.code, reason=str(exc))


@dataclass(frozen=True)
class CoupledBundleNetwork:
    bundle_ids: tuple[int, ...]
    terminal_junctions: tuple[int, ...]
    junctions: tuple[BundleSideJunction, ...]
    sections: tuple[IncidentSection, ...]
    seams: tuple[IncidentSeam, ...]
    member_spans: tuple[tuple[int, int, int], ...]  # Bundle, Street, ordinal.
    section_bundles: tuple[tuple[int, int], ...]
    deferred_streets: tuple[int, ...]
    free_ends: tuple[FreeBundleEnd, ...] = ()
    ordinary_junctions: tuple[int, ...] = ()

    def report(self):
        return dict(status='PLANNED_COUPLED_NETWORK', **asdict(self))


class CoupledBundlePlanner:
    """Connected components of explicit paired Bundle terminal incidences.

    Complete terminal patches and shared ordinary branch Sections are owned
    once. Plans never infer adjacency from nearby coordinates, silently drop a
    rejected member, or promote a partial Bundle to a complete one.
    """
    def __init__(self, categories=None, *, allow_free_ends=False, include_branch_junctions=False):
        self.member_planner=BoundedBundlePlanner(categories,allow_free_ends=allow_free_ends,
            include_branch_junctions=include_branch_junctions)

    @staticmethod
    def component(snapshot, seed):
        if seed not in snapshot.bundles:
            raise BundleJunctionError('MISSING_BUNDLE', 'A component references an absent Bundle')
        incidence={bid:set() for bid in snapshot.bundles}
        groups=[{c.bundle_source_id for c in j.connectors if c.bundle_source_id>=0}
                for j in snapshot.junctions.values()]
        # Retain connectivity even when a logical junction has no qualified
        # radial contract. Several such records can own incidences at one
        # original physical site; choosing separate Bundle components would
        # duplicate or disconnect that site.
        origins=getattr(snapshot,'junction_origins',{})
        logical={row['junction']['junction_id']:{c['item_id'] for c in row['junction']['ring']
                 if c['item_kind']=='Bundle'} for row in origins.get('logical',())}
        groups.extend(logical.values())
        for component in origins.get('components',()):
            groups.append({bid for jid in component['logical_junction_ids'] for bid in logical.get(jid,())})
        for ids in groups:
            for bid in ids:
                if bid not in incidence:
                    raise BundleJunctionError('MISSING_BUNDLE', 'A terminal references an absent Bundle')
                incidence[bid].update(ids-{bid})
        found=set();pending=[seed]
        while pending:
            bid=pending.pop()
            if bid in found:continue
            found.add(bid);pending.extend(incidence[bid]-found)
        return tuple(sorted(found))

    def plan(self, snapshot, seed):
        bids=self.component(snapshot,seed)
        origins=getattr(snapshot,'junction_origins',{})
        logical={row['junction']['junction_id']:{c['item_id'] for c in row['junction']['ring']
                 if c['item_kind']=='Bundle'} for row in origins.get('logical',())}
        shared={pid for component in origins.get('components',())
                if any(set(bids)&logical.get(jid,set()) for jid in component['logical_junction_ids'])
                for pid in component['shared_physical_junction_ids']}
        if shared:
            raise BundleJunctionError('SHARED_PHYSICAL_JUNCTION',
                'Logical Bundle junctions require one shared physical construction for sites '+
                ', '.join(map(str,sorted(shared))))
        terminal_ids={j.source_id for j in snapshot.junctions.values()
            if any(c.bundle_source_id in bids for c in j.connectors)}
        parts=[self.member_planner._plan(snapshot,bid,frozenset(bids),terminal_ids) for bid in bids]
        sections={};owners={};patches={};seams=set();spans=[];deferred=set();ordinary=set()
        for part in parts:
            deferred.update(part.deferred_streets)
            ordinary.update(part.ordinary_junctions)
            seams.update(part.seams)
            spans.extend((part.bundle_id,*pair) for pair in part.member_spans)
            for patch in part.junctions:
                if patch.junction_id in patches and patches[patch.junction_id]!=patch:
                    raise BundleJunctionError('JUNCTION_OWNER', 'A side patch has incompatible Bundle owners')
                patches[patch.junction_id]=patch
            for section in part.sections:
                old=sections.get(section.section_id)
                if old:
                    if (old.street_id,old.half_profile,old.inner_u)!=(section.street_id,section.half_profile,section.inner_u):
                        raise BundleJunctionError('SECTION_PROFILE_CONFLICT', 'A shared Section has conflicting ownership or profiles')
                    ends=[]
                    for a,b in ((old.start_junction,section.start_junction),(old.end_junction,section.end_junction)):
                        if a is not None and b is not None and a!=b:
                            raise BundleJunctionError('SECTION_END_OWNER', 'A shared Section end has multiple junction owners')
                        ends.append(b if a is None else a)
                    section=IncidentSection(section.section_id,section.street_id,section.half_profile,section.inner_u,*ends)
                sections[section.section_id]=section
                if section.half_profile:
                    previous=owners.setdefault(section.section_id,part.bundle_id)
                    if previous!=part.bundle_id:
                        raise BundleJunctionError('SECTION_BUNDLE_OWNER', 'A complete Section belongs to multiple Bundles')
        if terminal_ids & patches.keys():
            raise BundleJunctionError('JUNCTION_OWNER', 'A patch cannot be both a paired terminal and a side junction')
        if ordinary & (terminal_ids|patches.keys()):
            raise BundleJunctionError('JUNCTION_OWNER','An ordinary branch patch has a conflicting Bundle patch role')
        # Ordinary Streets may connect two included terminal patches. Their
        # two endpoint owners must be supplied after composing all members.
        for sid in {s.street_id for s in sections.values()}:
            street=snapshot.streets[sid]
            original=[s for s in street['items'] if s['kind']=='Section']
            for end,source in ((0,original[0]),(1,original[-1])):
                section=sections[source['id']]
                owner=section.start_junction if end==0 else section.end_junction
                declared=street.get('endpoints',(None,None))[end]
                if declared is not None and declared!=owner and snapshot.intersections.get(declared,{}).get('status')!='ENDPOINT':
                    raise BundleJunctionError('EXTERNAL_JUNCTION', 'A component still has a branch entering an unowned junction')
        return CoupledBundleNetwork(bids,tuple(sorted(terminal_ids)),
            tuple(patches[j] for j in sorted(patches)),tuple(sections[s] for s in sorted(sections)),
            tuple(sorted(seams,key=lambda s:(s.street_id,s.before_section,s.after_section))),
            tuple(sorted(spans)),tuple(sorted(owners.items())),tuple(sorted(deferred)),
            tuple(free for part in parts for free in part.free_ends),tuple(sorted(ordinary)))

    def review(self, snapshot):
        remaining=set(snapshot.bundles);result=[]
        while remaining:
            seed=min(remaining)
            bids=(seed,)
            try:
                bids=self.component(snapshot,seed)
                remaining.difference_update(bids)
                result.append(self.plan(snapshot,seed).report())
            except BundleJunctionError as exc:
                remaining.discard(seed)
                result.append(dict(status='UNSUPPORTED',bundle_ids=bids,
                    code=exc.code,reason=str(exc)))
        return result
