"""Map vehicle attachments inside a Blosm Street to the ordinary GN kernel.

Semantic Street ownership survives an internal junction. Its two through
Section ends become separate connectors and the kernel owns their join.
Non-carriageway attachments remain explicitly deferred feature records.
"""
import math
from .adapter import PreparedNetworkAdapter
from ..gn.kernel import inputs
from ..way_properties import wayCategoryProps


class MinorJunctionAdapter:
    def __init__(self, design):
        self.design = design

    def read(self, snapshot):
        sections = {i['id']: (street, i) for street in snapshot.streets.values()
                    for i in street['items'] if i['kind']=='Section'}
        next_connector = 1+max((c.source_id for j in snapshot.junctions.values()
                               for c in j.connectors), default=-1)
        # Avoid collisions with connectors belonging to deferred minor features.
        next_connector = max(next_connector, 1+max((a['connector_id']
            for s in snapshot.streets.values() for i in s['items']
            for a in i.get('attachments', ())), default=-1))
        for street in snapshot.streets.values():
            for item in street['items']:
                treatment = item.get('surface_treatment')
                if treatment not in ('NATIVE_JUNCTION', 'DEFERRED_ATTACHMENT'):
                    continue
                jid = item['id']
                row = dict(id=jid, location=item['location'], order=2, status='DEFERRED_ATTACHMENT', minor=True)
                snapshot.intersections[jid] = row
                if treatment == 'DEFERRED_ATTACHMENT':
                    continue
                try:
                    if item.get('attachment_error'):
                        raise ValueError(item['attachment_error'])
                    refs = list(item['through'])
                    refs.extend((a['section_id'], a['end']) for a in item['attachments']
                                if a['category'] in self.design.categories)
                    arms = []
                    for sid, end in refs:
                        owner, section = sections[sid]
                        if owner['status'] != 'SUPPORTED':
                            raise ValueError(f"Attached Street {owner['id']}: {owner['reason']}")
                        points = section['points']
                        a, b = (points[0],points[1]) if end==0 else (points[-1],points[-2])
                        arms.append((math.atan2(b[1]-a[1],b[0]-a[0]) % math.tau, owner, section, end))
                    arms.sort(key=lambda arm: arm[0])
                    sources, connectors, profiles = [], [], []
                    for index, (_angle, owner, section, end) in enumerate(arms):
                        points = section['points'][:2] if end==0 else section['points'][-2:]
                        sources.append(inputs.StreetInput(owner['id'], (inputs.SectionInput(section['id'],
                            points, section['width'], section['offset'] or 0., section['category']),)))
                        profiles.append(inputs.ProfileInput(self.design.profile_rise,
                            inputs.SMOOTH_CROWN if self.design.profile_rise else inputs.FLAT))
                        connectors.append(inputs.ConnectorInput(next_connector,
                            (inputs.MemberInput(index, end==0),), index,
                            radius=wayCategoryProps[section['category']]['radius']))
                        next_connector += 1
                    junction = inputs.PreparedJunction(jid, tuple(item['location']), tuple(sources),
                        tuple(connectors), tuple(profiles), 0., self.design.cell_size)
                    row.update(order=len(arms), status='UNSUPPORTED', streets=tuple(s.source_id for s in sources))
                    inputs.validate_initial_scope(junction)
                    snapshot.junctions[jid] = junction
                    row['status'] = 'SUPPORTED'
                except (ValueError, KeyError) as exc:
                    row.update(status='UNSUPPORTED', reason=str(exc))
                    snapshot.issue('MinorIntersection', jid, 'INTERSECTION_SCOPE', str(exc))
