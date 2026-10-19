"""Validate semantic Section ownership before junction/Bundle adaptation."""
from collections import defaultdict
import math


class PreparedOwnershipReview:
    def apply(self, snapshot):
        owners = defaultdict(list)
        for street in snapshot.streets.values():
            for item in street['items']:
                if item['kind'] == 'Section':
                    owners[item['id']].append(street['id'])
        duplicates = []
        for section_id, streets in sorted(owners.items()):
            if len(streets) < 2:
                continue
            reason = f'Section {section_id} occurs more than once in prepared Street ownership: {streets}'
            duplicates.append(dict(section_id=section_id, street_ids=streets))
            snapshot.issue('Section', section_id, 'DUPLICATE_SECTION_OWNERSHIP', reason)
            for sid in streets:
                snapshot.streets[sid].update(status='UNSUPPORTED', reason=reason)
        snapshot.ownership_audit = dict(section_count=len(owners),
            section_occurrence_count=sum(map(len, owners.values())), duplicate_sections=duplicates)
        self._attachments(snapshot, owners)

    @staticmethod
    def _attachments(snapshot, owners):
        """Resolve stale Street references through already identified Sections.

        Concatenation can leave a minor connector pointing to a superseded
        Street object. Its Section/end reference is the stable incidence; use
        only a unique owner, and verify the declared junction location. Never
        select a replacement by the nearest Street or create a new connection.
        """
        sections={i['id']:i for s in snapshot.streets.values()
                  for i in s['items'] if i['kind']=='Section' and len(owners[i['id']])==1}
        remapped,unresolved=[],[]
        for street in snapshot.streets.values():
            for item in street['items']:
                item.pop('attachment_error',None)
                for attachment in item.get('attachments', ()):
                    sid,end=attachment['section_id'],attachment['end']
                    section=sections.get(sid)
                    reason=None
                    if section is None or end not in (0,1):
                        reason='Attachment has no uniquely owned Section endpoint'
                    elif not section['points'] or not math.dist(
                            section['points'][0 if end==0 else -1][:2], item['location'][:2])<=1e-5:
                        reason='Attachment Section endpoint does not meet its declared junction'
                    if reason:
                        record=dict(junction_id=item['id'],connector_id=attachment['connector_id'],
                                    section_id=sid,end=end,reason=reason)
                        unresolved.append(record)
                        item['attachment_error']=reason
                        snapshot.issue('MinorIntersection',item['id'],'ATTACHMENT_OWNERSHIP',reason)
                        continue
                    actual=owners[sid][0]
                    if actual!=attachment['street_id']:
                        remapped.append(dict(junction_id=item['id'],connector_id=attachment['connector_id'],
                            section_id=sid,end=end,prepared_street_id=attachment['street_id'],street_id=actual))
                        attachment['street_id']=actual
        snapshot.ownership_audit.update(remapped_attachments=remapped,unresolved_attachments=unresolved)
