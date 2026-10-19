"""Review complete prepared Bundle members after minor-attachment adaptation.

Street and Section ownership remain distinct from the derived continuous sweep.
Future extent/transition providers can replace these bounded eligibility rules
without changing the renderer or the persistent source schema.
"""
import math


class PreparedBundleReview:
    def member_reason(self, street):
        if street['status'] != 'SUPPORTED':
            return street['reason'] or 'Bundle member is outside the current surface domain'
        sections = [item for item in street['items'] if item['kind'] == 'Section']
        if not sections:
            return 'Bundle member contains no Section'
        for item in street['items']:
            if item['kind'] == 'Section':
                continue
            if item.get('surface_treatment') == 'NATIVE_JUNCTION':
                return 'Internal Bundle member spans require a native junction surface patch'
            if item.get('surface_treatment') not in ('COMPATIBLE_SEAM', 'DEFERRED_ATTACHMENT', 'WIDTH_TRANSITION'):
                return 'Internal Bundle junction or transition requires explicit member spans'
        transitions={(i['transition']['before_section'],i['transition']['after_section'])
                     for i in street['items'] if i.get('transition')}
        for left, right in zip(sections, sections[1:]):
            if ((abs(left['width']-right['width']) > 1e-6
                    or abs((left['offset'] or 0.)-(right['offset'] or 0.)) > 1e-6)
                    and (left['id'],right['id']) not in transitions):
                return 'Changing Bundle member dimensions require an explicit transition provider'
            if math.dist(left['points'][-1], right['points'][0]) > 1e-5:
                return 'Bundle member Sections are not a continuous directed chain'
        return None

    def apply(self, snapshot):
        for bundle in snapshot.bundles.values():
            if bundle['status'] != 'SUPPORTED':
                continue
            topology=bundle.get('topology', {})
            reason = (topology.get('reason') if topology.get('status')=='UNRESOLVED' else
                      next((reason for sid in bundle['head']
                            if (reason := self.member_reason(snapshot.streets[sid]))), None))
            if reason:
                bundle.update(status='UNSUPPORTED', reason=reason)
                snapshot.issue('Bundle', bundle['id'], 'BUNDLE_SCOPE', reason)
