"""Adapt prepared lane-change items to explicit editable dimensional tapers.

The prepared full Section centerlines and constant dimensions remain intact.
Legacy trim polygons are not baked into output. GN measures the available
distance and calculates the taper from live widths and signed offsets.
"""
import math


class PreparedTransitionAdapter:
    kinds = ('SymLane', 'SideLane')

    def __init__(self, *, minimum_half_length=1.5, maximum_edge_slope=.3):
        if (not math.isfinite(minimum_half_length) or minimum_half_length <= 0
                or not math.isfinite(maximum_edge_slope) or not 0 < maximum_edge_slope <= 1):
            raise ValueError('Invalid dimensional transition policy')
        self.minimum_half_length = minimum_half_length
        self.maximum_edge_slope = maximum_edge_slope

    def read(self, entries):
        for index, item in enumerate(entries):
            deferred = item['kind'] == 'Intersection' and item.get('surface_treatment') == 'DEFERRED_ATTACHMENT'
            if (item['kind'] not in self.kinds and not deferred) or item.get('surface_treatment') == 'COMPATIBLE_SEAM':
                continue
            if not 0 < index < len(entries)-1:
                item['transition_reason'] = 'Lane transition has no two adjacent Sections'
                continue
            before, after = entries[index-1:index+2:2]
            if before['kind'] != 'Section' or after['kind'] != 'Section':
                item['transition_reason'] = 'Lane transition requires explicit adjacent Section ownership'
                continue
            if item.get('points'):
                item['transition_reason'] = 'A transition with its own prepared centerline needs an extent adapter'
                continue
            if not before['points'] or not after['points'] or math.dist(before['points'][-1],after['points'][0]) > 1e-6:
                item['transition_reason'] = 'Lane transition Section endpoints are disconnected'
                continue
            if deferred:
                through = tuple(tuple(end) for end in item.get('through', ()))
                if len(through) != 2 or set(through) != {(before['id'], 1), (after['id'], 0)}:
                    item['transition_reason'] = 'Deferred dimensional seam has different through Section owners'
                    continue
                location = item.get('location', ())
                if (not item.get('attachments') or len(location) < 2
                        or not all(math.isfinite(v) for v in location)
                        or math.dist(location[:2], before['points'][-1][:2]) > 1e-6):
                    item['transition_reason'] = 'Deferred dimensional seam needs its original attachment location'
                    continue
                # Keep an existing declaration dormant when dimensions become
                # equal. Otherwise equal-dimension attachments need no taper.
                if (before['width'], before['offset']) == (after['width'], after['offset']) and not item.get('transition'):
                    continue
                item['transition_policy'] = 'THROUGH_DIMENSIONS_AT_DEFERRED_ATTACHMENT'
            else:
                item['surface_treatment'] = 'WIDTH_TRANSITION'
            item.pop('transition_reason', None)
            item['transition'] = dict(before_section=before['id'],after_section=after['id'],
                minimum_half_length=self.minimum_half_length,maximum_edge_slope=self.maximum_edge_slope,
                source_kind=item['kind'],source_id=item['id'])
