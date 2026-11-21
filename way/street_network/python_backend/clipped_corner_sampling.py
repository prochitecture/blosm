"""Resample clipped circular fillets from their original arc definitions."""
from copy import copy
from dataclasses import replace
import numpy as np

from .compounds import rim_label
from .retained_arc_fillets import OriginalAnalyticalCircleCurve, _topology
from .sampling import AnalyticalBoundarySampler, CertifiedCurbSampler
from .source_circular_intervals import PreparedCircularInterval


class ClippedCircularCurb(OriginalAnalyticalCircleCurve):
    def certified_parameters(self, sampler, first, last):
        interval = self.original_prepared_curve.original_circle_interval
        angle = abs(self.corner.sweep) * max(np.diff(interval.circle_parameters))
        spacing = min(sampler.spacing, np.nextafter(self.corner.radius * angle, 0.))
        arc_sampler = CertifiedCurbSampler(min(sampler.tolerance, interval.maximum_chord_error),
                                           spacing, sampler.max_samples)
        low, high = self.circle_bounds
        breaks = [first, *(k for k in (low, high) if first < k < last), last]
        values = []
        for start, stop in zip(breaks, breaks[1:]):
            if low <= start and stop <= high:
                segment = super().certified_parameters(arc_sampler, start, stop)
            else:
                segment = self.original_prepared_curve.certified_parameters(sampler, start, stop)
            values.extend(segment[:-1])
        return np.r_[values, last]


def _mouth_edges(group):
    sources = {s.region.owner: s for s in group.footprints}
    result = {}
    for loop in group.loops:
        for i, span in enumerate(loop.spans):
            label = rim_label(span, sources)
            if label.kind == 'mouth':
                key = span.owner, label.index, label.first, label.last
                if key in result:
                    raise ValueError('Repeated physical profile edge in clipped-corner sampling')
                result[key] = loop.xy[i].tobytes(), loop.xy[(i+1) % len(loop.xy)].tobytes()
    return result


def resample(group, settings):
    sources = {s.region.owner: s for s in group.footprints}
    selected = set()
    for loop in group.loops:
        for i, span in enumerate(loop.spans):
            label = rim_label(span, sources)
            if label.kind != 'curb':
                continue
            curve = sources[span.owner].layout.boundaries[label.index]
            interval = getattr(curve, 'original_circle_interval', None)
            if not isinstance(interval, PreparedCircularInterval):
                continue
            low, high = interval.prepared_bounds
            if not low <= label.first < label.last <= high:
                continue
            start_clipped = loop.spans[i-1].owner != span.owner and low < label.first
            end_clipped = loop.spans[(i+1) % len(loop.spans)].owner != span.owner and label.last < high
            if not (start_clipped or end_clipped):
                continue
            first, last = interval.circle_parameter([label.first, label.last])
            step = max(np.diff(interval.circle_parameters))
            if 0. < last-first < step/8:
                selected.add((span.owner, label.index))
    if not selected:
        return None
    for owner, index in sorted(selected):
        source = sources[owner]
        layout = copy(source.layout)
        boundaries = list(layout.boundaries)
        original = boundaries[index]
        boundaries[index] = ClippedCircularCurb(owner, index, original,
                                                original.original_circle_interval.original_corner)
        layout.boundaries = tuple(boundaries)
        sources[owner] = replace(source, layout=layout)
    candidate = replace(group, footprints=tuple(sources[s.region.owner] for s in group.footprints))
    sampler = AnalyticalBoundarySampler(tolerance=settings.curve_tolerance, spacing=settings.spacing)
    candidate = sampler.build(candidate)
    _topology(group, candidate, settings.curve_tolerance)
    if _mouth_edges(candidate) != _mouth_edges(group):
        raise ValueError('Clipped-corner sampling changed an incoming profile edge')
    return candidate
