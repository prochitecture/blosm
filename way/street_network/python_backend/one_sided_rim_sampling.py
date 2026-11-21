"""Opt-in bounded chord resampling for a one-sided preparation's real holes.

This helper has no production hook. All incoming profile samples and the
exterior stay exact. Original hole polylines, source spans and source curves
remain available; only their derived sampling may project by at most 1 um.
"""
from copy import copy
from dataclasses import replace
import numpy as np
from .compounds import RimSegment, RimSpan, rim_label
from .regions import BoundaryLoop, Region, RegionUnion
from .source_rim_normalization import ResampledSourceCurb


def normalize_holes(group, *, maximum_projection=1e-6):
    from .updates import drain_steps
    return drain_steps(normalize_holes_steps(group, maximum_projection=maximum_projection))


def normalize_holes_steps(group, *, maximum_projection=1e-6):
    from .one_sided_paired_junctions import certify_steps
    if not np.isfinite(maximum_projection) or not 0. < maximum_projection <= 1e-6:
        raise ValueError('One-sided hole sampling allows at most one micrometre projection')
    yield from certify_steps(group)
    layout = group.physical_layout
    if getattr(layout, 'source_hole_normalization', None):
        if any(row['additional_error_allowance'] > maximum_projection for row in layout.source_hole_normalization):
            raise ValueError('A normalized physical hole requires its original projection allowance')
        return group
    if not any(loop.hole for loop in group.loops):
        return group
    physical = copy(layout)
    curves, reports, loops = list(layout.boundaries), [], []
    sources = {source.region.owner: source for source in group.footprints}
    for loop in group.loops:
        if not loop.hole:
            loops.append(loop)
            continue
        labels = [rim_label(span, sources) for span in loop.spans]
        centered = loop.xy-loop.xy.mean(axis=0)
        _, axes = np.linalg.eigh(centered.T@centered)
        projected = centered@axes[:, -1]
        # The original consolidation's half-by-index point is not a physical
        # tip. Unequal bank sampling can put it arbitrarily close to one tip,
        # pinning a microscopic derived interval. Keep the two actual extreme
        # vertices instead; all earlier samples survive as source provenance.
        starts = sorted((int(np.argmin(projected)), int(np.argmax(projected))))
        if starts[0] == starts[1]:
            raise ValueError('A real one-sided hole has no resolved finite bank endpoints')
        points, spans = [], []
        for first, last in zip(starts, (*starts[1:], starts[0]+len(loop.xy))):
            indices = np.arange(first, last+1) % len(loop.xy)
            old_spans = tuple(loop.spans[i % len(loop.xy)] for i in range(first, last))
            provenance = {(span.owner, labels[i % len(labels)].index):
                          sources[span.owner].layout.boundaries[labels[i % len(labels)].index]
                          for i, span in zip(range(first, last), old_spans)}
            curve = ResampledSourceCurb(loop.xy[indices], old_spans,
                [(owner, index, source) for (owner, index), source in sorted(provenance.items())], maximum_projection)
            curve.on_separator_rim = True
            index = len(curves)
            curves.append(curve)
            reports.append(curve.report())
            values = curve.initial_parameters
            points.extend(curve(values)[:-1])
            spans.extend(RimSpan(group.members[0], RimSegment('curb', index, float(a), float(b)))
                         for a, b in zip(values, values[1:]))
        boundary = BoundaryLoop(np.asarray(points), tuple(spans))
        if not boundary.hole:
            raise ValueError('One-sided hole resampling reversed its physical orientation')
        RegionUnion().build((Region(group.members[0], boundary.xy[::-1]),))
        loops.append(boundary)
    physical.boundaries = tuple(curves)
    physical.source_hole_normalization = tuple(reports)
    footprints = tuple(replace(source, layout=physical) if source.layout is layout else source
                       for source in group.footprints)
    result = replace(group, physical_layout=physical, footprints=footprints, loops=tuple(loops))
    yield from certify_steps(result)
    return result
