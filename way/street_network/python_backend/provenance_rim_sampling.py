"""Bounded derived-rim sampling with an explicit original-reference ledger.

This module changes sampling of prepared one-sided source unions only. Source
Paths, analytical local curbs, physical profiles and original hole tips remain
authoritative. Callers must try their existing chart first and opt in only
after a typed grading rejection; there is no discovery or publication hook.
"""
from copy import copy
from dataclasses import replace
import math
import numpy as np
from .compound_cycles import PolylineCurb
from .compounds import RimSegment, RimSpan, SampledFootprint, rim_label
from .model import Key
from .regions import BoundaryLoop, Region, RegionUnion
from .sampling import CertifiedCurbSampler


class BudgetedSourceCurb(PolylineCurb):
    def __init__(self, points, source_spans, source_curves, allowance):
        original = np.asarray(points, dtype=float)
        if (original.ndim != 2 or original.shape[1] != 2 or len(original) < 2
                or not np.isfinite(original).all() or not math.isfinite(allowance) or allowance <= 0.):
            raise ValueError('Budgeted source sampling needs a finite positive projection allowance')
        edge = np.diff(original, axis=0)
        lengths = np.linalg.norm(edge, axis=1)
        if np.any(lengths <= 1e-10):
            raise ValueError('Budgeted source sampling cannot conceal a collapsed original chord')
        unit = edge/lengths[:, None]
        turn = np.arctan2(unit[:-1, 0]*unit[1:, 1]-unit[:-1, 1]*unit[1:, 0],
                          np.sum(unit[:-1]*unit[1:], axis=1))
        keep = {0, len(original)-1, *(int(i+1) for i in np.flatnonzero(abs(turn) > .05))}
        protected = sorted(keep)
        pending = list(zip(protected, protected[1:]))
        while pending:
            first, last = pending.pop()
            if last-first < 2:
                continue
            delta = original[last]-original[first]
            fraction = np.clip((original[first+1:last]-original[first])@delta/(delta@delta), 0., 1.)
            error = np.linalg.norm(original[first+1:last]-original[first]-fraction[:, None]*delta, axis=1)
            at = int(np.argmax(error))
            if error[at] > allowance:
                index = first+1+at
                keep.add(index)
                pending.extend(((first, index), (index, last)))
        indices = sorted(keep)
        maximum = 0.
        for first, last in zip(indices, indices[1:]):
            delta = original[last]-original[first]
            fraction = np.clip((original[first:last+1]-original[first])@delta/(delta@delta), 0., 1.)
            maximum = max(maximum, float(np.linalg.norm(original[first:last+1]-original[first]
                -fraction[:, None]*delta, axis=1).max()))
        if maximum > allowance+1e-12:
            raise ValueError('Derived source curb exceeds its declared projection allowance')
        super().__init__(original[indices])
        self.original_polyline = original.copy()
        self.original_polyline.flags.writeable = False
        self.source_spans, self.source_curves = tuple(source_spans), tuple(source_curves)
        self.maximum_projection_error, self.projection_allowance = maximum, float(allowance)
        self.protected_original_indices = tuple(protected)
        parameters = self.stations/self.length
        if (len(parameters)-1) % 2:
            at = int(np.argmax(np.diff(parameters)))
            parameters = np.sort(np.r_[parameters, .5*(parameters[at]+parameters[at+1])])
        self.initial_parameters = parameters

    def report(self):
        return dict(reference='sampled_original_source_union_polyline',
            original_points=len(self.original_polyline), retained_points=len(self.points),
            maximum_projection_error=self.maximum_projection_error,
            additional_error_allowance=self.projection_allowance,
            original_analytical_primitives_retained=True, physical_profiles_exact=True,
            original_bank_tips_and_sharp_vertices_exact=True)


def original_local_chord_bounds(layout):
    """Bound the old prepared chords against their original analytical curbs.

    Fine certified source chords cover the analytical curve. Distance to a
    finite old chord is convex along each fine chord and 1-Lipschitz, so its
    maximum vertex distance plus the certified fine error bounds the complete
    original analytical interval, not only diagnostic sample points.
    """
    fine_error = 1e-7
    sampler = CertifiedCurbSampler(fine_error, .005, 20000)
    bounds = []
    for curve in layout.boundaries:
        parameters = getattr(curve, 'initial_parameters', None)
        if parameters is None:
            parameters = np.linspace(0., 1., layout.rows+1)
        maximum = 0.
        for first, last in zip(parameters, parameters[1:]):
            points = curve(sampler.sample(curve, float(first), float(last)))
            a, b = curve(first), curve(last)
            delta = b-a
            if float(delta@delta) <= 1e-20:
                raise ValueError('Original local analytical curb has a collapsed prepared chord')
            fraction = np.clip((points-a)@delta/(delta@delta), 0., 1.)
            maximum = max(maximum, float(np.linalg.norm(points-a-fraction[:, None]*delta, axis=1).max())+fine_error)
        bounds.append(maximum)
    return tuple(bounds)


def resample(group, context, *, maximum_projection=.0004):
    from .updates import drain_steps
    return drain_steps(resample_steps(group, context, maximum_projection=maximum_projection))


def resample_steps(group, context, *, maximum_projection=.0004):
    """Return a fresh prepared group with the same profiles and physical holes.

    The allowance is deducted from both the authored curve budget and the
    independent 2 mm complete-strip budget. Existing documented projection is
    never reset or hidden. No accepted mesh or physical mouth is modified.
    """
    from .one_sided_paired_junctions import certify_steps
    yield from certify_steps(group)
    layout = group.physical_layout
    if not math.isfinite(maximum_projection) or not 0. < maximum_projection <= .0004:
        raise ValueError('Budgeted source sampling exceeds its finite construction projection cap')
    if getattr(layout, 'provenance_rim_sampling', None):
        raise ValueError('A budgeted source rim cannot spend its reference allowance twice')
    nominal = context.network.settings.curve_tolerance
    old_bounds = original_local_chord_bounds(layout.original_junction_layout)
    old_outer = sum(max((row['maximum_projection_error'] for row in rows), default=0.) for rows in
                    (getattr(layout, 'source_rim_normalization', ()),))
    old_hole = max((row['maximum_projection_error'] for row in getattr(layout, 'source_hole_normalization', ())), default=0.)
    old_projection = old_outer+old_hole
    available = min(maximum_projection, nominal-max(old_bounds)-old_projection,
                    (.002-layout.source_curve_tolerance-old_projection)/2.)
    if not math.isfinite(available) or available <= 0.:
        raise ValueError('Original source/chord errors exhaust the remaining authored rim budget')
    sources = {source.region.owner: source for source in group.footprints}
    outer = next(loop for loop in group.loops if not loop.hole)
    labels = [rim_label(span, sources) for span in outer.spans]
    by_identity = {id(approach): i for i, approach in enumerate(layout.approaches)}
    tokens = []
    for span, label in zip(outer.spans, labels):
        approach = sources[span.owner].layout.approaches[label.index] if label.kind == 'mouth' else None
        if label.kind not in ('mouth', 'curb'):
            raise ValueError('Budgeted rim sampling cannot hide an internal source seam')
        tokens.append(by_identity.get(id(approach)) if approach is not None else None)
    starts = [i for i, token in enumerate(tokens) if token is not None and token != tokens[i-1]]
    if len(starts) != len(layout.approaches) or len({tokens[i] for i in starts}) != len(starts):
        raise ValueError('Budgeted source rim needs every complete original physical profile')
    owner = Key('budgeted_source_rim', group.members[0].id)
    if owner in sources:
        raise ValueError('Budgeted source rim repeats a derived physical owner')
    curves, points, segments, reports = [], [], [], []
    order = []
    for first, following in zip(starts, (*starts[1:], starts[0]+len(outer.xy))):
        token = tokens[first]
        order.append(token)
        at, u = first, 1.
        while at < following and tokens[at % len(tokens)] == token:
            label = labels[at % len(labels)]
            if abs(label.first-u) > 1e-8 or label.last >= label.first:
                raise ValueError('Budgeted source rim has a broken authoritative profile partition')
            points.append(outer.xy[at % len(outer.xy)])
            segments.append(RimSegment('mouth', len(order)-1, label.first, label.last))
            u = label.last
            at += 1
        if abs(u) > 1e-8 or at >= following:
            raise ValueError('Budgeted source rim clips an original physical profile')
        indices = np.arange(at, following+1) % len(outer.xy)
        old_spans = tuple(outer.spans[i % len(outer.xy)] for i in range(at, following))
        provenance = {(span.owner, labels[i % len(labels)].index):
                      sources[span.owner].layout.boundaries[labels[i % len(labels)].index]
                      for i, span in zip(range(at, following), old_spans)}
        curve = BudgetedSourceCurb(outer.xy[indices], old_spans,
            [(key, index, value) for (key, index), value in sorted(provenance.items())], available)
        index = len(curves)
        curves.append(curve)
        reports.append(curve.report())
        parameters = curve.initial_parameters
        points.extend(curve(parameters)[:-1])
        segments.extend(RimSegment('curb', index, float(a), float(b)) for a, b in zip(parameters, parameters[1:]))
    physical = copy(layout)
    physical.approaches = tuple(layout.approaches[i] for i in order)
    physical.outward_stations = tuple(layout.outward_stations[i] for i in order)
    checked = RegionUnion().build((Region(owner, np.asarray(points)),))
    if len(checked) != 1 or checked[0].hole:
        raise ValueError('Budgeted source sampling changed its physical exterior topology')
    # The union is a validation here. Reconstructing its endpoint arithmetic
    # can change a last profile coordinate by an ULP; retain the authoritative
    # copied profile rows and direct analytical labels exactly.
    loops = [BoundaryLoop(np.asarray(points), tuple(RimSpan(owner, label) for label in segments))]
    for loop in group.loops:
        if not loop.hole:
            continue
        centered = loop.xy-loop.xy.mean(axis=0)
        _, axes = np.linalg.eigh(centered.T@centered)
        along = centered@axes[:, -1]
        first, last = sorted((int(np.argmin(along)), int(np.argmax(along))))
        xy, spans = [], []
        for a, b in ((first, last), (last, first+len(loop.xy))):
            indices = np.arange(a, b+1) % len(loop.xy)
            old_spans = tuple(loop.spans[i % len(loop.xy)] for i in range(a, b))
            provenance = {}
            for span in old_spans:
                label = rim_label(span, sources)
                if label.kind != 'curb':
                    raise ValueError('A physical grass hole contains an unexpected incoming profile')
                provenance[span.owner, label.index] = sources[span.owner].layout.boundaries[label.index]
            curve = BudgetedSourceCurb(loop.xy[indices], old_spans,
                [(key, index, value) for (key, index), value in sorted(provenance.items())], available)
            curve.on_separator_rim = True
            index = len(curves)
            curves.append(curve)
            reports.append(curve.report())
            parameters = curve.initial_parameters
            xy.extend(curve(parameters)[:-1])
            spans.extend(RimSpan(owner, RimSegment('curb', index, float(x), float(y)))
                         for x, y in zip(parameters, parameters[1:]))
        hole = BoundaryLoop(np.asarray(xy), tuple(spans))
        if not hole.hole:
            raise ValueError('Budgeted source sampling changed a physical hole orientation')
        RegionUnion().build((Region(owner, hole.xy[::-1]),))
        loops.append(hole)
    physical.boundaries = tuple(curves)
    physical.rows = max(len(curve.points)-1 for curve in curves)
    maximum = max(row['maximum_projection_error'] for row in reports)
    ledger = dict(original_local_chord_bounds=old_bounds, previous_projection=old_projection,
        maximum_added_projection=maximum, local_curve_combined_bounds=[v+old_projection+maximum for v in old_bounds],
        authored_curve_budget=nominal, original_strip_chord_error=layout.source_curve_tolerance,
        independent_complete_width_budget=.002, curves=reports,
        exterior_projection=max(row['maximum_projection_error'] for row in reports[:len(order)]),
        hole_projection=max((row['maximum_projection_error'] for row in reports[len(order):]), default=0.))
    if max(ledger['local_curve_combined_bounds']) > nominal:
        raise ValueError('Budgeted source sampling exceeds the authored local analytical curve budget')
    physical.provenance_rim_sampling = ledger
    # Keep all old proof operands/layouts. Only this consolidated prepared rim
    # owns the newly sampled boundary; original source curves remain reachable.
    sampled = SampledFootprint(Region(owner, np.asarray(points)), tuple(segments), physical)
    result = replace(group, physical_layout=physical, footprints=(sampled, *group.footprints), loops=tuple(loops))
    yield from certify_steps(result)
    return result


def strict_sampling(group):
    from .updates import drain_steps
    return drain_steps(strict_sampling_steps(group))


def strict_sampling_steps(group):
    """Retain the established one-micrometre construction before a retry."""
    from .one_sided_paired_junctions import certify_steps
    from .one_sided_rim_sampling import normalize_holes_steps
    from .source_rim_normalization import normalize
    if any(loop.hole for loop in group.loops):
        group = yield from normalize_holes_steps(group)
        holes = tuple(loop for loop in group.loops if loop.hole)
        outer = normalize(replace(group, loops=tuple(loop for loop in group.loops if not loop.hole)))
        group = replace(outer, loops=(*outer.loops, *holes))
    else:
        yield from certify_steps(group)
        group = normalize(group)
    yield from certify_steps(group)
    return group


def grading_retry(chart_factory, group, context):
    from .updates import drain_steps
    return drain_steps(grading_retry_steps(chart_factory, group, context))


def grading_retry_steps(chart_factory, group, context):
    """Spend remaining reference budget only after a typed grading rejection."""
    from .partition import QuadChartGradingError

    def build(candidate):
        chart = chart_factory()
        if hasattr(chart, 'build_steps'):
            return (yield from chart.build_steps(candidate, context))
        return chart.build(candidate, context)

    strict = yield from strict_sampling_steps(group)
    yield
    try:
        return (yield from build(strict))
    except QuadChartGradingError as failure:
        yield
        candidate = yield from resample_steps(strict, context)
        yield
        output = yield from build(candidate)
        output.metadata = dict(output.metadata, provenance_rim_sampling=candidate.physical_layout.provenance_rim_sampling,
            previous_strict_chart_failure=str(failure))
        return output
