"""Opt-in bounded quad charts for certified one-sided paired preparation."""
from .conditions import BoundaryConditions
from .four_side_charts import quality
from .mesh import publication_coordinates
from .one_sided_paired_cuts import OneSidedPairedCuts
from .one_sided_paired_junctions import SourceEvaluationInterval, certify_steps
from .partition import CompoundPartitioner, CompoundQuadChart, QuadChartGradingError
from .provenance_rim_sampling import grading_retry_steps
from .updates import drain_steps


def near_intervals(layout):
    paired = []
    for approach, station in zip(layout.approaches, layout.outward_stations):
        if not hasattr(approach, 'near_role'):
            continue
        role = approach.near_role
        alias = approach.source_mouths(station)[role]
        if not alias.original_incidence:
            raise ValueError('One-sided chart lost its real original near-source incidence')
        path, visit = approach.paths[role], approach.visits[role]
        first, last = (0., alias.station) if alias.endpoint.end == 0 else (alias.station, path.length)
        paired.append((approach.end_index, SourceEvaluationInterval(path.section.id, first, last, visit.forward)))
    if len(paired) != 2 or {end for end, interval in paired} != {0, 1}:
        raise ValueError('One-sided chart needs both original near-source evaluation intervals')
    return tuple(interval for end, interval in sorted(paired, reverse=True))


class _SourceChart(CompoundQuadChart):
    def _repairable_leaf(self, group):
        return getattr(group.physical_layout, 'one_sided_source_preparation', False) or super()._repairable_leaf(group)


class _CutFamilies:
    def __init__(self, patterns, spacing, max_cells, quality_iterations):
        self.patterns, self.spacing = patterns, spacing
        self.max_cells, self.quality_iterations = max_cells, quality_iterations

    def build_steps(self, group, context):
        yield from certify_steps(group)
        layouts = {source.region.owner: source.layout for source in group.footprints}
        for loop in group.loops:
            BoundaryConditions().build(group, loop, layouts)
            yield
        failures, grading = [], None
        for existing in (False, True):
            policy = OneSidedPairedCuts(near_intervals(group.physical_layout), existing_vertices=existing)
            try:
                chart = _SourceChart(self.patterns, spacing=self.spacing, max_cells=self.max_cells,
                    quality_iterations=self.quality_iterations,
                    partitioner=CompoundPartitioner(policy=policy, require_visibility=False),
                    four_sided=True)
                result = yield from chart.build_steps(group, context)
                stored = quality(result.mesh.with_positions(publication_coordinates(result.mesh.xyz, 1., physical=True)))
                if stored['maximum_edge_ratio'] > 8. or stored['minimum_corner_sine'] < .03:
                    raise QuadChartGradingError('One-sided source atlas fails unchanged stored pavement quality')
            except ValueError as failure:
                failures.append(str(failure))
                if isinstance(failure, QuadChartGradingError):
                    grading = failure
                continue
            result.metadata = dict(result.metadata, one_sided_source_partition=True,
                existing_source_cut_vertices=existing, previous_one_sided_chart_failures=failures,
                stored_quality=stored)
            return result
        if grading is not None:
            raise QuadChartGradingError('One-sided computational cut families need better grading: '+str(grading),
                candidate=grading.candidate, details=dict(grading.details or {}, cut_failures=failures)) from grading
        raise ValueError('One-sided computational cut families failed: '+'; '.join(failures))


class OneSidedQuadChart:
    """Only explicitly prepared source groups select this bounded chart."""
    def __init__(self, patterns, *, spacing=.7, max_cells=50000, quality_iterations=300):
        self.patterns, self.spacing = patterns, spacing
        self.max_cells, self.quality_iterations = max_cells, quality_iterations

    def build(self, group, context):
        return drain_steps(self.build_steps(group, context))

    def build_steps(self, group, context):
        if not getattr(group.physical_layout, 'one_sided_source_preparation', False):
            raise ValueError('One-sided quad charts require their complete prepared-source proof')
        return (yield from grading_retry_steps(lambda: _CutFamilies(self.patterns, self.spacing,
            self.max_cells, self.quality_iterations), group, context))
