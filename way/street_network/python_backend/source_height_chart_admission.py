"""Necessary original-height feasibility for the failed ordinary prefix atlas.

This chooses no physical input and solves no heights. An optional Surface caller
can reject a candidate whose original immutable data already contradict its
configured entrance limit. All authoritative final guards still run normally.
"""
import numpy as np
from .entrance_height_recovery import entrance_slabs
from .source_body_references import BodyGuideLimits
from .source_body_heights import prepare_steps
from .updates import drain_steps
from .junction_bank_transitions import certify_body_domains


def entrance_interval_floor(slabs, heights, immutable, datum, ceiling):
    """Lower-bound |a.z| over the original height box and exact known nodes.

    Free heights vary independently in [datum,ceiling], giving a necessary,
    deliberately generous feasibility interval. Floating arithmetic is expanded
    outwards before rejecting any candidate; solver tolerances are unchanged.
    """
    known = immutable[slabs.vertices]
    coefficients = slabs.coefficients
    fixed = np.where(known, coefficients*heights[slabs.vertices], 0.)
    low = np.where(known, fixed, coefficients*np.where(coefficients >= 0., datum, ceiling))
    high = np.where(known, fixed, coefficients*np.where(coefficients >= 0., ceiling, datum))
    magnitude = np.maximum(abs(low), abs(high)).sum(axis=1)
    eps = np.finfo(float).eps
    error = (8.*eps/(1.-8.*eps))*magnitude
    first = np.nextafter(low.sum(axis=1)-error, -np.inf)
    last = np.nextafter(high.sum(axis=1)+error, np.inf)
    return np.maximum(0., np.maximum(first, -last))


def admission_steps(admission, group, chart):
    if isinstance(admission, SourceHeightChartAdmission):
        return (yield from admission.check_steps(group, chart))
    return admission(group, chart)


class SourceHeightChartAdmission:
    def __init__(self, builder, context):
        self.builder, self.context = builder, context

    def __call__(self, group, chart):
        return drain_steps(self.check_steps(group, chart))

    def check_steps(self, group, chart):
        builder, context = self.builder, self.context
        conditions, _, limit = builder.boundary_conditions(group, context)
        # Necessary row feasibility does not need an interior field or a jet
        # approximation. Give free nodes an arbitrary allowed height; only
        # actual physical profile/bank nodes and exact native body targets
        # enter the immutable part of each affine row.
        datum, ceiling = conditions[0].datum, conditions[0].ceiling
        if (not np.isfinite((datum, ceiling)).all() or ceiling < datum
                or any(abs(data.datum-datum)>1e-8 or abs(data.ceiling-ceiling)>1e-8 for data in conditions)):
            raise ValueError('Source-height chart admission needs the actual common elevation envelope')
        mesh = chart.mesh.with_positions(chart.mesh.xyz.copy())
        mesh.xyz[:, 2] = datum
        for node, ref in chart.boundary_refs.items():
            mesh.xyz[node, 2] = conditions[ref.loop].edge_samples(ref.edge, ref.fraction)[1]
        grass_cells = sum(len(hole.mesh.quads) for hole in chart.holes.values())
        guides = yield from prepare_steps(context, group, mesh, chart.boundary_refs,
            limits=BodyGuideLimits(maximum_mesh_cells=min(50000, builder.max_cells-grass_cells)))
        certify_body_domains(guides, conditions)
        fixed = np.zeros(len(mesh.xyz), bool)
        fixed[list(chart.boundary_refs)] = True
        immutable = fixed.copy()
        immutable[guides.core_nodes] = True
        rows = entrance_slabs(mesh, fixed, conditions, builder.maximum_join_slope, at_corners=False)
        floors = entrance_interval_floor(rows, guides.seed.xyz[:, 2], immutable, datum, ceiling)
        maximum = float(floors.max())
        if maximum > builder.maximum_join_slope:
            raise ValueError('Original immutable profile/body data make this chart entrance infeasible: '
                             f'{maximum!r} > {builder.maximum_join_slope!r}')
        return dict(construction='original-height-box entrance midpoint necessary feasibility',
            physical_limit=limit, entrance_limit=builder.maximum_join_slope,
            maximum_necessary_entrance_slope=maximum, entrance_rows=len(floors),
            original_body_nodes=len(guides.core_nodes), physical_boundary_nodes=int(fixed.sum()),
            source_profiles_unchanged=True, exact_native_body_targets=True,
            solver_ran=False, mandatory_final_surface_checks=True)
