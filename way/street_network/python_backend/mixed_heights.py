"""Bounded interior height reconciliation for otherwise valid mixed joins.

This is a fallback for matching full/half crowns. Source loci, direct quads,
complete mouths plus their neighboring rows, outer curbs and the whole crest
remain fixed. The interior may redistribute its profile while retaining a
downhill transverse route to a zero-height curb and the ordinary slope ceiling.
"""
from dataclasses import replace
import numpy as np
from .drainage import BoundaryDrainage
from .slope_projection import SlopeProjection


class MixedHeightCorrection:
    max_candidates = 3
    max_cells = 8192
    memory_bytes = 24 * 1024 * 1024
    max_iterations = 4000

    @staticmethod
    def compatible(plan):
        a, b = plan.approaches
        q = np.linspace(-1., 1., 5)
        return (abs(a.profile.datum - b.profile.datum) <= 1e-8
                and abs(a.profile.peak + b.profile.peak) <= 1e-12
                and np.max(abs(a.profile.value(-q) - b.profile.value(q))) <= 1e-12)

    @classmethod
    def score(cls, context, plan):
        """Keep only geometrically valid bounded candidates for the fallback."""
        if plan.rows * plan.across > cls.max_cells or not cls.compatible(plan):
            return None
        try:
            mesh = plan.mesh(context.patterns, 1.)
            if mesh.validate()['maximum_edge_ratio'] > 8.:
                return None
            return mesh.slope_statistics()['maximum']
        except ValueError:
            return None

    @classmethod
    def solve(cls, context, plan):
        if not cls.compatible(plan):
            raise ValueError('Mixed height correction requires matching oriented crown profiles')
        mesh = plan.mesh(context.patterns, 1.)
        grid = np.arange(len(mesh.xyz)).reshape(plan.rows + 1, plan.across + 1)
        peak = int(round(.5 * (plan.approaches[1].profile.peak + 1.) * plan.across))
        fixed = np.zeros(len(mesh.xyz), bool)
        fixed[grid[:2].ravel()] = True
        fixed[grid[-2:].ravel()] = True
        fixed[grid[:, (0, peak, plan.across)].ravel()] = True
        incoming = max(float(abs(a.profile.derivative(np.array((-1., 1.)))).max()) * 2. / a.width
                       for a in plan.approaches)
        limit = max(.03, 1.1 * incoming)
        projected, report = SlopeProjection(max_iterations=cls.max_iterations, max_cells=cls.max_cells,
            memory_bytes=cls.memory_bytes).solve(mesh, fixed, datum=plan.datum, ceiling=plan.ceiling,
                maximum_slope=.995 * limit, feasibility_slope=.999999 * limit)
        if (not np.array_equal(projected.xyz[:, :2], mesh.xyz[:, :2])
                or not np.array_equal(projected.quads, mesh.quads)
                or not np.array_equal(projected.xyz[fixed], mesh.xyz[fixed])):
            raise ValueError('Mixed height projection changed fixed geometry')
        heights = projected.xyz[:, 2].reshape(grid.shape)
        if ((peak and np.any(np.diff(heights[:, :peak + 1], axis=1) < -1e-10))
                or (peak < plan.across and np.any(np.diff(heights[:, peak:], axis=1) > 1e-10))):
            raise ValueError('Mixed height correction introduced a transverse depression')
        outlets = np.zeros(len(mesh.xyz), bool)
        for column in (0, plan.across):
            ids = grid[:, column]
            outlets[ids] = abs(projected.xyz[ids, 2] - plan.datum) <= 1e-10
        _, drainage = BoundaryDrainage().fill(projected, outlets)
        if drainage['maximum_fill'] > 1e-10:
            raise ValueError('Mixed height correction introduced a drainage basin')
        residual = (projected.xyz[:, 2] - mesh.xyz[:, 2]).reshape(grid.shape)
        residual.flags.writeable = False
        return replace(plan, height_corrections=residual,
            certification=dict(plan.certification, height_reconciliation=report,
                fixed_mouth_neighbor_rows=True, fixed_curbs_and_crest=True,
                transverse_drainage=True, drainage=drainage))
