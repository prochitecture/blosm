"""Geometry-selected chart strategies with explicit eligibility and fallback.

These strategies own neither source profiles nor terrain/UV policy. A strategy
may only change the interior chart after qualifying the prescribed physical rim.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class ChartEligibility:
    eligible: bool
    reason: str


class FlatPairedRimStrategy:
    """Cut a long plateau into visible ribbon blocks before chart construction.

    Only an already constant plateau is eligible: changing its interior chart
    cannot introduce a different prepared height field. End curves and every
    physical rim sample remain authoritative. Nonvisible/short/nonflat regions
    use the established generic chart search.
    """
    name = 'paired-rim'

    def __init__(self, patterns, *, spacing, max_cells):
        self.patterns, self.spacing, self.max_cells = patterns, spacing, max_cells

    def classify(self, boundary, *, datum, height):
        if len(boundary) < 64:
            return ChartEligibility(False, 'small rim')
        ceiling = datum + height
        tolerance = min(1e-12, 8. * np.spacing(max(1., abs(datum), abs(ceiling))))
        if np.max(abs(boundary[:, 2] - ceiling)) > tolerance:
            return ChartEligibility(False, 'nonconstant plateau')
        xy = boundary[:, :2] - boundary[0, :2]
        centered = xy - xy.mean(axis=0)
        _, vectors = np.linalg.eigh(centered.T @ centered)
        extents = np.ptp(centered @ vectors, axis=0)
        if extents[0] <= 1e-7 or extents[1] < 12. * extents[0]:
            return ChartEligibility(False, 'not a slender ribbon')
        from .quad_charts import visibility_center
        try:
            visibility_center(boundary[:, :2])
        except ValueError:
            return ChartEligibility(False, 'no strict visibility kernel')
        return ChartEligibility(True, 'long visible constant plateau')

    def construct(self, boundary):
        from .separator_surfaces import PartitionedSeparatorChart, SeparatorStripCuts

        class PairedRimChart(PartitionedSeparatorChart):
            def cut_policy(self):
                # Eligibility already established a strict global kernel.
                # The original initial visibility split would have no cuts.
                return SeparatorStripCuts(self.spacing, initial_visibility=False)

        return PairedRimChart(self.patterns, spacing=self.spacing,
            interior_spacing=2. * self.spacing, max_cells=self.max_cells).build(boundary)

    def validate(self, result, boundary):
        from .model import Key
        from .validation import validate_overlaps
        mesh = result.mesh
        quality = mesh.validate()
        if quality['maximum_edge_ratio'] > 8. or len(mesh.quads) > self.max_cells:
            raise ValueError('Paired-rim chart exceeds the existing quality/cell limits')
        if not np.array_equal(mesh.xyz[:len(boundary), :2], boundary):
            raise ValueError('Paired-rim chart changed the physical boundary')
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        expected = np.sort(np.column_stack((np.arange(len(boundary)),
                             np.roll(np.arange(len(boundary)), -1))), axis=1)
        expected = expected[np.lexsort((expected[:, 1], expected[:, 0]))]
        if np.any(counts > 2) or not np.array_equal(unique[counts == 1], expected):
            raise ValueError('Paired-rim chart lost conforming boundary incidence')
        # Positive local cells alone do not establish global non-overlap.
        validate_overlaps({Key('separator_chart', 0): mesh})
        return result

    def build(self, boundary):
        return self.validate(self.construct(boundary), boundary)
