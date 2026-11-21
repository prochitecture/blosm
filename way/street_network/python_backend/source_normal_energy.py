"""Extra energy seeds only inside the final source-certified normal atlas."""
import numpy as np
from .four_side_charts import quality
from .mesh import publication_coordinates
from .updates import drain_steps
from .partition import CompoundQuadChart, QuadChartGradingError
from .quad_charts import FixedRimQuality


def _energy_attempt(parent, original, rejected, power, retry_rotations):
    return drain_steps(_energy_attempt_steps(parent, original, rejected, power, retry_rotations))


def _energy_attempt_steps(parent, original, rejected, power, retry_rotations):
    original.validate()
    rejected.validate()
    if (max(len(original.quads), len(rejected.quads)) > parent.max_cells
            or max(float(np.ptp(original.xyz[:, 2])), float(np.ptp(rejected.xyz[:, 2]))) > 1e-10):
        raise ValueError('Source-normal energy repair requires a bounded planar candidate')
    seed, energy = yield from FixedRimQuality(iterations=500, power=power).improve_steps(rejected)
    correction_chart = CompoundQuadChart(parent.patterns, spacing=parent.spacing,
        max_cells=parent.max_cells, quality_iterations=300)
    result, correction = yield from correction_chart._correct_steps(seed, retry_rotations=retry_rotations)
    physical = quality(result)
    stored = quality(result.with_positions(publication_coordinates(result.xyz, 1., physical=True)))
    if (max(physical['maximum_edge_ratio'], stored['maximum_edge_ratio']) > 8.
            or min(physical['minimum_corner_sine'], stored['minimum_corner_sine']) < .03):
        raise ValueError('Source-normal energy seed fails unchanged physical or stored quality')
    if (set(result.seams) != set(original.seams)
            or any(not np.array_equal(result.seams[key], ids)
                   or not np.array_equal(result.xyz[ids], original.xyz[ids])
                   for key, ids in original.seams.items())):
        raise ValueError('Source-normal energy repair moved a physical rim or its source vertex IDs')
    return result, correction, energy, physical, stored


class EnergySourceNormalChart(CompoundQuadChart):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.energy_recoveries = []

    def _correct_steps(self, mesh, *, retry_rotations=False):
        try:
            return (yield from super()._correct_steps(mesh, retry_rotations=retry_rotations))
        except QuadChartGradingError as original:
            if mesh.material == 'grass' or original.candidate is None:
                raise
            failures = []
            # The ordinary source child keeps its cheaper search. This last
            # family reuses the general chart's finite 300-step relaxation.
            for power in (1, 2):
                try:
                    result, correction, energy, physical, stored = yield from _energy_attempt_steps(
                        self, mesh, original.candidate, power, retry_rotations)
                except ValueError as failure:
                    failures.append(dict(power=power, reason=str(failure)))
                    continue
                details = dict(power=power, energy_iterations=500, correction_energy_iterations=300, energy=energy,
                    previous_correction=original.details, previous_energy_failures=failures,
                    physical_quality=physical, stored_quality=stored, fixed_source_rim=True,
                    cell_budget=self.max_cells)
                self.energy_recoveries.append(details)
                return result, dict(correction, source_normal_energy_recovery=details)
            original.details = dict(original.details or {}, failed_source_normal_energy_seeds=failures)
            raise original
