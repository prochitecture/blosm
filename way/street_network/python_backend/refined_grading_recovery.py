"""Bounded energy seeds after fixed-topology refined shape correction fails."""
import numpy as np
from .four_side_charts import quality
from .mesh import publication_coordinates
from .quad_charts import FixedRimQuality
from .quality_correction import QuadShapeCorrection


def recover(mesh, fixed_vertices, max_cells, previous):
    mesh.validate()
    if (type(max_cells) is not int or not len(mesh.quads) <= max_cells <= 50000
            or np.ptp(mesh.xyz[:, 2]) > 1e-10):
        raise ValueError('Refined grading recovery requires a bounded planar mesh')
    fixed = np.asarray(tuple(fixed_vertices))
    if (fixed.ndim != 1 or not np.issubdtype(fixed.dtype, np.integer)
            or np.any(fixed < 0) or np.any(fixed >= len(mesh.xyz))):
        raise ValueError('Refined grading recovery lost its authoritative boundary vertices')
    failures = []
    for power in (1, 2):
        try:
            seed, energy = FixedRimQuality(iterations=500, power=power).improve(mesh)
            if not np.array_equal(seed.xyz[fixed], mesh.xyz[fixed]):
                raise ValueError('Refined energy seed moved an authoritative boundary vertex')
            candidate, correction = QuadShapeCorrection(max_cells=max_cells).improve(seed, fixed_vertices=fixed)
            if not correction['qualified']:
                raise ValueError('Refined energy seed did not meet the shape limits')
            double = quality(candidate)
            stored = quality(candidate.with_positions(publication_coordinates(candidate.xyz, 1., physical=True)))
            if (max(double['maximum_edge_ratio'], stored['maximum_edge_ratio']) > 8.
                    or min(double['minimum_corner_sine'], stored['minimum_corner_sine']) < .03):
                raise ValueError('Refined energy seed fails physical or stored pavement quality')
            if (not np.array_equal(candidate.quads, mesh.quads)
                    or not np.array_equal(candidate.xyz[fixed], mesh.xyz[fixed])
                    or set(candidate.seams) != set(mesh.seams)
                    or any(not np.array_equal(candidate.seams[key], ids) for key, ids in mesh.seams.items())):
                raise ValueError('Refined energy correction changed fixed topology or source boundary')
        except ValueError as failure:
            failures.append(dict(power=power, reason=str(failure)))
            continue
        return candidate, dict(correction, refined_energy_recovery=dict(
            previous_correction=previous, previous_energy_failures=failures,
            power=power, energy=energy, physical_quality=double, stored_quality=stored,
            fixed_topology=True, fixed_source_boundary=True, energy_iterations=500,
            cell_budget=max_cells))
    raise ValueError('Refined grading energy seeds exhausted their bounded shape alternatives: '+str(failures))
