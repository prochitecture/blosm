"""Keep physical slope bounds after Blender's float32 vertex storage.

Only an actual rounded-buffer violation triggers a correction. XY, topology,
boundary profiles and curb heights remain fixed; the numerical query uses the
same corrected physical heights as the mesh later sent for publication.
"""
import numpy as np
from .drainage import BoundaryDrainage
from .mesh import Mesh
from .slope_projection import SlopeProjection, ProjectionIterationBudgetError, corner_matrices
from .surface_query import QuadHeightField
from .updates import drain_steps, numerical_time_slices


class PublicationSlopeGuard:
    tolerance = 1e-7
    maximum_rounds = 8

    @staticmethod
    def fixed_gradient_floor(mesh, fixed, datum):
        """Maximum unavoidable corner gradient for this fixed-height mask."""
        return float(PublicationSlopeGuard.fixed_gradient_floors(mesh,fixed,datum).max())

    @staticmethod
    def fixed_gradient_floors(mesh, fixed, datum):
        """Per-corner minimum slopes allowed by immutable boundary heights.

        Free vertex heights can change a corner gradient only in the range of
        its two-by-four free-height matrix. The component in the orthogonal
        complement is fixed. A global roundoff allowance must not request a
        target below that component; doing so can reject a representable mesh
        before its actual stored slopes are tested.
        """
        coefficient = corner_matrices(mesh)
        free = ~np.asarray(fixed)[mesh.quads]
        movable = coefficient * free[:, None, None, :]
        gram = np.einsum('fcij,fckj->fcik', movable, movable)
        values, vectors = np.linalg.eigh(gram)
        gradient = np.einsum('fcij,fj->fci', coefficient,
                            mesh.xyz[mesh.quads, 2] - datum)
        component = np.einsum('fcji,fcj->fci', vectors, gradient)
        # Use the same rank convention as the authoritative projection. The
        # zero-rank case has two immutable gradient components.
        largest = np.maximum(values[..., 1], 0.)
        zero = largest <= 1e-20
        null = np.stack((zero | (values[..., 0] <= largest * 1e-12), zero), axis=-1)
        return np.linalg.norm(np.where(null, component, 0.), axis=-1)

    @staticmethod
    def physical_buffer(mesh, datum, exaggeration):
        xyz = mesh.xyz.copy()
        xyz[:, 2] = datum + exaggeration * (xyz[:, 2] - datum)
        xyz = xyz.astype(np.float32).astype(np.float64)
        xyz[:, 2] = datum + (xyz[:, 2] - datum) / exaggeration
        return Mesh(xyz, mesh.quads, mesh.material, dict(mesh.seams))

    def correct(self, mesh, field, fixed, *, exaggeration, limit, max_cells, check_joins,
                allow_projection=True,height_slabs=None):
        return drain_steps(self.correct_steps(mesh, field, fixed, exaggeration=exaggeration,
            limit=limit, max_cells=max_cells, check_joins=check_joins,
            allow_projection=allow_projection, height_slabs=height_slabs))

    def correct_steps(self, mesh, field, fixed, *, exaggeration, limit, max_cells, check_joins,
                      allow_projection=True,height_slabs=None):
        slab_tolerance=None
        if height_slabs is not None:
            from .height_slabs import HeightSlabs
            if not isinstance(height_slabs,HeightSlabs):
                raise ValueError('Publication directional constraints require HeightSlabs')
            slab_tolerance=SlopeProjection().slope_tolerance
        if exaggeration == 0.:
            return mesh, field, None
        rounded = self.physical_buffer(mesh, field.datum, exaggeration)
        initial = float(rounded.quad_slopes().max())
        if (initial <= limit + self.tolerance
                and (height_slabs is None or height_slabs.violation(mesh.xyz[:,2])<=slab_tolerance)):
            if height_slabs is not None:check_joins(mesh)
            return mesh, field, None
        if not allow_projection:
            raise ValueError('Publication slope correction is required but height projection is disabled')
        original = mesh
        fixed_floor = self.fixed_gradient_floor(mesh, fixed, field.datum)
        stored_fixed_floor = self.fixed_gradient_floor(rounded, fixed, field.datum)
        if max(fixed_floor, stored_fixed_floor) > limit + self.tolerance:
            raise ValueError('Fixed boundary gradients cannot meet the published physical slope limit')
        floors = np.full(len(mesh.xyz), field.datum)
        rounds = []
        continued_after_spill = False
        fallback_targets = 0
        previous_target = None
        for iteration in range(self.maximum_rounds):
            def gradients(candidate):
                z = candidate.xyz[candidate.quads, 2] - field.datum
                return np.einsum('fcij,fj->fci', corner_matrices(candidate), z)
            error = float(np.linalg.norm(gradients(rounded) - gradients(mesh), axis=-1).max())
            margin = max(1e-6, 4. * error)
            target = limit - margin
            # Retain the established target whenever it is feasible. Otherwise
            # reserve one quarter of the fixed-gradient headroom. This is
            # a search target, never an acceptance tolerance: the actual stored
            # buffer must still pass the unchanged limit below.
            if (target <= fixed_floor + 1e-9
                    or (fallback_targets and target >= previous_target)):
                # Keep the established first target exactly. Repeating that
                # infeasible stored target cannot improve rounding: spend
                # successively less of the same immutable-gradient headroom.
                fraction = .75 if fallback_targets == 0 else .5 ** fallback_targets
                target = fixed_floor + fraction * (limit - fixed_floor)
                fallback_targets += 1
                if fallback_targets > 1 and target >= previous_target:
                    target = fixed_floor + .5 * (previous_target - fixed_floor)
            if target <= 0. or target <= fixed_floor:
                raise ValueError('Published coordinate precision cannot represent the physical slope limit')
            previous_target = target
            exhausted = None
            try:
                directional={} if height_slabs is None else dict(height_slabs=height_slabs)
                mesh, projection = yield from numerical_time_slices(SlopeProjection(max_cells=max_cells).solve_steps(mesh, fixed,
                    datum=field.datum, ceiling=field.ceiling, maximum_slope=target,
                    minimum_heights=floors,**directional))
            except ProjectionIterationBudgetError as failure:
                # The conservative inner target is not an acceptance limit.
                # A last iterate can already satisfy the actual stored mesh
                # contract before nearest-displacement convergence. Preserve
                # the default solver's rejection unless every physical box,
                # publication, drainage and join condition can be certified.
                candidate=failure.candidate
                if (not np.array_equal(candidate.xyz[fixed],original.xyz[fixed])
                        or not np.array_equal(candidate.xyz[:,:2],original.xyz[:,:2])
                        or not np.array_equal(candidate.quads,original.quads)
                        or np.any(candidate.xyz[:,2]<floors-1e-10)
                        or np.any(candidate.xyz[:,2]>field.ceiling+1e-10)
                        or float(candidate.quad_slopes().max())>limit+self.tolerance
                        or float(self.physical_buffer(candidate,field.datum,exaggeration).quad_slopes().max())
                           >limit+self.tolerance
                        or (height_slabs is not None
                            and height_slabs.violation(candidate.xyz[:,2])>slab_tolerance)):
                    raise
                exhausted = failure
                mesh=candidate
                projection=dict(failure.report,convergence='independently_qualified_published_candidate',
                    conservative_target_converged=False)
            filled, drainage = BoundaryDrainage(max_cells=max_cells).fill(mesh, fixed)
            rounded = self.physical_buffer(mesh, field.datum, exaggeration)
            maximum = float(rounded.quad_slopes().max())
            rounds.append(dict(roundoff_gradient_error=error, target=target,
                               fixed_gradient_floor=fixed_floor,
                               stored_fixed_gradient_floor=stored_fixed_floor,
                               published_maximum=maximum, projection=projection, drainage=drainage))
            if ((exhausted is not None or iteration+1==self.maximum_rounds) and maximum<=limit+self.tolerance
                    and drainage['maximum_fill']>1e-6):
                # Reprojection can recreate a small depression even when the
                # actual spill surface itself satisfies every slope bound.
                # Test that surface before exhausting the existing budget.
                xyz=mesh.xyz.copy()
                xyz[:,2]=filled
                candidate=mesh.with_positions(xyz)
                stored=self.physical_buffer(candidate,field.datum,exaggeration)
                stored_max=float(stored.quad_slopes().max())
                if (np.array_equal(candidate.xyz[fixed],original.xyz[fixed])
                        and np.all((filled>=floors-1e-10)&(filled<=field.ceiling+1e-10))
                        and float(candidate.quad_slopes().max())<=limit+self.tolerance
                        and stored_max<=limit+self.tolerance
                        and (height_slabs is None or height_slabs.violation(filled)<=slab_tolerance)):
                    drainage_guard=BoundaryDrainage(max_cells=max_cells)
                    _,physical_drainage=drainage_guard.fill(candidate,fixed)
                    _,stored_drainage=drainage_guard.fill(stored,fixed)
                    if max(physical_drainage['maximum_fill'],stored_drainage['maximum_fill'])<=1e-6:
                        mesh,maximum,drainage=candidate,stored_max,physical_drainage
                        rounds[-1].update(accepted_spill_surface=True,
                            spill_surface_published_maximum=maximum,spill_surface_stored_drainage=stored_drainage)
                if exhausted is not None and not rounds[-1].get('accepted_spill_surface'):
                    # The iterate meets both slope bounds. Feed its spill
                    # floors into the next ordinary round if the height
                    # envelope and any authored slabs permit the fill.
                    if (iteration+1 == self.maximum_rounds or np.any(filled > field.ceiling + 1e-10)
                            or (height_slabs is not None and height_slabs.violation(filled)>slab_tolerance)):
                        raise exhausted
                    continued_after_spill = True
                    rounds[-1]['continued_with_spill_floors'] = True
                    projection['convergence'] = 'feasible_slopes_pending_drainage'
            if (maximum <= limit + self.tolerance and drainage['maximum_fill'] <= 1e-6
                    and (height_slabs is None or height_slabs.violation(mesh.xyz[:,2])<=slab_tolerance)):
                if continued_after_spill:
                    stored = self.physical_buffer(mesh, field.datum, exaggeration)
                    _, stored_drainage = BoundaryDrainage(max_cells=max_cells).fill(stored, fixed)
                    if stored_drainage['maximum_fill'] > 1e-6:
                        raise ValueError('Projected spill floors leave a depression after vertex storage')
                    rounds[-1]['stored_drainage'] = stored_drainage
                check_joins(mesh)
                if (not np.array_equal(mesh.xyz[fixed], original.xyz[fixed])
                        or not np.array_equal(mesh.xyz[:, :2], original.xyz[:, :2])
                        or not np.array_equal(mesh.quads, original.quads)):
                    raise ValueError('Publication slope correction changed its fixed geometry')
                details = dict(initial_published_maximum=initial, published_maximum=maximum,
                    maximum_height_change=float(abs(mesh.xyz[:, 2] - original.xyz[:, 2]).max()),
                    exact_boundary=True, exact_xy=True, exact_quads=True, rounds=rounds)
                if height_slabs is not None:
                    details['directional_height_slabs']=dict(rows=len(height_slabs.vertices),
                        maximum_violation=height_slabs.violation(mesh.xyz[:,2]))
                query = QuadHeightField(mesh, datum=field.datum, ceiling=field.ceiling,
                    metadata=dict(field.metadata, publication_roundoff=details), max_cells=max_cells)
                return mesh, query, details
            raised = filled > mesh.xyz[:, 2] + 1e-10
            floors[raised] = np.maximum(floors[raised], filled[raised])
        raise ValueError('Publication slope correction exceeded its iteration budget')
