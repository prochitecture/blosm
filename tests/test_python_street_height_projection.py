import gc
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import weakref
from dataclasses import replace
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.slope_projection import SlopeProjection, corner_matrices
from tested_python_streets.surface_query import QuadHeightField
from tested_python_streets.ownership import CompoundProduct


def peak_grid(datum=0.):
    pattern = PatternLibrary().grid(2, 2)
    xyz = np.column_stack((2. * pattern.uv, np.full(len(pattern.uv), datum)))
    xyz[4, 2] += .08
    fixed = np.ones(len(xyz), bool)
    fixed[4] = False
    return Mesh(xyz, pattern.quads), fixed


def tilted_peak_grid():
    pattern=PatternLibrary().grid(4,4)
    xyz=np.c_[4.*pattern.uv,.03*4.*pattern.uv[:,0]]
    fixed=np.any((pattern.uv==0.)|(pattern.uv==1.),axis=1)
    xyz[~fixed,2]+=.08
    return Mesh(xyz,pattern.quads),fixed


class PythonHeightProjectionTests(unittest.TestCase):
    def test_projection_can_be_cancelled_between_sweeps_without_changing_source_heights(self):
        from tested_python_streets.updates import drain_steps
        mesh, fixed = peak_grid()
        original = mesh.xyz.copy()
        solver = SlopeProjection()
        options = dict(datum=0., ceiling=.08, maximum_slope=.05)
        steps = solver.solve_steps(mesh, fixed, **options)
        next(steps)
        steps.close()
        np.testing.assert_array_equal(mesh.xyz, original)
        resumed, evidence = drain_steps(solver.solve_steps(mesh, fixed, **options))
        synchronous, expected = solver.solve(mesh, fixed, **options)
        np.testing.assert_array_equal(resumed.xyz, synchronous.xyz)
        self.assertEqual(evidence, expected)

    def test_coupled_successful_default_keeps_exact_result_without_retry(self):
        from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
        mesh,fixed=peak_grid();floors=np.zeros(len(mesh.xyz));solver=SlopeProjection()
        expected,expected_report=solver.solve(mesh,fixed,datum=0.,ceiling=.08,
            maximum_slope=.05,minimum_heights=floors)
        with patch('tested_python_streets.compound_surfaces.PublicationSlopeGuard.fixed_gradient_floor',
                side_effect=AssertionError('Successful default entered fallback')):
            result,report=CompoundSurfaceBuilder._coupled_projection(solver,mesh,fixed,
                datum=0.,ceiling=.08,limit=.05,floors=floors)
        np.testing.assert_array_equal(result.xyz,expected.xyz)
        np.testing.assert_array_equal(result.quads,expected.quads)
        self.assertEqual(report,expected_report)

    def test_coupled_budget_retry_uses_strict_feasible_mode_above_real_fixed_floor(self):
        from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
        from tested_python_streets.slope_projection import ProjectionIterationBudgetError
        mesh,fixed=tilted_peak_grid();before=mesh.xyz.copy();floors=np.zeros(len(mesh.xyz))
        solver=SlopeProjection(max_iterations=3)
        with self.assertRaises(ProjectionIterationBudgetError) as failed:
            solver.solve(mesh,fixed,datum=0.,ceiling=.2,maximum_slope=.05,minimum_heights=floors)
        self.assertGreater(float(failed.exception.candidate.quad_slopes().max()),.05+1e-8)
        result,report=CompoundSurfaceBuilder._coupled_projection(solver,mesh,fixed,
            datum=0.,ceiling=.2,limit=.05,floors=floors)
        retry=report['coupled_feasible_retry']
        self.assertAlmostEqual(retry['fixed_gradient_floor'],.03,places=13)
        self.assertGreater(retry['inner_target'],retry['fixed_gradient_floor'])
        self.assertLess(retry['inner_target'],.05)
        self.assertEqual(retry['physical_acceptance_limit'],.05)
        self.assertEqual(report['convergence'],'feasible_surface')
        self.assertEqual(retry['maximum_iterations'],3)
        self.assertEqual(retry['memory_bytes'],solver.memory_bytes)
        self.assertLessEqual(float(result.quad_slopes().max()),.05)
        np.testing.assert_array_equal(result.xyz[fixed],before[fixed])
        np.testing.assert_array_equal(result.xyz[:,:2],before[:,:2])
        np.testing.assert_array_equal(result.quads,mesh.quads)
        np.testing.assert_array_equal(mesh.xyz,before)
        self.assertTrue(np.all((result.xyz[:,2]>=floors)&(result.xyz[:,2]<=.2)))

    def test_coupled_fixed_floor_without_headroom_and_unrelated_errors_remain_fatal(self):
        from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
        from tested_python_streets.slope_projection import ProjectionIterationBudgetError
        mesh,fixed=peak_grid();fixed[:]=True;floors=np.zeros(len(mesh.xyz))
        solver=SlopeProjection(max_iterations=0)
        with patch.object(solver,'solve_steps',wraps=solver.solve_steps) as calls:
            with self.assertRaises(ProjectionIterationBudgetError):
                CompoundSurfaceBuilder._coupled_projection(solver,mesh,fixed,
                    datum=0.,ceiling=.08,limit=.05,floors=floors)
        self.assertEqual(calls.call_count,1)
        problem=ValueError('Unrelated source-height constraint')
        with patch.object(solver,'solve_steps',side_effect=problem) as calls:
            with self.assertRaises(ValueError) as caught:
                CompoundSurfaceBuilder._coupled_projection(solver,mesh,fixed,
                    datum=0.,ceiling=.08,limit=.05,floors=floors)
        self.assertIs(caught.exception,problem);self.assertEqual(calls.call_count,1)

    def test_coupled_retry_rejects_nonfinite_floor_and_moved_original_source(self):
        from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
        from tested_python_streets.slope_projection import ProjectionIterationBudgetError
        mesh,fixed=peak_grid();floors=np.zeros(len(mesh.xyz))
        failure=ProjectionIterationBudgetError('Budget exhausted',mesh,{})
        solver=SlopeProjection(max_iterations=0)
        with patch.object(solver,'solve_steps',side_effect=failure) as calls:
            with patch('tested_python_streets.compound_surfaces.PublicationSlopeGuard.fixed_gradient_floor',return_value=np.nan):
                with self.assertRaises(ProjectionIterationBudgetError) as caught:
                    CompoundSurfaceBuilder._coupled_projection(solver,mesh,fixed,
                        datum=0.,ceiling=.08,limit=.05,floors=floors)
        self.assertIs(caught.exception,failure);self.assertEqual(calls.call_count,1)
        moved=mesh.with_positions(mesh.xyz.copy());moved.xyz[0,2]+=.001
        failure=ProjectionIterationBudgetError('Invalid source candidate',moved,{})
        with patch.object(solver,'solve_steps',side_effect=failure) as calls:
            with self.assertRaises(ProjectionIterationBudgetError) as caught:
                CompoundSurfaceBuilder._coupled_projection(solver,mesh,fixed,
                    datum=0.,ceiling=.08,limit=.05,floors=floors)
        self.assertIs(caught.exception,failure);self.assertEqual(calls.call_count,1)

    def test_coupled_feasible_retry_cannot_accept_a_looser_slope_tolerance(self):
        from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
        from tested_python_streets.slope_projection import ProjectionIterationBudgetError
        mesh,fixed=peak_grid();floors=np.zeros(len(mesh.xyz))
        failure=ProjectionIterationBudgetError('Budget exhausted',mesh,{})
        almost=mesh.with_positions(mesh.xyz.copy());almost.xyz[4,2]=(.05+5e-9)/np.sqrt(2.)
        solver=SlopeProjection(max_iterations=1)
        def almost_feasible():
            yield
            return almost,{}
        with patch.object(solver,'solve_steps',side_effect=[failure,almost_feasible()]) as calls:
            with self.assertRaisesRegex(ValueError,'unchanged source or physical constraints'):
                CompoundSurfaceBuilder._coupled_projection(solver,mesh,fixed,
                    datum=0.,ceiling=.08,limit=.05,floors=floors)
        self.assertEqual(calls.call_count,2)

    def test_optional_feasibility_stop_keeps_physical_constraints_and_default_projection(self):
        mesh, fixed = peak_grid()
        result, report = SlopeProjection().solve(mesh, fixed, datum=0., ceiling=.08,
            maximum_slope=.049, feasibility_slope=.05)
        self.assertEqual(report['convergence'], 'feasible_surface')
        self.assertLessEqual(result.quad_slopes().max(), .05)
        np.testing.assert_array_equal(result.xyz[fixed], mesh.xyz[fixed])
        self.assertGreaterEqual(float(result.xyz[:, 2].min()), 0.)
        self.assertLessEqual(float(result.xyz[:, 2].max()), .08)
        for bound in (.048, float('nan'), float('inf')):
            with self.assertRaisesRegex(ValueError, 'Feasibility slope'):
                SlopeProjection().solve(mesh, fixed, datum=0., ceiling=.08,
                    maximum_slope=.049, feasibility_slope=bound)

    def test_disabled_projection_cannot_be_reenabled_by_the_publication_guard(self):
        from tested_python_streets.publication_slopes import PublicationSlopeGuard
        mesh, fixed = peak_grid()
        mesh.xyz[:, :2] = mesh.xyz[:, :2] * .1 + (12345., -12345.)
        mesh.xyz[4, 2] = .01 / np.sqrt(2.)
        field = SimpleNamespace(datum=0., ceiling=.08, metadata={})
        self.assertLessEqual(float(mesh.quad_slopes().max()), .1 + 1e-10)
        rounded = PublicationSlopeGuard.physical_buffer(mesh, 0., 1.)
        self.assertGreater(float(rounded.quad_slopes().max()), .1 + 1e-7)
        with patch('tested_python_streets.publication_slopes.SlopeProjection.solve_steps',
                   side_effect=AssertionError('Disabled projection was called')):
            with self.assertRaisesRegex(ValueError, 'projection is disabled'):
                PublicationSlopeGuard().correct(mesh, field, fixed, exaggeration=1., limit=.1,
                    max_cells=100, check_joins=lambda candidate: None, allow_projection=False)

    def test_final_spill_surface_is_checked_in_physical_and_stored_coordinates(self):
        from tested_python_streets.publication_slopes import PublicationSlopeGuard
        mesh,fixed=peak_grid()
        mesh.xyz[:,:2]=mesh.xyz[:,:2]*.1+(12345.,-12345.)
        mesh.xyz[:,2]=.04
        mesh.xyz[4,2]-=.01/np.sqrt(2.)
        guard=PublicationSlopeGuard()
        guard.maximum_rounds=1
        field=SimpleNamespace(datum=0.,ceiling=.08,metadata={})
        self.assertGreater(guard.physical_buffer(mesh,0.,1.).quad_slopes().max(),.1+guard.tolerance)
        checked=[]
        result,_,report=guard.correct(mesh,field,fixed,exaggeration=1.,limit=.1,max_cells=100,
                                      check_joins=checked.append)
        self.assertIs(checked[0],result)
        self.assertTrue(report['rounds'][-1]['accepted_spill_surface'])
        self.assertEqual(report['rounds'][-1]['spill_surface_stored_drainage']['maximum_fill'],0.)
        np.testing.assert_array_equal(result.xyz[:,2],np.full(len(mesh.xyz),.04))
        np.testing.assert_array_equal(result.xyz[fixed],mesh.xyz[fixed])
        np.testing.assert_array_equal(result.xyz[:,:2],mesh.xyz[:,:2])

    def test_budget_limited_spill_preserves_directional_constraints(self):
        from tested_python_streets.publication_slopes import PublicationSlopeGuard
        from tested_python_streets.slope_projection import ProjectionIterationBudgetError
        from tested_python_streets.height_slabs import HeightSlabs
        mesh,fixed=peak_grid()
        mesh.xyz[:,2]=.04
        mesh.xyz[4,2]=.12
        points=mesh.xyz.copy()
        points[4,2]=.039
        failure=ProjectionIterationBudgetError('bounded projection',mesh.with_positions(points),{})
        field=SimpleNamespace(datum=0.,ceiling=.12,metadata={})
        guard=PublicationSlopeGuard()
        options=dict(exaggeration=1.,limit=.10,max_cells=100,check_joins=lambda m:None)
        with patch.object(SlopeProjection,'solve_steps',side_effect=failure):
            result,_,report=guard.correct(mesh,field,fixed,**options)
            self.assertTrue(report['rounds'][0]['accepted_spill_surface'])
            np.testing.assert_array_equal(result.xyz[:,2],np.full(len(mesh.xyz),.04))
            slabs=HeightSlabs(np.array(((4,0,1,2),)),np.array(((1.,0.,0.,0.),)),np.array((.0395,)))
            with self.assertRaises(ProjectionIterationBudgetError):
                guard.correct(mesh,field,fixed,height_slabs=slabs,**options)

    def test_publication_target_respects_the_immutable_incoming_crossfall(self):
        from tested_python_streets.publication_slopes import PublicationSlopeGuard
        pattern = PatternLibrary().grid(2, 2)
        angle = .04
        xy = 1000. + .01 * np.c_[2. * pattern.uv[:, 0] + 2. * np.cos(angle) * pattern.uv[:, 1],
                                  2. * np.sin(angle) * pattern.uv[:, 1]]
        normal = np.array((-np.sin(angle / 2.), np.cos(angle / 2.)))
        heights = (xy - xy[0]) @ normal * .0857142838
        heights -= heights.min()
        mesh = Mesh(np.c_[xy, heights], pattern.quads)
        mesh.xyz[4, 2] += 2.83695410088279e-6
        fixed = np.ones(len(mesh.xyz), bool); fixed[4] = False
        field = SimpleNamespace(datum=0., ceiling=float(mesh.xyz[:, 2].max()), metadata={})
        guard = PublicationSlopeGuard()
        self.assertLess(float(mesh.quad_slopes().max()), .1)
        self.assertGreater(float(guard.physical_buffer(mesh, 0., 1.).quad_slopes().max()), .1)
        before = mesh.xyz.copy()
        checked = []
        result, query, report = guard.correct(mesh, field, fixed, exaggeration=1., limit=.1,
            max_cells=100, check_joins=lambda candidate: checked.append(candidate))
        self.assertLessEqual(float(guard.physical_buffer(result, 0., 1.).quad_slopes().max()), .1 + guard.tolerance)
        np.testing.assert_array_equal(result.xyz[fixed], before[fixed])
        np.testing.assert_array_equal(result.xyz[:, :2], before[:, :2])
        np.testing.assert_array_equal(result.quads, mesh.quads)
        np.testing.assert_array_equal(mesh.xyz, before)
        self.assertEqual(len(checked), 1)
        self.assertLess(.1 - 4. * report['rounds'][0]['roundoff_gradient_error'], .0857142838)
        self.assertGreater(report['rounds'][0]['target'], .0857142838)
        self.assertGreaterEqual(float(result.xyz[:, 2].min()), field.datum)
        self.assertLessEqual(float(result.xyz[:, 2].max()), field.ceiling)
        np.testing.assert_allclose(query(result.xyz[:, :2]), result.xyz[:, 2], atol=1e-10, rtol=0.)
        # The established first successful target and its numerical result
        # remain exact; the new schedule is never entered on this fixture.
        self.assertEqual(len(report['rounds']),1)
        expected,_=SlopeProjection(max_cells=100).solve(mesh,fixed,datum=0.,ceiling=field.ceiling,
            maximum_slope=guard.fixed_gradient_floor(mesh,fixed,0.)+.75*(.1-guard.fixed_gradient_floor(mesh,fixed,0.)),
            minimum_heights=np.zeros(len(mesh.xyz)))
        np.testing.assert_array_equal(result.xyz,expected.xyz)

    def test_publication_headroom_progresses_beyond_a_real_stored_target_fixed_point(self):
        from tested_python_streets.publication_slopes import PublicationSlopeGuard
        pattern=PatternLibrary().grid(2,2);angle=.06
        xy=1000.+.00004+.01*np.c_[2.*pattern.uv[:,0]+2.*np.cos(angle)*pattern.uv[:,1],
                                     2.*np.sin(angle)*pattern.uv[:,1]]
        normal=np.array((-np.sin(angle/2.),np.cos(angle/2.)))
        heights=(xy-xy[0])@normal*.0857142838;heights-=heights.min()
        mesh=Mesh(np.c_[xy,heights],pattern.quads);mesh.xyz[4,2]+=.000003
        fixed=np.ones(len(xy),bool);fixed[4]=False
        field=SimpleNamespace(datum=0.,ceiling=float(mesh.xyz[:,2].max()),metadata={})
        guard=PublicationSlopeGuard();guard.maximum_rounds=1
        # The old first target is already physically feasible and leaves the
        # same mesh, but that mesh fails after actual float32 coordinate storage.
        floor=guard.fixed_gradient_floor(mesh,fixed,0.)
        first_target=floor+.75*(.1-floor)
        first,details=SlopeProjection(max_cells=100).solve(mesh,fixed,datum=0.,ceiling=field.ceiling,
            maximum_slope=first_target,minimum_heights=np.zeros(len(xy)))
        np.testing.assert_array_equal(first.xyz,mesh.xyz)
        self.assertGreater(guard.physical_buffer(first,0.,1.).quad_slopes().max(),.1+guard.tolerance)
        with self.assertRaisesRegex(ValueError,'iteration budget'):
            guard.correct(mesh,field,fixed,exaggeration=1.,limit=.1,max_cells=100,check_joins=lambda m:None)
        guard.maximum_rounds=8;checked=[]
        result,query,report=guard.correct(mesh,field,fixed,exaggeration=1.,limit=.1,max_cells=100,
            check_joins=lambda m:checked.append(m))
        targets=np.array([row['target'] for row in report['rounds']])
        self.assertEqual(targets[0],first_target)
        self.assertGreater(len(targets),1);self.assertLessEqual(len(targets),8)
        self.assertTrue(np.all(np.diff(targets)<0.));self.assertTrue(np.all(targets>floor))
        self.assertLessEqual(guard.physical_buffer(result,0.,1.).quad_slopes().max(),.1+guard.tolerance)
        self.assertEqual(len(checked),1)
        np.testing.assert_array_equal(result.xyz[fixed],mesh.xyz[fixed])
        np.testing.assert_array_equal(result.xyz[:,:2],mesh.xyz[:,:2])
        np.testing.assert_array_equal(result.quads,mesh.quads)
        np.testing.assert_allclose(query(result.xyz[:,:2]),result.xyz[:,2],atol=1e-10,rtol=0.)
        with self.assertRaisesRegex(ValueError,'Fixed boundary gradients'):
            guard.correct(mesh,field,np.ones(len(fixed),bool),exaggeration=1.,limit=.1,
                max_cells=100,check_joins=lambda m:None)

    def test_fixed_gradient_floor_distinguishes_free_and_immutable_directions(self):
        from tested_python_streets.publication_slopes import PublicationSlopeGuard
        mesh, fixed = peak_grid(12.5)
        mesh.xyz[:, 2] = 12.5 + .08 * mesh.xyz[:, 0]
        mesh.xyz[4, 2] += .03
        guard = PublicationSlopeGuard()
        self.assertAlmostEqual(guard.fixed_gradient_floor(mesh, fixed, 12.5), .08, places=12)
        self.assertEqual(guard.fixed_gradient_floor(mesh, np.zeros(len(fixed), bool), 12.5), 0.)
        self.assertAlmostEqual(guard.fixed_gradient_floor(mesh, np.ones(len(fixed), bool), 12.5),
                               float(mesh.quad_slopes().max()), places=12)

    def test_publication_can_certify_a_bounded_last_iterate_without_accepting_a_failed_solve(self):
        from tested_python_streets.publication_slopes import PublicationSlopeGuard
        from tested_python_streets.slope_projection import ProjectionIterationBudgetError
        mesh,fixed=peak_grid()
        mesh.xyz[:,:2]=mesh.xyz[:,:2]*.1+(12345.,-12345.)
        mesh.xyz[4,2]=.01/np.sqrt(2.)
        field=SimpleNamespace(datum=0.,ceiling=.08,metadata={})
        solver=SlopeProjection(max_iterations=1)
        with self.assertRaises(ProjectionIterationBudgetError):
            solver.solve(mesh,fixed,datum=0.,ceiling=.08,maximum_slope=.049)
        guard=PublicationSlopeGuard();before=mesh.xyz.copy();checked=[]
        with patch('tested_python_streets.publication_slopes.SlopeProjection',return_value=solver):
            result,query,report=guard.correct(mesh,field,fixed,exaggeration=1.,limit=.1,
                max_cells=100,check_joins=lambda candidate:checked.append(candidate))
        self.assertEqual(report['rounds'][0]['projection']['convergence'],
                         'independently_qualified_published_candidate')
        self.assertLessEqual(float(guard.physical_buffer(result,0.,1.).quad_slopes().max()),.1+guard.tolerance)
        self.assertEqual(len(checked),1)
        np.testing.assert_array_equal(result.xyz[fixed],before[fixed])
        np.testing.assert_array_equal(result.xyz[:,:2],before[:,:2])
        np.testing.assert_array_equal(mesh.xyz,before)
        np.testing.assert_allclose(query(result.xyz[:,:2]),result.xyz[:,2],atol=1e-10,rtol=0.)
        with patch('tested_python_streets.publication_slopes.SlopeProjection',return_value=SlopeProjection(max_iterations=0)):
            with self.assertRaises(ProjectionIterationBudgetError):
                guard.correct(mesh,field,fixed,exaggeration=1.,limit=.1,max_cells=100,
                    check_joins=lambda candidate:None)

    def test_coupled_solver_still_checks_drainage_and_joins_after_feasibility_stop(self):
        from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
        from tested_python_streets.slope_projection import ProjectionIterationBudgetError
        mesh,fixed=peak_grid()
        field=QuadHeightField(mesh,datum=0.,ceiling=.08)
        refs={int(i):None for i in np.flatnonzero(fixed)}
        builder=CompoundSurfaceBuilder(PatternLibrary())
        with patch('tested_python_streets.compound_surfaces.SlopeProjection',return_value=SlopeProjection(max_iterations=1)):
            result,query,report=builder.correct_refined_grading(mesh,field,refs,(),.05,100)
        self.assertEqual(report['rounds'][0]['projection']['convergence'],
                         'physical_feasibility_after_iteration_budget')
        self.assertEqual(report['joins']['isolated_interior_minima'],0)
        self.assertLessEqual(report['rounds'][-1]['drainage']['maximum_fill'],1e-6)
        self.assertLessEqual(float(result.quad_slopes().max()),.05+1e-8)
        np.testing.assert_array_equal(result.xyz[fixed],mesh.xyz[fixed])
        np.testing.assert_allclose(query(result.xyz[:,:2]),result.xyz[:,2],atol=1e-10,rtol=0.)
        with patch('tested_python_streets.compound_surfaces.SlopeProjection',return_value=SlopeProjection(max_iterations=0)):
            with self.assertRaises(ProjectionIterationBudgetError):
                builder.correct_refined_grading(mesh,field,refs,(),.05,100)

    def test_stalled_spill_does_not_get_additional_projections(self):
        from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
        from tested_python_streets.drainage import BoundaryDrainage
        mesh, fixed = peak_grid()
        mesh.xyz[:, 2] = 0.
        field = QuadHeightField(mesh, datum=0., ceiling=.08)
        refs = dict.fromkeys(np.flatnonzero(fixed))

        def stationary_residual(candidate, mask):
            heights = candidate.xyz[:, 2].copy()
            heights[~mask] += .00001
            return heights, dict(maximum_fill=.00001, raised_vertices=1)

        with patch.object(BoundaryDrainage, 'fill', side_effect=stationary_residual) as fill:
            with self.assertRaisesRegex(ValueError, 'drainage exceeded its iteration budget'):
                CompoundSurfaceBuilder(PatternLibrary()).correct_refined_grading(mesh, field, refs, (), .05, 100)
        self.assertEqual(fill.call_count, 9)  # Initial measurement plus eight projections.
        np.testing.assert_array_equal(mesh.xyz[:, 2], 0.)

    def test_explicit_publication_directional_slabs_keep_profiles_and_default_search(self):
        from tested_python_streets.publication_slopes import PublicationSlopeGuard
        from tested_python_streets.height_slabs import HeightSlabs
        mesh,fixed=peak_grid()
        mesh.xyz[:,:2]=mesh.xyz[:,:2]*.1+(12345.,-12345.)
        mesh.xyz[4,2]=.01/np.sqrt(2.)
        field=SimpleNamespace(datum=0.,ceiling=.08,metadata={})
        matrices=corner_matrices(mesh)
        rows=mesh.quads.copy()
        coefficients=matrices[:,0,1,:].copy()
        slabs=HeightSlabs(rows,coefficients,np.full(len(rows),.03))
        guard=PublicationSlopeGuard()
        ordinary,_,_=guard.correct(mesh,field,fixed,exaggeration=1.,limit=.1,max_cells=100,
            check_joins=lambda candidate:None)
        self.assertGreater(slabs.violation(ordinary.xyz[:,2]),.001)
        checked=[]
        result,query,report=guard.correct(mesh,field,fixed,exaggeration=1.,limit=.1,max_cells=100,
            check_joins=lambda candidate:checked.append(candidate),height_slabs=slabs)
        self.assertLessEqual(slabs.violation(result.xyz[:,2]),1e-9)
        self.assertLessEqual(float(guard.physical_buffer(result,0.,1.).quad_slopes().max()),.1+guard.tolerance)
        self.assertEqual(report['directional_height_slabs']['rows'],len(rows))
        self.assertEqual(len(checked),1)
        np.testing.assert_array_equal(result.xyz[fixed],mesh.xyz[fixed])
        np.testing.assert_array_equal(result.xyz[:,:2],mesh.xyz[:,:2])
        np.testing.assert_array_equal(result.quads,mesh.quads)
        np.testing.assert_allclose(query(result.xyz[:,:2]),result.xyz[:,2],atol=1e-10,rtol=0.)

    def test_one_free_peak_has_the_analytical_minimum_displacement_solution(self):
        original, fixed = peak_grid()
        expected = .05 / np.sqrt(2.)
        result, report = SlopeProjection().solve(original, fixed, datum=0., ceiling=.08, maximum_slope=.05)
        self.assertAlmostEqual(result.xyz[4, 2], expected, places=10)
        np.testing.assert_array_equal(result.xyz[fixed], original.xyz[fixed])
        np.testing.assert_array_equal(result.xyz[:, :2], original.xyz[:, :2])
        np.testing.assert_array_equal(result.quads, original.quads)
        self.assertEqual(original.xyz[4, 2], .08)
        self.assertLessEqual(result.quad_slopes().max(), .05 + 1e-9)
        self.assertLess(report['iterations'], 1000)

    def test_corner_constraints_match_the_independent_nine_point_slope_audit(self):
        mesh, fixed = peak_grid()
        mesh.xyz[4, :2] = .8, .9
        matrices = corner_matrices(mesh)
        values = np.linalg.norm(np.einsum('fcij,fj->fci', matrices, mesh.xyz[mesh.quads, 2]), axis=-1).max(axis=1)
        np.testing.assert_allclose(values, mesh.quad_slopes(), atol=1e-14, rtol=0.)

    def test_repeated_runs_and_common_datum_translation_are_deterministic(self):
        mesh, fixed = peak_grid()
        solver = SlopeProjection()
        a, _ = solver.solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.05)
        b, _ = solver.solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.05)
        np.testing.assert_array_equal(a.xyz, b.xyz)
        raised, fixed = peak_grid(12.5)
        c, _ = solver.solve(raised, fixed, datum=12.5, ceiling=12.58, maximum_slope=.05)
        np.testing.assert_allclose(a.xyz[:, 2], c.xyz[:, 2] - 12.5, atol=1e-10, rtol=0.)

    def test_feasible_data_retains_its_buffers_and_bad_fixed_data_is_rejected(self):
        mesh, fixed = peak_grid()
        unchanged, report = SlopeProjection(max_iterations=0).solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.2)
        self.assertIs(unchanged, mesh)
        self.assertEqual(report['iterations'], 0)
        with self.assertRaisesRegex(ValueError, 'infeasible'):
            SlopeProjection().solve(mesh, np.ones(len(fixed), bool), datum=0., ceiling=.08, maximum_slope=.05)
        with self.assertRaisesRegex(ValueError, 'iteration budget'):
            SlopeProjection(max_iterations=0).solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.05)

    def test_projection_budgets_and_boundary_envelopes_are_checked(self):
        for kwargs in ({'max_cells': 0}, {'max_iterations': -1}, {'memory_bytes': 0}, {'max_colors': 0}):
            with self.assertRaises(ValueError):
                SlopeProjection(**kwargs)
        mesh, fixed = peak_grid()
        with self.assertRaisesRegex(ValueError, 'memory budget'):
            SlopeProjection(memory_bytes=4096).solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.05)
        with self.assertRaisesRegex(ValueError, 'Boolean'):
            SlopeProjection().solve(mesh, np.flatnonzero(fixed), datum=0., ceiling=.08, maximum_slope=.05)
        mesh.xyz[0, 2] = -.01
        with self.assertRaisesRegex(ValueError, 'crown envelope'):
            SlopeProjection().solve(mesh, fixed, datum=0., ceiling=.08, maximum_slope=.05)

    def test_a_boundary_gradient_equal_to_the_limit_is_not_mistaken_for_infeasibility(self):
        mesh, fixed = peak_grid()
        mesh.xyz[:, 2] = .05 * mesh.xyz[:, 0]
        mesh.xyz[4, 2] += .08
        result, _ = SlopeProjection().solve(mesh, fixed, datum=0., ceiling=.2, maximum_slope=.05)
        np.testing.assert_allclose(result.xyz[:, 2], .05 * result.xyz[:, 0], atol=1e-10, rtol=0.)


class PythonAdjustedSurfaceQueryTests(unittest.TestCase):
    def test_standalone_junction_anchor_queries_survive_publication_buffer_release(self):
        from _python_street_fixtures import fixture_networks
        from tested_python_streets import model
        from tested_python_streets.backend import PythonMeshBackend
        network = fixture_networks(model)[0]
        network = replace(network, furniture=(model.Placement(81, 0, 1., offset=-3.),))
        backend = PythonMeshBackend(network.settings)
        prepared = backend.prepare(network)
        before = prepared.products[model.Key('placement', 81)]
        self.assertGreater(before.position[2], .1)
        backend.commit(prepared)
        backend.release_published_meshes()
        changed = replace(network, furniture=(replace(network.furniture[0], station=1.2),))
        edited = backend.prepare(changed)
        fresh = PythonMeshBackend(changed.settings).prepare(changed)
        self.assertEqual(edited.rebuilt, (model.Key('placement', 81),))
        self.assertEqual(edited.products[model.Key('placement', 81)], fresh.products[model.Key('placement', 81)])

    def test_inverse_queries_match_bilinear_heights_on_varied_convex_quads(self):
        random = np.random.default_rng(731)
        for _ in range(12):
            xy = np.array(((0., 0.), (3., 0.), (3., 2.), (0., 2.))) + random.uniform(-.2, .2, (4, 2))
            height = random.uniform(0., .08, 4)
            mesh = Mesh(np.column_stack((xy, height)), np.array(((0, 1, 2, 3),)))
            field = QuadHeightField(mesh, datum=0., ceiling=.08)
            uv = np.vstack((random.uniform(0., 1., (40, 2)), ((0., 0.), (1., 0.), (1., 1.), (0., 1.))))
            u, v = uv.T
            weights = np.column_stack(((1. - u) * (1. - v), u * (1. - v), u * v, (1. - u) * v))
            np.testing.assert_allclose(field(weights @ xy), weights @ height, atol=1e-11, rtol=0.)
            offset = np.array((1_000_000., -2_000_000.))
            moved = QuadHeightField(Mesh(np.column_stack((xy + offset, height)), mesh.quads), datum=0., ceiling=.08)
            np.testing.assert_allclose(moved(weights @ xy + offset), weights @ height, atol=1e-9, rtol=0.)

    def test_queries_use_adjusted_heights_and_survive_release_and_exaggeration(self):
        original, fixed = peak_grid()
        adjusted, _ = SlopeProjection().solve(original, fixed, datum=0., ceiling=.08, maximum_slope=.05)
        field = QuadHeightField(adjusted, datum=0., ceiling=.08)
        expected = .05 / np.sqrt(2.)
        np.testing.assert_allclose(field(np.array(((1., 1.), (.75, .25)))), (expected, .75 * .25 * expected), atol=1e-10)
        reference = weakref.ref(adjusted)
        adjusted.xyz[:, 2] *= 10.
        del adjusted
        gc.collect()
        self.assertIsNone(reference())
        self.assertAlmostEqual(field(np.array(((1., 1.),)))[0], expected)

    def test_owner_excludes_separator_holes_before_strict_pavement_queries(self):
        pattern = PatternLibrary().grid(3, 3)
        mesh = Mesh(np.column_stack((pattern.uv, np.full(len(pattern.uv), .08))), np.delete(pattern.quads, 4, axis=0))
        field = QuadHeightField(mesh, datum=0., ceiling=.08)
        loops = (SimpleNamespace(hole=False, xy=np.array(((0., 0.), (1., 0.), (1., 1.), (0., 1.)))),
                 SimpleNamespace(hole=True, xy=np.array(((1/3, 1/3), (1/3, 2/3), (2/3, 2/3), (2/3, 1/3)))))
        from tested_python_streets.hole_surfaces import ConstantHoleHeight
        owner = CompoundProduct(None, SimpleNamespace(loops=loops), field, {}, {}, {}, {}, None, None, 10.,
                                hole_fields={1:ConstantHoleHeight(.08,0.,.08)})
        np.testing.assert_allclose(owner.height(np.array(((.1, .1), (.5, .5), (2., 2.)))), (.8, .8, 0.), atol=1e-14)
        with self.assertRaisesRegex(ValueError, 'outside pavement'):
            field(np.array(((.5, .5),)))

    def test_invalid_queries_and_spatial_budgets_are_explicit(self):
        mesh, fixed = peak_grid()
        field = QuadHeightField(mesh, datum=0., ceiling=.08)
        self.assertEqual(field(np.empty((0, 2))).shape, (0,))
        for query in (np.ones((2, 3)), np.array(((float('nan'), 0.),))):
            with self.assertRaises(ValueError):
                field(query)
        with self.assertRaisesRegex(ValueError, 'outside pavement'):
            field(np.array(((4., 4.),)))
        with self.assertRaisesRegex(ValueError, 'spatial/memory budget'):
            QuadHeightField(mesh, datum=0., ceiling=.08, max_entries=1)


if __name__ == '__main__':
    unittest.main()
