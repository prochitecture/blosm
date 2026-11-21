"""Failed primary charts get one candidate; original source errors stay fatal."""
from copy import copy
from dataclasses import replace
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.partition import CompoundQuadChart,QuadChartGradingError
from tested_python_streets.quad_charts import QuadChartCellBudgetError
from tested_python_streets.physical_hole_cuts import PhysicalHoleCutError
from tested_python_streets.source_strip_containment import SourceStripContainmentError
from tested_python_streets.source_exterior_sampling import (
    chart_sampling_scope,try_primary_sampling,eligible)
from tested_python_streets.model import Settings
from tested_python_streets.patterns import PatternLibrary
from test_python_street_source_rim_normalization import fixture as positive_exterior
from test_python_street_source_exterior_sampling import fixture as generated_contact


def protocol_fixture():
    # Source/phase geometry is genuine. These tests mock only the expensive
    # width/chart qualification to isolate ordering and error propagation;
    # they do not claim this local corner is an accepted complete compound.
    group,_=generated_contact()
    context=SimpleNamespace(network=SimpleNamespace(settings=Settings()))
    normalized=replace(group,physical_layout=copy(group.physical_layout))
    normalized.physical_layout.source_rim_normalization=({'maximum_projection_error':.5e-6},)
    normalized.physical_layout.source_exterior_sampling={'original_holes_exact':True}
    return group,context,normalized


class PrimarySamplingRecoveryTests(unittest.TestCase):
    def test_progressive_charts_own_their_retry_state_across_yields_and_cancellation(self):
        from tested_python_streets.source_exterior_sampling import chart_sampling_steps, _construction, _pool_key
        group, context = SimpleNamespace(), object()
        closed = []

        @chart_sampling_steps
        def chart(self, group, context, name):
            state = _construction.get()
            state['declined'][_pool_key(group, context)] = name
            try:
                yield state
                self.assertIs(_construction.get(), state)
                return SimpleNamespace(metadata={})
            finally:
                closed.append(_construction.get() is state)

        first, second = (chart(self, group, context, name) for name in ('first','second'))
        self.assertIsNot(next(first), next(second))
        self.assertIsNone(_construction.get())
        first.close()
        self.assertIsNone(_construction.get())
        with self.assertRaises(StopIteration) as done:
            next(second)
        self.assertEqual(done.exception.value.metadata['primary_generated_phase_sampling_declined'], 'second')
        self.assertEqual(closed, [True, True])
        self.assertIsNone(_construction.get())

    def test_positive_primary_retains_exact_geometry_and_never_enters_sampling(self):
        group=positive_exterior()
        context=SimpleNamespace(network=SimpleNamespace(settings=Settings()),patterns=PatternLibrary())
        chart=CompoundQuadChart(context.patterns)
        expected=chart._build(group,context,require_quality=True)
        with patch('tested_python_streets.partition.try_primary_sampling',
                side_effect=AssertionError('successful primary entered sampling')):
            result=chart.build(group,context)
        np.testing.assert_array_equal(result.mesh.xyz,expected.mesh.xyz)
        np.testing.assert_array_equal(result.mesh.quads,expected.mesh.quads)

    def test_primary_budget_recovery_returns_before_secondary_construction(self):
        group=positive_exterior()
        context=SimpleNamespace(network=SimpleNamespace(settings=Settings()),patterns=PatternLibrary())
        chart=CompoundQuadChart(context.patterns)
        sentinel=SimpleNamespace(metadata={})
        with patch.object(chart,'_build_steps',side_effect=QuadChartCellBudgetError('primary cell budget')):
            with patch('tested_python_streets.partition.try_primary_sampling',return_value=sentinel) as recover:
                with patch.object(chart,'_reflex',side_effect=AssertionError('secondary work ran')):
                    self.assertIs(chart.build(group,context),sentinel)
        self.assertEqual(recover.call_count,1)
        self.assertEqual(recover.call_args.args[-1],'primary chart budget')

    def test_ineligible_and_unrelated_failures_do_not_enter_source_or_candidate_work(self):
        group,context,_=protocol_fixture()
        with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',
                side_effect=AssertionError('unexpected width proof')):
            self.assertIsNone(try_primary_sampling(None,replace(group,claims=()),context,
                QuadChartGradingError('primary grading'),'primary'))
            self.assertIsNone(try_primary_sampling(None,group,context,ValueError('unrelated source error'),'primary'))

    def test_successful_candidate_uses_existing_proofs_and_avoids_raw_secondary(self):
        group,context,normalized=protocol_fixture()
        candidate=SimpleNamespace(metadata={},partition=SimpleNamespace(group=normalized))
        calls=[];parent=SimpleNamespace(max_cells=50000,build=lambda g,c:(calls.append(g) or candidate))
        @chart_sampling_scope
        def construction(self,g,c):
            result=try_primary_sampling(self,g,c,QuadChartGradingError('primary grading'),'first chart')
            if result is None:raise AssertionError('raw secondary ran after success')
            return result
        with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',return_value={'maximum_error':.002}) as proof:
            with patch('tested_python_streets.source_rim_normalization.normalize_exterior',return_value=normalized):
                result=construction(parent,group,context)
        self.assertEqual(calls,[normalized]);self.assertEqual(proof.call_count,3)
        self.assertEqual(result.metadata['primary_generated_phase_sampling_recovery']['maximum_candidate_rebuilds'],1)

    def test_declined_candidate_preserves_original_fallback_and_nested_retry_budget(self):
        group,context,normalized=protocol_fixture();calls=[]
        failure=QuadChartGradingError('normalized candidate cannot grade')
        def candidate(g,c):calls.append(g);raise failure
        parent=SimpleNamespace(max_cells=50000,build=candidate)
        fallback=SimpleNamespace(metadata={'original_raw_fallback':True})
        @chart_sampling_scope
        def nested(self,g,c):
            self_result=try_primary_sampling(self,g,c,QuadChartGradingError('later piecewise failure'),'nested')
            if self_result is not None:raise AssertionError('candidate repeated')
            return fallback
        @chart_sampling_scope
        def construction(self,g,c):
            result=try_primary_sampling(self,g,c,QuadChartGradingError('primary grading'),'first chart')
            if result is not None:raise AssertionError('failed candidate accepted')
            if eligible(g,c):raise AssertionError('attempt budget was not consumed')
            return nested(self,g,c)
        with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',return_value={'maximum_error':.002}):
            with patch('tested_python_streets.source_rim_normalization.normalize_exterior',return_value=normalized):
                result=construction(parent,group,context)
        self.assertIs(result,fallback);self.assertEqual(len(calls),1)
        self.assertEqual(result.metadata['primary_generated_phase_sampling_declined']['failure_type'],'QuadChartGradingError')
        self.assertTrue(eligible(group,context))  # independent constructions get a fresh budget

    def test_invalid_original_proof_is_fatal_before_any_candidate_or_fallback(self):
        group,context,normalized=protocol_fixture();source_error=SourceStripContainmentError('original is clipped')
        parent=SimpleNamespace(max_cells=50000,build=lambda g,c:(_ for _ in ()).throw(AssertionError('candidate ran')))
        @chart_sampling_scope
        def construction(self,g,c):
            result=try_primary_sampling(self,g,c,QuadChartGradingError('primary grading'),'first chart')
            raise AssertionError('invalid original fell through to secondary')
        with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',side_effect=source_error):
            with patch('tested_python_streets.source_rim_normalization.normalize_exterior') as normalize:
                with self.assertRaises(SourceStripContainmentError) as caught:construction(parent,group,context)
        self.assertIs(caught.exception,source_error);normalize.assert_not_called()

    def test_only_typed_candidate_clipping_can_decline_after_the_original_proof(self):
        group,context,normalized=protocol_fixture()
        parent=SimpleNamespace(max_cells=50000,build=lambda g,c:(_ for _ in ()).throw(AssertionError('chart ran')))
        @chart_sampling_scope
        def construction(self,g,c):
            return try_primary_sampling(self,g,c,QuadChartGradingError('primary grading'),'first chart')
        with patch('tested_python_streets.source_rim_normalization.normalize_exterior',return_value=normalized):
            with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',
                    side_effect=[{'maximum_error':.002},SourceStripContainmentError('candidate clips source')]):
                self.assertIsNone(construction(parent,group,context))
            budget=ValueError('original offset curvature certificate unavailable')
            with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',
                    side_effect=[{'maximum_error':.002},budget]):
                with self.assertRaises(ValueError) as caught:construction(parent,group,context)
            self.assertIs(caught.exception,budget)

    def test_unrelated_candidate_chart_error_is_not_a_typed_decline(self):
        group,context,normalized=protocol_fixture();error=ValueError('unrelated chart implementation failure')
        def broken(g,c):raise error
        parent=SimpleNamespace(max_cells=50000,build=broken)
        @chart_sampling_scope
        def construction(self,g,c):
            return try_primary_sampling(self,g,c,QuadChartGradingError('primary grading'),'first chart')
        with patch('tested_python_streets.source_rim_normalization.normalize_exterior',return_value=normalized):
            with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',return_value={'maximum_error':.002}):
                with self.assertRaises(ValueError) as caught:construction(parent,group,context)
        self.assertIs(caught.exception,error)

    def test_candidate_hole_partition_exhaustion_retains_original_chart_families(self):
        group,context,normalized=protocol_fixture();calls=[]
        failure=PhysicalHoleCutError('candidate has no compatible finite hole cuts')
        def candidate(g,c):calls.append(g);raise failure
        parent=SimpleNamespace(max_cells=50000,build=candidate)
        fallback=SimpleNamespace(metadata={'original_raw_fallback':True})
        @chart_sampling_scope
        def construction(self,g,c):
            result=try_primary_sampling(self,g,c,QuadChartGradingError('primary grading'),'first chart')
            self_outer.assertIsNone(result)
            return fallback
        self_outer=self
        with patch('tested_python_streets.source_rim_normalization.normalize_exterior',return_value=normalized):
            with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',
                    return_value={'maximum_error':.002}) as proof:
                result=construction(parent,group,context)
        self.assertIs(result,fallback);self.assertEqual(calls,[normalized])
        self.assertEqual(proof.call_count,2)
        decline=result.metadata['primary_generated_phase_sampling_declined']
        self.assertEqual(decline['failure_type'],'PhysicalHoleCutError')
        self.assertEqual(decline['stage'],'candidate chart')

    def test_original_source_error_stays_fatal_after_a_typed_hole_cut_failure(self):
        group,context,_=protocol_fixture()
        source_error=SourceStripContainmentError('original full-width proof failed')
        parent=SimpleNamespace(max_cells=50000)
        @chart_sampling_scope
        def construction(self,g,c):
            return try_primary_sampling(self,g,c,PhysicalHoleCutError('primary finite cut family'),'first chart')
        with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',side_effect=source_error):
            with patch('tested_python_streets.source_rim_normalization.normalize_exterior') as normalize:
                with self.assertRaises(SourceStripContainmentError) as caught:
                    construction(parent,group,context)
        self.assertIs(caught.exception,source_error);normalize.assert_not_called()


if __name__=='__main__':unittest.main()
