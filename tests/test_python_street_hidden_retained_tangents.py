"""Finite analytical insertion must prove every omitted tangent and operand."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_retained_arc_fillets import primitive_fixture, prepared_view
from tested_python_streets.compounds import SampledFootprint
from tested_python_streets.model import Key
from tested_python_streets.regions import Region
from tested_python_streets.retained_arc_fillets import RetainedArcReconciliationError
from tested_python_streets.retained_tangent_visibility import restore_hidden
from tested_python_streets.sampling import AnalyticalBoundarySampler


def checker():
    sampler = AnalyticalBoundarySampler(tolerance=.005, spacing=.7)
    sampler.retained_arc_contact_tokens = set()
    sampler.retained_arc_contact_handoffs = []
    return sampler


class HiddenRetainedTangentTests(unittest.TestCase):
    def test_exact_parent_insertion_preserves_all_other_incidence_under_rigid_motion(self):
        for reverse, angle, translation in ((False, 0., (0., 0.)),
                                            (True, .713, (1200., -1800.))):
            context, group = primitive_fixture(second_origin=(-.5, .5), reverse=reverse,
                                               angle=angle, translation=translation)
            fillet, _, view = prepared_view(group)
            # This primitive fixture retains the original two-parent union;
            # the proposed wedge is not yet an area-union qualification.
            before = tuple((loop.xy.copy(), loop.spans) for loop in view.loops)
            sampler = checker()
            result, proof = restore_hidden(view, fillet, sampler, context=context)
            self.assertIs(result.footprints, view.footprints)
            self.assertIs(result.physical_layout, view.physical_layout)
            self.assertLessEqual(proof['proof_points'], 4096)
            self.assertEqual(len(sampler.retained_arc_contact_tokens), 2)
            spans = [span for loop in result.loops for span in loop.spans
                     if span.owner == fillet.owner]
            self.assertEqual(len(spans), 1)
            self.assertEqual((spans[0].label.first, spans[0].label.last), (0., 1.))
            for loop, (points, old_spans) in zip(view.loops, before):
                np.testing.assert_array_equal(loop.xy, points)
                self.assertEqual(loop.spans, old_spans)
            with self.assertRaisesRegex(RetainedArcReconciliationError, 'already visible'):
                restore_hidden(result, fillet, checker(), context=context)

    def test_genuinely_covering_third_operand_cannot_be_replaced_by_tangent_insertion(self):
        context, group = primitive_fixture(second_origin=(-.5, .5))
        fillet, _, view = prepared_view(group)
        points = fillet(np.linspace(0., 1., 101))
        lo, hi = points.min(0)-.01, points.max(0)+.01
        polygon = np.array((lo, (hi[0], lo[1]), hi, (lo[0], hi[1])))
        covering = SampledFootprint(Region(Key('unrelated_source', 99), polygon), (),
                                    SimpleNamespace(boundaries=()))
        changed = replace(view, footprints=(*view.footprints, covering))
        with self.assertRaisesRegex(RetainedArcReconciliationError, 'complete original native strip'):
            restore_hidden(changed, fillet, checker(), context=context)

    def test_invalid_consumption_or_parent_identity_is_fatal(self):
        context, group = primitive_fixture(second_origin=(-.5, .5))
        fillet, _, view = prepared_view(group)
        parent = fillet.parent_arcs[0]
        for interval in ((np.nan, parent.consumed_circle_parameters[1]),
                         (-1., parent.consumed_circle_parameters[1]),
                         (parent.consumed_circle_parameters[1], parent.consumed_circle_parameters[0])):
            wrong = replace(fillet, parent_arcs=(replace(parent, consumed_circle_parameters=interval),
                                                 fillet.parent_arcs[1]))
            with self.assertRaisesRegex(RetainedArcReconciliationError, 'finite original circular interval'):
                restore_hidden(view, wrong, checker(), context=context)
        wrong = replace(fillet, parent_arcs=(replace(parent, original_curve=object()),
                                             fillet.parent_arcs[1]))
        with self.assertRaisesRegex(RetainedArcReconciliationError, 'exact original parent ancestry'):
            restore_hidden(view, wrong, checker(), context=context)

    def test_proof_refinement_allocation_is_finite_and_never_skipped(self):
        context, group = primitive_fixture(second_origin=(-.5, .5))
        fillet, _, view = prepared_view(group)
        with self.assertRaisesRegex(RetainedArcReconciliationError, 'proof-refinement budget'):
            restore_hidden(view, fillet, checker(), context=context,
                           budget=dict(maximum=64, used=0))
        for budget in (dict(maximum=4097, used=0), dict(maximum=64, used=65),
                       dict(maximum=True, used=0), dict(maximum=64, used=-1)):
            with self.assertRaisesRegex(RetainedArcReconciliationError, 'Invalid bounded'):
                restore_hidden(view, fillet, checker(), context=context, budget=budget)


if __name__ == '__main__':
    unittest.main()
