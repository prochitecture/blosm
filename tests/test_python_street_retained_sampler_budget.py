"""Retained analytical circles keep their actual prepared source budget."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_retained_arc_fillets import primitive_fixture, prepared_view
from tested_python_streets.sampling import CertifiedCurbSampler


class RetainedSamplerBudgetTests(unittest.TestCase):
    def record(self):
        context, group = primitive_fixture()
        _, _, view = prepared_view(group)
        return context, view.physical_layout.retained_arc_fillets[0]

    def test_coarser_request_keeps_every_source_knot_and_prepared_error_bound(self):
        for caller_tolerance in (.005, .00025, .00001):
            context, record = self.record()
            before = context.network.to_dict()
            sampler = CertifiedCurbSampler(caller_tolerance, .7, 4096)
            at = sampler.sample(record, 0., 1.)
            self.assertTrue(all(value in at for value in record.source_reference_parameters))
            error = record.radius * (1. - np.cos(.5 * abs(record.sweep) * np.diff(at)))
            self.assertLessEqual(float(error.max()), min(caller_tolerance, record.source_chord_tolerance))
            self.assertEqual((at[0], at[-1]), (0., 1.))
            self.assertEqual(context.network.to_dict(), before)

    def test_already_aligned_native_sampler_keeps_its_previous_output_exact(self):
        _, record = self.record()
        sampler = CertifiedCurbSampler(record.source_chord_tolerance, .2, 4096)
        step = min(sampler.spacing, np.sqrt(8. * sampler.tolerance * record.radius))
        old = sampler.regular(0., 1., step / record.length)
        aligned = replace(record, source_reference_parameters=tuple(old))
        np.testing.assert_array_equal(aligned.certified_parameters(sampler, 0., 1.), old)
        self.assertEqual(aligned.radius, record.radius)
        self.assertIs(aligned.parent_arcs, record.parent_arcs)

    def test_partial_and_reverse_queries_keep_exact_finite_bounds(self):
        _, record = self.record()
        sampler = CertifiedCurbSampler(.005, .7, 4096)
        forward = sampler.sample(record, .17, .83)
        reverse = sampler.sample(record, .83, .17)
        self.assertEqual((forward[0], forward[-1]), (.17, .83))
        self.assertTrue(all(value in forward for value in record.source_reference_parameters if .17 < value < .83))
        self.assertTrue(np.all((forward >= .17) & (forward <= .83)))
        np.testing.assert_array_equal(reverse, forward[::-1])

    def test_compulsory_knots_cannot_escape_the_existing_sample_allocation(self):
        _, record = self.record()
        sampler = CertifiedCurbSampler(.005, .7, 4)
        with self.assertRaisesRegex(ValueError, 'allocation|budget'):
            sampler.sample(record, 0., 1.)


if __name__ == '__main__':
    unittest.main()
