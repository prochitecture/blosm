"""Jointly prepare original curb phases before either physical mouth qualifies."""
from dataclasses import replace
from copy import copy
import unittest
import numpy as np
import _python_street_package
from test_python_street_one_sided_paired_junctions import fixture
from test_python_street_compounds import footprint_context
from tested_python_streets import model
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.one_sided_paired_junctions import OneSidedPairedSpecification
from tested_python_streets.paired_source_sampling import (OneSidedPreparationRequest,
    canonical_reprepare, shared_source_parameters, verify_shared_source_phase, sample_source_parameters)


def two_sites():
    context, site, _, visits, _ = fixture(3)
    sections = tuple(replace(s, points=((-123.456, -4.2, 0.), (0., -4.1, 0.),
                        (30., -4.05, 0.), (145.67, -3.9, 0.))) if s.id == 3 else s
                     for s in context.network.sections)
    sections += (model.Section(5, 5, ((30., 4., 0.), (45., 4.1, 0.), (60., 4., 0.)), width=10.),
                 model.Section(6, 6, ((30., 4., 0.), (30., 30., 0.)), width=2.))
    remote = model.Intersection(101, (model.End(2, 1), model.End(5, 0), model.End(6, 0)), radius=2.)
    near = model.BundleMember((*visits, model.Traversal(5)))
    bundle = replace(context.network.bundles[0], members=(near, context.network.bundles[0].members[1]))
    network = replace(context.network, sections=sections, intersections=(site, remote), bundles=(bundle,))
    context = footprint_context(network)
    requests = tuple(OneSidedPreparationRequest(junction, IntersectionProvider.footprint(context, junction),
        pair, model.Traversal(3), 0) for junction, pair in ((site, visits),
        (remote, (model.Traversal(2), model.Traversal(5)))))
    return context, requests


class PairedSourceSamplingTests(unittest.TestCase):
    def test_shared_original_prefixes_are_reprepared_before_any_mouth_is_hydrated(self):
        context, requests = two_sites()
        source = context.network.to_dict()
        initial = tuple(OneSidedPairedSpecification.build(context, r.site, r.local_layout,
            r.near_visits, r.opposite_visit, r.near_role) for r in requests)
        old = tuple(tuple(a.source_curve_stations[1].copy() for a in spec.entrances) for spec in initial)
        with self.assertRaisesRegex(ValueError, 'canonical shared-prefix'):
            verify_shared_source_phase(tuple(a for spec in initial for a in spec.entrances))
        prepared, parameters, proof = canonical_reprepare(context, requests)
        self.assertTrue(proof['overlapping_station_arrays_exact'])
        self.assertTrue(any(row['section'] == 3 for row in proof['shared_intervals']))
        self.assertIsInstance(parameters[3], tuple)
        self.assertEqual(context.network.to_dict(), source)
        for spec, old_arrays in zip(initial, old):
            for entrance, values in zip(spec.entrances, old_arrays):
                np.testing.assert_array_equal(entrance.source_curve_stations[1], values)
        for spec in prepared:
            self.assertTrue(spec.proof['shared_source_phase']['original_path_identity'])
        self.assertFalse(any(a is b for a, b in zip(initial, prepared)))

    def test_different_original_path_objects_and_excessive_allocations_are_rejected(self):
        context, requests = two_sites()
        request = requests[0]
        spec = OneSidedPairedSpecification.build(context, request.site, request.local_layout,
            request.near_visits, request.opposite_visit, request.near_role)
        with self.assertRaisesRegex(ValueError, 'finite sample budget'):
            shared_source_parameters(spec.entrances, maximum_per_source=2)
        with self.assertRaisesRegex(ValueError, 'repeats an original site'):
            canonical_reprepare(context, (request, request))
        a, b = spec.entrances
        other = copy(b)
        other.paths = (b.paths[0], copy(b.paths[1]))
        with self.assertRaisesRegex(ValueError, 'same original Path object'):
            shared_source_parameters((a, other))

    def test_supplied_station_sampling_is_idempotent_and_retains_real_knots(self):
        context, site, layout, visits, opposite = fixture()
        for sid in (1, 3):
            path = context[model.Key('path', sid)]
            compulsory = (0., .5, np.nextafter(.5, 1.), path.length)
            first = sample_source_parameters(path, 0., path.length, .5, .00025, compulsory)
            second = sample_source_parameters(path, 0., path.length, .5, .00025, first)
            np.testing.assert_array_equal(first, second)
            self.assertTrue(np.any(first == .5))
            for station in path.stations:
                self.assertTrue(np.any(first == station))
            with self.assertRaisesRegex(ValueError, 'finite sample budget'):
                sample_source_parameters(path, 0., path.length, .5, .00025, compulsory, maximum_samples=2)


if __name__ == '__main__':
    unittest.main()
