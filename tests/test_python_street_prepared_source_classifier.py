"""Exact source-proof classification, bounded cache and explicit lifetimes."""
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets import finite_containment, source_strip_containment
from tested_python_streets.regions import _classify
from tested_python_streets.prepared_source_classifier import SourceProofClassifier
from test_python_street_source_strip_containment import context, group


def parity(classifier, points, polygon, epsilon):
    for expected, actual in zip(_classify(points, polygon, epsilon), classifier(points, polygon, epsilon)):
        np.testing.assert_array_equal(expected, actual)


class PreparedSourceClassifierTests(unittest.TestCase):
    def test_original_edge_ties_and_buffer_views_share_one_preparation(self):
        polygon = np.array(((0., 0.), (1., 0.), (1., 1.), (0., 1.)))
        points = np.array(((0., 0.), (.5, .5), (1., 1.), (1.+1e-8, 1.+1e-8)))
        classifier = SourceProofClassifier()
        for _ in range(10):
            for epsilon in (1e-7, 1.):
                parity(classifier, points, polygon.view(), epsilon)
        self.assertEqual(classifier.preparations, 1)
        self.assertEqual(len(classifier._cache), 1)
        self.assertEqual(classifier(points, polygon, 1.)[1].tolist(), [0, 0, 1, 1])

    def test_rotations_translation_vertices_and_numerical_boundary_bands_are_exact(self):
        random = np.random.default_rng(483)
        angle = np.linspace(0., 2.*np.pi, 128, endpoint=False)
        radius = 5.+.6*np.cos(3.*angle)
        shape = np.stack((radius*np.cos(angle), radius*np.sin(angle)), axis=1)
        sampled = np.arange(0, len(shape), 8)
        for rotation_angle in (0., .713):
            rotation = np.array(((np.cos(rotation_angle), -np.sin(rotation_angle)),
                                 (np.sin(rotation_angle), np.cos(rotation_angle))))
            for translation in ((0., 0.), (1200., -1800.), (1e12, -1e12)):
                polygon = shape@rotation.T+translation
                edge = (np.roll(polygon, -1, axis=0)-polygon)[sampled]
                normal = np.stack((-edge[:, 1], edge[:, 0]), axis=1)/np.linalg.norm(edge, axis=1)[:, None]
                center = .5*(polygon+np.roll(polygon, -1, axis=0))[sampled]
                classifier = SourceProofClassifier()
                for epsilon in (1e-10, 1e-7, .002):
                    points = np.vstack((polygon[sampled], np.nextafter(polygon[sampled], np.inf),
                        np.nextafter(polygon[sampled], -np.inf),
                        *(center+amount*epsilon*normal for amount in (-1.001, -.999, 0., .999, 1.001)),
                        random.uniform(-7., 7., (32, 2))+translation))
                    parity(classifier, points, polygon, epsilon)
        long = np.array(((-1e16, -1e16), (1e16, 1e16), (1e16, 1e16+100.), (-1e16, -1e16+100.)))
        parity(SourceProofClassifier(), np.array(((0., 0.), (0., 50.), (1e16, 1e16))), long, 1e-7)

    def test_pathological_indices_and_full_cache_decline_without_repreparation(self):
        classifier = SourceProofClassifier()
        points = np.array(((0., 0.), (.5, .5)))
        # A source exterior and eight distinct holes fit without eviction.
        polygons = [np.array(((0., 0.), (1., 0.), (1., 1.), (0., 1.)))+i*2. for i in range(20)]
        for polygon in polygons[:9]:
            parity(classifier, points, polygon, 1e-7)
        self.assertEqual((len(classifier._cache), classifier.preparations), (9, 9))
        for polygon in polygons:
            parity(classifier, points, polygon, 1e-7)
        self.assertEqual((len(classifier._cache), classifier.preparations), (16, 16))
        for _ in range(10):
            parity(classifier, points, polygons[-1], 1e-7)
        self.assertEqual(classifier.preparations, 16)
        self.assertLessEqual(classifier.allocated_bytes, 8*1024*1024)
        pathological = np.stack((np.arange(2000, dtype=float), np.arange(2000) % 2), axis=1)
        declined = SourceProofClassifier()
        for _ in range(4):
            parity(declined, points, pathological, 1e-7)
        self.assertEqual(declined.preparations, 1)
        self.assertEqual(declined.allocated_bytes, 0)
        class TinyCache(SourceProofClassifier):
            maximum_cache_bytes = 1
        tiny = TinyCache()
        for _ in range(4):
            parity(tiny, points, polygons[0], 1e-7)
        self.assertEqual((tiny.preparations, tiny.allocated_bytes), (1, 0))

    def test_unsupported_arithmetic_domains_keep_baseline_results(self):
        polygon = np.array(((0., 0.), (1., 0.), (1., 2e-12), (0., 2e-12)))
        classifier = SourceProofClassifier()
        with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
            for points in (np.array(((1e308, 0.),)), np.array(((0, 0), (1, 0)), dtype=np.int64),
                           np.array(((0., 0.),), dtype=np.float32)):
                parity(classifier, points, polygon, 1e-10)
            degenerate = np.array(((0., 0.), (1., 0.), (1., 0.), (0., 1.)))
            parity(classifier, np.array(((.2, .2),)), degenerate, 1e-7)

    def test_optional_classifier_keeps_continuous_intersections_and_capsules(self):
        polygon = np.array(((-1., -2.), (11., -2.), (11., 2.), (8., 2.), (8., -1.),
            (5.0001, -1.), (5.0001, -.0001), (4.9999, -.0001), (4.9999, -1.), (2., -1.), (2., 2.), (-1., 2.)))
        rectangle = np.array(((-1., .001), (11., .001), (11., 2.), (-1., 2.)))
        for boundary in (polygon, rectangle):
            with patch.object(finite_containment, '_capsule_intervals',
                              wraps=finite_containment._capsule_intervals) as capsules:
                original = finite_containment.contains_segment((0., 0.), (10., 0.), boundary, tolerance=.002)
                original_calls = capsules.call_args_list[:]
                capsules.reset_mock()
                prepared = finite_containment.contains_segment((0., 0.), (10., 0.), boundary,
                    tolerance=.002, classifier=SourceProofClassifier())
                self.assertEqual(original, prepared)
                self.assertEqual(len(original_calls), len(capsules.call_args_list))
                for first, second in zip(original_calls, capsules.call_args_list):
                    for a, b in zip(first.args, second.args):
                        np.testing.assert_array_equal(a, b)

    def test_batch_matches_scalar_proof_across_an_inlet_and_tolerance_band(self):
        polygon = np.array(((0., 0.), (8., 0.), (8., 4.), (4.001, 4.),
                            (4.001, 1.), (3.999, 1.), (3.999, 4.), (0., 4.)))
        rng = np.random.default_rng(131)
        samples = rng.uniform((-1., -1.), (9., 5.), size=(2, 90, 2))
        first = np.vstack((samples[0], polygon, (0., 1.5), (0., -.001), (1., 1.)))
        last = np.vstack((samples[1], np.roll(polygon, 1, axis=0), (8., 1.5), (8., -.001), (1., 1.)))
        for offset in (0., 2351.):
            for inside in (True, False):
                for tolerance in (.00001, .002):
                    args = dict(inside=inside, tolerance=tolerance)
                    expected = [finite_containment.segment_on_polygon_side(a+offset, b+offset,
                        polygon+offset, **args) for a, b in zip(first, last)]
                    actual = finite_containment.segments_on_polygon_side(first+offset, last+offset,
                        polygon+offset, classifier=SourceProofClassifier(), **args)
                    np.testing.assert_array_equal(actual, expected)

    def test_new_source_proof_after_boundary_edit_builds_a_fresh_classifier(self):
        boundary = group(((-1., -2.), (11., -2.), (11., 2.), (-1., 2.)))
        instances = []
        def fresh():
            classifier = SourceProofClassifier(); instances.append(classifier)
            return classifier
        with patch.object(source_strip_containment, 'SourceProofClassifier', side_effect=fresh):
            source_strip_containment.certify(context(), boundary, {1})
            boundary.loops[0].xy[2:, 1] = .5
            with self.assertRaises(source_strip_containment.SourceStripContainmentError):
                source_strip_containment.certify(context(), boundary, {1})
        self.assertEqual(len(instances), 2)
        self.assertIsNot(instances[0], instances[1])


if __name__ == '__main__':
    unittest.main()
