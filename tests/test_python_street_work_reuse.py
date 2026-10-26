"""Numerical fidelity and bounded work for repeated boundary/dependency queries."""
from dataclasses import replace
import random
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from _python_street_fixtures import fixture_networks
from _python_work_reference import uncached_approach_frame, sorted_list_order as old_order
from tested_python_streets import intersections, model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.conditions import ProfileEdge
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.routes import Route, Piece
from tested_python_streets.updates import DependencyIndex, Task


class LegacyApproach(intersections.Approach):
    """Prior arithmetic, including Exact miter finite differences, without reuse."""
    frame = uncached_approach_frame


class CountedPath:
    def __init__(self, path):
        self.path, self.calls = path, 0

    def __getattr__(self, name):
        return getattr(self.path, name)

    def frame(self, station):
        self.calls += 1
        return self.path.frame(station)


class PythonWorkReuseTests(unittest.TestCase):
    def section(self, mode=3):
        return model.Section(1, 1, ((0., 0., 0.), (20., 2., 0.), (40., -1., 0.)),
                             width=5.5, offset=.4, interpolation=mode)

    def exact(self, first, second):
        first, second = np.asarray(first), np.asarray(second)
        self.assertEqual(first.dtype, second.dtype)
        self.assertEqual(first.shape, second.shape)
        self.assertEqual(first.tobytes(), second.tobytes())

    def test_cached_frames_and_profile_gradients_are_bit_exact_in_both_modes_and_directions(self):
        for mode in (0, 3):
            for end in (0, 1):
                with self.subTest(mode=mode, end=end):
                    path = Path(self.section(mode))
                    actual = intersections.Approach(path, end, Crown(.12, 37., 1))
                    reference = LegacyApproach(path, end, Crown(.12, 37., 1))
                    for station in (-.3, 0., .00005, .1 * path.length,
                                    .55 * path.length, path.length - .00001, path.length):
                        for _ in range(2):
                            for a, b in zip(actual.frame(station), reference.frame(station)):
                                self.exact(a, b)
                            for t in (0., np.array((0., .3, .5, 1.)), np.array(((0., .3), (.7, 1.)))):
                                for a, b in zip(ProfileEdge(actual, station, 0., 1.).sample(t),
                                                ProfileEdge(reference, station, 0., 1.).sample(t)):
                                    self.exact(a, b)
                    # Preserve the separate float32/scalar/vector arithmetic paths.
                    for station in (np.float32(.4), np.float64(.4), 1, np.array((.4,)), np.array(((.4, .8),))):
                        for a, b in zip(actual.frame(station), reference.frame(station)):
                            self.exact(a, b)

    def test_repeated_profile_samples_reuse_three_frames_without_mutable_aliases(self):
        path = CountedPath(Path(self.section()))
        approach = intersections.Approach(path, 0, Crown(.12, 0.))
        rule = ProfileEdge(approach, 12., 0., 1.)
        previous_path = CountedPath(Path(self.section()))
        previous = ProfileEdge(LegacyApproach(previous_path, 0, Crown(.12, 0.)), 12., 0., 1.)
        for count in range(1, 25):
            t = np.linspace(0., 1., count)
            for a, b in zip(rule.sample(t), previous.sample(t)):
                self.exact(a, b)
        self.assertEqual(path.calls, 4)  # Initial origin and three mouth stations.
        self.assertEqual(previous_path.calls, 97)
        expected = tuple(v.copy() for v in approach.frame(12.))
        received = approach.frame(12.)
        for value in received:
            value[:] = 999.
        for a, b in zip(approach.frame(12.), expected):
            self.exact(a, b)
        for station in np.linspace(1., path.length - 1., 80):
            approach.frame(float(station))
        self.assertLessEqual(len(approach._frames), approach.frame_cache_limit)
        # A fresh derived approach cannot retain another source revision's frame.
        moved = replace(self.section(), points=((0., 3., 0.), (20., 5., 0.), (40., 2., 0.)))
        changed = intersections.Approach(Path(moved), 0, Crown(.12, 0.))
        np.testing.assert_array_equal(changed.frame(12.)[0], expected[0] + (0., 3.))

    def test_multisection_outward_route_keeps_exact_frames_across_the_seam(self):
        a = Path(model.Section(1, 7, ((0., 0., 0.), (20., 0., 0.)), width=4.))
        b = Path(model.Section(2, 7, ((20., 0., 0.), (40., 1., 0.)), width=4.))
        path = Route((Piece(a, True), Piece(b, True)))
        actual, reference = (cls(path, 0, Crown(.08, 0.)) for cls in (intersections.Approach, LegacyApproach))
        for station in (-.2, 0., 19.9999, 20., 20.0001, path.length):
            for _ in range(2):
                for a, b in zip(actual.frame(station), reference.frame(station)):
                    self.exact(a, b)

    def test_heap_retains_the_previous_minimum_ready_order_under_input_permutations(self):
        rng = random.Random(173)
        for _ in range(10):
            keys = [model.Key('task', i) for i in range(240)]
            tasks = [Task(key, tuple(rng.sample(keys[:i], min(i, rng.randrange(5)))) +
                          (model.Key('source', i % 7),), None) for i, key in enumerate(keys)]
            reference = old_order(tasks)
            for _ in range(3):
                rng.shuffle(tasks)
                actual = DependencyIndex(tuple(tasks))
                self.assertEqual(actual.order, reference)
                self.assertEqual(actual.affected((model.Key('source', 0),)),
                                 DependencyIndex(tuple(reversed(tasks))).affected((model.Key('source', 0),)))

    def test_six_family_meshes_and_dependencies_equal_the_uncached_reference(self):
        for network in fixture_networks(model):
            with self.subTest(network=network.uuid):
                actual = PythonMeshBackend(network.settings).prepare(network)
                with patch.object(intersections, 'Approach', LegacyApproach):
                    reference = PythonMeshBackend(network.settings).prepare(network)
                self.assertEqual(actual.rebuilt, reference.rebuilt)
                self.assertEqual(actual.dependencies, reference.dependencies)
                self.assertEqual(actual.products.keys(), reference.products.keys())
                for key, product in actual.products.items():
                    if getattr(product, 'mesh', None) is not None:
                        self.exact(product.mesh.xyz, reference.products[key].mesh.xyz)
                        self.exact(product.mesh.quads, reference.products[key].mesh.quads)


if __name__ == '__main__':
    unittest.main()
