"""Stored precision must not manufacture overlaps between unchanged and new parts."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh, publication_coordinates
from tested_python_streets.model import Key
from tested_python_streets.validation import validate_overlaps


class PublicationCoordinatesTests(unittest.TestCase):
    def test_rebuilt_and_retained_seams_are_compared_at_the_same_precision(self):
        x = 1000.000035
        left = Mesh(np.array(((999., 0., 0.), (x, 0., 0.), (x, 1., 0.), (999., 1., 0.))), np.array(((0, 1, 2, 3),)))
        right = Mesh(np.array(((x, 0., 0.), (1001., 0., 0.), (1001., 1., 0.), (x, 1., 0.))), left.quads)
        stored = Mesh(publication_coordinates(left.xyz, 1., physical=True), left.quads)
        with self.assertRaisesRegex(ValueError, 'interiors overlap'):
            validate_overlaps({Key('left', 0): stored, Key('right', 0): right}, tolerance=2e-5)
        rebuilt = Mesh(publication_coordinates(right.xyz, 1., physical=True), right.quads)
        validate_overlaps({Key('left', 0): stored, Key('right', 0): rebuilt}, tolerance=1e-6)
        np.testing.assert_array_equal(stored.xyz[[1, 2]], rebuilt.xyz[[0, 3]])

    def test_rounded_bounds_match_stored_vertices_at_different_scene_scales(self):
        points = np.array(((1000.000035, -500.000017, 1.7), (1001.000011, -499.000009, 1.78)))
        for scale in (1., .01, 100.):
            with self.subTest(scale=scale):
                actual = publication_coordinates(points, scale, physical=True)
                bounds = np.array((points[:, :2].min(axis=0), points[:, :2].max(axis=0)))
                expected = np.array((actual[:, :2].min(axis=0), actual[:, :2].max(axis=0)))
                np.testing.assert_array_equal(publication_coordinates(bounds, scale, physical=True), expected)

    def test_unrepresentable_coordinates_and_invalid_scales_are_rejected(self):
        for scale in (0., -1., float('nan')):
            with self.assertRaisesRegex(ValueError, 'positive scene scale'):
                publication_coordinates(((0., 0., 0.),), scale)
        with self.assertRaisesRegex(ValueError, 'float32 publication range'):
            publication_coordinates(((1e100, 0., 0.),), 1.)
