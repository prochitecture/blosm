"""Unequal separator-bank samples remain exact in conforming capped disks."""
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.hole_fibre_chart import build_fibre_hole, _pair_rows
from tested_python_streets.mesh import publication_coordinates
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.regions import signed_area


def unequal_lens(second_count, angle):
    x = np.linspace(0., 12., 43)
    first = np.column_stack((x, -.075*4.*x/12.*(1.-x/12.)))
    x = 12.*np.linspace(0., 1., second_count)**1.01
    second = np.column_stack((x, .075*4.*x/12.*(1.-x/12.)))
    points = np.vstack((first[:-1], second[::-1][:-1]))
    c, s = math.cos(angle), math.sin(angle)
    return points @ np.array(((c, s), (-s, c))) + (1000., -1100.)


class HoleFibreTests(unittest.TestCase):
    def test_unequal_banks_retain_every_point_and_only_add_chord_samples(self):
        for count, angle in ((51, 0.), (52, .713)):
            with self.subTest(count=count, angle=angle):
                boundary = unequal_lens(count, angle)
                result = build_fibre_hole(PatternLibrary(), boundary, .7, 3000)
                mesh = result.mesh
                rim = mesh.xyz[mesh.seams[0], :2]
                positions = []
                for point in boundary:
                    at = np.flatnonzero(np.all(rim == point, axis=1))
                    self.assertEqual(len(at), 1)
                    positions.append(int(at[0]))
                self.assertEqual(positions, sorted(positions))
                self.assertGreater(len(rim), len(boundary))
                self.assertEqual(len(rim) % 2, 0)
                for i, start in enumerate(positions):
                    end = positions[(i+1) % len(positions)]
                    if i+1 == len(positions):
                        end += len(rim)
                    p, q = boundary[i], boundary[(i+1) % len(boundary)]
                    samples = rim[np.arange(start, end+1) % len(rim)]
                    d = q-p
                    fraction = (samples-p) @ d/(d @ d)
                    np.testing.assert_allclose(samples, p+fraction[:, None]*d, atol=1e-8, rtol=0.)
                    self.assertTrue(np.all(np.diff(fraction) > 0.))
                self.assertAlmostEqual(signed_area(rim), signed_area(boundary), places=8)
                self.assertLessEqual(mesh.validate()['maximum_edge_ratio'], 8.)
                stored = mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True))
                self.assertLessEqual(stored.validate()['maximum_edge_ratio'], 8.)
                edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
                unique, counts = np.unique(edges, axis=0, return_counts=True)
                self.assertEqual(len(mesh.xyz)-len(unique)+len(mesh.quads), 1)
                self.assertEqual(int(np.sum(counts == 1)), len(rim))
                self.assertTrue(set(counts).issubset({1, 2}))

    def test_nonmonotone_and_excessive_pairing_and_cell_allocation_reject(self):
        points = np.column_stack((np.linspace(0., 1., 800), np.ones(800)))
        with self.assertRaisesRegex(ValueError, 'pair budget'):
            _pair_rows(points, points-(0., 1.), np.zeros(2), np.array((1., 0.)), np.array((0., 1.)))
        with self.assertRaisesRegex(ValueError, 'monotone'):
            _pair_rows(points[[0, 2, 1, 3]], points[:4], np.zeros(2), np.array((1., 0.)), np.array((0., 1.)))
        with self.assertRaises(ValueError):
            build_fibre_hole(PatternLibrary(), unequal_lens(51, 0.), .7, 1)
        with self.assertRaisesRegex(ValueError, 'spacing or cell budget'):
            build_fibre_hole(PatternLibrary(), unequal_lens(51, 0.), 0., 3000)


if __name__ == '__main__':
    unittest.main()
