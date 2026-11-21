"""Thin separator disks retain the existing grass quality/publication contract."""
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.four_side_charts import FourSidedQuadChart
from tested_python_streets.hole_disk_chart import build_hole_disk
from tested_python_streets.mesh import publication_coordinates
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.regions import signed_area


def lens(angle=0.):
    t = np.linspace(0., 1., 61)
    x, y = 10.*t, .075*4.*t*(1.-t)
    xy = np.vstack((np.column_stack((x, -y))[:-1], np.column_stack((x, y))[::-1][:-1]))
    c, s = math.cos(angle), math.sin(angle)
    return xy @ np.array(((c, s), (-s, c))) + (1000., -1100.)


class HoleDiskTests(unittest.TestCase):
    def test_sharp_fixed_lens_is_positive_conforming_and_valid_after_float32_storage(self):
        for angle in (0., .713):
            with self.subTest(angle=angle):
                boundary = lens(angle)
                result = build_hole_disk(PatternLibrary(), boundary, .7, 3000, 100)
                mesh = result.mesh
                np.testing.assert_array_equal(mesh.xyz[mesh.seams[0], :2], boundary)
                report = result.metadata['quality']
                self.assertLessEqual(report['maximum_edge_ratio'], 8.)
                self.assertLess(report['minimum_corner_sine'], .03)
                # Pavement atlases retain their separate sine requirement.
                self.assertFalse(FourSidedQuadChart(PatternLibrary()).accepted(report))
                uv = mesh.xyz[mesh.quads, :2].copy()
                mesh = mesh.with_uv_layer('physical', uv)
                mesh.xyz[:, 2] = .08
                stored = mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True))
                self.assertLessEqual(stored.validate()['maximum_edge_ratio'], 8.)
                np.testing.assert_array_equal(stored.uv_layers['physical'], uv)
                np.testing.assert_array_equal(stored.quad_slopes(), np.zeros(len(stored.quads)))
                edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
                unique, counts = np.unique(edges, axis=0, return_counts=True)
                self.assertEqual(len(mesh.xyz)-len(unique)+len(mesh.quads), 1)
                self.assertEqual(int(np.sum(counts == 1)), len(boundary))
                self.assertTrue(set(counts).issubset({1, 2}))
                area = sum(signed_area(quad) for quad in mesh.xyz[mesh.quads, :2])
                self.assertAlmostEqual(area, signed_area(boundary), places=8)


if __name__ == '__main__':
    unittest.main()
