"""Source coverage cancels whole interfaces and preserves actual holes."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.gapped_source_validation import mesh_union_loops


def square(first, last):
    a, b = first, last
    return Mesh(np.array(((a, 0., 0.), (b, 0., 0.), (b, 1., 0.), (a, 1., 0.))), np.array(((0, 1, 2, 3),)))


class GappedSourceValidationTests(unittest.TestCase):
    def test_exact_whole_interfaces_cancel_and_real_hole_survives(self):
        loops = mesh_union_loops((square(0., 1.), square(1., 2.)))
        self.assertEqual(len(loops), 1); self.assertFalse(loops[0].hole)
        xyz = np.array(((0., 0., 0.), (3., 0., 0.), (3., 3., 0.), (0., 3., 0.),
                        (1., 1., 0.), (2., 1., 0.), (2., 2., 0.), (1., 2., 0.)))
        ring = Mesh(xyz, np.array(((0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))))
        ring.validate()
        loops = mesh_union_loops((ring,))
        self.assertEqual(len(loops), 2)
        self.assertEqual(sum(loop.hole for loop in loops), 1)
        np.testing.assert_array_equal(next(loop.xy for loop in loops if loop.hole), xyz[[4, 7, 6, 5], :2])

    def test_touching_or_overlapping_edges_cannot_be_repaired_into_coverage(self):
        first = square(0., 1.)
        with self.assertRaisesRegex(ValueError, 'equally directed'):
            mesh_union_loops((first, first))
        second = square(1., 2.)
        second.xyz[:, 1] += 1.
        with self.assertRaisesRegex(ValueError, 'boundary branch'):
            mesh_union_loops((first, second))
        # Splitting only one side of an otherwise identical interface creates
        # a T junction; source coverage must not infer a numerical weld.
        partial = Mesh(np.array(((1., 0., 0.), (2., 0., 0.), (2., .5, 0.), (1., .5, 0.))),
                       np.array(((0, 1, 2, 3),)))
        with self.assertRaisesRegex(ValueError, 'boundary branch'):
            mesh_union_loops((first, partial))


if __name__ == '__main__': unittest.main()
