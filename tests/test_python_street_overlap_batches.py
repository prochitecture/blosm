"""Batched overlap decisions against an independent all-pairs SAT oracle."""
import itertools
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.model import Key
from tested_python_streets.validation import validate_overlaps, _CellGrid, _PAIR_BATCH


def rectangle(x, y, w=1., h=1., angle=0.):
    points = np.array(((-w/2, -h/2), (w/2, -h/2), (w/2, h/2), (-w/2, h/2)))
    c, s = np.cos(angle), np.sin(angle)
    return points @ np.array(((c, s), (-s, c))) + (x, y)


def mesh(quads):
    xy = np.asarray(quads).reshape(-1, 2)
    return Mesh(np.column_stack((xy, np.zeros(len(xy)))), np.arange(len(xy)).reshape(-1, 4))


def oracle(meshes, changed, tolerance):
    faces = [(key, quad) for key, value in meshes.items() for quad in value.xyz[value.quads, :2]]
    hits, candidates = [], 0
    for j in range(len(faces)):
        key_b, b = faces[j]
        for i in range(j):
            key_a, a = faces[i]
            if key_a not in changed and key_b not in changed:
                continue
            if np.any(a.max(axis=0) - b.min(axis=0) <= tolerance) or np.any(b.max(axis=0) - a.min(axis=0) <= tolerance):
                continue
            candidates += 1
            # Independent scalar axis loop, with the same original face order
            # and local origin. No grid, batch kernel or triangulation is used.
            first, second = b - b[0], a - b[0]
            depths = []
            for polygon in (first, second):
                for p, q in zip(polygon, np.roll(polygon, -1, axis=0)):
                    delta = q - p
                    axis = np.array((delta[1], -delta[0])) / np.linalg.norm(delta)
                    pa, pb = first @ axis, second @ axis
                    depths.append(min(pa.max(), pb.max()) - max(pa.min(), pb.min()))
            if min(depths) > tolerance:
                hits.append((key_b, key_a))
    return hits, candidates


class PythonOverlapBatches(unittest.TestCase):
    def check_oracle(self, meshes, changed, cell_size=1., tolerance=1e-6):
        hits, pairs = oracle(meshes, changed, tolerance)
        if hits:
            with self.assertRaisesRegex(ValueError, 'interiors overlap'):
                validate_overlaps(meshes, changed=changed, cell_size=cell_size, tolerance=tolerance)
        else:
            result = validate_overlaps(meshes, changed=changed, cell_size=cell_size, tolerance=tolerance)
            self.assertEqual(result['candidate_pairs'], pairs)
            self.assertEqual(result['quads'], sum(len(m.quads) for m in meshes.values()))

    def test_random_convex_quads_and_changed_subsets_match_all_pairs(self):
        rng = np.random.default_rng(98173)
        for trial in range(12):
            meshes = {Key('road', k): mesh([rectangle(*rng.uniform(-8., 8., 2),
                *rng.uniform(.25, 1.7, 2), rng.uniform(-np.pi, np.pi)) for _ in range(4)])
                for k in range(6)}
            keys = list(meshes)
            for changed in (set(keys), set(keys[::2]), {keys[-1]}, set()):
                with self.subTest(trial=trial, changed=changed):
                    self.check_oracle(meshes, changed, cell_size=.7 if trial % 2 else 2.)

    def test_tolerance_edge_contacts_and_translations(self):
        for shift, penetration in itertools.product((0., -10000., 1000000.), (0., .5e-6, 2e-6)):
            meshes = {Key('road', 1): mesh([rectangle(shift, shift)]),
                      Key('road', 2): mesh([rectangle(shift + 1. - penetration, shift)])}
            for changed in ({Key('road', 1)}, {Key('road', 2)}):
                with self.subTest(shift=shift, penetration=penetration, changed=changed):
                    self.check_oracle(meshes, changed, cell_size=.5)

    def test_a_pair_sharing_many_cells_is_tested_only_once(self):
        diamond = np.array(((0., -2.), (2., 0.), (0., 2.), (-2., 0.)))
        meshes = {Key('road', 1): mesh([diamond]), Key('road', 2): mesh([diamond + (3., 3.)])}
        for changed in (set(meshes), {Key('road', 1)}, {Key('road', 2)}):
            result = validate_overlaps(meshes, changed=changed, cell_size=.25)
            self.assertEqual(result['candidate_pairs'], 1)

    def test_self_overlap_is_checked_but_unchanged_pairs_are_ignored(self):
        key = Key('road', 1)
        with self.assertRaisesRegex(ValueError, 'road:1 and road:1'):
            validate_overlaps({key: mesh([rectangle(0., 0.), rectangle(.1, .1)])})
        meshes = {key: mesh([rectangle(0., 0.)]), Key('road', 2): mesh([rectangle(.1, .1)]),
                  Key('road', 3): mesh([rectangle(10., 10.)])}
        self.check_oracle(meshes, {Key('road', 3)})

    def test_multiple_expansion_and_pair_batches_preserve_a_large_disjoint_grid(self):
        quads = [rectangle(x, y, .8, .8) for x in range(100) for y in range(100)]
        meshes = {Key('road', 1): mesh(quads)}
        result = validate_overlaps(meshes, cell_size=1.)
        self.assertEqual(result['quads'], 10000)
        self.assertEqual(result['candidate_pairs'], 0)

    def test_dense_cell_joins_stream_in_bounded_batches(self):
        quads = [rectangle(i % 40, i // 40, .8, .8) for i in range(400)]
        result = validate_overlaps({Key('road', 1): mesh(quads)}, cell_size=100.)
        self.assertEqual(result['candidate_pairs'], 0)

    def test_one_dense_query_row_cannot_exceed_the_pair_workspace(self):
        count = _PAIR_BATCH + 123
        first = np.zeros((count + 1, 2), dtype=np.int64)
        spans = np.ones_like(first)
        grid = _CellGrid(first, spans, np.arange(count), 2 * count)
        batches = list(grid.pairs(first, spans, np.array((count,))))
        self.assertEqual(len(batches), 2)
        self.assertLessEqual(max(len(target) for target, partner, slots in batches), _PAIR_BATCH)
        np.testing.assert_array_equal(np.concatenate([partner for target, partner, slots in batches]),
                                      np.arange(count))

    def test_explicit_entry_spatial_and_parameter_budgets(self):
        meshes = {Key('road', 1): mesh([rectangle(.5, .5)])}
        with self.assertRaisesRegex(ValueError, 'entry budget'):
            validate_overlaps(meshes, cell_size=.5, max_index_entries=8)
        with self.assertRaisesRegex(ValueError, 'spatial budget'):
            validate_overlaps({Key('road', 1): mesh([rectangle(0., 0., 65., 65.)])}, cell_size=1.)
        for options in ({'cell_size': 0.}, {'tolerance': -1.}, {'max_index_entries': 0}):
            with self.subTest(options=options), self.assertRaisesRegex(ValueError, 'limits'):
                validate_overlaps(meshes, **options)


if __name__ == '__main__':
    unittest.main()
