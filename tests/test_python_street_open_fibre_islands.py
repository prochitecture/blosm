"""A real narrow source-tip sector keeps its banks through direct quad grading."""
import unittest
from types import SimpleNamespace
import numpy as np
import _python_street_package
from tested_python_streets.bundle_angular_curves import AngularCoordinates, RayCurb
from tested_python_streets.bundle_open_fibre_islands import OpenFibreIsland
from tested_python_streets.four_side_charts import quality
from tested_python_streets.mesh import Mesh
from tested_python_streets.model import Key
from tested_python_streets.patterns import PatternLibrary


class Frame(AngularCoordinates):
    def __init__(self, angle, origin):
        self.scale, self.direction, self.angle = 10000., 1., angle-.5*np.pi
        self.center = np.array(origin)-self.scale*np.array((np.cos(self.angle), np.sin(self.angle)))


def fixture(*, end=0, angle=0., origin=(-170., -544.)):
    frame = Frame(angle, origin)
    first = np.unique(np.r_[np.linspace(0., 20., 41), .21])
    second = np.unique(np.r_[0., .015 if angle == 0. else .05,
                             np.arange(.4, 20., .5), 20.])
    if end:
        first, second = 20.-first[::-1], 20.-second[::-1]
    distance = first if end == 0 else 20.-first
    upper = frame.world(np.c_[first, .0128*distance+.0008*distance**2])
    lower = frame.world(np.c_[second, np.zeros(len(second))])
    rays = (None, RayCurb(frame, upper), RayCurb(frame, lower), None)
    context = SimpleNamespace(patterns=PatternLibrary(), network=SimpleNamespace(settings=SimpleNamespace(spacing=1.)))
    domain = SimpleNamespace(context=context, frame=frame, rays=rays,
                             field=SimpleNamespace(ceiling=.08), bundle=SimpleNamespace(key=Key('bundle', 1)))
    root = 20. if end else 0.
    return domain, {root: frame.world(np.array((root, 0.)))}, (upper, lower)


class OpenFibreIslandTests(unittest.TestCase):
    def test_long_narrow_banks_use_bounded_columns_after_the_original_family_fails(self):
        for end, angle, origin in ((0, 0., (0., 0.)), (1, .713, (1500., -1700.))):
            frame = Frame(angle, origin)
            length = 75.
            first = np.unique(np.r_[0., np.arange(.3, length, .7), length])
            second = np.unique(np.r_[0., np.arange(.21, length, .592), length])
            distance = first if end == 0 else length-first
            width = .0168*distance-.0003*distance**2+.000003*distance**3
            upper = frame.world(np.c_[first, width])
            lower = frame.world(np.c_[second, np.zeros(len(second))])
            rays = (None, RayCurb(frame, upper), RayCurb(frame, lower), None)
            context = SimpleNamespace(patterns=PatternLibrary(),
                network=SimpleNamespace(settings=SimpleNamespace(spacing=.7)))
            domain = SimpleNamespace(context=context, frame=frame, rays=rays,
                field=SimpleNamespace(ceiling=.08), bundle=SimpleNamespace(key=Key('bundle', 1)))
            root = length if end else 0.
            builder = OpenFibreIsland(domain)
            mesh, banks, report = builder.build(0., length,
                {root: frame.world(np.array((root, 0.)))}, 1)
            self.assertIn('width_adaptive_columns', report)
            self.assertEqual(len(report['fibre_candidate_failures']), 3)
            self.assertLessEqual(len(mesh.quads), 5000)
            self.assertLessEqual(quality(Mesh(mesh.xyz.astype(np.float32), mesh.quads,
                'grass'))['maximum_edge_ratio'], 8.)
            self.assertTrue(report['exact_physical_area'] and report['conforming_disk'])
            for source, (actual, keys), proof in zip((upper, lower), banks, report['bank_sampling']):
                np.testing.assert_array_equal(actual[proof['original_vertex_positions']], source)
                self.assertLessEqual(proof['maximum_chord_error'], 1e-8)

    def test_protected_unequal_banks_keep_sharp_tip_and_exact_physical_chords(self):
        for end, angle, origin in ((0, 0., (-170., -544.)), (1, .713, (1200., -1800.))):
            with self.subTest(end=end):
                domain, roots, originals = fixture(end=end, angle=angle, origin=origin)
                before = tuple(points.copy() for points in originals)
                builder = OpenFibreIsland(domain)
                mesh, banks, report = builder.build(0., 20., roots, 4)
                stored = Mesh(mesh.xyz.astype(np.float32), mesh.quads, 'grass')
                self.assertLessEqual(max(quality(mesh)['maximum_edge_ratio'], quality(stored)['maximum_edge_ratio']), 8.)
                self.assertLess(report['original_tip_sine'], .03)
                self.assertTrue(report['exact_physical_area'] and report['conforming_disk'])
                self.assertEqual(len(builder.keys), len(mesh.xyz))
                self.assertEqual(len(set(builder.keys)), len(mesh.xyz))
                self.assertLessEqual(len(mesh.quads), 5000)
                lookup = {key: i for i, key in enumerate(builder.keys)}
                for source, saved, (points, keys), proof in zip(originals, before, banks, report['bank_sampling']):
                    np.testing.assert_array_equal(source, saved)
                    np.testing.assert_array_equal(points[proof['original_vertex_positions']], source)
                    np.testing.assert_array_equal(mesh.xyz[[lookup[key] for key in keys], :2], points)
                    self.assertLessEqual(proof['maximum_chord_error'], 1e-8)
                    self.assertGreater(len(proof['added_chord_refs']), 0)
                points, keys = builder.cuts[1-end]
                np.testing.assert_array_equal(mesh.xyz[[lookup[key] for key in keys], :2], points)
                local = domain.frame.local(points)
                self.assertLess(np.ptp(local[:, 0]), 1e-8)
                self.assertEqual(report['flat_cut_intervals'][end], None)

    def test_finite_budget_missing_tip_and_reversed_lateral_roles_reject(self):
        domain, roots, _ = fixture()
        with self.assertRaisesRegex(ValueError, 'finite.*budget'):
            OpenFibreIsland(domain, max_cells=5001)
        with self.assertRaisesRegex(ValueError, 'quad budget'):
            OpenFibreIsland(domain, max_cells=10).build(0., 20., roots, 0)
        with self.assertRaisesRegex(ValueError, 'one finite physical tip'):
            OpenFibreIsland(domain).build(0., 20., {}, 0)
        with self.assertRaisesRegex(ValueError, 'one finite physical tip'):
            OpenFibreIsland(domain).build(0., 20., {**roots, 20.: np.zeros(2)}, 0)
        domain.rays = (None, domain.rays[2], domain.rays[1], None)
        with self.assertRaisesRegex(ValueError, 'lateral order'):
            OpenFibreIsland(domain).build(0., 20., roots, 0)


if __name__ == '__main__': unittest.main()
