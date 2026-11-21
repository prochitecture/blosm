"""A finite Bundle piece can join merged or separated actual physical mouths."""
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.bundle_angular_curves import AngularCoordinates, AngularCrown
from tested_python_streets.bundle_common_fibres import CommonFibreCuts
from tested_python_streets.bundle_merged_body import MergedRunBody
from tested_python_streets.mesh import Mesh
from tested_python_streets.model import Key, Settings
from tested_python_streets.patterns import PatternLibrary


def fixture(separated=False, *, angle=0., translation=(0., 0.)):
    frame = AngularCoordinates.__new__(AngularCoordinates)
    frame.center = np.asarray(translation); frame.scale = 50.; frame.angle = angle; frame.direction = 1.
    stations = np.linspace(-20., 20., 161)
    inner = .15 * (stations + .123) if separated == 'mixed' else np.full(len(stations), 1. if separated else -1.)
    curves = [frame.world(np.c_[stations, np.broadcast_to(y, stations.shape)])
              for y in (4., inner, -inner, -4.)]
    field = AngularCrown(frame, curves, datum=0., rise=.08)
    paths = (SimpleNamespace(exact=False), SimpleNamespace(exact=False))
    approach = SimpleNamespace(frame_chart=frame, field=field, paths=paths, sign=1.)
    cuts = {}
    for end, at in enumerate((-10., 10.)):
        edges, peaks = field.edges(at), field.peaks(at)
        for role in (0, 1):
            y = np.linspace(edges[0 if role == 0 else 3], peaks[role], 7)
            xy = frame.world(np.c_[np.full(len(y), at), y])
            cuts[role, end] = np.c_[xy, field.values(xy)]
    settings = Settings(resolution=1, exaggeration=1.)
    context = SimpleNamespace(network=SimpleNamespace(settings=settings), patterns=PatternLibrary())
    bundle = SimpleNamespace(key=Key('bundle', 1))
    return frame, field, approach, cuts, context, bundle


class CommonFibreBodyTests(unittest.TestCase):
    def test_merged_and_separated_both_end_profiles_are_actual_direct_quad_interfaces(self):
        for separated in (False, True, 'mixed'):
            with self.subTest(separated=separated):
                frame, field, approach, cuts, context, bundle = fixture(separated)
                before = {key: points.copy() for key, points in cuts.items()}
                body = MergedRunBody()
                pavement, grass, report = body.build_between(context, bundle, approach, cuts, namespace=('piece', 7))
                for key in cuts:
                    np.testing.assert_array_equal(cuts[key], before[key])
                    np.testing.assert_array_equal(body.cuts[key], cuts[key])
                    actual = {name: point for name, point in zip(body.pavement_keys, pavement.xyz)}
                    for name, point in zip(body.cut_keys[key], cuts[key]):
                        np.testing.assert_array_equal(actual[name], point)
                self.assertEqual(grass is not None, bool(separated))
                self.assertLess(len(pavement.quads) + (len(grass.quads) if grass is not None else 0), 2000)
                self.assertLessEqual(Mesh(pavement.xyz.astype(np.float32), pavement.quads).validate()['maximum_edge_ratio'], 8.)
                self.assertLessEqual(frame.longitudinal_bound(pavement), .03 + 1e-8)
                self.assertLessEqual(pavement.slope_statistics()['maximum'], .1 + 1e-8)
                self.assertTrue(report['all_direct_quads'])
                self.assertEqual(report['prepared_common_fibres']['merged_ends'],
                                 [0] if separated == 'mixed' else [] if separated else [0, 1])
                self.assertEqual(body.cut_keys[0, 0][-1] == body.cut_keys[1, 0][-1], separated is not True)

    def test_translated_rotated_piece_retains_the_same_physical_profiles_and_quality_limits(self):
        frame, field, approach, cuts, context, bundle = fixture('mixed', angle=.713, translation=(1800., -1200.))
        body = MergedRunBody()
        pavement, grass, report = body.build_between(context, bundle, approach, cuts, namespace=('piece', 9))
        self.assertIsNotNone(grass)
        self.assertEqual(report['prepared_common_fibres']['merged_ends'], [0])
        for key in cuts:np.testing.assert_array_equal(body.cuts[key], cuts[key])
        for details in report['prepared_common_fibres']['stored_quality']:
            self.assertLessEqual(details['maximum_edge_ratio'], 8.)
            self.assertGreaterEqual(details['minimum_corner_sine'], .03)

    def test_wrong_fibres_incomplete_profiles_and_extrapolation_are_rejected(self):
        frame, field, approach, cuts, context, bundle = fixture()
        mutations = []
        wrong = {key: points.copy() for key, points in cuts.items()}; wrong[0, 0] = wrong[0, 0][1:]
        mutations.append((wrong, 'curb or crown'))
        wrong = {key: points.copy() for key, points in cuts.items()}; wrong[1, 1][:, :2] = frame.world(
            frame.local(wrong[1, 1][:, :2]) + (1., 0.))
        mutations.append((wrong, 'same actual fibre'))
        wrong = {key: points.copy() for key, points in cuts.items()}; wrong[0, 1][3, 2] += .01
        mutations.append((wrong, 'physical profile'))
        wrong = {key: points.copy() for key, points in cuts.items()}
        for role in (0, 1):
            uv = frame.local(wrong[role, 1][:, :2]); uv[:, 0] = 21.
            wrong[role, 1][:, :2] = frame.world(uv)
        mutations.append((wrong, 'extrapolate'))
        for changed, message in mutations:
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                CommonFibreCuts(frame, field, changed, 'bad')


if __name__ == '__main__':unittest.main()
