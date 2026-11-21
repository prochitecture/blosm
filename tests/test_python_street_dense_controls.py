"""Certified smoothing retries for dense OSM controls, without source edits."""
from dataclasses import replace
import unittest
import numpy as np

import _python_street_package
from tested_python_streets.model import Section, Network, Key, Intersection, End, Connection
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.mesh import publication_coordinates
from tested_python_streets.paths import Path
from tested_python_streets.local_smoothing import coalesced_source_indices, within_polyline_envelope


def bergamo_dense_bend():
    # Original Bergamo Section997: tiny rounded controls before a long source
    # run previously demanded overlapping independent local corner setbacks.
    return Section(997, 997, (
        (1426.601806640625, -1321.603271484375, 0.),
        (1427.4986572265625, -1320.3194580078125, 0.),
        (1430.1204833984375, -1319.3392333984375, 0.),
        (1466.3984375, -1308.2322998046875, 0.),
        (1470.46728515625, -1307.0625, 0.)), width=7.)


def bergamo_tight_loop():
    return Network('bergamo_tight_loop', (
        Section(1087, 1087, ((1163.4005126953125, 826.5556640625, 0.),
            (1157.5948486328125, 761.5551147460938, 0.),
            (1157.651123046875, 751.2581176757812, 0.),
            (1159.13818359375, 743.844482421875, 0.)), width=7.),
        Section(1091, 1091, ((1159.13818359375, 743.844482421875, 0.),
            (1155.833251953125, 740.1814575195312, 0.),
            (1156.64306640625, 735.3058471679688, 0.),
            (1160.952392578125, 732.902099609375, 0.),
            (1165.5255126953125, 734.7731323242188, 0.),
            (1166.9090576171875, 739.515625, 0.),
            (1164.0616455078125, 743.5559692382812, 0.),
            (1159.13818359375, 743.844482421875, 0.)), width=7.)),
        (Intersection(2181, (End(1087, 1), End(1091, 0), End(1091, 1)), 2.,
                      (2., 2., 2.), ('STYLE',)*3),))


class DenseControlTests(unittest.TestCase):
    def test_real_dense_bend_is_certified_without_changing_authored_controls(self):
        source = bergamo_dense_bend()
        with self.assertRaisesRegex(ValueError, 'more source length'):
            Path(source, _coalescing=False)
        path = Path(source)
        self.assertIs(path.section, source)
        self.assertEqual(path.coalescing['source_control_count'], len(source.points))
        self.assertLess(len(path.coalescing['retained_indices']), len(source.points))
        np.testing.assert_array_equal(path.source_points, np.array(source.points)[:, :2])
        self.assertTrue(within_polyline_envelope(path.coefficients, path.source_points, path.deviation_limit))
        self.assertLess(path.curvature_bound * source.width*.5, .98)
        for end in (0, 1):
            original = np.array(source.points)[:, :2]
            delta = original[1]-original[0] if end == 0 else original[-1]-original[-2]
            np.testing.assert_allclose(path.frame(0. if end == 0 else path.length)[1],
                                       delta/np.linalg.norm(delta), atol=1e-12)

    def test_dense_retry_is_invariant_under_reversal_and_rigid_motion(self):
        source = bergamo_dense_bend()
        first = Path(source)
        reversed_path = Path(replace(source, points=source.points[::-1]))
        stations = np.linspace(0., first.length, 1001)
        np.testing.assert_allclose(first.point(stations),
            reversed_path.point(reversed_path.length-stations), atol=1e-8)
        rotation = np.array(((np.cos(.37), -np.sin(.37)), (np.sin(.37), np.cos(.37))))
        points = np.array(source.points)
        points[:, :2] = points[:, :2] @ rotation.T + (27., -14.)
        moved = Path(replace(source, points=tuple(map(tuple, points))))
        np.testing.assert_allclose(moved.point(stations), first.point(stations) @ rotation.T + (27., -14.), atol=1e-8)

    def test_exact_and_previously_accepted_curves_are_unchanged(self):
        source = bergamo_dense_bend()
        exact = Path(replace(source, interpolation=0))
        self.assertIsNone(exact.coalescing)
        np.testing.assert_array_equal(exact.source_points, exact.points)
        gentle = Section(1, 1, ((0., 0., 0.), (10., 1., 0.), (20., 2.5, 0.), (30., 4.5, 0.)))
        before, after = Path(gentle, _coalescing=False), Path(gentle)
        self.assertIsNone(after.coalescing)
        np.testing.assert_array_equal(after.coefficients, before.coefficients)
        np.testing.assert_array_equal(after.parameters, before.parameters)
        np.testing.assert_array_equal(after.stations, before.stations)

    def test_source_reversal_and_sharp_corner_cannot_be_erased(self):
        for points in (
            ((0., 0.), (0., -1.), (1., -1.)),
            ((0., 0.), (10., 0.), (2., .1), (20., .2)),
            ((0., 0.), (10., 0.), (10., 10.), (0., 10.), (0., 0.))):
            values = np.asarray(points)
            self.assertEqual(coalesced_source_indices(values, 100.), tuple(range(len(points))))
        source = Section(1, 1, ((0., 0., 0.), (0., -1., 0.), (1., -1., 0.)), width=5.)
        with self.assertRaises(ValueError):
            Path(source)

    def test_envelope_proof_rejects_an_excursion_between_valid_endpoints(self):
        # Cubic endpoints are on the source line, while its interior is far
        # outside. A control/sample endpoint-only check would falsely pass.
        coefficients = np.array((((0., 0.), (0., 300.), (30., -300.), (-20., 0.)),))
        source = np.array(((0., 0.), (10., 0.)))
        self.assertFalse(within_polyline_envelope(coefficients, source, 1.))
        self.assertTrue(within_polyline_envelope(coefficients, source, 100.1))

    def test_fixed_shared_endpoint_neighbor_is_protected(self):
        points = np.array(((0., 0.), (.01, 0.), (1., .001), (2., 0.)))
        self.assertEqual(coalesced_source_indices(points, .01), (0, 3))
        self.assertIn(1, coalesced_source_indices(points, .01, protected=(1,)))

    def test_short_valid_bend_uses_a_certified_radius_reserve_and_valid_quads(self):
        source = Section(2227, 2227, ((-1223.2373046875, -362.9956970214844, 0.),
            (-1223.6170654296875, -364.5513916015625, 0.),
            (-1122.204345703125, -390.962646484375, 0.)), width=2.)
        with self.assertRaisesRegex(ValueError, 'more source length'):
            Path(source, _coalescing=False)
        network = Network('short-valid-bend', (source,))
        update = PythonMeshBackend(network.settings).prepare(network)
        path = update.products[Key('path', 2227)]
        road = update.products[Key('road', 2227)]
        self.assertGreaterEqual(path.radius_reserve_fallback['minimum_radius_factor'], 1.35)
        self.assertIsNone(path.coalescing)
        self.assertLess(path.curvature_bound * source.width*.5, .98)
        self.assertTrue(within_polyline_envelope(path.coefficients, path.source_points, path.deviation_limit))
        self.assertLessEqual(road.mesh.validate()['maximum_edge_ratio'], 8.)
        self.assertEqual(network.sections[0], source)

    def test_tight_returning_loop_has_uniform_quads_at_published_precision(self):
        network = bergamo_tight_loop()
        update = PythonMeshBackend(network.settings).prepare(network)
        path = update.products[Key('path', 1091)]
        road = update.products[Key('road', 1091)]
        self.assertEqual(path.radius_reserve_fallback['minimum_radius_factor'], 1.3)
        published = road.mesh.with_positions(publication_coordinates(road.mesh.xyz, 1., physical=True))
        self.assertLessEqual(published.validate()['maximum_edge_ratio'], 8.)
        self.assertLess(len(published.quads), 1600)
        self.assertGreater(road.metadata['start'], 0.)
        self.assertLess(road.metadata['end'], path.length)

    def test_offset_chord_certificate_covers_both_curbs_and_interior_offsets(self):
        source = bergamo_tight_loop().section_index[1091]
        path = Path(source, returning=True)
        stations = path.samples(0., path.length, .7)
        for a, b in zip(stations, stations[1:]):
            bound = path._offset_chord_bound(a, b)
            self.assertLessEqual(bound, .005)
            query = np.linspace(a, b, 31)
            for offset in (-3.5, -1., 0., 2., 3.5):
                points = path.point(query, offset)
                delta = points[-1]-points[0]
                parameter = np.clip((points-points[0]) @ delta / (delta @ delta), 0., 1.)
                error = np.linalg.norm(points-points[0]-parameter[:, None]*delta, axis=1).max()
                self.assertLessEqual(error, bound+1e-9)
        with self.assertRaisesRegex(ValueError, 'allocation budget'):
            path.samples(0., path.length, .7, max_samples=5)

    def test_short_shared_leg_uses_available_length_without_losing_safe_radius(self):
        cases = ((7., ((0., 7.7, 0.), (0., 0., 0.)), ((0., 0., 0.), (175., 0., 0.))),
                 (5., ((0., 10., 0.), (0., 0., 0.)), ((0., 0., 0.), (.3, -1., 0.))))
        for width, first, last in cases:
            with self.subTest(width=width):
                network = Network('short-shared', (Section(1, 1, first, width=width),
                    Section(2, 1, last, width=width)),
                    connections=(Connection(1, (End(1, 1), End(2, 0))),))
                update = PythonMeshBackend(network.settings).prepare(network)
                a, b = (update.products[Key('path', sid)] for sid in (1, 2))
                self.assertTrue(a.endpoint_bends)
                self.assertTrue(b.endpoint_bends)
                np.testing.assert_allclose(a.frame(a.length)[0], b.frame(0.)[0], atol=1e-10)
                np.testing.assert_allclose(a.frame(a.length)[1], b.frame(0.)[1], atol=1e-10)
                self.assertGreaterEqual(a.endpoint_bends[1].radius, 1.5*width*.5-1e-10)
                for sid, path in ((1, a), (2, b)):
                    mesh = update.products[Key('road', sid)].mesh
                    published = mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True))
                    self.assertLessEqual(published.validate()['maximum_edge_ratio'], 8.)
                    self.assertLess(path.curvature_bound*width*.5, .98)

    def test_adjacent_dense_corner_keeps_original_corridor_and_published_quality(self):
        source = Section(1554, 1554, ((1708.05029296875, 857.9842529296875, 0.),
            (1793.4747314453125, 882.1200561523438, 0.),
            (1797.02099609375, 883.1007080078125, 0.),
            (1799.206787109375, 882.333251953125, 0.),
            (1799.2857666015625, 878.3257446289062, 0.),
            (1798.836181640625, 872.9043579101562, 0.)), width=5.)
        network = Network('dense-adjacent-corner', (source,))
        update = PythonMeshBackend(network.settings).prepare(network)
        path = update.products[Key('path', 1554)]
        self.assertEqual(path.coalescing['retained_indices'], (0, 3, 5))
        self.assertEqual(path.deviation_limit, 5.)
        self.assertTrue(within_polyline_envelope(path.coefficients, path.source_points, 5.))
        mesh = update.products[Key('road', 1554)].mesh
        published = mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True))
        self.assertLessEqual(published.validate()['maximum_edge_ratio'], 8.)


if __name__ == '__main__':
    unittest.main()
