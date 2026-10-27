"""Source locality, continuous offset safety and update contracts at sparse bends."""
from dataclasses import replace
import json
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section, Network, Settings, Placement, Key, Connection, End
from tested_python_streets.paths import Path
from tested_python_streets.backend import PythonMeshBackend


def distance_to_polyline(query, points):
    points = np.asarray(points)[:, :2]
    delta = np.diff(points, axis=0)
    t = np.clip(np.sum((query[:, None] - points[None, :-1]) * delta[None], axis=2)
                / np.sum(delta * delta, axis=1)[None], 0., 1.)
    return np.linalg.norm(query[:, None] - points[None, :-1] - t[..., None] * delta[None], axis=2).min(axis=1)


def courtyard():
    return Section(1, 1, ((0., 0., 0.), (0., -80., 0.), (120., -80., 0.), (120., 0., 0.)), width=5.)


class LocalSmoothingTests(unittest.TestCase):
    def test_sparse_courtyard_has_local_rounds_and_exact_straight_runs(self):
        section = courtyard()
        original = section.points
        path = Path(section)
        self.assertEqual(len(path.local_bends), 2)
        self.assertEqual(section.points, original)
        points = path.jet(np.linspace(0., path.knots[-1], 12001))[0]
        self.assertLessEqual(float(distance_to_polyline(points, section.points).max()), path.deviation_limit)
        for point in ((0., -40.), (60., -80.), (120., -40.)):
            self.assertLess(float(np.linalg.norm(points - point, axis=1).min()), .03)
        for index in range(len(path.chords)):
            first, last = path.knots[index:index + 2]
            q = path.jet(np.linspace(first, last, 15))[0]
            if abs(q[:, 0]).max() < 1e-8:
                self.assertLess(abs(q[:, 0]).max(), 1e-12)
        self.assertLess(path.length, sum(np.linalg.norm(np.diff(np.asarray(section.points)[:, :2], axis=0), axis=1)))

    def test_tangent_and_metric_station_continuity_survive_reverse_and_rotation(self):
        original = courtyard()
        path = Path(original)
        reversed_path = Path(replace(original, points=original.points[::-1]))
        stations = np.linspace(0., path.length, 251)
        np.testing.assert_allclose(path.point(stations), reversed_path.point(reversed_path.length - stations), atol=1e-8)
        parameters = np.linspace(0., path.knots[-1], 301)
        np.testing.assert_allclose(path.parameter(path.station(parameters)), parameters, atol=1e-8)
        for knot in path.knots[1:-1]:
            a, b = path.jet(np.array((knot - 1e-8, knot + 1e-8)))[1]
            np.testing.assert_allclose(a, b, atol=2e-8)
        angle = .731
        rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        source = np.asarray(original.points)
        source[:, :2] = source[:, :2] @ rotation.T + (913., -607.)
        moved = Path(replace(original, points=tuple(map(tuple, source))))
        np.testing.assert_allclose(moved.point(stations), path.point(stations) @ rotation.T + (913., -607.), atol=1e-8)

    def test_actual_short_long_source_spans_stay_inside_the_envelope(self):
        sections = (
            Section(284, 608, ((-126.37627410888672, -91.51412200927734, 0.),
                (-131.31402587890625, -72.43383026123047, 0.),
                (-133.0684051513672, -69.36136627197266, 0.),
                (-300.9466552734375, -40.76707077026367, 0.)), width=5.),
            Section(421, 615, ((337.8796081542969, -158.2178497314453, 0.),
                (337.0125427246094, -158.5184783935547, 0.),
                (239.9792938232422, -143.7743682861328, 0.),
                (235.23046875, -140.00086975097656, 0.)), width=2.))
        for section in sections:
            with self.subTest(section=section.id):
                path = Path(section)
                self.assertTrue(path.local_bends)
                xy = path.jet(np.linspace(0., path.knots[-1], 12001))[0]
                self.assertLessEqual(float(distance_to_polyline(xy, section.points).max()), path.deviation_limit)
                parameter = np.linspace(0., path.knots[-1], 12001)
                unused, velocity, acceleration = path.jet(parameter)
                curvature = abs(velocity[:, 0] * acceleration[:, 1] - velocity[:, 1] * acceleration[:, 0]) / np.linalg.norm(velocity, axis=1) ** 3
                self.assertLess(float(curvature.max()) * (abs(section.offset) + section.width * .5), .98)

    def test_short_bends_report_infeasibility_and_exact_mode_keeps_vertices(self):
        section = replace(courtyard(), points=((0., 0., 0.), (0., -1., 0.), (1., -1., 0.)))
        with self.assertRaisesRegex(ValueError, 'more source length|fold'):
            Path(section)
        exact = Path(replace(courtyard(), interpolation=0))
        np.testing.assert_allclose(exact.point(exact.knots), np.asarray(courtyard().points)[:, :2])

    def test_tiny_turn_next_to_a_large_bend_has_no_absolute_speed_cutoff(self):
        section=Section(1,1,((0.,0.,0.),(0.,40.,0.),(.0001,80.,0.),(80.,80.,0.)),width=2.)
        path=Path(section)
        self.assertTrue(path.localized)
        self.assertGreater(path.length,100.)
        self.assertLess(path.curvature_bound, .98)

    def test_endpoint_frames_and_native_quad_dimensions_remain_valid(self):
        section = courtyard()
        angle = .35
        endpoints = {0: np.array((np.sin(angle), -np.cos(angle))),
                     1: np.array((np.sin(angle), np.cos(angle)))}
        path = Path(section, endpoint_directions=endpoints)
        np.testing.assert_allclose(path.frame(0.)[1], endpoints[0], atol=1e-12)
        np.testing.assert_allclose(path.frame(path.length)[1], endpoints[1], atol=1e-12)
        stations = path.samples(0., path.length, .7)
        xy = path.point(stations[:, None], np.linspace(-2.5, 2.5, 9)[None]).astype(np.float32).astype(float)
        faces = np.stack((xy[:-1, :-1], xy[:-1, 1:], xy[1:, 1:], xy[1:, :-1]), axis=2)
        lengths = np.linalg.norm(np.roll(faces, -1, axis=2) - faces, axis=-1)
        self.assertLess(float((lengths.max(axis=2) / lengths.min(axis=2)).max()), 8.)

    def test_sampling_refines_bends_without_densifying_long_straights(self):
        path = Path(courtyard())
        stations = path.samples(0., path.length, .7)
        self.assertLess(len(stations), 800)
        self.assertGreater(np.count_nonzero(np.diff(stations) > .4), 300)
        query = np.linspace(0., path.length, 24001)
        at = np.clip(np.searchsorted(stations, query, side='right') - 1, 0, len(stations) - 2)
        for offset in (-2.5, 2.5):
            samples, actual = path.point(stations, offset), path.point(query, offset)
            a, b = samples[at], samples[at + 1]
            delta = b - a
            t = np.clip(np.sum((actual - a) * delta, axis=1) / np.sum(delta * delta, axis=1), 0., 1.)
            error = np.linalg.norm(actual - a - t[:, None] * delta, axis=1)
            self.assertLessEqual(float(error.max()), .005)
        with self.assertRaisesRegex(ValueError, 'allocation budget'):
            path.samples(0., path.length, .7, max_samples=5)

    def test_almost_collinear_endpoint_does_not_consume_a_short_source_leg(self):
        section = Section(1, 1, ((0., 0., 0.), (6., 0., 0.), (126., 120., 0.)), width=5.)
        path = Path(section, endpoint_directions={0: np.array((1., .00001))})
        self.assertTrue(path.local_bends)
        np.testing.assert_allclose(path.frame(0.)[1], np.array((1., .00001)) / np.hypot(1., .00001), atol=1e-12)

    def test_courtyard_split_into_two_control_sections_gets_the_same_envelope(self):
        source = Network('split-courtyard', (
            Section(1, 1, ((0., 0., 0.), (0., -80., 0.)), width=2.),
            Section(2, 1, ((0., -80., 0.), (120., -80., 0.)), width=2.),
            Section(3, 1, ((120., -80., 0.), (120., 0., 0.)), width=2.)),
            connections=(Connection(1, (End(1, 1), End(2, 0))), Connection(2, (End(2, 1), End(3, 0)))),
            settings=Settings(exaggeration=1.))
        result = PythonMeshBackend(source.settings).prepare(source)
        for section in source.sections:
            path = result.products[Key('path', section.id)]
            self.assertTrue(path.localized)
            xy = path.point(np.linspace(0., path.length, 2001))
            self.assertLessEqual(float(distance_to_polyline(xy, section.points).max()), path.deviation_limit)
            self.assertLess(result.products[Key('road', section.id)].mesh.validate()['maximum_edge_ratio'], 8.)
        for a, b in ((1, 2), (2, 3)):
            first, second = (result.products[Key('path', sid)] for sid in (a, b))
            np.testing.assert_allclose(first.frame(first.length)[0], second.frame(0.)[0], atol=1e-10)
            np.testing.assert_allclose(first.frame(first.length)[1], second.frame(0.)[1], atol=1e-10)

    def test_split_courtyard_has_one_signed_round_instead_of_endpoint_bulges(self):
        source = Network('shared-round', (
            Section(1, 1, ((0., 0., 0.), (0., -80., 0.)), width=2.),
            Section(2, 1, ((0., -80., 0.), (120., -80., 0.)), width=2.),
            Section(3, 1, ((120., -80., 0.), (120., 0., 0.)), width=2.)),
            connections=(Connection(1, (End(1, 1), End(2, 0))), Connection(2, (End(2, 1), End(3, 0)))),
            settings=Settings(exaggeration=1.))
        result = PythonMeshBackend(source.settings).prepare(source)
        paths = [result.products[Key('path', s.id)] for s in source.sections]
        whole = Path(replace(courtyard(), width=2.))
        shift = 0.
        for path in paths:
            stations = np.linspace(0., path.length, 2001)
            np.testing.assert_allclose(path.point(stations), whole.point(stations + shift), atol=2e-8)
            shift += path.length
            unused, d, dd = path.jet(np.linspace(0., path.knots[-1], 2001))
            signed = d[:, 0] * dd[:, 1] - d[:, 1] * dd[:, 0]
            self.assertGreaterEqual(float(signed.min()), -1e-9)
        self.assertAlmostEqual(shift, whole.length, places=9)
        for a, b in zip(paths, paths[1:]):
            # The derived seam lies inside the corner, with identical signed
            # curvature on both halves, instead of hitting the sharp vertex.
            self.assertGreater(np.linalg.norm(a.point(a.length) - np.array(a.section.points[-1])[:2]), .1)
            jets = [p.jet(np.asarray(t)) for p, t in ((a, a.knots[-1]), (b, 0.))]
            curvature = [float(d[0] * dd[1] - d[1] * dd[0]) / np.linalg.norm(d) ** 3
                         for unused, d, dd in jets]
            self.assertAlmostEqual(*curvature, places=10)

    def test_shared_rounds_are_invariant_to_section_order_and_direction(self):
        from test_python_street_dimensional_bundles import reverse_sections
        source = Network('shared-reversal', (
            Section(337, 595, ((-6.827916, 126.648186, 0.), (-19.589930, 48.479679, 0.)), width=2.),
            Section(359, 595, ((-19.589930, 48.479679, 0.), (99.155823, 27.007114, 0.)), width=2.),
            Section(431, 595, ((99.155823, 27.007114, 0.), (115.892105, 108.760513, 0.)), width=2.)),
            connections=(Connection(49, (End(337, 1), End(359, 0))),
                         Connection(50, (End(359, 1), End(431, 0)))))
        baseline = PythonMeshBackend(source.settings).prepare(source)
        for ids in ({337}, {359}, {431}, {337, 359, 431}):
            with self.subTest(reversed_ids=ids):
                reversed_source = reverse_sections(source, ids)
                reversed_source = replace(reversed_source, sections=reversed_source.sections[::-1])
                result = PythonMeshBackend(source.settings).prepare(reversed_source)
                for section in source.sections:
                    a, b = [r.products[Key('path', section.id)] for r in (baseline, result)]
                    stations = np.linspace(0., a.length, 1001)
                    np.testing.assert_allclose(a.point(stations),
                        b.point(b.length - stations if section.id in ids else stations), atol=2e-8)
                    unused, d, dd = a.jet(np.linspace(0., a.knots[-1], 2001))
                    self.assertGreaterEqual(float((d[:, 0] * dd[:, 1] - d[:, 1] * dd[:, 0]).min()), -1e-8)

    def test_shared_signed_turns_and_width_update_keep_one_coherent_round(self):
        from tested_python_streets.providers import PathProvider
        from types import SimpleNamespace
        for degrees in (30., 60., 90., 120., 150., -30., -90., -150.):
            with self.subTest(turn=degrees):
                angle = np.radians(degrees)
                source = Network('shared-angles', (
                    Section(1, 1, ((-200., 0., 0.), (0., 0., 0.)), width=2.),
                    Section(2, 1, ((0., 0., 0.), (200.*np.cos(angle), 200.*np.sin(angle), 0.)), width=2.1)),
                    connections=(Connection(1, (End(1, 1), End(2, 0))),))
                paths = [PathProvider.build(SimpleNamespace(network=source), s, source.connections)
                         for s in source.sections]
                self.assertTrue(all(path.endpoint_bends for path in paths))
                np.testing.assert_allclose(paths[0].frame(paths[0].length)[0], paths[1].frame(0.)[0], atol=1e-12)
                np.testing.assert_allclose(paths[0].frame(paths[0].length)[1], paths[1].frame(0.)[1], atol=1e-12)
                for path in paths:
                    unused, d, dd = path.jet(np.linspace(0., path.knots[-1], 2001))
                    self.assertGreaterEqual(float((np.sign(degrees) * (d[:, 0]*dd[:, 1]-d[:, 1]*dd[:, 0])).min()), -1e-8)

        source = Network('shared-width', (
            Section(1, 1, ((-80., 0., 0.), (0., 0., 0.)), width=2.),
            Section(2, 1, ((0., 0., 0.), (0., 80., 0.)), width=2.),
            Section(99, 99, ((200., 0., 0.), (280., 0., 0.)), width=2.)),
            connections=(Connection(1, (End(1, 1), End(2, 0))),),
            furniture=(Placement(81, 1, 40.),))
        backend = PythonMeshBackend(source.settings)
        first = backend.prepare(source)
        backend.commit(first)
        edited = replace(source, sections=(replace(source.sections[0], width=2.1), *source.sections[1:]))
        local = backend.prepare(edited)
        full = PythonMeshBackend(edited.settings).prepare(edited)
        self.assertIs(local.products[Key('road', 99)], first.products[Key('road', 99)])
        self.assertTrue(all(local.products[Key('path', sid)].endpoint_bends for sid in (1, 2)))
        for key, part in full.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(local.products[key].mesh.xyz, part.mesh.xyz)
                np.testing.assert_array_equal(local.products[key].mesh.quads, part.mesh.quads)

    def test_local_update_matches_full_and_preserves_remote_and_furniture_dependencies(self):
        source = Network('local-smoothing', (courtyard(), Section(99, 99, ((0., 100., 0.), (120., 100., 0.)))),
            settings=Settings(exaggeration=1.), furniture=(Placement(81, 1, 45.),))
        backend = PythonMeshBackend(source.settings)
        first = backend.prepare(source)
        backend.commit(first)
        original = first.products[Key('road', 1)].mesh.xyz.copy()
        source_copy = json.dumps(source.to_dict(), sort_keys=True)
        remote = first.products[Key('road', 99)]
        section = replace(source.sections[0], width=5.5,
            points=((0., 0., 0.), (0., -82., 0.), (120., -80., 0.), (120., 0., 0.)))
        edited = replace(source, sections=(section, source.sections[1]))
        local = backend.prepare(edited)
        full = PythonMeshBackend(edited.settings).prepare(edited)
        self.assertIs(local.products[Key('road', 99)], remote)
        self.assertIn(Key('placement', 81), local.rebuilt)
        for key, value in full.products.items():
            if getattr(value, 'mesh', None) is not None:
                np.testing.assert_array_equal(local.products[key].mesh.xyz, value.mesh.xyz)
                np.testing.assert_array_equal(local.products[key].mesh.quads, value.mesh.quads)
        backend.commit(local)
        restored = backend.prepare(source)
        np.testing.assert_array_equal(restored.products[Key('road', 1)].mesh.xyz, original)
        self.assertEqual(json.dumps(source.to_dict(), sort_keys=True), source_copy)


if __name__ == '__main__':
    unittest.main()
