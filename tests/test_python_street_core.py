from dataclasses import replace
import json
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section, Network, Settings, Placement, Key, fingerprint
from tested_python_streets.model import Connection, End
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.paths import Path
from tested_python_streets.providers import PathProvider, CarriagewayProvider, FurnitureProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.updates import UpdateCoordinator, DependencyIndex, Task


def source():
    return Network('test-network',
        (Section(1, 100, ((0., 0., 0.), (20., 0., 0.), (40., 0., 0.))),
         Section(2, 200, ((0., 30., 0.), (20., 30., 0.)))),
        furniture=(Placement(81, 1, 18.),))


def updater():
    return UpdateCoordinator((PathProvider(), ProfileProvider(), CarriagewayProvider(), FurnitureProvider()), PatternLibrary())


class PythonStreetCoreTests(unittest.TestCase):
    def test_source_roundtrip_preserves_identity_and_fingerprints(self):
        network = source()
        reopened = Network.from_dict(json.loads(json.dumps(network.to_dict())))
        self.assertEqual(network, reopened)
        self.assertEqual(fingerprint(network), fingerprint(reopened))
        self.assertEqual(Key.parse('section:1'), Key('section', 1))

    def test_invalid_or_duplicate_source_identity_is_rejected(self):
        network = source()
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            replace(network, sections=(network.sections[0], network.sections[0]))
        for width in (0., -1., float('nan')):
            with self.assertRaises(ValueError):
                replace(network.sections[0], width=width)
        with self.assertRaises(ValueError):
            Key('section', True)

    def test_pattern_connectivity_is_direct_immutable_and_reused(self):
        library = PatternLibrary(512)
        grid = library.grid(2, 2)
        self.assertEqual(grid.quads.shape, (4, 4))
        self.assertIs(library.grid(2, 2), grid)
        with self.assertRaises(ValueError):
            grid.quads[0, 0] = 99
        library.grid(3, 3)
        self.assertLessEqual(library.bytes, 512)
        self.assertEqual(library.stats()['hits'], 1)
        self.assertEqual(library.stats()['patterns'], 1)

    def test_exact_path_retains_controls_and_mitered_offset_lines(self):
        section = Section(4, 4, ((0., 0., 0.), (10., 0., 0.), (10., 10., 0.)), interpolation=0)
        path = Path(section)
        self.assertEqual(path.length, 20.)
        np.testing.assert_allclose(path.frame(path.knots)[0], np.array(section.points)[:, :2])
        np.testing.assert_allclose(path.point(10., 3.), (13., -3.))
        np.testing.assert_allclose(path.point(5., 3.), (6.5, -3.))
        np.testing.assert_allclose(path.point(15., 3.), (13., 3.5))

    def test_smooth_cubic_passes_through_controls_and_has_metric_stations(self):
        section = Section(1, 1, ((0., 0., 0.), (20., 3., 0.), (40., 0., 0.)))
        path = Path(section)
        np.testing.assert_allclose(path.jet(path.knots)[0], np.array(section.points)[:, :2], atol=1e-12)
        self.assertGreater(path.length, path.knots[-1])
        stations = np.linspace(0., path.length, 17)
        parameters = path.parameter(stations)
        lengths = [sum(path._integral(np.asarray(path.knots[i]),
                   np.asarray(min(t, path.knots[i + 1]))).item()
                   for i in range(len(path.chords)) if t > path.knots[i]) for t in parameters]
        np.testing.assert_allclose(lengths, stations, atol=1e-7)
        self.assertIn(path.certificate_intervals, (8, 16, 32, 64))

    def test_offset_fold_or_reversal_is_rejected_without_changing_inputs(self):
        for points, mode in ((((0., 0., 0.), (1., 0., 0.), (1., 1., 0.)), 3),
                             (((0., 0., 0.), (10., 0., 0.), (0., .1, 0.)), 0)):
            section = Section(1, 1, points, width=8., interpolation=mode)
            original = fingerprint(section)
            with self.assertRaises(ValueError):
                Path(section)
            self.assertEqual(fingerprint(section), original)

    def test_update_rebuilds_only_changed_section_and_its_furniture(self):
        network, engine = source(), updater()
        initial = engine.prepare(network)
        engine.commit(initial)
        untouched = engine.products[Key('road', 2)]
        changed = replace(network, sections=(replace(network.sections[0], width=8.), network.sections[1]))
        update = engine.prepare(changed)
        self.assertEqual(set(update.rebuilt), {Key('path', 1), Key('road', 1), Key('placement', 81)})
        self.assertIs(update.products[Key('road', 2)], untouched)
        before = initial.products[Key('placement', 81)]
        after = update.products[Key('placement', 81)]
        self.assertAlmostEqual(after.position[0], before.position[0])
        self.assertAlmostEqual(after.position[1] - before.position[1], -1.)
        fresh = updater().prepare(changed)
        for key in (Key('road', 1), Key('road', 2)):
            np.testing.assert_array_equal(update.products[key].mesh.xyz, fresh.products[key].mesh.xyz)
            np.testing.assert_array_equal(update.products[key].mesh.quads, fresh.products[key].mesh.quads)

    def test_noop_does_not_rebuild_and_source_movement_repositions_anchor(self):
        network, engine = source(), updater()
        engine.commit(engine.prepare(network))
        self.assertEqual(engine.prepare(network).rebuilt, ())
        moved = replace(network.sections[0], points=((0., 4., 0.), (20., 4., 0.), (40., 4., 0.)))
        update = engine.prepare(replace(network, sections=(moved, network.sections[1])))
        self.assertAlmostEqual(update.products[Key('placement', 81)].position[1], .4)

    def test_failed_prepare_and_stale_commit_do_not_advance_built_state(self):
        network, engine = source(), updater()
        engine.commit(engine.prepare(network))
        manifest = engine.manifest()
        invalid = replace(network, furniture=(replace(network.furniture[0], station=60.),))
        with self.assertRaisesRegex(ValueError, 'outside'):
            engine.prepare(invalid)
        self.assertEqual(manifest, engine.manifest())
        pending = engine.prepare(network)
        engine.commit(engine.prepare(network))
        with self.assertRaisesRegex(ValueError, 'superseded'):
            engine.commit(pending)

    def test_persistence_and_removed_output_have_explicit_lifecycles(self):
        network, engine = source(), updater()
        engine.commit(engine.prepare(network))
        reopened = updater()
        reopened.restore_manifest(json.loads(json.dumps(engine.manifest())))
        update = reopened.prepare(network)
        self.assertEqual(len(update.rebuilt), 7)
        removed = engine.prepare(replace(network, furniture=()))
        self.assertEqual(removed.removed, (Key('placement', 81),))
        self.assertIn(Key('furniture', 81), removed.changed)

    def test_dependency_cycle_or_missing_input_is_a_diagnostic(self):
        a, b = Key('a', 1), Key('b', 2)
        with self.assertRaisesRegex(ValueError, 'cycle'):
            DependencyIndex((Task(a, (b,), None), Task(b, (a,), None)))
        graph = DependencyIndex((Task(b, (a,), None),))
        self.assertEqual(graph.affected((a,)), {b})

    def test_flat_profile_has_zero_edges_and_preserves_crown_bound(self):
        result = updater().prepare(source())
        part = result.products[Key('road', 1)]
        self.assertLess(part.metadata['quality']['maximum_edge_ratio'], 2.)
        self.assertAlmostEqual(part.mesh.xyz[:, 2].max(), .8)
        self.assertTrue(np.all(part.mesh.xyz[part.mesh.seams[0][[0, -1]], 2] == 0.))

    def test_new_spatial_interaction_is_detected_without_an_old_dependency(self):
        network = source()
        engine = PythonMeshBackend(network.settings)
        engine.commit(engine.prepare(network))
        manifest = engine.manifest()
        moved = replace(network.sections[1], points=((0., .5, 0.), (20., .5, 0.)))
        with self.assertRaisesRegex(ValueError, 'interiors overlap'):
            engine.prepare(replace(network, sections=(network.sections[0], moved)))
        self.assertEqual(engine.manifest(), manifest)

    def test_smooth_output_sampling_respects_physical_offset_curve_accuracy(self):
        section = Section(1, 1, ((0., 0., 0.), (20., 3., 0.), (40., 0., 0.)))
        path = Path(section, tolerance=.005)
        stations = path.samples(2.123, path.length - 1.234, .7)
        query = np.linspace(stations[0], stations[-1], 2001)
        for offset in (-3., 3.):
            vertices = path.point(stations, offset)
            actual = path.point(query, offset)
            at = np.clip(np.searchsorted(stations, query, side='right') - 1, 0, len(stations) - 2)
            a, b = vertices[at], vertices[at + 1]
            d = b - a
            parameter = np.clip(np.sum((actual - a) * d, axis=1) / np.sum(d * d, axis=1), 0., 1.)
            distance = np.linalg.norm(actual - a - parameter[:, None] * d, axis=1)
            self.assertLessEqual(float(distance.max()), .005)

    def test_explicit_section_seam_has_matching_profiles_without_duplicate_surface(self):
        first = Section(1, 1, ((-20., 0., 0.), (0., 0., 0.)))
        second = Section(2, 1, ((0., 0., 0.), (20., 3., 0.)))
        network = Network('seam', (first, second), connections=(Connection(1, (End(1, 1), End(2, 0))),))
        update = PythonMeshBackend(network.settings).prepare(network)
        a, b = update.products[Key('road', 1)].mesh, update.products[Key('road', 2)].mesh
        np.testing.assert_allclose(a.xyz[a.seams[1]], b.xyz[b.seams[0]], atol=1e-8)
        self.assertEqual({k.kind for k in update.products}, {'path', 'profile', 'road', 'seam_layout'})

    def test_published_geometry_cache_can_be_released_and_pattern_budget_reduced(self):
        network = source()
        engine = PythonMeshBackend(network.settings)
        engine.commit(engine.prepare(network))
        engine.release_published_meshes()
        self.assertIsNone(engine.updates.products[Key('road', 1)].mesh)
        changed = replace(network, settings=replace(network.settings, pattern_cache_bytes=0))
        update = engine.prepare(changed)
        self.assertLessEqual(engine.patterns.bytes, 0)
        self.assertIsNotNone(update.products[Key('road', 1)].mesh)


if __name__ == '__main__':
    unittest.main()
