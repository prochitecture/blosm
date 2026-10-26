"""Actual provider/update contracts for curved overlapping Bundle members."""
from dataclasses import replace
from functools import partial
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from _python_angular_fixtures import angular_network as make_angular_network


angular_network = partial(make_angular_network, model)


def buffers(update):
    return {key:(part.mesh.xyz.copy(), part.mesh.quads.copy()) for key, part in update.products.items()
            if getattr(part, 'mesh', None) is not None}


class PythonAngularBundleTests(unittest.TestCase):
    def qualify(self, network, *, backend=None):
        source = network.to_dict()
        backend = backend or PythonMeshBackend(network.settings)
        update = backend.prepare(network)
        self.assertEqual(source, network.to_dict())
        plan = update.products[model.Key('bundle_plan', 13)]
        self.assertEqual(type(plan).__name__, 'AngularBundlePlan')
        self.assertEqual(plan.source_ids, (17, 29))
        part = update.products[model.Key('bundle_pavement', 13)]
        self.assertEqual(part.metadata['source_sections'], [17, 29])
        self.assertTrue(part.metadata['all_direct_quads'])
        self.assertLessEqual(part.metadata['physical_slopes']['maximum'], plan.maximum_slope+1e-8)
        self.assertLessEqual(part.metadata['angular_longitudinal_bound'], .03+1e-8)
        self.assertTrue(part.metadata['physical_boundary_contract']['degree_two'])
        self.assertLessEqual(part.metadata['quality']['maximum_edge_ratio'], 8.)
        separator = update.products[model.Key('separator', 13)]
        if separator is not None:
            self.assertLessEqual(separator.mesh.validate()['maximum_edge_ratio'], 8.)
        self.assertEqual(update.products[model.Key('road', 17)].metadata['composite_owner'], 'bundle:13')
        for mesh_key in (model.Key('road', 17), model.Key('road', 29)):
            self.assertFalse(hasattr(update.products[mesh_key], 'mesh'))
        mesh = part.mesh
        self.assertGreaterEqual(float(mesh.xyz[:, 2].min()), plan.datum-1e-8)
        self.assertLessEqual(float(mesh.xyz[:, 2].max()), plan.datum+plan.rise*network.settings.exaggeration+1e-8)
        placement = update.products[model.Key('placement', 901)]
        self.assertAlmostEqual(placement.position[2], float(plan.height(np.asarray(placement.position[:2]))), places=10)
        return backend, update, plan

    def test_all_twelve_source_cases_use_real_provider(self):
        for sweep in (180., 225., 270.):
            for variable in (False, True):
                for mode in (0, 3):
                    with self.subTest(sweep=sweep, variable=variable, mode=mode):
                        network = angular_network(sweep, variable, mode)
                        _, update, _ = self.qualify(network)
                        pavement = update.products[model.Key('bundle_pavement', 13)]
                        self.assertEqual(pavement.metadata['islands'], int(variable))
                        separator = update.products[model.Key('separator', 13)]
                        self.assertEqual(separator is not None, variable)
                        count = len(pavement.mesh.quads)+(len(separator.mesh.quads) if variable else 0)
                        self.assertLess(count, 6000)
                        if variable:
                            self.assertTrue(pavement.metadata['physical_boundary_contract']['grass_conformity']['whole_edges_exact'])

    def test_widths_offsets_shared_rise_and_translated_datum(self):
        self.qualify(angular_network(225., True, widths=(5., 7.), offsets=(.25, -.15),
            rises=(.04, .1), datum=37.))

    def test_reversed_source_storage_and_rigid_world_transform(self):
        self.qualify(angular_network(270., True, reversed_sources=(17, 29),
            rotation=.713, translation=(1234., -5678.)))

    def test_resolution_and_zero_inspection_exaggeration(self):
        for resolution, exaggeration in ((2, 0.), (4, 10.)):
            with self.subTest(resolution=resolution):
                self.qualify(angular_network(180., True,
                    settings=model.Settings(resolution=resolution, exaggeration=exaggeration)))

    def test_explicit_updates_are_local_equal_full_and_restore(self):
        original = angular_network(225., True)
        backend, initial, _ = self.qualify(original)
        backend.commit(initial)
        saved = buffers(initial)
        remote = initial.products[model.Key('road', 99)]
        for changed in (angular_network(225., True, displacement=.1),
                        angular_network(225., True, widths=(6.1, 6.))):
            _, local, _ = self.qualify(changed, backend=backend)
            full = PythonMeshBackend(changed.settings).prepare(changed)
            self.assertEqual(buffers(local).keys(), buffers(full).keys())
            for key, pair in buffers(local).items():
                np.testing.assert_array_equal(pair[0], buffers(full)[key][0])
                np.testing.assert_array_equal(pair[1], buffers(full)[key][1])
            self.assertIs(local.products[model.Key('road', 99)], remote)
            self.assertNotIn(model.Key('placement', 999), local.rebuilt)
            self.assertIn(model.Key('placement', 901), local.rebuilt)
            backend.commit(local)
            restored = backend.prepare(original)
            for key, pair in buffers(restored).items():
                np.testing.assert_array_equal(pair[0], saved[key][0])
                np.testing.assert_array_equal(pair[1], saved[key][1])
            backend.commit(restored)

    def test_published_mesh_release_retains_source_attachment_query(self):
        network = angular_network(180., True)
        backend, first, plan = self.qualify(network)
        backend.commit(first)
        position = np.asarray(first.products[model.Key('placement', 901)].position[:2])
        before = float(plan.height(position))
        backend.release_published_meshes()
        self.assertIsNone(first.products[model.Key('bundle_pavement', 13)].mesh)
        self.assertAlmostEqual(float(plan.height(position)), before, places=12)
        furniture = tuple(replace(p, station=p.station+.1) if p.id==901 else p for p in network.furniture)
        changed = replace(network, furniture=furniture)
        local = backend.prepare(changed)
        full = PythonMeshBackend(changed.settings).prepare(changed)
        self.assertEqual(local.products[model.Key('placement', 901)], full.products[model.Key('placement', 901)])
        self.assertNotIn(model.Key('bundle_plan', 13), local.rebuilt)

    def test_different_datums_and_crossing_centerlines_are_not_silently_merged(self):
        network = angular_network(180., True)
        changed = replace(network, sections=tuple(replace(s, datum=2.) if s.id==29 else s for s in network.sections))
        with self.assertRaises(ValueError):
            PythonMeshBackend(changed.settings).prepare(changed)
        # Each individual offset remains regular; this input specifically
        # exercises the continuous cross-member contact certificate.
        points = np.asarray(network.sections[0].points).copy()
        crossing = replace(network, sections=tuple(replace(s, points=tuple(map(tuple, points))) if s.id==29 else s
            for s in network.sections))
        with self.assertRaisesRegex(ValueError, 'centerlines cross or touch'):
            PythonMeshBackend(crossing.settings).prepare(crossing)


if __name__ == '__main__':
    unittest.main(verbosity=2)
