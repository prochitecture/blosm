"""Regression for a valid translated Exact junction rejected by the old chart."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
import test_python_street_merged_bundles as merged
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.drainage import BoundaryDrainage
from tested_python_streets.entrance_slopes import EntranceSlopeGuard
from tested_python_streets.mesh import Mesh
from tested_python_streets.model import Key
from tested_python_streets.ownership import surface_key
from tested_python_streets.patterns import PatternLibrary


def source(exaggeration=1.):
    original = merged.source_for(0, exaggeration)
    return replace(original, sections=tuple(replace(section,
        points=tuple(map(tuple, (np.asarray(section.points) + (1000., -500., 0.)).astype(np.float32).astype(float))),
        width=float(np.float32(section.width)), rise=float(np.float32(section.rise)), datum=1.7)
        for section in original.sections))


class CornerSeedTests(unittest.TestCase):
    def qualify_surface(self, network, result):
        merged.MergedBundleTests().qualify(network, result)
        part = result.products[Key('junction_mesh', 210)]
        self.assertTrue(part.metadata['planar_chart']['corner_band_seed'])
        self.assertTrue(part.metadata['corner_seed_reconciliation']['exact_boundary'])
        product = result.products[surface_key(Key('intersection', 210))]
        datum = network.sections[0].datum
        physical = Mesh(part.mesh.xyz.astype(np.float32), part.mesh.quads)
        if network.settings.exaggeration == 0.:
            self.assertTrue(all(np.all(rim[:, 2] >= datum - 1e-8) for rim in product.physical_rims.values()))
            self.assertGreater(max(rim[:, 2].max() for rim in product.physical_rims.values()), datum + .07)
            return
        physical.xyz[:, 2] = datum + (physical.xyz[:, 2] - datum) / network.settings.exaggeration
        self.assertLessEqual(physical.slope_statistics()['maximum'], .1 + 1e-6)
        edges = np.sort(np.stack((physical.quads, np.roll(physical.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        boundary = np.unique(unique[counts == 1])
        fixed = np.zeros(len(physical.xyz), dtype=bool)
        fixed[boundary] = True
        _, drainage = BoundaryDrainage().fill(physical, fixed)
        self.assertLessEqual(drainage['maximum_fill'], 1e-6)
        group = product.group
        layouts = {footprint.region.owner: footprint.layout for footprint in group.footprints}
        for loop in group.loops:
            condition = BoundaryConditions().build(group, loop, layouts)
            for port in condition.ports:
                points = condition.xy[list(port.vertices)].astype(np.float32).astype(float)
                ids = []
                for point in points:
                    distances = np.linalg.norm(physical.xyz[boundary, :2] - point, axis=1)
                    at = int(np.argmin(distances))
                    self.assertEqual(distances[at], 0.)
                    ids.append(int(boundary[at]))
                self.assertEqual(len(set(ids)), len(ids))
                approach = layouts[port.owner].approaches[port.approach_index]
                q = 2. * np.asarray(port.outward_u) - 1.
                velocities = np.broadcast_to(approach.frame(port.station)[1], (len(physical.xyz), 2)).copy()
                velocities[ids] = approach.jet(port.station, approach.offset + .5 * approach.width * q)[1]
                bound = EntranceSlopeGuard().validate(physical, np.asarray(ids), velocities, physical.xyz[:, 2] - datum)
                self.assertLessEqual(bound, .03 + 1e-8)

    def test_translated_exact_network_preserves_whole_interfaces_profiles_and_drainage(self):
        network = source()
        original = network.to_dict()
        result = PythonMeshBackend(network.settings).prepare(network)
        self.qualify_surface(network, result)
        self.assertEqual(network.to_dict(), original)

    def test_zero_and_exaggerated_views_keep_physical_boundary_conditions(self):
        for exaggeration in (0., 10.):
            with self.subTest(exaggeration=exaggeration):
                network = source(exaggeration)
                result = PythonMeshBackend(network.settings).prepare(network)
                self.qualify_surface(network, result)

    def test_width_edit_matches_full_build_and_exact_restoration(self):
        network = source(10.)
        backend = PythonMeshBackend(network.settings)
        first = backend.prepare(network)
        original = merged.snapshot(first)
        backend.commit(first)
        edited = replace(network, sections=tuple(replace(section,
            width=float(np.float32(section.width + .02))) if section.id == 117 else section for section in network.sections))
        local = backend.prepare(edited)
        self.qualify_surface(edited, local)
        self.assertIs(local.products[Key('road', 99)], first.products[Key('road', 99)])
        after = merged.snapshot(local)
        backend.commit(local)
        full = backend.prepare(edited, force=True)
        current = merged.snapshot(full)
        self.assertEqual(after.keys(), current.keys())
        for key in after:
            np.testing.assert_array_equal(after[key][0], current[key][0])
            np.testing.assert_array_equal(after[key][1], current[key][1])
        backend.commit(full)
        restored = merged.snapshot(backend.prepare(network))
        self.assertEqual(original.keys(), restored.keys())
        for key in original:
            np.testing.assert_array_equal(original[key][0], restored[key][0])
            np.testing.assert_array_equal(original[key][1], restored[key][1])

    def test_staged_pattern_retains_complete_boundary_profiles_and_bounded_cache(self):
        library = PatternLibrary(max_bytes=128 * 1024)
        pattern = library.staged_strip(40, 4, 8, middle=24, end_rows=4)
        for value, count in ((0., 4), (1., 8)):
            np.testing.assert_allclose(np.sort(pattern.uv[pattern.uv[:, 0] == value, 1]), np.linspace(0., 1., count + 1), atol=1e-14)
        for value in (0., 1.):
            np.testing.assert_allclose(np.sort(pattern.uv[pattern.uv[:, 1] == value, 0]), np.linspace(0., 1., 41), atol=1e-14)
        self.assertIs(library.staged_strip(40, 4, 8, middle=24, end_rows=4), pattern)
        self.assertFalse(pattern.uv.flags.writeable)
        self.assertFalse(pattern.quads.flags.writeable)
        self.assertLessEqual(library.bytes, library.max_bytes)
        with self.assertRaisesRegex(ValueError, 'positive counts'):
            library.staged_strip(40, 4, 8, middle=24, end_rows=20)
        with self.assertRaisesRegex(ValueError, 'vertex budget'):
            library.staged_strip(100000, 4, 8, middle=24, end_rows=4)
