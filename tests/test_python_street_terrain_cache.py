"""Exact terrain constraints, bounded retention and lowering/invalidation tests."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.terrain.alignment import (
    CorridorComponent, ElevationReference, PlacementSettings, fit)
from tested_python_streets.terrain.cache import TerrainConstraintCache
from tested_python_streets.patterns import PatternLibrary
from _python_terrain_fit_reference import fit as reference_fit


class SampleCounter:
    def __init__(self, offset=0.):
        self.samples, self.offset = 0, offset

    def __call__(self, xy):
        self.samples += len(xy)
        return 3. + self.offset + .01 * xy[:, 0] + .08 * np.exp(-np.sum((xy-(7., 1.))**2,axis=1)/10.)


def components():
    pattern = PatternLibrary().grid(24, 4)
    refs, comps = [], []
    for i in range(2):
        xy = pattern.uv * (12., 3.) + (12. * i, 0.)
        z = .08 * (1. - (2. * pattern.uv[:, 1] - 1.) ** 2)
        name = f'road:{i}'
        refs.append(ElevationReference(name, np.arange(0., 12.1, 2.)))
        comps.append(CorridorComponent(name, name, np.column_stack((xy, z)),
            pattern.quads, 12. * pattern.uv[:, 0], contacts=((12. if i == 0 else 0., 1),)))
    return refs, comps


class TerrainCacheTests(unittest.TestCase):
    def exact(self, expected, actual):
        self.assertEqual(expected.report, actual.report)
        self.assertEqual(expected.coordinates.keys(), actual.coordinates.keys())
        for key in expected.coordinates:
            np.testing.assert_array_equal(expected.coordinates[key].view(np.uint64), actual.coordinates[key].view(np.uint64))
        for key in expected.curves:
            for name in ('stations', 'heights', 'derivatives'):
                np.testing.assert_array_equal(getattr(expected.curves[key], name), getattr(actual.curves[key], name))

    def test_cold_and_warm_constraints_match_the_original_fit_without_any_warm_rays(self):
        refs, comps = components()
        expected = reference_fit(refs, comps, SampleCounter())
        cache, terrain = TerrainConstraintCache(), SampleCounter()
        self.exact(expected, fit(refs, comps, terrain, cache=cache, terrain_token='terrain-1'))
        before = terrain.samples
        self.exact(expected, fit(refs, comps, terrain, cache=cache, terrain_token='terrain-1'))
        self.assertEqual(terrain.samples, before)
        self.assertEqual(cache.hits, len(comps))

    def test_geometry_controls_sampling_and_terrain_revisions_invalidate_exact_inputs(self):
        refs, comps = components()
        cache, terrain = TerrainConstraintCache(), SampleCounter()
        base = fit(refs, comps, terrain, cache=cache, terrain_token='terrain-1')
        xyz = comps[0].xyz.copy()
        xyz[20:23, 1] += .05
        changed = [replace(comps[0], xyz=xyz), comps[1]]
        self.exact(reference_fit(refs, changed, SampleCounter()),
                   fit(refs, changed, terrain, cache=cache, terrain_token='terrain-1'))
        self.assertEqual(cache.hits, 1)
        closer = [replace(r, stations=np.arange(0., 12.1, 1.)) for r in refs]
        self.exact(reference_fit(closer, comps, SampleCounter()),
                   fit(closer, comps, terrain, cache=cache, terrain_token='terrain-1'))
        self.assertEqual(cache.hits, 1)
        settings = PlacementSettings(sample_spacing=.4, clearance=.03)
        self.exact(reference_fit(refs, comps, SampleCounter(), settings),
                   fit(refs, comps, terrain, settings, cache=cache, terrain_token='terrain-1'))
        lowered = SampleCounter(-1.)
        actual = fit(refs, comps, lowered, cache=cache, terrain_token='terrain-2')
        self.exact(reference_fit(refs, comps, SampleCounter(-1.)), actual)
        self.assertGreater(lowered.samples, 0)
        for key in base.coordinates:
            np.testing.assert_allclose(actual.coordinates[key][:, 2], base.coordinates[key][:, 2] - 1., atol=1e-12)

    def test_grade_contacts_and_lift_are_recomputed_even_on_cache_hits(self):
        refs, comps = components()
        cache = TerrainConstraintCache()
        fit(refs, comps, SampleCounter(), cache=cache, terrain_token='terrain-1')
        settings = PlacementSettings(maximum_grade=.02)
        terrain = SampleCounter()
        changed = [replace(c, contacts=()) for c in comps]
        self.exact(reference_fit(refs, changed, SampleCounter(), settings),
                   fit(refs, changed, terrain, settings, cache=cache, terrain_token='terrain-1'))
        self.assertEqual(terrain.samples, 0)
        with self.assertRaisesRegex(ValueError, 'lifts'):
            fit(refs, comps, terrain, PlacementSettings(maximum_lift=.001),
                cache=cache, terrain_token='terrain-1')
        with self.assertRaisesRegex(ValueError, 'cover'):
            fit(refs, comps, lambda xy: np.full(len(xy), np.nan), cache=cache, terrain_token='missing')

    def test_unknown_terrain_revision_and_zero_budget_keep_uncached_behavior(self):
        refs, comps = components()
        for cache, token in ((TerrainConstraintCache(), None), (TerrainConstraintCache(max_bytes=0), 'terrain-1')):
            for offset in (0., -1.):
                terrain = SampleCounter(offset)
                self.exact(reference_fit(refs, comps, SampleCounter(offset)),
                           fit(refs, comps, terrain, cache=cache, terrain_token=token))
                self.assertGreater(terrain.samples, 0)
                self.assertEqual(cache.bytes, 0)

    def test_cache_budget_eviction_and_mutable_alias_isolation(self):
        cache = TerrainConstraintCache(max_bytes=1800, max_entries=2)
        floor, ground = np.array((1., 2.)), np.arange(10.)
        for i in range(5):
            cache.put(bytes((i,)), floor, ground, 90)
            self.assertLessEqual(cache.bytes, cache.max_bytes)
            self.assertLessEqual(len(cache.entries), 2)
        self.assertEqual(cache.evictions, 3)
        floor[:] = ground[:] = -99.
        saved = cache.get(bytes((4,)))
        np.testing.assert_array_equal(saved.floors, (1., 2.))
        np.testing.assert_array_equal(saved.ground, np.arange(10.))
        self.assertFalse(saved.ground.flags.writeable)
        self.assertFalse(saved.floors.flags.writeable)


if __name__ == '__main__':
    unittest.main()
