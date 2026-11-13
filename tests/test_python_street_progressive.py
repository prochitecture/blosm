"""Numerical checkpoints preserve results and never accept pending work."""
from copy import deepcopy
from dataclasses import replace
import unittest
from unittest.mock import patch

import numpy as np
import _python_street_package
from _python_street_fixtures import clipped_fixture_source
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.mesh import Mesh
from tested_python_streets.model import Key, Network, Section
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.providers import Part
from tested_python_streets.source_support import ClippedApproachSupport
from tested_python_streets.updates import UpdateCoordinator, Task, drain_steps
from tested_python_streets.uv import DEFAULT_UV_MAP
from tested_python_streets.validation import validate_overlaps, validate_overlaps_steps


def network():
    return Network('progressive-numerical', (
        Section(2, 2, ((0., 20., 0.), (8., 20., 0.))),
        Section(1, 1, ((0., 0., 0.), (8., 0., 0.)))))


def quad(x=0.):
    return Mesh(np.array(((x, 0., 0.), (x + 1., 0., 0.),
                          (x + 1., 1., 0.), (x, 1., 0.))), np.array(((0, 1, 2, 3),)))


class Parts:
    source_kinds = frozenset(('section',))

    def __init__(self, build=None):
        self.build = build or (lambda context, section:
            Part(Key('road', section.id), quad(section.id * 3.), {'source': section.id}))

    def tasks(self, source):
        for section in source.sections:
            yield Task(Key('road', section.id), (section.key,),
                       lambda context, section=section: self.build(context, section))


class PythonStreetProgressiveTests(unittest.TestCase):
    def collect(self, steps, pause=None):
        count = 0
        try:
            while True:
                try:
                    self.assertIsNone(next(steps))
                except StopIteration as finished:
                    return finished.value, count
                count += 1
                if pause is not None:
                    pause()
        finally:
            steps.close()

    def test_real_backend_blocking_and_checkpointed_meshes_are_exactly_equal(self):
        source = network()
        blocking = PythonMeshBackend(source.settings).prepare(source)
        stepped, checkpoints = self.collect(PythonMeshBackend(source.settings).prepare_steps(source))
        self.assertGreater(checkpoints, len(blocking.rebuilt))
        for name in ('network_id', 'source_fingerprints', 'task_fingerprints', 'rebuilt',
                     'removed', 'changed', 'dependencies', 'base_revision'):
            self.assertEqual(getattr(blocking, name), getattr(stepped, name))
        self.assertEqual(blocking.timings['overlap_check'], stepped.timings['overlap_check'])
        for key, product in blocking.products.items():
            if isinstance(product, Part):
                actual = stepped.products[key]
                self.assertEqual(product.metadata, actual.metadata)
                np.testing.assert_array_equal(product.mesh.xyz, actual.mesh.xyz)
                np.testing.assert_array_equal(product.mesh.quads, actual.mesh.quads)
                self.assertEqual(product.mesh.material, actual.mesh.material)
                self.assertEqual(set(product.mesh.uv_layers), set(actual.mesh.uv_layers))
                for name, uv in product.mesh.uv_layers.items():
                    np.testing.assert_array_equal(uv, actual.mesh.uv_layers[name])

    def test_cancel_at_every_checkpoint_keeps_last_accepted_products_and_manifest(self):
        source = network()
        engine = PythonMeshBackend(source.settings, providers=(Parts(),))
        engine.commit(engine.prepare(source))
        previous = dict(engine.updates.products)
        manifest = engine.manifest()
        changed = replace(source, sections=tuple(replace(s, width=7.) for s in source.sections))
        _, count = self.collect(engine.prepare_steps(changed))
        for stop in range(count + 1):
            with self.subTest(checkpoint=stop):
                steps = engine.prepare_steps(changed)
                for _ in range(stop):
                    self.assertIsNone(next(steps))
                steps.close()
                self.assertEqual(engine.manifest(), manifest)
                for key, product in previous.items():
                    self.assertIs(engine.updates.products[key], product)
                    self.assertIsNotNone(product.mesh)
        accepted = drain_steps(engine.prepare_steps(changed))
        engine.commit(accepted)
        self.assertEqual(engine.updates.revision, 2)

    def test_task_product_and_digest_become_available_together_in_dependency_order(self):
        events = []
        a, b = Key('first', 0), Key('second', 0)

        def build(context, key):
            self.assertNotIn(key, context.task_digests)
            if key == b:
                self.assertEqual(context[a], a.text())
                self.assertEqual(context.digest(a), context.task_digests[a])
            events.append(key)
            return key.text()

        class Provider:
            def tasks(self, source):
                yield Task(b, (a,), lambda context: build(context, b))
                yield Task(a, (), lambda context: build(context, a))

        engine = UpdateCoordinator((Provider(),), PatternLibrary())
        result = drain_steps(engine.prepare_steps(network()))
        self.assertEqual(events, [a, b])
        self.assertEqual(result.rebuilt, (a, b))
        self.assertEqual(engine.products, {})
        self.assertEqual(engine.task_fingerprints, {})

    def test_intervening_commit_cannot_promote_a_pending_revision(self):
        source = network()
        engine = UpdateCoordinator((Parts(),), PatternLibrary())
        steps = engine.prepare_steps(source)
        self.assertIsNone(next(steps))
        engine.commit(engine.prepare(source))
        pending = drain_steps(steps)
        self.assertEqual(pending.base_revision, 0)
        with self.assertRaisesRegex(ValueError, 'superseded'):
            engine.commit(pending)
        self.assertEqual(engine.revision, 1)

    def test_failed_builder_retains_accepted_state_and_key_diagnostic(self):
        source = network()
        provider = Parts()
        engine = PythonMeshBackend(source.settings, providers=(provider,))
        engine.commit(engine.prepare(source))
        manifest = engine.manifest()
        products = dict(engine.updates.products)

        def fail(context, section):
            raise ValueError('injected failure')

        provider.build = fail
        with self.assertRaisesRegex(ValueError, 'road:1: injected failure'):
            drain_steps(engine.prepare_steps(source, force=True))
        self.assertEqual(engine.manifest(), manifest)
        for key in products:
            self.assertIs(engine.updates.products[key], products[key])

    def test_uv_fallback_and_overlap_phases_are_visible_without_mutating_raw_products(self):
        source = network()
        engine = PythonMeshBackend(source.settings, providers=(Parts(),))
        engine.updates.commit(engine.updates.prepare(source))
        raw = dict(engine.updates.products)
        phases = []
        prepared = drain_steps(engine.prepare_steps(source, on_progress=phases.append))
        self.assertIn('Mapping street UVs for road:1', phases)
        self.assertIn('Validating street overlaps', phases)
        self.assertIn('Street numerical validation complete', phases)
        for key, product in raw.items():
            self.assertNotIn(DEFAULT_UV_MAP, product.mesh.uv_layers)
            self.assertIn(DEFAULT_UV_MAP, prepared.products[key].mesh.uv_layers)
        self.assertEqual(prepared.rebuilt, ())

    def test_nested_timings_exclude_pauses_between_checkpoints(self):
        clock = [0.]

        def build(context, section):
            clock[0] += .25
            return Part(Key('road', section.id), quad(section.id * 3.), {})

        def pause():
            clock[0] += 10.

        source = network()
        engine = PythonMeshBackend(source.settings, providers=(Parts(build),))
        with patch('tested_python_streets.updates.time.perf_counter', side_effect=lambda: clock[0]):
            prepared, count = self.collect(engine.prepare_steps(source), pause)
        self.assertGreater(count, 10)
        self.assertEqual(prepared.timings['total'], .5)
        self.assertEqual(prepared.timings['geometry_seconds'], .5)
        self.assertEqual(sum(prepared.timings['tasks'].values()), .5)

    def test_overlap_steps_preserve_batch_counts_and_inputs_when_closed(self):
        meshes = {Key('road', index): quad(3. * index) for index in range(4)}
        expected = validate_overlaps(meshes)
        phases = []
        actual, count = self.collect(validate_overlaps_steps(meshes, on_progress=phases.append))
        self.assertEqual(actual, expected)
        self.assertIn('Sorting street overlap cells', phases)
        original = {key: value.xyz.copy() for key, value in meshes.items()}
        for stop in range(count + 1):
            steps = validate_overlaps_steps(meshes)
            for _ in range(stop):
                next(steps)
            steps.close()
            for key, value in meshes.items():
                np.testing.assert_array_equal(value.xyz, original[key])
        with self.assertRaisesRegex(ValueError, 'interiors overlap'):
            drain_steps(validate_overlaps_steps({Key('road', 0): quad(), Key('road', 1): quad(.5)}))

    def test_clipped_source_steps_match_blocking_recovery_and_bounded_failure(self):
        for reverse, interpolation, budget in ((False, 0, 40.), (True, 0, 40.),
                                               (False, 3, 40.), (True, 3, 40.),
                                               (False, 3, 1.)):
            with self.subTest(reverse=reverse, interpolation=interpolation, budget=budget):
                source, snapshot, context = clipped_fixture_source(model, reverse=reverse,
                                                                   interpolation=interpolation)
                support = ClippedApproachSupport(snapshot, context, maximum_extension=budget)
                expected = support.prepare(source)
                phases = []
                actual, count = self.collect(support.prepare_steps(source, on_progress=phases.append))
                self.assertEqual(actual, expected)
                self.assertGreater(count, 10)
                self.assertIn('Fitting clipped footprint for intersection:1', phases)
                self.assertIn('Extending original street source section:2', phases)

    def test_cancelled_clipped_source_recovery_preserves_original_input_and_context(self):
        source, snapshot, context = clipped_fixture_source(model, reverse=True)
        before = source.to_dict()
        original_streets, original_context = deepcopy(snapshot.streets), deepcopy(context)
        support = ClippedApproachSupport(snapshot, context)
        expected, count = self.collect(support.prepare_steps(source))
        for stop in range(count + 1):
            steps = support.prepare_steps(source)
            for _ in range(stop):
                next(steps)
            steps.close()
            self.assertEqual(source.to_dict(), before)
            self.assertEqual(snapshot.streets, original_streets)
            self.assertEqual(context, original_context)
        self.assertEqual(drain_steps(support.prepare_steps(source)), expected)


if __name__ == '__main__':
    unittest.main()
