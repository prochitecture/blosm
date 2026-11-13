"""Early numerical product notifications remain optional and transactional."""
from dataclasses import FrozenInstanceError, replace
import unittest
from unittest.mock import patch

import numpy as np
import _python_street_package
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.mesh import Mesh
from tested_python_streets.model import Key, Network, Section
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.providers import Part
from tested_python_streets.updates import (
    ProductReady, Task, UpdateCoordinator, _PreparationProgress, drain_steps)
from tested_python_streets.uv import DEFAULT_UV_MAP


def network():
    return Network('early-product-events', (
        Section(2, 2, ((0., 20., 0.), (8., 20., 0.))),
        Section(1, 1, ((0., 0., 0.), (8., 0., 0.)))))


def part(section, *, overlap=False):
    x = 0. if overlap else section.id * 3.
    mesh = Mesh(np.array(((x, 0., 0.), (x + 1., 0., 0.),
                          (x + 1., 1., 0.), (x, 1., 0.))), np.array(((0, 1, 2, 3),)))
    return Part(Key('road', section.id), mesh, {'source': section.id})


class Parts:
    source_kinds = frozenset(('section',))

    def __init__(self, build=None):
        self.build = build or (lambda context, section: part(section))

    def tasks(self, source):
        for section in source.sections:
            yield Task(Key('road', section.id), (section.key,),
                       lambda context, section=section: self.build(context, section))


def collect(steps, inspect=None):
    commands = []
    try:
        while True:
            try:
                command = next(steps)
            except StopIteration as finished:
                return finished.value, commands
            commands.append(command)
            if inspect is not None:
                inspect(command)
    finally:
        steps.close()


class PythonStreetProductEventTests(unittest.TestCase):
    def test_default_callers_keep_identical_none_checkpoints_and_progress(self):
        source = network()
        for make_engine in (
            lambda: UpdateCoordinator((Parts(),), PatternLibrary()),
            lambda: PythonMeshBackend(source.settings, providers=(Parts(),)),
        ):
            with self.subTest(engine=make_engine().__class__.__name__):
                default_phases, explicit_phases = [], []
                default, commands = collect(make_engine().prepare_steps(
                    source, on_progress=default_phases.append))
                explicit, explicit_commands = collect(make_engine().prepare_steps(
                    source, on_progress=explicit_phases.append, emit_products=False))
                self.assertTrue(all(command is None for command in commands))
                self.assertEqual(commands, explicit_commands)
                self.assertEqual(default_phases, explicit_phases)
                self.assertEqual(default.task_fingerprints, explicit.task_fingerprints)
                _, opted_in = collect(make_engine().prepare_steps(source, emit_products=True))
                self.assertEqual(commands, [item for item in opted_in if item is None])

    def test_event_arrives_after_product_and_digest_before_dependent_builder(self):
        first, second = Key('first', 0), Key('second', 0)
        built, contexts, events = [], {}, []

        def build(context, key):
            if key == second:
                self.assertEqual(events, [first])
                self.assertEqual(context[first], first.text())
                self.assertEqual(context.digest(first), context.task_digests[first])
            contexts[key] = context
            built.append(key)
            return key.text()

        class Provider:
            def tasks(self, source):
                yield Task(second, (first,), lambda context: build(context, second))
                yield Task(first, (), lambda context: build(context, first))

        def inspect(command):
            if isinstance(command, ProductReady):
                events.append(command.key)
                self.assertEqual(built, events)
                self.assertIs(command.product, contexts[command.key][command.key])
                self.assertEqual(contexts[command.key].digest(command.key),
                                 contexts[command.key].task_digests[command.key])
                with self.assertRaises(FrozenInstanceError):
                    command.key = second

        engine = UpdateCoordinator((Provider(),), PatternLibrary())
        prepared, _ = collect(engine.prepare_steps(network(), emit_products=True), inspect)
        self.assertEqual(tuple(events), prepared.rebuilt)
        self.assertEqual(engine.products, {})
        self.assertEqual(engine.task_fingerprints, {})

    def test_only_rebuilt_products_emit_in_dependency_order(self):
        source = network()
        engine = UpdateCoordinator((Parts(),), PatternLibrary())
        engine.commit(engine.prepare(source))
        unchanged, commands = collect(engine.prepare_steps(source, emit_products=True))
        self.assertEqual(unchanged.rebuilt, ())
        self.assertTrue(all(command is None for command in commands))
        changed = replace(source, sections=(source.sections[0], replace(source.sections[1], width=8.)))
        prepared, commands = collect(engine.prepare_steps(changed, emit_products=True))
        self.assertEqual([command.key for command in commands if isinstance(command, ProductReady)],
                         list(prepared.rebuilt))
        self.assertEqual(prepared.rebuilt, (Key('road', 1),))
        self.assertIs(prepared.products[Key('road', 2)], engine.products[Key('road', 2)])
        forced, commands = collect(engine.prepare_steps(source, force=True, emit_products=True))
        self.assertEqual([command.key for command in commands if isinstance(command, ProductReady)],
                         list(forced.rebuilt))
        self.assertEqual(forced.rebuilt, (Key('road', 1), Key('road', 2)))

    def test_backend_events_precede_uv_mapping_and_whole_network_validation(self):
        source, phases, raw = network(), [], {}
        engine = PythonMeshBackend(source.settings, providers=(Parts(),))

        def inspect(command):
            if isinstance(command, ProductReady):
                self.assertNotIn('Street geometry prepared', phases)
                self.assertNotIn('Checking street UV maps', phases)
                self.assertNotIn('Validating street overlaps', phases)
                self.assertNotIn(DEFAULT_UV_MAP, command.product.mesh.uv_layers)
                raw[command.key] = command.product

        prepared, _ = collect(engine.prepare_steps(
            source, emit_products=True, on_progress=phases.append), inspect)
        self.assertEqual(tuple(raw), prepared.rebuilt)
        self.assertIn('Street numerical validation complete', phases)
        for key, original in raw.items():
            self.assertNotIn(DEFAULT_UV_MAP, original.mesh.uv_layers)
            self.assertIn(DEFAULT_UV_MAP, prepared.products[key].mesh.uv_layers)
            self.assertIsNot(original, prepared.products[key])
        self.assertEqual(engine.updates.products, {})

    def test_close_at_every_yield_preserves_accepted_products_and_manifest(self):
        source = network()
        engine = PythonMeshBackend(source.settings, providers=(Parts(),))
        engine.commit(engine.prepare(source))
        accepted, manifest = dict(engine.updates.products), engine.manifest()
        changed = replace(source, sections=tuple(replace(section, width=8.) for section in source.sections))
        _, commands = collect(engine.prepare_steps(changed, emit_products=True))
        for stop in range(len(commands) + 1):
            with self.subTest(stop=stop):
                steps = engine.prepare_steps(changed, emit_products=True)
                for _ in range(stop):
                    next(steps)
                steps.close()
                self.assertEqual(engine.manifest(), manifest)
                for key, product in accepted.items():
                    self.assertIs(engine.updates.products[key], product)
                    self.assertIsNotNone(product.mesh)

    def test_late_overlap_failure_follows_events_without_accepting_products(self):
        source, events = network(), []
        provider = Parts()
        engine = PythonMeshBackend(source.settings, providers=(provider,))
        engine.commit(engine.prepare(source))
        accepted, manifest = dict(engine.updates.products), engine.manifest()
        provider.build = lambda context, section: part(section, overlap=True)
        with self.assertRaisesRegex(ValueError, 'interiors overlap'):
            collect(engine.prepare_steps(source, force=True, emit_products=True),
                    lambda command: events.append(command) if isinstance(command, ProductReady) else None)
        self.assertEqual([event.key for event in events], [Key('road', 1), Key('road', 2)])
        self.assertEqual(engine.manifest(), manifest)
        for key, product in accepted.items():
            self.assertIs(engine.updates.products[key], product)

    def test_forward_preserves_commands_return_value_and_close(self):
        event, closed, clock = ProductReady(Key('road', 1), object()), [], [0.]

        def nested():
            try:
                clock[0] += .5
                yield event
                clock[0] += .25
                yield None
                return 'finished'
            finally:
                closed.append(True)

        with patch('tested_python_streets.updates.time.perf_counter', side_effect=lambda: clock[0]):
            progress = _PreparationProgress()
            steps = progress.forward(nested())
            self.assertIs(next(steps), event)
            clock[0] += 10.
            self.assertIsNone(next(steps))
            clock[0] += 20.
            self.assertEqual(drain_steps(steps), 'finished')
            self.assertEqual(progress.elapsed(), .75)
            self.assertEqual(closed, [True])
            steps = progress.forward(nested())
            self.assertIs(next(steps), event)
            clock[0] += 30.
            steps.close()
            self.assertEqual(progress.elapsed(), 1.25)
            self.assertEqual(closed, [True, True])

    def test_nested_backend_timings_exclude_pauses_on_product_events(self):
        clock = [0.]

        def build(context, section):
            clock[0] += .25
            return part(section)

        def pause(command):
            clock[0] += 20. if isinstance(command, ProductReady) else 10.

        source = network()
        engine = PythonMeshBackend(source.settings, providers=(Parts(build),))
        with patch('tested_python_streets.updates.time.perf_counter', side_effect=lambda: clock[0]):
            prepared, commands = collect(engine.prepare_steps(source, emit_products=True), pause)
        self.assertEqual(sum(isinstance(item, ProductReady) for item in commands), 2)
        self.assertEqual(prepared.timings['total'], .5)
        self.assertEqual(prepared.timings['geometry_seconds'], .5)
        self.assertEqual(sum(prepared.timings['tasks'].values()), .5)

    def test_real_backend_final_geometry_uvs_and_manifest_match_blocking_exactly(self):
        source = network()
        blocking_engine, event_engine = PythonMeshBackend(source.settings), PythonMeshBackend(source.settings)
        blocking = blocking_engine.prepare(source)
        prepared, commands = collect(event_engine.prepare_steps(source, emit_products=True))
        self.assertEqual([event.key for event in commands if isinstance(event, ProductReady)],
                         list(prepared.rebuilt))
        for name in ('network_id', 'source_fingerprints', 'task_fingerprints', 'rebuilt',
                     'removed', 'changed', 'dependencies', 'base_revision'):
            self.assertEqual(getattr(blocking, name), getattr(prepared, name))
        self.assertEqual(blocking.timings['overlap_check'], prepared.timings['overlap_check'])
        for key, product in blocking.products.items():
            if isinstance(product, Part):
                actual = prepared.products[key]
                self.assertEqual(product.metadata, actual.metadata)
                np.testing.assert_array_equal(product.mesh.xyz, actual.mesh.xyz)
                np.testing.assert_array_equal(product.mesh.quads, actual.mesh.quads)
                self.assertEqual(product.mesh.material, actual.mesh.material)
                self.assertEqual(set(product.mesh.uv_layers), set(actual.mesh.uv_layers))
                for name, uv in product.mesh.uv_layers.items():
                    np.testing.assert_array_equal(uv, actual.mesh.uv_layers[name])
        blocking_engine.commit(blocking)
        event_engine.commit(prepared)
        self.assertEqual(blocking_engine.manifest(), event_engine.manifest())


if __name__ == '__main__':
    unittest.main()
