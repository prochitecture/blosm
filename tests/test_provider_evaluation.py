"""Native consumers must stay disabled until their providers are ready."""
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch


class Modifier:
    type = 'NODES'
    def __init__(self, schema, enabled=True):
        self.node_group = {'blosm_street_input_schema': schema}
        self.show_viewport = enabled


class Object:
    def __init__(self, name, schema, parent=None, enabled=True):
        self.name = name; self.parent = parent
        self.modifiers = [Modifier(schema, enabled)]
        self.data = SimpleNamespace(vertices=(), polygons=())
    def evaluated_get(self, graph): return self


class ProviderEvaluationTests(unittest.TestCase):
    def setUp(self):
        self.controller = Object('controller', 'original_spans.1')
        schemas = ('original_patch_provider.1', 'original_strip_provider.1',
            'free_bundle_provider.1', 'separator_ownership_provider.1',
            'separator_fill_provider.1', 'original_intake_provider.1', 'aggregate_provider.1')
        self.providers = [Object(str(i), s, self.controller) for i, s in enumerate(schemas)]
        self.display = Object('labels', 'original_failure_display.1', self.controller)
        self.disabled = Object('disabled', 'original_strip_provider.1', self.controller, False)
        self.unrelated = Object('unrelated', 'original_intake_provider.1')
        objects = [self.controller, self.display, self.disabled, self.unrelated, *reversed(self.providers)]
        self.updates = []
        fake = SimpleNamespace(data=SimpleNamespace(objects=objects), context=SimpleNamespace(
            view_layer=SimpleNamespace(update=lambda: self.updates.append(self.states())),
            evaluated_depsgraph_get=lambda: None))
        path = Path(__file__).resolve().parents[1]/'way/street_network/provider_evaluation.py'
        spec = importlib.util.spec_from_file_location('_provider_evaluation_test', path)
        self.module = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, {'bpy': fake}): spec.loader.exec_module(self.module)
    def states(self):
        return tuple(o.modifiers[0].show_viewport for o in
                     [self.controller, self.display, self.disabled, self.unrelated, *self.providers])

    def test_complete_native_order_keeps_consumers_and_labels_deferred(self):
        before = self.states()
        rows = self.module.evaluate_prepared_network(self.controller)
        self.assertEqual([r['name'] for r in rows], [o.name for o in self.providers])
        for i, state in enumerate(self.updates[:len(self.providers)]):
            self.assertEqual(state[:4], (False, False, False, True))
            self.assertEqual(state[4:], (True,)*(i+1)+(False,)*(len(self.providers)-i-1))
        self.assertEqual(self.updates[-2][:2], (True, False))
        self.assertEqual(self.updates[-1], before)
        self.assertEqual(self.states(), before)

    def test_failure_restores_initial_flags_including_disabled_provider(self):
        before = self.states()
        def fail(event):
            if event['phase'] == 'starting_provider' and event['name'] == '2':
                raise RuntimeError('test failure')
        with self.assertRaisesRegex(RuntimeError, 'test failure'):
            self.module.evaluate_prepared_network(self.controller, observe=fail)
        self.assertEqual(self.states(), before)

    def test_disabled_controller_is_not_enabled_by_schedule(self):
        self.controller.modifiers[0].show_viewport = False
        self.display.modifiers[0].show_viewport = False
        before = self.states()
        self.module.evaluate_prepared_network(self.controller)
        self.assertEqual(self.states(), before)
        self.assertTrue(all(not state[0] and not state[1] for state in self.updates))


if __name__ == '__main__': unittest.main()
