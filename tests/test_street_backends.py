"""Host checks of lazy street backends and the importer lifecycle facade."""
import importlib
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[1]
PREFIX = '_street_backends_tests'
for suffix, directory in (('', 'way'), ('.street_network', 'way/street_network')):
    package = ModuleType(PREFIX+suffix)
    package.__path__ = [str(ROOT/directory)]
    sys.modules[package.__name__] = package
backends = importlib.import_module(PREFIX+'.street_network.backends')
facade = importlib.import_module(PREFIX+'.renderer_streets')


class RecordingBackend(backends.StreetRenderBackend):
    def __init__(self, app, **options):
        self.app, self.options, self.calls = app, options, []

    def prepare(self):
        self.calls.append('prepare')
        return 'prepared'

    def render(self, manager, data):
        self.calls.append(('render', manager, data))
        return 'rendered'

    def finalize(self):
        self.calls.append('finalize')
        return 'finalized'

    def cleanup(self):
        self.calls.append('cleanup')
        return 'cleaned'


def definition(key='python', number=0, factory=RecordingBackend, **kwargs):
    return backends.StreetBackendDefinition(key, key.title(), 'Test provider', number, factory, **kwargs)


class StreetBackendTests(unittest.TestCase):
    def test_backend_lookup_and_menu_do_not_construct_an_engine(self):
        calls = []
        def factory(app, **options):
            calls.append((app, options))
            return RecordingBackend(app, **options)
        registry = backends.StreetBackendRegistry()
        registry.register(definition(factory=factory))
        items = registry.enum_items()
        self.assertIs(items, registry.enum_items())
        self.assertEqual(items, (('python', 'Python', 'Test provider', 0),))
        registry.get('python')
        self.assertEqual(calls, [])
        app = object()
        result = registry.create('python', app)
        self.assertIs(result.app, app)
        self.assertEqual(len(calls), 1)

    def test_duplicate_keys_values_and_retired_values_are_rejected_atomically(self):
        registry = backends.StreetBackendRegistry()
        registry.register(definition())
        for candidate in (definition(), definition('other', 0), definition('old_gn', 1), definition('old_legacy', 2)):
            with self.subTest(candidate=candidate.key), self.assertRaises(ValueError):
                registry.register(candidate)
        self.assertEqual(len(registry), 1)

    def test_future_gn_identity_uses_a_new_stable_number_and_can_be_removed(self):
        registry = backends.StreetBackendRegistry()
        registry.register(definition())
        registry.register(definition('geometry_nodes', 3))
        self.assertEqual(registry.enum_items()[-1][::3], ('geometry_nodes', 3))
        registry.unregister('geometry_nodes')
        with self.assertRaises(ValueError):
            registry.get('geometry_nodes')
        with self.assertRaises(ValueError):
            registry.unregister('python')

    def test_facade_forwards_options_results_and_lifecycle_to_another_backend(self):
        backends.registry.register(definition('alternate', 3))
        try:
            app = SimpleNamespace(streetBackend='alternate')
            design, loader, features, manager, data = (object() for _ in range(5))
            renderer = facade.StreetRenderer(app, design=design, features=features, asset_loader=loader)
            self.assertIs(renderer.implementation.app, app)
            self.assertEqual(renderer.implementation.options,
                             dict(design=design, features=features, asset_loader=loader))
            self.assertEqual(renderer.prepare(), 'prepared')
            self.assertEqual(renderer.render(manager, data), 'rendered')
            self.assertEqual(renderer.finalize(), 'finalized')
            self.assertEqual(renderer.cleanup(), 'cleaned')
            self.assertEqual(renderer.implementation.calls,
                             ['prepare', ('render', manager, data), 'finalize', 'cleanup'])
        finally:
            backends.registry.unregister('alternate')

    def test_unknown_backend_and_incomplete_implementation_fail_explicitly(self):
        with self.assertRaisesRegex(ValueError, 'Unknown street backend'):
            facade.StreetRenderer(SimpleNamespace(), backend='missing')
        registry = backends.StreetBackendRegistry()
        registry.register(definition(factory=lambda app, **kwargs: object()))
        with self.assertRaisesRegex(TypeError, 'StreetRenderBackend'):
            registry.create('python', object())

    def test_settings_hook_does_not_evaluate_the_backend(self):
        draws = []
        registry = backends.StreetBackendRegistry()
        registry.register(definition(draw_settings=lambda layout, settings: draws.append((layout, settings))))
        layout, settings = object(), object()
        registry.get('python').draw_settings(layout, settings)
        self.assertEqual(draws, [(layout, settings)])

    def test_registry_and_facade_have_not_imported_a_geometry_implementation(self):
        self.assertFalse(any(name.startswith(PREFIX+'.street_network.python_backend') for name in sys.modules))


if __name__ == '__main__':
    unittest.main()
