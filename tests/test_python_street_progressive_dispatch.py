"""Legacy subclass overrides survive optional cooperative street dispatch."""
import importlib
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import unittest


ROOT = Path(__file__).resolve().parents[1]
PREFIX = '_street_progressive_dispatch_tests'


def package(suffix, path):
    module = ModuleType(PREFIX + suffix)
    module.__path__ = [str(ROOT / path)]
    sys.modules[module.__name__] = module
    return module


def stub(suffix, **attributes):
    module = ModuleType(PREFIX + suffix)
    module.__dict__.update(attributes)
    sys.modules[module.__name__] = module


# Import the real manager/facade while replacing unrelated geometry/style
# dependencies. Scheduling itself must remain usable without Blender.
package('', '')
way = package('.way', 'way')
way.Way = way.Railway = object
package('.defs', 'defs')
package('.way.waymap', 'way/waymap')
package('.way.item', 'way/item')
package('.way.street_network', 'way/street_network')
stub('.style', StyleStore=object)
stub('.way.waymap.waymap', WayMap=object)
stub('.way.item.street', Street=object)
stub('.way.item.connectors', IntConnector=object)
WayManager = importlib.import_module(PREFIX + '.way.manager').WayManager
StreetRenderer = importlib.import_module(PREFIX + '.way.renderer_streets').StreetRenderer


class Session:
    def __init__(self):
        self.calls = self.units = 0

    def call(self, callback):
        self.calls += 1
        return callback(), []

    def unit(self, callback):
        self.units += 1
        callback()
        yield None


class BaseAction:
    def __init__(self, events):
        self.events = events

    def do(self, manager):
        self.events.append('base action')

    def do_steps(self, manager, session=None):
        self.events.append('base action producer')
        yield None


class LegacyAction(BaseAction):
    def do(self, manager):
        self.events.append('custom action')


class CooperativeAction(LegacyAction):
    def do_steps(self, manager, session=None):
        self.events.append('custom action producer')
        yield None


class BaseRenderer:
    def __init__(self, events):
        self.events = events

    def render(self, manager, data):
        self.events.append('base renderer')

    def render_steps(self, manager, data, session):
        self.events.append('base renderer producer')
        yield None


class LegacyRenderer(BaseRenderer):
    def render(self, manager, data):
        self.events.append('custom renderer')


class CooperativeRenderer(LegacyRenderer):
    def render_steps(self, manager, data, session):
        self.events.append('custom renderer producer')
        yield None


class ProgressiveStreetDispatchTests(unittest.TestCase):
    @staticmethod
    def manager():
        manager = WayManager.__new__(WayManager)
        manager.app = SimpleNamespace(trafficSide='rht')
        manager.data = object()
        manager.getAllWays = lambda: iter(())
        manager.actions, manager.renderers = [], []
        return manager

    def test_legacy_action_override_runs_once_in_both_drivers(self):
        for progressive in (False, True):
            events, session, manager = [], Session(), self.manager()
            manager.actions = [LegacyAction(events)]
            if progressive:
                list(manager.process_steps(session))
            else:
                manager.process()
            self.assertEqual(events, ['custom action'])
            self.assertEqual(session.calls, int(progressive))

    def test_action_can_opt_back_into_its_own_producer(self):
        for action, expected in ((BaseAction, 'base action producer'),
                                 (CooperativeAction, 'custom action producer')):
            events, session, manager = [], Session(), self.manager()
            manager.actions = [action(events)]
            list(manager.process_steps(session))
            self.assertEqual(events, [expected])
            self.assertEqual(session.calls, 0)

    def test_manager_preserves_sync_renderer_override_as_one_unit(self):
        events, session, manager = [], Session(), self.manager()
        manager.renderers = [LegacyRenderer(events)]
        list(manager.render_steps(session))
        self.assertEqual(events, ['custom renderer'])
        self.assertEqual(session.units, 1)

    def test_facade_preserves_sync_backend_override_as_one_unit(self):
        events, session = [], Session()
        renderer = StreetRenderer.__new__(StreetRenderer)
        renderer.implementation = LegacyRenderer(events)
        list(renderer.render_steps('manager', 'data', session))
        self.assertEqual(events, ['custom renderer'])
        self.assertEqual(session.units, 1)

    def test_explicit_backend_producer_remains_cooperative_through_both_layers(self):
        events, session, manager = [], Session(), self.manager()
        renderer = StreetRenderer.__new__(StreetRenderer)
        renderer.implementation = CooperativeRenderer(events)
        manager.renderers = [renderer]
        list(manager.render_steps(session))
        self.assertEqual(events, ['custom renderer producer'])
        self.assertEqual(session.units, 0)


if __name__ == '__main__':
    unittest.main()
