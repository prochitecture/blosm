"""Published source identities survive a changed compound approach inventory."""
from copy import copy
from dataclasses import replace
import unittest
import _python_street_package
from test_python_street_hydrated_external_runs import fixture
from tested_python_streets.model import End, Key
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.ownership import surface_key
from tested_python_streets.separator_ports import (
    PORTS, PHYSICAL_PORT_BINDING_POLICY, VisiblePort, bind_physical_port)
from tested_python_streets.source_mouths import source_mouths
from tested_python_streets.updates import BuildContext


class PhysicalPortBindingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original, _, _, cls.mouth, *_ = fixture()

    def setUp(self):
        self.context = BuildContext(self.original.network, dict(self.original.products), self.original.patterns)
        self.product = copy(self.context[surface_key(self.mouth.owner)])
        self.layout = copy(self.product.group.physical_layout)
        self.product.group = replace(self.product.group, physical_layout=self.layout)
        self.context.products[surface_key(self.mouth.owner)] = self.product
        self.index, self.alias = next((index, alias)
            for index, (arm, station) in enumerate(zip(self.layout.approaches, self.layout.outward_stations))
            for alias in source_mouths(arm, station)
            if alias.original_incidence and self.context[PORTS].resolve(alias.endpoint) is not None)
        self.visible = self.context[PORTS].resolve(self.alias.endpoint)

    def test_reordered_completed_layout_keeps_the_exact_published_source_profile(self):
        original_arms = tuple(self.layout.approaches)
        order = tuple(range(len(original_arms)-1, -1, -1))
        self.layout.approaches = tuple(original_arms[i] for i in order)
        self.layout.outward_stations = tuple(self.layout.outward_stations[i] for i in order)
        binding = bind_physical_port(self.context, self.visible, self.alias.endpoint)
        self.assertEqual(binding.approach_index, order.index(self.index))
        self.assertIs(binding.layout.approaches[binding.approach_index], original_arms[self.index])
        self.assertEqual(binding.alias, self.alias)
        self.assertEqual(binding.owner, self.product.owner)
        self.assertIs(binding.product, self.product)
        self.assertEqual(self.visible.approach_index,
                         self.original[PORTS].resolve(self.alias.endpoint).approach_index)

    def test_old_layout_index_cannot_redirect_an_actual_source_identity(self):
        visible = VisiblePort(self.visible.site, len(self.layout.approaches)+9)
        binding = bind_physical_port(self.context, visible, self.alias.endpoint)
        self.assertEqual(binding.approach_index, self.index)
        self.assertEqual(binding.alias.endpoint, self.alias.endpoint)

    def test_hidden_or_unpublished_source_profile_is_refused(self):
        self.product.peaks = dict(self.product.peaks)
        self.product.peaks.pop((self.alias.endpoint.section_id, self.alias.endpoint.end))
        with self.assertRaisesRegex(ValueError, 'published source profiles'):
            bind_physical_port(self.context, self.visible, self.alias.endpoint)

    def test_changed_published_native_station_is_refused(self):
        self.product.peaks = dict(self.product.peaks)
        key = self.alias.endpoint.section_id, self.alias.endpoint.end
        loop, index, station = self.product.peaks[key]
        self.product.peaks[key] = loop, index, station+.01
        with self.assertRaisesRegex(ValueError, 'source station'):
            bind_physical_port(self.context, self.visible, self.alias.endpoint)

    def test_duplicate_physical_alias_is_refused(self):
        self.layout.approaches = (*self.layout.approaches, self.layout.approaches[self.index])
        self.layout.outward_stations = (*self.layout.outward_stations, self.layout.outward_stations[self.index])
        with self.assertRaisesRegex(ValueError, 'unique published source identity'):
            bind_physical_port(self.context, self.visible, self.alias.endpoint)

    def test_missing_physical_alias_is_refused(self):
        self.layout.approaches = tuple(a for i, a in enumerate(self.layout.approaches) if i != self.index)
        self.layout.outward_stations = tuple(s for i, s in enumerate(self.layout.outward_stations) if i != self.index)
        with self.assertRaisesRegex(ValueError, 'unique published source identity'):
            bind_physical_port(self.context, self.visible, self.alias.endpoint)

    def test_derived_port_retains_its_separate_owner_scoped_identity(self):
        registry = self.context[PORTS]
        owner, sid, end = next(key for key in registry.derived_ports if key[0] == self.mouth.owner)
        visible = registry.resolve_derived(owner, sid, end)
        binding = bind_physical_port(self.context, visible, End(sid, end))
        self.assertFalse(binding.alias.original_incidence)
        self.assertEqual(binding.alias, visible.alias)
        self.assertEqual(binding.owner, owner)
        with self.assertRaises(ValueError):
            bind_physical_port(self.context, self.visible, End(sid, end))

    def test_port_task_records_the_physical_binding_policy(self):
        backend = PythonMeshBackend(self.context.network.settings)
        task, = (task for provider in backend.providers for task in provider.tasks(self.context.network)
                 if task.key == PORTS)
        self.assertEqual(task.version, PHYSICAL_PORT_BINDING_POLICY)
        plans = [task for provider in backend.providers for task in provider.tasks(self.context.network)
                 if task.key.kind == 'bundle_plan']
        self.assertTrue(plans)
        self.assertTrue(all(task.version == PHYSICAL_PORT_BINDING_POLICY for task in plans))


if __name__ == '__main__':
    unittest.main()
