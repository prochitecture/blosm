"""Resolved source-cycle roots prepare real side owners before surfaces."""
from dataclasses import replace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from test_python_street_merged_run_inventory import side_junction, prepared
from tested_python_streets.bundle_merged_entrances import PairedAngularEntrance, preparation_key
from tested_python_streets.bundles import BundleProvider
from tested_python_streets.model import BundleMember, Key
from tested_python_streets.ownership import PLAN, JunctionOwnership, bundle_signature
from tested_python_streets.separator_ports import PORTS
from tested_python_streets.source_companions import PreparedSourceRun, prepare


class SourceCompanionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context, _ = prepared(side_junction(segment=6))
        context = cls.context
        original = context.network.bundles[0]
        cls.specification = context[preparation_key(original)]
        cls.approach = next(a for a in cls.specification.layout.approaches
            if isinstance(a, PairedAngularEntrance)
            and any(v.section_id == 18 for m in a.run.members for v in m.visits))
        cls.bundle = replace(original, members=tuple(BundleMember(m.visits) for m in cls.approach.run.members),
                             interior_section_ids=())
        context.network = replace(context.network, bundles=(cls.bundle,))
        cls.before = tuple((a, tuple(s.copy() for s in a.source_curve_stations))
            for a in cls.specification.layout.approaches if isinstance(a, PairedAngularEntrance))
        cls.groups, cls.resources = prepare(context, (cls.specification.group,))

    def test_side_incidence_has_its_own_original_owner_and_complete_run_resource(self):
        resource = self.resources[self.bundle.id]
        self.assertIsInstance(resource, PreparedSourceRun)
        self.assertEqual(resource.expected_run_count, 1)
        self.assertIs(resource.group, self.specification.group)
        self.assertIs(resource.layout.approaches[0], self.approach)
        self.assertEqual([side.owner.id for side in resource.paired_sources], [600])
        self.assertEqual(len(self.groups), 2)
        self.assertTrue(resource.source_phase_proof['overlapping_station_arrays_exact'])
        for approach, arrays in self.before:
            for actual, old in zip(approach.source_curve_stations, arrays):
                np.testing.assert_array_equal(actual, old)

    def test_opposite_source_invalidates_side_surface_without_a_fake_incidence(self):
        resource = self.resources[self.bundle.id]
        side = resource.paired_sources[0]
        original_ends = {end.section_id for site in self.context.network.intersections
                         if site.key == side.owner for end in site.ends}
        self.assertNotIn(29, original_ends)
        self.assertIn(29, side.group.physical_layout.paired_source_dependencies)
        context = self.context
        saved = dict(context.task_digests)
        try:
            context.task_digests = {**saved, Key('path', 29): 'source before'}
            before = bundle_signature(context, side.group)
            context.task_digests[Key('path', 29)] = 'source edited'
            self.assertNotEqual(before, bundle_signature(context, side.group))
        finally:
            context.task_digests = saved

    def test_declared_inputs_include_opposite_profiles_and_actual_remainder_attachments(self):
        network = self.context.network
        ownership = next(task for task in JunctionOwnership.tasks(network) if task.key == PLAN)
        for sid in self.bundle.sections:
            self.assertIn(Key('path', sid), ownership.dependencies)
            self.assertIn(Key('profile', sid), ownership.dependencies)
        plan = next(task for task in BundleProvider().tasks(network)
                    if task.key == Key('bundle_plan', self.bundle.id))
        self.assertIn(PORTS, plan.dependencies)
        self.assertIn(Key('path', 401), plan.dependencies)
        self.assertIn(Key('profile', 401), plan.dependencies)

    def test_existing_prepared_resource_is_not_repeated_and_owners_cannot_overlap(self):
        with patch('tested_python_streets.source_companions.canonical_reprepare') as rebuild:
            groups, resources = prepare(self.context, (self.specification.group,),
                                         prepared_bundle_ids=(self.bundle.id,))
            rebuild.assert_not_called()
        self.assertEqual(groups, (self.specification.group,))
        self.assertEqual(resources, {})
        with self.assertRaisesRegex(ValueError, 'existing original physical owner'):
            prepare(self.context, self.groups)


if __name__ == '__main__':
    unittest.main()
