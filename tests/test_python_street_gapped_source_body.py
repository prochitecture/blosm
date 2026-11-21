"""Actual separated junction profiles bound opposite-body-only pavement tiles."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_merged_run_inventory import side_junction, prepared
from tested_python_streets.gapped_source_inventory import GappedSourceInventory, GappedSourceUnavailable
from tested_python_streets.gapped_source_body import GappedSourceBodyPiece
from tested_python_streets.gapped_source_plan import GappedSourcePlan
from tested_python_streets.merged_run_inventory import inventory_frame
from tested_python_streets.separator_ports import PORTS
from tested_python_streets.model import Key
from tested_python_streets.uv import apply_mapping, PlanarMapping


class GappedSourceBodyTests(unittest.TestCase):
    def test_private_bank_terminal_sampling_preserves_tips_and_respects_chord_error(self):
        bank = np.array([(0., 0., .08), (.004, 2e-7, .08), (.704, 0., .08),
                         (1.404, 0., .08), (2.104, 0., .08)])
        for source in (bank, bank[::-1]):
            before = source.copy()
            result = GappedSourceBodyPiece._balance_bank_ends(source, 1e-6)
            np.testing.assert_array_equal(source, before)
            np.testing.assert_array_equal(result[[0,-1]], source[[0,-1]])
            np.testing.assert_array_equal(result[:,2], source[:,2])
            self.assertGreater(np.min(np.linalg.norm(np.diff(result, axis=0), axis=1)), .35)
        bank[1,1] = .001
        self.assertIsNone(GappedSourceBodyPiece._balance_bank_ends(bank, 1e-6))

    @classmethod
    def setUpClass(cls):
        cls.fixtures = []
        for reverse in (False, True):
            context, approach = prepared(side_junction(reverse))
            bundle = context.network.bundles[0]
            inventory = GappedSourceInventory.build(context, bundle, run=approach.run,
                frame=inventory_frame(context, bundle, approach))
            cls.fixtures.append((context, inventory))
        source = side_junction(True)
        angle = .713
        matrix = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        sections = []
        for section in source.sections:
            xyz = np.asarray(section.points).copy()
            xyz[:, :2] = xyz[:, :2]@matrix.T+np.array((1200., -1800.))
            sections.append(replace(section, points=tuple(map(tuple, xyz))))
        context, approach = prepared(replace(source, sections=tuple(sections)))
        bundle = context.network.bundles[0]
        cls.fixtures.append((context, GappedSourceInventory.build(context, bundle, run=approach.run,
            frame=inventory_frame(context, bundle, approach))))

    def test_real_side_junction_keeps_all_profiles_whole_curbs_and_opposite_source_ownership(self):
        for context, inventory in self.fixtures:
            with self.subTest(reversed=not inventory.run.members[0].visits[0].forward):
                before = context.network.to_dict()
                registry = context[PORTS]
                original_ports, derived_ports = dict(registry.ports), dict(registry.derived_ports)
                paths = {sid: context[Key('path', sid)] for sid in inventory.source_sections}
                piece = GappedSourceBodyPiece()
                pavement, grass, report = piece.build(context, inventory, namespace=('test gapped', id(context)))
                self.assertFalse(report['junction_area_remeshed'])
                self.assertFalse(report['runtime_discovery_enabled'])
                self.assertEqual(report['new_source_aliases'], 0)
                self.assertEqual(report['exact_source_profiles'], 6)
                self.assertTrue(report['boundary_manifold'])
                self.assertGreater(report['actual_owner_curb_edges'], 0)
                self.assertGreater(report['shared_whole_grass_edges'], 0)
                self.assertEqual(report['complete_source_width']['maximum_error'], .002)
                self.assertEqual(set(report['complete_source_width']['complete_original_source_ids']), {18, 17, 29})
                self.assertLess(report['total_quads'], 50000)
                for role, spans in enumerate(inventory.spans_by_role):
                    for index, span in enumerate(spans):
                        for end, node in enumerate((span.head, span.tail)):
                            actual = node.physical_xyz if node.prepared_mouth is None else node.prepared_mouth.xyz
                            np.testing.assert_array_equal(pavement.xyz[piece.profile_ids[role, index, end]], actual)
                for mesh, row in zip((pavement, grass), report['stored_checks']):
                    self.assertLessEqual(row['stored_quality']['maximum_edge_ratio'], 8.)
                    if mesh.material != 'grass':
                        self.assertGreaterEqual(row['stored_quality']['minimum_corner_sine'], .03)
                        self.assertLessEqual(row['longitudinal_bound'], .03+1e-8)
                    self.assertLessEqual(row['maximum_slope'], .1+1e-6)
                    self.assertLessEqual(row['drainage']['maximum_fill'], 1e-6)
                    self.assertTrue(np.isfinite(apply_mapping(mesh, PlanarMapping()).uv_layers['UVMap']).all())
                self.assertEqual(registry.ports, original_ports)
                self.assertEqual(registry.derived_ports, derived_ports)
                self.assertEqual(context.network.to_dict(), before)
                for sid, path in paths.items(): self.assertIs(context[Key('path', sid)], path)

    def test_changed_source_interval_and_actual_junction_curb_reject_before_meshing(self):
        context, inventory = self.fixtures[0]
        span = inventory.spans_by_role[0][0]
        rows = list(span.source_intervals)
        rows[0] = dict(rows[0], last=rows[0]['last']-.01)
        wrong = replace(span, source_intervals=tuple(rows))
        own = ((wrong, *inventory.spans_by_role[0][1:]), inventory.spans_by_role[1])
        with self.assertRaisesRegex(ValueError, 'actual complete profile'):
            GappedSourceBodyPiece().build(context, replace(inventory, spans_by_role=own), namespace='changed')
        gap = inventory.gaps[0]
        points = gap.bridge.physical_xyz.copy(); points[0, 0] += .001
        changed = replace(gap, bridge=replace(gap.bridge, physical_xyz=points))
        with self.assertRaisesRegex(ValueError, 'real near-side junction curb'):
            GappedSourceBodyPiece().build(context, replace(inventory, gaps=(changed,)), namespace='changed')

    def test_incomplete_original_member_stays_ineligible_and_no_registry_is_extended(self):
        context, inventory = self.fixtures[0]
        bundle = inventory.bundle
        incomplete = replace(bundle, members=(replace(bundle.members[0], complete=False), bundle.members[1]))
        with self.assertRaises(GappedSourceUnavailable):
            GappedSourcePlan.build(context, incomplete)


if __name__ == '__main__': unittest.main()
