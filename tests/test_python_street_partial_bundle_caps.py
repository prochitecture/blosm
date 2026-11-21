"""Finite partial Bundle caps follow actual junction curbs to free/scope ends."""
from dataclasses import replace
import unittest
from unittest.mock import patch
import _python_street_package
from tested_python_streets.model import Section, Network, Bundle, BundleMember, Traversal, Intersection, End, Key, Settings
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_partial_caps import terminal_chain, dependencies


def bergamo_ground_ramps():
    # Original controls from the airport Bergamo extract, with the removed
    # bridge span represented by its original scope incidences.
    sections = (
        Section(93, 93, ((-1696.84228515625, 773.0778198242188, 0.0), (-1698.4686279296875, 770.3397827148438, 0.0), (-1707.9224853515625, 756.8281860351562, 0.0)), width=14.0, rise=0.08, scope_ends=()),
        Section(95, 95, ((-1707.9224853515625, 756.8281860351562, 0.0), (-1711.3333740234375, 741.6785278320312, 0.0), (-1710.27099609375, 729.8672485351562, 0.0), (-1707.3258056640625, 720.3375244140625, 0.0), (-1701.363037109375, 709.704833984375, 0.0), (-1694.7850341796875, 701.8995361328125, 0.0), (-1687.51416015625, 696.0198974609375, 0.0), (-1682.023681640625, 692.9794311523438, 0.0)), width=14.0, rise=0.08, scope_ends=()),
        Section(96, 96, ((-1707.9224853515625, 756.8281860351562, 0.0), (-1719.23583984375, 741.547119140625, 0.0), (-1738.945068359375, 714.0121459960938, 0.0)), width=3.5, rise=0.08, scope_ends=()),
        Section(105, 2752, ((-1738.945068359375, 714.0121459960938, 0.0), (-1773.569091796875, 671.6314697265625, 0.0), (-1797.8270263671875, 650.4431762695312, 0.0)), width=3.5, rise=0.08, scope_ends=(1,)),
        Section(106, 2753, ((-1738.945068359375, 714.0121459960938, 0.0), (-1760.6337890625, 700.5819702148438, 0.0), (-1781.3486328125, 693.2964477539062, 0.0), (-1801.807373046875, 684.675048828125, 0.0)), width=2.5, rise=0.08, scope_ends=(1,)),
        Section(190, 2759, ((-1702.1434326171875, 785.702880859375, 0.0), (-1721.97216796875, 775.0662231445312, 0.0), (-1735.4620361328125, 764.2496948242188, 0.0), (-1764.38232421875, 729.9268798828125, 0.0), (-1787.223876953125, 710.4747314453125, 0.0), (-1809.3255615234375, 695.4974975585938, 0.0)), width=3.5, rise=0.08, scope_ends=(1,)),
        Section(189, 2759, ((-1833.9014892578125, 677.749267578125, 0.0), (-1845.401123046875, 665.1846923828125, 0.0)), width=3.5, rise=0.08, scope_ends=(0,)),
    )
    sites = (Intersection(190, (End(93, 1), End(96, 0), End(95, 0)), 10.,
                          (10., 10., 10.), ('STYLE', 'STYLE', 'STYLE')),
             Intersection(210, (End(96, 1), End(106, 0), End(105, 0)), 10.,
                          (10., 3., 10.), ('STYLE', 'STYLE', 'STYLE')))
    bundle = Bundle(4, (BundleMember((Traversal(190, False),)),
                        BundleMember((Traversal(96, False),))))
    return Network('partial_ground_bundle', sections, sites, (bundle,), settings=Settings())


class PartialBundleCapTests(unittest.TestCase):
    def test_conservative_dependencies_include_an_unselected_clipped_branch(self):
        network = bergamo_ground_ramps()
        network = replace(network, sections=tuple(replace(s, selection_ends=(1,)) if s.id == 95 else s
                                                  for s in network.sections))
        # The actual inside cap follows Section 93. An unrelated outside
        # branch reaching the import rectangle must still be a dependency,
        # without being mistaken for permission to close a cap on that cut.
        self.assertIn(95, dependencies(network, network.bundles[0]))
        engine = PythonMeshBackend(network.settings)
        from tested_python_streets.updates import DependencyIndex
        DependencyIndex(tuple(t for provider in engine.providers for t in provider.tasks(network)))
        path, last = terminal_chain(network, End(95, 0), require_terminal=False)
        self.assertEqual([t.section_id for t in path], [95])
        self.assertEqual(last, End(95, 1))
        with self.assertRaisesRegex(ValueError, 'selection cut'):
            terminal_chain(network, End(95, 0))

    def test_actual_unequal_ground_ramps_keep_nearby_roads_and_rounded_cap_curbs(self):
        network = bergamo_ground_ramps()
        with patch('tested_python_streets.finite_free_closure.prepare',
                   side_effect=AssertionError('A successful original terminal walk must stay first')):
            update = PythonMeshBackend(network.settings).prepare(network)
        separator = update.products[Key('separator', 4)]
        self.assertGreater(len(separator.mesh.quads), 0)
        self.assertLess(separator.mesh.validate()['maximum_edge_ratio'], 8.)
        plan = update.products[Key('bundle_plan', 4)]
        self.assertEqual(plan.mode, 'separated')
        self.assertIsNone(update.products[Key('bundle_pavement', 4)])
        for key in (Key('road', 96), Key('road', 190), Key('junction_mesh', 190), Key('junction_mesh', 210)):
            pavement = update.products[key].mesh
            self.assertGreater(len(pavement.quads), 0)
            self.assertLess(pavement.validate()['maximum_edge_ratio'], 8.)
        run, = separator.metadata['runs']
        records = run['cap_continuations']
        self.assertNotIn('finite_free_reconciliation', run)
        self.assertTrue(all('finite_free_closure' not in r for r in records))
        self.assertEqual({r['terminal']['section_id'] for r in records}, {93, 106})
        self.assertEqual({r['terminal']['scope'] for r in records}, {False, True})
        self.assertTrue(all(not r['membership_changed'] for r in records))
        for sid in (93, 95, 105, 106, 189):
            self.assertGreater(len(update.products[Key('road', sid)].mesh.quads), 0)
            self.assertEqual(update.products[Key('profile', sid)].half, 0)
        for sid in (96, 190):
            self.assertNotEqual(update.products[Key('profile', sid)].half, 0)
        self.assertTrue({93, 95, 105, 106} <= dependencies(network, network.bundles[0]))

    def test_terminal_walk_keeps_original_scope_and_rejects_other_junctions_or_selection_cuts(self):
        a = Section(1, 9, ((0., 0., 0.), (10., 0., 0.)))
        b = Section(2, 9, ((10., 0., 0.), (20., 0., 0.)), scope_ends=(1,))
        seam = Intersection(7, (End(1, 1), End(2, 0)))
        network = Network('walk', (a, b), (), connections=(seam,))
        path, terminal = terminal_chain(network, End(1, 0))
        self.assertEqual([t.section_id for t in path], [1, 2])
        self.assertEqual(terminal, End(2, 1))
        invalid = replace(network, sections=(a, replace(b, scope_ends=(), selection_ends=(1,))))
        with self.assertRaisesRegex(ValueError, 'selection cut'):
            terminal_chain(invalid, End(1, 0))
        c = Section(3, 10, ((10., 0., 0.), (10., 10., 0.)))
        fork = Intersection(7, (End(1, 1), End(2, 0), End(3, 0)))
        invalid = Network('fork', (a, b, c), (fork,))
        with self.assertRaisesRegex(ValueError, 'another physical junction'):
            terminal_chain(invalid, End(1, 0))


if __name__ == '__main__':
    unittest.main()
