"""Regression: installed open strips retain source sampling; corridors keep row repair.

Run with factory-startup Blender. This reads the bundled declarations and uses
the common asset loader; no candidate override or graph generation is allowed.
"""
from pathlib import Path
import argparse
from dataclasses import replace
import sys
import unittest
import bpy

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--blosm', type=Path, default=Path(__file__).resolve().parents[1])
parser.add_argument('--library', type=Path)
parser.add_argument('--revision')
OPTIONS = parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
if bool(OPTIONS.library) != bool(OPTIONS.revision):
    parser.error('A candidate path and revision are required together')
sys.path.insert(0, str(OPTIONS.blosm.parent))
from blosm.app.asset_library import BlenderAssetLoader
from blosm.way.street_network import assets
from blosm.way.street_network.patch_programs import OriginalPatchCatalog
from blosm.way.gn.editable_input.patch_program_selection import KINDS


class StreetSamplingPolicyTests(unittest.TestCase):
    def test_open_and_closed_families_use_separate_shared_charts(self):
        if OPTIONS.library:
            assets.ORIGINAL_LIBRARY = replace(assets.ORIGINAL_LIBRARY,
                path=str(OPTIONS.library.resolve()), source='external', revision=OPTIONS.revision)
        loader = BlenderAssetLoader(OPTIONS.blosm/'assets', assets.catalog())
        references = (assets.ORIGINAL_SOLVER, assets.ORIGINAL_PATCH_PROVIDER, assets.ORIGINAL_STRIP_PROVIDER)
        descriptor = OriginalPatchCatalog.load_optional(loader)
        self.assertIsNotNone(descriptor)
        references += tuple(descriptor.references(set(KINDS)).values())
        roots = loader.load_many(references).values()
        pending = list(roots)
        trees = {}
        while pending:
            tree = pending.pop()
            if tree.as_pointer() in trees:
                continue
            trees[tree.as_pointer()] = tree
            pending.extend(n.node_tree for n in tree.nodes if n.type == 'GROUP' and n.node_tree)
        prefix = 'Prochitecture Streets / '
        chart_name = prefix+'Original corridor density chart v2'
        open_name = chart_name+' / Preserved open rows v1'
        bridge_name = prefix+'Conforming boundary quad bridge v3'
        open_families = ('Original curved one-ended approach v1', 'Unattached profiled Street v1',
                         'Category-bounded profiled run v1')
        closed_families = ('Original curved corridor with density transitions v1',
                           'Original Street dimensions physical corridor v1',
                           'Finite original seam transition v1')
        charts = {True: set(), False: set()}
        seen = set()
        for tree in trees.values():
            calls = [n for n in tree.nodes if n.type == 'GROUP' and n.node_tree
                     and n.node_tree.name in (chart_name, open_name)]
            if not calls:
                continue
            family = next((f for f in open_families+closed_families
                           if tree.name.startswith(prefix+f)), None)
            self.assertIsNotNone(family, tree.name)
            opened = family in open_families
            seen.add(family)
            for node in calls:
                chart = node.node_tree
                charts[opened].add(chart.as_pointer())
                self.assertEqual(chart.name, open_name if opened else chart_name)
                self.assertEqual(chart.get('st_uniform_equal_corridor_rows', 0), 0 if opened else 1)
                bridges = [n.node_tree for n in chart.nodes if n.type == 'GROUP' and n.node_tree
                           and n.node_tree.name.startswith(bridge_name)]
                self.assertEqual(len(bridges), 1, chart.name)
                self.assertEqual(bridges[0].name, bridge_name if opened else
                                 bridge_name+' / Uniform equal-count rows v1')
        self.assertEqual(seen, set(open_families+closed_families))
        self.assertEqual([len(charts[k]) for k in (True, False)], [1, 1])
        self.assertTrue(charts[True].isdisjoint(charts[False]))


if __name__ == '__main__':
    unittest.main(argv=[sys.argv[0]], verbosity=2)
