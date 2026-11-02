"""Provider-specific ground-passage admission; run in Blender, no asset loading."""
from copy import deepcopy
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from blosm.way.street_network.adapter import PreparedNetworkAdapter, StreetDesign
from blosm.way.street_network.python_backend.renderer import PythonStreetRenderer


class GroundSurfaceScopeTests(unittest.TestCase):
    def test_python_provider_admits_ground_passages_without_rewriting_tags(self):
        renderer = PythonStreetRenderer(SimpleNamespace(), design=StreetDesign(), asset_loader=object())
        self.assertTrue(renderer.adapter.ground_passages)
        section = dict(kind='Section', category='residential', tags={'tunnel': 'building_passage'},
                       points=((0., 0.), (20., 0.)), width=6.)
        original = deepcopy(section)
        self.assertIsNone(renderer.adapter.street_support([section]))
        self.assertEqual(section, original)

    def test_default_and_ordinary_tunnel_scope_is_unchanged(self):
        section = dict(kind='Section', category='residential', tags={}, points=((0., 0.), (20., 0.)), width=6.)
        self.assertIsNone(PreparedNetworkAdapter().street_support([section]))
        section['tags']['tunnel'] = 'building_passage'
        self.assertIsNotNone(PreparedNetworkAdapter().street_support([section]))
        for tags in ({'tunnel': 'yes'}, {'tunnel': 'building_passage', 'layer': '-1'},
                     {'tunnel': 'building_passage', 'level': '1'}, {'bridge': 'yes'}):
            section['tags'] = tags
            self.assertIsNotNone(PreparedNetworkAdapter(ground_passages=True).street_support([section]))


if __name__ == '__main__':
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(GroundSurfaceScopeTests))
    if not result.wasSuccessful():
        raise RuntimeError('Ground passage scope regression failed')
