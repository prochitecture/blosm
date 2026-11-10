"""Host-only imports and graph operations under legacy/extension namespaces.

No Blender installation, pip dependency or global sys.path change is needed.
Root package shells deliberately avoid importing bpy-dependent addon startup.
"""
import importlib
import pickle
from pathlib import Path
import sys
from types import ModuleType
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


class VendorImportTests(unittest.TestCase):
    def load(self, package):
        before = set(sys.modules)
        names = package.split('.')
        for i in range(1, len(names) + 1):
            name = '.'.join(names[:i])
            if name not in sys.modules:
                shell = ModuleType(name)
                shell.__path__ = [str(ROOT)] if name == package else []
                sys.modules[name] = shell
        def cleanup():
            for name in set(sys.modules) - before:
                if name == package or name.startswith(package + '.') or name in (
                        '.'.join(names[:i]) for i in range(1, len(names))):
                    sys.modules.pop(name, None)
        self.addCleanup(cleanup)
        return importlib.import_module(package + '.lib.blosm_networkx')

    def exercise_graphs(self, nx):
        graph = nx.MultiDiGraph()
        graph.add_edge(1, 2, width=7.)
        graph.add_edge(1, 2, width=4.)
        graph.add_edge(2, 3)
        self.assertEqual(graph.number_of_edges(), 3)
        self.assertEqual(nx.shortest_path(graph, 1, 3), [1, 2, 3])
        self.assertEqual(list(nx.weakly_connected_components(graph)), [{1, 2, 3}])
        self.assertEqual(list(nx.topological_sort(graph)), [1, 2, 3])
        restored = pickle.loads(pickle.dumps(graph))
        self.assertEqual(list(restored.edges(keys=True, data=True)), list(graph.edges(keys=True, data=True)))
        data = nx.node_link_data(graph)
        self.assertEqual(nx.node_link_graph(data).number_of_edges(), 3)
        self.assertEqual(len(nx.graph_atlas_g()), 1253)

    def test_legacy_package(self):
        self.exercise_graphs(self.load('blosm_vendor_legacy_test'))

    def test_nested_extension_with_foreign_lib(self):
        foreign = ModuleType('lib')
        foreign.__path__ = []
        paths = list(sys.path)
        with patch.dict(sys.modules, lib=foreign):
            nx = self.load('bl_ext.vendor_test.blosm')
            self.exercise_graphs(nx)
            self.assertTrue(nx.Graph.__module__.startswith('bl_ext.vendor_test.blosm.'))
            self.assertIs(sys.modules['lib'], foreign)
            self.assertEqual(paths, sys.path)
            self.assertNotIn('lib.blosm_networkx', sys.modules)

    def test_two_extension_repositories_do_not_share_graph_types(self):
        first = self.load('bl_ext.vendor_first.blosm')
        second = self.load('bl_ext.vendor_second.blosm')
        self.assertIsNot(first.Graph, second.Graph)
        self.assertIsNot(first.NetworkXError, second.NetworkXError)
        self.exercise_graphs(first)
        self.exercise_graphs(second)


if __name__ == '__main__':
    unittest.main(verbosity=2)
