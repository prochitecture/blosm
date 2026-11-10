"""Extension packaging and network-policy checks that run without Blender."""
import ast
from pathlib import Path
import sys
import tomllib
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from release_target import ADDON_ROOT as ROOT
sys.path.insert(0, str(ROOT))
from downloads.access import OnlineAccessError, require_online_access
from downloads.transport import DownloadTransport


class ExtensionContractChecks(unittest.TestCase):
    def test_manifest_metadata(self):
        manifest = tomllib.loads((ROOT / 'blender_manifest.toml').read_text(encoding='utf-8'))
        self.assertEqual(manifest['id'], 'blosm')
        self.assertEqual(manifest['blender_version_min'], '5.2.0')
        self.assertLessEqual(len(manifest['tagline']), 64)
        self.assertIn('SPDX:GPL-3.0-or-later', manifest['license'])
        self.assertEqual(set(manifest['permissions']), {'network', 'files', 'clipboard'})
        tree = ast.parse((ROOT / '__init__.py').read_text(encoding='utf-8'))
        self.assertFalse(any(isinstance(node, ast.Name) and node.id == 'bl_info' for node in ast.walk(tree)))

    def test_offline_blocks_helper_before_creating_resources(self):
        fake_bpy = SimpleNamespace(app=SimpleNamespace(online_access=False))
        with patch.dict(sys.modules, bpy=fake_bpy), \
                patch('downloads.transport.socket.socket') as socket, \
                patch('downloads.transport.subprocess.Popen') as process:
            with self.assertRaises(OnlineAccessError):
                DownloadTransport()
            socket.assert_not_called()
            process.assert_not_called()

    def test_runtime_imports_stay_in_the_extension_namespace(self):
        # Include namespace packages and nested modules, not only root folders
        # with __init__.py. The old check missed lib.blosm_networkx, debug and
        # the generated parser's top-level fallback. Stdlib imports such as
        # random/typing remain valid despite same-named internal files.
        runtime_paths = [path for path in ROOT.rglob('*.py')
                         if path.relative_to(ROOT).parts[0] not in ('tests', 'tools', 'docs', '.git')]
        roots = {path.relative_to(ROOT).parts[0] for path in runtime_paths
                 if len(path.relative_to(ROOT).parts) > 1}
        roots.update(path.stem for path in runtime_paths)
        roots.difference_update(sys.stdlib_module_names)
        for path in runtime_paths:
            relative = path.relative_to(ROOT)
            if relative.parts[0] in ('tests', 'tools', 'docs', '.git') or relative.as_posix() == 'downloads/worker.py':
                continue  # The download helper deliberately runs outside Blender.
            tree = ast.parse(path.read_text(encoding='utf-8-sig'), filename=str(relative))
            for node in ast.walk(tree):
                names = []
                if isinstance(node, ast.ImportFrom) and not node.level:
                    names = [node.module]
                elif isinstance(node, ast.Import):
                    names = [alias.name for alias in node.names]
                for name in names:
                    self.assertNotIn(name.split('.')[0], roots, f'{relative}:{node.lineno}: {name}')
                if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                        and node.func.attr in ('files', 'import_module') and node.args
                        and isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str)):
                    name = node.args[0].value
                    self.assertNotIn(name.split('.')[0], roots,
                                     f'{relative}:{node.lineno}: dynamic package reference {name}')

    def test_revoking_online_access_closes_existing_transport(self):
        fake_bpy = SimpleNamespace(app=SimpleNamespace(online_access=False))
        for action in (lambda t: t.submit('https://example.invalid'), lambda t: t.poll()):
            transport = object.__new__(DownloadTransport)
            transport.close = Mock()
            with patch.dict(sys.modules, bpy=fake_bpy), self.assertRaises(OnlineAccessError):
                action(transport)
            transport.close.assert_called_once_with()

    def test_online_and_standalone_access(self):
        with patch.dict(sys.modules, bpy=SimpleNamespace(app=SimpleNamespace(online_access=True))):
            require_online_access()
        with patch.dict(sys.modules, bpy=None):
            require_online_access()


if __name__ == '__main__':
    unittest.main()
