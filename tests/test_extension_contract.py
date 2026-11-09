"""Extension packaging and network-policy checks that run without Blender."""
import ast
from pathlib import Path
import sys
import tomllib
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from downloads.access import OnlineAccessError, require_online_access
from downloads.transport import DownloadTransport


class ExtensionContractChecks(unittest.TestCase):
    def test_manifest_metadata(self):
        manifest = tomllib.loads((ROOT / 'blender_manifest.toml').read_text(encoding='utf-8'))
        self.assertEqual(manifest['id'], 'blosm')
        self.assertEqual(manifest['blender_version_min'], '4.5.0')
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
