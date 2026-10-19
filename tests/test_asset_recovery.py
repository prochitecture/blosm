"""Filesystem and reload-boundary tests; no Blender process required."""
import hashlib
import importlib
import os
from pathlib import Path
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import Mock, patch

package = ModuleType('blosm_asset_recovery_tests')
package.__path__ = [str(Path(__file__).resolve().parents[1]/'app')]
sys.modules[package.__name__] = package
recovery = importlib.import_module(package.__name__+'.asset_recovery')


class AssetRecoveryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='blosm_recovery_', dir=os.environ.get('BLOSM_TEST_TMP'))
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root/'backup.blend'
        self.payload = b'Exact saved self-contained library revision'*30000
        self.source.write_bytes(self.payload)
        self.digest = hashlib.sha256(self.payload).hexdigest()
        self.target = self.root/'retained'/self.digest/'library.blend'
        self.library = SimpleNamespace(name='Saved library', filepath=str(self.target),
                                       parent=None, reload=Mock(), users_id=[])
        fake = SimpleNamespace(path=SimpleNamespace(abspath=lambda path, library=None: path))
        context = patch.dict(sys.modules, bpy=fake)
        context.start(); self.addCleanup(context.stop)

    def test_missing_revision_restores_exact_bytes_and_reloads_same_reference(self):
        result = recovery.restore_library(self.library, self.source, expected_digest=self.digest)
        self.assertEqual(result, self.target)
        self.assertEqual(result.read_bytes(), self.payload)
        self.assertEqual(self.library.filepath, str(self.target))
        self.library.reload.assert_called_once_with()

    def test_wrong_revision_cannot_publish_or_reload(self):
        wrong = self.root/'other.blend'; wrong.write_bytes(b'New incompatible revision')
        with self.assertRaisesRegex(recovery.AssetStoreError, 'fingerprint mismatch'):
            recovery.restore_library(self.library, wrong)
        self.assertFalse(self.target.exists())
        self.library.reload.assert_not_called()
        self.assertEqual(list(self.root.rglob('*.partial')), [])

    def test_changed_reference_while_file_dialog_is_open_is_rejected(self):
        with self.assertRaisesRegex(recovery.AssetStoreError, 'revision changed'):
            recovery.restore_library(self.library, self.source, expected_digest='0'*64)
        self.assertFalse(self.target.parent.exists())
        self.library.reload.assert_not_called()

    def test_corrupt_original_is_preserved_and_an_explicit_new_store_can_recover(self):
        self.target.parent.mkdir(parents=True); self.target.write_bytes(b'damaged')
        with self.assertRaisesRegex(recovery.AssetStoreError, 'damaged or was replaced'):
            recovery.restore_library(self.library, self.source)
        self.library.reload.assert_not_called()
        moved = recovery.restore_library(self.library, self.source, root=self.root/'new store')
        self.assertEqual(moved, self.root/'new store'/self.digest/'library.blend')
        self.assertEqual(moved.read_bytes(), self.payload)
        self.assertEqual(self.target.read_bytes(), b'damaged')
        self.assertEqual(self.library.filepath, str(moved))

    def test_reload_failure_preserves_verified_copy_and_previous_link(self):
        self.library.reload.side_effect = RuntimeError('native reload failed')
        with self.assertRaisesRegex(recovery.AssetStoreError, 'could not reload'):
            recovery.restore_library(self.library, self.source, root=self.root/'new store')
        self.assertEqual(self.library.filepath, str(self.target))
        self.assertEqual((self.root/'new store'/self.digest/'library.blend').read_bytes(), self.payload)

    def test_unresolved_datablocks_are_not_reported_as_success(self):
        self.library.users_id = [SimpleNamespace(name='Missing group', is_missing=True)]
        with self.assertRaisesRegex(recovery.AssetStoreError, 'still missing'):
            recovery.restore_library(self.library, self.source)

    def test_only_retained_layout_is_offered_for_exact_revision_recovery(self):
        for path in (self.root/'other.blend', self.root/'invalid digest'/'library.blend',
                     Path(self.digest)/'library.blend'):
            self.assertIsNone(recovery.RetainedLibraryReference.inspect('unmanaged', path))
        self.assertEqual(recovery.inspect_library(self.library).digest, self.digest)
        self.library.filepath = str(self.root/'unmanaged.blend')
        with self.assertRaisesRegex(recovery.AssetStoreError, 'not a retained'):
            recovery.restore_library(self.library, self.source)
        self.library.reload.assert_not_called()


if __name__ == '__main__': unittest.main()
