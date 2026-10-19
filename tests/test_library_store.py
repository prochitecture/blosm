"""Filesystem-only library retention tests; no Blender or addon import needed."""
import hashlib
import importlib.util
import os
from pathlib import Path
import tempfile
from concurrent.futures import ThreadPoolExecutor
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('asset_store',
    Path(__file__).resolve().parents[1] / 'app' / 'library_store.py')
store_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(store_module)


class AssetStoreTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='blosm_store_tests_',
                                                     dir=os.environ.get('BLOSM_TEST_TMP'))
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / 'source.blend'
        self.payload = b'first library revision' * 100000
        self.digest = hashlib.sha256(self.payload).hexdigest()
        self.source.write_bytes(self.payload)
        self.store = store_module.ImmutableLibraryStore(self.root / 'retained')

    def test_retained_revision_survives_source_upgrade_and_removal(self):
        first = self.store.install(self.source, self.digest)
        replacement = b'new library revision'
        self.source.write_bytes(replacement)
        second_digest = hashlib.sha256(replacement).hexdigest()
        second = self.store.install(self.source, second_digest)
        self.source.unlink()
        self.assertNotEqual(first, second)
        self.assertEqual(self.store.install(self.source, self.digest).read_bytes(), self.payload)
        self.assertEqual(self.store.existing(second_digest).read_bytes(), replacement)

    def test_mismatch_never_publishes_a_partial_revision(self):
        with self.assertRaisesRegex(store_module.AssetStoreError, 'fingerprint mismatch'):
            self.store.install(self.source, '0' * 64)
        self.assertFalse(self.store.path('0' * 64).exists())
        self.assertEqual(list(self.store.root.rglob('*.partial')), [])

    def test_corrupt_retained_revision_is_not_silently_replaced(self):
        target = self.store.install(self.source, self.digest)
        target.write_bytes(b'corrupted')
        with self.assertRaisesRegex(store_module.AssetStoreError, 'damaged or was replaced'):
            self.store.install(self.source, self.digest)
        self.assertEqual(target.read_bytes(), b'corrupted')

    def test_concurrent_installers_reuse_one_complete_revision(self):
        with ThreadPoolExecutor(max_workers=4) as executor:
            paths = list(executor.map(lambda _: self.store.install(self.source, self.digest), range(4)))
        self.assertEqual(len(set(paths)), 1)
        self.assertEqual(paths[0].read_bytes(), self.payload)
        self.assertEqual(list(self.store.root.rglob('*.partial')), [])

    def test_failed_atomic_publication_leaves_no_partial_file(self):
        with patch.object(store_module, '_publish', side_effect=OSError('unsupported filesystem')):
            with self.assertRaisesRegex(store_module.AssetStoreError, 'Cannot retain library'):
                self.store.install(self.source, self.digest)
        self.assertFalse(self.store.path(self.digest).exists())
        self.assertEqual(list(self.store.root.rglob('*.partial')), [])

    def test_store_can_be_relocated_without_changing_revision_layout(self):
        installed = self.store.install(self.source, self.digest)
        relocated = self.root / 'moved'
        self.store.root.rename(relocated)
        moved_store = store_module.ImmutableLibraryStore(relocated)
        self.assertEqual(moved_store.existing(self.digest).relative_to(relocated),
                         installed.relative_to(self.store.root))

    def test_invalid_digest_cannot_become_a_path(self):
        for value in ('../outside', '', 'A'*64, '1'*63, None):
            with self.assertRaises(ValueError):
                self.store.path(value)


if __name__ == '__main__':
    unittest.main()
