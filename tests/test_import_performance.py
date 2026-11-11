"""Host-only failure and policy tests; no Windows or Blender changes are made."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('blosm_import_performance_test', ROOT / 'app/import_performance.py')
performance = importlib.util.module_from_spec(spec)
spec.loader.exec_module(performance)


class Policy:
    def __init__(self, original=(0, 0)):
        self.current = original
        self.writes = []
        self.read_error = None
        self.write_error = None

    def read(self):
        if self.read_error:
            error, self.read_error = self.read_error, None
            raise error
        return self.current

    def write(self, value):
        self.writes.append(value)
        self.current = value
        if self.write_error:
            error, self.write_error = self.write_error, None
            raise error


class ImportPerformanceTests(unittest.TestCase):
    def setUp(self):
        self.policy = Policy()
        self.platform = patch.object(performance.os, 'name', 'nt')
        self.factory = patch.object(performance, '_WindowsProcessPolicy', return_value=self.policy)
        self.platform.start()
        self.api = self.factory.start()
        self.addCleanup(self.platform.stop)
        self.addCleanup(self.factory.stop)

    def test_request_restores_and_preserves_other_bits(self):
        self.policy.current = (2, 2)
        with performance.ImportPerformance() as request:
            self.assertTrue(request.applied)
            self.assertEqual(self.policy.current, (3, 2))
        self.assertEqual(self.policy.current, (2, 2))
        request.close()
        self.assertEqual(self.policy.writes, [(3, 2), (2, 2)])

    def test_disabled_and_unsupported_never_load_windows_api(self):
        with performance.ImportPerformance(False) as request:
            self.assertEqual(request.status, 'disabled')
        with patch.object(performance.os, 'name', 'posix'), performance.ImportPerformance() as request:
            self.assertEqual(request.status, 'unsupported')
        self.api.assert_not_called()

    def test_explicit_eco_and_existing_high_are_untouched(self):
        for original, expected in (((3, 3), 'power_saving'), ((3, 2), 'already_high')):
            self.policy.current = original
            with performance.ImportPerformance() as request:
                self.assertEqual(request.status, expected)
                self.assertFalse(request.applied)
            self.assertEqual(self.policy.current, original)
        self.assertEqual(self.policy.writes, [])

    def test_nested_requests_restore_only_at_outer_exit(self):
        with performance.ImportPerformance():
            with performance.ImportPerformance() as inner:
                self.assertEqual(inner.status, 'already_high')
            self.assertEqual(self.policy.current, (1, 0))
        self.assertEqual(self.policy.writes, [(1, 0), (0, 0)])

    def test_body_error_and_cancellation_propagate_after_restoration(self):
        for error in (RuntimeError('failed import'), KeyboardInterrupt()):
            with self.assertRaises(type(error)):
                with performance.ImportPerformance():
                    raise error
            self.assertEqual(self.policy.current, (0, 0))

    def test_unavailable_api_and_read_failure_continue(self):
        for error in (AttributeError('missing API'), OSError('unsupported policy')):
            self.api.side_effect = error
            with performance.ImportPerformance() as request:
                self.assertEqual(request.status, 'unavailable')
        self.api.side_effect = None
        self.policy.read_error = OSError('read denied')
        with performance.ImportPerformance() as request:
            self.assertEqual(request.status, 'unavailable')
        self.assertEqual(self.policy.writes, [])

    def test_failed_request_rolls_back_even_if_write_changed_policy(self):
        self.policy.write_error = OSError('request denied')
        with performance.ImportPerformance() as request:
            self.assertFalse(request.applied)
            self.assertEqual(request.status, 'unavailable')
            self.assertEqual(self.policy.current, (0, 0))
        self.assertEqual(self.policy.writes, [(1, 0), (0, 0)])

    def test_interruption_while_entering_rolls_back_and_propagates(self):
        self.policy.write_error = KeyboardInterrupt()
        with self.assertRaises(KeyboardInterrupt):
            with performance.ImportPerformance():
                self.fail('interrupted request entered')
        self.assertEqual(self.policy.current, (0, 0))

    def test_unconfirmed_request_is_rolled_back(self):
        with patch.object(self.policy, 'read', side_effect=[(0, 0), (0, 0), (0, 0)]):
            with performance.ImportPerformance() as request:
                self.assertFalse(request.applied)
        self.assertEqual(self.policy.writes, [(1, 0), (0, 0)])

    def test_restore_failure_does_not_mask_import_error_and_can_retry(self):
        with self.assertLogs(performance._logger, level='WARNING'), self.assertRaisesRegex(RuntimeError, 'import failed'):
            with performance.ImportPerformance() as request:
                self.policy.write_error = OSError('restore denied')
                raise RuntimeError('import failed')
        self.assertEqual(request.status, 'restore_failed')
        request.close()
        self.assertEqual(request.status, 'restored')
        self.assertFalse(request._restore_needed)

    def test_unconfirmed_restoration_is_reported(self):
        request = performance.ImportPerformance()
        request.__enter__()
        with patch.object(self.policy, 'read', return_value=(1, 0)), \
                self.assertLogs(performance._logger, level='WARNING'):
            request.close()
        self.assertEqual(request.status, 'restore_failed')
        request.close()
        self.assertFalse(request._restore_needed)

    def test_shared_preference_and_missing_entry_default(self):
        context = SimpleNamespace(preferences=SimpleNamespace(addons={}))
        self.assertTrue(performance.import_performance(context, 'bl_ext.repo.blosm').enabled)
        prefs = SimpleNamespace(preferPerformanceDuringImport=False)
        context.preferences.addons['bl_ext.repo.blosm'] = SimpleNamespace(preferences=prefs)
        self.assertFalse(performance.import_performance(context, 'bl_ext.repo.blosm').enabled)


if __name__ == '__main__':
    unittest.main(verbosity=2)
