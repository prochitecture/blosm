"""Standalone date/model filtering checks; no Blender or network dependency."""
from pathlib import Path
import sys
import unittest

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mapillary.filter_data import CameraFilter, photo_category, capture_timestamp, date_bound, DAY_MS


class FilterTests(unittest.TestCase):
    def test_inclusive_utc_days(self):
        start = date_bound('2024-02-29', 'From')
        spec = CameraFilter.build('2024-02-29', '2024-02-29')
        for stamp, expected in ((start-1, False), (start, True),
                                (start+DAY_MS-1, True), (start+DAY_MS, False)):
            self.assertEqual(spec.matches(str(stamp), 'perspective'), expected)

    def test_open_bounds_and_calendar_extremes(self):
        self.assertTrue(CameraFilter.build(date_to='1970-01-01').matches(0, 'spherical'))
        self.assertTrue(CameraFilter.build(date_from='1970-01-01').matches(0, 'spherical'))
        self.assertFalse(CameraFilter.build(date_from='1970-01-02').matches(0, 'spherical'))
        self.assertTrue(CameraFilter.build('0001-01-01', '9999-12-31').matches(-1, 'fisheye'))
        self.assertTrue(CameraFilter.build(date_to='9999-12-31').matches(
            date_bound('9999-12-31', 'To')+DAY_MS-1, 'perspective'))

    def test_invalid_dates_rejected(self):
        for bad in ('2023-02-29', '2024-02-30', '0000-01-01', '01/02/2024', '2024-1-01', '2024-01-01T12:00'):
            with self.assertRaises(ValueError):
                CameraFilter.build(bad)
        with self.assertRaisesRegex(ValueError, 'From date'):
            CameraFilter.build('2024-01-02', '2024-01-01')

    def test_undated_policy(self):
        for bad in (None, '', 'unknown', True, [], {}, float('nan'), float('inf'), '1e9999'):
            self.assertIsNone(capture_timestamp(bad))
            self.assertTrue(CameraFilter.build().matches(bad, 'unknown'))
            self.assertFalse(CameraFilter.build('2020-01-01').matches(bad, 'perspective'))
            self.assertTrue(CameraFilter.build('2020-01-01', include_undated=True).matches(bad, 'perspective'))
        self.assertEqual(capture_timestamp('1700000000000.0'), 1700000000000)

    def test_model_categories_are_disjoint(self):
        expected = {'perspective': 'ordinary', 'spherical': 'panoramic',
                    'equirectangular': 'panoramic', 'fisheye': 'fisheye', 'brown': 'other',
                    'unknown': 'other', '': 'other', None: 'other', ' SPHERICAL ': 'panoramic'}
        for model, category in expected.items():
            self.assertEqual(photo_category(model), category)
            for selected in ('ordinary', 'panoramic', 'fisheye', 'other'):
                self.assertEqual(CameraFilter.build(categories=[selected]).matches(None, model), selected == category)
        self.assertEqual(photo_category([]), 'other')
        self.assertFalse(CameraFilter.build(categories=[]).matches(0, 'perspective'))

    def test_date_and_type_are_combined(self):
        spec = CameraFilter.build('2023-11-14', '2023-11-14', ['panoramic'])
        self.assertTrue(spec.matches(1700000000000, 'spherical'))
        self.assertFalse(spec.matches(1700000000000, 'perspective'))
        self.assertFalse(spec.matches(1700000000000-DAY_MS, 'spherical'))


if __name__ == '__main__':
    unittest.main()
