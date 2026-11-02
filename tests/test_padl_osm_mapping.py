"""OSM mapping semantics without registering Blender or importing a backend."""
import importlib.util
from pathlib import Path
import sys
import unittest

PATH = Path(__file__).resolve().parents[1] / 'building/padl_integration/mapping.py'
spec = importlib.util.spec_from_file_location('_padl_osm_mapping_tested', PATH)
mapping = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mapping
spec.loader.exec_module(mapping)


class OsmMappingTests(unittest.TestCase):
    def setUp(self):
        self.mapper = mapping.OsmFaceAttributeMapper()

    def test_osm_heights_are_separate_tmin_offsets(self):
        result = self.mapper.map({'height': '24 m', 'min_height': '3', 'roof:height': '4',
                                  'building:levels': '6', 'roof:shape': 'hipped'})
        self.assertTrue(result.accepted)
        self.assertEqual(result.attributes, {'building_use': 'unknown', 'min_elevation': 3,
            'max_elevation': 24, 'roof_height': 4, 'levels_wall_count': 6, 'roof_type': 'hipped'})

    def test_units(self):
        for raw, expected in [('12.5', 12.5), ('12.5 m', 12.5), ('20 ft', 6.096),
                              ('20 feet', 6.096), ("20'", 6.096), ('20\' 6"', 6.2484), ('.5', .5)]:
            with self.subTest(raw=raw):
                self.assertAlmostEqual(mapping.length_metres(raw), expected)

    def test_invalid_values_do_not_become_authoritative_overrides(self):
        for raw in ('12;14', '12-14', '~12', 'nan', 'inf', '-2', '10 storeys', '', '12,5', '5\' 12"'):
            with self.subTest(raw=raw):
                result = self.mapper.map({'height': raw})
                self.assertNotIn('max_elevation', result.attributes)
                self.assertEqual(result.issues[0].code, 'invalid_tag')

    def test_zero_and_absent_are_distinct(self):
        self.assertEqual(self.mapper.map({'roof:height': '0'}).attributes['roof_height'], 0)
        self.assertNotIn('roof_height', self.mapper.map({}).attributes)
        self.assertNotIn('max_elevation', self.mapper.map({'height': '0'}).attributes)

    def test_no_default_policy_values_are_sourced_from_osm(self):
        attrs = self.mapper.map({'building': 'apartments', 'ele': '80', 'levels': '7', 'level': '2'}).attributes
        self.assertEqual(attrs, {'building_use': 'apartments', 'min_elevation': 0})

    def test_supported_uses_and_house_forms(self):
        for form in ('house', 'detached', 'semidetached_house', 'terrace', 'bungalow'):
            self.assertEqual(self.mapper.map({'building': form}).attributes['building_use'], 'house')
        for form in ('apartments', 'office'):
            self.assertEqual(self.mapper.map({'building': form}).attributes['building_use'], form)

    def test_residential_is_not_assumed_to_be_apartments(self):
        self.assertEqual(self.mapper.map({'building': 'residential'}).attributes['building_use'], 'unknown')

    def test_current_use_takes_priority(self):
        self.assertEqual(self.mapper.map({'building': 'warehouse', 'building:use': 'apartments'}).attributes['building_use'], 'apartments')
        self.assertEqual(self.mapper.map({'building': 'office', 'building:use': 'school'}).attributes['building_use'], 'unknown')

    def test_office_poi_does_not_reclassify_apartment_block(self):
        self.assertEqual(self.mapper.map({'building': 'apartments', 'office': 'lawyer'}).attributes['building_use'], 'apartments')
        self.assertEqual(self.mapper.map({'building': 'yes', 'office': 'company'}).attributes['building_use'], 'office')
        self.assertEqual(self.mapper.map({'building': 'yes', 'office': 'no'}).attributes['building_use'], 'unknown')

    def test_unmapped_functions_stay_unknown(self):
        self.assertEqual(self.mapper.map({'building': 'yes', 'amenity': 'school'}).attributes['building_use'], 'unknown')
        self.assertEqual(self.mapper.map({'building': 'yes', 'man_made': 'water_tower'}).attributes['building_use'], 'unknown')

    def test_flat_roof_ignores_roof_height(self):
        self.assertNotIn('roof_height', self.mapper.map({'roof:shape': 'flat', 'roof:height': '2'}).attributes)
        self.assertEqual(self.mapper.map({'roof:shape': 'flat', 'roof:height': 'invalid'}).issues, [])

    def test_unsupported_roof_is_not_silently_flat(self):
        result = self.mapper.map({'roof:shape': 'dome'})
        self.assertFalse(result.accepted)
        self.assertEqual(result.issues[0].code, 'unsupported_roof')
        self.assertNotIn('roof_type', result.attributes)

    def test_only_positive_integral_wall_counts(self):
        for raw in ('0', '-1', '3.5', '3;4', '2147483648'):
            self.assertNotIn('levels_wall_count', self.mapper.map({'building:levels': raw}).attributes)

    def test_roof_and_basement_counts_are_not_wall_counts(self):
        self.assertNotIn('levels_wall_count', self.mapper.map({'roof:levels': '2', 'building:levels:underground': '3'}).attributes)

    def test_skipped_storeys_require_an_adapter(self):
        self.assertFalse(self.mapper.map({'building:min_level': '2', 'building:levels': '5'}).accepted)
        self.assertTrue(self.mapper.map({'building:min_level': '0', 'building:levels': '5'}).accepted)

    def test_rules_can_be_extended_without_changing_a_backend(self):
        class TowerRule(mapping.TagRule):
            def apply(self, tags, result):
                if tags.get('man_made') == 'tower':
                    result.attributes['wall_min_height'] = 18.0
        mapper = mapping.OsmFaceAttributeMapper((*self.mapper.rules, TowerRule()))
        self.assertEqual(mapper.map({'man_made': 'tower'}).attributes['wall_min_height'], 18)


if __name__ == '__main__':
    unittest.main()
