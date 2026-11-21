"""OSM placement specifies the way reference, not an invented source shift."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('street_source_placement',
    Path(__file__).resolve().parents[1] / 'way/source_placement.py')
policy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(policy)


class SourcePlacementTests(unittest.TestCase):
    def resolve(self, value, *, width=10.5, lanes=3, **tags):
        return policy.constant_placement(dict(oneway='yes', placement=value, **tags), width, lanes)

    def test_lane_reference_is_measured_from_the_left_in_original_way_direction(self):
        expected = {'left_of:1': 5.25, 'middle_of:1': 3.5, 'right_of:1': 1.75,
                    'left_of:2': 1.75, 'middle_of:2': 0., 'right_of:2': -1.75,
                    'left_of:3': -1.75, 'middle_of:3': -3.5, 'right_of:3': -5.25}
        for value, offset in expected.items():
            with self.subTest(value=value):
                result = self.resolve(value)
                self.assertEqual(result['status'], 'RESOLVED')
                self.assertEqual(result['offset'], offset)

    def test_explicit_width_and_nonuniform_lane_widths_control_the_reference(self):
        self.assertEqual(self.resolve('right_of:1', width=12.)['offset'], 2.)
        result = self.resolve('right_of:1', width=10., **{'width:lanes': '2|3 m|5'})
        self.assertEqual(result['offset'], 3.)
        self.assertEqual(result['lane_widths'], (2., 3., 5.))

    def test_reversing_prepared_traversal_reverses_offset_only(self):
        tags = dict(oneway='yes', placement='right_of:1', lanes='3')
        before = deepcopy(tags)
        forward = policy.constant_placement(tags, 10.5, 3, direction=1)
        backward = policy.constant_placement(tags, 10.5, 3, direction=-1)
        self.assertEqual(backward['offset'], -forward['offset'])
        self.assertEqual(tags, before)
        tags['oneway'] = '-1'  # OSM-way direction still defines lane numbering.
        self.assertEqual(policy.constant_placement(tags, 10.5, 3)['offset'], 1.75)

    def test_equal_endpoint_declarations_are_constant(self):
        result = policy.constant_placement(dict(oneway='yes', **{
            'placement:start': 'right_of:1', 'placement:end': 'right_of:1'}), 10.5, 3)
        self.assertEqual(result['offset'], 1.75)

    def test_variable_directional_and_invalid_placement_stay_unresolved(self):
        variants = [dict(placement='transition'), dict(placement='right_of:4'),
                    dict(placement='middle_of:0'), dict(placement='right_of:1.5'),
                    dict(placement='right_of:1', **{'placement:end': 'middle_of:2'}),
                    {'placement:start': 'middle_of:1'}, {'placement:forward': 'right_of:1'},
                    dict(placement='right_of:1', oneway='no'),
                    dict(placement='right_of:1', **{'width:lanes': '2||3'}),
                    dict(placement='right_of:1', **{'width:lanes': '2|3|4'})]
        for extra in variants:
            with self.subTest(tags=extra):
                result = policy.constant_placement(dict(dict(oneway='yes'), **extra), 10.5, 3)
                self.assertEqual(result['status'], 'UNRESOLVED')
                self.assertNotIn('offset', result)

    def test_untagged_and_ambiguous_source_directions_are_distinct(self):
        self.assertEqual(policy.constant_placement({}, 6., 2), dict(status='NOT_TAGGED', offset=0.))
        for direction in (0, True, 1., -2):
            with self.subTest(direction=direction):
                self.assertEqual(policy.constant_placement(dict(oneway='yes', placement='right_of:1'),
                                                          10.5, 3, direction=direction)['status'], 'UNRESOLVED')


if __name__ == '__main__':
    unittest.main()
