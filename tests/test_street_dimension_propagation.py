"""Lane corrections must not overwrite a real source change at a junction."""
from types import SimpleNamespace
import unittest
from mathutils import Vector
from test_street_preparation import source
from blosm.way.item import Intersection
from blosm.way.dimension_propagation import TransitionDimensionPropagation
from blosm.way.item.side_lane import SideLane


def minor(before, after):
    node = Intersection(Vector(before.centerline[-1][:2]).freeze())
    node.isMinor = True
    node.arriving = SimpleNamespace(tail=before)
    node.leaving = SimpleNamespace(head=after)
    before.succ = node
    after.pred = node
    return node


class StreetDimensionPropagationTests(unittest.TestCase):
    def test_explicit_source_placement_survives_lane_inference_and_pixel_correction(self):
        before = source((0., 0.), (20., 0.), 'NN').head
        following = source((20., 0.), (50., 0.), 'NNR').head
        following.tags = dict(oneway='yes', placement='right_of:1', lanes='3')
        following.setSectionAttributes(True, 'NNR', '', 0, dict(laneWidth=2.5))
        controls = tuple(following.centerline)
        SideLane(True, following.src, before, following)
        self.assertEqual(following.width, 7.5)
        self.assertEqual(following.offset, 1.25)
        self.assertEqual(following.sourcePlacement['status'], 'RESOLVED')
        self.assertEqual(tuple(following.centerline), controls)

    def test_placement_boundary_stops_inferred_dimension_propagation(self):
        first = source((0., 0.), (20., 0.), 'NNN').head
        following = source((20., 0.), (50., 0.), 'NNN').head
        first.tags = dict(oneway='yes', placement='right_of:1', lanes='3')
        first.setSectionAttributes(True, 'NNN', '', 0, dict(laneWidth=2.5))
        minor(first, following)
        first.width = 7.6
        TransitionDimensionPropagation().apply(SimpleNamespace(succ=first), offset=True)
        self.assertEqual(following.width, 7.5)
        self.assertEqual(following.offset, 0.)

    def test_nonuniform_lane_widths_keep_their_distinct_source_offsets(self):
        first = source((0., 0.), (20., 0.), 'NNN').head
        following = source((20., 0.), (50., 0.), 'NNN').head
        for section, widths in ((first, '2|3|5'), (following, '3|3|4')):
            section.tags = dict(oneway='yes', placement='right_of:1', width='10',
                                **{'width:lanes': widths})
            section.setSectionAttributes(True, 'NNN', '', 0, dict(laneWidth=3.5))
        minor(first, following)
        TransitionDimensionPropagation().apply(SimpleNamespace(succ=first), offset=True)
        self.assertEqual((first.offset, following.offset), (3., 2.))

    def test_minor_lane_loss_keeps_the_original_two_lane_width(self):
        wide = source((0., 0.), (20., 0.), 'NNN').head
        narrow = source((20., 0.), (50., 0.), 'NN').head
        minor(wide, narrow)
        before = narrow.width
        TransitionDimensionPropagation().apply(SimpleNamespace(succ=wide))
        self.assertEqual(wide.sourceWidth, 7.5)
        self.assertEqual(narrow.sourceWidth, 5.)
        self.assertEqual(narrow.width, before)
        self.assertEqual(narrow.widthOrigin, 'Lane count and style/category lane width')

    def test_pixel_and_offset_correction_reaches_equal_source_width_only(self):
        first = source((0., 0.), (20., 0.), 'NN').head
        same = source((20., 0.), (40., 0.), 'NN').head
        changed = source((40., 0.), (60., 0.), 'NNN').head
        later = source((60., 0.), (80., 0.), 'NN').head
        minor(first, same)
        minor(same, changed)
        minor(changed, later)
        first.width = 5.2  # Optional prepared pixel correction; original is 5 m.
        first.offset = 1.75
        TransitionDimensionPropagation().apply(SimpleNamespace(succ=first), offset=True)
        self.assertEqual(same.width, 5.2)
        self.assertEqual(same.offset, 1.75)
        self.assertEqual(same.sourceWidth, 5.)
        self.assertEqual(changed.width, 7.5)
        self.assertEqual(changed.offset, 0.)
        self.assertEqual(later.width, 5.)

    def test_explicit_width_change_stops_even_with_the_same_lane_count(self):
        first = source((0., 0.), (20., 0.), 'NN').head
        following = source((20., 0.), (40., 0.), 'NN').head
        following.tags = dict(following.tags, width='6')
        following.setSectionAttributes(True, 'NN', '', 0, dict(laneWidth=2.5))
        minor(first, following)
        TransitionDimensionPropagation().apply(SimpleNamespace(succ=first), offset=True)
        self.assertEqual(following.width, 6.)
        self.assertEqual(following.sourceWidth, 6.)

    def test_symmetric_correction_does_not_propagate_offsets(self):
        first = source((0., 0.), (20., 0.), 'NN').head
        following = source((20., 0.), (40., 0.), 'NN').head
        minor(first, following)
        first.width, first.offset, following.offset = 5.2, 1.75, .3
        TransitionDimensionPropagation().apply(SimpleNamespace(succ=first))
        self.assertEqual(following.width, 5.2)
        self.assertEqual(following.offset, .3)

    def test_corrupt_item_cycle_is_explicit(self):
        first = source((0., 0.), (20., 0.)).head
        first.succ = first
        with self.assertRaisesRegex(ValueError, 'cyclic original item chain'):
            TransitionDimensionPropagation().apply(SimpleNamespace(succ=first))


if __name__ == '__main__':
    unittest.main()
