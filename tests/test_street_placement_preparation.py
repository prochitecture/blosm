"""Native Section preparation must preserve source-authored lateral placement."""
import unittest
from mathutils import Vector
from blosm.lib.CompGeom.PolyLine import PolyLine
from blosm.way.item import Section
from blosm.way.way_network import NetSection, WayNetwork
from blosm.way.way_algorithms import createSectionNetwork


def segment(a, b, placement='right_of:1'):
    a, b = Vector(a).freeze(), Vector(b).freeze()
    tags = dict(highway='trunk', oneway='yes', lanes='3')
    if placement is not None:
        tags['placement'] = placement
    result = NetSection(a, b, 'trunk', tags)
    result.sourceSegments = (('source', str(tuple(a)), str(tuple(b)), 0., 1.),)
    return result


class PlacementPreparationTests(unittest.TestCase):
    def prepare(self, segments):
        network = WayNetwork(False)
        for source in segments:
            network.addSegment(source, False)
        return tuple(createSectionNetwork(network, False).iterAllForwardSegments())

    def test_same_constant_placement_retains_a_continuous_section(self):
        rows = self.prepare((segment((0., 0.), (10., 0.)), segment((10., 0.), (20., 0.))))
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(rows[0].path), 3)
        self.assertEqual(len(rows[0].sourceSegments), 2)

    def test_equal_lane_counts_do_not_erase_a_placement_boundary(self):
        for after in ('middle_of:2', None):
            with self.subTest(after=after):
                rows = self.prepare((segment((0., 0.), (10., 0.)),
                                     segment((10., 0.), (20., 0.), after)))
                self.assertEqual(len(rows), 2)
                self.assertEqual(sum(len(r.sourceSegments) for r in rows), 2)

    def test_opposed_original_directions_cannot_share_one_constant_offset(self):
        rows = self.prepare((segment((0., 0.), (10., 0.)), segment((20., 0.), (10., 0.))))
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(len(r.sourceSegments) == 1 for r in rows))

    def test_inverse_traversal_preserves_physical_placement_and_original_controls(self):
        original = segment((0., 0.), (20., 0.))
        forward = Section(original, PolyLine(original.path), None)
        reversed_source = ~original
        backward = Section(reversed_source, PolyLine(reversed_source.path), None)
        for section in (forward, backward):
            before = tuple(section.centerline)
            section.setSectionAttributes(True, 'NNN', '', 0, dict(laneWidth=3.5))
            self.assertEqual(tuple(section.centerline), before)
            self.assertEqual(section.width, 10.5)
        self.assertEqual((forward.offset, backward.offset), (1.75, -1.75))
        self.assertEqual(backward.sourcePlacement['original_direction'], -1)

    def test_reassigning_untagged_dimensions_cannot_retain_a_stale_source_offset(self):
        source = segment((0., 0.), (20., 0.))
        section = Section(source, PolyLine(source.path), None)
        section.setSectionAttributes(True, 'NNN', '', 0, dict(laneWidth=3.5))
        self.assertEqual(section.offset, 1.75)
        section.tags = dict(oneway='yes', lanes='3')
        section.setSectionAttributes(True, 'NNN', '', 0, dict(laneWidth=3.5))
        self.assertEqual(section.offset, 0.)


if __name__ == '__main__':
    unittest.main()
