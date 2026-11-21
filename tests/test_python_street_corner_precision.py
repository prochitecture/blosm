"""Accepted circular contacts refine their roots without rejecting a stall."""
import math
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.corner_patch import corner_patch
from tested_python_streets.intersections import Approach, _local_corner
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.regions import Region, RegionUnion, signed_area


def original_contacts(angle=0., offset=(0., 0.)):
    # Actual finite source controls and endpoint directions at Bergamo1198.
    controls = (((1664.4317626953125, 947.7738647460938),
                 (1671.4512939453125, 932.781005859375)),
                ((1664.4317626953125, 947.7738647460938),
                 (1663.8153076171875, 955.5772094726562)))
    directions = (((.42401924954921716, -.9056531764487544),
                   (.43045749436560093, -.9026108494497995)),
                  ((-.07875346607576854, .9968941225531691),
                   (-.07875346607576854, .9968941225531691)))
    c, s = math.cos(angle), math.sin(angle)
    rotation = np.array(((c, -s), (s, c)))
    result = []
    for identity, points, tangent in zip((600, 599), controls, directions):
        xy = np.asarray(points) @ rotation.T+offset
        section = model.Section(identity, identity, tuple((x, y, 0.) for x, y in xy), width=5.)
        endpoint = {i: rotation @ np.asarray(value) for i, value in enumerate(tangent)}
        result.append(Approach(Path(section, endpoint_directions=endpoint), 0, Crown(.08, 0.)))
    return tuple(result)


def residual(corner):
    a, b = corner.first, corner.second
    distances = (a.offset-.5*a.width-corner.radius, b.offset+.5*b.width+corner.radius)
    centers = [arm.point(station, distance) for arm, station, distance in
               zip((a, b), corner.stations, distances)]
    return float(np.linalg.norm(centers[0]-centers[1]))


class CornerPrecisionTests(unittest.TestCase):
    def test_actual_accepted_contact_refines_before_exact_source_union(self):
        for angle, offset in ((0., (0., 0.)), (.713, (1000., -1100.))):
            with self.subTest(angle=angle):
                first, second = original_contacts(angle, offset)
                original = [a.source.points for a in (first, second)]
                corner = _local_corner(first, second, 2.)
                self.assertLess(residual(corner), 1e-10)
                self.assertEqual(corner.radius, 2.)
                self.assertEqual(original, [a.source.points for a in (first, second)])
                ends = corner.center+corner.radius*np.column_stack((
                    np.cos((corner.start_angle, corner.start_angle+corner.sweep)),
                    np.sin((corner.start_angle, corner.start_angle+corner.sweep))))
                exact = (first.point(corner.stations[0], -.5*first.width),
                         second.point(corner.stations[1], .5*second.width))
                np.testing.assert_allclose(ends, exact, atol=1e-10, rtol=0.)
                regions = [corner_patch(model.Key('cycle_corner', 11), corner, .5, .00025).region]
                for arm, contact in zip((first, second), corner.stations):
                    at = np.unique(np.r_[np.linspace(0., arm.length, 51), contact])
                    xy = np.vstack((arm.point(at, .5*arm.width), arm.point(at, -.5*arm.width)[::-1]))
                    if signed_area(xy) < 0.:
                        xy = xy[::-1]
                    regions.append(Region(model.Key('cycle_strip', arm.source.id), xy))
                loop, = RegionUnion().build(regions)
                published = loop.xy.astype(np.float32)
                self.assertFalse(np.any(np.all(published == np.roll(published, -1, axis=0), axis=1)))

    def test_a_previously_accepted_root_survives_a_floating_point_stall(self):
        first, second = original_contacts()
        seed = np.array((.8144450597185512, .8158651539997414))
        # A singular numerical derivative prevents further improvement after
        # the established contact criterion has already accepted this root.
        def stalled(arm, station, distance):
            return arm.point(station, distance), np.zeros(2)
        with patch.object(Approach, 'jet', stalled):
            corner = _local_corner(first, second, 2., seed=seed)
        np.testing.assert_array_equal(corner.stations, seed)
        self.assertLess(residual(corner), 1e-8)


if __name__ == '__main__':
    unittest.main()
