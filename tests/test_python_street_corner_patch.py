"""Finite source fillets close through strips rather than radial arc spokes."""
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.corner_patch import corner_patch
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.intersections import Approach
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.regions import Region, RegionUnion, signed_area
from tested_python_streets.routes import Piece, Route


def approaches(case, angle=0., offset=(0., 0.), routed=False):
    # Original local metre-space controls at Beijing source site193. The
    # opposite through arms need negative auxiliary circular contacts.
    controls = {
        87: (10.8, ((865.008544921875, -513.74609375), (870.5442504882812, -528.9962158203125))),
        90: (10.8, ((870.5442504882812, -528.9962158203125), (873.9628295898438, -551.0148315429688),
                    (877.1314697265625, -701.8857421875))),
        97: (2.7, ((862.883056640625, -522.5739135742188), (870.5442504882812, -528.9962158203125)))}
    if case == 'positive':
        controls = {1: (6., ((0., 0.), (40., 0.))), 2: (6., ((0., 0.), (0., 40.)))}
        ends = ((1, 0), (2, 0))
    else:
        ends = ((90, 0), (87, 1)) if case == 'negative' else ((97, 1), (90, 0))
    c, s = math.cos(angle), math.sin(angle)
    rotation = np.array(((c, -s), (s, c)))
    result = []
    for identity, end in ends:
        width, points = controls[identity]
        xy = np.array(points) @ rotation.T + offset
        source = model.Section(identity, identity, tuple((float(x), float(y), 0.) for x, y in xy), width=width)
        path = Path(source)
        if routed:
            path = Route((Piece(path, end == 0),))
        result.append(Approach(path, end, Crown(.08, 0., -1)))
    return tuple(result)


def strip(approach, identity):
    at = np.linspace(0., min(20., approach.length), 201)
    xy = np.vstack((approach.point(at, .5*approach.width), approach.point(at, -.5*approach.width)[::-1]))
    if signed_area(xy) < 0.:
        xy = xy[::-1]
    return Region(model.Key('strip', identity), xy)


class CornerPatchTests(unittest.TestCase):
    def test_actual_positive_negative_and_mixed_contacts_have_simple_connected_source_coverage(self):
        for case in ('positive', 'negative', 'mixed'):
            for routed in (False, True):
                for angle, offset in ((0., (0., 0.)), (.713, (1200., -1800.))):
                    with self.subTest(case=case, routed=routed, angle=angle):
                        first, second = approaches(case, angle, offset, routed)
                        corner = CornerRadiusPolicy.solve(first, second, 4., 'STRICT')
                        if case == 'negative':
                            self.assertTrue(all(value < 0. for value in corner.stations))
                        elif case == 'mixed':
                            self.assertLess(min(corner.stations), 0.)
                            self.assertGreater(max(corner.stations), 0.)
                        source = corner_patch(model.Key('patch', 1), corner, .5, .00025)
                        self.assertGreater(source.region.xy.shape[0], 3)
                        self.assertEqual(len(source.segments), len(source.region.xy))
                        RegionUnion().build((source.region,))
                        loops = RegionUnion().build((source.region, strip(first, 1), strip(second, 2)))
                        self.assertEqual(len(loops), 1)
                        retained = [span for span in loops[0].spans if span.owner == source.region.owner]
                        self.assertTrue(retained)
                        self.assertTrue(all(source.segments[span.segment].kind == 'curb' for span in retained))
                        for i, label in enumerate(source.segments):
                            curve = source.layout.boundaries[label.index]
                            np.testing.assert_allclose(source.region.xy[i], curve(label.first), atol=1e-10, rtol=0.)
                            np.testing.assert_allclose(source.region.xy[(i+1) % len(source.region.xy)],
                                                       curve(label.last), atol=1e-10, rtol=0.)

    def test_invalid_requests_and_sampling_budget_fail_explicitly(self):
        first, second = approaches('positive')
        corner = CornerRadiusPolicy.solve(first, second, 4., 'STRICT')
        with self.assertRaisesRegex(ValueError, 'budget'):
            corner_patch(model.Key('patch', 1), corner, .1, .00001, max_samples=4)
        with self.assertRaisesRegex(ValueError, 'sampling request'):
            corner_patch(model.Key('patch', 1), corner, 0., .00025)


if __name__ == '__main__':
    unittest.main()
