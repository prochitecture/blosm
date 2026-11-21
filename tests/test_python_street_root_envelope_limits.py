"""Explicit root-projection bounds leave established envelope calls unchanged."""
import unittest
import numpy as np
import _python_street_package
from test_python_street_root_envelopes import actual274
from tested_python_streets.model import Settings,Key
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.compound_root_envelopes import _search


class RootEnvelopeLimitTests(unittest.TestCase):
    def test_default_and_explicit_none_keep_actual_source_output_byte_exact(self):
        a,b=actual274();settings=Settings(curve_tolerance=.00025)
        corner=CornerRadiusPolicy.solve(a,b,10.,'STYLE')
        stations=[IntersectionProvider.mouth_station(p,max(0.,s)+min(.25,.05*p.width),settings)
                  for p,s in zip((a,b),corner.stations)]
        default=_search(Key('unchanged',1),corner,settings,stations)
        explicit=_search(Key('unchanged',1),corner,settings,stations,outward_limit=None)
        self.assertEqual(default.report,explicit.report)
        self.assertEqual(default.stations,explicit.stations)
        np.testing.assert_array_equal(default.source.region.xy,explicit.source.region.xy)
        np.testing.assert_array_equal(default.curve(np.linspace(0.,1.,101)),explicit.curve(np.linspace(0.,1.,101)))
        self.assertNotIn('outward_limit_policy',default.report)

    def test_explicit_limit_requires_original_root_points_and_source_width_bound(self):
        a,b=actual274();settings=Settings(curve_tolerance=.00025)
        corner=CornerRadiusPolicy.solve(a,b,10.,'STYLE');stations=(.25,.25)
        with self.assertRaisesRegex(ValueError,'actual root projections'):
            _search(Key('invalid',1),corner,settings,stations,outward_limit=.1)
        points=(a.point(0.,-.5*a.width),b.point(0.,.5*b.width))
        # A declared seven-metre support permits at most half its width
        # difference from the 3.5-metre donor, rather than an arbitrary bell.
        for limit in (-.01,float('nan'),float('inf'),1.750001):
            with self.subTest(limit=limit),self.assertRaisesRegex(ValueError,'source-width bound'):
                _search(Key('invalid',1),corner,settings,stations,root_points=points,
                        support_width=7.,outward_limit=limit)


if __name__=='__main__':unittest.main()
