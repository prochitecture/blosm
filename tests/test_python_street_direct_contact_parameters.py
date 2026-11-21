"""Direct original stations survive reversed native aliases without an epsilon."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route,Piece
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.compound_contact_routes import outward_station


def approach():
    section=Section(1,1,((0.,0.,0.),(330.5250496738632,0.,0.)),width=5.,interpolation=0)
    return Approach(Route((Piece(Path(section),False),)),1,Crown(.08,0.))


class DirectContactParameterTests(unittest.TestCase):
    def test_own_reversed_contact_preserves_exact_original_station_and_native_profile_cut(self):
        arm=approach();station=32.344042464820596
        cuts={key:{value} for key,value in arm.path.trims(station).items()}
        rounded=outward_station(arm,cuts)
        self.assertNotEqual(station,rounded)
        exact=outward_station(arm,cuts,direct_contacts=((arm,station),))
        self.assertEqual(exact,station)
        self.assertEqual(arm.path.trims(exact+.25),arm.path.trims(station+.25))
        np.testing.assert_array_equal(arm.point(exact+.25,0.),arm.point(station+.25,0.))

    def test_another_actual_native_contact_still_extends_the_required_finite_bound(self):
        arm=approach();station=32.344042464820596;far=40.
        cuts={key:{value} for key,value in arm.path.trims(station).items()}
        for key,value in arm.path.trims(far).items():cuts.setdefault(key,set()).add(value)
        self.assertEqual(outward_station(arm,cuts,direct_contacts=((arm,station),)),far)
        other=approach()
        self.assertEqual(outward_station(arm,cuts,direct_contacts=((other,station),)),outward_station(arm,cuts))

    def test_nonfinite_or_extrapolated_direct_coordinates_cannot_replace_a_native_bound(self):
        arm=approach()
        for station in (float('nan'),float('inf'),-1.,arm.length+1.):
            with self.assertRaises(ValueError):outward_station(arm,{},direct_contacts=((arm,station),))


if __name__=='__main__':unittest.main()
