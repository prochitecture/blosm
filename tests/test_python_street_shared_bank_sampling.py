"""One curved carriageway bank belongs to a dimensional corner and a Section port."""
import unittest

import numpy as np

import _python_street_package
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.model import End, Intersection, Network, Section, Settings
from tested_python_streets.ownership import PLAN
from tested_python_streets.regions import BoundaryIncidenceError, Region
from tested_python_streets.sampling import CertifiedCurbSampler
from tested_python_streets.shared_bank_sampling import SharedOffsetCurb
from tested_python_streets.shallow_bank_fillets import candidate_union


def curved_fork():
    # Three original centrelines suffice to reproduce the two independently
    # sampled representations of the 7 m carriageway at a 5 m approach.
    origin = (1238.843017578125, 586.776611328125)
    controls = (
        (3, 5., ((1249.1884765625, 660.5831298828125), (1243.8709716796875, 647.7025756835938),
                 (1237., 630.7474975585938), (1237.9185791015625, 600.056884765625), origin)),
        (1, 7., (origin, (1239.185302734375, 582.2236938476562), (1235.107177734375, 560.4823608398438))),
        (2, 5., (origin, (1253.9571533203125, 584.4859008789062), (1359.053466796875, 584.392578125))),
    )
    sections = tuple(Section(sid, sid, tuple((x, y, 0.) for x, y in xy),
                             width=width, rise=float(np.float32(.08))) for sid, width, xy in controls)
    site = Intersection(1, (End(3, 1), End(1, 0), End(2, 0)), 2., (2.,)*3, ('STYLE',)*3)
    return Network('curved-shared-bank', sections, (site,),
                   settings=Settings(cell_size=float(np.float32(.7)), exaggeration=1.))


class SharedBankTests(unittest.TestCase):
    def test_point_only_contact_without_a_fitted_arc_is_rejected(self):
        from types import SimpleNamespace
        from tested_python_streets.model import Key
        polygons = (((0., 0.), (1., 0.), (1., 1.), (0., 1.)),
                    ((1., 1.), (2., 1.), (2., 2.), (1., 2.)))
        sources = [SimpleNamespace(region=Region(Key('port', i), np.array(xy)))
                   for i, xy in enumerate(polygons)]
        with self.assertRaises(BoundaryIncidenceError):
            candidate_union(None, None, sources, ())

    def test_shallow_native_fillet_with_a_third_approach(self):
        junction = (2067.45751953125, -120.63502502441406)
        xy = (((2061.772705078125, -66.44623565673828), junction),
              (junction, (2068.468017578125, -130.3084259033203), (2078.72216796875, -186.12123107910156)),
              ((2066.24267578125, -175.29322814941406), (2062.802978515625, -150.3251953125),
               (2061.764404296875, -130.85569763183594), junction))
        sections = tuple(Section(sid, 1 if sid != 3 else 2, tuple((x, y, 0.) for x, y in path),
                                 width=width, rise=float(np.float32(.08)))
                         for sid, path, width in zip((1, 2, 3), xy, (5., 5., 2.)))
        site = Intersection(1, (End(1, 1), End(2, 0), End(3, 1)), 2., (2., 2., 1.), ('STYLE',)*3)
        network = Network('shallow-native-fillet', sections, (site,),
                          settings=Settings(cell_size=float(np.float32(.7)), exaggeration=1.))
        update = PythonMeshBackend(network.settings).prepare(network)
        layout = update.products[PLAN].by_site[site.key].physical_layout
        record = layout.shallow_bank_fillet
        self.assertLess(record['sagitta'], record['base_welding_distance'])
        self.assertEqual(record['replaced_bank_spans'], 2)
        self.assertTrue(layout.source_exterior_sampling['original_profiles_exact'])
        self.assertTrue(layout.source_exterior_sampling['original_area_operands_unchanged'])
        self.assertIn('overlap_check', update.timings)

    def test_shared_edge_matches_both_polygons_and_reaches_mesh_output(self):
        network = curved_fork()
        engine = PythonMeshBackend(network.settings)
        update = engine.prepare(network)
        group = update.products[PLAN].by_site[network.intersections[0].key]
        self.assertTrue(group.physical_layout.shared_bank_signature)
        banks = [(curve, i, piece) for source in group.area_footprints
                 for curve in source.layout.boundaries for i, piece in enumerate(getattr(curve, 'pieces', ()))
                 if isinstance(piece, SharedOffsetCurb)]
        self.assertTrue(banks)
        for curve, i, bank in banks:
            knots = (curve.lengths[i] + bank.knots * bank.length) / curve.length
            np.testing.assert_array_equal(curve(knots), bank.points)
            port, = (s for s in group.area_footprints if len(s.layout.approaches) == 1
                     and s.layout.approaches[0] is bank.approach)
            port_vertices = set(map(tuple, port.region.xy))
            self.assertTrue(all(tuple(p) in port_vertices for p in bank.points))
            with self.assertRaisesRegex(ValueError, 'finer sampling'):
                bank.sample(CertifiedCurbSampler(bank.error * .25, .7, 8192), 0., 1.)
        self.assertIn('overlap_check', update.timings)
        self.assertTrue(any(getattr(part, 'mesh', None) is not None and part.owner == network.intersections[0].key
                            and np.isfinite(part.mesh.uv_layers['UVMap']).all() for part in update.products.values()))


if __name__ == '__main__':
    unittest.main()
