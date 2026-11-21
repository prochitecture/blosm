"""A divided approach meets two full-profile roads at a shallow angle."""
from dataclasses import replace
import unittest

import numpy as np
import _python_street_package
from test_python_street_compounds import footprint_context
from tested_python_streets import model
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.ordinary_mouth_clearance import prepare
from tested_python_streets.mixed_bank_continuation import _mixed_side


def divided_approach():
    root = (-1611.645263671875, 224.0830841064453)
    half = (root, (-1623.4234619140625, 308.0433044433594), (-1630.0091552734375, 377.8534851074219))
    full = ((-1382.116943359375, -311.95318603515625), (-1402.9962158203125, -276.2594909667969),
            (-1479.528076171875, -145.4191131591797), (-1521.5689697265625, -61.62992858886719),
            (-1550.6126708984375, 12.07073974609375), (-1577.9051513671875, 87.51881408691406),
            (-1600.5189208984375, 174.71002197265625), root)
    link = ((-1575.1964111328125, 125.09959411621094), (-1595.9095458984375, 192.263916015625),
            (-1601.9722900390625, 210.0877227783203), root)
    partner = [(x+20., y) for x, y in half]
    sections = tuple(model.Section(i, i, tuple((x, y, 0.) for x, y in xy), width=w, offset=delta,
                                   rise=float(np.float32(.08)))
                     for i, xy, w, delta in ((1, half, 10.5, 1.75), (2, full, 7., 0.),
                                              (3, link, 3.5, 0.), (4, partner, 10.5, 0.)))
    site = model.Intersection(1, (model.End(1, 0), model.End(2, 1), model.End(3, 1)),
                              10., (10.,)*3, ('STYLE',)*3)
    bundle = model.Bundle(1, tuple(model.BundleMember((model.Traversal(i, True),)) for i in (4, 1)))
    return model.Network('mixed-grade', sections, (site,), (bundle,),
                         settings=model.Settings(cell_size=float(np.float32(.7)), exaggeration=1.))


class MixedBankGradeTests(unittest.TestCase):
    def test_style_fit_uses_original_profiles_and_certifies_width_radius_and_grade(self):
        network = divided_approach()
        ctx = footprint_context(network)
        site, = network.intersections
        raw = IntersectionProvider.footprint(ctx, site)
        effective, _ = prepare(ctx, (site,), {site.key: raw}, targets={site.key}, source_singletons={})
        layout = effective[site.key]
        proof = layout.ordinary_complete_mouth_clearance
        record = proof['mixed_bank_grade']
        self.assertLess(record['grade_bound'], record['slope_limit'])
        self.assertLessEqual(record['curvature_bound'], 1./record['original_radius'])
        self.assertTrue(proof['complete_source_width']['sections'])
        original = raw.corners[record['corner']]
        self.assertIsNone(_mixed_side(replace(original, radius_mode='STRICT')))
        for before, after in zip(raw.approaches, layout.approaches):
            self.assertIs(before, after)


if __name__ == '__main__':
    unittest.main()
