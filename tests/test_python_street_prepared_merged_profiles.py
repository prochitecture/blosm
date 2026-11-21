"""Physical body interfaces use both real and owner-scoped derived samples."""
import unittest
import numpy as np
import _python_street_package
from test_python_street_one_sided_paired_junctions import fixture
from tested_python_streets.one_sided_paired_junctions import OneSidedPairedSpecification
from tested_python_streets.ownership import CompoundProduct
from tested_python_streets.prepared_merged_profiles import PreparedMergedMouth


class PreparedMergedProfileTests(unittest.TestCase):
    def test_both_actual_one_sided_mouths_include_their_complete_shared_profile(self):
        context, site, layout, visits, opposite = fixture()
        specification = OneSidedPairedSpecification.build(context, site, layout, visits, opposite, 0)
        product = CompoundProduct.build(context, site.key, specification.group)
        before = context.network.to_dict()
        count = 0
        for index, approach in enumerate(product.group.physical_layout.approaches):
            if not hasattr(approach, 'source_mouths'):continue
            mouth = PreparedMergedMouth.from_product(product, index)
            self.assertFalse(mouth.xyz.flags.writeable)
            first, last = mouth.split_body_cut(mouth.xyz)
            np.testing.assert_array_equal(first[-1], last[-1])
            self.assertGreaterEqual(len(first), 3); self.assertGreaterEqual(len(last), 3)
            for point, (loop, vertex) in zip(mouth.xyz, mouth.rim_references):
                np.testing.assert_array_equal(point, product.physical_rims[loop][vertex])
            for role, points in enumerate((first, last)):
                local = approach.frame_chart.local(points[:, :2])
                self.assertTrue(np.all(np.diff(local[:, 1]) * (1. if role == 0 else -1.) < 0.))
            self.assertEqual(sum(a.original_incidence for a in approach.source_mouths(mouth.station)), 1)
            with self.assertRaisesRegex(ValueError, 'whole-profile partition'):
                mouth.split_body_cut(mouth.xyz[:-1])
            count += 1
        self.assertEqual(count, 2)
        self.assertEqual(context.network.to_dict(), before)


if __name__ == '__main__':unittest.main()
