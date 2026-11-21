"""Merged pavement never admits crossing or touching visible centerlines."""
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.gapped_source_validity import certify_visible_centerlines
from tested_python_streets.model import Section, Traversal, Key
from tested_python_streets.paths import Path


def fixture(second, *, first_interval=None, second_interval=None,
            angle=0., shift=(0., 0.)):
    rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
    products = {};spans = []
    for role, (points, interval) in enumerate((([(-5., 0.), (5., 0.)], first_interval),
                                              (second, second_interval))):
        xy = np.asarray(points) @ rotation.T + shift
        section = Section(role+1, role+1, tuple((x, y, 0.) for x, y in xy),
                          width=1., interpolation=0)
        path = Path(section)
        products[Key('path', section.id)] = path
        first, last = (0., path.length) if interval is None else interval
        spans.append((SimpleNamespace(role=role, visits=(Traversal(section.id),),
            source_intervals=({'section': section.id, 'forward': True,
                               'first': first, 'last': last},)),))
    return products, tuple(spans)


class GappedSourceValidityTests(unittest.TestCase):
    def test_crossing_touching_and_overlapping_visible_original_curves_are_rejected(self):
        for second in (((0., -5.), (0., 5.)), ((0., 0.), (0., 5.)),
                       ((-2., 0.), (3., 0.))):
            for angle, shift in ((0., (0., 0.)), (.713, (1400., -1900.))):
                products, spans = fixture(second, angle=angle, shift=shift)
                with self.assertRaisesRegex(ValueError, 'cross or touch'):
                    certify_visible_centerlines(products, spans)

    def test_actual_owned_prefix_trim_excludes_a_junction_contact_without_editing_source(self):
        products, spans = fixture(((0., -5.), (0., 5.)), first_interval=(0., 4.5))
        before = {key: path.points.copy() for key, path in products.items()}
        report = certify_visible_centerlines(products, spans)
        self.assertGreater(report['proof_checks'], 0)
        self.assertEqual(report['source_intervals'][0][0]['last'], 4.5)
        self.assertFalse(report['mesh_samples_used'])
        self.assertFalse(report['new_source_connectivity'])
        for key, path in products.items():np.testing.assert_array_equal(path.points, before[key])
        # Certifying the complete source instead would include the owned
        # junction contact. The actual visible interval is authoritative.
        spans[0][0].source_intervals[0]['last'] = products[Key('path', 1)].length
        with self.assertRaisesRegex(ValueError, 'cross or touch'):
            certify_visible_centerlines(products, spans)

    def test_unmatched_traversals_and_extrapolated_intervals_cannot_change_the_proof(self):
        products, spans = fixture(((-5., 2.), (5., 2.)))
        spans[0][0].source_intervals[0]['forward'] = False
        with self.assertRaisesRegex(ValueError, 'original source traversal'):
            certify_visible_centerlines(products, spans)
        spans[0][0].source_intervals[0]['forward'] = True
        spans[0][0].source_intervals[0]['last'] = 11.
        with self.assertRaisesRegex(ValueError, 'visible station span'):
            certify_visible_centerlines(products, spans)


if __name__ == '__main__':unittest.main()
