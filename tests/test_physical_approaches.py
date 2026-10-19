"""External approach ownership must not invent a free end or omit a transition."""
from copy import deepcopy
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest

from test_physical_junctions import contract as k

spec = importlib.util.spec_from_file_location('tested_physical_approaches',
    Path(__file__).resolve().parents[1] / 'way/street_network/physical_approaches.py')
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def fixture(*, leaving=True, dimensions=False, far=None):
    sections = (k.SectionInput(10, ((0., 0.), (10., 0.)), 6.),
                k.SectionInput(11, ((10., 0.), (20., 1.)), 8. if dimensions else 6.))
    street = k.StreetInput(7, sections)
    connector = k.ConnectorInput(3, (k.MemberInput(0, leaving),), 0, radius=3.)
    junction = k.PreparedJunction(1, (0., 0.), (street,), (connector,), (k.ProfileInput(.08),))
    snapshot = NS(streets={7: dict(status='SUPPORTED', endpoints=(1, far) if leaving else (far, 1))})
    return snapshot, (junction,)


class PhysicalApproachTests(unittest.TestCase):
    def test_true_free_and_remote_are_distinct_in_both_directions(self):
        for leaving in (False, True):
            for far, role in ((None, module.ApproachEndRole.FREE), (25, module.ApproachEndRole.REMOTE_JUNCTION)):
                snapshot, junctions = fixture(leaving=leaving, far=far)
                original = deepcopy((snapshot, junctions))
                row, = module.PreparedPhysicalApproaches().prepare(snapshot, junctions)
                self.assertEqual((row.street_id, row.leaving, row.section_ids), (7, leaving, (10, 11)))
                self.assertEqual((row.end_role, row.owner_id), (role, -1 if far is None else far))
                self.assertEqual((snapshot, junctions), original)

    def test_dimension_handoff_precedes_even_an_explicit_free_end(self):
        for leaving, selected, following in ((True, (10,), 11), (False, (11,), 10)):
            snapshot, junctions = fixture(leaving=leaving, dimensions=True)
            row, = module.PreparedPhysicalApproaches().prepare(snapshot, junctions)
            self.assertEqual((row.end_role, row.section_ids, row.next_section_id),
                             (module.ApproachEndRole.DIMENSION_CHANGE, selected, following))

    def test_missing_roles_never_default_to_free(self):
        for value in (None, (), (1,), (1, -1), (1, '25')):
            snapshot, junctions = fixture()
            snapshot.streets[7]['endpoints'] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                module.PreparedPhysicalApproaches().prepare(snapshot, junctions)

    def test_both_ends_owned_are_internal_and_duplicates_reject(self):
        snapshot, junctions = fixture()
        first, = junctions
        second = replace(first, source_id=2, connectors=(replace(first.connectors[0],
            members=(k.MemberInput(0, False),)),))
        self.assertEqual(module.PreparedPhysicalApproaches().prepare(snapshot, (first, second)), ())
        with self.assertRaisesRegex(ValueError, 'multiple physical'):
            module.PreparedPhysicalApproaches().prepare(snapshot, (first, first))


if __name__ == '__main__':
    unittest.main()
