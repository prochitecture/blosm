"""Coordinate alternatives align actual finite profiles without moving sources."""
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_terminal_ports import fixture
from tested_python_streets.finite_terminal_entrances import FinitePairedEntrance
from tested_python_streets.bundle_mouth_frames import candidate_frames
from tested_python_streets.model import Key


class MouthFrameTests(unittest.TestCase):
    def test_declared_mouth_normals_align_both_whole_profiles_and_certify_every_original_offset(self):
        context, pair, sites = fixture(3, reverse=True, reflection=True, angle=.713,
                                       translation=(1800., -1200.))
        first = FinitePairedEntrance(context, pair.run, 0, source_tolerance=.00025)
        last = FinitePairedEntrance(context, pair.run, 1, source_tolerance=.00025)
        a, b = first.minimum_station+3., last.minimum_station+3.
        layout = SimpleNamespace(approaches=(first,), outward_stations=(a,))
        side_layout = SimpleNamespace(approaches=(last,), outward_stations=(b,))
        specification = SimpleNamespace(group=SimpleNamespace(members=(sites[0].key,), physical_layout=layout),
            paired_sources=(SimpleNamespace(owner=Key('intersection', 900), entrances=(last,),
                            group=SimpleNamespace(physical_layout=side_layout)),))
        before = context.network.to_dict()
        choices = list(candidate_frames(context, context.network.bundles[0], first, specification))
        self.assertEqual(choices[0][1], {'family': 'established'})
        self.assertGreaterEqual(len(choices), 2)
        for frame, proof in choices[1:]:
            for member in pair.run.members:
                for visit in member.visits:
                    self.assertGreater(frame.certify(context[Key('path', visit.section_id)], forward=visit.forward), 0.)
            for approach, station in ((first, a), (last, b)):
                xy = approach.point(station, approach.width*np.linspace(-.5, .5, 17))
                self.assertLess(float(np.ptp(frame.local(xy)[:, 0])), 1e-8)
            self.assertFalse(proof['physical_geometry_changed'])
            self.assertEqual(set(proof['original_sections']), set(context.network.bundles[0].sections))
        self.assertEqual(context.network.to_dict(), before)

    def test_unrelated_approach_cannot_supply_a_normal_line_and_no_sides_keep_original_candidate(self):
        context, pair, sites = fixture()
        first = FinitePairedEntrance(context, pair.run, 0, source_tolerance=.00025)
        layout = SimpleNamespace(approaches=(first,), outward_stations=(first.minimum_station+3.,))
        specification = SimpleNamespace(group=SimpleNamespace(members=(sites[0].key,), physical_layout=layout), paired_sources=())
        self.assertEqual(len(list(candidate_frames(context, context.network.bundles[0], first, specification))), 1)
        wrong = SimpleNamespace(frame_chart=first.frame_chart)
        specification.group.physical_layout = SimpleNamespace(approaches=(wrong,), outward_stations=(3.,))
        with self.assertRaisesRegex(ValueError, 'actual original prepared approach'):
            list(candidate_frames(context, context.network.bundles[0], first, specification))


if __name__ == '__main__': unittest.main()
