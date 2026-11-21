"""Mixed full and Bundle-half junctions with exposed carriageway ends."""
import math
import unittest

import numpy as np
import _python_street_package
from test_python_street_finite_scope_compounds import OwnershipOnly
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext, UpdateCoordinator
from tested_python_streets.model import Bundle, BundleMember, End, Intersection, Network, Section, Settings, Traversal
from tested_python_streets.ownership import PLAN, footprint_key, bundle_signature

CAP_CHAIN = (
 ((11.0802001953125,16.2496337890625),(9.453857421875,13.5115966796875),(0.,0.)),
 ((0.,0.),(-3.410888671875,-15.149658203125),(-2.3485107421875,-26.9609375),
  (.5966796875,-36.49066162109375),(6.5594482421875,-47.12335205078125),
  (13.137451171875,-54.92864990234375),(20.4083251953125,-60.80828857421875),
  (25.8988037109375,-63.8487548828125)),
 ((0.,0.),(-11.3133544921875,-15.28106689453125),(-31.0225830078125,-42.8160400390625)))



ENDPOINT = (
 ((31.0225830078125,42.8160400390625),(19.709228515625,27.53497314453125),(0.,0.)),
 ((0.,0.),(-21.688720703125,-13.43017578125),(-42.403564453125,-20.7156982421875),(-62.8623046875,-29.33709716796875)),
 ((0.,0.),(-34.6240234375,-42.38067626953125),(-58.8819580078125,-63.5689697265625)))




def network_for(controls, widths, member, radii, reverse):
    angle = .713 if reverse else 0.
    rotation = np.array(((math.cos(angle), -math.sin(angle)), (math.sin(angle), math.cos(angle))))
    sections, ends = [], []
    for sid, (points, width, end) in enumerate(zip(controls, widths, (1, 0, 0)), 1):
        xy = np.asarray(points) @ rotation.T
        if reverse:
            xy = xy[::-1] + (1200., -1800.)
            end = 1 - end
        sections.append(Section(sid, sid, tuple((x, y, 0.) for x, y in xy),
                                width=width, interpolation=3, scope_ends=(1-end,)))
        ends.append(End(sid, end))
    member_xy = np.array(sections[member-1].points)[:, :2]
    travel = member_xy[-1] - member_xy[0] if reverse else member_xy[0] - member_xy[-1]
    left = np.array((-travel[1], travel[0])) / np.linalg.norm(travel)
    xy = member_xy + 12. * left
    sections.append(Section(4, 4, tuple((x, y, 0.) for x, y in xy),
                            width=3.5, interpolation=3, scope_ends=(0, 1)))
    bundle = Bundle(1, tuple(BundleMember((Traversal(sid, reverse),)) for sid in (4, member)))
    site = Intersection(1, tuple(ends), 10., radii, ('STYLE',)*3)
    return Network('native-root-joins', tuple(sections), (site,), (bundle,),
                   settings=Settings(cell_size=.7, curve_tolerance=.005, exaggeration=1.))


class NativeRootJoinTests(unittest.TestCase):
    def check_join(self, controls, widths, member, radii, records):
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                network = network_for(controls, widths, member, radii, reverse)
                patterns = PatternLibrary()
                coordinator = UpdateCoordinator((PathProvider(), ProfileProvider(), OwnershipOnly()), patterns)
                update = coordinator.prepare(network)
                site = network.intersections[0].key
                group = update.products[PLAN].by_site[site]
                original = update.products[footprint_key(site)]
                self.assertEqual(group.physical_layout.mouths, original.mouths)
                self.assertEqual(group.physical_layout.outward_stations, original.outward_stations)
                record, = getattr(group.physical_layout, records)
                self.assertAlmostEqual(record.report['full_curve_exposure'], 1., places=8)
                context = BuildContext(network, update.products, patterns, update.task_fingerprints,
                                       update.source_fingerprints)
                record.reference.validate(context)
                self.assertTrue(bundle_signature(context, group))
                coordinator.commit(update)
                self.assertFalse(coordinator.prepare(network).rebuilt)

    def test_two_exposed_carriageway_ends(self):
        self.check_join(CAP_CHAIN, (14., 14., 3.5), 3, (), 'source_cycle_owned_native_cap_chains')

    def test_full_crown_endpoint_joining_half_crown(self):
        self.check_join(ENDPOINT, (3.5, 2.5, 3.5), 1, (10., 3., 10.), 'source_cycle_native_zero_bank_roots')


if __name__ == '__main__':
    unittest.main()
