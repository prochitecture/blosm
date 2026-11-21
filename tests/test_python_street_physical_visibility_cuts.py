"""A contained source owner needs no invented ownership-switch seam."""
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.compounds import CompoundGroup, RimSegment, SampledFootprint, SharedClaim
from tested_python_streets.partition import BoundaryOwnershipCuts, CompoundPartitioner
from tested_python_streets.paths import Path
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.regions import Region, RegionUnion
from tested_python_streets.updates import BuildContext


def contained_owner(*, outside=False, partial=False):
    owners = (model.Key('intersection', 1), model.Key('intersection', 2))
    polygons = (((0., 0.), (10., 0.), (10., 10.), (0., 10.)),
                ((4., 4.), (6., 4.), (6., 6.), (4., 6.)))
    sources = []
    for owner, xy in zip(owners, polygons):
        labels = [RimSegment('curb', i, 0., 1.) for i in range(4)]
        if partial and owner == owners[0]:
            labels[0] = RimSegment('mouth', 0, .8, 0.)
        sources.append(SampledFootprint(Region(owner, np.asarray(xy)), tuple(labels),
                                       SimpleNamespace(approaches=())))
    points = ((-3., 5., 0.), (6., 5., 0.)) if outside else ((4., 5., 0.), (6., 5., 0.))
    section = model.Section(1, 1, points, width=2., interpolation=0)
    network = model.Network('contained-source-owner', (section,), ())
    context = BuildContext(network, {model.Key('path', 1): Path(section)}, PatternLibrary())
    group = CompoundGroup(owners, (SharedClaim(1, *owners, 1.5, .5),), tuple(sources),
                          RegionUnion().build(f.region for f in sources))
    return group, context


class PhysicalVisibilityCutsTests(unittest.TestCase):
    def test_contained_owner_uses_qualified_physical_rim_without_changing_claims(self):
        group, context = contained_owner()
        policy = BoundaryOwnershipCuts()
        partition = CompoundPartitioner(policy=policy).build(group, context)
        self.assertEqual(partition.group.claims, group.claims)
        self.assertEqual(partition.group.members, group.members)
        self.assertEqual(len(partition.blocks), 1)
        self.assertEqual(partition.cuts, ())
        self.assertEqual(partition.area_error, 0.)
        np.testing.assert_array_equal(partition.xy[list(partition.boundary_loops[0])], group.loops[0].xy)
        self.assertEqual(policy.visibility_source_proof['total_cells'], 1)

    def test_a_clipped_original_strip_cannot_use_the_cut_fallback(self):
        group, context = contained_owner(outside=True)
        with self.assertRaisesRegex(ValueError, 'complete source strip'):
            CompoundPartitioner().build(group, context)

    def test_a_partial_external_profile_cannot_use_the_cut_fallback(self):
        group, context = contained_owner(partial=True)
        with self.assertRaisesRegex(ValueError, 'no unique pair'):
            CompoundPartitioner().build(group, context)


if __name__ == '__main__':
    unittest.main()
