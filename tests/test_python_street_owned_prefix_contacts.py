"""Displaced contact closures follow the real consumed dimensional prefix."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_owned_prefixes import fixture
from test_python_street_dimensional_bundles import reverse_sections
from tested_python_streets.model import Key,End
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import UpdateCoordinator,BuildContext
from tested_python_streets.routes import Route,RouteTopology
from tested_python_streets.intersections import Approach
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.owned_prefixes import approaches,contact_walks
from tested_python_streets.compound_displaced_contacts import patch


class OwnedPrefixContactTests(unittest.TestCase):
    @staticmethod
    def source(reverse=False):
        n=reverse_sections(fixture(),{1,2,3,4}) if reverse else fixture()
        patterns=PatternLibrary()
        p=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(n).products
        ctx=BuildContext(n,p,patterns);site=n.intersections[0]
        original=[Approach(Route.build(ctx,RouteTopology(n).trace(end)),end.end,p[Key('profile',end.section_id)])
                  for end in site.ends]
        resolved,claims=approaches(ctx,site,original)
        by_id={a.source.id:a for a in resolved}
        corner=CornerRadiusPolicy.solve(by_id[4],by_id[1],4.,'STYLE')
        return ctx,site,corner

    def test_actual_frames_profiles_radius_and_original_prefix_survive_reversal(self):
        for reverse in (False,True):
            ctx,site,corner=self.source(reverse)
            before=ctx.network.to_dict()
            result=contact_walks(ctx,site.key,corner)
            self.assertIsNot(result,corner)
            np.testing.assert_array_equal(result.center,corner.center)
            self.assertEqual(result.stations,corner.stations)
            self.assertEqual(result.radius,4.)
            for old,new in zip((corner.first,corner.second),(result.first,result.second)):
                self.assertIs(old,new)
            source=patch(Key('contact',1),result,.3,.00025)
            for label in source.segments:
                if label.kind=='curb':
                    curve=source.layout.boundaries[label.index]
                    # The physical curve keeps the original sampled circle.
                    points=curve.points
                    np.testing.assert_allclose(np.linalg.norm(points-result.center,axis=1),4.,atol=1e-8)
            walk=result.source_contact_walks[1]
            self.assertEqual(walk.prefix[0].section.id,1)
            self.assertEqual(ctx.network.to_dict(),before)
            self.assertFalse(hasattr(corner,'source_contact_walks'))

    def test_claims_cannot_be_attached_to_another_source_junction(self):
        ctx,site,corner=self.source()
        with self.assertRaisesRegex(ValueError,'explicit owned-prefix incidence'):
            contact_walks(ctx,Key('intersection',999),corner)
        path=corner.second.path
        path.claim=replace(path.claim,external=End(3,0))
        with self.assertRaisesRegex(ValueError,'explicit owned-prefix incidence'):
            contact_walks(ctx,site.key,corner)


if __name__=='__main__':unittest.main()
