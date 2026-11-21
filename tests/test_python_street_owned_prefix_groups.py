"""Source unions certify owned prefixes using their actual loops, including holes."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_owned_prefixes import fixture
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.model import Key
from tested_python_streets.updates import BuildContext
from tested_python_streets.compounds import CompoundGroup,sample_footprint
from tested_python_streets.regions import RegionUnion,BoundaryLoop
from tested_python_streets.owned_prefixes import certify


class OwnedPrefixGroupTests(unittest.TestCase):
    @staticmethod
    def source():
        n=fixture();engine=PythonMeshBackend(n.settings);p=engine.prepare(n).products
        site=n.intersections[0];layout=p[site.layout_key]
        footprint=sample_footprint(site.key,layout)
        group=CompoundGroup((site.key,),(),(footprint,),RegionUnion().build((footprint.region,)),layout)
        return BuildContext(n,p,engine.patterns),site,layout,group,p[Key('path',1)].owned_prefix

    def test_actual_group_loops_replace_the_single_boundary_sampling_protocol(self):
        context,site,layout,group,claim=self.source()
        # Source-cycle layouts have no single concatenated boundary protocol.
        layout.boundaries=()
        certify(context,site,layout,(claim,),group=group)
        proof=layout.owned_prefix_width_proof
        self.assertEqual([row['section'] for row in proof['sections']],[1])
        self.assertLessEqual(proof['maximum_error'],.002)

    def test_a_foreign_group_or_claim_cannot_authorize_the_prefix(self):
        context,site,layout,group,claim=self.source()
        with self.assertRaisesRegex(ValueError,'actual prepared source group'):
            certify(context,site,replace(layout), (claim,),group=group)
        with self.assertRaisesRegex(ValueError,'real original junction claim'):
            certify(context,site,layout,(replace(claim,junction=Key('intersection',999)),),group=group)

    def test_a_real_hole_inside_the_original_prefix_still_rejects(self):
        context,site,layout,group,claim=self.source()
        hole=BoundaryLoop(np.array(((.10,.32),(.10,.38),(.16,.38),(.16,.32))),())
        group=replace(group,loops=(*group.loops,hole))
        with self.assertRaisesRegex(ValueError,'hole'):
            certify(context,site,layout,(claim,),group=group)


if __name__=='__main__':unittest.main()
