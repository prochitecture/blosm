"""Explicitly qualified source preparation can use the complete Bundle plan."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_merged_bundles import source_for
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_merged_entrances import PairedAngularEntrance,preparation_key
from tested_python_streets.bundle_merged_plan import MergedBundlePlan
from tested_python_streets.paired_source_sampling import verify_shared_source_phase
from tested_python_streets.updates import BuildContext


class HydratedBundlePlanTests(unittest.TestCase):
    def test_complete_run_scope_query_ownership_and_release_in_both_modes(self):
        for mode,exaggeration in ((3,1.),(0,0.)):
            with self.subTest(mode=mode,exaggeration=exaggeration):
                network=source_for(mode,exaggeration)
                engine=PythonMeshBackend(network.settings)
                prepared=engine.prepare(network)
                context=BuildContext(network,prepared.products,engine.patterns)
                bundle=network.bundles[0]
                specification=context[preparation_key(bundle)]
                phase=verify_shared_source_phase(tuple(a for a in specification.layout.approaches
                                                       if isinstance(a,PairedAngularEntrance)))
                plan=MergedBundlePlan(context,bundle,replace(specification,source_phase_proof=phase))
                metadata=plan._pavement.metadata
                self.assertTrue(metadata['complete_selected_source_coverage'])
                self.assertFalse(metadata['junction_geometry_republished'])
                self.assertEqual(len(metadata['runs']),2)
                self.assertTrue(all(run['coverage']['complete_selected_scope'] for run in metadata['runs']))
                self.assertEqual(set(plan.source_ids),{17,29,117,129})
                self.assertGreater(metadata['overlap']['quads'],len(plan._pavement.mesh.quads))
                xy=plan._pavement.mesh.xyz[::13,:2].copy()
                before=plan.height(xy)
                physical=plan.query.evaluate(xy)
                self.assertGreater(float(physical.max()),.07)
                np.testing.assert_allclose(before,physical*exaggeration,atol=1e-12,rtol=0.)
                plan.release_geometry()
                self.assertIsNone(plan._pavement.mesh)
                if plan._separator is not None:self.assertIsNone(plan._separator.mesh)
                np.testing.assert_array_equal(before,plan.height(xy))


if __name__=='__main__':unittest.main()
