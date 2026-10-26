"""Finite Bundle mouths retain source profiles and shared whole-edge seams."""
from dataclasses import replace
import gc
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import skew_bundle_fixture, height_fixture_networks
from test_python_street_attached_unions import meshes
from test_python_street_bundles import boundary_edges
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_caps import FiniteAxialPlan
from tested_python_streets.providers import Part


class PythonSkewCapTests(unittest.TestCase):
    def equal(self,first,second):
        self.assertEqual(first.keys(),second.keys())
        for key in first:
            np.testing.assert_array_equal(first[key][0],second[key][0])
            np.testing.assert_array_equal(first[key][1],second[key][1])

    def test_bounded_cap_selection_preserves_profiles_curbs_and_source_controls(self):
        for move in (.2,-.2,1.,-1.,3.):
            with self.subTest(move=move):
                network=skew_bundle_fixture(model,move=move)
                before=network.to_dict()
                result=PythonMeshBackend(network.settings).prepare(network)
                plan=result.products[model.Key('bundle_plan',7)]
                self.assertIsInstance(plan,FiniteAxialPlan)
                self.assertEqual(network.to_dict(),before)
                part=result.products[model.Key('bundle_pavement',7)]
                self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'],8.)
                self.assertLessEqual(part.mesh.slope_statistics(10.)['maximum'],.1+1e-8)
                self.assertLessEqual(part.mesh.xyz[:,2].max(),1.2+1e-8)
                self.assertGreaterEqual(part.mesh.xyz[:,2].min(),-1e-8)
                self.assertLessEqual(part.metadata['maximum_entrance_longitudinal_slope_bound'], .03+1e-8)
                self.assertEqual(len(part.metadata['cap_entrance_longitudinal_slope_bounds']), 4)
                junctions=np.concatenate([boundary_edges(p.mesh) for key,p in result.products.items()
                    if key.kind=='junction_mesh' and getattr(p,'mesh',None) is not None])
                matches=0
                for edge in boundary_edges(part.mesh):
                    error=np.minimum(np.max(abs(junctions-edge),axis=(1,2)),np.max(abs(junctions[:,::-1]-edge),axis=(1,2)))
                    matches+=float(error.min())<1e-8
                self.assertEqual(matches,64)
                self.assertEqual(result.products[model.Key('separator',7)].metadata['islands'],2)
                self.assertTrue(all(not p.metadata['fitted_corner_radii'] for key,p in result.products.items()
                    if key.kind=='junction_mesh' and isinstance(p,Part)))

    def test_free_skew_caps_and_consumed_section_chains_use_the_same_solver(self):
        for split,free in ((False,True),(True,False)):
            with self.subTest(split=split,free=free):
                network=skew_bundle_fixture(model,split=split,free=free)
                result=PythonMeshBackend(network.settings).prepare(network)
                self.assertIsInstance(result.products[model.Key('bundle_plan',7)],FiniteAxialPlan)
                self.assertEqual(result.products[model.Key('separator',7)].metadata['attached_caps'],[not free]*2)
                self.assertEqual(result.products[model.Key('separator',7)].metadata['islands'],2)

    def test_finer_resolutions_exact_smooth_and_translated_zero_height(self):
        for mode in (0,3):
            for resolution in (1,2,3,4):
                with self.subTest(mode=mode,resolution=resolution):
                    network=skew_bundle_fixture(model)
                    network=replace(network,sections=tuple(replace(s,interpolation=mode,datum=37.) for s in network.sections),
                        settings=replace(network.settings,resolution=resolution,exaggeration=0. if resolution==2 else 10.))
                    result=PythonMeshBackend(network.settings).prepare(network)
                    self.assertIsInstance(result.products[model.Key('bundle_plan',7)],FiniteAxialPlan)
                    self.assertLessEqual(result.products[model.Key('bundle_pavement',7)].metadata[
                        'maximum_entrance_longitudinal_slope_bound'], .03+1e-8)
                    for kind in ('bundle_pavement','separator'):
                        mesh=result.products[model.Key(kind,7)].mesh
                        self.assertLessEqual(mesh.validate()['maximum_edge_ratio'],8.)
                        self.assertGreaterEqual(mesh.xyz[:,2].min(),37.-1e-8)
                        self.assertLessEqual(mesh.xyz[:,2].max(),37.+.12*network.settings.exaggeration+1e-8)
                    del mesh,result
                    gc.collect()

    def test_staggered_free_ends_unequal_widths_and_offsets(self):
        network=skew_bundle_fixture(model,free=True)
        lower=network.section_index[2]
        points=((-38.,-8.,0.),*lower.points[1:])
        network=replace(network,sections=tuple(replace(s,width=8.,offset=.5,points=points) if s.id==2 else s
                                               for s in network.sections))
        result=PythonMeshBackend(network.settings).prepare(network)
        self.assertIsInstance(result.products[model.Key('bundle_plan',7)],FiniteAxialPlan)
        self.assertEqual(result.products[model.Key('separator',7)].metadata['attached_caps'],[False,False])
        self.assertLessEqual(result.products[model.Key('bundle_pavement',7)].mesh.validate()['maximum_edge_ratio'],8.)

    def test_compound_cropped_ports_retain_their_existing_region_owner(self):
        from types import SimpleNamespace
        from tested_python_streets.bundle_compounds import compound_owners
        network, = height_fixture_networks(model)
        result = PythonMeshBackend(network.settings).prepare(network)
        found = False
        for bundle in network.bundles:
            if not compound_owners(SimpleNamespace(products=result.products), bundle):
                continue
            found = True
            plan = result.products[model.Key('bundle_plan', bundle.id)]
            self.assertNotIsInstance(plan, FiniteAxialPlan)
            self.assertEqual(plan.mode, 'separated')
            self.assertIsNone(result.products[model.Key('bundle_pavement', bundle.id)])
            self.assertIsNotNone(result.products[model.Key('separator', bundle.id)].mesh)
        self.assertTrue(found)

    def test_surface_queries_after_release_and_local_edits_equal_full_rebuilds(self):
        network=skew_bundle_fixture(model)
        engine=PythonMeshBackend(network.settings)
        first=engine.prepare(network)
        saved=meshes(first)
        plan=first.products[model.Key('bundle_plan',7)]
        centers=np.concatenate([m.xyz[m.quads].mean(axis=1) for m in plan.cap_meshes.values()])
        np.testing.assert_allclose(plan.height(centers[:,:2]),centers[:,2],atol=1e-10,rtol=0.)
        engine.commit(first); engine.release_published_meshes()
        self.assertIsNone(plan.combined)
        self.assertFalse(plan.cap_meshes)
        np.testing.assert_allclose(plan.height(centers[:,:2]),centers[:,2],atol=1e-10,rtol=0.)
        placed=replace(network,furniture=tuple(replace(p,station=5.2) if p.id==84 else p for p in network.furniture))
        local=engine.prepare(placed)
        fresh=PythonMeshBackend(placed.settings).prepare(placed)
        self.assertEqual(local.rebuilt,(model.Key('placement',84),))
        self.assertEqual(local.products[model.Key('placement',84)],fresh.products[model.Key('placement',84)])
        engine.commit(local)
        first_section=network.sections[0]
        points=tuple((x,y+.05 if i==1 else y,z) for i,(x,y,z) in enumerate(first_section.points))
        edited=replace(placed,sections=(replace(first_section,points=points),*network.sections[1:]))
        local=engine.prepare(edited)
        fresh=PythonMeshBackend(edited.settings).prepare(edited)
        published=dict(saved);published.update(meshes(local))
        for key in local.removed:published.pop(key,None)
        self.equal(published,meshes(fresh))
        self.assertNotIn(model.Key('road',99),local.rebuilt)
        self.assertNotIn(model.Key('placement',99),local.rebuilt)
        engine.commit(local)
        restored=engine.prepare(network)
        published.update(meshes(restored))
        for key in restored.removed:published.pop(key,None)
        self.equal(published,saved)


if __name__=='__main__':unittest.main(verbosity=2)
