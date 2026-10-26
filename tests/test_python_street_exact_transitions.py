"""Exact dimensional seams keep source bends and full physical profiles."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.exact_transitions import PiecewiseCurb
from _python_street_fixtures import two_merged_terminal_fixture, exact_bent_seam_fixture
from test_python_street_dimensional_bundles import meshes, reverse_sections


def source_network(*, offsets=(0.,0.), rises=(.08,.08), datum=0.):
    network=exact_bent_seam_fixture(model,offsets=offsets,rises=rises,datum=datum)
    return replace(network,sections=tuple(replace(s,interpolation=0) if s.id!=99 else s for s in network.sections))


class PythonExactTransitionTests(unittest.TestCase):
    def equal(self,a,b):
        self.assertEqual(a.keys(),b.keys())
        for key in a:
            np.testing.assert_array_equal(a[key][0],b[key][0])
            np.testing.assert_array_equal(a[key][1],b[key][1])

    def build(self,network):
        authored = network.to_dict()
        result = PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(network.to_dict(),authored)
        plan = result.products[model.Key('seam_layout',8)].transition
        self.assertTrue(plan.certification['exact_miter'])
        self.assertTrue(all(isinstance(c,PiecewiseCurb) for c in plan.curbs))
        self.assertLessEqual(result.products[model.Key('seam_mesh',8)].mesh.validate()['maximum_edge_ratio'],8.)
        self.assertLessEqual(plan.certification['maximum_entrance_longitudinal_slope'],.03+1e-8)
        return result,plan

    def test_extra_authored_bends_are_compulsory_and_centerline_remains_exact(self):
        result,plan = self.build(source_network())
        self.assertGreater(plan.certification['polynomial_pieces'],2)
        knots = plan.curbs[0].knots
        self.assertTrue(np.isin(knots,plan.parameters).all())
        t = np.unique(np.r_[np.linspace(0.,1.,2001),knots])
        distance = t*sum(plan.distances)
        first = distance < plan.distances[0]
        expected = np.empty((len(t),2))
        expected[first] = plan.approaches[0].frame(plan.distances[0]-distance[first])[0]
        expected[~first] = plan.approaches[1].frame(distance[~first]-plan.distances[0])[0]
        xy,_ = plan.surface(t,np.zeros_like(t))
        np.testing.assert_allclose(xy,expected,atol=1e-10,rtol=0.)

    def test_complete_mouth_profiles_and_piecewise_sampling_bound(self):
        result,plan = self.build(source_network())
        mesh = plan.mesh(PythonMeshBackend(source_network().settings).patterns,exaggeration=1.)
        q = np.linspace(-1.,1.,plan.across+1)
        for i,approach in enumerate(plan.approaches):
            sign = -1. if i==0 else 1.
            xy = approach.point(plan.distances[i],approach.offset+sign*.5*approach.width*q)
            expected = np.column_stack((xy,approach.profile.datum+approach.profile.value(sign*q)))
            row = mesh.xyz[:plan.across+1] if i==0 else mesh.xyz[-plan.across-1:]
            np.testing.assert_allclose(row,expected,atol=1e-10,rtol=0.)
        for curb in plan.curbs:
            parameters = curb.sampling_parameters(.7,.005,max_intervals=10000)
            query = np.linspace(0.,1.,4001)
            actual = curb(query);samples = curb(parameters)
            chord = np.column_stack([np.interp(query,parameters,samples[:,i]) for i in range(2)])
            self.assertLessEqual(float(np.linalg.norm(actual-chord,axis=1).max()),.005+1e-10)
            self.assertFalse(curb.knots.flags.writeable)
            self.assertTrue(all(not p.controls.flags.writeable for p in curb.segments))
            with self.assertRaisesRegex(ValueError,'allocation budget'):
                curb.sampling_parameters(.7,.005,max_intervals=1)

    def test_reversed_authored_sections_rotation_offsets_and_common_datum(self):
        network = source_network(offsets=(.5,-.3),rises=(.06,.1),datum=37.)
        for reverse in (False,True):
            with self.subTest(reverse=reverse):
                candidate = reverse_sections(network,{1,2}) if reverse else network
                self.build(candidate)
        angle=.713;rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
        sources=[]
        for s in network.sections:
            p=np.asarray(s.points).copy();p[:,:2]=p[:,:2]@rotation.T+(1234.,-5678.)
            sources.append(replace(s,points=tuple(map(tuple,p))))
        self.build(replace(network,sections=tuple(sources)))

    def test_finer_mesh_and_zero_display_height_keep_physical_checks(self):
        network=source_network(datum=37.)
        for resolution in (1,2,4):
            with self.subTest(resolution=resolution):
                candidate=replace(network,settings=replace(network.settings,resolution=resolution,
                    exaggeration=0. if resolution==2 else 10.))
                result,plan=self.build(candidate)
                self.assertGreater(plan.certification['physical_slopes']['maximum'],0.)
                self.assertLessEqual(plan.certification['physical_slopes']['maximum'],plan.certification['maximum_slope']+1e-8)

    def test_original_bent_bundle_has_two_merged_caps_and_one_separator_island(self):
        network=two_merged_terminal_fixture(model,bent_seams=True)
        network=replace(network,sections=tuple(replace(s,interpolation=0) for s in network.sections))
        authored=network.to_dict();result=PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(network.to_dict(),authored)
        plan=result.products[model.Key('bundle_plan',7)]
        self.assertEqual(plan.merged_ends,(0,1))
        self.assertEqual(result.products[model.Key('separator',7)].metadata['islands'],1)
        self.assertLessEqual(result.products[model.Key('bundle_pavement',7)].mesh.validate()['maximum_edge_ratio'],8.)
        for connection in (2,6):
            transition=result.products[model.Key('seam_layout',connection)].transition
            self.assertTrue(transition.certification['exact_miter'])
            self.assertIsNone(result.products[model.Key('seam_mesh',connection)])

    def test_bent_bundle_reversed_sources_and_flat_display_keep_physical_ownership(self):
        original=two_merged_terminal_fixture(model,bent_seams=True)
        original=replace(original,sections=tuple(replace(s,interpolation=0,datum=37.) for s in original.sections))
        for reverse,resolution in ((False,1),(True,2)):
            with self.subTest(reverse=reverse,resolution=resolution):
                network=reverse_sections(original,{1,4,6}) if reverse else original
                network=replace(network,settings=replace(network.settings,resolution=resolution,exaggeration=0.))
                result=PythonMeshBackend(network.settings).prepare(network)
                self.assertEqual(result.products[model.Key('bundle_plan',7)].merged_ends,(0,1))
                self.assertEqual(result.products[model.Key('separator',7)].metadata['islands'],1)
                for connection in (2,6):
                    plan=result.products[model.Key('seam_layout',connection)].transition
                    self.assertTrue(plan.certification['exact_miter'])
                    self.assertGreater(plan.certification['physical_slopes']['maximum'],0.)
                self.assertTrue(np.isfinite(result.products[model.Key('placement',87)].position).all())

    def test_local_edit_full_rebuild_release_and_restoration_keep_dependencies(self):
        network=source_network();backend=PythonMeshBackend(network.settings)
        initial=backend.prepare(network);backend.commit(initial);saved=meshes(initial)
        pose=initial.products[model.Key('placement',85)].position;remote=initial.products[model.Key('road',99)]
        section=network.sections[0];points=list(section.points);points[-2]=(-6.,.3,0.)
        changed=replace(network,sections=(replace(section,points=tuple(points),width=6.1),*network.sections[1:]))
        local=backend.prepare(changed)
        self.equal(meshes(local),meshes(PythonMeshBackend(network.settings).prepare(changed)))
        self.assertIs(local.products[model.Key('road',99)],remote)
        self.assertNotEqual(local.products[model.Key('placement',85)].position,pose)
        backend.commit(local);restored=backend.prepare(network)
        self.equal(meshes(restored),saved)
        backend.commit(restored);backend.release_published_meshes()
        self.equal(meshes(backend.prepare(network,force=True)),saved)


if __name__=='__main__':
    unittest.main(verbosity=2)
