"""One separated/skew cap opposite merged pavement retains physical inputs."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import one_sided_free_cap_fixture
from test_python_street_dimensional_bundles import meshes,reverse_sections
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend


class PythonOneSidedFreeCapTests(unittest.TestCase):
    def equal(self,a,b):
        self.assertEqual(a.keys(),b.keys())
        for key in a:
            np.testing.assert_array_equal(a[key][0],b[key][0])
            np.testing.assert_array_equal(a[key][1],b[key][1])

    def qualify(self,network,merged_end=1):
        source=network.to_dict()
        update=PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(network.to_dict(),source)
        plan=update.products[model.Key('bundle_plan',7)]
        self.assertEqual(type(plan).__name__,'MergedFreeCapPlan')
        self.assertEqual(plan.mode,'union')
        part=update.products[model.Key('bundle_pavement',7)]
        self.assertEqual(part.metadata['finite_caps'],2)
        self.assertEqual(part.metadata['merged_free_end'],merged_end)
        self.assertEqual(part.metadata['finite_strip_ends'],[1-merged_end])
        self.assertLessEqual(part.metadata['maximum_entrance_longitudinal_slope_bound'],.03+1e-8)
        for kind in ('bundle_pavement','separator'):
            mesh=update.products[model.Key(kind,7)].mesh
            self.assertLessEqual(mesh.validate()['maximum_edge_ratio'],8.)
            self.assertGreaterEqual(float(mesh.xyz[:,2].min()),plan.datum-1e-8)
            self.assertLessEqual(float(mesh.xyz[:,2].max()),plan.datum+plan.rise*network.settings.exaggeration+1e-8)
        self.assertEqual(update.products[model.Key('separator',7)].metadata['islands'],1)
        mesh=part.mesh
        edges=np.sort(np.stack((mesh.quads,np.roll(mesh.quads,-1,axis=1)),axis=-1).reshape(-1,2),axis=1)
        unique,counts=np.unique(edges,axis=0,return_counts=True)
        self.assertLessEqual(int(counts.max()),2)
        self.assertEqual(len(mesh.xyz)-len(unique)+len(mesh.quads),1)
        x=(mesh.xyz[:,:2]-plan.origin)@plan.axis
        plane=plan.first if merged_end==1 else plan.last
        self.assertFalse(np.any(np.all(abs(x[unique[counts==1]]-plane)<1e-8,axis=1)))
        self.assertFalse(any(key.kind in ('road','seam_mesh') and key.id!=99 for key in meshes(update)))
        return update

    def test_captured_edit_and_both_common_traversal_directions(self):
        for reverse in (False,True):
            with self.subTest(reverse=reverse):
                network=one_sided_free_cap_fixture(model,reversed_direction=reverse)
                update=self.qualify(network,merged_end=0 if reverse else 1)
                self.assertLess(len(update.products[model.Key('bundle_pavement',7)].mesh.quads),10000)

    def test_actual_separated_profile_and_merged_end_curve_remain_complete(self):
        network=one_sided_free_cap_fixture(model)
        update=self.qualify(network)
        plan=update.products[model.Key('bundle_plan',7)]
        for role,member in enumerate(network.bundles[0].members):
            visit=member.path[0]
            path=update.products[model.Key('path',visit.section_id)]
            q=np.linspace(-1.,1.,plan.port_intervals+1)
            expected=path.point(0.,path.section.offset+.5*path.section.width*q)
            cap=plan.cap_meshes[role,0]
            distances=np.linalg.norm(cap.xyz[:,None,:2]-expected[None,:,:],axis=-1)
            self.assertLess(float(distances.min(axis=0).max()),1e-8)
            profile=update.products[model.Key('profile',visit.section_id)]
            np.testing.assert_allclose(cap.xyz[distances.argmin(axis=0),2],
                profile.datum+network.settings.exaggeration*profile.value(q),atol=1e-10)
        mesh=update.products[model.Key('bundle_pavement',7)].mesh
        x=(mesh.xyz[:,:2]-plan.origin)@plan.axis
        end=mesh.xyz[abs(x-plan.last)<1e-8]
        self.assertAlmostEqual(float(end[:,1].min()),-9.,places=8)
        self.assertAlmostEqual(float(end[:,1].max()),9.,places=8)
        self.assertAlmostEqual(float(end[:,2].max()),.8,places=7)

    def test_both_modes_displacements_and_physical_slope_limits(self):
        for mode in (0,3):
            for move in (-.2,.05,.2,1.):
                with self.subTest(mode=mode,move=move):
                    network=one_sided_free_cap_fixture(model,move=move)
                    network=replace(network,sections=tuple(replace(s,interpolation=mode) for s in network.sections))
                    update=self.qualify(network)
                    metadata=update.products[model.Key('bundle_pavement',7)].metadata
                    self.assertLessEqual(metadata['physical_slopes']['maximum'],metadata['measured_slope_limit']+1e-8)

    def test_reversed_sections_fine_resolutions_and_zero_exaggeration(self):
        for mode in (0,3):
            for resolution in (1,2,3,4):
                with self.subTest(mode=mode,resolution=resolution):
                    network=one_sided_free_cap_fixture(model)
                    network=replace(network,sections=tuple(replace(s,interpolation=mode,datum=37.) for s in network.sections),
                        settings=replace(network.settings,resolution=resolution,exaggeration=0. if resolution==2 else 10.))
                    self.qualify(reverse_sections(network,{1,4}))

    def test_local_edits_match_full_rebuild_and_restore_queries_after_release(self):
        network=one_sided_free_cap_fixture(model)
        backend=PythonMeshBackend(network.settings)
        initial=backend.prepare(network);backend.commit(initial)
        saved=meshes(initial)
        remote=initial.products[model.Key('road',99)]
        changed=one_sided_free_cap_fixture(model,move=.1)
        for candidate in (changed,replace(network,sections=tuple(replace(s,width=9.5) if s.id in (2,4) else s for s in network.sections))):
            local=backend.prepare(candidate)
            self.equal(meshes(local),meshes(PythonMeshBackend(candidate.settings).prepare(candidate)))
            self.assertIs(local.products[model.Key('road',99)],remote)
            self.assertNotIn(model.Key('placement',99),local.rebuilt)
            backend.commit(local)
            restored=backend.prepare(network);self.equal(meshes(restored),saved);backend.commit(restored)
        plan=backend.updates.products[model.Key('bundle_plan',7)]
        position=np.asarray(initial.products[model.Key('placement',84)].position[:2])
        height=float(plan.height(position))
        backend.release_published_meshes()
        self.assertFalse(plan.cap_meshes)
        self.assertAlmostEqual(float(plan.height(position)),height,places=12)
        self.equal(meshes(backend.prepare(network,force=True)),saved)

    def test_rotated_coordinates_preserve_the_profile_and_curbs(self):
        network=one_sided_free_cap_fixture(model)
        angle=.731
        rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
        sections=[]
        for source in network.sections:
            p=np.asarray(source.points).copy();p[:,:2]=p[:,:2]@rotation.T+np.array((1234.,-5678.))
            sections.append(replace(source,points=tuple(map(tuple,p))))
        self.qualify(replace(network,sections=tuple(sections)))

    def test_constant_width_members_use_the_same_cap_ownership(self):
        upper=((-40.,8.,0.),(-30.,8.05,0.),(-20.,8.,0.),(0.,2.,0.),(20.,2.,0.),(40.,2.,0.))
        lower=tuple((x,-y if i!=1 else -8.,z) for i,(x,y,z) in enumerate(upper))
        for mode in (0,3):
            with self.subTest(mode=mode):
                sections=(model.Section(1,1,upper,interpolation=mode),model.Section(3,2,lower,interpolation=mode),
                          model.Section(99,99,((-40.,35.,0.),(40.,35.,0.))))
                bundle=model.Bundle(7,(model.BundleMember((model.Traversal(1),)),model.BundleMember((model.Traversal(3),))))
                network=model.Network('one_sided_constant_width',sections,bundles=(bundle,))
                update=self.qualify(network)
                self.assertEqual(update.products[model.Key('bundle_plan',7)].transition_owners,())


if __name__=='__main__':unittest.main(verbosity=2)
