"""Equal-profile mixed seams use certified geometry and persistent ownership."""
from dataclasses import replace
import math
import unittest
from types import SimpleNamespace
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.providers import PathProvider
from tested_python_streets.intersections import Approach
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.transitions import required
from test_python_street_mixed_transitions import source_network
from test_python_street_dimensional_bundles import meshes,reverse_sections
from _python_mixed_fixtures import uniform_mixed_networks,uniform_mixed_bundle_network


def uniform_source(**options):
    source=source_network(**options)
    return replace(source,uuid='uniform_mixed_seam',sections=(source.sections[0],
        replace(source.sections[1],width=source.sections[0].width),source.sections[2]))


class PythonUniformMixedSeamTests(unittest.TestCase):
    def build(self,source):
        authored=source.to_dict()
        self.assertTrue(required(source,source.connections[0]))
        result=PythonMeshBackend(source.settings).prepare(source)
        self.assertEqual(source.to_dict(),authored)
        plan=result.products[model.Key('seam_layout',8)].transition
        self.assertTrue(plan.certification['mixed_miter'])
        self.assertIsNotNone(result.products[model.Key('seam_mesh',8)].mesh)
        self.assertLessEqual(plan.certification['quality']['maximum_edge_ratio'],8.)
        self.assertLessEqual(plan.certification['physical_slopes']['maximum'],plan.certification['maximum_slope']+1e-8)
        self.assertLessEqual(plan.certification['maximum_source_entrance_longitudinal_bound'],.03+1e-8)
        return result,plan

    def equal(self,a,b):
        self.assertEqual(a.keys(),b.keys())
        for key in a:
            np.testing.assert_array_equal(a[key][0],b[key][0])
            np.testing.assert_array_equal(a[key][1],b[key][1])

    def test_source_dispatch_matches_actual_mouth_frames(self):
        for modes in ((0,3),(3,0),(3,3),(0,0)):
            for angle in (0.,.0001,.0012,5.,65.):
                source=uniform_source(modes=modes,angle=angle)
                site=source.connections[0]
                context=SimpleNamespace(network=source)
                approaches=[Approach(PathProvider.build(context,source.section_index[e.section_id],(site,)),
                    e.end,ProfileProvider.build(source.section_index[e.section_id],None,source.section_index))
                    for e in site.ends]
                mismatch=np.linalg.norm(approaches[0].direction+approaches[1].direction)>1e-5
                expected=modes[0]!=modes[1] and mismatch
                self.assertEqual(required(source,site),expected,(modes,angle))

    def test_normal_dispatch_preserves_complete_profiles_for_both_orders_and_signed_turns(self):
        for length,angles in ((30.,(5.,35.,65.,110.)),(60.,(5.,35.,65.,110.))):
            for angle in angles:
                for sign in (-1.,1.):
                    for modes in ((0,3),(3,0)):
                        with self.subTest(length=length,angle=sign*angle,modes=modes):
                            source=uniform_source(length=length,angle=sign*angle,modes=modes)
                            result,plan=self.build(source)
                            mesh=plan.mesh(PythonMeshBackend(source.settings).patterns,1.)
                            q=np.linspace(-1.,1.,plan.across+1)
                            for side,a in enumerate(plan.approaches):
                                outward_q=-q if side==0 else q
                                expected=np.column_stack((a.point(plan.distances[side],a.offset+.5*a.width*outward_q),
                                    a.profile.datum+a.profile.value(outward_q)))
                                np.testing.assert_array_equal(mesh.xyz[mesh.seams[side]],expected)
                            self.assertGreaterEqual(float(mesh.xyz[:,2].min()),-1e-12)
                            self.assertLessEqual(float(mesh.xyz[:,2].max()),.08+1e-12)

    def test_uniform_profile_and_original_polynomial_loci(self):
        source=uniform_source()
        _,plan=self.build(source)
        t=np.linspace(0.,1.,601)
        for q in (-1.,-.5,0.,.5,1.):
            np.testing.assert_allclose(plan.surface(t,q)[1],.08*(1.-q*q),atol=1e-14,rtol=0.)
        paths={a.source.id:a.path for a in plan.approaches}
        for index,(identity,first,last) in enumerate(plan.certification['source_parameter_pieces']):
            lo,hi=plan.curbs[0].knots[index:index+2]
            u=np.linspace(0.,1.,501)
            expected=paths[identity].jet(first+(last-first)*u)[0]
            np.testing.assert_allclose(plan.surface(lo+(hi-lo)*u,0.)[0],expected,atol=3e-10,rtol=0.)

    def test_source_orientation_rotation_datum_and_resolution(self):
        for resolution in (1,2,4):
            source=reverse_sections(uniform_source(resolution=resolution,datum=37.),{1,2})
            angle=.713
            matrix=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
            sections=[]
            for section in source.sections:
                points=np.asarray(section.points);points[:,:2]=points[:,:2]@matrix.T+(1234.,-5678.)
                sections.append(replace(section,points=tuple(map(tuple,points))))
            source=replace(source,sections=tuple(sections),settings=replace(source.settings,
                exaggeration=0. if resolution==2 else 10.))
            self.build(source)

    def test_mode_and_dimension_edits_remove_restore_owned_surface_and_anchors(self):
        for source in uniform_mixed_networks(model):
            backend=PythonMeshBackend(source.settings)
            initial=backend.prepare(source);backend.commit(initial)
            saved=meshes(initial)
            poses={key:value for key,value in initial.products.items() if key.kind=='placement'}
            remote=initial.products[model.Key('road',99)]
            smooth=replace(source,sections=tuple(replace(section,interpolation=3) if section.id in (1,2)
                else section for section in source.sections))
            result=backend.prepare(smooth)
            self.assertIsNone(result.products[model.Key('seam_layout',8)].transition)
            self.assertIn(model.Key('seam_mesh',8),result.removed)
            self.assertNotIn(model.Key('seam_surface',8),result.products)
            self.equal(meshes(result),meshes(PythonMeshBackend(smooth.settings).prepare(smooth)))
            self.assertIs(result.products[model.Key('road',99)],remote)
            backend.commit(result)
            restored=backend.prepare(source)
            self.equal(meshes(restored),saved)
            for key,value in poses.items():self.assertEqual(restored.products[key],value)
            backend.commit(restored)
            second=source.sections[1];points=list(second.points)
            points[-1]=(points[-1][0]+.025,points[-1][1],points[-1][2])
            changed=replace(source,sections=(replace(source.sections[0],width=6.02),
                replace(second,points=tuple(points)),source.sections[2]))
            edited=backend.prepare(changed);full=PythonMeshBackend(changed.settings).prepare(changed)
            self.equal(meshes(edited),meshes(full))
            for key in poses:self.assertEqual(edited.products[key],full.products[key])
            self.assertIs(edited.products[model.Key('road',99)],remote)
            backend.commit(edited)
            self.equal(meshes(backend.prepare(source)),saved)

    def test_neighboring_junction_remains_a_physical_clearance_constraint(self):
        source=uniform_source(angle=35.)
        point=np.asarray(source.sections[1].points[-1])
        added=tuple(model.Section(identity,identity,(tuple(point),tuple(point+delta)),width=6.)
            for identity,delta in ((3,np.array((30.,0.,0.))),(4,np.array((0.,30.,0.)))))
        junction=model.Intersection(12,(model.End(2,1),model.End(3,0),model.End(4,0)))
        source=replace(source,sections=(*source.sections,*added),intersections=(junction,))
        result,plan=self.build(source)
        self.assertTrue(any(row['junction']==junction.key.text() for row in plan.certification['clearance_constraints']))
        from tested_python_streets.ownership import footprint_key
        self.assertIn(footprint_key(junction.key),result.dependencies[model.Key('seam_layout',8)])

    def test_bundle_members_keep_shared_profile_separator_and_single_surface_owner(self):
        for modes in ((0,3),(3,0)):
            for overlap in (False,True):
                with self.subTest(modes=modes,overlap=overlap):
                    source=uniform_mixed_bundle_network(model,modes,overlap=overlap)
                    authored=source.to_dict()
                    result=PythonMeshBackend(source.settings).prepare(source)
                    self.assertEqual(source.to_dict(),authored)
                    plan=result.products[model.Key('bundle_plan',7)]
                    self.assertEqual(plan.mode,'union' if overlap else 'separated')
                    if overlap:
                        self.assertEqual({site.key for site in source.connections},set(plan.transition_owners))
                    for site in source.connections:
                        self.assertTrue(required(source,site))
                        part=result.products[model.Key('seam_mesh',site.id)]
                        if overlap:self.assertIsNone(part)
                        else:self.assertIsNotNone(part.mesh)
                    for section in source.sections[:4]:
                        part=result.products[model.Key('road',section.id)]
                        self.assertEqual(getattr(part,'mesh',None) is None,overlap)
                    pavement=result.products[model.Key('bundle_pavement',7)]
                    separator=result.products[model.Key('separator',7)]
                    if overlap:
                        self.assertIsNotNone(pavement.mesh)
                        self.assertIsNone(separator)
                        checked=(pavement,)
                    else:
                        self.assertIsNone(pavement)
                        self.assertIsNotNone(separator.mesh)
                        checked=(separator,)
                    for part in checked:
                        self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'],8.)
                        self.assertGreaterEqual(float(part.mesh.xyz[:,2].min()),-1e-8)
                        self.assertLessEqual(float(part.mesh.xyz[:,2].max()),.8+1e-8)

    def test_invalid_common_datum_preserves_last_committed_state(self):
        source=uniform_source()
        backend=PythonMeshBackend(source.settings)
        initial=backend.prepare(source);backend.commit(initial)
        saved=meshes(initial);manifest=backend.manifest()
        for modes in ((0,3),(3,0)):
            for angle in (-110.,110.):
                bad=uniform_source(length=30.,angle=angle,modes=modes)
                bad=replace(bad,sections=tuple(replace(section,datum=.01) if section.id==2
                    else section for section in bad.sections))
                authored=bad.to_dict()
                with self.assertRaisesRegex(ValueError,'common minimum height'):
                    backend.prepare(bad)
                self.assertEqual(bad.to_dict(),authored)
                self.assertEqual(backend.manifest(),manifest)
        self.equal(meshes(backend.prepare(source)),saved)


if __name__=='__main__':unittest.main(verbosity=2)
