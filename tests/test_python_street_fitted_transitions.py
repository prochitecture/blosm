"""Sharp Exact seams retain physical profiles, certified curves and ownership."""
from dataclasses import replace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.paths import Path
from tested_python_streets.intersections import Approach
from tested_python_streets.profiles import Crown
from tested_python_streets.transitions import TransitionClearance
from tested_python_streets.transition_fitting import ExactTransitionFitter, ExactRowSampler
from tested_python_streets.exact_transitions import PolynomialCurb, PiecewiseCurb
from _python_street_fixtures import two_merged_terminal_fixture
from test_python_street_dimensional_bundles import meshes, reverse_sections


def network(*, length=30., angle=150., resolution=1, datum=0., offsets=(0.,0.)):
    theta=math.radians(angle)
    sources=(model.Section(1,10,((-length,0.,0.),(0.,0.,0.)),width=6.,offset=offsets[0],interpolation=0,datum=datum),
        model.Section(2,10,((0.,0.,0.),(length*math.cos(theta),length*math.sin(theta),0.)),
                      width=10.,offset=offsets[1],interpolation=0,datum=datum),
        model.Section(99,99,((1000.,-100.,0.),(1040.,-100.,0.)),datum=datum))
    return model.Network('sharp_exact_dimensional_seam',sources,
        connections=(model.Connection(8,(model.End(1,1),model.End(2,0))),),
        settings=model.Settings(resolution=resolution))


class PythonFittedTransitionTests(unittest.TestCase):
    def build(self,source):
        authored=source.to_dict();result=PythonMeshBackend(source.settings).prepare(source)
        self.assertEqual(source.to_dict(),authored)
        plan=result.products[model.Key('seam_layout',8)].transition
        self.assertLessEqual(plan.certification['quality']['maximum_edge_ratio'],8.)
        self.assertLessEqual(plan.certification['physical_slopes']['maximum'],plan.certification['maximum_slope']+1e-8)
        self.assertLessEqual(plan.certification['maximum_entrance_longitudinal_slope'],.03+1e-8)
        return result,plan

    def test_original_twenty_angle_cases_including_both_sharp_turns(self):
        for length in (30.,60.):
            for angle in (5.,35.,65.,110.,150.):
                for sign in (-1.,1.):
                    with self.subTest(length=length,angle=angle*sign):
                        result,plan=self.build(network(length=length,angle=angle*sign))
                        if angle==150.:
                            self.assertEqual(plan.certification['sampling_policy'],'certified_curve_error')
                            self.assertEqual(len(plan.certification['extended_free_ends']),2)
                            self.assertLessEqual(plan.rows*plan.across,1000)
                            self.assertLessEqual(plan.certification['certified_curve_error'],.005)
                        else:
                            self.assertNotIn('sampling_policy',plan.certification)

    def test_dense_independent_curve_error_and_exact_centerline(self):
        result,plan=self.build(network())
        t=np.linspace(0.,1.,12001)
        for curve in plan.curbs:
            points=curve(plan.parameters)
            linear=np.column_stack([np.interp(t,plan.parameters,points[:,i]) for i in (0,1)])
            error=float(np.linalg.norm(curve(t)-linear,axis=1).max())
            self.assertLessEqual(error,plan.certification['certified_curve_error'])
        self.assertTrue(np.isin(plan.curbs[0].knots,plan.parameters).all())
        for a,b in zip(plan.curbs[0].knots,plan.curbs[0].knots[1:]):
            self.assertEqual((np.count_nonzero((plan.parameters>=a)&(plan.parameters<=b))-1)%2,0)
        d=t*sum(plan.distances);first=d<plan.distances[0];expected=np.empty((len(t),2))
        expected[first]=plan.approaches[0].frame(plan.distances[0]-d[first])[0]
        expected[~first]=plan.approaches[1].frame(d[~first]-plan.distances[0])[0]
        np.testing.assert_allclose(plan.surface(t,np.zeros_like(t))[0],expected,atol=1e-10,rtol=0.)

    def test_complete_mouth_profiles_and_remaining_source_spans(self):
        result,plan=self.build(network())
        mesh=result.products[model.Key('seam_mesh',8)].mesh
        q=np.linspace(-1.,1.,plan.across+1)
        for i,approach in enumerate(plan.approaches):
            row=mesh.xyz[mesh.seams[i]];side=-q if i==0 else q
            expected=approach.point(plan.distances[i],approach.offset+.5*approach.width*side)
            np.testing.assert_allclose(row[:,:2],expected,atol=1e-12)
            np.testing.assert_allclose(row[:,2],10.*approach.profile.value(side),atol=1e-12)
            road=result.products[model.Key('road',approach.source.id)]
            self.assertAlmostEqual(road.metadata['end']-road.metadata['start'],1.5,places=9)

    def test_resolution_reverse_source_orientation_datum_and_flat_display(self):
        for resolution in (1,2,4):
            with self.subTest(resolution=resolution):
                source=network(length=60.,resolution=resolution,datum=37.)
                source=replace(source,settings=replace(source.settings,exaggeration=0.))
                source=reverse_sections(source,{1,2})
                self.build(source)

    def test_collinear_control_density_does_not_create_a_short_miter_span(self):
        source=network()
        first=replace(source.sections[0],points=((-30.,0.,0.),(-15.,0.,0.),(-.001,0.,0.),(0.,0.,0.)))
        original,_=self.build(source)
        resampled,plan=self.build(replace(source,sections=(first,*source.sections[1:])))
        self.assertTrue(plan.certification['coherent_collinear_frames'])
        for key,part in meshes(original).items():
            np.testing.assert_allclose(part[0],meshes(resampled)[key][0],atol=1e-10,rtol=0.)
            np.testing.assert_array_equal(part[1],meshes(resampled)[key][1])

    def test_free_terminal_budget_respects_selection_incidence_and_bundle_ownership(self):
        source=network()
        approaches=tuple(Approach(Path(s),e,Crown(s.rise,s.datum)) for s,e in zip(source.sections,(1,0)))
        clearance=TransitionClearance(np.array((12.,12.)),())
        limits,records=ExactTransitionFitter.free_limits(source,approaches,clearance)
        np.testing.assert_allclose(limits,(28.5,28.5));self.assertEqual(len(records),2)
        selected=replace(source,sections=(replace(source.sections[0],selection_ends=(0,)),
            replace(source.sections[1],selection_ends=(1,)),source.sections[2]))
        selected_approaches=tuple(Approach(Path(s),e,Crown(s.rise,s.datum)) for s,e in zip(selected.sections,(1,0)))
        limits,records=ExactTransitionFitter.free_limits(selected,selected_approaches,clearance)
        np.testing.assert_array_equal(limits,clearance.limits);self.assertFalse(records)
        with self.assertRaisesRegex(ValueError,'No bounded dimensional'):
            PythonMeshBackend(selected.settings).prepare(selected)
        attached=replace(source,connections=(*source.connections,model.Connection(9,(model.End(1,0),model.End(99,0)))))
        limits,records=ExactTransitionFitter.free_limits(attached,approaches,clearance)
        np.testing.assert_allclose(limits,(12.,28.5));self.assertEqual(len(records),1)
        bundle=two_merged_terminal_fixture(model,bent_seams=True)
        bundle=replace(bundle,sections=tuple(replace(s,interpolation=0) for s in bundle.sections))
        result=PythonMeshBackend(bundle.settings).prepare(bundle)
        for connection in (2,6):
            plan=result.products[model.Key('seam_layout',connection)].transition
            cap=TransitionClearance(.4*np.array([a.length for a in plan.approaches]),())
            limits,records=ExactTransitionFitter.free_limits(bundle,plan.approaches,cap)
            np.testing.assert_array_equal(limits,cap.limits);self.assertFalse(records)

    def test_sampler_rejects_bad_tolerance_and_preallocation_budget(self):
        _,plan=self.build(network())
        with self.assertRaisesRegex(ValueError,'allocation budget'):
            ExactRowSampler(plan.curbs,.005,1)
        with self.assertRaisesRegex(ValueError,'positive finite'):
            ExactRowSampler(plan.curbs,float('nan'),100)

    def test_sampling_search_expands_its_seed_for_conservative_control_hulls(self):
        curve=PiecewiseCurb((0.,1.),(PolynomialCurb(((0.,0.),(.5,0.),(1.,.08))),))
        sampler=ExactRowSampler((curve,),.005,100)
        self.assertEqual(sampler.counts,(4,))
        for parameters,bound in sampler.candidates(1):
            query=np.linspace(0.,1.,1001);points=curve(parameters)
            linear=np.column_stack([np.interp(query,parameters,points[:,i]) for i in (0,1)])
            self.assertLessEqual(float(np.linalg.norm(curve(query)-linear,axis=1).max()),bound)
            self.assertLessEqual(bound,.005)

    def test_local_edit_full_build_furniture_release_and_restoration(self):
        source=network(length=60.)
        source=replace(source,furniture=(model.Placement(85,1,59.5,'RIGHT',-1.),model.Placement(99,99,20.,'RIGHT',-1.)))
        backend=PythonMeshBackend(source.settings);initial=backend.prepare(source);backend.commit(initial)
        saved=meshes(initial);pose=initial.products[model.Key('placement',85)].position;remote=initial.products[model.Key('road',99)]
        source2=source.sections[1];points=list(source2.points);points[-1]=(points[-1][0]+.02,points[-1][1],0.)
        changed=replace(source,sections=(replace(source.sections[0],width=6.05),replace(source2,points=tuple(points)),source.sections[2]))
        edited=backend.prepare(changed);full=PythonMeshBackend(source.settings).prepare(changed)
        self.assertIs(edited.products[model.Key('road',99)],remote)
        self.assertNotEqual(edited.products[model.Key('placement',85)].position,pose)
        for key,part in meshes(edited).items():
            np.testing.assert_array_equal(part[0],meshes(full)[key][0]);np.testing.assert_array_equal(part[1],meshes(full)[key][1])
        backend.commit(edited);restored=backend.prepare(source)
        for key,part in meshes(restored).items():
            np.testing.assert_array_equal(part[0],saved[key][0]);np.testing.assert_array_equal(part[1],saved[key][1])
        backend.commit(restored);backend.release_published_meshes()
        rebuilt=backend.prepare(source,force=True)
        self.assertEqual(rebuilt.products[model.Key('placement',85)].position,initial.products[model.Key('placement',85)].position)


if __name__=='__main__':
    unittest.main(verbosity=2)
