"""Mixed dimensional seams preserve both source curves and physical contracts."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from test_python_street_fitted_transitions import network
from test_python_street_dimensional_bundles import meshes,reverse_sections


def source_network(*,length=60.,angle=65.,modes=(0,3),resolution=1,datum=0.):
    source=network(length=length,angle=angle,resolution=resolution,datum=datum)
    return replace(source,uuid='mixed_dimensional_seam',sections=(
        replace(source.sections[0],interpolation=modes[0]),
        replace(source.sections[1],interpolation=modes[1]),source.sections[2]))


class PythonMixedTransitionTests(unittest.TestCase):
    def build(self,source):
        authored=source.to_dict()
        result=PythonMeshBackend(source.settings).prepare(source)
        self.assertEqual(source.to_dict(),authored)
        plan=result.products[model.Key('seam_layout',8)].transition
        self.assertTrue(plan.certification['mixed_miter'])
        self.assertLessEqual(plan.certification['quality']['maximum_edge_ratio'],8.)
        self.assertLessEqual(plan.certification['physical_slopes']['maximum'],plan.certification['maximum_slope']+1e-8)
        self.assertLessEqual(plan.certification['maximum_entrance_longitudinal_slope'],.03+1e-8)
        self.assertLessEqual(plan.certification['maximum_source_entrance_longitudinal_bound'],.03+1e-8)
        self.assertLessEqual(plan.certification['certified_curve_error'],source.settings.curve_tolerance)
        return result,plan

    def equal(self,first,last):
        self.assertEqual(first.keys(),last.keys())
        for key in first:
            np.testing.assert_array_equal(first[key][0],last[key][0])
            np.testing.assert_array_equal(first[key][1],last[key][1])

    def test_both_modes_and_four_turn_sizes_with_both_signed_turns(self):
        for length,angles in ((30.,(5.,35.,65.)),(60.,(5.,35.,65.,110.))):
            for angle in angles:
                for sign in (-1.,1.):
                    for modes in ((0,3),(3,0)):
                        with self.subTest(length=length,angle=sign*angle,modes=modes):
                            self.build(source_network(length=length,angle=sign*angle,modes=modes))

    def test_source_centerline_loci_and_mandatory_authored_piece_knots(self):
        for modes in ((0,3),(3,0)):
            _,plan=self.build(source_network(modes=modes))
            by_id={a.source.id:a.path for a in plan.approaches}
            knots=plan.curbs[0].knots
            self.assertTrue(np.isin(knots,plan.parameters).all())
            for (identity,start,end),lo,hi in zip(plan.certification['source_parameter_pieces'],knots,knots[1:]):
                u=np.linspace(0.,1.,501)
                expected=by_id[identity].jet(start+(end-start)*u)[0]
                actual,_=plan.surface(lo+(hi-lo)*u,np.zeros_like(u))
                np.testing.assert_allclose(actual,expected,atol=2e-11,rtol=0.)

    def test_complete_mouth_profiles_and_independent_dense_chord_error(self):
        source=source_network()
        _,plan=self.build(source)
        mesh=plan.mesh(PythonMeshBackend(source.settings).patterns,1.)
        q=np.linspace(-1.,1.,plan.across+1)
        for side,a in enumerate(plan.approaches):
            sign=-1. if side==0 else 1.
            expected=np.column_stack((a.point(plan.distances[side],a.offset+sign*.5*a.width*q),
                                      a.profile.datum+a.profile.value(sign*q)))
            np.testing.assert_array_equal(mesh.xyz[mesh.seams[side]],expected)
        query=np.linspace(0.,1.,12001)
        for curve in plan.curbs:
            points=curve(plan.parameters)
            linear=np.column_stack([np.interp(query,plan.parameters,points[:,axis]) for axis in (0,1)])
            error=float(np.linalg.norm(curve(query)-linear,axis=1).max())
            self.assertLessEqual(error,plan.certification['certified_curve_error']+1e-11)
            self.assertFalse(curve.knots.flags.writeable)
            self.assertTrue(all(not p.controls.flags.writeable for p in curve.segments))

    def test_resolution_reversed_source_storage_rotation_and_common_datum(self):
        for resolution in (1,2,4):
            source=source_network(resolution=resolution,datum=37.)
            source=reverse_sections(source,{1,2})
            source=replace(source,settings=replace(source.settings,exaggeration=0. if resolution==2 else 10.))
            angle=.713
            matrix=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
            sections=[]
            for section in source.sections:
                points=np.array(section.points);points[:,:2]=points[:,:2]@matrix.T+(1234.,-5678.)
                sections.append(replace(section,points=tuple(map(tuple,points))))
            self.build(replace(source,sections=tuple(sections)))

    def test_offsets_and_unequal_rises_keep_authored_parameters(self):
        source=source_network(angle=35.)
        source=replace(source,sections=(replace(source.sections[0],rise=.06,offset=.2),
            replace(source.sections[1],rise=.1,offset=-.1),source.sections[2]))
        _,plan=self.build(source)
        mesh=plan.mesh(PythonMeshBackend(source.settings).patterns,1.)
        self.assertGreaterEqual(float(mesh.xyz[:,2].min()),-1e-12)
        self.assertLessEqual(float(mesh.xyz[:,2].max()),.1+1e-12)

    def test_furniture_metric_stations_follow_native_polynomial_parameters(self):
        for modes in ((0,3),(3,0)):
            for reverse in (False,True):
                source=source_network(modes=modes)
                if reverse:
                    source=reverse_sections(source,{1,2})
                backend=PythonMeshBackend(source.settings)
                initial=backend.prepare(source)
                backend.commit(initial)
                plan=initial.products[model.Key('seam_layout',8)].transition
                original=meshes(initial)
                placements=[]
                expected_positions={}
                for side,approach in enumerate(plan.approaches):
                    for fraction in (0.,.25,.5,.75,1.):
                        distance=fraction*plan.distances[side]
                        station=distance if approach.end==0 else approach.length-distance
                        parameter=float(approach.path.parameter(np.asarray(station)))
                        pieces=plan.certification['source_parameter_pieces']
                        index,first,last=next((i,a,b) for i,(identity,a,b) in enumerate(pieces)
                            if identity==approach.source.id and min(a,b)-1e-9<=parameter<=max(a,b)+1e-9)
                        lo,hi=plan.curbs[0].knots[index:index+2]
                        expected_t=float(lo+(hi-lo)*(parameter-first)/(last-first))
                        centerline=plan.surface(np.asarray(expected_t),np.asarray(0.))[0]
                        np.testing.assert_allclose(centerline,approach.path.frame(station)[0],atol=3e-10,rtol=0.)
                        for placement_side in ('LEFT','RIGHT'):
                            with self.subTest(modes=modes,reverse=reverse,source=approach.source.id,
                                              fraction=fraction,side=placement_side):
                                placement=model.Placement(7085+len(placements),approach.source.id,
                                                          station,placement_side,-1.)
                                placements.append(placement)
                                direction=(-1 if side==0 else 1)*(1 if approach.end==0 else -1)
                                sign=direction*(1 if placement_side=='RIGHT' else -1)
                                left,right=(curb(np.asarray(expected_t)) for curb in plan.curbs)
                                transverse=right-left
                                transverse/=np.linalg.norm(transverse)
                                expected=(right if sign==1 else left)+sign*placement.offset*transverse
                                expected_positions[placement.id]=expected
                changed=replace(source,furniture=tuple(placements))
                result=backend.prepare(changed)
                self.equal(meshes(result),original)
                full=PythonMeshBackend(source.settings).prepare(changed)
                for placement in placements:
                    key=model.Key('placement',placement.id)
                    np.testing.assert_allclose(result.products[key].position[:2],
                                               expected_positions[placement.id],atol=3e-10,rtol=0.)
                    self.assertEqual(result.products[key],full.products[key])
                backend.release_published_meshes()
                released=backend.prepare(changed,force=True)
                for placement in placements:
                    key=model.Key('placement',placement.id)
                    self.assertEqual(released.products[key],result.products[key])

    def test_explicit_update_full_restoration_remote_identity_and_released_queries(self):
        source=source_network(angle=35.)
        source=replace(source,furniture=(model.Placement(85,1,59.5,'RIGHT',-1.),
                                        model.Placement(99,99,20.,'RIGHT',-1.)))
        backend=PythonMeshBackend(source.settings)
        initial=backend.prepare(source);backend.commit(initial)
        saved=meshes(initial);pose=initial.products[model.Key('placement',85)].position
        remote=initial.products[model.Key('road',99)]
        second=source.sections[1];points=list(second.points)
        points[-1]=(points[-1][0]+.02,points[-1][1],0.)
        changed=replace(source,sections=(replace(source.sections[0],width=6.05),
            replace(second,points=tuple(points)),source.sections[2]))
        local=backend.prepare(changed);full=PythonMeshBackend(source.settings).prepare(changed)
        self.assertIs(local.products[model.Key('road',99)],remote)
        self.assertNotEqual(local.products[model.Key('placement',85)].position,pose)
        self.equal(meshes(local),meshes(full))
        backend.commit(local);restored=backend.prepare(source)
        self.equal(meshes(restored),saved)
        backend.commit(restored);backend.release_published_meshes()
        rebuilt=backend.prepare(source,force=True)
        self.equal(meshes(rebuilt),saved)
        self.assertEqual(rebuilt.products[model.Key('placement',85)].position,pose)

    def test_invalid_smooth_offset_is_rejected_before_transition_construction(self):
        with self.assertRaisesRegex(ValueError,'Smooth offset curvature'):
            PythonMeshBackend(source_network(length=30.,angle=150.,modes=(0,3)).settings).prepare(
                source_network(length=30.,angle=150.,modes=(0,3)))

    def test_authored_gallery_factories_keep_the_previous_translated_case_exact(self):
        from _python_mixed_fixtures import mixed_transition_networks,translated_corner_seed_network
        translated=translated_corner_seed_network(model)
        # Captured from the accepted native scene. A double-rounding factory
        # changes these controls by an ulp and can redistribute interior quads.
        points={17:(6,(1033.2587890625,-522.2228393554688,0.)),
                117:(5,(1003.9207153320312,-460.1925964355469,0.)),
                102:(1,(944.8456420898438,-555.1543579101562,0.)),
                112:(1,(944.8456420898438,-444.8456726074219,0.))}
        for identity,(index,expected) in points.items():
            self.assertEqual(translated.section_index[identity].points[index],expected)
        self.assertTrue({904,905}<={p.id for p in translated.furniture})
        for network in mixed_transition_networks(model):
            self.build(network)


if __name__=='__main__':
    unittest.main(verbosity=2)
