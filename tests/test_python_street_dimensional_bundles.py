"""Certified width transitions participate in Bundle union and Update ownership."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import dimensional_bundle_fixture, dimensional_bundle_fixture_networks
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_boundaries import dimensional_sites
from tested_python_streets.transitions import output_key
from tested_python_streets.paths import Path


def meshes(update):
    return {k: (p.mesh.xyz.copy(), p.mesh.quads.copy()) for k,p in update.products.items()
            if getattr(p, 'mesh', None) is not None}


def reverse_sections(network, ids):
    lengths = {sid: Path(network.section_index[sid]).length for sid in ids}
    swap = lambda e: model.End(e.section_id, 1-e.end) if e.section_id in ids else e
    return replace(network,
        sections=tuple(replace(s, points=s.points[::-1], offset=-s.offset) if s.id in ids else s for s in network.sections),
        intersections=tuple(replace(j, ends=tuple(map(swap,j.ends))) for j in network.intersections),
        connections=tuple(replace(j, ends=tuple(map(swap,j.ends))) for j in network.connections),
        bundles=tuple(replace(b, members=tuple(replace(m, path=tuple(
            replace(t, forward=not t.forward) if t.section_id in ids else t for t in m.path)) for m in b.members)) for b in network.bundles),
        furniture=tuple(replace(p, station=lengths[p.section_id]-p.station,
            side='RIGHT' if p.side == 'LEFT' else 'LEFT') if p.section_id in ids else p for p in network.furniture))


class PythonDimensionalBundleTests(unittest.TestCase):
    def assert_meshes_equal(self, a, b):
        self.assertEqual(a.keys(), b.keys())
        for key in a:
            np.testing.assert_array_equal(a[key][0], b[key][0], err_msg=key.text())
            np.testing.assert_array_equal(a[key][1], b[key][1], err_msg=key.text())

    def qualify(self, network):
        source = network.to_dict()
        update = PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(network.to_dict(), source)
        plan = update.products[model.Key('bundle_plan', 7)]
        self.assertEqual(plan.mode, 'union')
        sites = dimensional_sites(network, network.bundles[0])
        self.assertEqual(set(plan.transition_owners), {j.key for j in sites})
        for site in sites:
            self.assertIsNone(update.products[output_key(site.key)])
            self.assertIn(site.layout_key, update.dependencies[model.Key('bundle_plan',7)])
            self.assertIn(model.Key('bundle_plan',7), update.dependencies[output_key(site.key)])
        for sid in network.bundles[0].sections:
            self.assertIsNone(getattr(update.products[model.Key('road',sid)], 'mesh', None))
        for kind in ('bundle_pavement','separator'):
            part=update.products[model.Key(kind,7)]
            self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
            self.assertGreaterEqual(float(part.mesh.xyz[:,2].min()), plan.datum-1e-8)
            self.assertLessEqual(float(part.mesh.xyz[:,2].max()), plan.datum+plan.rise*network.settings.exaggeration+1e-8)
        self.assertLess(plan.minimum_gap,0.)
        self.assertGreater(plan.maximum_gap,0.)
        for p in network.furniture:
            spans=update.products[model.Key('road',p.section_id)].metadata.get('transition_spans',())
            if not any(s['start']-1e-8 <= p.station <= s['end']+1e-8 for s in spans):
                continue
            transform=update.products[model.Key('placement',p.id)]
            self.assertAlmostEqual(transform.position[2], float(plan.height(np.asarray(transform.position[:2])))+p.height, places=10)
        return update

    def test_straight_bent_and_rounded_attached_width_transitions_both_modes(self):
        for base in dimensional_bundle_fixture_networks(model):
            for mode in (0,3):
                with self.subTest(case=base.uuid, mode=mode):
                    network=replace(base, sections=tuple(replace(s,interpolation=mode) for s in base.sections))
                    update=self.qualify(network)
                    pavement=update.products[model.Key('bundle_pavement',7)]
                    self.assertLess(len(pavement.mesh.quads), 10000)
                    self.assertEqual(update.products[model.Key('separator',7)].metadata['islands'], 2 if network.intersections else 1)

    def test_transition_curbs_match_independent_dense_analytical_samples(self):
        network=dimensional_bundle_fixture(model,bent=True)
        update=self.qualify(network)
        plan=update.products[model.Key('bundle_plan',7)]
        for site in dimensional_sites(network,network.bundles[0]):
            transition=update.products[site.layout_key].transition
            role=0 if transition.approaches[0].source.id in (1,2) else 1
            for curb in transition.curbs:
                reference=curb(np.linspace(0.,1.,2001))-plan.origin
                x,y=reference@plan.axis,reference@plan.left
                if x[-1] < x[0]:
                    x,y=x[::-1],y[::-1]
                curves=plan.curves[2*role:2*role+2]
                errors=[np.max(abs(np.interp(x,c[:,0],c[:,1])-y)) for c in curves]
                self.assertLessEqual(min(errors), network.settings.curve_tolerance*.11)
        # Find the physical inner-curb contact independently, on much denser
        # analytical samples. Unequal bends can move it away from x=40.
        inner=[]
        for role,site in enumerate(dimensional_sites(network,network.bundles[0])):
            transition=update.products[site.layout_key].transition
            curves=[curb(np.linspace(0.,1.,20001))-plan.origin for curb in transition.curbs]
            curve=sorted(curves,key=lambda c:float((c@plan.left).mean()))[0 if role==0 else -1]
            xy=np.column_stack((curve@plan.axis,curve@plan.left))
            inner.append(xy if xy[-1,0]>xy[0,0] else xy[::-1])
        lo=max(c[0,0] for c in inner); hi=min(c[-1,0] for c in inner)
        for _ in range(40):
            mid=.5*(lo+hi)
            gap=np.interp(mid,inner[0][:,0],inner[0][:,1])-np.interp(mid,inner[1][:,0],inner[1][:,1])
            if gap>0.: lo=mid
            else: hi=mid
        root=.5*(lo+hi)
        grass=update.products[model.Key('separator',7)].mesh
        self.assertAlmostEqual(float(((grass.xyz[:,:2]-plan.origin)@plan.axis).max()),root,delta=.005)

    def test_reversed_sources_resolutions_and_zero_display_exaggeration(self):
        for mode in (0,3):
            for resolution in (1,2,3,4):
                with self.subTest(mode=mode,resolution=resolution):
                    network=dimensional_bundle_fixture(model,bent=True)
                    network=replace(network,sections=tuple(replace(s,interpolation=mode,datum=37.) for s in network.sections),
                        settings=replace(network.settings,resolution=resolution,exaggeration=0. if resolution==2 else 10.))
                    self.qualify(reverse_sections(network,{2,3}))

    def test_local_edits_match_clean_builds_preserve_remote_outputs_and_restore(self):
        network=dimensional_bundle_fixture(model,bent=True)
        backend=PythonMeshBackend(network.settings)
        original=backend.prepare(network); backend.commit(original)
        saved=meshes(original)
        remote=original.products[model.Key('road',99)]
        edits=(replace(network,sections=tuple(replace(s,width=9.5) if s.id in (2,4) else s for s in network.sections)),
               replace(network,sections=tuple(replace(s,points=tuple((x,y+.05 if i==2 else y,z) for i,(x,y,z) in enumerate(s.points))) if s.id==1 else s for s in network.sections)))
        for changed in edits:
            local=backend.prepare(changed)
            clean=PythonMeshBackend(changed.settings).prepare(changed)
            self.assert_meshes_equal(meshes(local),meshes(clean))
            self.assertIs(local.products[model.Key('road',99)],remote)
            self.assertNotIn(model.Key('placement',99),local.rebuilt)
            self.assertIn(model.Key('bundle_plan',7),local.rebuilt)
            self.assertEqual(local.products[model.Key('placement',81)],clean.products[model.Key('placement',81)])
            backend.commit(local)
            restored=backend.prepare(network)
            self.assert_meshes_equal(meshes(restored),saved); backend.commit(restored)
        backend.release_published_meshes()
        self.assert_meshes_equal(meshes(backend.prepare(network,force=True)),saved)

    def test_skew_end_strips_keep_each_terminal_sections_actual_width_and_profile(self):
        network=dimensional_bundle_fixture(model,skew=True)
        update=self.qualify(network)
        plan=update.products[model.Key('bundle_plan',7)]
        self.assertEqual(type(plan).__name__,'FiniteAxialPlan')
        metadata=update.products[model.Key('bundle_pavement',7)].metadata
        self.assertLessEqual(metadata['maximum_entrance_longitudinal_slope_bound'],.03+1e-8)
        for role,member in enumerate(network.bundles[0].members):
            for end in (0,1):
                t=member.path[0 if end==0 else -1]
                source=network.section_index[t.section_id]
                path=update.products[model.Key('path',source.id)]
                q=np.linspace(-1.,1.,plan.port_intervals+1)
                expected=path.point(0. if end==0 else path.length,source.offset+.5*source.width*q)
                pattern_mesh=plan.cap_meshes[role,end]
                # The quad atlas need not store a boundary contiguously.
                distances=np.linalg.norm(pattern_mesh.xyz[:,None,:2]-expected[None,:,:],axis=-1)
                self.assertLess(float(distances.min(axis=0).max()),1e-8)
                ids=distances.argmin(axis=0)
                profile=update.products[model.Key('profile',source.id)]
                np.testing.assert_allclose(pattern_mesh.xyz[ids,2],profile.datum+network.settings.exaggeration*profile.value(q),atol=1e-10)

    def test_attached_higher_resolutions_and_rotated_coordinates_keep_source_contracts(self):
        angle=.731
        rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
        for mode in (0,3):
            for resolution in (2,4):
                with self.subTest(mode=mode,resolution=resolution):
                    network=dimensional_bundle_fixture(model,attached=True)
                    sections=[]
                    for source in network.sections:
                        p=np.asarray(source.points).copy()
                        p[:,:2]=p[:,:2]@rotation.T+np.array((1234.,-5678.))
                        sections.append(replace(source,points=tuple(map(tuple,p)),datum=37.,interpolation=mode))
                    network=replace(network,sections=tuple(sections),settings=replace(network.settings,resolution=resolution))
                    self.qualify(reverse_sections(network,{11,2}))

    def test_separated_to_union_switch_retires_and_restores_transition_pavements(self):
        network=dimensional_bundle_fixture(model)
        separated=replace(network,sections=tuple(replace(s,width=7.) if s.id in (2,4) else s for s in network.sections))
        backend=PythonMeshBackend(network.settings)
        before=backend.prepare(separated); backend.commit(before)
        self.assertEqual(before.products[model.Key('bundle_plan',7)].mode,'separated')
        for site in separated.connections:
            self.assertIsNotNone(before.products[output_key(site.key)].mesh)
        union=backend.prepare(network)
        self.assert_meshes_equal(meshes(union),meshes(PythonMeshBackend(network.settings).prepare(network)))
        for site in network.connections:
            self.assertIsNone(union.products[output_key(site.key)])
        backend.commit(union)
        after=backend.prepare(separated)
        self.assert_meshes_equal(meshes(before),meshes(after))

    def test_crossing_or_inconsistent_member_controls_leave_committed_build_unchanged(self):
        network=dimensional_bundle_fixture(model)
        backend=PythonMeshBackend(network.settings)
        backend.commit(backend.prepare(network)); manifest=backend.manifest()
        broken=replace(network,sections=tuple(replace(s,points=tuple((x,4.,z) for x,y,z in s.points)) if s.id in (3,4) else s for s in network.sections))
        source=broken.to_dict()
        with self.assertRaisesRegex(ValueError,'cross or touch'):
            backend.prepare(broken)
        self.assertEqual(backend.manifest(),manifest)
        self.assertEqual(broken.to_dict(),source)

    def test_one_sided_and_staggered_width_changes_keep_offsets_and_remove_obsolete_tasks(self):
        base=dimensional_bundle_fixture(model)
        for mode in (0,3):
            sections=[]
            for s in base.sections:
                width=12. if s.id==2 else 6.
                offset=.5 if s.id==2 else -.3 if s.id==4 else 0.
                seam=-3. if s.id<3 else 3.
                points=tuple((seam if x==0. else x,y,z) for x,y,z in s.points)
                sections.append(replace(s,width=width if s.id!=99 else s.width,offset=offset,
                                        points=points,interpolation=mode))
            network=replace(base,sections=tuple(sections),furniture=())
            result=self.qualify(network)
            self.assertEqual(len(result.products[model.Key('bundle_plan',7)].transition_owners),2)
        # One member changes width; the other keeps a compatible area-free seam.
        one=replace(base,sections=tuple(replace(s,width=12.) if s.id==2 else
            replace(s,width=6.) if s.id==4 else s for s in base.sections))
        result=self.qualify(one)
        self.assertEqual(result.products[model.Key('bundle_plan',7)].transition_owners,(model.Key('connection',40),))
        backend=PythonMeshBackend(one.settings); backend.commit(result)
        compatible=replace(one,sections=tuple(replace(s,width=6.) if s.id==2 else s for s in one.sections))
        local=backend.prepare(compatible)
        self.assertNotIn(model.Key('seam_surface',40),local.products)
        self.assert_meshes_equal(meshes(local),meshes(PythonMeshBackend(compatible.settings).prepare(compatible)))

    def test_ownership_switch_after_buffer_release_reconstructs_an_unchanged_transition(self):
        network=dimensional_bundle_fixture(model)
        backend=PythonMeshBackend(network.settings)
        initial=backend.prepare(network); backend.commit(initial)
        original_product=initial.products[model.Key('seam_surface',40)]
        backend.release_published_meshes()
        self.assertIsNone(original_product.pavement.mesh)
        narrower=replace(network,sections=tuple(replace(s,width=2.) if s.id in (3,4) else s for s in network.sections))
        local=backend.prepare(narrower)
        self.assertEqual(local.products[model.Key('bundle_plan',7)].mode,'separated')
        self.assertIs(local.products[model.Key('seam_surface',40)],original_product)
        self.assertIsNone(original_product.pavement.mesh)  # Committed query remains immutable.
        mesh=local.products[model.Key('seam_mesh',40)].mesh
        fresh=PythonMeshBackend(narrower.settings).prepare(narrower).products[model.Key('seam_mesh',40)].mesh
        np.testing.assert_array_equal(mesh.xyz,fresh.xyz)
        np.testing.assert_array_equal(mesh.quads,fresh.quads)


if __name__=='__main__': unittest.main(verbosity=2)
