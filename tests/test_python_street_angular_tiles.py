"""Physical island interfaces, real edit recovery and high-resolution grading."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_angular_fixtures import angular_dimensional_network
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_angular_curves import polyline_distance
from test_python_street_angular_bundles import buffers


def source_network(sweep,mode,*,float32=False,resolution=1):
    source=angular_dimensional_network(model,sweep,True,mode)
    if float32:
        source=replace(source,sections=tuple(replace(s,points=tuple(map(tuple,
            np.asarray(s.points,dtype=np.float32).astype(float)))) for s in source.sections))
    return replace(source,settings=replace(source.settings,resolution=resolution))


class AngularPhysicalTileTests(unittest.TestCase):
    @staticmethod
    def equal(first,second):
        if first.keys()!=second.keys(): raise AssertionError('Different output keys')
        for key in first:
            np.testing.assert_array_equal(first[key][0],second[key][0])
            np.testing.assert_array_equal(first[key][1],second[key][1])

    def qualify(self,source,update):
        plan=update.products[model.Key('bundle_plan',13)]
        self.assertEqual(type(plan).__name__,'AngularBundlePlan')
        self.assertEqual(len(plan.transition_owners),4)
        self.assertEqual(len(plan.islands),1)
        pavement=update.products[model.Key('bundle_pavement',13)]
        grass=update.products[model.Key('separator',13)]
        for part in (pavement,grass):
            self.assertEqual(part.mesh.quads.shape[1],4)
            self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'],8.)
        self.assertLessEqual(pavement.metadata['angular_longitudinal_bound'],.03+1e-8)
        self.assertLessEqual(pavement.metadata['physical_slopes']['maximum'],plan.maximum_slope+1e-8)
        self.assertGreaterEqual(float(pavement.mesh.xyz[:,2].min()),-1e-10)
        self.assertLessEqual(float(pavement.mesh.xyz[:,2].max()),.08*source.settings.exaggeration+1e-8)
        # Audit the actual published interfaces independently of patch keys.
        def edges(mesh):
            all_edges=np.sort(np.stack((mesh.quads,np.roll(mesh.quads,-1,axis=1)),axis=-1).reshape(-1,2),axis=1)
            unique,count=np.unique(all_edges,axis=0,return_counts=True)
            self.assertLessEqual(int(count.max()),2)
            return unique[count==1]
        boundary=edges(grass.mesh)
        road_boundary=edges(pavement.mesh)
        road_ids=np.unique(road_boundary)
        mapped={}
        for index in np.unique(boundary):
            point=grass.mesh.xyz[index,:2]
            distance=np.linalg.norm(pavement.mesh.xyz[road_ids,:2]-point,axis=1)
            at=int(distance.argmin())
            self.assertLessEqual(float(distance[at]),1e-8)
            mapped[int(index)]=int(road_ids[at])
        available={tuple(pair) for pair in road_boundary}
        for a,b in boundary:
            self.assertIn(tuple(sorted((mapped[int(a)],mapped[int(b)]))),available)
        # Every sampled physical grass edge stays on the original source curbs
        # to the retained curve budget; no interior computational cuts appear.
        points=grass.mesh.xyz[boundary,:2]
        query=np.concatenate([(1.-t)*points[:,0]+t*points[:,1] for t in np.linspace(0.,1.,9)])
        error=np.min(np.stack([polyline_distance(query,curve.points) for curve in plan.field.curves]),axis=0)
        self.assertLessEqual(float(error.max()),source.settings.curve_tolerance+1e-8)
        return plan

    def test_exact_islands_all_sweeps_float64_and_blender_float32(self):
        for sweep in (180.,225.,270.):
            for precision in (False,True):
                with self.subTest(sweep=sweep,float32=precision):
                    source=source_network(sweep,0,float32=precision)
                    before=source.to_dict()
                    update=PythonMeshBackend(source.settings).prepare(source)
                    self.qualify(source,update)
                    self.assertEqual(source.to_dict(),before)
                    total=sum(len(update.products[key].mesh.quads) for key in
                              (model.Key('bundle_pavement',13),model.Key('separator',13)))
                    self.assertLess(total,10000)

    def test_native_first_member_width_edit_equal_full_restore_and_keep_remote(self):
        source=source_network(180.,3,float32=True)
        backend=PythonMeshBackend(source.settings)
        initial=backend.prepare(source)
        backend.commit(initial)
        saved=buffers(initial)
        remote=initial.products[model.Key('road',99)]
        changed=replace(source,sections=tuple(replace(s,width=float(np.float32(s.width+.1)))
                        if s.id in (17,117,217) else s for s in source.sections))
        local=backend.prepare(changed)
        self.qualify(changed,local)
        full=PythonMeshBackend(changed.settings).prepare(changed)
        self.equal(buffers(local),buffers(full))
        self.assertIs(local.products[model.Key('road',99)],remote)
        self.assertNotIn(model.Key('placement',999),local.rebuilt)
        for placement in changed.furniture:
            self.assertEqual(local.products[model.Key('placement',placement.id)],
                             full.products[model.Key('placement',placement.id)])
        backend.commit(local)
        restored=backend.prepare(source)
        self.equal(saved,buffers(restored))

    def test_actual_float32_r4_island_keeps_world_ratio_crown_and_slope(self):
        source=source_network(225.,3,float32=True,resolution=4)
        before=source.to_dict()
        update=PythonMeshBackend(source.settings).prepare(source)
        self.qualify(source,update)
        self.assertEqual(source.to_dict(),before)

    def test_actual_float32_exact_180_r4_separator_keeps_complete_interfaces(self):
        from tested_python_streets.mesh import Mesh
        source=source_network(180.,0,float32=True,resolution=4)
        source=replace(source,sections=tuple(replace(s,rise=float(np.float32(s.rise))) for s in source.sections))
        before=source.to_dict()
        update=PythonMeshBackend(source.settings).prepare(source)
        self.qualify(source,update)
        grass=update.products[model.Key('separator',13)]
        published=Mesh(grass.mesh.xyz.astype(np.float32),grass.mesh.quads)
        self.assertLessEqual(published.validate()['maximum_edge_ratio'],8.)
        self.assertTrue(np.all(grass.mesh.xyz[:,2]==grass.mesh.xyz[0,2]))
        correction=grass.metadata.get('refinement_correction')
        if correction is not None:
            self.assertTrue(correction['exact_fixed_xyz'])
            self.assertEqual(correction['physical_slopes']['maximum'],0.)
        self.assertEqual(source.to_dict(),before)

    def test_exact_270_float32_control_move_equal_full_and_restores(self):
        source=source_network(270.,0,float32=True)
        source=replace(source,sections=tuple(replace(s,rise=float(np.float32(s.rise))) for s in source.sections))
        backend=PythonMeshBackend(source.settings)
        initial=backend.prepare(source);backend.commit(initial)
        saved=buffers(initial);sections=[]
        for section in source.sections:
            if section.id==229:
                points=list(section.points);index=len(points)//2
                x,y,z=points[index];points[index]=(float(np.float32(x+.1)),y,z)
                section=replace(section,points=tuple(points))
            sections.append(section)
        changed=replace(source,sections=tuple(sections))
        local=backend.prepare(changed);self.qualify(changed,local)
        full=PythonMeshBackend(changed.settings).prepare(changed)
        self.equal(buffers(local),buffers(full))
        self.assertIs(local.products[model.Key('road',99)],initial.products[model.Key('road',99)])
        for placement in changed.furniture:
            self.assertEqual(local.products[model.Key('placement',placement.id)],
                             full.products[model.Key('placement',placement.id)])
        backend.commit(local)
        self.equal(saved,buffers(backend.prepare(source)))


if __name__=='__main__': unittest.main()
