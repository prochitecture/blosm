"""Junction ownership of near dimensional seams, profile ports and Updates."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.ownership import footprint_key, surface_key
from tested_python_streets.providers import FurnitureProvider
from _python_street_fixtures import absorbed_transition_fixture_networks
from test_python_street_dimensional_bundles import meshes, reverse_sections


class PythonJunctionTransitionTests(unittest.TestCase):
    def test_extra_width_certification_does_not_change_cached_source_sampling(self):
        from tested_python_streets.paths import Path
        path=Path(model.Section(1,1,((0.,0.,0.),(10.,4.,0.),(15.,-2.,0.),(30.,5.,0.)),width=.1))
        saved=path.certificate_intervals,path.curvature_bound
        samples=path.samples(0.,path.length,.5)
        # This wider reach needs a finer proof subdivision than the original
        # source. It must not overwrite the original path's sampling bound.
        self.assertTrue(path._certified(.98/saved[1]*1.15))
        self.assertEqual((path.certificate_intervals,path.curvature_bound),saved)
        np.testing.assert_array_equal(path.samples(0.,path.length,.5),samples)
        self.assertFalse(path._certified(1e6))
        self.assertEqual((path.certificate_intervals,path.curvature_bound),saved)

    def fixture(self):
        return absorbed_transition_fixture_networks(model)[0]

    def qualify(self, network):
        saved = network.to_dict()
        result = PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(saved, network.to_dict())
        layout = result.products[footprint_key(model.Key('intersection', 1))]
        self.assertEqual(layout.absorbed_transitions, (model.Key('connection', 6),))
        self.assertEqual(result.products[model.Key('seam_layout', 6)].absorbed_by, model.Key('intersection', 1))
        self.assertIsNone(result.products[model.Key('seam_surface', 6)])
        self.assertIsNone(result.products[model.Key('seam_mesh', 6)])
        self.assertIsNone(getattr(result.products[model.Key('road', 1)], 'mesh', None))
        pavement = result.products[model.Key('junction_mesh', 1)].mesh
        for part in result.products.values():
            if getattr(part, 'mesh', None) is not None:
                self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
                self.assertTrue(np.all(part.mesh.xyz[:, 2] >= network.settings.exaggeration * 0. - 1e-8))
                self.assertLessEqual(float(part.mesh.xyz[:, 2].max()), .08 * network.settings.exaggeration + 1e-8)
        # The real external port belongs to Section 2, including its complete
        # width/profile and original path normal, regardless of root direction.
        road = result.products[model.Key('road', 2)]
        port_end = next(end for (sid, end) in layout.mouths if sid == 2)
        profile = result.products[model.Key('profile', 2)]
        row = road.mesh.xyz[road.mesh.seams[port_end]]
        q = np.linspace(-1., 1., len(row))
        np.testing.assert_allclose(row[:, 2], profile.datum + network.settings.exaggeration * profile.value(q), atol=1e-12)
        self.assertAlmostEqual(float(np.linalg.norm(row[-1, :2] - row[0, :2])), network.section_index[2].width, places=9)
        distance = np.linalg.norm(row[:, None] - pavement.xyz[None], axis=-1)
        self.assertLess(float(distance.min(axis=1).max()), 1e-8)
        context = type('Context', (), {'__getitem__':lambda self,k: result.products[k]})()
        for p in network.furniture[:2]:
            self.assertEqual(FurnitureProvider.compound_owner(context,p), model.Key('intersection',1))
            transform=result.products[model.Key('placement',p.id)]
            expected=result.products[surface_key(model.Key('intersection',1))].height(np.asarray([transform.position[:2]]))[0]
            self.assertAlmostEqual(transform.position[2], expected+p.height, places=10)
        return result

    def test_complete_profile_ports_and_furniture_in_both_directions_and_modes(self):
        for mode in (0,3):
            for reversed_ids in (set(), {1}, {2}, {1,2}):
                with self.subTest(mode=mode, reverse=reversed_ids):
                    n=self.fixture()
                    n=replace(n, sections=tuple(replace(s,interpolation=mode) for s in n.sections))
                    self.qualify(reverse_sections(n,reversed_ids))

    def test_resolution_and_rigid_transform_preserve_ownership_and_profiles(self):
        for resolution in (1,2,4):
            with self.subTest(resolution=resolution):
                n=self.fixture()
                angle=.731
                rot=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
                sources=[]
                for s in n.sections:
                    p=np.asarray(s.points).copy();p[:,:2]=p[:,:2]@rot.T+(1234.,-5678.)
                    sources.append(replace(s,points=tuple(map(tuple,p))))
                self.qualify(replace(n,sections=tuple(sources),settings=replace(n.settings,resolution=resolution)))

    def test_local_edits_and_ownership_reversal_match_full_after_buffer_release(self):
        n=self.fixture();backend=PythonMeshBackend(n.settings)
        first=backend.prepare(n);saved=meshes(first);digests=first.task_fingerprints.copy();backend.commit(first)
        remote=first.products[model.Key('road',99)]
        backend.release_published_meshes()
        changed=replace(n,sections=tuple(replace(s,width=6.1) if s.id==2 else s for s in n.sections))
        update=backend.prepare(changed);self.assertIs(update.products[model.Key('road',99)],remote)
        full=PythonMeshBackend(changed.settings).prepare(changed)
        for key in update.rebuilt:
            a,b=update.products[key],full.products[key]
            if getattr(a,'mesh',None) is not None:
                np.testing.assert_array_equal(a.mesh.xyz,b.mesh.xyz)
                np.testing.assert_array_equal(a.mesh.quads,b.mesh.quads)
        backend.commit(update);backend.release_published_meshes()
        apart=replace(n,sections=tuple(replace(s,points=((0.,0.,0.),(6.,0.,0.),(12.,0.,0.))) if s.id==1
            else replace(s,points=((12.,0.,0.),*s.points[1:])) if s.id==2 else s for s in n.sections))
        separate=backend.prepare(apart)
        self.assertIsNotNone(separate.products[model.Key('seam_mesh',6)])
        self.assertIsNone(separate.products[model.Key('seam_layout',6)].absorbed_by)
        self.assertIs(separate.products[model.Key('road',99)],remote)
        backend.commit(separate);backend.release_published_meshes()
        restored=backend.prepare(n)
        for key, buffers in saved.items():
            if key==model.Key('road',99):continue
            part=restored.products[key]
            if part.mesh is None:
                self.assertNotIn(key,restored.rebuilt)
                self.assertEqual(restored.task_fingerprints[key],digests[key])
            else:
                np.testing.assert_array_equal(part.mesh.xyz,buffers[0])
                np.testing.assert_array_equal(part.mesh.quads,buffers[1])
        for p in n.furniture:
            self.assertEqual(first.products[model.Key('placement',p.id)],restored.products[model.Key('placement',p.id)])

    def test_datum_street_or_mode_changes_cannot_be_hidden_and_failed_update_keeps_manifest(self):
        n=self.fixture();backend=PythonMeshBackend(n.settings);backend.commit(backend.prepare(n))
        manifest=backend.manifest()
        for fields in (dict(datum=1.),dict(street_id=2),dict(interpolation=0)):
            with self.subTest(fields=fields):
                changed=replace(n,sections=tuple(replace(s,**fields) if s.id==2 else s for s in n.sections))
                with self.assertRaises(ValueError):backend.prepare(changed)
                self.assertEqual(backend.manifest(),manifest)


if __name__=='__main__':unittest.main()
