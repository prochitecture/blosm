"""A short forced-tangent failure can be owned by its real incident junction."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_dimensional_bundles import reverse_sections
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.mesh import Mesh, publication_coordinates
from tested_python_streets.ownership import surface_key
from tested_python_streets.providers import RoadSpan
from tested_python_streets.validation import validate_overlaps


def fixture():
    sections=(model.Section(1,1,((0.,0.,0.),(.25,.70,0.)),width=10.8),
        model.Section(2,2,((-30.,.70,0.),(.25,.70,0.)),width=5.4),
        model.Section(3,3,((0.,0.,0.),(30.,0.,0.)),width=10.8),
        model.Section(4,4,((0.,0.,0.),(0.,30.,0.)),width=10.8))
    return model.Network('owned_prefix_T',sections,
        (model.Intersection(1,(model.End(1,0),model.End(3,0),model.End(4,0))),),
        connections=(model.Connection(1,(model.End(1,1),model.End(2,1))),))


class OwnedPrefixTests(unittest.TestCase):
    def qualify(self,n):
        original=n.to_dict()
        update=PythonMeshBackend(n.settings).prepare(n)
        self.assertEqual(original,n.to_dict())
        path=update.products[model.Key('path',1)]
        self.assertTrue(path._certified())
        self.assertEqual(path.owned_prefix.endpoint.section_id,1)
        np.testing.assert_array_equal(path.source_points,np.asarray(n.section_index[1].points)[:,:2])
        self.assertIsInstance(update.products[model.Key('road',1)],RoadSpan)
        self.assertEqual(update.products[model.Key('seam_layout',1)].absorbed_by,model.Key('intersection',1))
        self.assertIsNone(update.products[model.Key('seam_mesh',1)])
        meshes={key:Mesh(publication_coordinates(p.mesh.xyz,1.,physical=True),p.mesh.quads)
                for key,p in update.products.items() if getattr(p,'mesh',None) is not None}
        self.assertTrue(meshes)
        for mesh in meshes.values():self.assertLessEqual(mesh.validate()['maximum_edge_ratio'],8.)
        validate_overlaps(meshes,tolerance=2e-5)
        self.assertEqual(set(meshes),{model.Key('junction_mesh',1),*(model.Key('road',sid) for sid in (2,3,4))})
        actual=update.products[model.Key('junction_mesh',1)].mesh.xyz
        layout=update.products[model.Key('mouths',1)]
        claim=path.owned_prefix
        exterior=update.products[model.Key('road',2)].mesh
        end=claim.external.end
        row=exterior.xyz[exterior.seams[end]]
        self.assertLess(float(np.linalg.norm(row[:,None]-actual[None],axis=2).min(axis=1).max()),1e-8)
        self.assertIn((2,end),layout.mouths)
        return update

    def test_real_external_frame_and_profile_replace_impossible_short_seam(self):
        self.qualify(fixture())

    def test_authored_directions_and_source_street_grouping_are_preserved(self):
        for reverse in ({1},{2},{1,2,3,4}):
            with self.subTest(reverse=reverse):self.qualify(reverse_sections(fixture(),reverse))
        n=fixture();n=replace(n,sections=tuple(replace(s,street_id=1) if s.id==2 else s for s in n.sections))
        self.qualify(n)

    def test_no_real_junction_or_known_scope_boundary_cannot_defer(self):
        n=replace(fixture(),intersections=())
        with self.assertRaisesRegex(ValueError,'Smooth'):
            PythonMeshBackend(n.settings).prepare(n)
        with self.assertRaisesRegex(ValueError,'scope boundary'):
            replace(fixture(),sections=tuple(replace(s,scope_ends=(0,)) if s.id==1 else s
                                             for s in fixture().sections))

    def test_uncertified_visible_road_cannot_use_prefix_as_a_general_escape(self):
        n=fixture()
        engine=PythonMeshBackend(n.settings);accepted=engine.prepare(n);engine.commit(accepted)
        manifest=engine.manifest()
        # Remove the owner after a successful build: cached prefix approval
        # must not survive a topology change to an ordinary short bent seam.
        with self.assertRaisesRegex(ValueError,'Smooth'):
            engine.prepare(replace(n,intersections=()))
        self.assertEqual(engine.manifest(),manifest)

    def test_local_update_restoration_and_consumed_source_anchor_use_physical_surface(self):
        n=replace(fixture(),furniture=(model.Placement(1,1,.3,offset=-5.4),))
        engine=PythonMeshBackend(n.settings);initial=engine.prepare(n);engine.commit(initial)
        changed=replace(n,sections=tuple(replace(s,width=5.7) if s.id==2 else s for s in n.sections))
        local=engine.prepare(changed);fresh=PythonMeshBackend(changed.settings).prepare(changed)
        for key,part in local.products.items():
            if getattr(part,'mesh',None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz,fresh.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads,fresh.products[key].mesh.quads)
        anchor=local.products[model.Key('placement',1)]
        product=local.products[surface_key(model.Key('intersection',1))]
        self.assertAlmostEqual(anchor.position[2],product.height(np.asarray([anchor.position[:2]]))[0])
        engine.commit(local)
        restored=engine.prepare(n)
        np.testing.assert_array_equal(restored.products[model.Key('junction_mesh',1)].mesh.xyz,
                                      initial.products[model.Key('junction_mesh',1)].mesh.xyz)


if __name__=='__main__':unittest.main()
