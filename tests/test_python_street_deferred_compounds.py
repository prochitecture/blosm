"""Short connected junctions can own one rim before isolated mouths exist."""
from dataclasses import replace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from _python_street_fixtures import short_junction_pair
from test_python_street_dimensional_bundles import reverse_sections
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.compound_preparation import DeferredCompoundPlanner, DeferredCompoundFootprint
from tested_python_streets.ownership import PLAN, footprint_key
from tested_python_streets.mesh import Mesh, publication_coordinates


def fixture(mode=3):
    n=short_junction_pair(model,length=2.527)
    return replace(n,sections=tuple(replace(s,interpolation=mode) for s in n.sections))


class DeferredCompoundTests(unittest.TestCase):
    def qualify(self,n):
        original=n.to_dict()
        update=PythonMeshBackend(n.settings).prepare(n)
        for product in update.products.values():
            mesh=getattr(product,'mesh',None)
            if mesh is not None:
                stored=Mesh(publication_coordinates(mesh.xyz,1.,physical=True),mesh.quads)
                self.assertLessEqual(stored.validate()['maximum_edge_ratio'],8.)
        self.assertEqual(original,n.to_dict())
        root=model.Key('intersection',1)
        group=update.products[PLAN].by_site[root]
        self.assertEqual(group.members,(root,model.Key('intersection',2)))
        self.assertIsNotNone(group.physical_layout)
        self.assertEqual({c.section_id for c in group.claims},{1})
        self.assertTrue(any(isinstance(update.products[footprint_key(site)],DeferredCompoundFootprint)
                            for site in group.members))
        self.assertIsNone(update.products[model.Key('junction_mesh',2)])
        self.assertIsNone(getattr(update.products[model.Key('road',1)],'mesh',None))
        part=update.products[model.Key('junction_mesh',1)]
        self.assertEqual(part.owner,model.Key('compound',1))
        self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'],8.)
        physical=group.physical_layout
        self.assertGreater(max(np.linalg.norm(a.origin-physical.approaches[0].origin)
                               for a in physical.approaches),2.)
        self.assertEqual(set(physical.mouths),{(1,0),(1,1),(2,0),(3,0),(4,0),(5,0)})
        xyz=part.mesh.xyz
        self.assertGreaterEqual(float(xyz[:,2].min()),-1e-8)
        self.assertLessEqual(float(xyz[:,2].max()),n.settings.exaggeration*max(s.rise for s in n.sections)+1e-8)
        for sid in (2,3,4,5):
            road=update.products[model.Key('road',sid)].mesh
            seam=road.xyz[road.seams[0]]
            self.assertLess(float(np.linalg.norm(seam[:,None]-xyz[None],axis=2).min(axis=1).max()),1e-8)
        return update

    def test_both_interpolation_modes_keep_full_profiles_and_original_source_origins(self):
        for mode in (0,3):
            with self.subTest(mode=mode):self.qualify(fixture(mode))

    def test_rotated_translated_unequal_streets_preserve_requested_radii(self):
        n=fixture()
        angle=.731
        matrix=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
        streets=[]
        for s in n.sections:
            p=np.asarray(s.points).copy();p[:,:2]=p[:,:2]@matrix.T+(245.,-119.)
            streets.append(replace(s,points=tuple(map(tuple,p)),width=5.4 if s.id in (1,2) else 5.))
        n=replace(n,sections=tuple(streets))
        result=self.qualify(n)
        physical=result.products[PLAN].groups[0].physical_layout
        self.assertTrue(all(c.radius==2. for c in physical.corners if c.center is not None))

    def test_original_directions_do_not_change_physical_ownership(self):
        n=reverse_sections(fixture(),{1,2,4})
        result=PythonMeshBackend(n.settings).prepare(n)
        group=result.products[PLAN].groups[0]
        self.assertEqual({c.section_id for c in group.claims},{1})
        self.assertEqual(set(group.physical_layout.mouths),{(1,0),(1,1),(2,1),(3,0),(4,1),(5,0)})
        self.assertIsNone(getattr(result.products[model.Key('road',1)],'mesh',None))

    def test_edit_and_restoration_match_fresh_builds_without_rebuilding_remote_meshes(self):
        n=replace(fixture(),sections=(*fixture().sections,
            model.Section(99,99,((100.,100.,0.),(140.,100.,0.)))))
        engine=PythonMeshBackend(n.settings);initial=engine.prepare(n);engine.commit(initial)
        changed=replace(n,sections=tuple(replace(s,width=6.3) if s.id==3 else s for s in n.sections))
        local=engine.prepare(changed);fresh=PythonMeshBackend(changed.settings).prepare(changed)
        self.assertIs(local.products[model.Key('road',99)],initial.products[model.Key('road',99)])
        for key,part in local.products.items():
            if getattr(part,'mesh',None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz,fresh.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads,fresh.products[key].mesh.quads)
        engine.commit(local)
        restored=engine.prepare(n)
        np.testing.assert_array_equal(restored.products[model.Key('junction_mesh',1)].mesh.xyz,
                                      initial.products[model.Key('junction_mesh',1)].mesh.xyz)

    def test_connected_member_budget_is_a_hard_guard(self):
        with patch.object(DeferredCompoundPlanner,'maximum_members',1):
            with self.assertRaisesRegex(ValueError,'member budget'):
                PythonMeshBackend(fixture().settings).prepare(fixture())

    def test_unresolvable_exterior_contacts_are_not_published(self):
        n=fixture()
        n=replace(n,sections=tuple(replace(s,points=((2.527,0.,0.),(2.527,30.,0.)))
                                   if s.id==4 else s for s in n.sections))
        with self.assertRaisesRegex(ValueError,'indistinguishable outward directions'):
            PythonMeshBackend(n.settings).prepare(n)


if __name__=='__main__':unittest.main()
