"""Bundle ownership permissions invalidate cached plans without geometry edits."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import _python_street_package
from tested_python_streets import model
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.ownership import PLAN,JunctionOwnership,footprint_key,surface_key
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.updates import UpdateCoordinator,BuildContext


class OwnershipOnly(IntersectionProvider):
    def tasks(self,network):
        # Test preparation permission and caching without publishing an
        # independently unqualified Bundle pavement or separator surface.
        for task in super().tasks(network):
            if task.key.kind in ('approach','physical_intersection','compound_plan'):
                yield task


def fixture():
    sections=(model.Section(1,1,((0.,-20.,0.),(0.,0.,0.)),width=4.),
        model.Section(2,1,((0.,0.,0.),(0.,.1,0.)),width=4.),
        model.Section(3,1,((0.,.1,0.),(0.,20.1,0.)),width=4.),
        model.Section(4,2,((10.,-20.,0.),(10.,22.,0.)),width=4.),
        model.Section(5,3,((0.,0.,0.),(-20.,-10.,0.)),width=4.),
        model.Section(6,4,((0.,.1,0.),(-20.,10.1,0.)),width=4.))
    sites=(model.Intersection(10,(model.End(1,1),model.End(2,0),model.End(5,0)),radius_modes=('STYLE',)*3),
        model.Intersection(11,(model.End(2,1),model.End(3,0),model.End(6,0)),radius_modes=('STYLE',)*3))
    bundle=model.Bundle(1,(model.BundleMember(tuple(model.Traversal(i) for i in (1,2,3))),
        model.BundleMember((model.Traversal(4),))))
    return model.Network('partial-bundle-cache',sections,sites,(bundle,),settings=model.Settings(exaggeration=1.))


class PartialBundleCacheTests(unittest.TestCase):
    def coordinator(self):
        return UpdateCoordinator((PathProvider(),ProfileProvider(),OwnershipOnly()),PatternLibrary())

    def test_permission_edit_rejects_cached_plan_just_as_fresh_preparation(self):
        n=fixture();engine=self.coordinator()
        initial=engine.prepare(n)
        self.assertTrue(initial.products[PLAN].groups)
        self.assertTrue(initial.products[PLAN].groups[0].physical_layout.partial_bundle_claims)
        for bundle in n.bundles:self.assertIn(bundle.key,initial.dependencies[PLAN])
        engine.commit(initial);manifest=engine.manifest()
        b=n.bundles[0]
        path=b.members[0].path
        reordered=replace(n,bundles=(replace(b,members=(replace(b.members[0],
            path=(path[1],path[0],path[2])),b.members[1])),))
        # The old half/rise geometry is exactly unchanged. Only the original
        # member ordering (and hence its interval-ownership permission) moved.
        profile_update=UpdateCoordinator((PathProvider(),ProfileProvider()),PatternLibrary()).prepare(reordered)
        for key in initial.task_fingerprints:
            if key.kind=='profile':
                self.assertEqual(initial.task_fingerprints[key],profile_update.task_fingerprints[key])
        for coordinator in (engine,self.coordinator()):
            with self.assertRaisesRegex(ValueError,'Partial Bundle'):
                coordinator.prepare(reordered)
        self.assertEqual(engine.manifest(),manifest)

    def test_surface_and_mouth_signatures_follow_only_relevant_bundle_semantics(self):
        n=fixture();result=self.coordinator().prepare(n)
        context=BuildContext(n,result.products,PatternLibrary(),result.task_fingerprints,result.source_fingerprints)
        site=n.intersections[0];owner=result.products[PLAN].owner(site.key)
        before_surface=JunctionOwnership.surface_signature(context,owner)
        before_layout=JunctionOwnership.layout_signature(context,site)
        digests=dict(result.source_fingerprints);digests[n.bundles[0].key]='changed bundle permission'
        changed=BuildContext(n,result.products,PatternLibrary(),result.task_fingerprints,digests)
        self.assertNotEqual(before_surface,JunctionOwnership.surface_signature(changed,owner))
        self.assertNotEqual(before_layout,JunctionOwnership.layout_signature(changed,site))
        unrelated=model.Bundle(99,(model.BundleMember((model.Traversal(99),)),
            model.BundleMember((model.Traversal(100),))))
        extra=(model.Section(99,99,((100.,0.,0.),(100.,20.,0.))),
               model.Section(100,100,((110.,0.,0.),(110.,20.,0.))))
        extended=replace(n,sections=(*n.sections,*extra),bundles=(*n.bundles,unrelated))
        digests=dict(result.source_fingerprints);digests[unrelated.key]='unrelated bundle'
        remote=BuildContext(extended,result.products,PatternLibrary(),result.task_fingerprints,digests)
        self.assertEqual(before_surface,JunctionOwnership.surface_signature(remote,owner))
        self.assertEqual(before_layout,JunctionOwnership.layout_signature(remote,site))

    def test_joint_source_phase_invalidates_surface_and_mouth_without_local_geometry_edit(self):
        from tested_python_streets.bundle_merged_entrances import preparation_key
        n=fixture();result=self.coordinator().prepare(n)
        site=n.intersections[0];owner=result.products[PLAN].owner(site.key)
        phase=preparation_key(n.bundles[0])
        products=dict(result.products)
        products[phase]=SimpleNamespace(source_phase_proof={'qualified':True})
        digests=dict(result.task_fingerprints);digests[phase]='original shared source partition'
        context=BuildContext(n,products,PatternLibrary(),digests,result.source_fingerprints)
        before_surface=JunctionOwnership.surface_signature(context,owner)
        before_layout=JunctionOwnership.layout_signature(context,site)
        changed=dict(digests);changed[phase]='opposite member changes shared source partition'
        edited=BuildContext(n,products,PatternLibrary(),changed,result.source_fingerprints)
        self.assertNotEqual(before_surface,JunctionOwnership.surface_signature(edited,owner))
        self.assertNotEqual(before_layout,JunctionOwnership.layout_signature(edited,site))
        # Raw local physical tasks and Bundle membership did not change.
        self.assertEqual(result.task_fingerprints[footprint_key(site.key)],digests[footprint_key(site.key)])
        self.assertEqual(context.digest(n.bundles[0].key),edited.digest(n.bundles[0].key))
        # A prepared resource for an unrelated Bundle cannot invalidate it.
        unrelated=preparation_key(SimpleNamespace(id=99))
        products=dict(products);products[unrelated]=SimpleNamespace(source_phase_proof={'qualified':True})
        changed=dict(digests);changed[unrelated]='remote source partition'
        remote=BuildContext(n,products,PatternLibrary(),changed,result.source_fingerprints)
        self.assertEqual(before_surface,JunctionOwnership.surface_signature(remote,owner))
        self.assertEqual(before_layout,JunctionOwnership.layout_signature(remote,site))


if __name__=='__main__':unittest.main()
