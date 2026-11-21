"""A real short-prefix sector keeps its later native contact/profile separate."""
from copy import copy
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_owned_prefixes import fixture as base_network
from test_python_street_dimensional_bundles import reverse_sections
from tested_python_streets import model
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import UpdateCoordinator,BuildContext
from tested_python_streets.routes import Route,RouteTopology
from tested_python_streets.intersections import Approach,IntersectionProvider,Mouth
from tested_python_streets.owned_prefixes import approaches as prefix_approaches,contact_walks
from tested_python_streets.ownership import footprint_key
from tested_python_streets.compounds import CompoundGroup,sample_footprint
from tested_python_streets.regions import RegionUnion
from tested_python_streets.junction_bank_transitions import original_corner
from tested_python_streets.owned_prefix_corner_authority import (validate_owned_prefix_corner,
    verified_owned_prefix_ranges,OwnedPrefixCornerAuthorityError)


def case(*,reverse=False,angle=0.,translation=(0.,0.),mode='STRICT'):
    n=base_network()
    fifth=model.Section(5,5,((-30.,-3.,0.),(.25,-3.,0.)),width=5.4)
    bundle=model.Bundle(7,(model.BundleMember((model.Traversal(2,False),)),model.BundleMember((model.Traversal(5,False),))))
    n=replace(n,sections=(*n.sections,fifth),bundles=(bundle,),
        intersections=tuple(replace(node,radius_modes=(mode,)*len(node.ends)) for node in n.intersections))
    if reverse:n=reverse_sections(n,{1,2,3,4,5})
    rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
    n=replace(n,sections=tuple(replace(section,points=tuple(map(tuple,np.c_[np.asarray(section.points)[:,:2]@rotation.T+translation,np.zeros(len(section.points))]))) for section in n.sections))
    patterns=PatternLibrary();update=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(n)
    ctx=BuildContext(n,update.products,patterns,update.task_fingerprints,update.source_fingerprints)
    prefix=ctx[model.Key('path',1)]
    if not hasattr(prefix,'owned_prefix'):raise AssertionError('The actual Path producer did not select the short-prefix claim')
    node=n.intersections[0];topology=RouteTopology(n)
    original=[Approach(Route.build(ctx,topology.trace(end)),end.end,ctx[model.Key('profile',end.section_id)]) for end in node.ends]
    arms,claims=prefix_approaches(ctx,node,original)
    layout=IntersectionProvider.curved_boundary(ctx,node,arms);ctx.products[footprint_key(node.key)]=layout
    source=sample_footprint(node.key,layout)
    group=CompoundGroup((node.key,),(),(source,),RegionUnion().build((source.region,)),layout)
    corner=next(contact_walks(ctx,node.key,c) for c in layout.corners if c.center is not None and any(hasattr(a.path,'claim') for a in (c.first,c.second)))
    return ctx,group,node,corner


class OwnedPrefixCornerAuthorityTests(unittest.TestCase):
    def test_real_prefix_and_half_profile_contact_keep_distinct_logical_ends(self):
        for mode in ('STRICT','STYLE'):
            ctx,group,node,corner=case(mode=mode);original=ctx.network.to_dict()
            self.assertFalse(original_corner(type('Resolver',(),dict(context=ctx))(),node,corner))
            proof=validate_owned_prefix_corner(ctx,group,node,corner)
            self.assertIn(model.End(1,0),proof.logical_ends);self.assertIn(model.End(2,1),proof.contact_ends)
            self.assertNotIn(model.End(2,1),node.ends)
            self.assertTrue(all(p.final_mouth_width_proof['maximum_error']<=.002 for p in proof.prefix_range_proofs))
            if mode=='STRICT':self.assertEqual(corner.radius,node.corner_radius(proof.logical_ends[0])[0])
            self.assertEqual(ctx.network.to_dict(),original)

    def test_reversed_rotated_and_equal_reconstructed_inputs_retain_canonical_paths(self):
        for reverse,angle,translation in ((False,.713,(1200.,-1800.)),(True,0.,(0.,0.))):
            ctx,group,node,corner=case(reverse=reverse,angle=angle,translation=translation)
            self.assertEqual(validate_owned_prefix_corner(ctx,group,node,corner).junction,node.key)
            clone=model.Network.from_dict(ctx.network.to_dict())
            reconstructed=BuildContext(clone,ctx.products,ctx.patterns,ctx.task_digests,ctx.source_digests)
            proof=validate_owned_prefix_corner(reconstructed,group,clone.intersections[0],corner)
            self.assertEqual(proof.logical_ends,tuple(model.End(arm.source.id,arm.end) for arm in (corner.first,corner.second)))

    def test_public_parent_authority_and_owned_ranges_use_only_the_typed_proof(self):
        from tested_python_streets.displaced_corner_authority import original_parent_corner
        from tested_python_streets.source_bank_reconciliation import owned_ranges
        ctx,group,node,corner=case()
        resolver=type('Resolver',(),dict(context=ctx,group=group))()
        self.assertFalse(original_corner(resolver,node,corner))
        self.assertTrue(original_parent_corner(resolver,node,corner))
        actual=owned_ranges(ctx,group)
        arm=next(arm for arm in (corner.first,corner.second) if hasattr(arm.path,'claim'))
        mouth=next(at for candidate,at in zip(group.physical_layout.approaches,group.physical_layout.outward_stations)
                   if getattr(candidate.path,'claim',None) is arm.path.claim)
        expected=dict(verified_owned_prefix_ranges(ctx,group,arm,mouth).ranges)
        self.assertTrue(expected)
        for sid,ranges in expected.items():self.assertEqual(actual[sid],ranges)

    def test_copied_claim_and_borrowed_native_path_do_not_inherit_authority(self):
        ctx,group,node,corner=case();walks=corner.source_contact_walks
        index=next(i for i,w in enumerate(walks) if hasattr(w.approach.path,'claim'))
        for borrowed in (False,True):
            arm=copy(walks[index].approach);arm.path=copy(arm.path)
            if borrowed:
                pieces=list(arm.path.pieces);pieces[0]=replace(pieces[0],path=copy(pieces[0].path));arm.path.pieces=tuple(pieces)
            else:arm.path.claim=replace(arm.path.claim)
            changed=copy(corner);changed.first=arm if index==0 else corner.first;changed.second=arm if index==1 else corner.second
            changed.source_contact_walks=tuple(replace(walk,approach=arm) if i==index else walk for i,walk in enumerate(walks))
            with self.assertRaises(OwnedPrefixCornerAuthorityError):validate_owned_prefix_corner(ctx,group,node,changed)

    def test_missing_alias_and_contact_beyond_unchanged_mouth_are_refused(self):
        ctx,group,node,corner=case();index=next(i for i,a in enumerate((corner.first,corner.second)) if hasattr(a.path,'claim'))
        arm=(corner.first,corner.second)[index];claim=arm.path.claim
        saved=group.physical_layout.mouths.pop((claim.endpoint.section_id,claim.endpoint.end))
        try:
            with self.assertRaises(OwnedPrefixCornerAuthorityError):verified_owned_prefix_ranges(ctx,group,arm,corner.stations[index])
        finally:group.physical_layout.mouths[claim.endpoint.section_id,claim.endpoint.end]=saved
        mouth=next(s for a,s in zip(group.physical_layout.approaches,group.physical_layout.outward_stations)
                   if getattr(a.path,'claim',None) is claim)
        with self.assertRaises(OwnedPrefixCornerAuthorityError):verified_owned_prefix_ranges(ctx,group,arm,mouth+.001)

    def test_radius_and_actual_contact_geometry_cannot_be_relabelled(self):
        ctx,group,node,corner=case()
        for changes in (dict(requested_radius=corner.requested_radius+.1),dict(radius_mode='STYLE'),
                        dict(radius_attempts=(corner.radius,.5)),dict(center=corner.center+(.001,0.)),
                        dict(gap=corner.gap+.01),dict(sweep=-corner.sweep)):
            changed=copy(corner)
            for key,value in changes.items():setattr(changed,key,value)
            with self.assertRaises(OwnedPrefixCornerAuthorityError):validate_owned_prefix_corner(ctx,group,node,changed)


if __name__=='__main__':unittest.main()
