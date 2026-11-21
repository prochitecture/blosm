"""A compound may consume an interior Bundle interval without losing its cuts."""
from dataclasses import replace
import unittest
import _python_street_package
from tested_python_streets.model import (Section,End,Intersection,Connection,Bundle,BundleMember,
    Traversal,Network,Key)
from tested_python_streets.partial_bundle_claims import prepare
from tested_python_streets.providers import PathProvider
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.updates import BuildContext,UpdateCoordinator
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.routes import Route,RouteTopology
from tested_python_streets.intersections import Approach


def sources(reverse=False):
    sections=(Section(1,1,((0.,0.,0.),(0.,10.,0.)),width=4.),
        Section(2,1,((0.,10.,0.),(0.,20.,0.)),width=4.),
        Section(3,1,((0.,20.,0.),(0.,30.,0.)),width=4.),
        Section(4,2,((10.,0.,0.),(10.,30.,0.)),width=4.),
        Section(5,3,((0.,10.,0.),(-20.,10.,0.)),width=4.),
        Section(6,4,((0.,20.,0.),(-20.,20.,0.)),width=4.))
    sites=(Intersection(10,(End(1,1),End(2,0),End(5,0))),
           Intersection(11,(End(2,1),End(3,0),End(6,0))))
    visits=tuple(Traversal(sid,True) for sid in (1,2,3))
    if reverse:
        sections=tuple(replace(s,points=s.points[::-1]) if s.id in (1,2,3) else s for s in sections)
        sites=tuple(replace(j,ends=tuple(End(e.section_id,1-e.end) if e.section_id in (1,2,3) else e for e in j.ends)) for j in sites)
        visits=tuple(Traversal(sid,False) for sid in (1,2,3))
    bundle=Bundle(1,(BundleMember(visits),BundleMember((Traversal(4,True),))))
    network=Network('partial_original_bundle',sections,sites,(bundle,))
    products=UpdateCoordinator((PathProvider(),ProfileProvider()),PatternLibrary()).prepare(network).products
    context=BuildContext(network,products,PatternLibrary())
    cuts=(End(1,0 if reverse else 1),End(3,1 if reverse else 0))
    approaches=[Approach(Route.build(context,RouteTopology(network).trace(e)),e.end,
        products[Key('profile',e.section_id)]) for e in cuts]
    return context,bundle,approaches,{j.key for j in sites}


def paired_ends(whole=False):
    sections=[Section(sid,street,((x,float(y),0.),(x,float(y+10),0.)),width=4.)
        for street,x,ids in ((1,0.,(1,2,3)),(2,10.,(4,5,6)))
        for sid,y in zip(ids,(0,10,20))]
    sections.extend((Section(7,3,((0.,20.,0.),(-20.,20.,0.)),width=4.),
        Section(8,4,((10.,20.,0.),(30.,20.,0.)),width=4.),
        Section(9,5,((0.,30.,0.),(10.,30.,0.)),width=4.),
        Section(10,6,((0.,30.,0.),(0.,50.,0.)),width=4.),
        Section(11,7,((10.,30.,0.),(30.,30.,0.)),width=4.)))
    sites=[Intersection(10,(End(2,1),End(3,0),End(7,0))),
        Intersection(11,(End(5,1),End(6,0),End(8,0))),
        Intersection(12,(End(3,1),End(9,0),End(10,0))),
        Intersection(13,(End(6,1),End(9,1),End(11,0)))]
    if whole:
        sections.extend((Section(12,8,((0.,0.,0.),(0.,-20.,0.)),width=4.),
            Section(13,9,((10.,0.,0.),(30.,0.,0.)),width=4.),
            Section(14,10,((0.,0.,0.),(10.,0.,0.)),width=4.)))
        sites.extend((Intersection(14,(End(1,0),End(12,0),End(14,0))),
            Intersection(15,(End(4,0),End(13,0),End(14,1)))))
    bundle=Bundle(1,tuple(BundleMember(tuple(Traversal(sid,True) for sid in ids))
        for ids in ((1,2,3),(4,5,6))))
    network=Network('paired_original_ends',tuple(sections),tuple(sites),(bundle,),
        connections=(Connection(0,(End(1,1),End(2,0))),Connection(1,(End(4,1),End(5,0)))))
    products=UpdateCoordinator((PathProvider(),ProfileProvider()),PatternLibrary()).prepare(network).products
    context=BuildContext(network,products,PatternLibrary())
    cuts=() if whole else (End(2,1),End(5,1))
    approaches=[Approach(Route.build(context,RouteTopology(network).trace(e)),e.end,
        products[Key('profile',e.section_id)]) for e in cuts]
    return context,bundle,approaches,{j.key for j in sites}


class PartialBundleClaimTests(unittest.TestCase):
    def test_both_terminal_intervals_keep_original_junction_owners_and_cut_profiles(self):
        context,bundle,approaches,members=paired_ends()
        before=context.network.to_dict()
        claims=prepare(context,bundle,{3,6},approaches,members)
        self.assertEqual([c.sections for c in claims],[(3,),(6,)])
        self.assertTrue(all(c.bundle_end=='TAIL' for c in claims))
        self.assertEqual([c.cuts for c in claims],[(End(2,1),),(End(5,1),)])
        self.assertEqual(context.network.to_dict(),before)
        with self.assertRaisesRegex(ValueError,'terminal junction'):
            prepare(context,bundle,{3,6},approaches,members-{Key('intersection',13)})
        with self.assertRaises(ValueError):
            prepare(context,bundle,{3},approaches,members)

    def test_asymmetric_terminal_consumption_keeps_the_opposite_original_half_profile(self):
        context,bundle,_,members=paired_ends()
        members=members-{Key('intersection',11)}
        cuts=(End(2,1),End(6,1))
        approaches=[Approach(Route.build(context,RouteTopology(context.network).trace(e)),e.end,
            context[Key('profile',e.section_id)]) for e in cuts]
        claims=prepare(context,bundle,{3},approaches,members)
        self.assertEqual([(c.role,c.sections,c.cuts) for c in claims],
            [(0,(3,),(End(2,1),)),(1,(),(End(6,1),))])
        self.assertEqual(claims[1].report()['policy'],'RETAINED_TAIL_HALF_PROFILE')
        with self.assertRaises(ValueError):
            prepare(context,bundle,{3},approaches[:1],members)
        with self.assertRaisesRegex(ValueError,'actually belongs'):
            prepare(context,bundle,{3},approaches,members|{Key('intersection',11)})

    def test_whole_pair_requires_all_four_original_terminal_owners(self):
        context,bundle,approaches,members=paired_ends(whole=True)
        claims=prepare(context,bundle,{1,2,3,4,5,6},approaches,members)
        self.assertTrue(all(c.bundle_end=='WHOLE' and not c.cuts for c in claims))
        with self.assertRaisesRegex(ValueError,'all original terminal'):
            prepare(context,bundle,{1,2,3,4,5,6},approaches,members-{Key('intersection',15)})

    def test_original_member_interval_and_both_half_profiles_survive_reversal(self):
        for reverse in (False,True):
            with self.subTest(reverse=reverse):
                context,bundle,approaches,members=sources(reverse)
                before=context.network.to_dict()
                claim,=prepare(context,bundle,{2,5},approaches,members)
                self.assertEqual(claim.sections,(2,))
                self.assertEqual(claim.role,0)
                self.assertEqual(claim.cuts,tuple(End(a.source.id,a.end) for a in approaches))
                self.assertEqual(context.network.to_dict(),before)
                self.assertEqual(claim.report()['policy'],'CONTIGUOUS_INTERNAL_MEMBER_INTERVAL')

    def test_missing_profile_or_wrong_original_site_is_rejected(self):
        context,bundle,approaches,members=sources()
        with self.assertRaisesRegex(ValueError,'half-profile'):
            prepare(context,bundle,{2},approaches[:1],members)
        with self.assertRaisesRegex(ValueError,'half-profile'):
            prepare(context,bundle,{2},approaches,{Key('intersection',10)})
        approaches[0].profile=replace(approaches[0].profile,half=0)
        with self.assertRaisesRegex(ValueError,'half-profile'):
            prepare(context,bundle,{2},approaches,members)

    def test_terminal_whole_body_and_two_lateral_claims_need_another_policy(self):
        context,bundle,approaches,members=sources()
        for consumed in ({1},{1,2,3},{2,4},{1,3}):
            with self.subTest(consumed=consumed),self.assertRaises(ValueError):
                prepare(context,bundle,consumed,approaches,members)
        incomplete=replace(bundle,members=(replace(bundle.members[0],complete=False),bundle.members[1]))
        with self.assertRaisesRegex(ValueError,'complete'):
            prepare(context,incomplete,{2},approaches,members)


if __name__=='__main__':unittest.main()
