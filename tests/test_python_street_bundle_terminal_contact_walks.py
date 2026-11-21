"""A consumed original link can close a circle on real opposed half banks."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets import model
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.routes import Piece,Route
from tested_python_streets.intersections import Approach
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.ownership import footprint_key
from tested_python_streets.bundle_terminal_contact_walks import contact_pair,retry
from tested_python_streets.compound_displaced_contacts import patch
from tested_python_streets.regions import RegionUnion


def fixture(reverse=False,reflection=False,angle=0.,translation=(0.,0.)):
    # Original controls around a short transverse terminal link. Each contact
    # keeps its own half profile; the link remains a full original road.
    controls={
        29:((825.7032470703125,309.91363525390625),(823.9221801757812,355.5766906738281),
            (821.952880859375,404.82421875),(820.5657958984375,439.94537353515625),
            (820.283203125,447.203369140625),(820.1602172851562,479.3301696777344),
            (820.0282592773438,515.3976440429688),(830.8951416015625,544.7649536132812)),
        190:((836.8914184570312,448.5966796875),(836.1804809570312,468.67864990234375),
             (834.613037109375,512.9613647460938),(833.64599609375,532.2640380859375)),
        397:((833.64599609375,532.2640380859375),(830.8951416015625,544.7649536132812)),
        398:((833.64599609375,532.2640380859375),(843.1497802734375,532.5433959960938)),
        402:((830.8951416015625,544.7649536132812),(831.457763671875,553.30322265625))}
    matrix=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))@np.diag((1.,-1. if reflection else 1.))
    sections=[]
    for sid,points in controls.items():
        xy=(np.asarray(points)-np.array((830.,540.)))@matrix.T+translation
        if reverse:xy=xy[::-1]
        sections.append(model.Section(sid,sid,tuple((x,y,0.) for x,y in xy),
            width=5. if sid==398 else 10.800000190734863,interpolation=3))
    end=lambda sid,e:model.End(sid,e^int(reverse))
    far=model.Intersection(794,tuple(end(s,e) for s,e in ((398,0),(397,0),(190,1))),
        radius=4.,radius_modes=('STYLE',)*3)
    root=model.Intersection(802,tuple(end(s,e) for s,e in ((402,0),(29,1),(397,1))),
        radius=4.,radius_modes=('STYLE',)*3)
    members=tuple(model.BundleMember((model.Traversal(sid,not reverse),)) for sid in (29,190))
    if reflection:members=members[::-1]
    bundle=model.Bundle(12,members)
    n=model.Network('terminal-bank-contact',tuple(sections),intersections=(far,root),bundles=(bundle,))
    p={model.Key('path',s.id):Path(s) for s in sections}
    for s in sections:
        half=next(((1 if role==0 else -1)*(1 if member.path[0].forward else -1)
            for role,member in enumerate(members) if member.path[0].section_id==s.id),0)
        p[model.Key('profile',s.id)]=Crown(.08,0.,half)
    approaches={e:Approach(Route((Piece(p[model.Key('path',e.section_id)],e.end==0),)),e.end,
        p[model.Key('profile',e.section_id)]) for site in (far,root) for e in site.ends}
    p[footprint_key(far.key)]=SimpleNamespace(approaches=tuple(approaches[e] for e in far.ends))
    p[footprint_key(root.key)]=SimpleNamespace(approaches=tuple(approaches[e] for e in root.ends))
    context=BuildContext(n,p,PatternLibrary())
    pair=(approaches[end(29,1)],approaches[end(397,1)])
    if reflection:pair=pair[::-1]
    return context,pair,{far.key,root.key},{397}


class BundleTerminalContactWalkTests(unittest.TestCase):
    def test_original_half_banks_and_consumed_link_survive_storage_and_rigid_changes(self):
        for reverse,reflection,angle,translation in ((False,False,0.,(0.,0.)),
                (True,False,.713,(1800.,-1200.)),(True,True,-.417,(1200.,900.))):
            with self.subTest(reverse=reverse,reflection=reflection):
                c,pair,members,consumed=fixture(reverse,reflection,angle,translation)
                before=c.network.to_dict();prepared=contact_pair(c,*pair,members,consumed)
                self.assertIsNotNone(prepared)
                corner=retry(c,*pair,4.,'STYLE',members,consumed)
                self.assertIsNotNone(corner)
                self.assertAlmostEqual(corner.radius,2.56)
                self.assertEqual(corner.radius_attempts,(4.,3.2,2.5600000000000005))
                self.assertEqual({corner.first.source.id,corner.second.source.id},{29,190})
                self.assertAlmostEqual(float(corner.first.profile.value(-1.)),.08)
                self.assertAlmostEqual(float(corner.second.profile.value(1.)),.08)
                prefix=[w for w in corner.source_contact_walks if w.prefix]
                self.assertEqual(len(prefix),1)
                self.assertEqual([p.section.id for p in prefix[0].prefix],[397])
                self.assertFalse(prefix[0].original.profile.half)
                source=patch(model.Key('proof',1),corner,c.network.settings.spacing,.00025)
                self.assertEqual(len(RegionUnion().build((source.region,))),1)
                self.assertEqual(c.network.to_dict(),before)

    def test_no_unowned_site_profile_substitution_wrong_bank_or_strict_radius_admission(self):
        c,pair,members,consumed=fixture()
        self.assertIsNone(retry(c,*pair,4.,'STRICT',members,consumed))
        self.assertIsNone(contact_pair(c,*pair,{model.Key('intersection',802)},consumed))
        self.assertIsNone(contact_pair(c,*pair,members,set()))
        original=pair[0].profile
        for profile in (Crown(.08,0.,1),Crown(.09,0.,-1),Crown(.08,1.,-1),Crown(.08,0.,0)):
            pair[0].profile=profile
            self.assertIsNone(contact_pair(c,*pair,members,consumed))
        pair[0].profile=original
        pair[0].width-=.1
        self.assertIsNone(contact_pair(c,*pair,members,consumed))
        pair[0].width+=.1
        # Same coordinates without the real far source incidence cannot own
        # a contact walk, even when its local source curves would fit a circle.
        far=c.network.intersections[0]
        changed=replace(c.network,intersections=(replace(far,ends=(far.ends[0],far.ends[1]),
            radius_modes=('STYLE',)*2),c.network.intersections[1]))
        self.assertIsNone(contact_pair(BuildContext(changed,c.products,c.patterns),*pair,members,consumed))


if __name__=='__main__':unittest.main()
