"""Actual source dimensions and explicit original centreline patch closure."""
from dataclasses import replace
import copy
import math
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_cycle_overlap_expansion import fixture as base_fixture
from tested_python_streets import model
from tested_python_streets.compound_displaced_contacts import retry,patch
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.intersections import Approach
from tested_python_streets.ownership import footprint_key
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import ProfileProvider
from tested_python_streets.providers import PathProvider
from tested_python_streets.routes import Route,RouteTopology
from tested_python_streets.updates import UpdateCoordinator,BuildContext


def fixture():
    n=base_fixture();origin=np.array((-1249.9444580078125,-685.5804443359375))
    controls=np.array(((-1247.052734375,-658.1073608398438),(-1243.987060546875,-655.7591552734375),
        (-1241.652587890625,-652.6871948242188),(-1240.212646484375,-649.1141357421875)))-origin
    sections=[]
    for s in n.sections:
        if s.id==4:s=replace(s,width=5.,points=tuple((*p,0.) for p in controls))
        if s.id==5:
            b=np.asarray(s.points[0]);s=replace(s,points=(s.points[0],tuple(b+(-20.,-5.,0.))))
        sections.append(s)
    end=controls[-1]
    for sid,delta in ((8,(-20.,20.)),(9,(20.,20.))):
        sections.append(model.Section(sid,sid,(tuple((*end,0.)),tuple((*(end+delta),0.))),width=5.))
    site=model.Intersection(13,(model.End(4,1),model.End(8,0),model.End(9,0)))
    return replace(n,sections=tuple(sections),intersections=tuple(replace(j,radius_modes=('STYLE',)*len(j.ends))
        for j in (*n.intersections,site)))


def context(n):
    patterns=PatternLibrary();products=UpdateCoordinator((PathProvider(),ProfileProvider()),patterns).prepare(n).products
    result=BuildContext(n,products,patterns);topology=RouteTopology(n)
    for site in n.intersections:
        approaches=tuple(Approach(Route.build(result,topology.trace(end)),end.end,
            result[model.Key('profile',end.section_id)]) for end in site.ends)
        products[footprint_key(site.key)]=SimpleNamespace(approaches=approaches)
    root=products[footprint_key(model.Key('intersection',10))]
    a=next(a for a in root.approaches if a.source.id==2)
    b=next(a for a in root.approaches if a.source.id==1)
    return result,a,b


class DisplacedSourceContactTests(unittest.TestCase):
    def test_future_real_contact_keeps_the_same_source_expansion_in_either_member_order(self):
        from tested_python_streets.compound_preparation import DeferredCompoundPlanner,CompoundExpansionRequired
        from tested_python_streets.compound_cycles import source_cycle_group
        requested=[]
        for ids in ((10,11),(11,10)):
            n=fixture();before=n.to_dict();c,a,b=context(n)
            # This actual finite circle contacts source4 at site11 while the
            # first member's source2/source1 corner is being processed.
            contact=retry(c,a,b,10.,'STYLE',set(map(lambda i:model.Key('intersection',i),ids)),{1})
            self.assertEqual(contact.second.source.id,4)
            layouts={j.key:c[footprint_key(j.key)] for j in n.intersections}
            planner=DeferredCompoundPlanner(c,n.intersections,layouts)
            with self.assertRaises(CompoundExpansionRequired) as caught:
                source_cycle_group(planner,tuple(model.Key('intersection',i) for i in ids),{1},
                    {model.End(1,0),model.End(1,1)},
                    {1:{0:model.Key('intersection',10),1:model.Key('intersection',11)}},set())
            requested.append(caught.exception.neighbors)
            self.assertEqual(n.to_dict(),before)
        self.assertEqual(requested[0],requested[1])

    def test_legacy_support_list_does_not_admit_a_self_intersecting_circle_patch(self):
        c,a,b=context(fixture())
        corner=retry(c,a,b,10.,'STYLE',{model.Key('intersection',10),model.Key('intersection',11)},{1})
        self.assertIsNotNone(corner)
        self.assertIsInstance(corner.source_contact_support,list)
        valid=patch(model.Key('proof',1),corner,.7,.00025)
        invalid=copy.copy(corner)
        invalid.sweep+=math.copysign(math.tau,invalid.sweep)
        with self.assertRaisesRegex(ValueError,'boundary intersects or touches itself'):
            patch(model.Key('invalid',1),invalid,.7,.00025)
        again=patch(model.Key('proof',1),corner,.7,.00025)
        np.testing.assert_array_equal(valid.region.xy,again.region.xy)
        self.assertEqual(valid.segments,again.segments)

    def qualify(self,n):
        before=n.to_dict();c,a,b=context(n)
        from tested_python_streets.compound_displaced_contacts import retry_steps
        args=(c,a,b,10.,'STYLE',{model.Key('intersection',10),model.Key('intersection',11)},{1})
        pending=retry_steps(*args)
        products=dict(c.products)
        next(pending)
        pending.close()
        self.assertEqual(products.keys(),c.products.keys())
        self.assertTrue(all(c.products[key] is product for key,product in products.items()))
        contact=retry(*args)
        self.assertIsNotNone(contact)
        self.assertIn(contact.radius,tuple(CornerRadiusPolicy.contact_candidates(a,b,10.,'STYLE')))
        self.assertGreaterEqual(contact.radius,.5)
        self.assertEqual(contact.first.source.id,2);self.assertEqual(contact.first.width,3.5)
        self.assertEqual(contact.second.source.id,4);self.assertEqual(contact.second.width,5.)
        self.assertGreater(contact.stations[0],0.);self.assertGreater(contact.stations[1],0.)
        self.assertLess(contact.stations[1],contact.second.length)
        self.assertEqual(contact.fit_policy,'DISPLACED_ORIGINAL_DIMENSIONAL_SOURCE_CURBS')
        self.assertEqual([p.section.id for p in contact.source_contact_walks[1].prefix],[1])
        alias=contact.source_contact_support[1]['actual_contact_aliases'][-1]
        self.assertEqual(alias['section'],4)
        self.assertAlmostEqual(alias['station'],contact.stations[1])
        f=patch(model.Key('proof',1),contact,.7,.00025)
        inventory={f.region.owner:f}
        for i,label in enumerate(f.segments):
            if label.kind=='curb':
                np.testing.assert_allclose(np.linalg.norm(f.region.xy[i]-contact.center),contact.radius,atol=1e-10)
            else:self.assertEqual(label.kind,'internal')
        self.assertEqual(before,n.to_dict())
        invalid=copy.copy(contact);invalid.stations=(contact.first.length+.1,contact.stations[1])
        with self.assertRaisesRegex(ValueError,'extrapolate'):
            patch(model.Key('invalid',1),invalid,.7,.00025)

    def test_different_actual_widths_and_original_closure_aliases(self):self.qualify(fixture())

    def test_reversed_storage_and_rotated_translated_original_controls(self):
        n=fixture();swap=lambda e:model.End(e.section_id,1-e.end)
        reversed=replace(n,sections=tuple(replace(s,points=s.points[::-1]) for s in n.sections),
            intersections=tuple(replace(j,ends=tuple(map(swap,j.ends))) for j in n.intersections))
        # Explicit source aliases follow the original stored orientation.
        c,a,b=context(reversed)
        contact=retry(c,a,b,10.,'STYLE',{model.Key('intersection',10),model.Key('intersection',11)},{1})
        self.assertIsNotNone(contact);self.assertEqual(contact.second.end,1)
        alias=contact.source_contact_support[1]['actual_contact_aliases'][-1]
        self.assertEqual(alias['end'],1)
        self.assertAlmostEqual(alias['station'],contact.second.length-contact.stations[1])
        patch(model.Key('proof',2),contact,.7,.00025)
        rotation=np.array(((np.cos(.7),-np.sin(.7)),(np.sin(.7),np.cos(.7))))
        transformed=[]
        for s in n.sections:
            p=np.asarray(s.points).copy();p[:,:2]=p[:,:2]@rotation.T+(1200.,-800.)
            transformed.append(replace(s,points=tuple(map(tuple,p))))
        self.qualify(replace(n,sections=tuple(transformed)))
        reflected=replace(n,sections=tuple(replace(s,points=tuple((x,-y,z) for x,y,z in s.points)) for s in n.sections))
        c,outside,inside=context(reflected)
        contact=retry(c,inside,outside,10.,'STYLE',{model.Key('intersection',10),model.Key('intersection',11)},{1})
        self.assertIsNotNone(contact)
        self.assertEqual(contact.first.source.id,4);self.assertEqual(contact.second.source.id,2)
        patch(model.Key('proof',3),contact,.7,.00025)

    def test_strict_unowned_connection_scope_and_missing_footprint_are_rejected(self):
        n=fixture();c,a,b=context(n);members={model.Key('intersection',10),model.Key('intersection',11)}
        self.assertIsNone(retry(c,a,b,10.,'STRICT',members,{1}))
        self.assertIsNone(retry(c,a,b,10.,'STYLE',{model.Key('intersection',10)},{1}))
        self.assertIsNone(retry(c,a,b,10.,'STYLE',members,set()))
        strict=replace(n,intersections=tuple(replace(j,radius_modes=('STRICT',)*len(j.ends)) if j.id==11 else j for j in n.intersections))
        sc,sa,sb=context(strict)
        self.assertIsNone(retry(sc,sa,sb,10.,'STYLE',members,{1}))
        c.products.pop(footprint_key(model.Key('intersection',11)))
        self.assertIsNone(retry(c,a,b,10.,'STYLE',members,{1}))
        cut=replace(n,sections=tuple(replace(s,scope_ends=(1,)) if s.id==1 else s for s in n.sections),
            intersections=tuple(j for j in n.intersections if j.id!=11))
        c,a,b=context(cut)
        self.assertIsNone(retry(c,a,b,10.,'STYLE',members,{1}))


if __name__=='__main__':unittest.main()
