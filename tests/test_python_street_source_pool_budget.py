"""Nine original sites need measured source proof, not a raised global cap."""
from dataclasses import replace
from types import SimpleNamespace
import math
import unittest
import _python_street_package
from tested_python_streets.model import Section,End,Intersection,Network,Key
from tested_python_streets.source_pool_budget import certify, certify_steps
from tested_python_streets.source_cycle_atlas import eligible,OriginalCorridorCuts
from test_python_street_source_cycles import prepare


def ring():
    points=[(20.*math.cos(i*math.tau/9.),20.*math.sin(i*math.tau/9.),0.) for i in range(9)]
    sections=[];sites=[]
    for i,p in enumerate(points):
        sections.append(Section(i+1,i+1,(p,points[(i+1)%9]),width=2.))
        sections.append(Section(i+10,i+10,(p,tuple(1.75*x for x in p)),width=2.))
        sites.append(Intersection(i+1,(End(i+1,0),End((i-1)%9+1,1),End(i+10,0)),
            radius=1.,radius_modes=('STYLE',)*3))
    return Network('nine-original-source-sites',tuple(sections),tuple(sites))


class SourcePoolBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context,cls.group=prepare(ring())

    def test_nine_real_original_sites_pass_fresh_complete_width_and_measured_work(self):
        context,group=self.context,self.group;before=context.network.to_dict()
        report=certify(context,group)
        self.assertEqual(report['ordinary_site_budget'],8)
        self.assertEqual(report['original_sites'],9)
        self.assertEqual(report['complete_external_profiles'],9)
        self.assertEqual(report['physical_holes'],1)
        self.assertLessEqual(report['source_segments'],20000)
        self.assertEqual(report['maximum_source_error'],.002)
        self.assertTrue(report['downstream_surface_and_stored_checks_required'])
        self.assertEqual(before,context.network.to_dict())
        self.assertTrue(eligible(group,context))

    def test_uncertified_oversized_or_partial_source_pools_reject(self):
        context,group=self.context,self.group
        with self.assertRaisesRegex(ValueError,'nine to sixteen'):
            certify(context,replace(group,members=(*group.members,*(Key('intersection',i) for i in range(20,28)))))
        with self.assertRaisesRegex(ValueError,'complete original source-cycle'):
            certify(context,replace(group,physical_layout=SimpleNamespace()))
        with self.assertRaisesRegex(ValueError,'operand budget'):
            certify(context,replace(group,footprints=group.footprints*256))
        claim=replace(group.claims[0],first_trim=1.)
        with self.assertRaisesRegex(ValueError,'complete original source ownership'):
            certify(context,replace(group,claims=(claim,*group.claims[1:])))
        with self.assertRaisesRegex(ValueError,'physical-hole budget'):
            certify(context,replace(group,loops=(group.loops[0],*(group.loops[1:]*9))))

    def test_cancellation_during_width_check_preserves_published_proof(self):
        layout=self.group.physical_layout
        previous=getattr(layout,'extended_source_pool',None)
        previous_width=getattr(layout,'source_cycle_full_width_proof',None)
        steps=certify_steps(self.context,self.group)
        next(steps)
        steps.close()
        self.assertIs(getattr(layout,'extended_source_pool',None),previous)
        self.assertIs(getattr(layout,'source_cycle_full_width_proof',None),previous_width)

    def test_original_corridor_axes_do_not_infer_a_bundle_or_change_physical_claims(self):
        context,group=self.context,self.group
        certify(context,group)
        claims=group.claims;sid=claims[0].section_id
        selected,paired=OriginalCorridorCuts(sid)._source_inventory(group,context)
        self.assertEqual(selected,{sid});self.assertFalse(paired)
        self.assertIs(group.claims,claims)
        with self.assertRaisesRegex(ValueError,'consumed original'):
            OriginalCorridorCuts(999)._source_inventory(group,context)
        with_bundle=SimpleNamespace(network=SimpleNamespace(bundles=(SimpleNamespace(sections=(sid,)),)))
        self.assertFalse(eligible(group,with_bundle))


if __name__=='__main__':unittest.main()
