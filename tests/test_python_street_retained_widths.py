"""A complete visible profile does not prove every consumed original width."""
from dataclasses import replace
from types import SimpleNamespace,MethodType
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,End,Intersection,Network,Key
from tested_python_streets.paths import Path
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.compounds import CompoundGroup,SharedClaim,SampledFootprint,RimSpan,RimSegment
from tested_python_streets.regions import BoundaryLoop,Region
from tested_python_streets.compound_ports import CompoundPortPlanner
from tested_python_streets.compound_preparation import DeferredCompoundPlanner
from tested_python_streets.compound_curbs import RetainedSourceWidthError


def fixture(reverse=False,angle=0.,offset=(0.,0.)):
    rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
    def points(xy):return tuple((float(x),float(y),0.) for x,y in np.asarray(xy)@rotation.T+offset)
    source=Section(1,1,points(((0.,0.),(10.,0.))),width=2.)
    if reverse:source=replace(source,points=source.points[::-1])
    sections=(source,Section(2,2,points(((0.,0.),(-10.,4.)))),
        Section(3,3,points(((0.,0.),(-10.,-4.)))),
        Section(4,4,points(((10.,0.),(20.,4.)))),Section(5,5,points(((10.,0.),(20.,-4.)))))
    sites=(Intersection(1,(End(1,int(reverse)),End(2,0),End(3,0))),
           Intersection(2,(End(1,1-int(reverse)),End(4,0),End(5,0))))
    network=Network('retained-full-width',sections,sites)
    context=BuildContext(network,{Key('path',1):Path(source)},PatternLibrary())
    xy=np.asarray(points(((0.,-1.),(10.,-1.),(10.,.8),(0.,.8))))[:,:2]
    owner=sites[0].key
    labels=tuple(RimSegment('curb',0,float(i)/4.,float(i+1)/4.) for i in range(4))
    layout=SimpleNamespace(mouths={},absorbed_transitions=())
    footprint=SampledFootprint(Region(owner,xy),labels,layout)
    loop=BoundaryLoop(xy,tuple(RimSpan(owner,label) for label in labels))
    claim=SharedClaim(1,sites[int(reverse)].key,sites[1-int(reverse)].key,10.,0.)
    group=CompoundGroup(tuple(site.key for site in sites),(claim,),(footprint,),(loop,))
    return context,group,{site.key:layout for site in sites}


class RetainedWidthTests(unittest.TestCase):
    def test_complete_profiles_cannot_hide_a_clipped_whole_original_strip(self):
        for reverse,angle,offset in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            context,group,layouts=fixture(reverse,angle,offset);before=context.network.to_dict()
            ports=CompoundPortPlanner(context,context.network.intersections,layouts)
            self.assertFalse(ports.incomplete_ports(group))
            self.assertTrue(ports.incomplete_source_widths(group))
            planner=DeferredCompoundPlanner(context,context.network.intersections,layouts)
            calls=[]
            def exterior(self,members):
                yield
                calls.append(members)
                return replace(group,physical_layout=SimpleNamespace(mouths={}))
            planner.exterior_steps=MethodType(exterior,planner)
            result,=planner.build((group,))
            self.assertEqual(calls,[group.members]);self.assertEqual(result.members,group.members)
            self.assertEqual(context.network.to_dict(),before)

    def test_already_contained_width_and_existing_physical_policy_keep_old_group(self):
        context,group,layouts=fixture();ports=CompoundPortPlanner(context,context.network.intersections,layouts)
        xy=np.array(((0.,-1.),(10.,-1.),(10.,1.),(0.,1.)))
        full=replace(group,loops=(replace(group.loops[0],xy=xy),))
        self.assertFalse(ports.incomplete_source_widths(full))
        self.assertFalse(ports.incomplete_source_widths(replace(group,claims=())))
        self.assertFalse(ports.incomplete_source_widths(replace(group,physical_layout=SimpleNamespace())))

    def test_typed_failure_retains_the_exact_original_source_inventory(self):
        error=RetainedSourceWidthError((9,2,7),ValueError('source strip is clipped'))
        self.assertEqual(error.sections,(2,7,9))
        self.assertIn('complete original source widths',str(error))
        self.assertIn('source strip is clipped',str(error))

    def test_unavailable_certificate_is_not_an_actual_width_clipping_trigger(self):
        context,group,layouts=fixture()
        ports=CompoundPortPlanner(context,context.network.intersections,layouts)
        with patch('tested_python_streets.source_strip_containment.certify',
                side_effect=ValueError('Complete source-strip proof exceeds its finite cell budget')):
            with self.assertRaisesRegex(ValueError,'finite cell budget'):
                ports.incomplete_source_widths(group)


if __name__=='__main__':unittest.main()
