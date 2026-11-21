"""A certified width-continuation recipient keeps its full native profile."""
from copy import copy
from dataclasses import replace
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,End,Intersection,Network,Settings,Key
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.profiles import Crown
from tested_python_streets.providers import PathProvider
from tested_python_streets.routes import Route,RouteTopology,route_key
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.updates import BuildContext
from tested_python_streets.ordinary_mouth_clearance import correct
from tested_python_streets.ordinary_singleton_sources import build,_ends,_dimensional_recipient_lead
from tested_python_streets.source_strip_containment import SourceStripContainmentError


def fixture(*,reverse=False,angle=0.,shift=(0.,0.),mode='STYLE'):
    controls=(
        ((1344.16552734375,212.04165649414062),(1328.467529296875,215.96786499023438),
         (1313.400146484375,217.07786560058594),(1294.078369140625,215.4930877685547),
         (1276.6083984375,211.86045837402344),(1237.2900390625,199.44041442871094)),
        ((1344.16552734375,212.04165649414062),(1329.7847900390625,203.3667755126953),
         (1322.747802734375,190.93089294433594),(1321.0694580078125,182.0138397216797)),
        ((1390.1976318359375,184.65602111816406),(1381.6396484375,191.4334716796875),
         (1371.52587890625,198.2885284423828),(1360.9302978515625,203.57388305664062),
         (1344.16552734375,212.04165649414062)))
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    sections=[]
    for sid,(points,width) in enumerate(zip(controls,(7.,5.400000095367432,5.400000095367432)),1):
        xy=np.asarray(points)@rotation.T+shift
        if reverse:xy=xy[::-1]
        sections.append(Section(sid,sid,tuple((x,y,0.) for x,y in xy),width=width,interpolation=3))
    ends=tuple(End(sid,(0 if sid!=3 else 1)^int(reverse)) for sid in (1,2,3))
    site=Intersection(1,ends,2.,(2.,)*3,(mode,)*3)
    network=Network('dimensional-recipient-source',tuple(sections),(site,),settings=Settings(exaggeration=1.))
    products={Key('profile',s.id):Crown(.08,0.) for s in sections}
    context=BuildContext(network,products,PatternLibrary())
    for section in sections:products[Key('path',section.id)]=PathProvider.build(context,section,())
    topology=RouteTopology(network)
    for end in ends:products[route_key(end)]=Route.build(context,topology.trace(end))
    layout=IntersectionProvider.footprint(context,site)
    planner=SimpleNamespace(context=context,layouts={site.key:layout},topology=topology)
    try:correct(planner,site.key,layout)
    except SourceStripContainmentError as error:failure=error
    else:raise AssertionError('The original ordinary layout must exercise its real width failure')
    group=build(context,site.key,layout,planner.layouts,(site,),failure)
    return planner,site,layout,group


class DimensionalRecipientLeadTests(unittest.TestCase):
    def test_real_recipient_advance_preserves_profiles_radius_and_owned_native_cuts(self):
        for reverse,angle,shift,mode in ((False,0.,(0.,0.),'STYLE'),
                (True,0.,(0.,0.),'STYLE'),(False,.713,(600.,-1200.),'STYLE')):
            with self.subTest(reverse=reverse,angle=angle,mode=mode):
                planner,site,layout,group=fixture(reverse=reverse,angle=angle,shift=shift,mode=mode)
                record=group.physical_layout.ordinary_complete_mouth_clearance
                lead,=record['dimensional_recipient_lead_advances']
                self.assertEqual(lead['section'],3)
                self.assertGreater(lead['advance'],0.)
                self.assertLess(lead['actual_old_exposed_profile_measure'],1.-1e-8)
                self.assertGreaterEqual(lead['actual_new_exposed_profile_measure'],1.-1e-8)
                self.assertEqual(lead['minimum_radius'],2.)
                self.assertLessEqual(lead['curvature_bound'],.5)
                self.assertFalse(record['original_mouth_planes_preserved'])
                self.assertTrue(record['original_paths_profiles_incidence_preserved'])
                old,new=_ends(layout),_ends(group.physical_layout)
                for key,(a,s) in old.items():
                    b,t=new[key]
                    self.assertIs(a,b);self.assertIs(a.path,b.path);self.assertIs(a.profile,b.profile)
                    if key[0]!=3:self.assertEqual(s,t)
                self.assertEqual(group.members,(site.key,));self.assertEqual(group.claims,())
                self.assertLessEqual(record['complete_source_width']['maximum_error'],.002)

    def test_an_unsupported_strict_circle_does_not_gain_width_continuation_authority(self):
        # The pre-existing circular-to-dimensional retry is STYLE-only here.
        # Recipient admission must not reinterpret an authored STRICT circle.
        for reverse in (False,True):
            with self.assertRaisesRegex(ValueError,'internal fillet preparation seam'):
                fixture(reverse=reverse,mode='STRICT')

    def test_only_the_actual_narrow_producer_and_exact_normal_reserve_are_admitted(self):
        planner,_,layout,group=fixture();old,new=_ends(layout),_ends(group.physical_layout)
        key=next(key for key in old if key[0]==3);a,s=old[key];b,t=new[key]
        self.assertIsNotNone(_dimensional_recipient_lead(planner.context,group,a,s,b,t))
        for before,after in ((copy(a),b),(a,copy(b))):
            self.assertIsNone(_dimensional_recipient_lead(planner.context,group,before,s,after,t))
        self.assertIsNone(_dimensional_recipient_lead(planner.context,group,a,s,b,t+.125))
        self.assertIsNone(_dimensional_recipient_lead(planner.context,group,a,t-.1,b,t))
        self.assertIsNone(_dimensional_recipient_lead(planner.context,replace(group,claims=(object(),)),a,s,b,t))
        donor=copy(group.physical_layout);donor.source_cycle_root_support_envelopes=(object(),)
        self.assertIsNone(_dimensional_recipient_lead(planner.context,replace(group,physical_layout=donor),a,s,b,t))
        corner=next(c for c in group.physical_layout.corners if hasattr(c,'dimensional_corner'))
        without=replace(group,footprints=tuple(source for source in group.footprints if source is not corner.dimensional_corner.source))
        self.assertIsNone(_dimensional_recipient_lead(planner.context,without,a,s,b,t))

    def test_changed_native_profiles_or_fabricated_curvature_reports_are_refused(self):
        planner,_,layout,group=fixture();old,new=_ends(layout),_ends(group.physical_layout)
        key=next(key for key in old if key[0]==3);a,s=old[key];b,t=new[key]
        products=dict(planner.context.products);products[Key('profile',3)]=Crown(.09,0.)
        context=BuildContext(planner.context.network,products,planner.context.patterns)
        self.assertIsNone(_dimensional_recipient_lead(context,group,a,s,b,t))
        corner=next(c for c in group.physical_layout.corners if hasattr(c,'dimensional_corner'))
        report=corner.dimensional_corner.report
        original=report['curvature_bound']
        try:
            report['curvature_bound']=float('nan')
            self.assertIsNone(_dimensional_recipient_lead(planner.context,group,a,s,b,t))
        finally:report['curvature_bound']=original


if __name__=='__main__':unittest.main()
