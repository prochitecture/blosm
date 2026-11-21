"""Preferred radii retain complete ordinary source widths and fixed mouths."""
import math
from types import SimpleNamespace
from unittest import mock
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
from tested_python_streets.ordinary_mouth_clearance import correct,_group,_canonical_candidate
from tested_python_streets.ordinary_source_width_fit import fit
from tested_python_streets.cycle_dimensional_retry import certify_layout
from tested_python_streets.source_strip_containment import SourceStripContainmentError


def fixture(*,reverse=False,angle=0.,shift=(0.,0.),strict=False):
    # Renumbered original smooth controls: the R10 contact erodes the root of
    # source2 by9.48mm analytically. Complete mouth planes are already clear.
    controls=(
      ((-58.6639404296875,-75.82923889160156),(-62.433837890625,-62.35871887207031),
       (-76.3367919921875,-34.08045959472656),(-78.25634765625,-25.820144653320312),
       (-78.37109375,-17.215118408203125),(-75.4830322265625,-8.187744140625),
       (-70.9610595703125,1.796600341796875),(-61.3519287109375,10.221343994140625),
       (-49.35546875,14.916366577148438),(-37.8497314453125,16.617034912109375),
       (-24.8433837890625,14.71063232421875),(-13.176025390625,8.874984741210938),(0.,0.)),
      ((0.,0.),(8.5697021484375,-14.039199829101562),(13.4288330078125,-27.543258666992188),
       (19.8660888671875,-46.11271667480469)),
      ((0.,0.),(-13.4542236328125,17.435516357421875),(-29.258544921875,30.4410400390625),
       (-73.1048583984375,50.933502197265625),(-85.751953125,58.4615478515625),
       (-98.7799072265625,67.79307556152344),(-115.8204345703125,82.54682922363281),
       (-124.96630859375,92.10015869140625),(-132.6107177734375,102.73295593261719),
       (-219.2685546875,262.00747871398926),(-232.4613037109375,293.90379333496094)))
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    sections=[]
    for sid,points in enumerate(controls,1):
        points=np.asarray(points)
        if reverse:points=points[::-1]
        points=points@rotation.T+shift
        sections.append(Section(sid,sid,tuple((x,y,0.) for x,y in points),width=3.5,interpolation=3))
    ends=(End(3,0),End(1,1),End(2,0))
    if reverse:ends=tuple(End(e.section_id,1-e.end) for e in ends)
    site=Intersection(1,ends,10.,(10.,10.,10.),('STRICT',)*3 if strict else ('STYLE',)*3)
    n=Network('ordinary-style-width',tuple(sections),(site,),settings=Settings(exaggeration=1.))
    p={Key('profile',s.id):Crown(.08,0.) for s in sections};c=BuildContext(n,p,PatternLibrary())
    for s in sections:p[Key('path',s.id)]=PathProvider.build(c,s,())
    topology=RouteTopology(n)
    for end in ends:p[route_key(end)]=Route.build(c,topology.trace(end))
    layout=IntersectionProvider.footprint(c,site)
    return SimpleNamespace(context=c,topology=topology,layouts={site.key:layout}),site,layout


class OrdinarySourceWidthFitTests(unittest.TestCase):
    def test_existing_style_candidates_keep_original_mouths_and_other_corners(self):
        for reverse,angle,shift in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            planner,site,raw=fixture(reverse=reverse,angle=angle,shift=shift)
            before=planner.context.network.to_dict()
            with self.assertRaises(SourceStripContainmentError):
                _canonical_candidate(planner.context,site.key,raw)
            result=correct(planner,site.key,raw)
            changed=[i for i,c in enumerate(result.corners) if c is not raw.corners[i]]
            self.assertEqual(len(changed),1)
            corner=result.corners[changed[0]]
            self.assertAlmostEqual(corner.radius,3.2768,places=10)
            self.assertEqual(corner.requested_radius,10.)
            self.assertEqual(corner.radius_attempts,(10.,8.,6.400000000000001,5.120000000000001,
                                                    4.096000000000001,3.2768000000000006))
            self.assertIs(result.approaches,raw.approaches)
            self.assertIs(result.mouths,raw.mouths)
            self.assertIs(result.outward_stations,raw.outward_stations)
            record=result.ordinary_complete_mouth_clearance
            self.assertFalse(record['original_corners_contacts_radii_preserved'])
            self.assertEqual(record['style_source_width_fit']['actual_radius'],corner.radius)
            trials=record['style_source_width_fit']['candidate_trials']
            self.assertEqual([trial['radius'] for trial in trials],list(corner.radius_attempts[1:]))
            self.assertTrue(all(trial['result']=='complete_source_width_rejected' for trial in trials[:-1]))
            self.assertEqual(trials[-1]['result'],'complete_initial_and_canonical_source_width_accepted')
            self.assertTrue(record['original_paths_profiles_incidence_preserved'])
            self.assertEqual(planner.context.network.to_dict(),before)
            proof=certify_layout(planner.context,_group(site.key,result),())
            self.assertEqual({r['section'] for r in proof['sections']},{1,2,3})
            self.assertIs(correct(planner,site.key,result),result)

    def test_explicit_radii_and_unrelated_failures_remain_strict(self):
        planner,site,raw=fixture(strict=True)
        with self.assertRaises(SourceStripContainmentError):correct(planner,site.key,raw)
        self.assertTrue(all(c.radius==10. and c.radius_mode=='STRICT' for c in raw.corners))
        planner,site,raw=fixture()
        failure=SourceStripContainmentError('real source-width witness',section_id=2)
        with self.assertRaisesRegex(ValueError,'allocation witness'):
            fit(planner.context,site.key,raw,failure,
                mock.Mock(side_effect=ValueError('allocation witness')))
        with self.assertRaisesRegex(ValueError,'finite candidate budget'):
            fit(planner.context,site.key,raw,failure,
                mock.Mock(side_effect=failure),maximum_candidates=1)
        self.assertIsNone(fit(planner.context,site.key,raw,
            SourceStripContainmentError('no original source identity'),mock.Mock()))


if __name__=='__main__':unittest.main()
