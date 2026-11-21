"""A real narrow source lens admits bounded STYLE noses before meshing."""
from dataclasses import replace
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section, End, Intersection, Network, Key
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route, Piece
from tested_python_streets.profiles import Crown
from tested_python_streets.intersections import Approach
from tested_python_streets.regions import Region, RegionUnion
from tested_python_streets.corner_patch import corner_patch
from tested_python_streets.thin_median_corners import prepare, certify, radius_candidates


def fixture(reverse=False, angle=0., translation=(0., 0.), *, sealed=False, scope=False, mode='STYLE'):
    controls = (
        ((-242.05604553222656,-829.281005859375),(-242.4990692138672,-818.0377197265625),
         (-243.66567993164062,-808.4752807617188),(-246.66822814941406,-796.263427734375)),
        ((-246.66822814941406,-796.263427734375),(-245.75074768066406,-810.0114135742188),
         (-244.37416076660156,-821.5219116210938),(-242.05604553222656,-829.281005859375)))
    rotation = np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    sections=[]; products={}; arms={}
    for sid, xy in enumerate(controls, 1):
        xy=np.array(xy) @ rotation.T + translation
        if reverse: xy=xy[::-1]
        source=Section(sid,sid,tuple((x,y,0.) for x,y in xy),width=4. if sealed else 2.,scope_ends=(0,) if scope and sid==1 else ())
        sections.append(source); path=Path(source); products[Key('path',sid)]=path
        for end in (0,1): arms[End(sid,end)]=Approach(Route((Piece(path,end==0),)),end,Crown(.08,0.))
    bottom=(End(1,1 if reverse else 0),End(2,0 if reverse else 1))
    top=(End(2,1 if reverse else 0),End(1,0 if reverse else 1))
    # Both endpoints remain explicit physical junctions, each with its actual
    # third incidence. No finite cut is relabelled as a physical free end.
    sites=[]
    for i,pair in enumerate((bottom,top),1):
        origin=arms[pair[0]].origin
        destination=origin+np.array((10.,-10. if i==1 else 10.)) @ rotation.T
        source=Section(2+i,2+i,tuple((x,y,0.) for x,y in (origin,destination)),width=2.)
        sections.append(source); path=Path(source); products[Key('path',source.id)]=path
        endpoint=End(source.id,0)
        arms[endpoint]=Approach(Route((Piece(path,True),)),0,Crown(.08,0.))
        sites.append(Intersection(i,(*pair,endpoint),radius=1.,radius_modes=(mode,)*3))
    network=Network('thin-median',tuple(sections),tuple(sites))
    class Context:
        def __init__(self):self.network=network
        def __getitem__(self,key):return products[key]
    planner=SimpleNamespace(context=Context(),sites={s.key:s for s in sites},
        layouts={s.key:SimpleNamespace(approaches=tuple(arms[e] for e in s.ends)) for s in sites})
    return planner,(bottom,top),arms


class ThinMedianTests(unittest.TestCase):
    def test_actual_original_banks_reversed_rotated_and_whole_inner_arcs(self):
        for reverse,angle,translation in ((False,0.,(0.,0.)),(True,.713,(1200.,-1800.))):
            planner,pairs,arms=fixture(reverse,angle,translation)
            before=planner.context.network.to_dict();members=tuple(planner.sites)
            corners=[prepare(planner,pair,{1,2},members,exhausted=pairs) for pair in pairs]
            regions=[];sources=[]
            from tested_python_streets.compound_cycles import _canonical_samples
            cuts={}
            for c in corners:
                for a,s in zip((c.first,c.second),c.stations):
                    for (sid,end),at in a.path.trims(s).items():cuts.setdefault(sid,set()).add(at)
            sampled=_canonical_samples({sid:planner.context[Key('path',sid)] for sid in (1,2)},cuts,.5,.00025)
            for sid in (1,2):
                path=planner.context[Key('path',sid)]
                at=sampled[sid]
                regions.append(Region(Key('section',sid),np.vstack((path.point(at,1.),path.point(at,-1.)[::-1]))))
            for i,c in enumerate(corners):
                source=corner_patch(Key('nose',i),c,1.,.00025)
                sources.append(source);regions.append(source.region)
            loops=RegionUnion().build(regions)
            certify(corners,sources,loops)
            with self.assertRaisesRegex(ValueError,'whole original tangent-to-tangent'):
                certify(corners,sources,tuple(loop for loop in loops if not loop.hole))
            self.assertEqual(sum(loop.hole for loop in loops),1)
            for c in corners:
                report=c.thin_median_gore.report
                self.assertEqual(c.radius,.125)
                self.assertLessEqual(report['candidate_attempts'],2)
                self.assertLess(report['characteristic_radius_cap'],.25)
                self.assertGreater(report['old_circle_required_projected_width'],report['available_bank_width'])
                self.assertGreater(report['retained_hole_fraction'],.9)
                self.assertTrue(report['whole_inner_circle_exposed'])
                self.assertTrue(all(0. < s < a.length for a,s in zip((c.first,c.second),c.stations)))
            self.assertEqual(before,planner.context.network.to_dict())

    def test_strict_scope_missing_owner_or_false_hole_cannot_qualify(self):
        for options,message in ((dict(mode='STRICT'),'STYLE'),(dict(sealed=True),'genuine median')):
            planner,pairs,arms=fixture(**options)
            with self.subTest(options=options),self.assertRaisesRegex(ValueError,message):
                prepare(planner,pairs[0],{1,2},tuple(planner.sites),exhausted=pairs)
        with self.assertRaisesRegex(ValueError,'scope boundary'):
            fixture(scope=True)
        planner,pairs,arms=fixture()
        with self.assertRaisesRegex(ValueError,'authored physical'):
            prepare(planner,pairs[0],{1,2},(Key('intersection',1),),exhausted=pairs)
        with self.assertRaisesRegex(ValueError,'exhausted'):
            prepare(planner,pairs[0],{1,2},tuple(planner.sites),exhausted=())
        with self.assertRaisesRegex(ValueError,'complete consumed'):
            prepare(planner,pairs[0],{1},tuple(planner.sites),exhausted=pairs)
        a=arms[pairs[0][0]];original=a.profile
        for change in (dict(rise=.04),dict(datum=.01),dict(half=1)):
            a.profile=replace(original,**change)
            with self.subTest(change=change),self.assertRaisesRegex(ValueError,'full original'):
                prepare(planner,pairs[0],{1,2},tuple(planner.sites),exhausted=pairs)
        a.profile=original

    def test_characteristic_width_prefers_largest_fitting_bounded_style_nose(self):
        self.assertEqual(radius_candidates(5.360716,30.154124,.5),(.125,))
        self.assertEqual(radius_candidates(12.,40.,.5),(.25,.125))
        self.assertEqual(radius_candidates(1.,40.,.5),())
        self.assertEqual(radius_candidates(12.,40.,.125),())
        for values in ((-1.,40.,.5),(1.,float('inf'),.5),(1.,40.,0.)):
            with self.assertRaises(ValueError):radius_candidates(*values)

    def test_previously_successful_style_cycle_never_enters_thin_median_family(self):
        from unittest.mock import patch
        from test_python_street_source_cycles import triangle, prepare as prepare_cycle
        network=triangle(2.)
        network=replace(network,intersections=tuple(replace(site,radius_modes=('STYLE',)*len(site.ends))
                                                   for site in network.intersections))
        with patch('tested_python_streets.thin_median_corners.prepare') as fallback:
            context,group=prepare_cycle(network)
        fallback.assert_not_called()
        self.assertFalse(getattr(group.physical_layout,'source_cycle_thin_median_gores',()))
        self.assertEqual(sum(loop.hole for loop in group.loops),1)


if __name__=='__main__':unittest.main()
