"""An exact exposed source cap may activate its original internal contact."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np

import _python_street_package
from tested_python_streets.model import Section,End,Key,Intersection
from tested_python_streets.paths import Path
from tested_python_streets.intersections import Approach
from tested_python_streets.profiles import Crown
from tested_python_streets.compound_cycles import _footprint
from tested_python_streets.cycle_exposed_root_contacts import required
from tested_python_streets.regions import Region,RegionUnion


def fixture(reverse=False,angle=0.,shift=(0.,0.)):
    rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
    approaches=[];sources=[];ends=[]
    for sid,direction,width in ((1,(1.,0.),4.),(2,(-1.,0.),8.),(3,(0.,-1.),2.)):
        xy=np.array(((0.,0.),np.asarray(direction)*12.))@rotation.T+shift
        points=tuple((*p,0.) for p in xy)
        if reverse:points=points[::-1]
        a=Approach(Path(Section(sid,sid,points,width=width,interpolation=0)),
                   int(reverse),Crown(.08,0.))
        end=End(sid,int(reverse));ends.append(end);approaches.append(a)
        at=np.array((0.,a.length))
        polygon=np.vstack((a.point(at,.5*width),a.point(at,-.5*width)[::-1]))
        source=_footprint(Key('original_strip',sid),polygon,internal_edges=(3,))
        source.layout.original_root_caps={3:end};sources.append(source)
    site=Intersection(100,tuple(ends),radius_modes=('STYLE',)*3)
    planner=SimpleNamespace(sites={site.key:site},
        layouts={site.key:SimpleNamespace(approaches=tuple(approaches))})
    loops=RegionUnion().build(tuple(s.region for s in sources))
    pair=tuple(ends[:2]);unresolved=((*pair,'finite original circle failed'),)
    return planner,loops,sources,pair,unresolved


class ExposedRootContactTests(unittest.TestCase):
    def test_real_union_cap_activates_only_its_original_consumed_style_pair(self):
        for reverse,angle,shift in ((False,0.,(0.,0.)),(True,.71,(1200.,-800.))):
            planner,loops,sources,pair,unresolved=fixture(reverse,angle,shift)
            found=required(planner,loops,sources,unresolved,set(pair),{1,2},())
            self.assertIsNotNone(found)
            selected,proof=found
            self.assertEqual(selected,pair)
            self.assertTrue(proof['exact_facing_root_cap_exposed'])
            self.assertTrue(proof['global_full_cap_and_source_width_guards_unchanged'])
            self.assertIsNone(required(planner,loops,sources,unresolved,set(pair),{1,2},(pair,)))
            self.assertIsNone(required(planner,loops,sources,unresolved,{pair[0]},{1,2},()))
            self.assertIsNone(required(planner,loops,sources,unresolved,set(pair),{1},()))
            site=next(iter(planner.sites.values()))
            planner.sites[site.key]=replace(site,radius_modes=('STRICT',)*3)
            self.assertIsNone(required(planner,loops,sources,unresolved,set(pair),{1,2},()))

    def test_root_identity_and_geometry_cannot_be_replaced_by_a_nearby_tag(self):
        planner,loops,sources,pair,unresolved=fixture()
        sources[1].layout.original_root_caps={3:End(99,0)}
        self.assertIsNone(required(planner,loops,sources,unresolved,set(pair),{1,2},()))
        sources[1].layout.original_root_caps={3:pair[1]}
        moved=replace(sources[1],region=Region(sources[1].region.owner,sources[1].region.xy+(1e-4,0.)))
        self.assertIsNone(required(planner,loops,(sources[0],moved,sources[2]),unresolved,set(pair),{1,2},()))
        a=planner.layouts[next(iter(planner.sites))].approaches[1]
        a.profile=replace(a.profile,half=-1)
        self.assertIsNone(required(planner,loops,sources,unresolved,set(pair),{1,2},()))


if __name__=='__main__':unittest.main()
