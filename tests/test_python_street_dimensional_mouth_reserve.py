"""A finite width transition retains a complete exposed incoming profile."""
from copy import copy
from dataclasses import replace
import math
from types import SimpleNamespace
import unittest
from unittest.mock import patch
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
from tested_python_streets.ordinary_singleton_sources import (
    build,_dimensional_lead_reserve,OrdinarySourceSingletonUnavailable)
from tested_python_streets.source_strip_containment import SourceStripContainmentError


def fixture(reverse=False,angle=0.,shift=(0.,0.),mode='STYLE'):
    # Renumbered real three-arm controls: the wide through road meets two
    # curved narrow branches. No graph cycle or lateral Bundle is introduced.
    points=(
        ((-758.5760498046875,-547.478759765625),(-767.6460571289062,-540.0748901367188),
         (-772.6166381835938,-535.5213012695312),(-777.9837646484375,-529.6986083984375),
         (-781.3129272460938,-525.6572875976562),(-782.930419921875,-520.3026733398438)),
        ((-794.3450317382812,-536.5761108398438),(-786.5499877929688,-536.9889526367188),
         (-780.0230712890625,-538.04736328125),(-772.259521484375,-541.4769287109375),
         (-768.4243774414062,-542.9468383789062),(-758.5760498046875,-547.478759765625)),
        ((-758.5760498046875,-547.478759765625),(-736.4286499023438,-556.03076171875),
         (-711.6204833984375,-563.3807373046875),(-421.5722961425781,-640.7741088867188),
         (-412.4468688964844,-642.80078125)))
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    sections=[];rise=float(np.float32(.08))
    for sid,(controls,width) in enumerate(zip(points,(2.5,2.5,5.)),1):
        controls=np.asarray(controls)-points[0][0]
        if reverse:controls=controls[::-1]
        controls=controls@rotation.T+shift
        sections.append(Section(sid,sid,tuple((x,y,0.) for x,y in controls),
                                width=width,rise=rise,interpolation=3))
    ends=(End(1,0),End(2,1),End(3,0))
    if reverse:ends=tuple(End(e.section_id,1-e.end) for e in ends)
    site=Intersection(1,ends,3.,(3.,)*3,(mode,)*3)
    network=Network('dimensional-mouth-reserve',tuple(sections),(site,),
        settings=Settings(cell_size=float(np.float32(.7)),exaggeration=1.))
    products={Key('profile',s.id):Crown(s.rise,s.datum) for s in sections}
    context=BuildContext(network,products,PatternLibrary())
    for section in sections:products[Key('path',section.id)]=PathProvider.build(context,section,())
    topology=RouteTopology(network)
    for end in ends:products[route_key(end)]=Route.build(context,topology.trace(end))
    raw=IntersectionProvider.footprint(context,site)
    planner=SimpleNamespace(context=context,topology=topology,layouts={site.key:raw})
    return planner,site,raw


def candidate(planner,site,raw):
    try:correct(planner,site.key,raw)
    except SourceStripContainmentError as failure:
        return build(planner.context,site.key,raw,planner.layouts,(site,),failure)
    raise AssertionError('Fixture must exercise the complete-source-width failure')


class DimensionalMouthReserveTests(unittest.TestCase):
    def test_normal_update_retains_full_profiles_and_uses_only_the_normal_wide_reserve(self):
        from tested_python_streets.backend import PythonMeshBackend
        from tested_python_streets.ownership import PLAN,footprint_key,surface_key
        planner,site,raw=fixture();network=planner.context.network;before=network.to_dict()
        backend=PythonMeshBackend(network.settings);update=backend.prepare(network)
        group=update.products[PLAN].by_site[site.key]
        initial=update.products[footprint_key(site.key)]
        record,=group.physical_layout.ordinary_complete_mouth_clearance['corrections']
        self.assertEqual(len(record['dimensional_lead_reserve_advances']),1)
        advance,=record['dimensional_lead_reserve_advances']
        self.assertEqual(advance['section'],3)
        self.assertEqual(advance['candidate_route_station']-advance['original_route_station'],.25)
        self.assertFalse(record['original_mouth_planes_preserved'])
        self.assertEqual(group.members,(site.key,));self.assertEqual(group.claims,())
        old={(a.source.id,a.end):(a,s) for a,s in zip(initial.approaches,initial.outward_stations)}
        product=update.products[surface_key(site.key)]
        for a,s in zip(group.physical_layout.approaches,group.physical_layout.outward_stations):
            prior,t=old[a.source.id,a.end]
            self.assertIs(a.path,prior.path);self.assertIs(a.profile,prior.profile)
            if a.source.id!=3:self.assertEqual(a.path.trims(s),prior.path.trims(t))
            for u in (0.,.5,1.):
                loop,index,station=product.port_points[a.source.id,a.end,u]
                xyz=product.physical_rims[loop][index]
                np.testing.assert_allclose(xyz[:2],a.point(s,a.offset+a.width*(u-.5)),atol=1e-9,rtol=0.)
                self.assertAlmostEqual(xyz[2],a.profile.datum+float(a.profile.value(2*u-1.)),places=10)
        self.assertLessEqual(record['complete_source_width']['maximum_error'],.002)
        self.assertLessEqual(product.pavement.mesh.validate()['maximum_edge_ratio'],8.)
        self.assertEqual(network.to_dict(),before)
        backend.commit(update);self.assertEqual(backend.prepare(network).rebuilt,())

    def test_reversed_rotated_translated_sources_keep_the_same_bounded_rule(self):
        for reverse,angle,shift in ((False,0.,(0.,0.)),(True,.513,(1200.,-1800.))):
            planner,site,raw=fixture(reverse,angle,shift)
            before=planner.context.network.to_dict();group=candidate(planner,site,raw)
            records=group.physical_layout.ordinary_complete_mouth_clearance['dimensional_lead_reserve_advances']
            self.assertEqual(len(records),1);record,=records
            self.assertEqual(record['section'],3);self.assertEqual(record['end'],int(reverse))
            self.assertEqual(record['normal_reserve'],.25)
            self.assertLessEqual(record['candidate_route_station'],record['normal_mouth_station'])
            self.assertEqual(planner.context.network.to_dict(),before)

    def test_no_extra_movement_narrow_root_or_unsupported_primitive_is_admitted(self):
        planner,site,raw=fixture();group=candidate(planner,site,raw)
        old={(a.source.id,a.end):(a,s) for a,s in zip(raw.approaches,raw.outward_stations)}
        wide=next((a,t) for a,t in zip(group.physical_layout.approaches,
            group.physical_layout.outward_stations) if a.source.id==3)
        b,t=wide;a,s=old[b.source.id,b.end]
        self.assertIsNotNone(_dimensional_lead_reserve(planner.context,group,a,s,b,t))
        for changed in (s-.01,t+.01,a.length,float('nan')):
            self.assertIsNone(_dimensional_lead_reserve(planner.context,group,a,s,b,changed))
        narrow=next((arm,at) for arm,at in zip(group.physical_layout.approaches,
            group.physical_layout.outward_stations) if arm.source.id==2)
        self.assertIsNone(_dimensional_lead_reserve(planner.context,group,*old[narrow[0].source.id,narrow[0].end],
                                                  narrow[0],narrow[1]+.125))
        for name,value in (('profile',Crown(.04,0.)),('width',b.width+1.),('offset',.1),('end',1-b.end)):
            changed=copy(b);setattr(changed,name,value)
            self.assertIsNone(_dimensional_lead_reserve(planner.context,group,a,s,changed,t))
        stripped=copy(group.physical_layout)
        stripped.corners=tuple(replace(c) for c in group.physical_layout.corners)
        self.assertIsNone(_dimensional_lead_reserve(planner.context,replace(group,physical_layout=stripped),a,s,b,t))
        borrowed=copy(group.physical_layout)
        borrowed.corners=tuple(copy(c) for c in group.physical_layout.corners)
        self.assertTrue(any(hasattr(c,'dimensional_corner') for c in borrowed.corners))
        self.assertIsNone(_dimensional_lead_reserve(planner.context,replace(group,physical_layout=borrowed),a,s,b,t))
        without_operand=replace(group,footprints=tuple(source for source in group.footprints
            if not any(getattr(c,'dimensional_corner',None) is not None
                       and source is c.dimensional_corner.source for c in group.physical_layout.corners)))
        self.assertIsNone(_dimensional_lead_reserve(planner.context,without_operand,a,s,b,t))
        with patch('tested_python_streets.cycle_dimensional_retry.certify_layout',
                   side_effect=SourceStripContainmentError('fresh width proof failed',section_id=2)):
            with self.assertRaisesRegex(SourceStripContainmentError,'fresh width proof failed'):
                candidate(planner,site,raw)

    def test_explicit_radius_remains_strict_on_the_real_width_transition(self):
        planner,site,raw=fixture(mode='STRICT')
        group=candidate(planner,site,raw)
        self.assertTrue(all(c.radius_mode=='STRICT' and c.radius==3. and c.requested_radius==3.
                            for c in group.physical_layout.corners))
        with self.assertRaises(SourceStripContainmentError) as failure:
            correct(planner,site.key,raw)
        physical=copy(group.physical_layout)
        physical.corners=(replace(physical.corners[0],radius=2.),*physical.corners[1:])
        with patch('tested_python_streets.compound_cycles.source_cycle_group_steps',
                   return_value=_python_street_package.prepared_steps(replace(group,physical_layout=physical))):
            with self.assertRaisesRegex(OrdinarySourceSingletonUnavailable,'STRICT'):
                build(planner.context,site.key,raw,planner.layouts,(site,),failure.exception)


if __name__=='__main__':unittest.main()
