"""Outside recovery preserves source identity and rejects incomplete node plans."""
from copy import deepcopy
import importlib
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace as NS
import unittest

root = Path(__file__).resolve().parents[1]
for name,path in (
        ('_boundary_adapter_tests',root),
        ('_boundary_adapter_tests.way',root/'way'),
        ('_boundary_adapter_tests.way.street_network',root/'way/street_network'),
        ('_boundary_adapter_tests.way.gn',root/'way/gn'),
        ('_boundary_adapter_tests.way.gn.kernel',root/'way/gn/kernel')):
    package = ModuleType(name); package.__path__ = [str(path)]; sys.modules[name] = package
m = importlib.import_module('_boundary_adapter_tests.way.street_network.boundary_junctions')
c = importlib.import_module('_boundary_adapter_tests.way.clipped_context')
p = importlib.import_module('_boundary_adapter_tests.way.street_network.physical_ports')


class Snapshot(NS):
    def issue(self,kind,source,code,message):
        self.issues.append(dict(kind=kind,id=source,code=code,message=message))


def fixture(reverse=False):
    node = (0.,-1.)
    def seg(a,v,b,w):return NS(id1=a,v1=v,id2=b,v2=w)
    road = NS(source_id='road',category='residential',element=NS(tags={}),segments=[
        seg('a',(0.,30.),'node',node),seg('node',node,'b',(0.,-50.))])
    side = NS(source_id='side',category='residential',element=NS(tags={}),segments=[
        seg('c',(-30.,10.),'node',node)])
    context = c.ClippedSourceContext((-50.,50.,0.,50.),source_selection='SURFACE_VEHICLES')
    cuts = (30./31.,10./11.)
    points = (((0.,30.,0.),(0.,0.,0.)),((-30.,10.,0.),(-30.+30.*cuts[1],0.,0.)))
    streets = {}
    free = []
    for sid,way,t,coords in zip((1,2),(road,side),cuts,points):
        segment = way.segments[0]
        context.observe(segment,(way.source_id,segment.id1,segment.id2,0.,t),coords)
        interval = (t,0.) if reverse and sid==1 else (0.,t)
        if reverse and sid==1:coords=coords[::-1]
        section = dict(kind='Section',id=sid*10,category=way.category,width=6.,
            offset=-.5 if reverse and sid==1 else .5,points=coords,tags={},
            source_way_ids=[way.source_id],source_segments=[(way.source_id,segment.id1,segment.id2,*interval)])
        streets[sid] = dict(id=sid,items=[section],status='SUPPORTED',reason=None,
                            endpoints=(None,None),bundle_id=None)
        for end in (0,1):
            free.append(dict(endpoint=dict(street_id=sid,end=end,terminal_id=section['id'],
                             terminal_kind='Section'),location=coords[-1 if end else 0]))
    context.finalize((road,side))
    return Snapshot(streets=streets,retained_physical_streets={},bundles={},
        intersections={},junctions={},issues=[],preparation_audit={},
        junction_origins=dict(physical=[],logical=[],components=[],free_endpoints=free,seams=[]),
        clipped_source_context=c.resolve_open_ends(context.report(),streets.values()))


class BoundaryJunctionTests(unittest.TestCase):
    def setUp(self):
        self.design = NS(categories=('residential',),profile_rise=.08,cell_size=.7)
        self.adapter = m.BoundaryJunctionAdapter(self.design,{'residential':{'radius':2.}})

    def test_restores_all_legs_with_original_controls_dimensions_and_source_intervals(self):
        source = fixture(); before=deepcopy(source)
        self.adapter.apply(source)
        report = source.preparation_audit['boundary_junction_recovery']
        self.assertEqual(report['nodes'][0]['status'],'RESTORED')
        self.assertFalse(report['geometry_qualified'])
        jid, = source.junctions
        self.assertEqual(len(source.junction_origins['physical'][0]['ring']),3)
        self.assertEqual(set(source.streets),{1,2,3})
        for sid in (1,2):
            old=before.streets[sid]['items'][0];new=source.streets[sid]['items'][0]
            self.assertEqual(new['id'],old['id'])
            self.assertEqual(new['points'][:-1],old['points'])
            self.assertEqual(new['points'][-1],(0.,-1.,0.))
            self.assertEqual(new['source_segments'][:-1],tuple(old['source_segments']))
            self.assertEqual((new['width'],new['offset']),(old['width'],old['offset']))
        support, = source.streets[3]['items']
        record, = support['source_segments']
        self.assertEqual(record[:3],('road','node','b'))
        self.assertEqual(record[3:],(0.,1.))
        self.assertAlmostEqual(support['points'][1][1],-1.-49.*record[4])
        self.assertEqual(len(source.junction_origins['free_endpoints']),2)
        boundary, = source.junction_origins['selection_boundaries']
        self.assertEqual(boundary['endpoint'],dict(street_id=3,section_id=support['id'],end=1))
        self.assertGreater(boundary['trim_m'],0.)
        plan = p.PreparedPhysicalPorts(('residential',)).prepare(source,physical_ids=(jid,))
        self.assertEqual(set(plan.source_streets),{1,2,3})
        self.assertEqual(len(plan.sites[0].ports),3)
        span, = [span for span in plan.spans if span.street_id==3]
        self.assertEqual(span.end.owner_kind,'SELECTION_BOUNDARY')

    def test_reversed_source_retains_offset_orientation_and_stable_ids(self):
        forward=fixture();reverse=fixture(True)
        self.adapter.apply(forward);self.adapter.apply(reverse)
        self.assertEqual(set(forward.junctions),set(reverse.junctions))
        self.assertEqual(reverse.streets[1]['items'][0]['points'][0],(0.,-1.,0.))
        self.assertEqual(forward.streets[3]['items'][0]['offset'],reverse.streets[3]['items'][0]['offset'])
        self.assertEqual(forward.streets[3]['items'][0]['source_segments'],reverse.streets[3]['items'][0]['source_segments'])

    def rejected(self,change,pattern):
        source=fixture();change(source);before=deepcopy(source)
        self.adapter.apply(source)
        self.assertEqual(source.streets,before.streets)
        self.assertEqual(source.junction_origins,before.junction_origins)
        self.assertEqual(source.junctions,{})
        row,=source.preparation_audit['boundary_junction_recovery']['nodes']
        self.assertEqual(row['status'],'UNSUPPORTED');self.assertRegex(row['reason'],pattern)
        return source

    def test_missing_way_dimensions_do_not_create_half_a_junction(self):
        def change(s):s.clipped_source_context['nodes'][0]['legs'][1]['source_way_id']='new-road'
        self.rejected(change,'dimensions')

    def test_short_original_support_is_not_extrapolated(self):
        def change(s):
            leg=next(l for l in s.clipped_source_context['nodes'][0]['legs'] if l['other_node_id']=='b')
            leg['other_position']=(0.,-3.)
        self.rejected(change,'more source context')

    def test_grade_separated_leg_does_not_become_a_ground_connection(self):
        def change(s):s.clipped_source_context['nodes'][0]['legs'][0]['tags']['bridge']='yes'
        self.rejected(change,'elevation')

    def test_bundle_member_requires_joint_recovery(self):
        self.rejected(lambda s:s.streets[1].update(bundle_id=4),'unbundled')

    def test_explicit_trim_is_not_erased(self):
        self.rejected(lambda s:s.streets[1]['items'][0].update(clip_end=1.),'explicit selection trim')

    def test_missing_direction_provenance_is_not_guessed(self):
        self.rejected(lambda s:s.clipped_source_context['nodes'][0]['legs'][0].pop('source_nodes'),
                      'direction provenance')

    def test_inconsistent_original_leg_direction_is_not_guessed(self):
        def change(s):
            leg=s.clipped_source_context['nodes'][0]['legs'][0]
            leg['source_end']=1-leg['source_end']
        self.rejected(change,'original node incidence')

    def test_nonzero_level_is_reported(self):
        def change(s):s.clipped_source_context['nodes'][0]['legs'][0]['tags']['level']='1'
        self.rejected(change,'elevation')

    def test_unequal_widths_are_kept_and_new_support_inherits_only_its_own_way(self):
        source=fixture()
        source.streets[1]['items'][0]['width']=7.
        source.streets[2]['items'][0]['width']=2.
        self.adapter.apply(source)
        self.assertEqual([source.streets[i]['items'][0]['width'] for i in (1,2,3)],[7.,2.,7.])
        self.assertEqual(source.preparation_audit['boundary_junction_recovery']['nodes'][0]['status'],'RESTORED')

    def test_source_incidence_and_location_are_verified(self):
        self.rejected(lambda s:s.streets[1]['items'][0].update(points=((0.,30.,0.),(1.,0.,0.))),
                      'retained source interval')

    def test_distant_outside_junction_does_not_expand_import(self):
        source=fixture();source.clipped_source_context['nodes'][0]['boundary_distance_m']=1000.
        before=deepcopy(source.streets);self.adapter.apply(source)
        self.assertEqual(source.streets,before)
        self.assertEqual(source.preparation_audit['boundary_junction_recovery']['nodes'][0]['status'],
                         'OUTSIDE_INFLUENCE')

    def test_repeated_recovery_is_an_explicit_error(self):
        source=fixture();self.adapter.apply(source)
        with self.assertRaises(m.BoundaryContextError):self.adapter.apply(source)

    def test_support_responds_to_angles_and_rejects_coincident_legs(self):
        arms=[dict(width=6.,offset=0.,radius=2.,angle=a) for a in (0.,2.,4.)]
        baseline=self.adapter.support.lengths(arms)
        arms[1]['angle']=.1
        self.assertGreater(max(self.adapter.support.lengths(arms)),max(baseline))
        arms[1]['angle']=0.
        with self.assertRaises(m.BoundaryContextError):self.adapter.support.lengths(arms)

    def test_four_and_five_retained_legs_use_the_same_identity_and_port_contract(self):
        import math
        for count,span in ((4,160.),(5,160.),(4,135.)):
            source=fixture();source.streets={};source.junction_origins['free_endpoints']=[]
            context=c.ClippedSourceContext((-50.,50.,0.,50.));ways=[]
            for i in range(count):
                angle=math.radians((180.-span)*.5+span*i/(count-1))
                far=(30.*math.cos(angle),30.*math.sin(angle)-1.)
                cut=(math.cos(angle)/math.sin(angle),0.)
                t=1.-1./(30.*math.sin(angle))
                segment=NS(id1=str(i),v1=far,id2='node',v2=(0.,-1.))
                way=NS(source_id='way'+str(i),category='residential',element=NS(tags={}),segments=[segment])
                ways.append(way);context.observe(segment,(way.source_id,str(i),'node',0.,t),(far,cut))
                section=dict(kind='Section',id=100+i,category='residential',width=6.,offset=0.,
                    points=(far+(0.,),cut+(0.,)),tags={},source_way_ids=[way.source_id],
                    source_segments=[(way.source_id,str(i),'node',0.,t)])
                source.streets[i+1]=dict(id=i+1,items=[section],status='SUPPORTED',reason=None,
                    endpoints=(None,None),bundle_id=None)
                for end in (0,1):
                    source.junction_origins['free_endpoints'].append(dict(endpoint=dict(street_id=i+1,
                        end=end,terminal_id=100+i,terminal_kind='Section'),location=section['points'][end]))
            context.finalize(ways)
            source.clipped_source_context=c.resolve_open_ends(context.report(),source.streets.values())
            before=deepcopy(source.streets);self.adapter.apply(source)
            site,=source.junction_origins['physical'];jid=site['junction_id']
            if span==135.:
                self.assertFalse(source.junctions)
                self.assertIn('210 degrees',source.intersections[jid]['reason'])
            else:self.assertIn(jid,source.junctions)
            plan=p.PreparedPhysicalPorts(('residential',)).prepare(source,physical_ids=(jid,))
            self.assertEqual(len(plan.sites[0].ports),count)
            self.assertEqual(set(source.streets),set(before))
            for sid,old in before.items():
                self.assertEqual(source.streets[sid]['items'][0]['points'][:-1],old['items'][0]['points'])


if __name__=='__main__':unittest.main()
