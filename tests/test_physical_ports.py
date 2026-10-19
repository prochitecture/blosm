"""Original Section-end addressing must survive geometry edits and regrouping."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest

spec = importlib.util.spec_from_file_location('tested_physical_ports',
    Path(__file__).resolve().parents[1]/'way/street_network/physical_ports.py')
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)


class GroundSurfacePolicyTests(unittest.TestCase):
    def test_ground_passages_preserve_original_tags(self):
        for tags in ({}, {'tunnel': 'no'}, {'tunnel': 'building_passage'},
                     {'tunnel': 'building_passage', 'layer': '0', 'level': '0',
                      'bridge': 'no', 'indoor': 'no', 'location': 'surface'},
                     {'tunnel': 'building_passage', 'layer': 0, 'level': 0}):
            original = deepcopy(tags)
            self.assertIsNone(m.GroundSurfacePolicy.reason(tags))
            self.assertEqual(tags, original)

    def test_passage_does_not_bypass_grade_or_indoor_treatment(self):
        for key, values in {'tunnel': ('yes', 'culvert'), 'bridge': ('yes', 'viaduct'),
                            'layer': ('-1', '1', '0;1'), 'level': ('-1', '1', '0;1'),
                            'indoor': ('yes', 'room'), 'location': ('underground', 'roof')}.items():
            for value in values:
                with self.subTest(key=key, value=value):
                    tags = {'tunnel': 'building_passage', key: value}
                    original = deepcopy(tags)
                    self.assertIsNotNone(m.GroundSurfacePolicy.reason(tags))
                    self.assertEqual(tags, original)

    def test_mixed_category_scope_applies_same_passage_policy(self):
        road = dict(kind='Section', id=1, category='residential', width=6., offset=0.,
                    points=((0., 0.), (20., 0.)), tags={'tunnel': 'building_passage'})
        foot = dict(road, id=2, category='footway', points=((20., 0.), (30., 0.)), tags={})
        source = dict(id=3, items=[road, foot], status='UNSUPPORTED',
                      reason='Reserved component/category: footway')
        original = deepcopy(source)
        self.assertEqual(m.StreetSectionScope.prepare(source, ('residential',)).included_sections, (1,))
        self.assertEqual(source, original)
        road['tags']['level'] = '1'
        with self.assertRaisesRegex(ValueError, 'Grade/elevation'):
            m.StreetSectionScope.prepare(source, ('residential',))


class ThroughConnectorTests(unittest.TestCase):
    def test_authored_through_connectors_survive_original_span_planning(self):
        data=fixture(); item=data.streets[1]['items'][1]
        item['through_connectors']=[dict(section_id=10,end=1,connector_id=901),
                                    dict(section_id=11,end=0,connector_id=902)]
        before=deepcopy(data)
        plan=m.PreparedPhysicalPorts(('residential',)).prepare(data,interior_ids=(21,))
        site=plan.sites[0]
        self.assertEqual({(p.endpoint.section_id,p.endpoint.end,p.source_kind,p.connector_id) for p in site.ports},
                         {(10,1,'THROUGH',901),(11,0,'THROUGH',902),(20,1,'ATTACHMENT',201)})
        self.assertEqual(data,before)
        for port in site.ports:
            span=plan.span_for(port)
            self.assertIn(port.endpoint,(span.start.endpoint,span.end.endpoint))

    def test_incomplete_duplicate_or_foreign_through_connector_records_reject(self):
        valid=[dict(section_id=10,end=1,connector_id=901),dict(section_id=11,end=0,connector_id=902)]
        for change in ('empty','missing','duplicate_end','foreign','bad_end','duplicate_id','negative','bool','attachment_id'):
            data=fixture(); rows=deepcopy(valid)
            if change=='empty':rows=[]
            elif change=='missing':rows.pop()
            elif change=='duplicate_end':rows[1]=dict(rows[0],connector_id=902)
            elif change=='foreign':rows[1]['section_id']=12
            elif change=='bad_end':rows[1]['end']=1
            elif change=='duplicate_id':rows[1]['connector_id']=901
            elif change=='negative':rows[1]['connector_id']=-1
            elif change=='bool':rows[1]['connector_id']=True
            else:rows[1]['connector_id']=201
            data.streets[1]['items'][1]['through_connectors']=rows
            with self.subTest(change=change),self.assertRaises(ValueError):
                m.PreparedPhysicalPorts(('residential',)).prepare(data,interior_ids=(21,))


def section(sid, points, category='residential'):
    return dict(kind='Section', id=sid, points=points, width=6., offset=0., category=category)


def fixture():
    def attachment(sid, sec, end, cid, category='residential'):
        return dict(street_id=sid, section_id=sec, end=end, connector_id=cid, category=category)
    a = dict(kind='Intersection', id=21, location=(0., 0.), through=((10, 1), (11, 0)),
             surface_treatment='NATIVE_JUNCTION', attachments=[attachment(2, 20, 1, 201)])
    b = dict(kind='Intersection', id=22, location=(40., 0.), through=((12, 1), (13, 0)),
             surface_treatment='NATIVE_JUNCTION', attachments=[attachment(3, 30, 0, 202)])
    items = [section(10, ((-40., 0.), (-20., .5), (0., 0.))), a,
             section(11, ((0., 0.), (15., 1.))),
             dict(kind='Intersection', id=23, surface_treatment='DEFERRED_ATTACHMENT'),
             section(12, ((15., 1.), (40., 0.))), b,
             section(13, ((40., 0.), (60., -.5), (80., 0.)))]
    streets = {1: dict(id=1, items=items, status='SUPPORTED', endpoints=(1001, 1002), bundle_id=None)}
    for sid, points, ends in (
        (2, ((0., 30.), (0., 0.)), (None, 21)),
        (3, ((40., 0.), (40., -30.)), (22, None)),
        (4, ((-40., 0.), (-60., -20.)), (1001, None)),
        (5, ((-40., 0.), (-60., 20.)), (1001, None)),
        (6, ((80., 0.), (100., 20.)), (1002, None)),
        (7, ((80., 0.), (100., -20.)), (1002, None)),
    ):
        streets[sid] = dict(id=sid, items=[section(sid*10, points)], status='SUPPORTED', endpoints=ends, bundle_id=None)
    sites = []
    for jid, location, refs in ((1001, (-40., 0.), ((1, 10, 0), (4, 40, 0), (5, 50, 0))),
                                (1002, (80., 0.), ((1, 13, 1), (6, 60, 0), (7, 70, 0)))):
        sites.append(dict(junction_id=jid, location=location, ring=[dict(
            item_kind='Street', connector_id=10*jid+i,
            endpoints=[dict(street_id=s, terminal_id=sec, terminal_kind='Section', end=end)])
            for i, (s, sec, end) in enumerate(refs)]))
    return NS(streets=streets, bundles={}, junction_origins={'physical': sites})


def mixed_fixture():
    snapshot = fixture()
    snapshot.streets[99] = dict(id=99, status='UNSUPPORTED',
        reason='Category footway is outside the carriageway scope',
        items=[section(990, ((-40., 0.), (-40., 20.)), 'footway')],
        endpoints=(1001, None), bundle_id=None)
    snapshot.junction_origins['physical'][0]['ring'].insert(1, dict(
        item_kind='Street', connector_id=999,
        endpoints=[dict(street_id=99, terminal_id=990, terminal_kind='Section', end=0)]))
    return snapshot


def scoped_fixture():
    data = fixture()
    source = data.streets[4]
    source.update(status='UNSUPPORTED', reason='Reserved component/category: footway')
    source['items'].extend([
        dict(kind='SymLane', id=5, points=(), surface_treatment='WIDTH_TRANSITION',
             transition=dict(source_kind='SymLane', source_id=5, before_section=40,
                             after_section=41, minimum_half_length=1.5, maximum_edge_slope=.3)),
        dict(section(41, ((-60., -20.), (-75., -25.)), 'footway'), width=2.)])
    return data


def seam_fixture():
    data=fixture();location=(-60.,-20.)
    data.streets[8]=dict(id=8,items=[section(80,(location,(-80.,-20.)))],
                        status='SUPPORTED',endpoints=(None,None),bundle_id=None)
    data.junction_origins['physical'].append(dict(junction_id=2000,location=location,ring=[]))
    data.junction_origins['seams']=[dict(junction_id=2000,location=location,endpoints=[
        dict(street_id=4,terminal_id=40,terminal_kind='Section',end=1),
        dict(street_id=8,terminal_id=80,terminal_kind='Section',end=0)])]
    return data


class PhysicalPortTests(unittest.TestCase):
    def setUp(self):
        self.reader = m.PreparedPhysicalPorts(('residential',))

    def test_cross_street_category_boundary_retains_original_provenance(self):
        for own_end in (0,1):
            for peer_end in (0,1):
                data=seam_fixture();data.streets[8]['items'][0]['category']='footway'
                data.streets[8].update(status='UNSUPPORTED',reason='Reserved component/category: footway')
                if own_end==0:
                    data.streets[4]['items'][0]['points']=data.streets[4]['items'][0]['points'][::-1]
                    data.streets[4]['endpoints']=data.streets[4]['endpoints'][::-1]
                    data.junction_origins['physical'][0]['ring'][1]['endpoints'][0]['end']=1
                    data.junction_origins['seams'][0]['endpoints'][0]['end']=0
                if peer_end==1:
                    data.streets[8]['items'][0]['points']=data.streets[8]['items'][0]['points'][::-1]
                    data.junction_origins['seams'][0]['endpoints'][1]['end']=1
                before=deepcopy(data)
                reader=m.PreparedPhysicalPorts(('residential',),section_scopes=True,seam_connections=True,category_connections=True)
                plan=reader.prepare(data,physical_ids=(1001,))
                self.assertEqual(reader.category_seam_ids(data),{2000})
                span=next(s for s in plan.spans if s.street_id==4)
                end=span.start if own_end==0 else span.end
                self.assertEqual((end.owner_kind,end.junction_id),('CATEGORY_BOUNDARY',None))
                boundary=plan.category_boundary(end.endpoint)
                self.assertEqual(boundary.source_id,2000)
                self.assertEqual(boundary.endpoints,(m.SectionEnd(4,40,own_end),m.SectionEnd(8,80,peer_end)))
                self.assertEqual(boundary.categories,('residential','footway'))
                self.assertNotIn(8,plan.source_streets)
                self.assertEqual(data,before)

    def test_cross_street_category_boundary_does_not_bypass_bad_provenance(self):
        reader=m.PreparedPhysicalPorts(('residential',),section_scopes=True,seam_connections=True,category_connections=True)
        for mutation in ('missing','location','terminal','category','duplicate','physical','grade'):
            data=seam_fixture();data.streets[8]['items'][0]['category']='footway'
            seam=data.junction_origins['seams'][0]
            if mutation=='missing':del data.streets[8]
            elif mutation=='location':seam['location']=(0.,0.)
            elif mutation=='terminal':seam['endpoints'][1]['end']=1
            elif mutation=='category':data.streets[8]['items'][0]['category']=None
            elif mutation=='duplicate':data.junction_origins['seams'].append(deepcopy(seam))
            elif mutation=='physical':data.junction_origins['physical'][-1]['ring']=[dict(endpoints=[])]
            elif mutation=='grade':data.streets[4].update(status='UNSUPPORTED',reason='Grade/elevation treatment is not yet implemented')
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):
                reader.prepare(data,physical_ids=(1001,))

    def test_cross_street_category_boundaries_require_capability(self):
        for scopes,seams in ((False,False),(True,False),(False,True)):
            with self.assertRaises(ValueError):
                m.PreparedPhysicalPorts(('residential',),section_scopes=scopes,seam_connections=seams,category_connections=True)
        data=seam_fixture();data.streets[8]['items'][0]['category']='footway'
        plan=m.PreparedPhysicalPorts(('residential',),section_scopes=True,seam_connections=True).prepare(data,physical_ids=(1001,))
        self.assertEqual(next(s for s in plan.spans if s.street_id==4).end.owner_kind,'SEAM')

    def test_selected_seam_needs_capability_and_retains_two_exact_original_ports(self):
        data=seam_fixture();before=deepcopy(data)
        with self.assertRaisesRegex(ValueError,'native connection provider'):
            self.reader.prepare(data,seam_ids=(2000,))
        reader=m.PreparedPhysicalPorts(('residential',),seam_connections=True)
        plan=reader.prepare(data,physical_ids=(1001,),seam_ids=(2000,))
        site=next(s for s in plan.sites if s.junction_id==2000)
        self.assertEqual([(p.endpoint.street_id,p.endpoint.end,p.source_kind,p.connector_id) for p in site.ports],
                         [(4,1,'SEAM',None),(8,0,'SEAM',None)])
        self.assertEqual(data,before)
        self.assertEqual(plan.span_for(site.ports[0]).end.owner_kind,'SEAM')
        self.assertEqual(plan.span_for(site.ports[1]).start.owner_kind,'SEAM')
        with self.assertRaisesRegex(ValueError,'ordinary junction'):
            reader.prepare(data,physical_ids=(2000,),seam_ids=(2000,))
        with self.assertRaisesRegex(ValueError,'at least three'):
            reader.prepare(data,physical_ids=(2000,))
        data.streets[8]['status']='UNSUPPORTED'
        with self.assertRaisesRegex(ValueError,'unsupported'):reader.prepare(data,seam_ids=(2000,))

    def test_original_two_ended_seam_is_never_an_independent_free_end(self):
        data=fixture();location=(-60.,-20.)
        data.streets[8]=dict(id=8,items=[section(80,(location,(-80.,-20.)))],
                            status='UNSUPPORTED',reason='Grade/elevation treatment',endpoints=(None,None))
        data.junction_origins['physical'].append(dict(junction_id=2000,location=location,ring=[]))
        data.junction_origins['seams']=[dict(junction_id=2000,location=location,endpoints=[
            dict(street_id=4,terminal_id=40,terminal_kind='Section',end=1),
            dict(street_id=8,terminal_id=80,terminal_kind='Section',end=0)])]
        before=deepcopy(data)
        plan=self.reader.prepare(data,physical_ids=(1001,))
        span=next(s for s in plan.spans if s.street_id==4)
        self.assertEqual((span.end.owner_kind,span.end.junction_id),('SEAM',2000))
        self.assertEqual(data,before)
        self.assertNotIn(8,plan.source_streets)
        for damage in ('terminal','location','free','repeated','physical'):
            with self.subTest(damage=damage):
                bad=deepcopy(data);seam=bad.junction_origins['seams'][0]
                if damage=='terminal':seam['endpoints'][1]['terminal_id']=999
                elif damage=='location':bad.streets[8]['items'][0]['points']=((0.,0.),(-80.,-20.))
                elif damage=='free':bad.junction_origins['free_endpoints']=[dict(
                    endpoint=seam['endpoints'][0],location=location)]
                elif damage=='repeated':seam['endpoints'][1]=seam['endpoints'][0]
                else:bad.junction_origins['physical'].pop()
                with self.assertRaises(ValueError):self.reader.prepare(bad,physical_ids=(1001,))

    def test_proven_original_free_end_overrides_later_logical_bundle_link(self):
        data=fixture()
        data.streets[4]['endpoints']=(1001,999)
        original=self.reader.prepare(data,physical_ids=(1001,))
        self.assertEqual(next(s for s in original.spans if s.street_id==4).end.owner_kind,'REMOTE_LOGICAL')
        record=dict(endpoint=dict(street_id=4,terminal_id=40,terminal_kind='Section',end=1),
                    location=(-60.,-20.))
        data.junction_origins['free_endpoints']=[record]
        before=deepcopy(data)
        plan=self.reader.prepare(data,physical_ids=(1001,))
        self.assertEqual(data,before)
        span=next(s for s in plan.spans if s.street_id==4)
        self.assertEqual(span.start.junction_id,1001)
        self.assertEqual((span.end.owner_kind,span.end.junction_id),('FREE',None))
        self.assertEqual(span.section_ids,(40,))
        self.assertEqual([s for s in plan.spans if s.street_id!=4],
                         [s for s in original.spans if s.street_id!=4])

    def test_original_free_end_cannot_change_identity_location_or_existing_owner(self):
        for change in ('duplicate','section','end','location','kind','physical','interior'):
            with self.subTest(change=change):
                data=fixture()
                record=dict(endpoint=dict(street_id=4,terminal_id=40,terminal_kind='Section',end=1),
                            location=(-60.,-20.))
                data.junction_origins['free_endpoints']=[record]
                if change=='duplicate':data.junction_origins['free_endpoints'].append(deepcopy(record))
                elif change=='section':record['endpoint']['terminal_id']=999
                elif change=='end':record['endpoint']['end']=2
                elif change=='location':record['location']=(-60.,-21.)
                elif change=='kind':record['endpoint']['terminal_kind']='SideLane'
                elif change=='physical':record['endpoint']['end']=0;record['location']=(-40.,0.)
                else:
                    record['endpoint'].update(street_id=2,terminal_id=20)
                    record['location']=(0.,0.)
                with self.assertRaises(ValueError):
                    self.reader.prepare(data,physical_ids=(1001,),interior_ids=(21,))

    def test_scoped_road_keeps_complete_street_and_declared_category_seam(self):
        data = scoped_fixture(); before = deepcopy(data)
        reader = m.PreparedPhysicalPorts(('residential',), section_scopes=True)
        plan = reader.prepare(data, physical_ids=(1001,))
        self.assertEqual(data, before)
        road, path = [s for s in plan.spans if s.street_id == 4]
        self.assertTrue(road.in_scope); self.assertFalse(path.in_scope)
        self.assertEqual((road.end.owner_kind, road.end.junction_id), ('CATEGORY_BOUNDARY', None))
        self.assertEqual(path.start.owner_kind, 'CATEGORY_BOUNDARY')
        self.assertEqual(path.end.owner_kind, 'FREE')
        self.assertEqual((road.first_control, road.last_control, road.source_control_count), (0, 1, 3))
        self.assertEqual(plan.category_boundary(road.end.endpoint),
                         m.CategoryBoundary(4, 40, 41, 'SymLane', 5))
        self.assertTrue(plan.source_supported(data.streets[4]))
        self.assertEqual(plan.section_scopes[0].curve_groups,(0,1))

    def test_curve_groups_join_consecutive_owned_sections_across_junction_spans(self):
        scope=m.StreetSectionScope(1,('residential',),(10,11,12,13,14,15),
                                  (10,11,14,15),())
        self.assertEqual(scope.curve_groups,(0,0,1,1,2,2))

    def test_scoped_category_policy_is_explicit_and_does_not_hide_other_rejections(self):
        with self.assertRaises(ValueError):self.reader.prepare(scoped_fixture(), physical_ids=(1001,))
        for change in ('grade', 'width', 'transition_id', 'missing_transition', 'extent', 'other_reason'):
            with self.subTest(change=change):
                data = scoped_fixture(); road, seam, path = data.streets[4]['items']
                if change == 'grade':road['tags'] = {'bridge':'yes'}
                elif change == 'width':road['width'] = float('nan')
                elif change == 'transition_id':seam['transition']['source_id'] += 1
                elif change == 'missing_transition':del seam['transition']
                elif change == 'extent':seam['points'] = [(1., 2.), (2., 3.)]
                else:data.streets[4]['reason'] = 'Disconnected source'
                with self.assertRaises(ValueError):
                    m.PreparedPhysicalPorts(('residential',), section_scopes=True).prepare(data, physical_ids=(1001,))

    def test_scope_preserves_opposite_authored_direction_and_original_ids(self):
        data = scoped_fixture(); source = data.streets[4]
        source['items'].reverse(); source['endpoints'] = source['endpoints'][::-1]
        for item in source['items']:
            if item['kind'] == 'Section':item['points'] = tuple(reversed(item['points']))
            else:item['transition'].update(before_section=41, after_section=40)
        data.junction_origins['physical'][0]['ring'][1]['endpoints'][0]['end'] = 1
        plan = m.PreparedPhysicalPorts(('residential',), section_scopes=True).prepare(data, physical_ids=(1001,))
        path, road = [s for s in plan.spans if s.street_id == 4]
        self.assertFalse(path.in_scope); self.assertTrue(road.in_scope)
        self.assertEqual((road.start.owner_kind, road.end.junction_id), ('CATEGORY_BOUNDARY', 1001))
        self.assertEqual(plan.section_scopes[0].included_sections, (40,))

    def test_two_interiors_partition_complete_street_once(self):
        snapshot = fixture()
        before = deepcopy(snapshot)
        plan = self.reader.prepare(snapshot, interior_ids=(21, 22))
        self.assertEqual(snapshot, before)
        self.assertEqual(plan.source_streets, (1, 2, 3))
        main = [s for s in plan.spans if s.street_id == 1]
        self.assertEqual([s.section_ids for s in main], [(10,), (11, 12), (13,)])
        self.assertEqual([(s.start.junction_id, s.end.junction_id) for s in main], [(1001, 21), (21, 22), (22, 1002)])
        self.assertEqual([s.identity for s in main], [(1, 10, 10), (1, 11, 12), (1, 13, 13)])
        self.assertEqual([(s.first_control, s.last_control, s.source_control_count) for s in main],
                         [(0, 2, 7), (2, 4, 7), (4, 6, 7)])
        self.assertEqual(len(plan.sites[0].ports), 3)
        self.assertEqual({p.endpoint for p in plan.sites[0].ports},
                         {m.SectionEnd(1, 10, 1), m.SectionEnd(1, 11, 0), m.SectionEnd(2, 20, 1)})
        for site in plan.sites:
            for port in site.ports:
                self.assertEqual(plan.span_for(port).street_id, port.endpoint.street_id)
        through = [p for s in plan.sites for p in s.ports if p.source_kind == 'THROUGH']
        self.assertTrue(all(p.connector_id is None for p in through))
        self.assertEqual({p.connector_id for s in plan.sites for p in s.ports if p.source_kind == 'ATTACHMENT'}, {201, 202})

    def test_original_sites_and_native_sites_share_section_identity(self):
        plan = self.reader.prepare(fixture(), physical_ids=(1001, 1002), interior_ids=(21, 22))
        self.assertEqual(plan.source_streets, tuple(range(1, 8)))
        self.assertEqual([p.connector_id for p in plan.sites[0].ports], [10010, 10011, 10012])
        head = next(p for p in plan.sites[0].ports if p.endpoint.street_id == 1)
        tail = next(p for p in plan.sites[2].ports if p.endpoint == m.SectionEnd(1, 10, 1))
        self.assertEqual(plan.span_for(head), plan.span_for(tail))
        self.assertEqual(len([s for s in plan.spans if s.start.owner_kind == 'FREE' or s.end.owner_kind == 'FREE']), 6)

    def test_unselected_interior_stays_an_explicit_span_boundary(self):
        plan = self.reader.prepare(fixture(), physical_ids=(1001,))
        spans = [s for s in plan.spans if s.street_id == 1]
        self.assertEqual(len(spans), 3)
        self.assertEqual(spans[0].end.owner_kind, 'INTERIOR')
        self.assertEqual(spans[0].end.junction_id, 21)
        self.assertEqual(spans[-1].end.owner_kind, 'PHYSICAL')
        self.assertEqual(spans[-1].end.junction_id, 1002)

    def test_geometry_edits_and_object_metadata_do_not_reassign_ports(self):
        snapshot = fixture()
        before = self.reader.prepare(snapshot, physical_ids=(1001, 1002), interior_ids=(21, 22))
        snapshot.streets[1]['items'][0]['points'] = ((-40., 0.), (-20., 1.1), (0., 0.))
        for sid, street in snapshot.streets.items():
            street['source_object'] = 'Grouped source' if sid % 2 else 'Whole Street '+str(sid)
            for section in street['items']:
                if section['kind'] == 'Section':
                    section['width'] += .3
        self.assertEqual(before, self.reader.prepare(snapshot, physical_ids=(1001, 1002), interior_ids=(21, 22)))

    def test_reversed_authoring_preserves_all_original_section_owners(self):
        snapshot = fixture()
        source = snapshot.streets[1]
        source['items'].reverse()
        source['endpoints'] = source['endpoints'][::-1]
        for item in source['items']:
            if item['kind'] == 'Section':
                item['points'] = tuple(reversed(item['points']))
            elif item.get('through'):
                item['through'] = tuple((sid, 1-end) for sid, end in reversed(item['through']))
        for site in snapshot.junction_origins['physical']:
            site['ring'][0]['endpoints'][0]['end'] = 1-site['ring'][0]['endpoints'][0]['end']
        plan = self.reader.prepare(snapshot, physical_ids=(1001, 1002), interior_ids=(21, 22))
        spans = [s for s in plan.spans if s.street_id == 1]
        self.assertEqual([s.section_ids for s in spans], [(13,), (12, 11), (10,)])
        self.assertEqual([(s.start.junction_id, s.end.junction_id) for s in spans], [(1002, 22), (22, 21), (21, 1001)])

    def test_out_of_scope_attachment_is_retained_as_deferred(self):
        snapshot = fixture()
        item = snapshot.streets[1]['items'][1]
        item['attachments'].append(dict(street_id=99, section_id=990, end=0, connector_id=999, category='footway'))
        plan = self.reader.prepare(snapshot, interior_ids=(21,))
        self.assertEqual(len(plan.sites[0].ports), 3)
        self.assertEqual(plan.deferred, (m.DeferredPort(21, m.SectionEnd(99, 990, 0), 999, 'footway'),))

    def test_physical_category_exclusion_retains_identity_order_and_sources(self):
        for status in ('SUPPORTED', 'UNSUPPORTED'):
            snapshot = mixed_fixture()
            snapshot.streets[99]['status'] = status
            before = deepcopy(snapshot)
            plan = self.reader.prepare(snapshot, physical_ids=(1001,))
            expected = self.reader.prepare(fixture(), physical_ids=(1001,))
            self.assertEqual(snapshot, before)
            self.assertEqual(plan.sites, expected.sites)
            self.assertEqual(plan.spans, expected.spans)
            self.assertEqual(plan.source_streets, expected.source_streets)
            self.assertEqual(plan.deferred,
                (m.DeferredPort(1001, m.SectionEnd(99, 990, 0), 999, 'footway'),))
        # An expanded scope must require a supported source for the same arm.
        reader = m.PreparedPhysicalPorts(('residential', 'footway'))
        with self.assertRaisesRegex(ValueError, 'Street 99 is unsupported'):
            reader.prepare(mixed_fixture(), physical_ids=(1001,))
        snapshot.streets[99]['status'] = 'SUPPORTED'
        plan = reader.prepare(snapshot, physical_ids=(1001,))
        self.assertEqual([p.connector_id for p in plan.sites[0].ports], [10010, 999, 10011, 10012])
        self.assertFalse(plan.deferred)

    def test_deferred_physical_port_still_requires_authentic_source_incidence(self):
        for mutation in ('missing', 'location', 'terminal', 'nonfinite', 'id', 'category',
                         'duplicate_connector', 'duplicate_endpoint'):
            snapshot = mixed_fixture()
            source = snapshot.streets[99]
            ring = snapshot.junction_origins['physical'][0]['ring']
            if mutation == 'missing':
                del snapshot.streets[99]
            elif mutation == 'location':
                source['items'][0]['points'] = ((-39., 0.), (-40., 20.))
            elif mutation == 'terminal':
                ring[1]['endpoints'][0]['terminal_id'] = 991
            elif mutation == 'nonfinite':
                source['items'][0]['points'] = ((-40., 0.), (float('nan'), 20.))
            elif mutation == 'id':
                source['id'] = 98
            elif mutation == 'category':
                source['items'][0]['category'] = None
            elif mutation == 'duplicate_connector':
                ring[1]['connector_id'] = ring[0]['connector_id']
            else:
                extra = deepcopy(ring[1]); extra['connector_id'] = 998
                ring.append(extra)
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                self.reader.prepare(snapshot, physical_ids=(1001,))

    def test_category_exclusion_does_not_invent_two_arm_junction_or_hide_grade(self):
        snapshot = mixed_fixture()
        snapshot.streets[5]['items'][0]['category'] = 'footway'
        with self.assertRaisesRegex(ValueError, 'at least three in-scope ports'):
            self.reader.prepare(snapshot, physical_ids=(1001,))
        snapshot = mixed_fixture()
        snapshot.streets[5].update(status='UNSUPPORTED', reason='Grade treatment is not implemented')
        with self.assertRaisesRegex(ValueError, 'Street 5 is unsupported: Grade'):
            self.reader.prepare(snapshot, physical_ids=(1001,))

    def test_missing_source_or_malformed_through_references_reject(self):
        for mutation in ('absent', 'unsupported', 'nonconsecutive', 'duplicate', 'location', 'direction', 'terminal'):
            snapshot = fixture()
            item = snapshot.streets[1]['items'][1]
            if mutation == 'absent':
                del snapshot.streets[2]
            elif mutation == 'unsupported':
                snapshot.streets[2]['status'] = 'UNSUPPORTED'
            elif mutation == 'nonconsecutive':
                item['through'] = ((10, 1), (12, 0))
            elif mutation == 'duplicate':
                item['through'] = ((10, 1), (10, 1))
            elif mutation == 'location':
                item['location'] = (.01, 0.)
            elif mutation == 'direction':
                snapshot.streets[2]['items'][0]['points'] = ((0., 0.), (0., 0.))
            else:
                item['attachments'][0]['section_id'] = 10
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                self.reader.prepare(snapshot, interior_ids=(21,))

    def test_duplicate_site_record_and_conflicting_ownership_reject(self):
        snapshot = fixture()
        snapshot.junction_origins['physical'].append(deepcopy(snapshot.junction_origins['physical'][0]))
        with self.assertRaisesRegex(ValueError, 'one original ring'):
            self.reader.prepare(snapshot, physical_ids=(1001,))
        snapshot = fixture()
        fake = deepcopy(snapshot.junction_origins['physical'][0])
        fake['junction_id'] = 9999
        snapshot.junction_origins['physical'].append(fake)
        with self.assertRaisesRegex(ValueError, 'conflicting'):
            self.reader.prepare(snapshot, interior_ids=(21,))

    def test_duplicate_section_owner_and_disconnected_sources_reject(self):
        snapshot = fixture()
        snapshot.streets[1]['items'][4]['points'] = ((16., 1.), (40., 0.))
        with self.assertRaisesRegex(ValueError, 'disconnected'):
            self.reader.prepare(snapshot, interior_ids=(21,))
        snapshot = fixture()
        snapshot.streets[2]['items'].append(deepcopy(snapshot.streets[2]['items'][0]))
        with self.assertRaisesRegex(ValueError, 'duplicate Sections'):
            self.reader.prepare(snapshot, interior_ids=(21,))

    def test_selection_and_provider_scope_are_explicit(self):
        for kwargs in ({}, {'interior_ids': (21, 21)}, {'physical_ids': (True,)},
                       {'physical_ids': (21,), 'interior_ids': (21,)}, {'interior_ids': (999,)}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.reader.prepare(fixture(), **kwargs)
        with self.assertRaises(ValueError):
            m.PreparedPhysicalPorts(())

    def test_ids_and_section_ends_are_not_coerced(self):
        for args in ((True, 10, 0), (1, -10, 0), (1, 10, True), (1, 10, 1.0)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                m.SectionEnd(*args)
        for value in (True, 1.0, -1):
            snapshot = fixture()
            snapshot.streets[1]['items'][1]['through'] = ((10, value), (11, 0))
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.reader.prepare(snapshot, interior_ids=(21,))
        snapshot = fixture()
        snapshot.streets[1]['id'] = 99
        with self.assertRaisesRegex(ValueError, 'inventory key'):
            self.reader.prepare(snapshot, interior_ids=(21,))

    def test_duplicate_remote_owners_and_connectors_are_rejected(self):
        snapshot = fixture()
        snapshot.junction_origins['physical'].append(deepcopy(snapshot.junction_origins['physical'][0]))
        with self.assertRaisesRegex(ValueError, 'duplicate or conflicting'):
            self.reader.prepare(snapshot, interior_ids=(21,))
        snapshot = fixture()
        ring = snapshot.junction_origins['physical'][0]['ring']
        ring[1]['connector_id'] = ring[0]['connector_id']
        with self.assertRaisesRegex(ValueError, 'repeats an original Connector'):
            self.reader.prepare(snapshot, physical_ids=(1001,))

    def test_managed_control_insertion_rebinds_indices_without_new_source_ids(self):
        snapshot = fixture()
        before = self.reader.prepare(snapshot, interior_ids=(21, 22))
        snapshot.streets[1]['items'][0]['points'] = ((-40., 0.), (-30., .25), (-20., .5), (0., 0.))
        after = self.reader.prepare(snapshot, interior_ids=(21, 22))
        self.assertEqual([s.identity for s in before.spans], [s.identity for s in after.spans])
        self.assertEqual(before.sites, after.sites)
        main = [s for s in after.spans if s.street_id == 1]
        self.assertEqual([(s.first_control, s.last_control, s.source_control_count) for s in main],
                         [(0, 3, 8), (3, 5, 8), (5, 7, 8)])

    def test_port_designs_preserve_section_radii_and_bundle_side(self):
        from test_physical_junctions import module as design
        snapshot = fixture()
        snapshot.streets[1]['bundle_id'] = 7
        snapshot.bundles[7] = dict(topology=dict(status='MAPPED', members=[
            dict(street_id=1, complete=True, head_order=0, tail_order=0, forward=True),
            dict(street_id=99, complete=True, head_order=1, tail_order=1, forward=True)]))
        snapshot.streets[1]['items'][2]['category'] = 'tertiary'
        reader = m.PreparedPhysicalPorts(('residential', 'tertiary'))
        plan = reader.prepare(snapshot, interior_ids=(21, 22))
        designer = design.PreparedPhysicalJunctions(None, None,
            dict(residential=dict(radius=2.), tertiary=dict(radius=7.)))
        result = designer.prepare_port_designs(snapshot, plan)
        self.assertEqual([(p.section_id, p.end, p.half_profile) for p in result if p.street_id==1],
                         [(10, 1, -1), (11, 0, 1), (12, 1, -1), (13, 0, 1)])
        self.assertEqual(next(p.radius for p in result if p.section_id==11), 7.)
        self.assertTrue(all(p.rise==.08 for p in result))
        self.assertTrue(all(p.bundle_id == (7 if p.street_id == 1 else -1) for p in result))
        del designer.category_properties['tertiary']
        with self.assertRaisesRegex(ValueError, 'declared radius'):
            designer.prepare_port_designs(snapshot, plan)

    def test_default_carriageway_scope_includes_living_streets(self):
        spec=importlib.util.spec_from_file_location('tested_street_design',
            Path(__file__).resolve().parents[1]/'way/street_network/design.py')
        design=importlib.util.module_from_spec(spec);sys.modules[spec.name]=design;spec.loader.exec_module(design)
        source=fixture()
        source.streets[2]['items'][0]['category']='living_street'
        source.streets[1]['items'][1]['attachments'][0]['category']='living_street'
        before=deepcopy(source)
        plan=m.PreparedPhysicalPorts(design.StreetDesign().categories).prepare(source,interior_ids=(21,))
        self.assertEqual(source,before)
        self.assertEqual(len(plan.sites[0].ports),3)
        self.assertFalse(plan.deferred)
        # A caller can still deliberately defer this category for a narrower provider.
        with self.assertRaisesRegex(ValueError,'at least three'):
            self.reader.prepare(source,interior_ids=(21,))


if __name__ == '__main__':
    unittest.main()
