"""Packaged integration smoke/regression tests; run in background Blender 5.2."""
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
import bpy

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import blosm
from blosm.app.asset_library import BlenderAssetLoader
from blosm.way.street_network import assets
from blosm.way.street_network.adapter import NetworkSnapshot, StreetDesign, PreparedNetworkAdapter
from blosm.way.street_network.authoring import EditableNetworkBuilder
from blosm.way.street_network.minor import MinorJunctionAdapter
from blosm.way.street_network.transaction import ImportTransaction
from blosm.way.gn.editable_input import model, runtime
from blosm.way.gn.kernel import inputs


def sample():
    snapshot = NetworkSnapshot()
    points = [((0.,0.,0.),(0.,12.,0.),(0.,30.,0.)),
              ((0.,0.,0.),(-12.,0.,0.),(-30.,0.,0.)),
              ((0.,0.,0.),(12.,0.,0.),(30.,0.,0.))]
    sources, connectors = [], []
    for i,p in enumerate(points):
        section = dict(kind='Section',id=i,points=p,width=6.,offset=0.,category='residential',
                       tags={},source_way_ids=[str(90000000000+i)])
        snapshot.streets[i] = dict(id=i, items=[section],status='SUPPORTED',reason=None,
                                   endpoints=(0,None),bundle_id=None)
        sources.append(inputs.StreetInput(i, (inputs.SectionInput(i,p,6.),)))
        connectors.append(inputs.ConnectorInput(i,(inputs.MemberInput(i,True),),i,radius=2.))
    # Connectors must be ordered counterclockwise from +X.
    connectors = [connectors[2],connectors[0],connectors[1]]
    snapshot.junctions[0] = inputs.PreparedJunction(0,(0.,0.),tuple(sources),tuple(connectors),
        (inputs.ProfileInput(.08),)*3,0.,.7)
    snapshot.intersections[0] = dict(id=0,order=3,status='SUPPORTED')
    return snapshot


def divided_sample(widths=(4.,8.), reverse_second=False):
    """Synthetic T: independently sized parallel members join two side roads."""
    snapshot=NetworkSnapshot()
    paths=[((-5.,0.,0.),(-5.,15.,0.),(-5.,35.,0.)),
           ((5.,0.,0.),(5.,15.,0.),(5.,35.,0.)),
           ((0.,0.,0.),(-20.,0.,0.),(-40.,0.,0.)),
           ((0.,0.,0.),(20.,0.,0.),(40.,0.,0.))]
    if reverse_second:paths[1]=tuple(reversed(paths[1]))
    streets=[]
    for sid,(points,width) in enumerate(zip(paths,(*widths,6.,6.))):
        section=dict(kind='Section',id=sid,points=points,width=width,offset=0.,
                     category='residential',tags={},source_way_ids=[str(sid)])
        snapshot.streets[sid]=dict(id=sid,items=[section],status='SUPPORTED',reason=None,
            endpoints=(None,0) if sid==1 and reverse_second else (0,None),
            bundle_id=7 if sid<2 else None)
        streets.append(inputs.StreetInput(sid,(inputs.SectionInput(sid,points,width),)))
    members=(inputs.MemberInput(0,True,0,0,0,0),
             inputs.MemberInput(1,not reverse_second,1,0,1,1))
    profiles=(inputs.ProfileInput(.08,separator=True),inputs.ProfileInput(.08),inputs.ProfileInput(.08))
    connectors=(inputs.ConnectorInput(2,(inputs.MemberInput(3,True),),2),
                inputs.ConnectorInput(0,members,0,item_kind=inputs.BUNDLE_ITEM,bundle_source_id=7),
                inputs.ConnectorInput(1,(inputs.MemberInput(2,True),),1))
    snapshot.junctions[0]=inputs.PreparedJunction(0,(0.,0.),tuple(streets),connectors,profiles,0.,.7)
    snapshot.intersections[0]=dict(id=0,order=3,status='SUPPORTED')
    return snapshot


class StreetIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        model.register()
        cls.loader = BlenderAssetLoader(Path(blosm.__file__).parent/'assets', assets.catalog())
        loaded = cls.loader.load_many((assets.SOLVER,assets.PAVEMENT,assets.GRASS))
        cls.blocks = dict(solver=loaded[assets.SOLVER],pavement=loaded[assets.PAVEMENT],grass=loaded[assets.GRASS])

    def test_registration_is_repeatable(self):
        model.register(); model.register()
        self.assertTrue(model.ST_NetworkRecord.is_registered)

    def test_closed_multi_section_street_is_explicitly_outside_path_scope(self):
        first=dict(sample().streets[0]['items'][0],points=((0.,0.,0.),(30.,0.,0.)))
        second=dict(first,id=99,points=((30.,0.,0.),(15.,20.,0.),(0.,0.,0.)))
        self.assertEqual(PreparedNetworkAdapter().street_support([first,second]),
                         'Closed Street requires a dedicated path treatment')

    def test_duplicate_section_ownership_is_diagnosed_before_adaptation(self):
        from copy import deepcopy
        from blosm.way.street_network.ownership import PreparedOwnershipReview
        for same_street in (False,True):
            snapshot=sample()
            if same_street:
                snapshot.streets[0]['items'].append(deepcopy(snapshot.streets[0]['items'][0]))
                owners=[0,0]
            else:
                snapshot.streets[40]=deepcopy(snapshot.streets[0]);snapshot.streets[40]['id']=40
                owners=[0,40]
            PreparedOwnershipReview().apply(snapshot)
            self.assertEqual(snapshot.ownership_audit['duplicate_sections'],[dict(section_id=0,street_ids=owners)])
            self.assertEqual(snapshot.issues[-1]['code'],'DUPLICATE_SECTION_OWNERSHIP')
            self.assertTrue(all(snapshot.streets[sid]['status']=='UNSUPPORTED' for sid in owners))
            self.assertEqual(snapshot.streets[1]['status'],'SUPPORTED')

    def test_bundle_member_review_requires_declared_dimension_transitions(self):
        from copy import deepcopy
        from blosm.way.street_network.bundles import PreparedBundleReview
        review=PreparedBundleReview()
        row=deepcopy(divided_sample().streets[0])
        first=row['items'][0];second=dict(first,id=100,points=first['points'][1:])
        first['points']=first['points'][:2]
        row['items']=[first,dict(kind='Intersection',surface_treatment='DEFERRED_ATTACHMENT'),second]
        self.assertIsNone(review.member_reason(row))
        second['width']+=1.
        self.assertIn('transition provider',review.member_reason(row))
        from blosm.way.street_network.transitions import PreparedTransitionAdapter
        row['items'][1]=dict(kind='SideLane',id=4,points=())
        PreparedTransitionAdapter().read(row['items'])
        self.assertIsNone(review.member_reason(row))
        row['items'][1]['transition']['after_section']=999
        self.assertIn('transition provider',review.member_reason(row))
        row['items'][1]['transition']['after_section']=second['id']
        points=second['points'];second['points']=((999.,999.,0.),*points[1:])
        self.assertIn('continuous directed chain',review.member_reason(row))
        second['points']=points
        second['width']=first['width']
        row['items'][1]['surface_treatment']='NATIVE_JUNCTION'
        self.assertIn('member spans',review.member_reason(row))

    def test_minor_attachment_street_is_resolved_by_unique_section_ownership(self):
        from blosm.way.street_network.ownership import PreparedOwnershipReview
        snapshot=sample()
        item=dict(kind='Intersection',id=9,location=(0.,0.),attachments=[
            dict(street_id=999,section_id=1,end=0,connector_id=77,category='residential')])
        snapshot.streets[0]['items'].insert(0,item)
        original=snapshot.streets[1]['items'][0]['points']
        PreparedOwnershipReview().apply(snapshot)
        self.assertEqual(item['attachments'][0]['street_id'],1)
        self.assertEqual(snapshot.streets[1]['items'][0]['points'],original)
        self.assertEqual(snapshot.ownership_audit['remapped_attachments'],[
            dict(junction_id=9,connector_id=77,section_id=1,end=0,prepared_street_id=999,street_id=1)])
        self.assertFalse(snapshot.ownership_audit['unresolved_attachments'])

    def test_minor_attachment_remap_does_not_guess_a_nearby_or_ambiguous_owner(self):
        from copy import deepcopy
        from blosm.way.street_network.ownership import PreparedOwnershipReview
        for mismatch in ('location','duplicate','missing'):
            snapshot=sample()
            item=dict(kind='Intersection',id=9,location=(0.,0.),attachments=[
                dict(street_id=999,section_id=1,end=0,connector_id=77,category='residential')])
            snapshot.streets[0]['items'].insert(0,item)
            if mismatch=='location':item['location']=(.01,0.)
            if mismatch=='duplicate':snapshot.streets[3]=deepcopy(snapshot.streets[1]);snapshot.streets[3]['id']=3
            if mismatch=='missing':item['attachments'][0]['section_id']=999
            PreparedOwnershipReview().apply(snapshot)
            self.assertEqual(item['attachments'][0]['street_id'],999)
            self.assertTrue(item['attachment_error'])
            self.assertEqual(len(snapshot.ownership_audit['unresolved_attachments']),1)
            self.assertFalse(snapshot.ownership_audit['remapped_attachments'])

    def test_multi_section_bundle_retains_semantic_faces_and_live_widths(self):
        from blosm.way.street_network.evidence import GeometryEvidence
        snapshot=divided_sample()
        snapshot.streets={sid:snapshot.streets[sid] for sid in (0,1)}
        snapshot.junctions.clear();snapshot.intersections.clear()
        for sid,street in snapshot.streets.items():
            first=street['items'][0]
            second=dict(first,id=100+sid,points=first['points'][1:])
            first['points']=first['points'][:2]
            street.update(items=[first,second],endpoints=(None,None))
        snapshot.bundles[7]=dict(id=7,head=(0,1),tail=(0,1),status='SUPPORTED',members=[
            dict(street_id=sid,forward=True,head_order=sid,tail_order=sid) for sid in (0,1)])
        collection=bpy.data.collections.new('Multi-Section Bundle regression');bpy.context.scene.collection.children.link(collection)
        obj=EditableNetworkBuilder(StreetDesign(),self.blocks).create('Complete members',snapshot,{0,1},set(),collection)
        reader=GeometryEvidence();baseline=reader.read(obj)
        self.assertEqual(baseline['diagnostics'],[])
        self.assertEqual(baseline['section_ids'],[0,1,100,101])
        source=next(o for o in obj.st_network.source_collection.objects if o.name.endswith('Street 0'))
        for row in source.data.attributes['st_width'].data:row.value=5.
        source.data.update();source.update_tag()
        self.assertEqual(reader.read(obj)['diagnostics'],[])
        source.data.attributes['st_width'].data[-1].value=4.
        source.data.update();source.update_tag()
        self.assertTrue(reader.read(obj)['diagnostics'])
        for row in source.data.attributes['st_width'].data:row.value=4.
        source.data.update();source.update_tag()
        self.assertEqual(reader.read(obj),baseline)

    def test_bundle_width_transition_updates_without_moving_the_end_cuts(self):
        from blosm.way.street_network.evidence import GeometryEvidence
        from blosm.way.street_network.transitions import PreparedTransitionAdapter
        from blosm.way.street_network.bundles import PreparedBundleReview
        snapshot=divided_sample()
        snapshot.streets={sid:snapshot.streets[sid] for sid in (0,1)}
        snapshot.junctions.clear();snapshot.intersections.clear()
        for street in snapshot.streets.values():
            street['endpoints']=(None,None)
            section=street['items'][0]
            section['points']=tuple((x,y*3.,z) for x,y,z in section['points'])
        first=snapshot.streets[0]['items'][0]
        second=dict(first,id=100,points=first['points'][1:],width=8.,offset=2.)
        first['points']=first['points'][:2]
        items=[first,dict(kind='SideLane',id=10,points=()),second]
        PreparedTransitionAdapter().read(items);snapshot.streets[0]['items']=items
        snapshot.bundles[7]=dict(id=7,head=(0,1),tail=(0,1),status='SUPPORTED',members=[
            dict(street_id=sid,forward=True,head_order=sid,tail_order=sid) for sid in (0,1)])
        PreparedBundleReview().apply(snapshot)
        self.assertEqual(snapshot.bundles[7]['status'],'SUPPORTED')
        collection=bpy.data.collections.new('Bundle dimension transition');bpy.context.scene.collection.children.link(collection)
        obj=EditableNetworkBuilder(StreetDesign(),self.blocks).create('Declared Bundle taper',snapshot,{0,1},set(),collection)
        reader=GeometryEvidence();baseline=reader.read(obj)
        self.assertFalse(baseline['diagnostics'],baseline)
        self.assertEqual(baseline['section_ids'],[0,1,100])
        self.assertEqual(baseline['transition_ids'],[0])
        source=next(o for o in obj.st_network.source_collection.objects if o.name.endswith('Street 0'))
        original=source.data.attributes['st_width'].data[-1].value
        for changed in (9.,200.,original):
            source.data.attributes['st_width'].data[-1].value=changed
            source.data.update();source.update_tag();current=reader.read(obj)
            if changed==200.:self.assertIn(65,{d['code'] for d in current['diagnostics']})
            else:self.assertFalse(current['diagnostics'],current)
        self.assertEqual(reader.read(obj),baseline)

    def test_oblique_bundle_uses_live_curved_leads_and_rejects_unqualified_exact_bends(self):
        from dataclasses import replace
        from blosm.way.street_network.evidence import GeometryEvidence
        snapshot=divided_sample()
        streets=list(snapshot.junctions[0].streets)
        for sid,slope in ((0,-.035),(1,.052)):
            row=snapshot.streets[sid]['items'][0]
            points=tuple((p[0]+slope*p[1],p[1],p[2]) for p in row['points'])
            row.update(points=points,clip_end=5.)
            streets[sid]=replace(streets[sid],sections=(replace(streets[sid].sections[0],centerline=points),))
        snapshot.junctions[0]=replace(snapshot.junctions[0],streets=tuple(streets))
        collection=bpy.data.collections.new('Oblique Bundle regression')
        bpy.context.scene.collection.children.link(collection)
        obj=EditableNetworkBuilder(StreetDesign(),self.blocks).create('Oblique Bundle',snapshot,{0,1,2,3},{0},collection)
        reader=GeometryEvidence()
        before=reader.read(obj)
        self.assertEqual(before['diagnostics'],[])
        self.assertEqual(before['section_ids'],[0,1,2,3])
        self.assertEqual(before['intersection_ids'],[0])
        source=next(o for o in obj.st_network.source_collection.objects if o.name.endswith('Street 0'))
        original=[v.co.copy() for v in source.data.vertices]
        for vertex in source.data.vertices:
            vertex.co.x+=.005*(source.matrix_world@vertex.co).y
        source.data.update();source.update_tag()
        self.assertEqual(reader.read(obj)['diagnostics'],[])
        for vertex,co in zip(source.data.vertices,original):vertex.co=co
        source.data.vertices[1].co.x+=.1
        source.data.update();source.update_tag()
        curved=reader.read(obj)
        self.assertEqual(curved['diagnostics'],[])
        self.assertEqual(curved['section_ids'],[0,1,2,3])
        self.assertEqual(curved['intersection_ids'],[0])
        for value in source.data.attributes['st_interpolation'].data:value.value=0
        source.data.update();source.update_tag()
        self.assertIn(62,{d['code'] for d in reader.read(obj)['diagnostics']})
        for value in source.data.attributes['st_interpolation'].data:value.value=3
        for vertex,co in zip(source.data.vertices,original):vertex.co=co
        source.data.update();source.update_tag()
        self.assertEqual(reader.read(obj),before)

    def test_import_preserves_metres_with_nondefault_scene_units(self):
        scale=bpy.context.scene.unit_settings.scale_length
        collection=bpy.data.collections.new('Metre conversion');bpy.context.scene.collection.children.link(collection)
        try:
            bpy.context.scene.unit_settings.scale_length=.1
            obj=EditableNetworkBuilder(StreetDesign(),self.blocks).create('Decimetre scene',sample(),{0,1,2},{0},collection)
            self.assertEqual(obj.st_network.status,'READY',obj.st_network.diagnostic)
            source=next(o for o in obj.st_network.source_collection.objects if o.name.endswith('Street 0'))
            self.assertAlmostEqual((source.matrix_world@source.data.vertices[-1].co).y*.1,30.,places=5)
            self.assertEqual(source.data.attributes['st_width'].data[0].value,6.)
            ev=obj.evaluated_get(bpy.context.evaluated_depsgraph_get());mesh=ev.to_mesh()
            try:
                self.assertTrue(all(len(f.vertices)==4 for f in mesh.polygons))
                self.assertAlmostEqual(max(v.co.z for v in mesh.vertices)*.1,.8,places=5)
            finally:ev.to_mesh_clear()
        finally:bpy.context.scene.unit_settings.scale_length=scale

    def test_transaction_preserves_existing_data(self):
        old = bpy.data.objects.new('Keep existing', None)
        before = set(bpy.data.objects)
        with self.assertRaises(ValueError):
            with ImportTransaction():
                mesh=bpy.data.meshes.new('Rollback mesh')
                obj=bpy.data.objects.new('Rollback object',mesh)
                bpy.context.scene.collection.objects.link(obj)
                raise ValueError('Deliberate import failure')
        self.assertEqual(before, set(bpy.data.objects))
        self.assertIsNone(bpy.data.meshes.get('Rollback mesh'))

    def test_builtin_styles_without_asset_directory(self):
        from blosm.way.street_network.defaults import styles
        from blosm.style import StyleStore
        from blosm.way.item.street import Street
        from blosm.way.dimension_policy import StreetDimensionPolicy
        store=StyleStore(None,'',styles())
        from mathutils import Vector
        s=Street(Vector((0.,0.)),Vector((1.,0.))); s.setStyleBlockFromTop(store.get('primary'))
        self.assertEqual(StreetDimensionPolicy().resolve(s,'primary',False)['nrLanes'],4)

    def test_transition_is_not_discarded(self):
        row = sample().streets[0]
        row['items'].insert(1,dict(kind='SideLane'))
        self.assertIn('SideLane',PreparedNetworkAdapter().street_support(row['items']))

    def test_native_shared_output_and_edit(self):
        collection=bpy.data.collections.new('Integration regression')
        bpy.context.scene.collection.children.link(collection)
        builder=EditableNetworkBuilder(StreetDesign(),self.blocks)
        outputs=[builder.create('Regression '+str(i),sample(),{0,1,2},{0},collection) for i in range(2)]
        self.assertIs(outputs[0].modifiers[0].node_group,outputs[1].modifiers[0].node_group)
        self.assertIs(outputs[0].modifiers[0].node_group,self.blocks['solver'])
        def capture(obj):
            bpy.context.view_layer.update()
            ev=obj.evaluated_get(bpy.context.evaluated_depsgraph_get()); mesh=ev.to_mesh()
            try:
                self.assertGreater(len(mesh.polygons),100)
                self.assertTrue(all(len(f.vertices)==4 for f in mesh.polygons))
                z=[v.co.z for v in mesh.vertices]
                self.assertGreaterEqual(min(z),-3e-5)
                self.assertLessEqual(max(z),.80003)
                return tuple(tuple(v.co) for v in mesh.vertices)
            finally: ev.to_mesh_clear()
        original=capture(outputs[0])
        source=next(o for o in outputs[0].st_network.source_collection.objects if o.name.endswith('Street 0'))
        self.assertEqual(set(range(len(source.data.vertices))), {i for e in source.data.edges for i in e.vertices})
        source.data.vertices[-1].co.x += 1.
        source.data.update(); source.update_tag()
        self.assertNotEqual(original,capture(outputs[0]))
        # No runtime callback has been registered and no refresh was called.
        self.assertFalse(runtime._timer_active)
        self.assertNotIn(runtime._changed,bpy.app.handlers.depsgraph_update_post)

    def test_native_edits_are_not_watched_as_legacy_after_registration(self):
        collection=bpy.data.collections.new('Registered native editing');bpy.context.scene.collection.children.link(collection)
        obj=EditableNetworkBuilder(StreetDesign(),self.blocks).create('Registered network',sample(),{0,1,2},{0},collection)
        source=next(iter(obj.st_network.source_collection.objects))
        runtime.register()
        try:
            self.assertFalse(any(obj.name in names for names in runtime._legacy_watch.values()))
            widths=source.data.attributes['st_width']
            for row in widths.data:row.value=100.
            source.data.update();source.update_tag()
            runtime.refresh(obj,force=True,raise_errors=True)
            self.assertEqual(obj.st_network.status,'INVALID')
            self.assertFalse(runtime._needs_binding_recovery(obj),'GN rejection must retain native recovery')
            runtime._changed(None,SimpleNamespace(updates=[SimpleNamespace(id=source)]))
            self.assertNotIn(obj.name,runtime._pending)
            for row in widths.data:row.value=6.
            source.data.update();source.update_tag()
            # A real binding failure still needs a one-shot Python recompile.
            roles=source.data.attributes['st_endpoint_role'];saved=[v.value for v in roles.data]
            for row in roles.data:row.value=0
            source.data.update();runtime.refresh(obj,force=True)
            self.assertTrue(runtime._needs_binding_recovery(obj))
            for row,value in zip(roles.data,saved):row.value=value
            source.data.update();source.update_tag()
            runtime._changed(None,SimpleNamespace(updates=[SimpleNamespace(id=source)]))
            self.assertIn(obj.name,runtime._pending)
            runtime._tick()
            self.assertEqual(obj.st_network.status,'READY',obj.st_network.diagnostic)
            self.assertFalse(runtime._needs_binding_recovery(obj))
        finally:runtime.unregister()

    def test_smooth_sections_share_a_profile_on_a_bent_seam(self):
        snapshot=NetworkSnapshot()
        sections=[]
        for i,points in enumerate((((-30.,-3.,0.),(-15.,0.,0.),(0.,0.,0.)),
                                  ((0.,0.,0.),(15.,3.,0.),(30.,4.,0.)))):
            sections.append(dict(kind='Section',id=i,points=points,width=6.,offset=0.,
                category='residential',tags={},source_way_ids=['123']))
        snapshot.streets[0]=dict(id=0,items=sections,status='SUPPORTED',reason=None,endpoints=(None,None),bundle_id=None)
        collection=bpy.data.collections.new('Bent seam regression'); bpy.context.scene.collection.children.link(collection)
        obj=EditableNetworkBuilder(StreetDesign(),self.blocks).create('Bent seam',snapshot,{0},set(),collection)
        self.assertEqual(obj.st_network.status,'READY',obj.st_network.diagnostic)
        ev=obj.evaluated_get(bpy.context.evaluated_depsgraph_get()); mesh=ev.to_mesh()
        try:
            from mathutils import Vector
            tangent=(Vector((1.,0.))+Vector((15.,3.)).normalized()).normalized()
            # Interior seam edges must be shared, not two unmatched open rows.
            users={tuple(sorted(e.vertices)):0 for e in mesh.edges}
            for face in mesh.polygons:
                corners=list(face.vertices)
                for a,b in zip(corners,corners[1:]+corners[:1]): users[tuple(sorted((a,b)))]+=1
            seam_edges=[e for e in mesh.edges if all(abs(Vector((obj.matrix_world@mesh.vertices[v].co)[:2]).dot(tangent))<2e-5 for v in e.vertices)]
            self.assertGreater(len(seam_edges),2)
            self.assertTrue(all(users[tuple(sorted(e.vertices))]==2 for e in seam_edges))
        finally: ev.to_mesh_clear()

    def test_native_endpoint_parameters_and_invalid_edit_recovery(self):
        from blosm.way.street_network.evidence import GeometryEvidence
        collection=bpy.data.collections.new('Native parameter edits');bpy.context.scene.collection.children.link(collection)
        obj=EditableNetworkBuilder(StreetDesign(),self.blocks).create('Native parameters',sample(),{0,1,2},{0},collection)
        source=next(o for o in obj.st_network.source_collection.objects if o.name.endswith('Street 0'))
        original=[v.co.copy() for v in source.data.vertices]
        def evidence():
            source.data.update();source.update_tag()
            return GeometryEvidence().read(obj)
        def valid():
            result=evidence()
            self.assertFalse(result['diagnostics'],result)
            self.assertEqual(result['section_ids'],[0,1,2])
            self.assertEqual(result['intersection_ids'],[0])
            self.assertEqual(result['other_faces'],0)
        valid()
        # The first control is at a junction; the next controls its approach.
        source.data.vertices[0].co.x+=.05
        valid()
        source.data.vertices[0].co=original[0]
        source.data.vertices[1].co.x+=.4
        for interpolation in (0,3):
            for row in source.data.attributes['st_interpolation'].data:row.value=interpolation
            valid()
        for row in source.data.attributes['st_width'].data:row.value=6.5
        for row in source.data.attributes['st_offset'].data:row.value=.2
        valid()
        source.data.vertices[1].co=source.data.vertices[0].co
        self.assertTrue(evidence()['diagnostics'],'A collapsed input edge must be diagnosed')
        for vertex,coordinate in zip(source.data.vertices,original):vertex.co=coordinate
        for row in source.data.attributes['st_width'].data:row.value=6.
        for row in source.data.attributes['st_offset'].data:row.value=0.
        valid()
        self.assertFalse(runtime._timer_active)

    def test_free_bundle_is_not_split_into_independent_crowned_streets(self):
        snapshot=NetworkSnapshot()
        for sid,x in enumerate((-4.,4.)):
            section=dict(kind='Section',id=sid,points=((x,0.,0.),(x,20.,0.),(x,40.,0.)),
                width=6.,offset=0.,category='residential',tags={},source_way_ids=[str(sid)])
            snapshot.streets[sid]=dict(id=sid,items=[section],status='SUPPORTED',reason=None,
                endpoints=(None,None),bundle_id=7)
        snapshot.bundles[7]=dict(id=7,head=(0,1),tail=(0,1),status='SUPPORTED',reason=None,
            members=[dict(street_id=i,forward=True,head_order=i,tail_order=i) for i in range(2)])
        collection=bpy.data.collections.new('Free Bundle regression');bpy.context.scene.collection.children.link(collection)
        from blosm.way.street_network.renderer import CarriagewayFeature, RenderContext
        context=RenderContext(None,snapshot,collection,StreetDesign(),self.loader)
        CarriagewayFeature().render(context)
        self.assertEqual(len(context.outputs),1)
        obj=context.outputs[0]
        self.assertEqual(obj.st_network.status,'READY',obj.st_network.diagnostic)
        self.assertEqual(len(obj.st_network.bundles),1)
        ev=obj.evaluated_get(bpy.context.evaluated_depsgraph_get());mesh=ev.to_mesh()
        try:
            self.assertGreater(len(mesh.polygons),100)
            self.assertTrue(all(len(f.vertices)==4 for f in mesh.polygons))
            self.assertEqual(len({f.material_index for f in mesh.polygons}),2,'Separator grass was not produced')
        finally:ev.to_mesh_clear()

    def test_unequal_bundle_mouth_widths_update_natively(self):
        from blosm.way.street_network.evidence import GeometryEvidence
        collection=bpy.data.collections.new('Unequal Bundle junction');bpy.context.scene.collection.children.link(collection)
        obj=EditableNetworkBuilder(StreetDesign(),self.blocks).create('Unequal joined Bundle',
            divided_sample(reverse_second=True),{0,1,2,3},{0},collection)
        self.assertEqual(obj.st_network.status,'READY',obj.st_network.diagnostic)
        sources=[next(o for o in obj.st_network.source_collection.objects if o.name.endswith('Street '+str(sid))) for sid in (0,1)]
        for widths in ((4.,8.),(8.,4.),(5.4,8.1)):
            for source,width in zip(sources,widths):
                for row in source.data.attributes['st_width'].data:row.value=width
                source.data.update();source.update_tag()
            evidence=GeometryEvidence().read(obj)
            self.assertFalse(evidence['diagnostics'],evidence)
            self.assertEqual(evidence['section_ids'],[0,1,2,3])
            self.assertEqual(evidence['intersection_ids'],[0])
            self.assertEqual(evidence['other_faces'],0)
        self.assertFalse(runtime._timer_active)

    def test_explicit_selection_clip_preserves_authored_points(self):
        snapshot=NetworkSnapshot()
        section=dict(kind='Section',id=0,points=((0.,0.,0.),(0.,20.,0.),(0.,40.,0.)),
            width=6.,offset=0.,category='residential',tags={},source_way_ids=['1'],clip_end=5.)
        snapshot.streets[0]=dict(id=0,items=[section],status='SUPPORTED',reason=None,
            endpoints=(None,None),bundle_id=None)
        collection=bpy.data.collections.new('Selection clip regression');bpy.context.scene.collection.children.link(collection)
        obj=EditableNetworkBuilder(StreetDesign(),self.blocks).create('Clipped extent',snapshot,{0},set(),collection)
        self.assertEqual(obj.st_network.status,'READY',obj.st_network.diagnostic)
        source=next(iter(obj.st_network.source_collection.objects))
        self.assertAlmostEqual((source.matrix_world@source.data.vertices[-1].co).y,40.)
        ev=obj.evaluated_get(bpy.context.evaluated_depsgraph_get());mesh=ev.to_mesh()
        try:self.assertAlmostEqual(max((obj.matrix_world@v.co).y for v in mesh.vertices),35.,places=4)
        finally:ev.to_mesh_clear()

    def test_open_bundle_boundary_clip_is_explicit_and_preserves_sources(self):
        from blosm.way.street_network.renderer import CarriagewayFeature, RenderContext
        for clearance in (0., 5.):
            snapshot=divided_sample()
            snapshot.bundles[7]=dict(id=7,head=(0,1),tail=(0,1),status='SUPPORTED',reason=None)
            original=tuple(snapshot.streets[s]['items'][0]['points'] for s in (0,1))
            collection=bpy.data.collections.new('Explicit open Bundle clipping '+str(clearance))
            bpy.context.scene.collection.children.link(collection)
            context=RenderContext(SimpleNamespace(),snapshot,collection,StreetDesign(),self.loader)
            CarriagewayFeature(junction_ids={0},open_bundle_clearance=clearance).render(context)
            self.assertEqual(context.results[0]['status'],'READY',context.results)
            for sid in (0,1):
                self.assertEqual(snapshot.streets[sid]['items'][0]['points'],original[sid])
                self.assertEqual(snapshot.streets[sid]['items'][0].get('clip_end',0.),clearance)
                source=next(o for o in context.outputs[0].st_network.source_collection.objects
                            if o.data.attributes['st_section_id'].data[0].value==sid)
                self.assertAlmostEqual((source.matrix_world@source.data.vertices[-1].co).y,35.)
            explicit=[i for i in snapshot.issues if i['code']=='EXPLICIT_OPEN_BUNDLE_CLIP']
            self.assertEqual(len(explicit),2 if clearance else 0)

    def test_coverage_does_not_count_partial_mesh_as_generated_length(self):
        from blosm.way.street_network.coverage import CoverageReport
        snapshot=sample()
        for i,status in enumerate(('GENERATED','PARTIAL','DELIBERATELY_EXCLUDED')):
            snapshot.streets[i]['surface_status']=status
        report=CoverageReport().build(snapshot)
        self.assertEqual(report['street_counts'],dict(GENERATED=1,PARTIAL=1,DELIBERATELY_EXCLUDED=1))
        self.assertEqual(report['centerline_length_m']['GENERATED'],30.)

    def test_connectivity_audit_distinguishes_geometry_from_osm_topology(self):
        from blosm.way.connectivity_audit import ConnectivityAudit
        def way(sid,points,nodes,tags):
            return SimpleNamespace(source_id=sid,category='residential',element=SimpleNamespace(tags=tags),
                segments=[SimpleNamespace(v1=points[0],v2=points[1],id1=nodes[0],id2=nodes[1])])
        a=way('1',((-1.,0.),(1.,0.)),('a','b'),{})
        b=way('2',((0.,-1.),(0.,1.)),('c','d'),{'layer':'1'})
        splits={(a.segments[0].v1,a.segments[0].v2):((-1.,0.),(0.,0.),(1.,0.)),
                (b.segments[0].v1,b.segments[0].v2):((0.,-1.),(0.,0.),(0.,1.))}
        self.assertEqual(ConnectivityAudit().read([a,b],splits)['counts'],{'GRADE_SEPARATED_XY':1})
        a=way('1',((-1.,0.),(0.,0.)),('a','common'),{})
        b=way('2',((0.,0.),(0.,1.)),('common','d'),{})
        self.assertEqual(ConnectivityAudit().read([a,b],{})['events'],[])

    def test_unsupported_junction_keeps_neighbors_with_an_explicit_gap(self):
        from blosm.way.street_network.renderer import CarriagewayFeature, RenderContext
        snapshot=sample()
        snapshot.streets[1]['status']='UNSUPPORTED'
        snapshot.streets[1]['reason']='Unsupported transition regression'
        snapshot.streets[1]['items'][0]['width']=20.
        snapshot.junctions.clear()
        snapshot.intersections[0].update(status='UNSUPPORTED',location=(0.,0.),reason='Unsupported transition regression')
        collection=bpy.data.collections.new('Local failure regression');bpy.context.scene.collection.children.link(collection)
        context=RenderContext(None,snapshot,collection,StreetDesign(),self.loader)
        CarriagewayFeature().render(context)
        self.assertEqual(len(context.outputs),2)
        for obj in context.outputs:
            self.assertEqual(obj.st_network.status,'READY',obj.st_network.diagnostic)
            ev=obj.evaluated_get(bpy.context.evaluated_depsgraph_get());mesh=ev.to_mesh()
            try:self.assertGreaterEqual(min((obj.matrix_world@v.co).xy.length for v in mesh.vertices),11.999)
            finally:ev.to_mesh_clear()
        self.assertEqual(snapshot.streets[0]['surface_status'],'PARTIAL')
        self.assertEqual(snapshot.streets[2]['surface_status'],'PARTIAL')

    def test_unverified_xy_connectivity_is_not_published_as_a_junction(self):
        from blosm.way.street_network.connectivity import SurfaceConnectivityReview
        snapshot=sample()
        snapshot.preparation_audit=dict(events=[dict(osm_way_ids=['90000000000','missing'],
            code='INFERRED_XY_CONNECTION',location=(0.,0.),categories=('residential','footway'))])
        SurfaceConnectivityReview().apply(snapshot,StreetDesign().categories)
        self.assertEqual(snapshot.streets[0]['status'],'UNSUPPORTED')
        self.assertEqual(snapshot.streets[1]['status'],'SUPPORTED')
        self.assertEqual(snapshot.preparation_audit['surface_review_counts'],dict(unverified=1,deferred=0,unmapped=0))
        self.assertEqual(snapshot.issues[0]['code'],'UNVERIFIED_CONNECTIVITY')
        snapshot=sample()
        snapshot.preparation_audit=dict(events=[dict(osm_way_ids=['foot','cycle'],
            code='INFERRED_XY_CONNECTION',location=(10.,10.),categories=('footway','cycleway'))])
        SurfaceConnectivityReview().apply(snapshot,StreetDesign().categories)
        self.assertTrue(all(s['status']=='SUPPORTED' for s in snapshot.streets.values()))
        self.assertEqual(snapshot.preparation_audit['surface_review_counts'],dict(unverified=0,deferred=1,unmapped=0))
        snapshot.preparation_audit['events'][0]['categories']=('residential','cycleway')
        SurfaceConnectivityReview().apply(snapshot,StreetDesign().categories)
        self.assertEqual(snapshot.preparation_audit['surface_review_counts'],dict(unverified=0,deferred=0,unmapped=1))
        self.assertEqual(snapshot.issues[-1]['code'],'UNMAPPED_SOURCE')

    def test_provenance_survives_network_copy_join_and_inversion(self):
        from blosm.way.way_network import NetSection
        a=NetSection((0.,0.),(1.,0.),'residential',{})
        b=NetSection((1.,0.),(2.,0.),'residential',{})
        a.sourceWayIds=frozenset(('123',)); b.sourceWayIds=frozenset(('456',))
        a.sourceSegments=(('123','1','2',.2,.6),)
        b.sourceSegments=(('456','3','4',0.,1.),)
        copy=NetSection(a);copy.join(b);reverse=~copy
        self.assertEqual(reverse.sourceWayIds,frozenset(('123','456')))
        self.assertEqual(reverse.sourceSegments,(('456','3','4',1.,0.),('123','1','2',.6,.2)))

    def test_transition_adapter_preserves_prepared_section_dimensions(self):
        from blosm.way.street_network.transitions import PreparedTransitionAdapter
        for kind in ('SymLane','SideLane'):
            a=dict(kind='Section',id=10,points=((-40.,0.,0.),(0.,0.,0.)),width=6.,offset=0.)
            b=dict(kind='Section',id=11,points=((0.,0.,0.),(40.,0.,0.)),width=10.,offset=2.)
            transition=dict(kind=kind,id=4,points=())
            PreparedTransitionAdapter().read([a,transition,b])
            self.assertEqual(transition['surface_treatment'],'WIDTH_TRANSITION')
            self.assertEqual((a['width'],a['offset'],b['width'],b['offset']),(6.,0.,10.,2.))
            b['points']=((1.,0.,0.),(40.,0.,0.))
            rejected=dict(kind=kind,id=4,points=())
            PreparedTransitionAdapter().read([a,rejected,b])
            self.assertNotIn('surface_treatment',rejected)
            self.assertIn('disconnected',rejected['transition_reason'])

    def test_native_transition_and_coverage_require_generated_transition_faces(self):
        from blosm.way.street_network.adapter import NetworkSnapshot
        from blosm.way.street_network.transitions import PreparedTransitionAdapter
        from blosm.way.street_network.renderer import CarriagewayFeature,RenderContext
        from blosm.way.street_network.evidence import GeometryEvidence
        base=sample().streets[0]['items'][0]
        a=dict(base,id=10,points=((-40.,0.,0.),(-20.,0.,0.),(0.,0.,0.)),width=6.,offset=0.)
        b=dict(base,id=11,points=((0.,0.,0.),(20.,0.,0.),(40.,0.,0.)),width=10.,offset=2.)
        transition=dict(kind='SideLane',id=3,points=())
        PreparedTransitionAdapter().read([a,transition,b])
        snapshot=NetworkSnapshot(streets={5:dict(id=5,items=[a,transition,b],status='SUPPORTED',reason=None,
                                                endpoints=(None,None),bundle_id=None)})
        collection=bpy.data.collections.new('Native transition regression');bpy.context.scene.collection.children.link(collection)
        context=RenderContext(None,snapshot,collection,StreetDesign(),self.loader)
        CarriagewayFeature().render(context)
        self.assertEqual(snapshot.streets[5]['surface_status'],'GENERATED',context.results)
        obj=context.outputs[0]
        self.assertEqual(GeometryEvidence().read(obj)['transition_ids'],[0])
        source=next(iter(obj.st_network.source_collection.objects));widths=source.data.attributes['st_width'].data
        widths[2].value=widths[3].value=100.;source.data.update();source.update_tag()
        evidence=GeometryEvidence().read(obj)
        self.assertTrue(any(d['code']==65 for d in evidence['diagnostics']),evidence)
        widths[2].value=widths[3].value=10.;source.data.update();source.update_tag()
        self.assertFalse(GeometryEvidence().read(obj)['diagnostics'])

    def test_deferred_attachment_taper_uses_shared_asset_and_keeps_crossing(self):
        from copy import deepcopy
        from blosm.way.street_network.transitions import PreparedTransitionAdapter
        from blosm.way.street_network.renderer import CarriagewayFeature,RenderContext
        from blosm.way.street_network.evidence import GeometryEvidence
        base=sample().streets[0]['items'][0]
        a=dict(base,id=10,points=((-40.,0.,0.),(-20.,0.,0.),(0.,0.,0.)),width=6.,offset=0.)
        b=dict(base,id=11,points=((0.,0.,0.),(20.,0.,0.),(40.,0.,0.)),width=9.,offset=-.4)
        crossing=dict(kind='Intersection',id=33,points=(),location=(0.,0.),through=((10,1),(11,0)),
            attachments=[dict(street_id=99,section_id=99,connector_id=99,end=0,side='LEFT',category='footway')],
            surface_treatment='DEFERRED_ATTACHMENT')
        original=deepcopy(crossing)
        PreparedTransitionAdapter().read([a,crossing,b])
        self.assertIsNone(PreparedNetworkAdapter().street_support([a,crossing,b]))
        snapshot=NetworkSnapshot(streets={5:dict(id=5,items=[a,crossing,b],status='SUPPORTED',reason=None,
                                                endpoints=(None,None),bundle_id=None)})
        collection=bpy.data.collections.new('Deferred attachment taper');bpy.context.scene.collection.children.link(collection)
        context=RenderContext(None,snapshot,collection,StreetDesign(),self.loader)
        CarriagewayFeature().render(context)
        self.assertEqual(snapshot.streets[5]['surface_status'],'GENERATED',context.results)
        self.assertEqual({k:crossing[k] for k in original},original)
        obj=context.outputs[0];record=obj.st_network.transitions[0]
        self.assertEqual((record.source_kind,record.source_id,record.before_section,record.after_section),
                         ('Intersection',33,10,11))
        self.assertEqual(GeometryEvidence().read(obj)['transition_ids'],[0])
        source=next(iter(obj.st_network.source_collection.objects));mesh=source.data
        self.assertEqual((len(mesh.vertices),len(mesh.edges)),(6,4))
        widths=mesh.attributes['st_width'].data
        widths[2].value=widths[3].value=10.;mesh.update();source.update_tag()
        self.assertFalse(GeometryEvidence().read(obj)['diagnostics'])
        # A declaration on a malformed deferred attachment must not remain
        # eligible merely because its surface role still says deferred.
        crossing['through']=((10,0),(11,0));PreparedTransitionAdapter().read([a,crossing,b])
        self.assertIn('owners',PreparedNetworkAdapter().street_support([a,crossing,b]))

    def test_two_transition_extents_cannot_overlap(self):
        from blosm.way.street_network.adapter import NetworkSnapshot
        from blosm.way.street_network.transitions import PreparedTransitionAdapter
        from blosm.way.street_network.renderer import CarriagewayFeature,RenderContext
        from blosm.way.street_network.evidence import GeometryEvidence
        base=sample().streets[0]['items'][0]
        sections=[dict(base,id=10+k,points=((a,0.,0.),(b,0.,0.)),width=w,offset=0.)
                  for k,(a,b,w) in enumerate(((-50.,-10.,6.),(-10.,10.,10.),(10.,50.,6.)))]
        items=[sections[0],dict(kind='SymLane',id=0,points=()),sections[1],
               dict(kind='SymLane',id=1,points=()),sections[2]]
        PreparedTransitionAdapter().read(items)
        snapshot=NetworkSnapshot(streets={5:dict(id=5,items=items,status='SUPPORTED',reason=None,
                                                endpoints=(None,None),bundle_id=None)})
        collection=bpy.data.collections.new('Two taper extents');bpy.context.scene.collection.children.link(collection)
        context=RenderContext(None,snapshot,collection,StreetDesign(),self.loader)
        CarriagewayFeature().render(context);obj=context.outputs[0]
        self.assertEqual(GeometryEvidence().read(obj)['transition_ids'],[0,1])
        source=next(iter(obj.st_network.source_collection.objects))
        source.data.attributes['st_width'].data[1].value=40.;source.data.update();source.update_tag()
        evidence=GeometryEvidence().read(obj)
        self.assertEqual({d['section_id'] for d in evidence['diagnostics'] if d['code']==65},{10,11,12})
        source.data.attributes['st_width'].data[1].value=10.;source.data.update();source.update_tag()
        self.assertFalse(GeometryEvidence().read(obj)['diagnostics'])

    def test_bundle_outlier_with_two_coincident_opposite_ends(self):
        from mathutils import Vector
        from blosm.way.item.bundle import orderHeadTail
        class Street:
            pred=succ=None
            def __init__(self,a,b):
                self.src=Vector(a).freeze();self.dst=Vector(b).freeze()
                self.head=self.tail=SimpleNamespace(centerline=[self.src,self.dst])
            def endVectors(self):return self.dst-self.src,self.src-self.dst
        class Group:
            def __init__(self,streets):self.group=streets
            def __iter__(self):return iter(self.group)
        for head_outlier in (True,False):
            outlier=Street((-1000.,0.) if head_outlier else (0.,0.),
                           (20.,0.) if head_outlier else (1000.,0.))
            group=Group([outlier,*[Street((0.,float(y)),(20.,float(y))) for y in (0,1,2)]])
            head,tail=orderHeadTail(group)
            self.assertNotIn(outlier,group.group)
            self.assertTrue(head and tail)
            self.assertNotIn(outlier,[end['street'] for end in head+tail])

    def test_bundle_centerline_preserves_loops_and_section_seams(self):
        from mathutils import Vector
        from blosm.way.item import Section, Street
        from blosm.way.item.bundle import bundleCenterline, parallelToSample
        def street(parts):
            obj=Street(Vector(parts[0][0]),Vector(parts[-1][-1]))
            for controls in parts:
                section=Section.__new__(Section)
                section.centerline=[Vector(p) for p in controls]
                section.category='residential'
                obj.insertEnd(section)
            return obj
        a,b,c=(0.,0.),(40.,0.),(80.,0.)
        open_street=street(((a,b),(b,c)))
        self.assertEqual([tuple(v) for v in open_street.iterCenterlineVertices()],[a,b,c])
        self.assertEqual(len(bundleCenterline(open_street)),3)
        for parts in (((a,b,c,a),),((a,b),(b,c,a)),((a,b,c,b,(90.,0.)),)):
            loop=street(parts)
            original=[[tuple(p) for p in s.centerline] for s in loop.iterItems()]
            self.assertIsNone(bundleCenterline(loop))
            self.assertFalse(parallelToSample(open_street,loop))
            self.assertFalse(parallelToSample(loop,open_street))
            self.assertEqual([[tuple(p) for p in s.centerline] for s in loop.iterItems()],original)
        collapsed=street(((a,a),))
        self.assertIsNone(bundleCenterline(collapsed))

    def test_parallel_preparation_keeps_a_nearby_loop_outside_the_bundle(self):
        from math import cos, sin, tau
        from mathutils import Vector
        from blosm.action.generate_streets import StreetGenerator
        from blosm.way.item import Section, Street
        def street(points):
            controls=[Vector(p).freeze() for p in points]
            obj=Street(controls[0],controls[-1])
            section=Section.__new__(Section)
            section.centerline=controls;section.category='residential'
            obj.insertEnd(section)
            return obj
        a=street(((0.,0.),(40.,0.),(80.,0.)))
        b=street(((0.,8.),(40.,8.),(80.,8.)))
        controls=[(40.+8.*cos(tau*i/32),8.+8.*sin(tau*i/32)) for i in range(32)]
        loop=street(controls+[controls[0]])
        population=(a,b,loop)
        generator=StreetGenerator(None,None)
        generator.wayManager=SimpleNamespace(iterStreets=lambda:iter(population))
        generator.createParallelStreets()
        groups=[set(group) for group in generator.parallelStreets]
        self.assertEqual(groups,[{a,b}])
        self.assertEqual(tuple(generator.wayManager.iterStreets()),population)
        self.assertEqual(loop.head.centerline[0],loop.head.centerline[-1])


if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(StreetIntegrationTests))
    if not result.wasSuccessful():
        raise RuntimeError('Street integration regression failed')
