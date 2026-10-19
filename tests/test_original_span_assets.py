"""Validate the distributed original-span library and required feature inputs.

Run in factory-startup background Blender after building the current asset.
Malformed copies test preflight rejection before source objects are authored.
"""
from pathlib import Path
import argparse
import json
import re
from dataclasses import replace
from functools import partial
import sys
import unittest

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import blosm
from blosm.app.asset_library import BlenderAssetLoader
from blosm.way.street_network import assets


class OriginalSpanAssetTests(unittest.TestCase):
    linked_library = False
    require_local_failures=False
    require_standalone=False
    require_filtered_ports=False
    require_section_scopes=False
    require_internal_ends=False
    require_shared_owners=False
    require_seams=False
    require_category_connections=False
    require_stable_quad_energy=False
    require_free_ended=False
    require_boundary_runs=False
    require_chained=False
    require_free_chained=False
    require_crossing_regions=False
    require_continuation_regions=False
    require_internal_regions=False
    require_assembly_payload_cleanup=False
    require_interpolation=False
    require_cycle_cores=False
    require_shared_kernels=False
    require_strip_refinement=False
    require_guide_correction=False
    require_guarded_profile=False
    require_core_homotopy=False
    require_separator_parity=False
    require_checked_chain_contacts=False
    require_scoped_separator_parity=False
    require_selection_boundaries=False
    require_authored_through_connectors=False
    require_authored_dimension_runs=False
    require_authored_dimension_approaches=False
    require_interior_selection_cuts=False
    require_rejected_profiles=False
    require_longitudinal_separator_chains=False

    def check_original_port_reader(self, names, suffix):
        """Keep exact constructor names while supporting the guarded revision."""
        prefix='Prochitecture Streets / Original Section-end port path '
        old=prefix+'v3 / Feature curves and seams'+suffix
        new=prefix+'v4 / Profile admission'+suffix
        self.assertTrue(set(names).intersection((old,new)),(old,new))
        if new in names:
            self.assertEqual(bpy.data.node_groups[new].get('st_port_profile_admission'),1)

    def test_rejected_profiles_require_consistent_native_consumers(self):
        roots=tuple(self.blocks[key] for key in
                    (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER))
        enabled=assets.validate_original_rejected_profiles(*roots)
        if self.require_rejected_profiles:self.assertTrue(enabled)
        if not enabled:self.skipTest('No explicit unresolved-profile support in this revision')
        readers={tree.as_pointer():tree for root in roots for tree in self.reachable_groups(root).values()
                 if tree.name.startswith('Prochitecture Streets / Original Section-end port path ')
                 and ' / Shared feature producer v1' not in tree.name}
        self.assertTrue(readers)
        for tree in readers.values():self.check_profile_reader_admission(tree)

    def check_profile_reader_admission(self,tree):
        self.assertEqual(tree.get('st_port_profile_admission'),1,tree.name)
        sockets={s.name:s.socket_type for s in tree.interface.items_tree
                 if s.item_type=='SOCKET' and s.in_out=='OUTPUT'}
        output,=[n for n in tree.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output]
        for key,kind in (('Geometry','NodeSocketGeometry'),('Valid','NodeSocketBool')):
            self.assertEqual(sockets.get(key),kind,tree.name)
            self.assertTrue(output.inputs[key].is_linked,(tree.name,key))
            pending=[link.from_node for link in output.inputs[key].links];live=set()
            while pending:
                node=pending.pop()
                if node in live:continue
                live.add(node)
                pending.extend(link.from_node for socket in node.inputs for link in socket.links)
            for field,datatype in (('st_cj_profile_error','BOOLEAN'),('st_cj_half','INT')):
                self.assertTrue(any(n.bl_idname=='GeometryNodeInputNamedAttribute' and n.data_type==datatype
                                    and n.inputs['Name'].default_value==field for n in live),
                                (tree.name,key,field))

    def test_interior_selection_cuts_have_reciprocal_native_readers(self):
        roots=[self.blocks[ref] for ref in (assets.ORIGINAL_SOLVER,
            assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER)]
        flags=[root.get('blosm_interior_selection_cuts',0) for root in roots]
        self.assertEqual(len(set(flags)),1);self.assertIn(flags[0],(0,1))
        if self.require_interior_selection_cuts:self.assertEqual(flags[0],1)
        if not flags[0]:self.skipTest('No paired interior selection cuts in this revision')
        groups={name:tree for root in roots for name,tree in self.reachable_groups(root).items()}
        consumers=[tree for tree in groups.values() if tree.get('st_interior_selection_cuts')==1]
        self.assertEqual(len(consumers),4)
        readers=[tree for name,tree in groups.items() if name.startswith(
            'Prochitecture Streets / Category-bounded original span path v1')]
        self.assertEqual(len(readers),2)
        for tree in consumers+readers:
            names={n.inputs['Name'].default_value for n in tree.nodes if n.bl_idname=='GeometryNodeInputNamedAttribute'}
            self.assertTrue({'st_ps_start_clip_peer','st_ps_end_clip_peer'}<=names)
        self.assertTrue(all('Paired selection spans' in tree.name for tree in readers))

    def test_missing_interior_cut_capability_rejects_before_authoring(self):
        from types import SimpleNamespace as NS
        from blosm.way.street_network.original_span_authoring import OriginalSpanNetworkBuilder
        from blosm.way.street_network.design import StreetDesign
        root=self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            if 'blosm_interior_selection_cuts' in root:del root['blosm_interior_selection_cuts']
            root['blosm_selection_boundaries']=1
            builder=OriginalSpanNetworkBuilder(StreetDesign(),root,self.blocks[assets.ORIGINAL_PAVEMENT],
                patch_solver=self.blocks[assets.ORIGINAL_PATCH_PROVIDER],strip_solver=self.blocks[assets.ORIGINAL_STRIP_PROVIDER])
            plan=NS(sites=[],selection_boundaries=[NS(endpoint=1)],selection_peer=lambda endpoint:2)
            before={o.as_pointer() for o in bpy.data.objects}
            with self.assertRaisesRegex(ValueError,'paired interior selection cuts'):
                builder.create('Should not create interior-cut objects',None,NS(ports=plan,network=None,free_bundles=()),bpy.context.scene.collection)
            self.assertEqual(before,{o.as_pointer() for o in bpy.data.objects})
        finally:bpy.data.node_groups.remove(root)

    def test_authored_dimension_approaches_use_shared_declared_dispatch(self):
        roots=[self.blocks[ref] for ref in (assets.ORIGINAL_SOLVER,
            assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER)]
        flags=[root.get('blosm_authored_dimension_approaches',0) for root in roots]
        self.assertEqual(len(set(flags)),1);self.assertIn(flags[0],(0,1))
        if self.require_authored_dimension_approaches:self.assertEqual(flags[0],1)
        if not flags[0]:self.skipTest('Authored dimensional approaches are not in this revision')
        groups={name:tree for root in roots for name,tree in self.reachable_groups(root).items()}
        tree,=[g for g in groups.values() if g.get('st_authored_dimension_approaches')]
        target,=[n for n in tree.nodes if n.type=='GROUP' and n.label=='Declared original dimensional approach']
        self.assertTrue(target.inputs['Transitions'].is_linked)
        self.assertIn('Declared dimensional transitions v1',target.node_tree.name)
        for caller in groups.values():
            for node in caller.nodes:
                if node.type=='GROUP' and node.node_tree==tree:self.assertTrue(node.inputs['Transitions'].is_linked)
        source=[g for g in groups.values() if 'Authored transition identity v1' in g.name]
        self.assertTrue(source)

    def test_authored_dimension_runs_use_shared_declared_dispatch(self):
        roots=[self.blocks[ref] for ref in (assets.ORIGINAL_SOLVER,
            assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER)]
        flags=[root.get('blosm_authored_dimension_runs',0) for root in roots]
        self.assertEqual(len(set(flags)),1)
        self.assertIn(flags[0],(0,1))
        if self.require_authored_dimension_runs:self.assertEqual(flags[0],1)
        if not flags[0]:self.skipTest('Authored dimensional runs are not in this revision')
        groups=self.reachable_groups(roots[0])
        dispatch=[g for g in groups.values() if g.get('st_authored_dimension_runs')]
        self.assertEqual(len(dispatch),2)
        for tree in dispatch:
            self.assertIn('Transitions',{s.name for s in tree.interface.items_tree if s.item_type=='SOCKET' and s.in_out=='INPUT'})
            target,=[n for n in tree.nodes if n.type=='GROUP' and n.label=='Declared original dimensional run']
            self.assertTrue(target.inputs['Transitions'].is_linked)
            self.assertIn('Declared dimensional transitions v1',target.node_tree.name)
            for caller in groups.values():
                for node in caller.nodes:
                    if node.type=='GROUP' and node.node_tree==tree:
                        self.assertTrue(node.inputs['Transitions'].is_linked)

    def test_authored_through_connector_capability_reaches_all_port_readers(self):
        roots=[self.blocks[ref] for ref in (assets.ORIGINAL_SOLVER,
            assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER)]
        roots.append(self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY])
        flags=[root.get('blosm_authored_through_connectors',0) for root in roots]
        self.assertEqual(len(set(flags)),1)
        self.assertIn(flags[0],(0,1))
        if self.require_authored_through_connectors:self.assertEqual(flags[0],1)
        if not flags[0]:self.skipTest('Authored through Connectors are not in this revision')
        groups={name:tree for root in roots for name,tree in self.reachable_groups(root).items()}
        readers=[g for g in groups.values() if g.name.startswith('Prochitecture Streets / Original Section-end port path ')
                 and 'Shared feature producer' not in g.name]
        self.assertTrue(readers)
        self.assertTrue(all(g.get('st_authored_through_connectors')==1 for g in readers))

    def test_missing_through_connector_capability_rejects_before_authoring(self):
        from types import SimpleNamespace as NS
        from blosm.way.street_network.original_span_authoring import OriginalSpanNetworkBuilder
        from blosm.way.street_network.design import StreetDesign
        root=self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            if 'blosm_authored_through_connectors' in root:del root['blosm_authored_through_connectors']
            builder=OriginalSpanNetworkBuilder(StreetDesign(),root,self.blocks[assets.ORIGINAL_PAVEMENT],
                patch_solver=self.blocks[assets.ORIGINAL_PATCH_PROVIDER],strip_solver=self.blocks[assets.ORIGINAL_STRIP_PROVIDER])
            port=NS(source_kind='THROUGH',connector_id=123)
            component=NS(ports=NS(sites=[NS(ports=[port])]),network=None,free_bundles=())
            before={o.as_pointer() for o in bpy.data.objects}
            with self.assertRaisesRegex(ValueError,'authored through Connectors'):
                builder.create('Should not be authored',None,component,bpy.context.scene.collection)
            self.assertEqual(before,{o.as_pointer() for o in bpy.data.objects})
        finally:bpy.data.node_groups.remove(root)

    def test_selection_boundaries_share_owned_live_path_consumers(self):
        roots=[self.blocks[ref] for ref in (assets.ORIGINAL_SOLVER,
            assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER)]
        flags=[root.get('blosm_selection_boundaries',0) for root in roots]
        self.assertEqual(len(set(flags)),1)
        self.assertIn(flags[0],(0,1))
        if self.require_selection_boundaries:self.assertEqual(flags[0],1)
        if not flags[0]:self.skipTest('Selection trims are not in this revision')
        groups={name:tree for root in roots for name,tree in self.reachable_groups(root).items()}
        for prefix in ('Prochitecture Streets / Original curved one-ended approach v1',
                       'Prochitecture Streets / Category-bounded profiled run v1'):
            consumers=[tree for tree in groups.values() if tree.name.startswith(prefix)
                       and not tree.get('st_authored_dimension_runs') and not tree.get('st_authored_dimension_approaches')]
            capability='blosm_authored_dimension_runs' if 'profiled run' in prefix else 'blosm_authored_dimension_approaches'
            expected=2 if roots[0].get(capability) else 1
            self.assertEqual(len(consumers),expected)
            for consumer in consumers:
                self.assertEqual(consumer.get('st_selection_boundaries'),1)
                names={n.inputs['Name'].default_value for n in consumer.nodes
                       if n.bl_idname=='GeometryNodeInputNamedAttribute'}
                self.assertTrue({'st_ps_start_clip','st_ps_end_clip'}<=names)

    def test_depth_recovery_reserve_preserves_inactive_geometry(self):
        roots=[self.blocks[ref] for ref in (assets.ORIGINAL_SOLVER,
            assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER)]
        roots.append(self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY])
        depths=[tree.get('blosm_corner_minimum_profile_depth',0.) for tree in roots]
        reserves=[tree.get('blosm_corner_profile_slope_reserve',0.) for tree in roots]
        self.assertEqual(len(set(depths)),1)
        self.assertEqual(len(set(reserves)),1)
        if not depths[0]:
            self.assertFalse(reserves[0])
            self.skipTest('No minimum-depth corner recovery in this revision')
        self.assertGreater(depths[0],.5);self.assertLess(depths[0],1.)
        self.assertGreaterEqual(reserves[0],0.);self.assertLessEqual(reserves[0],.1)
        # Patch providers evaluate independently; the assembly root receives
        # their results through object/collection inputs, not group calls.
        groups=list(self.reachable_groups(roots[1]).values())
        layouts=[tree for tree in groups if tree.get('st_capped_profile_depth')
                 and tree.get('st_corner_depth_floor')]
        self.assertTrue(layouts,'Depth capability has no reachable layout')
        for tree in layouts:
            self.assertAlmostEqual(tree['st_corner_depth_floor'],depths[0])
            if not reserves[0]:continue
            proposal,=[n for n in tree.nodes if n.type=='GROUP'
                       and n.label=='Recover profile room after infeasible guide compression']
            stored,=[n for n in tree.nodes if n.bl_idname=='GeometryNodeStoreNamedAttribute'
                     and n.inputs['Name'].default_value=='st_cq_profile_slope_reserve']
            self.assertEqual(stored.data_type,'FLOAT');self.assertEqual(stored.domain,'POINT')
            self.assertFalse(stored.inputs['Selection'].is_linked)
            self.assertTrue(stored.inputs['Selection'].default_value)
            self.assertFalse(stored.inputs['Value'].is_linked)
            self.assertAlmostEqual(stored.inputs['Value'].default_value,reserves[0])
            original=stored.inputs['Geometry'].links[0].from_socket
            switches=[n for n in tree.nodes if n.bl_idname=='GeometryNodeSwitch' and n.input_type=='GEOMETRY'
                      and any(link.from_node==stored for link in n.inputs['True'].links)]
            self.assertEqual(len(switches),1)
            self.assertEqual(switches[0].inputs['False'].links[0].from_socket,original)
            self.assertEqual(switches[0].inputs['Switch'].links[0].from_socket,proposal.outputs['Applied'])
        heights=[tree for tree in groups if tree.get('st_height_profile_experiment')=='physical_depth_capped_quadratic_v1']
        self.assertTrue(heights)
        if reserves[0]:
            self.assertTrue(all(tree.get('st_corner_profile_slope_reserve') for tree in heights))

    def test_standalone_separator_boundary_welds_are_complete(self):
        roots=[self.blocks[ref] for ref in (assets.ORIGINAL_SOLVER,
            assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER)]
        roots.append(self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY])
        flags=tuple(tree.get('blosm_separator_boundary_weld',0) for tree in roots)
        self.assertIn(flags,((0,)*4,(1,)*4))
        if not flags[0]:self.skipTest('No standalone boundary-weld upgrade in this revision')
        groups=self.reachable_groups(roots[0]).values()
        standalone=[tree for tree in groups
            if tree.name.startswith('Prochitecture Streets / Curved capped separator quads v1')
            and 'Shared standalone publication' in tree.name]
        self.assertTrue(standalone,'Boundary-weld capability has no standalone callers')
        for tree in standalone:
            with self.subTest(tree=tree.name):
                self.assertEqual(tree.get('st_boundary_weld_candidate'),1)
                merges=[n for n in tree.nodes if n.bl_idname=='GeometryNodeMergeByDistance'
                        and n.label=='Weld exact integer identities only']
                self.assertEqual(len(merges),1)
                self.assertTrue(merges[0].inputs['Selection'].is_linked)
                self.assertAlmostEqual(merges[0].inputs['Distance'].default_value,.1)
                # Keep the original shared construction and explicit weld and
                # orientation diagnostics, rather than replacing checked output
                # with an unchecked local fill.
                self.assertTrue(any(n.type=='GROUP' and n.node_tree.get('st_unchecked_separator_construction')
                                    for n in tree.nodes))
                fields={n.inputs['Name'].default_value for n in tree.nodes
                        if n.bl_idname=='GeometryNodeStoreNamedAttribute' and not n.inputs['Name'].is_linked}
                self.assertTrue({'st_separator_oriented','st_separator_agreement',
                    'st_separator_expected_weld','st_separator_actual_weld','st_separator_valid'} <= fields)

    def test_typed_records_cover_both_patch_branches_and_preserve_strip_edges(self):
        patch=self.blocks[assets.ORIGINAL_PATCH_PROVIDER]
        strip=self.blocks[assets.ORIGINAL_STRIP_PROVIDER]
        root=self.blocks[assets.ORIGINAL_SOLVER]
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        flags=tuple(tree.get('blosm_typed_provider_records',0) for tree in (root,patch,strip,display))
        self.assertIn(flags,((0,)*4,(1,)*4))
        if not flags[0]:self.skipTest('No separate provider record components in this revision')

        def roles(tree):
            found=set()
            for node in tree.nodes:
                if node.bl_idname=='GeometryNodeMeshToPoints' and node.label.startswith('Isolated record transport role '):
                    self.assertEqual(node.mode,'VERTICES')
                    self.assertTrue(node.inputs['Mesh'].is_linked)
                    self.assertTrue(node.outputs['Points'].is_linked)
                    found.add(int(node.label.rsplit(' ',1)[1]))
            return found

        patch_groups=list(self.reachable_groups(patch).values())
        dispatches=[g for g in patch_groups if g.get('st_complete_constructor_dispatch') and
                    any(n.type=='GROUP' and n.label=='Complete ordinary or paired constructor' for n in g.nodes)]
        self.assertTrue(dispatches,'Typed complete-constructor capability has no active dispatch')
        for dispatch in dispatches:
            for label in ('Unchanged original patch provider','Complete ordinary or paired constructor'):
                call,=[n for n in dispatch.nodes if n.type=='GROUP' and n.label==label]
                self.assertTrue(call.node_tree.get('st_typed_provider_records'))
                self.assertEqual(roles(call.node_tree),{2,3},'Both live patch constructors must isolate their records')
            self.assertEqual(roles(dispatch),{4},'Internal-region status must use the point-cloud component')
        strip_groups=list(self.reachable_groups(strip).values())
        packers=[g for g in strip_groups if g.get('st_typed_provider_records') and
                 g.name.startswith('Prochitecture Streets / Cached original strip owner v2')]
        self.assertTrue(packers)
        for packer in packers:
            self.assertEqual(roles(packer),{2,3},'Role 4 carries end-profile edges and must remain a mesh')

    def test_declared_separator_parity_reaches_every_chained_assembly(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        flags=tuple(root.get('blosm_separator_parity_preparation',0) for root in (*roots,display))
        self.assertIn(flags,((0,)*4,(1,)*4))
        if self.require_separator_parity:self.assertEqual(flags,(1,)*4)
        if not flags[0]:self.skipTest('No native separator parity preparation in this revision')
        groups=list(self.reachable_groups(roots[0]).values())
        assemblies=[tree for tree in groups
                    if tree.name.startswith('Prochitecture Streets / Original closed separator assembly v3')
                    and 'Chained boundaries' in tree.name]
        self.assertTrue(assemblies,'Capability metadata has no reachable chained assembly')
        preparations=set()
        for assembly in assemblies:
            self.assertEqual(assembly.get('st_separator_parity_preparation'),1)
            prepare,=[n for n in assembly.nodes
                      if n.label=='Conforming separator parity before unchanged fill']
            preparations.add(prepare.node_tree)
            inp,=[n for n in assembly.nodes if n.type=='GROUP_INPUT']
            links=list(inp.outputs['Pavement'].links)
            self.assertEqual(len(links),1,'Some pavement consumers bypass the parity preparation')
            self.assertEqual(links[0].to_socket,prepare.inputs['Pavement'])
            boundary,=[n for n in assembly.nodes if n.label=='Recover the actual declared pavement opening']
            fill,=[n for n in assembly.nodes if n.label=='Direct quad fill with exact source boundary']
            for call in (boundary,fill):
                link,=call.inputs['Pavement'].links
                self.assertEqual(link.from_socket,prepare.outputs['Geometry'])
        self.assertEqual(len(preparations),1,'Every network must use the same preparation group')
        preparation=preparations.pop()
        self.assertEqual(preparation.get('st_separator_parity_policy'),
                         'sequential_declared_holes_protect_all_other_openings_v1')
        strips=[tree for tree in self.reachable_groups(preparation).values()
                if tree.get('st_separator_parity_policy')=='one_same_owner_quad_strip_to_unprotected_boundary_v1']
        self.assertEqual(len(strips),1)
        scope_flags=tuple(root.get('blosm_separator_parity_scope',0) for root in (*roots,display))
        self.assertIn(scope_flags,((0,)*4,(2,)*4))
        if self.require_scoped_separator_parity:self.assertEqual(scope_flags,(2,)*4)
        if scope_flags[0]:
            self.assertEqual(strips[0].get('st_separator_parity_scope'),
                             'all_faces_of_side_incident_carriageway_owners_v2')
            labels={n.label for n in strips[0].nodes}
            self.assertIn('Owners touching this candidate side',labels)
            self.assertIn('Compact graph slot for this original face',labels)
        prohibited={'GeometryNodeTriangulate','GeometryNodeSubdivideMesh','GeometryNodeMeshBoolean'}
        self.assertFalse(any(n.bl_idname in prohibited for tree in strips for n in tree.nodes),
                         'The local parity correction must retain direct quads')

    def test_checked_chain_contacts_preserve_primary_and_repeat_original_guards(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        flags=tuple(root.get('blosm_checked_chain_contacts',0) for root in (*roots,display))
        self.assertIn(flags,((0,)*4,(1,)*4))
        if self.require_checked_chain_contacts:self.assertEqual(flags,(1,)*4)
        if not flags[0]:self.skipTest('No checked logical contact alternative in this revision')
        groups=self.reachable_groups(roots[1]).values()
        choices=[g for g in groups if g.get('st_checked_circle_contact_alternative')==1]
        self.assertEqual(len(choices),1,'All logical surfaces must share one checked selection')
        choice=choices[0]
        primary,=[n for n in choice.nodes if n.label=='Retain the complete valid original outline']
        circle,=[n for n in choice.nodes if n.label=='Same radius and guards, existing circular-contact endpoint']
        signature=lambda tree:[(s.identifier,s.name,s.in_out,s.socket_type)
                               for s in tree.interface.items_tree if s.item_type=='SOCKET']
        self.assertEqual(signature(choice),signature(primary.node_tree))
        self.assertEqual(signature(choice),signature(circle.node_tree))
        for group in (primary.node_tree,circle.node_tree):
            guards=[n for n in group.nodes if n.label=='Complete original curved carriageway footprints']
            self.assertEqual(len(guards),1,'Both alternatives must retain the original ribbon guard')
        selector,=[n for n in choice.nodes if n.label.startswith('Only replace an intruding outline with a fully checked')]
        link,=selector.inputs['Switch'].links
        self.assertEqual(link.from_socket,primary.outputs['Valid'])
        self.assertFalse(selector.inputs['True'].default_value)
        gate=selector.inputs['False'].links[0].from_node
        self.assertEqual(gate.inputs['Switch'].links[0].from_socket,primary.outputs['Footprints Checked'])
        self.assertEqual(gate.inputs['True'].links[0].from_socket,circle.outputs['Valid'])
        out,=[n for n in choice.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output]
        for socket in out.inputs:
            if socket.is_unavailable or not socket.is_linked:continue
            switch=socket.links[0].from_node
            self.assertEqual(switch.inputs['Switch'].links[0].from_socket,selector.outputs['Output'])
            self.assertEqual(switch.inputs['False'].links[0].from_socket,primary.outputs[socket.name])
            self.assertEqual(switch.inputs['True'].links[0].from_socket,circle.outputs[socket.name])

    def test_continuous_core_deformation_reaches_profiles_and_stays_ordinary(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        flags=tuple(root.get('blosm_corner_core_homotopy',0) for root in (*roots,display))
        self.assertIn(flags,((0,)*4,(2,)*4))
        if self.require_core_homotopy:self.assertEqual(flags,(2,)*4)
        if not flags[0]:self.skipTest('No continuous ordinary core deformation in this revision')
        recipes=[tree for tree in self.reachable_groups(roots[1]).values()
                 if tree.get('st_corner_pattern_recipe')]
        self.assertEqual({tree['st_corner_pattern_family'] for tree in recipes},{'ordinary','multi_site'})
        for recipe in recipes:
            with self.subTest(recipe=recipe.name):
                primitives=[tree for tree in self.reachable_groups(recipe).values()
                            if tree.get('st_ordinary_core_homotopy')]
                if recipe['st_corner_pattern_family']=='multi_site':
                    self.assertFalse(primitives,'Ordinary core deformation leaked into multi-site construction')
                    self.assertFalse(recipe.get('st_corner_core_homotopy',0))
                    continue
                self.assertEqual(recipe.get('st_corner_core_homotopy'),2)
                self.assertTrue(primitives,'Declared core deformation has no live constructor')
                for tree in primitives:
                    self.assertEqual(tree.get('st_ordinary_core_homotopy'),2)
                    self.assertEqual(tree.get('st_ordinary_guide_policy'),'guarded_cell_retry_v1')
                    deformation,=[n for n in tree.nodes if n.label=='Continuous convex-core guide deformation']
                    profile,=[n for n in tree.nodes if n.label=='Normal depth for equal-crown profile budget']
                    transport,=[n for n in tree.nodes if n.label=='Hub follows continuous core deformation']
                    link,=profile.inputs['Samples'].links
                    self.assertEqual(link.from_socket,deformation.outputs['Samples'],
                                     'Physical profile depth must use the deformed guides')
                    link,=transport.inputs['Hub'].links
                    self.assertEqual(link.from_socket,deformation.outputs['Hub'])
                    for key in ('Samples','Mouths','Rows','Applied'):
                        self.assertTrue(deformation.inputs[key].is_linked)
                    self.assertTrue(transport.outputs['Geometry'].is_linked)
                    self.assertTrue(any(n.label=='Corner-pattern shape admission' for n in tree.nodes))
                    self.assertTrue(any(n.bl_idname=='GeometryNodeInputNamedAttribute' and
                        n.inputs['Name'].default_value=='st_cj_half' for n in deformation.node_tree.nodes),
                        'Separator boundary eligibility was lost')

    def test_guide_correction_is_live_and_limited_to_ordinary_recipes(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        weights=tuple(root.get('blosm_corner_guide_blend',0.) for root in (*roots,display))
        self.assertEqual(len(set(weights)),1,'Inconsistent guide policy across exported roots')
        if self.require_guide_correction:self.assertGreater(weights[0],0.)
        if not weights[0]:self.skipTest('No reversed-mouth guide correction in this revision')
        self.assertLess(weights[0],1.)
        policy=roots[0].get('blosm_corner_guide_policy')
        self.assertIn(policy,('negative_mouth_turn_v1','guarded_cell_retry_v1'))
        for root in (*roots,display):
            self.assertEqual(root.get('blosm_corner_guide_policy'),policy)
        recipes=[tree for tree in self.reachable_groups(roots[1]).values() if tree.get('st_corner_pattern_recipe')]
        self.assertEqual({tree['st_corner_pattern_family'] for tree in recipes},{'ordinary','multi_site'})
        for recipe in recipes:
            with self.subTest(recipe=recipe.name):
                groups=list(self.reachable_groups(recipe).values())
                calls=[(tree,node) for tree in groups for node in tree.nodes
                       if node.type=='GROUP' and node.label in ('Correct only reversed guide turns at mouths',
                                                              'Correct full-band guides after rejected cells')]
                if recipe['st_corner_pattern_family']=='multi_site':
                    self.assertFalse(calls,'Ordinary guide correction leaked into multi-site construction')
                    self.assertFalse(recipe.get('st_corner_guide_blend',0.))
                    continue
                self.assertEqual(recipe.get('st_corner_guide_blend'),weights[0])
                self.assertEqual(recipe.get('st_corner_guide_policy'),policy)
                self.assertTrue(calls,'Declared guide policy has no live construction')
                for tree,node in calls:
                    self.assertEqual(tree.get('st_ordinary_guide_policy'),policy)
                    for key in ('Samples','Mouths','Rows'):self.assertTrue(node.inputs[key].is_linked)
                    self.assertTrue(node.outputs['Samples'].is_linked)
                    self.assertTrue(any(n.label=='Corner-pattern shape admission' for n in tree.nodes),
                                    'Guide correction lost its caller shape checks')
                    self.assertTrue(any(n.bl_idname=='GeometryNodeInputNamedAttribute'
                        and n.inputs['Name'].default_value=='st_cj_half' for n in node.node_tree.nodes),
                        'Half-profile boundary eligibility was lost')

    def test_guarded_profile_alternative_keeps_original_admission(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        flags=tuple(root.get('blosm_corner_capped_profile',0) for root in (*roots,display))
        self.assertIn(flags,((0,)*4,(1,)*4))
        if self.require_guarded_profile:self.assertEqual(flags,(1,)*4)
        if not flags[0]:self.skipTest('No capped-profile candidate in this revision')
        policy=roots[0].get('blosm_corner_guide_policy')
        self.assertIn(policy,('negative_mouth_turn_v1','guarded_cell_retry_v1'))
        contact_retry=roots[0].get('blosm_corner_contact_retry',0)
        if self.require_guarded_profile:
            self.assertEqual(policy,'guarded_cell_retry_v1')
            self.assertEqual(contact_retry,1)
        for root in (*roots,display):
            self.assertEqual(root.get('blosm_corner_guide_policy'),policy)
            self.assertEqual(root.get('blosm_corner_contact_retry',0),contact_retry)
        recipes=[g for g in self.reachable_groups(roots[1]).values() if g.get('st_corner_pattern_recipe')]
        self.assertEqual({g['st_corner_pattern_family'] for g in recipes},{'ordinary','multi_site'})
        for recipe in recipes:
            with self.subTest(recipe=recipe.name):
                groups=list(self.reachable_groups(recipe).values())
                retries=[g for g in groups if g.get('st_corner_guarded_retry')]
                contacts=[g for g in groups if g.get('st_contact_corner_weight_policy')]
                profiles=[g for g in groups if g.get('st_height_profile_experiment')]
                if recipe['st_corner_pattern_family']=='multi_site':
                    self.assertFalse(retries or contacts or profiles,'Ordinary corrections leaked into a multi-site recipe')
                    continue
                self.assertEqual(bool(retries),policy=='guarded_cell_retry_v1')
                self.assertEqual(bool(contacts),bool(contact_retry))
                self.assertTrue(profiles,'Declared profile correction has no live dependencies')
                for tree in retries:
                    original,=[n for n in tree.nodes if n.type=='GROUP' and n.label=='Original accepted corner layout']
                    alternate,=[n for n in tree.nodes if n.type=='GROUP' and n.label=='Bounded retry after rejected corner cells']
                    out,=[n for n in tree.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output]
                    # Verify the actual selector wiring, not just capability tags.
                    for item in tree.interface.items_tree:
                        if item.item_type!='SOCKET' or item.in_out!='OUTPUT':continue
                        link,=out.inputs[item.name].links; select=link.from_node
                        self.assertEqual(select.bl_idname,'GeometryNodeSwitch')
                        for key,expected in (('Switch',original.outputs['Valid']),
                                             ('True',original.outputs[item.name]),
                                             ('False',alternate.outputs[item.name])):
                            incoming,=select.inputs[key].links
                            self.assertEqual(incoming.from_socket,expected)
                for tree in contacts:
                    checks=[n for n in tree.nodes if n.type=='GROUP' and n.label in (
                        'Unchanged complete 50 micrometre polygon bound','Retry unchanged boundary check per corner')]
                    self.assertEqual(len(checks),2)
                    self.assertTrue(all(n.outputs['Valid'].is_linked for n in checks))
                    if recipe.get('st_direct_corner_records'):
                        self.assertEqual(recipe.get('st_corner_direct_source_before_contact_retry'),1)
                        original,=[n for n in tree.nodes if n.type=='GROUP' and
                                   n.label=='Same direct quad construction with mesh-only knots']
                        alternate,=[n for n in tree.nodes if n.type=='GROUP' and
                                    n.label=='Retry direct layout with independent corner weights']
                        self.assertTrue(original.node_tree.get('st_direct_corner_records'))
                        self.assertTrue(alternate.node_tree.get('st_direct_corner_records'))
                        self.assertTrue(alternate.inputs['Corner Weights'].is_linked)
                for tree in profiles:
                    self.assertEqual(tree.get('st_height_profile_experiment'),'physical_depth_capped_quadratic_v1')
                    valid,=[n for n in tree.nodes if n.bl_idname=='GeometryNodeStoreNamedAttribute'
                            and n.inputs['Name'].default_value=='st_ch_valid']
                    self.assertTrue(valid.inputs['Value'].is_linked)

    def test_directional_refinement_preserves_ordinary_publication_guards(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        flags=tuple(root.get('blosm_corner_strip_refinement',0) for root in roots)
        self.assertIn(flags,((0,0,0),(1,1,1)))
        if self.require_strip_refinement:self.assertEqual(flags,(1,1,1))
        if not flags[0]:self.skipTest('No directional ordinary-quad refinement in this revision')
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        self.assertEqual(display.get('blosm_corner_strip_refinement'),1)
        spacing=tuple(roots[0].get(key,default) for key,default in (
            ('blosm_corner_radial_edge_budget',6.5),('blosm_corner_refinement_limit',7.)))
        for value in spacing:
            self.assertGreater(value,1.)
            self.assertLess(value,8.,'Construction budget must not relax the acceptance ceiling')
        for root in (*roots,display):
            self.assertEqual(root.get('blosm_corner_ordinary_core_repair'),1)
            self.assertEqual(tuple(root.get(key,default) for key,default in (
                ('blosm_corner_radial_edge_budget',6.5),('blosm_corner_refinement_limit',7.))),spacing)
        recipes=[tree for tree in self.reachable_groups(roots[1]).values()
                 if tree.get('st_corner_pattern_recipe')]
        self.assertEqual({tree['st_corner_pattern_family'] for tree in recipes},{'ordinary','multi_site'})
        for recipe in recipes:
            with self.subTest(recipe=recipe.name):
                groups=list(self.reachable_groups(recipe).values())
                layouts=[tree for tree in groups if tree.get('st_quad_strip_refinement')]
                if recipe['st_corner_pattern_family']=='multi_site':
                    self.assertFalse(recipe.get('st_corner_strip_refinement',0))
                    self.assertFalse(layouts,'Ordinary refinement leaked into multi-site construction')
                    continue
                self.assertEqual(recipe.get('st_corner_strip_refinement'),1)
                self.assertEqual(recipe.get('st_corner_ordinary_core_repair'),1)
                self.assertEqual((recipe.get('st_corner_radial_edge_budget',6.5),
                                  recipe.get('st_corner_refinement_limit',7.)),spacing)
                grading=[tree for tree in groups if tree.name.startswith(
                    'Prochitecture Streets / Boundary-budgeted local corner band grading v2')
                    and ' / Geometric radial spacing v1' in tree.name]
                self.assertTrue(grading,'Ordinary recipe lost its geometric radial spacing')
                for tree in grading:
                    match=re.search(r' / First-edge budget ([0-9]+(?:\.[0-9]+)?)',tree.name)
                    self.assertEqual(float(match[1]) if match else 6.5,spacing[0])
                self.assertTrue(layouts,'Declared capability has no live refinement')
                for layout in layouts:
                    checks=[node for node in layout.nodes if node.label=='Corner-pattern shape admission']
                    self.assertEqual(len(checks),1)
                    links=list(checks[0].inputs['Geometry'].links)
                    self.assertEqual(len(links),1)
                    refinement=links[0].from_node
                    if layout.get('st_corner_depth_floor') and roots[0].get('blosm_corner_profile_slope_reserve'):
                        # The reserve adds only an attribute after refinement.
                        # Both branches must carry that same refined geometry,
                        # and only the actual recovery may add the attribute.
                        self.assertEqual(refinement.bl_idname,'GeometryNodeSwitch')
                        self.assertEqual(refinement.input_type,'GEOMETRY')
                        proposal,=[node for node in layout.nodes if node.type=='GROUP'
                                   and node.label=='Recover profile room after infeasible guide compression']
                        self.assertEqual(refinement.inputs['Switch'].links[0].from_socket,
                                         proposal.outputs['Applied'])
                        original=refinement.inputs['False'].links[0].from_socket
                        stored=refinement.inputs['True'].links[0].from_node
                        self.assertEqual(stored.bl_idname,'GeometryNodeStoreNamedAttribute')
                        self.assertEqual(stored.data_type,'FLOAT');self.assertEqual(stored.domain,'POINT')
                        self.assertEqual(stored.inputs['Name'].default_value,'st_cq_profile_slope_reserve')
                        self.assertFalse(stored.inputs['Name'].is_linked)
                        self.assertFalse(stored.inputs['Value'].is_linked)
                        self.assertAlmostEqual(stored.inputs['Value'].default_value,
                                               roots[0]['blosm_corner_profile_slope_reserve'])
                        self.assertFalse(stored.inputs['Selection'].is_linked)
                        self.assertTrue(stored.inputs['Selection'].default_value)
                        self.assertEqual(stored.inputs['Geometry'].links[0].from_socket,original)
                        refinement=original.node
                    self.assertEqual(refinement.type,'GROUP')
                    self.assertEqual(refinement.label,'Conforming directional quads with exact entrance edges')
                    self.assertTrue(refinement.inputs['Geometry'].is_linked)
                    self.assertIn('Profile records',refinement.node_tree.name)
                    match=re.search(r'Direct quad strip refinement v1 / [0-9]+ / ([0-9.]+) /',
                                    refinement.node_tree.name)
                    self.assertIsNotNone(match,'Unrecognized refinement kernel revision')
                    self.assertEqual(float(match[1]),spacing[1])
                    live=self.reachable_groups(refinement.node_tree).values()
                    self.assertFalse(any(node.bl_idname=='GeometryNodeTriangulate'
                                         for tree in live for node in tree.nodes))
                    self.assertTrue(any(node.bl_idname=='GeometryNodeMeshGrid'
                                        for node in refinement.node_tree.nodes))

    def test_shared_kernels_preserve_checked_callers(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        flags=tuple((root.get('blosm_shared_separator_construction',0),
                     root.get('blosm_shared_sector_coordinates',0)) for root in roots)
        for column in zip(*flags):self.assertIn(column,((0,0,0),(1,1,1)))
        if self.require_shared_kernels:self.assertEqual(flags,((1,1),)*3)
        if flags[0]==(0,0):self.skipTest('No shared construction kernels in this revision')
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        for key,flag in zip(('blosm_shared_separator_construction','blosm_shared_sector_coordinates'),flags[0]):
            self.assertEqual(display.get(key,0),flag)
        groups={tree.as_pointer():tree for root in roots for tree in self.reachable_groups(root).values()}
        separators=[tree for tree in groups.values() if tree.get('st_shared_separator_construction')]
        coordinates=[tree for tree in groups.values() if tree.get('st_shared_coordinates')]
        if flags[0][0]:self.assertTrue(separators)
        if flags[0][1]:self.assertTrue(coordinates)
        for tree in coordinates:
            calls=[node for node in tree.nodes if node.type=='GROUP' and node.node_tree and
                   node.label=='Shared original coordinate arithmetic']
            self.assertEqual(len(calls),1)
            self.assertEqual(calls[0].node_tree.get('st_shared_coordinate_core'),1)
            if roots[0].get('blosm_original_interpolation'):
                self.assertIn(' / Exact entrance rays v1',calls[0].node_tree.name)
        for tree in groups.values():
            for node in tree.nodes:
                if node.type=='GROUP' and node.node_tree and node.node_tree.get('st_unchecked_separator_construction'):
                    self.assertEqual(tree.get('st_shared_separator_construction'),1,
                                     'Raw separator core bypasses its checked publication wrapper')
            if flags[0][0] and tree.name.startswith('Prochitecture Streets / Curved capped separator quads v1'):
                self.assertTrue(tree.get('st_shared_separator_construction') or
                                tree.get('st_unchecked_separator_construction'),tree.name)

    def test_cycle_cores_reach_only_multi_site_recipes(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        flags=tuple(root.get('blosm_corner_cycle_cores',0) for root in roots)
        self.assertIn(flags,((0,0,0),(1,1,1)))
        if self.require_cycle_cores:self.assertEqual(flags,(1,1,1))
        if not flags[0]:self.skipTest('No cycle-core capability')
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        self.assertEqual(display.get('blosm_corner_cycle_cores'),1)
        recipes=[tree for tree in self.reachable_groups(roots[1]).values()
                 if tree.get('st_corner_pattern_recipe')]
        self.assertTrue(recipes)
        self.assertEqual({tree['st_corner_pattern_family'] for tree in recipes},{'ordinary','multi_site'})
        for recipe in recipes:
            with self.subTest(recipe=recipe.name):
                groups=list(self.reachable_groups(recipe).values())
                cores=[tree for tree in groups if tree.name.startswith('Prochitecture Streets / Physical-site corner core v1')]
                if recipe['st_corner_pattern_family']=='ordinary':
                    self.assertFalse(recipe.get('st_corner_cycle_cores',0))
                    self.assertFalse(cores)
                    continue
                self.assertEqual(recipe.get('st_corner_cycle_cores'),1)
                self.assertTrue(cores)
                for core in cores:
                    calls=[node for node in core.nodes if node.type=='GROUP' and node.node_tree and
                           node.label=='Geometric spacing next to short cycle-boundary segments']
                    self.assertEqual(len(calls),1,'A multi-site core lost its graded-cycle branch')
                    self.assertEqual(calls[0].node_tree.name,'Prochitecture Streets / Simple-cycle graded core v2')
                    self.assertTrue(calls[0].inputs['Input Valid'].is_linked)
                    iterations=[node for node in calls[0].node_tree.nodes
                                if node.bl_idname=='GeometryNodeRepeatInput' and
                                node.label=='Twelve bounded geometric-spacing steps']
                    self.assertEqual(len(iterations),1)
                    self.assertFalse(iterations[0].inputs['Iterations'].is_linked)
                    self.assertEqual(iterations[0].inputs['Iterations'].default_value,12)

    def test_direct_corner_sources_omit_preliminary_centers(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        flags=tuple(root.get('blosm_direct_corner_records',0) for root in roots)
        self.assertIn(flags,((0,0,0),(1,1,1)))
        if not flags[0]:self.skipTest('No direct corner-source capability')
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        self.assertEqual(display.get('blosm_direct_corner_records'),1)
        sources={node.node_tree for tree in self.reachable_groups(roots[1]).values()
                 for node in tree.nodes if node.type=='GROUP' and node.node_tree and
                 node.label=='Original analytic boundary and ray samples; no interior descent'}
        self.assertTrue(sources,'Declared direct records have no active local corner source')
        for source in sources:
            with self.subTest(source=source.name):
                groups=self.reachable_groups(source).values()
                layouts=[tree for tree in groups if any(node.bl_idname=='GeometryNodeMeshGrid'
                         and node.label=='Direct regular half-profile band' for node in tree.nodes)]
                self.assertTrue(layouts,'Missing direct entrance/boundary layout')
                for tree in layouts:
                    self.assertEqual(tree.get('st_direct_corner_records'),1)
                    self.assertFalse(any(node.label=='Direct center block between guide halves and shared spokes'
                                         for node in tree.nodes),'A source still builds an unused center')

    def test_batched_corner_guides_reach_every_exported_recipe(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        flags=tuple(root.get('blosm_corner_batched_guides',0) for root in roots)
        self.assertIn(flags,((0,0,0),(1,1,1)))
        if not flags[0]:self.skipTest('No batched corner-guide capability')
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        self.assertEqual(display.get('blosm_corner_batched_guides'),1)
        recipes=[tree for tree in self.reachable_groups(roots[1]).values()
                 if tree.get('st_corner_pattern_recipe')]
        self.assertTrue(recipes,'Declared guide batching has no live corner recipes')
        for recipe in recipes:
            with self.subTest(recipe=recipe.name):
                self.assertEqual(recipe.get('st_corner_batched_guides'),1)
                primitives=[tree for tree in self.reachable_groups(recipe).values()
                            if tree.name.startswith('Prochitecture Streets / Local corner pattern v1')]
                self.assertTrue(primitives,'Corner recipe does not call a local pattern')
                for primitive in primitives:
                    self.assertEqual(primitive.get('st_corner_batched_guides'),1,
                                     'A live pattern branch still uses separate guide preparation: '+primitive.name)

    def test_corner_builder_entry_points_share_the_published_recipes(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        flags=tuple(root.get('blosm_corner_builder_sources_reconciled',0) for root in roots)
        self.assertIn(flags,((0,0,0),(1,1,1)))
        if not flags[0]:self.skipTest('No reconciled builder-source capability')
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        self.assertEqual(display.get('blosm_corner_builder_sources_reconciled'),1)
        live=self.reachable_groups(roots[1])
        recipes={g.as_pointer() for g in live.values() if g.get('st_corner_pattern_recipe')}
        sources={g['st_complete_constructor_source'] for g in live.values()
                 if 'st_complete_constructor_source' in g}
        self.assertTrue(sources,'The runtime export lost its canonical builder entry point')
        self.assertTrue(recipes)
        for source in sources:
            with self.subTest(source=source.name):
                self.assertTrue(source.name.startswith('Prochitecture Streets / Cached original patch owner v2'))
                shared={g.as_pointer() for g in self.reachable_groups(source).values()
                        if g.get('st_corner_pattern_recipe')}
                self.assertTrue(shared,'Canonical source still uses superseded frames')
                self.assertLessEqual(shared,recipes,'Canonical source duplicates a live recipe')

    def test_initial_corner_sources_expose_only_sample_and_control_records(self):
        roots=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                                            assets.ORIGINAL_STRIP_PROVIDER))
        flags=tuple(root.get('blosm_corner_initial_samples_only',0) for root in roots)
        self.assertIn(flags,((0,0,0),(1,1,1)))
        if not flags[0]:self.skipTest('No initial-only corner source capability')
        display=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,))[assets.ORIGINAL_FAILURE_DISPLAY]
        self.assertEqual(display.get('blosm_corner_initial_samples_only'),1)
        sources={node.node_tree for tree in self.reachable_groups(roots[1]).values()
                 for node in tree.nodes if node.type=='GROUP' and node.node_tree and
                 node.label=='Original analytic boundary and ray samples; no interior descent'}
        self.assertTrue(sources,'Declared sample-only recipe has no active corner sources')
        # A checked identity-station wrapper publishes the same record contract.
        # Verify its gate, both unchanged source interfaces and the explicit
        # identity diagnostics before applying the existing source checks to
        # both branches. It must never turn into a surface-generating shortcut.
        expanded=set()
        for tree in sources:
            if not tree.get('st_contact_identity_fallback'):
                expanded.add(tree)
                continue
            self.assertEqual(tuple(root.get('blosm_corner_contact_identity_fallback',0) for root in roots),(1,1,1))
            self.assertEqual(tree.get('st_initial_samples_only'),1)
            self.assertGreater(tree.get('st_contact_identity_station_calls',0),0)
            first,=[node for node in tree.nodes if node.type=='GROUP' and node.label=='Existing checked contact stations']
            second,=[node for node in tree.nodes if node.type=='GROUP' and node.label=='Original stations with unchanged complete boundary checks']
            signature=lambda group:[(s.name,s.socket_type,s.in_out) for s in group.interface.items_tree if s.item_type=='SOCKET']
            self.assertEqual(signature(tree),signature(first.node_tree))
            self.assertEqual(signature(tree),signature(second.node_tree))
            self.assertNotEqual(first.node_tree.as_pointer(),second.node_tree.as_pointer())
            out,=[node for node in tree.nodes if node.type=='GROUP_OUTPUT' and node.is_active_output]
            gates=set()
            for item in tree.interface.items_tree:
                if item.item_type!='SOCKET' or item.in_out!='OUTPUT':continue
                link,=out.inputs[item.name].links
                choice=link.from_node
                self.assertEqual(choice.bl_idname,'GeometryNodeSwitch')
                original,=choice.inputs['False'].links
                self.assertEqual(original.from_socket,first.outputs[item.name])
                proposed,=choice.inputs['True'].links
                if item.socket_type=='NodeSocketGeometry':
                    stored=proposed.from_node
                    self.assertEqual(stored.bl_idname,'GeometryNodeStoreNamedAttribute')
                    self.assertEqual(stored.inputs['Name'].default_value,
                                     '_cw_weight' if item.name=='Corner Weights' else 'st_cq_contact_blend')
                    self.assertFalse(stored.inputs['Value'].is_linked)
                    self.assertEqual(stored.inputs['Value'].default_value,0.)
                    proposed,=stored.inputs['Geometry'].links
                self.assertEqual(proposed.from_socket,second.outputs[item.name])
                gate,=choice.inputs['Switch'].links;gates.add(gate.from_socket)
            gate,=gates
            selector=gate.node
            self.assertEqual((selector.bl_idname,selector.input_type),('GeometryNodeSwitch','BOOLEAN'))
            self.assertFalse(selector.inputs['True'].is_linked)
            self.assertFalse(selector.inputs['True'].default_value)
            for socket,source in ((selector.inputs['Switch'],first),(selector.inputs['False'],second)):
                link,=socket.links
                product=link.from_node
                self.assertEqual((product.bl_idname,product.operation),('ShaderNodeMath','MULTIPLY'))
                self.assertEqual({port.links[0].from_socket for port in list(product.inputs)[:2]},
                                 {source.outputs['Inputs Valid'],source.outputs['Seams Valid']})
            expanded.update((first.node_tree,second.node_tree))
        sources=expanded
        expected={'Initial Geometry','Plan','Inputs Valid','Seams Valid','Radial Count'}
        for tree in sources:
            with self.subTest(source=tree.name):
                self.assertEqual(tree.get('st_initial_samples_only'),1)
                outputs=set(expected)
                if tree.get('st_contact_corner_weight_policy'):
                    # The retry adds a control-table output to the five sample
                    # outputs. It is the very table already consumed by both
                    # retry branches, not an additional surface construction.
                    self.assertEqual(tree['st_contact_corner_weight_policy'],
                                     'retry_rejected_boundary_per_corner_v1')
                    self.assertEqual(roots[0].get('blosm_corner_contact_retry'),1)
                    outputs.add('Corner Weights')
                    out,=[n for n in tree.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output]
                    link,=out.inputs['Corner Weights'].links
                    stored=link.from_node
                    self.assertEqual(stored.bl_idname,'GeometryNodeStoreNamedAttribute')
                    self.assertEqual(stored.inputs['Name'].default_value,'_cw_weight')
                    self.assertEqual((stored.domain,stored.data_type),('POINT','FLOAT'))
                    callers=[n for n in tree.nodes if n.type=='GROUP' and n.label in (
                        'Retry direct layout with independent corner weights',
                        'Retry unchanged boundary check per corner')]
                    self.assertEqual(len(callers),2)
                    for caller in callers:
                        binding,=caller.inputs['Corner Weights'].links
                        self.assertEqual(binding.from_socket,link.from_socket)
                self.assertEqual({s.name for s in tree.interface.items_tree
                                  if s.item_type=='SOCKET' and s.in_out=='OUTPUT'},outputs)
                self.assertFalse(any(n.label=='Pinned direct-quad quality correction' for n in tree.nodes))

    def test_original_interpolation_agrees_across_assets_and_source_bridges(self):
        from itertools import product
        blocks=tuple(self.blocks[key] for key in
                     (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER))
        enabled=assets.validate_original_interpolation(*blocks)
        if self.require_interpolation:self.assertTrue(enabled)
        copies=tuple(block.copy() for block in blocks)
        try:
            for flags in product((0,1),repeat=3):
                for block,flag in zip(copies,flags):block['blosm_original_interpolation']=flag
                if len(set(flags))==1:self.assertEqual(assets.validate_original_interpolation(*copies),bool(flags[0]))
                else:
                    with self.assertRaisesRegex(ValueError,'interpolation'):assets.validate_original_interpolation(*copies)
            for block in copies:block['blosm_original_interpolation']=1
            for block,validate in zip(copies,(self.check_root,assets.validate_original_patch_provider,
                                              assets.validate_original_strip_provider)):
                block['blosm_original_interpolation']=2
                with self.assertRaisesRegex(ValueError,'interpolation'):validate(block)
                with self.assertRaisesRegex(ValueError,'interpolation'):assets.validate_original_interpolation(*copies)
                block['blosm_original_interpolation']=1
        finally:
            for block in copies:bpy.data.node_groups.remove(block)
        prefix='Prochitecture Streets / Canonical Section source views v1'
        for block in blocks:
            names=set(self.reachable_groups(block))
            expected=prefix+(' / Exact or smooth v1' if enabled else '')
            self.assertEqual({name for name in names if name.startswith(prefix)},{expected})
            if enabled:
                self.assertTrue(any(name.startswith('Prochitecture Streets / Original exact polyline table v1') for name in names))

    def test_interpolated_patch_constructor_reuses_the_packaged_group(self):
        root=self.blocks[assets.ORIGINAL_PATCH_PROVIDER]
        if not root.get('blosm_original_interpolation'):
            self.skipTest('Legacy library does not contain the combined interpolation constructor')
        from blosm.way.gn.editable_input import curved_patch_provider_nodes as patch
        from blosm.way.gn.editable_input import curved_internal_patch_nodes as internal
        options=dict(internal_connectors=bool(root.get('blosm_internal_logical_ends')),
                     shared_logical=bool(root.get('blosm_shared_logical_owners')),
                     seam_connections=bool(root.get('blosm_original_seams')),
                     seam_transitions=bool(root.get('blosm_finite_seam_transitions')),
                     interpolation=True)
        before={tree.as_pointer() for tree in bpy.data.node_groups}
        core=patch.build(**options)
        self.assertIn(core.as_pointer(),before,'Constructor rebuilt the packaged shared group')
        self.assertLessEqual(len(core.name.encode('utf8')),200)
        self.assertIn(' / Exact or smooth v3',core['st_construction_label'])
        self.assertEqual(core.as_pointer(),patch.build(**options).as_pointer())
        if root.get('blosm_internal_cycle_regions'):
            extension=internal.extend(core)
            self.assertIn(extension.as_pointer(),before)
            self.assertEqual(extension.as_pointer(),internal.extend(core).as_pointer())
        self.assertEqual({tree.as_pointer() for tree in bpy.data.node_groups},before)

    def test_interpolated_failure_wrapper_keeps_its_own_strict_core(self):
        root=self.blocks[assets.ORIGINAL_SOLVER]
        if not root.get('blosm_original_interpolation'):
            self.skipTest('Legacy library predates core-specific failure-wrapper identities')
        from blosm.way.gn.editable_input import curved_failure_surface_nodes as failure
        prefix='Prochitecture Streets / Original local failure assembly candidate v6'
        # The standalone-Street extension inherits this prefix but is not the
        # failure constructor. Require the exact bounded constructor identity.
        candidates=[tree for tree in self.reachable_groups(root).values()
                    if re.fullmatch(re.escape(prefix)+r' / [0-9a-f]{16}',tree.name)]
        self.assertEqual(len(candidates),1)
        wrapper=candidates[0]
        accepted=[node for node in wrapper.nodes if node.type=='GROUP' and
                  node.label=='Keep the accepted complete construction verbatim']
        self.assertEqual(len(accepted),1)
        strict=accepted[0].node_tree
        self.assertEqual(failure.construction_identity(strict)[0],wrapper.name,
                         'Saved wrapper identity must match its upgraded strict core')
        before={tree.as_pointer() for tree in bpy.data.node_groups}
        self.assertEqual(failure.build(strict).as_pointer(),wrapper.as_pointer())
        self.assertEqual({tree.as_pointer() for tree in bpy.data.node_groups},before)
        clone=strict.copy();clone.name='Asset contract / distinct strict core'
        other=None
        try:
            other=failure.build(clone)
            self.assertNotEqual(other.as_pointer(),wrapper.as_pointer())
            calls=[node for node in other.nodes if node.type=='GROUP' and
                   node.label=='Keep the accepted complete construction verbatim']
            self.assertEqual(len(calls),1)
            self.assertEqual(calls[0].node_tree.as_pointer(),clone.as_pointer())
            self.assertEqual(failure.build(clone).as_pointer(),other.as_pointer())
            self.assertEqual(accepted[0].node_tree.as_pointer(),strict.as_pointer())
        finally:
            if other is not None:bpy.data.node_groups.remove(other)
            bpy.data.node_groups.remove(clone)

    def test_cleanup_contract_is_confined_to_controller_intake(self):
        root=self.blocks[assets.ORIGINAL_SOLVER]
        enabled=root.get('blosm_assembly_payload_cleanup',0)
        if self.require_assembly_payload_cleanup:self.assertIn(enabled,(1,2))
        if not enabled:self.skipTest('No assembly-intake cleanup in this revision')
        from blosm.way.gn.editable_input import curved_network_input_nodes as intake
        contract=json.loads(root['st_assembly_payload_discard'])
        self.assertIn(enabled,(1,2))
        fields=(intake.ASSEMBLY_PATTERN_CONSTRUCTION_FIELDS if enabled==2
                else intake.ASSEMBLY_CONSTRUCTION_FIELDS)
        self.assertEqual(contract,list(fields))
        self.assertEqual(tuple(contract),intake.validate_payload_discard(root,contract))
        # Provider-to-provider mouth metadata must survive; only the controller
        # may contain this optimization, including through nested group calls.
        pending=[self.blocks[assets.ORIGINAL_PATCH_PROVIDER],self.blocks[assets.ORIGINAL_STRIP_PROVIDER]]
        seen=set()
        while pending:
            tree=pending.pop()
            if tree.as_pointer() in seen:continue
            seen.add(tree.as_pointer())
            self.assertNotIn('st_assembly_payload_discard',tree)
            pending.extend(n.node_tree for n in tree.nodes if n.type=='GROUP' and n.node_tree)

    def test_cleanup_refuses_consumed_public_and_wildcard_attributes(self):
        from blosm.way.gn.editable_input import curved_network_input_nodes as intake
        root=self.blocks[assets.ORIGINAL_SOLVER]
        for names in (('st_cn_height',),('_qr_*',),('position',)):
            with self.subTest(names=names),self.assertRaises(ValueError):
                intake.validate_payload_discard(root,names)
        # This private name is consumed by the current exact identity weld.
        with self.assertRaisesRegex(ValueError,'consumed'):
            intake.validate_payload_discard(root,('_bridge_key',))

    def test_cleanup_refuses_dynamic_reader_names(self):
        from blosm.way.gn.editable_input import curved_network_input_nodes as intake
        tree=self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            reader=tree.nodes.new('GeometryNodeInputNamedAttribute')
            name=tree.nodes.new('FunctionNodeInputString')
            tree.links.new(name.outputs[0],reader.inputs['Name'])
            with self.assertRaisesRegex(ValueError,'Dynamic'):
                intake.validate_payload_discard(tree,('_qr_gradient',))
        finally:bpy.data.node_groups.remove(tree)

    def test_cleanup_protects_material_attribute_consumers(self):
        from blosm.way.gn.editable_input import curved_network_input_nodes as intake
        material=bpy.data.materials.new('Payload contract material consumer')
        try:
            material.use_nodes=True
            reader=material.node_tree.nodes.new('ShaderNodeAttribute')
            reader.attribute_name='_qr_gradient'
            with self.assertRaisesRegex(ValueError,'consumed'):
                intake.validate_payload_discard(self.blocks[assets.ORIGINAL_SOLVER],('_qr_gradient',))
        finally:bpy.data.materials.remove(material)

    def test_internal_region_assets_have_matching_capability_and_constructors(self):
        root=self.blocks[assets.ORIGINAL_SOLVER];patch=self.blocks[assets.ORIGINAL_PATCH_PROVIDER]
        enabled=assets.validate_original_internal_regions(root,patch)
        if self.require_internal_regions:self.assertTrue(enabled)
        if not enabled:self.skipTest('No internal-region capability')
        self.check_root(root);assets.validate_original_patch_provider(patch)
        names=set();pending=[root,patch]
        while pending:
            tree=pending.pop()
            if tree.name in names:continue
            names.add(tree.name);pending.extend(n.node_tree for n in tree.nodes if n.type=='GROUP' and n.node_tree)
        self.assertTrue(any(' / Shared internal separator ownership v1' in name for name in names))
        self.assertIn('Prochitecture Streets / Owned internal plateau partition v1',names)
        self.assertIn('Prochitecture Streets / Canonical internal cycle view v1',names)
        for block in (root,patch):
            copy=block.copy()
            try:
                socket=next(s for s in copy.interface.items_tree if s.item_type=='SOCKET' and s.in_out=='INPUT' and s.name=='Internal Cycles')
                copy.interface.remove(socket)
                with self.assertRaises(ValueError):
                    assets.validate_original_internal_regions(copy,patch) if block is root else assets.validate_original_internal_regions(root,copy)
            finally:bpy.data.node_groups.remove(copy)

    def test_internal_region_capability_cannot_be_partial_or_unknown(self):
        root=self.blocks[assets.ORIGINAL_SOLVER].copy();patch=self.blocks[assets.ORIGINAL_PATCH_PROVIDER].copy()
        try:
            for a,b in ((0,1),(1,0),(2,2)):
                root['blosm_internal_cycle_regions']=a;patch['blosm_internal_cycle_regions']=b
                with self.assertRaisesRegex(ValueError,'internal-region'):
                    assets.validate_original_internal_regions(root,patch)
            with self.assertRaisesRegex(ValueError,'internal-region'):self.check_root(root)
            with self.assertRaisesRegex(ValueError,'internal-region'):assets.validate_original_patch_provider(patch)
        finally:
            bpy.data.node_groups.remove(root);bpy.data.node_groups.remove(patch)

    def test_continuation_capability_requires_inputs_and_reachable_constructor(self):
        root=self.blocks[assets.ORIGINAL_SOLVER]
        if self.require_continuation_regions:self.assertEqual(root.get('blosm_continuation_separator_regions'),1)
        if not root.get('blosm_continuation_separator_regions'):self.skipTest('No continuation-region capability')
        self.check_root(root)
        names=set();pending=[root]
        while pending:
            tree=pending.pop()
            if tree.name in names:continue
            names.add(tree.name);pending.extend(n.node_tree for n in tree.nodes if n.type=='GROUP' and n.node_tree)
        self.assertTrue(any(' / Continuation regions v1 / Metric 32' in name for name in names))
        self.assertTrue(any('Crossing-to-free separator candidate v1 / Explicit region namespace / Metric heights 32' in name for name in names))
        for name in ('Continuation Regions','Continuation Owners','Continuation Count'):
            copy=root.copy()
            try:
                socket=next(s for s in copy.interface.items_tree if s.item_type=='SOCKET' and s.in_out=='INPUT' and s.name==name)
                copy.interface.remove(socket)
                with self.assertRaisesRegex(ValueError,name):self.check_root(copy)
            finally:bpy.data.node_groups.remove(copy)

    def test_unknown_continuation_capability_rejected(self):
        root=self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            root['blosm_continuation_separator_regions']=2
            with self.assertRaisesRegex(ValueError,'continuation-region'):self.check_root(root)
        finally:bpy.data.node_groups.remove(root)

    def test_continuations_require_crossing_regions(self):
        root=self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            root['blosm_continuation_separator_regions']=1;root['blosm_crossing_separator_regions']=0
            with self.assertRaisesRegex(ValueError,'Continuation regions require'):self.check_root(root)
        finally:bpy.data.node_groups.remove(root)

    @classmethod
    def setUpClass(cls):
        loader = BlenderAssetLoader(Path(blosm.__file__).parent/'assets', assets.catalog())
        # Exercise the actual universal loader's explicit mode consistently,
        # including later role requests made by individual contract tests.
        loader.load_many = partial(loader.load_many, link=cls.linked_library)
        cls.loader = loader
        cls.refs = (assets.ORIGINAL_SOLVER, assets.ORIGINAL_PATCH_PROVIDER,
                    assets.ORIGINAL_STRIP_PROVIDER, assets.ORIGINAL_PAVEMENT, assets.ORIGINAL_GRASS)
        cls.blocks = loader.load_many(cls.refs)

    def check_root(self, root):
        assets.validate_original_solver(root, cached=True, cached_strips=True, separators=True)

    def test_independent_aggregate_provider_contract(self):
        root = self.blocks[assets.ORIGINAL_SOLVER]
        if not root.get('blosm_aggregate_provider_input'):
            self.skipTest('This revision does not use independent aggregate providers')
        provider = self.loader.load_many((assets.ORIGINAL_AGGREGATE_PROVIDER,))[assets.ORIGINAL_AGGREGATE_PROVIDER]
        self.check_root(root)
        assets.validate_original_aggregate_provider(provider)
        # Ordinary-only entry evaluation must not retain the aggregate solver
        # through a nested group call. Its native Object Info dependency is
        # supplied only when authoring an aggregate network.
        self.assertNotIn(provider.name, self.reachable_groups(root))
        self.assertTrue(any(n.bl_idname == 'GeometryNodeObjectInfo' and
                            n.inputs['Object'].is_linked for n in root.nodes))
        for property_name, invalid in (('blosm_aggregate_provider', 0),
                ('blosm_street_asset_version', 'wrong'),
                ('blosm_street_input_schema', 'original_spans.1'),
                ('st_aggregate_collection_inputs', 0)):
            copy = provider.copy()
            try:
                copy[property_name] = invalid
                with self.assertRaisesRegex(ValueError, 'Incompatible native aggregate provider'):
                    assets.validate_original_aggregate_provider(copy)
            finally:
                bpy.data.node_groups.remove(copy)
        for name in ('Aggregate Provider Key', 'Source Objects', 'Mapping', 'Mesh Resolution'):
            copy = provider.copy()
            try:
                item, = [s for s in copy.interface.items_tree if s.item_type == 'SOCKET'
                         and s.in_out == 'INPUT' and s.name == name]
                copy.interface.remove(item)
                with self.assertRaisesRegex(ValueError, 'Aggregate provider input contract mismatch'):
                    assets.validate_original_aggregate_provider(copy)
            finally:
                bpy.data.node_groups.remove(copy)

    def test_independent_aggregate_provider_requires_matching_authoring_asset(self):
        from blosm.way.street_network.aggregate_authoring import AggregateNetworkBuilder
        root = self.blocks[assets.ORIGINAL_SOLVER]
        if not root.get('blosm_aggregate_provider_input'):
            self.skipTest('This revision does not use independent aggregate providers')
        provider = self.loader.load_many((assets.ORIGINAL_AGGREGATE_PROVIDER,))[assets.ORIGINAL_AGGREGATE_PROVIDER]
        before = {o.as_pointer() for o in bpy.data.objects}
        with self.assertRaisesRegex(ValueError, 'Incompatible native aggregate provider'):
            AggregateNetworkBuilder(root)
        copy = root.copy()
        try:
            copy['blosm_street_asset_version'] = 'wrong'
            with self.assertRaisesRegex(ValueError, 'revisions must match'):
                AggregateNetworkBuilder(copy, provider=provider)
        finally:
            bpy.data.node_groups.remove(copy)
        self.assertEqual(before, {o.as_pointer() for o in bpy.data.objects})

    def reachable_groups(self, root):
        pending, result = [root], {}
        while pending:
            tree = pending.pop()
            if tree.name in result: continue
            result[tree.name] = tree
            pending.extend(n.node_tree for n in tree.nodes if n.type == 'GROUP' and n.node_tree)
        return result

    def check_connected_outputs(self, tree, expected):
        sockets = {s.name: s.socket_type for s in tree.interface.items_tree
                   if s.item_type == 'SOCKET' and s.in_out == 'OUTPUT'}
        self.assertEqual(sockets, expected, tree.name)
        output, = [n for n in tree.nodes if n.type == 'GROUP_OUTPUT' and n.is_active_output]
        for name in expected:
            self.assertTrue(output.inputs[name].is_linked, (tree.name, name))

    def check_separator_dispatch(self, root, suffix):
        reachable = self.reachable_groups(root)
        for policy in ('', ' / Retain valid surfaces'):
            base = 'Prochitecture Streets / Declared separator dispatch v1'+policy+suffix
            # Batching changes evaluation structure, but keeps both failure
            # policies and the complete public dispatch interface.
            candidates=[base]
            if root.get('blosm_longitudinal_separator_chains'):
                token=' / Free chained separators v1' if root.get('blosm_free_chained_separators') else ' / Chained separators v1'
                candidates.append(base.replace(token,token+' / Original terminal owner rows v1'))
            candidates += [name+' / Independent free-fill batching candidate' for name in tuple(candidates)]
            matches = [reachable[name] for name in candidates if name in reachable]
            self.assertEqual(len(matches), 1, base)
            tree = matches[0]
            if root.get('blosm_longitudinal_separator_chains'):
                self.assertEqual(tree.get('st_longitudinal_separator_chains'),1)
            inputs = {s.name: s.socket_type for s in tree.interface.items_tree
                      if s.item_type == 'SOCKET' and s.in_out == 'INPUT'}
            expected = {name: 'NodeSocketGeometry' for name in
                        ('Pavement', 'Corridors', 'External Approaches', 'Mouths', 'Separators')}
            if root.get('blosm_chained_separators'):
                expected['Separator Owners'] = 'NodeSocketGeometry'
            expected.update({'Expected Count': 'NodeSocketInt', 'Cell Size': 'NodeSocketFloat',
                             'Grass': 'NodeSocketMaterial'})
            self.assertEqual(inputs, expected, tree.name)
            self.check_connected_outputs(tree, {'Geometry': 'NodeSocketGeometry',
                'Status': 'NodeSocketGeometry', 'Valid': 'NodeSocketBool', 'Output Valid': 'NodeSocketBool'})

    def test_distributed_library_has_all_shared_roots_and_separator_construction(self):
        self.check_root(self.blocks[assets.ORIGINAL_SOLVER])
        assets.validate_original_patch_provider(self.blocks[assets.ORIGINAL_PATCH_PROVIDER])
        assets.validate_original_strip_provider(self.blocks[assets.ORIGINAL_STRIP_PROVIDER])
        pending = [self.blocks[ref] for ref in self.refs[:3]]; names = set()
        while pending:
            tree = pending.pop()
            if tree.name in names:continue
            names.add(tree.name)
            pending.extend(n.node_tree for n in tree.nodes if n.type == 'GROUP' and n.node_tree)
        self.assertIn('Prochitecture Streets / Original closed separator assembly v3 / Unequal cap counts', names)
        # Diagnostic performance alternatives remain opt-in development code.
        # Packaging their source must not silently enable a failed experiment
        # in the normal shared library.
        for experimental in ('Native registry cache v1', 'Direct ordinary incidence addresses v1',
                             'Independent mouth identities v1', ' / Mouth IDs v1'):
            self.assertFalse(any(experimental in name for name in names), experimental)

    def test_missing_separator_input_is_rejected(self):
        for name in ('Separators', 'Separator Count', 'Grass'):
            with self.subTest(name=name):
                root = self.blocks[assets.ORIGINAL_SOLVER].copy()
                try:
                    socket = next(s for s in root.interface.items_tree
                                  if s.item_type == 'SOCKET' and s.in_out == 'INPUT' and s.name == name)
                    root.interface.remove(socket)
                    with self.assertRaisesRegex(ValueError, name):self.check_root(root)
                finally:bpy.data.node_groups.remove(root)

    def test_free_ended_capability_has_shared_strict_and_partial_dispatch(self):
        root = self.blocks[assets.ORIGINAL_SOLVER]
        enabled = root.get('blosm_free_ended_separators', 0)
        if self.require_free_ended: self.assertEqual(enabled, 1)
        if not enabled: self.skipTest('No free-ended separator capability')
        pending = [root]; names = set()
        while pending:
            tree = pending.pop()
            if tree.name in names: continue
            names.add(tree.name)
            pending.extend(n.node_tree for n in tree.nodes if n.type == 'GROUP' and n.node_tree)
        suffix=' / Chained separators v1' if root.get('blosm_chained_separators') else ''
        if root.get('blosm_free_chained_separators'):suffix += ' / Free chained separators v1'
        if root.get('blosm_crossing_separator_regions'):suffix += ' / Crossing regions v1'
        self.check_separator_dispatch(root, suffix)
        self.assertIn('Prochitecture Streets / Free-ended separator boundary v1', names)

    def test_unknown_free_ended_capability_is_rejected(self):
        root = self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            root['blosm_free_ended_separators'] = 2
            with self.assertRaisesRegex(ValueError, 'free-ended'): self.check_root(root)
        finally: bpy.data.node_groups.remove(root)

    def test_chained_capability_has_shared_strict_and_partial_construction(self):
        root=self.blocks[assets.ORIGINAL_SOLVER]
        if self.require_chained:self.assertEqual(root.get('blosm_chained_separators'),1)
        if not root.get('blosm_chained_separators'):self.skipTest('No chained separator capability')
        pending,names=[root],set()
        while pending:
            tree=pending.pop()
            if tree.name in names:continue
            names.add(tree.name)
            pending.extend(n.node_tree for n in tree.nodes if n.type=='GROUP' and n.node_tree)
        name='Prochitecture Streets / Original closed separator assembly v3 / Unequal cap counts / Chained boundaries'
        parity=root.get('blosm_separator_parity_preparation',0)
        self.assertIn(parity,(0,1))
        suffix=' / Local parity repair v1' if parity else ''
        for branch in ('',' / Retain valid separators v1'):
            expected=name+branch+suffix
            self.assertIn(expected,names)
            self.assertEqual(bpy.data.node_groups[expected].get('st_separator_parity_preparation',0),parity)

    def test_unknown_chained_capability_is_rejected(self):
        root=self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            root['blosm_chained_separators']=2
            with self.assertRaisesRegex(ValueError,'chained'):self.check_root(root)
        finally:bpy.data.node_groups.remove(root)

    def test_free_chained_capability_has_shared_boundary_and_quad_fill(self):
        root=self.blocks[assets.ORIGINAL_SOLVER]
        if self.require_free_chained:self.assertEqual(root.get('blosm_free_chained_separators'),1)
        if not root.get('blosm_free_chained_separators'):self.skipTest('No free chained separator capability')
        pending,names=[root],set()
        while pending:
            tree=pending.pop()
            if tree.name in names:continue
            names.add(tree.name)
            pending.extend(n.node_tree for n in tree.nodes if n.type=='GROUP' and n.node_tree)
        self.assertTrue(any('Free-ended separator boundary' in name and 'Chained' in name for name in names))
        self.assertTrue(any('Free cap' in name for name in names))

    def test_unknown_free_chained_capability_is_rejected(self):
        root=self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            root['blosm_free_chained_separators']=2
            with self.assertRaisesRegex(ValueError,'free-chain'):self.check_root(root)
        finally:bpy.data.node_groups.remove(root)

    def test_longitudinal_capability_reaches_complete_and_partial_ownership(self):
        root=self.blocks[assets.ORIGINAL_SOLVER]
        if self.require_longitudinal_separator_chains:
            self.assertEqual(root.get('blosm_longitudinal_separator_chains'),1)
        if not root.get('blosm_longitudinal_separator_chains'):
            self.skipTest('No multi-Street separator side capability')
        self.check_root(root)
        groups=self.reachable_groups(root)
        policies=set()
        for tree in groups.values():
            calls=[n for n in tree.nodes if n.type=='GROUP' and n.node_tree
                   and n.node_tree.name.startswith('Prochitecture Streets / Declared separator dispatch v1')
                   and ' / Chained separators v1' in n.node_tree.name]
            if not calls:continue
            self.assertEqual(tree.get('st_longitudinal_separator_chains'),1,tree.name)
            registries=[n for n in tree.nodes if n.type=='GROUP' and n.node_tree
                        and n.node_tree.get('st_longitudinal_separator_registry_group')==1]
            self.assertLessEqual(len(registries),1,tree.name)
            registry=registries[0].node_tree if registries else tree
            if registries:
                self.assertEqual(registry.get('st_separator_patch_status_transport'),1,registry.name)
                self.assertTrue(any(n.bl_idname=='GeometryNodePointsToVertices' and
                    n.label=='Read point-cloud patch status as vertices' for n in registry.nodes),registry.name)
                from blosm.way.gn.editable_input.curved_longitudinal_upgrade import normalize_saved_status
                counts=(len(registry.nodes),len(registry.links))
                self.assertFalse(normalize_saved_status(registry))
                self.assertEqual((len(registry.nodes),len(registry.links)),counts)
                self.assertTrue(registries[0].outputs['Valid'].is_linked,tree.name)
                for socket in registries[0].inputs:
                    self.assertTrue(socket.is_linked,(tree.name,socket.name))
                self.check_connected_outputs(registry,{'Valid':'NodeSocketBool'})
            loops=[n for n in registry.nodes if n.bl_idname=='GeometryNodeRepeatOutput'
                   and n.get('st_longitudinal_separator_registry')==1]
            self.assertEqual(len(loops),1,tree.name)
            self.assertTrue(loops[0].outputs['Valid'].is_linked,tree.name)
            for call in calls:
                self.assertEqual(call.node_tree.get('st_longitudinal_separator_chains'),1)
                policies.add(' / Retain valid surfaces' in call.node_tree.name)
        self.assertEqual(policies,{False,True})
        if root.get('blosm_free_chained_separators'):
            name='Prochitecture Streets / Free-ended separator boundary v1 / Chained v1 / Original terminal owner rows v1'
            self.assertIn(name,groups)
            self.assertEqual(groups[name].get('st_original_terminal_owner_rows'),1)
            self.assertNotIn('Prochitecture Streets / Free-ended separator boundary v1 / Chained v1',groups)

    def test_crossing_region_capability_has_shared_boundary_and_graded_fill(self):
        root=self.blocks[assets.ORIGINAL_SOLVER]
        if self.require_crossing_regions:self.assertEqual(root.get('blosm_crossing_separator_regions'),1)
        if not root.get('blosm_crossing_separator_regions'):self.skipTest('No crossing-region capability')
        pending,names=[root],set()
        while pending:
            tree=pending.pop()
            if tree.name in names:continue
            names.add(tree.name)
            pending.extend(n.node_tree for n in tree.nodes if n.type=='GROUP' and n.node_tree)
        self.assertIn('Prochitecture Streets / Chained separator boundary v1',names)
        base = 'Prochitecture Streets / Curved capped separator quads v1 / Unequal cap counts / Graded tail cap'
        reachable = self.reachable_groups(root)
        matches = [reachable[name] for name in (base, base+' / Shared standalone publication v1')
                   if name in reachable]
        self.assertEqual(len(matches), 1, base)
        fill = matches[0]
        self.check_connected_outputs(fill, {'Geometry': 'NodeSocketGeometry',
            'Separator': 'NodeSocketGeometry', 'Diagnostics': 'NodeSocketGeometry', 'Valid': 'NodeSocketBool'})
        if fill.name != base:
            self.assertEqual(fill.get('st_shared_separator_construction'), 1)
            core_name = base.replace(' / Unequal cap counts', ' / Local assembly surface / Unequal cap counts')
            core_name += ' / Shared construction core v1'
            calls = [n for n in fill.nodes if n.type == 'GROUP' and n.node_tree and n.node_tree.name == core_name]
            self.assertEqual(len(calls), 1, core_name)
            self.assertTrue(calls[0].outputs['Separator'].is_linked)
            self.assertTrue(calls[0].outputs['Valid'].is_linked)
        suffix=' / Chained separators v1'
        if root.get('blosm_free_chained_separators'):suffix+=' / Free chained separators v1'
        suffix+=' / Crossing regions v1'
        self.check_separator_dispatch(root, suffix)

    def test_disconnected_separator_publication_is_rejected(self):
        root = self.blocks[assets.ORIGINAL_SOLVER]
        dispatch = next((tree for name, tree in self.reachable_groups(root).items()
                         if name.startswith('Prochitecture Streets / Declared separator dispatch v1')), None)
        if dispatch is None: self.skipTest('No separator dispatch in this revision')
        copy = dispatch.copy()
        try:
            output, = [n for n in copy.nodes if n.type == 'GROUP_OUTPUT' and n.is_active_output]
            for link in tuple(output.inputs['Valid'].links): copy.links.remove(link)
            with self.assertRaises(AssertionError):
                self.check_connected_outputs(copy, {'Geometry': 'NodeSocketGeometry',
                    'Status': 'NodeSocketGeometry', 'Valid': 'NodeSocketBool', 'Output Valid': 'NodeSocketBool'})
        finally: bpy.data.node_groups.remove(copy)

    def test_unknown_crossing_region_capability_is_rejected(self):
        root=self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            root['blosm_crossing_separator_regions']=2
            with self.assertRaisesRegex(ValueError,'crossing-region'):self.check_root(root)
        finally:bpy.data.node_groups.remove(root)

    def test_crossing_region_capability_requires_chained_support(self):
        root=self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            root['blosm_crossing_separator_regions']=1
            root['blosm_free_chained_separators']=0
            # Isolate this dependency from the independent longitudinal one.
            root['blosm_longitudinal_separator_chains']=0
            root['blosm_chained_separators']=0
            with self.assertRaisesRegex(ValueError,'Crossing regions require chained'):self.check_root(root)
        finally:bpy.data.node_groups.remove(root)

    def test_longitudinal_capability_requires_chained_support(self):
        root=self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            root['blosm_longitudinal_separator_chains']=1
            root['blosm_free_chained_separators']=0
            root['blosm_chained_separators']=0
            with self.assertRaisesRegex(ValueError,'Longitudinal separators require chained'):
                self.check_root(root)
        finally:bpy.data.node_groups.remove(root)

    def test_free_chained_capability_requires_chained_support(self):
        root=self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            root['blosm_free_chained_separators']=1;root['blosm_chained_separators']=0
            with self.assertRaisesRegex(ValueError,'chained'):self.check_root(root)
        finally:bpy.data.node_groups.remove(root)

    def test_chained_capability_requires_its_inputs_and_dispatch(self):
        original=self.blocks[assets.ORIGINAL_SOLVER]
        if not original.get('blosm_chained_separators'):self.skipTest('No chained separator capability')
        for key in ('Separator Owners','blosm_free_ended_separators','blosm_local_failures'):
            with self.subTest(key=key):
                root=original.copy()
                try:
                    if key=='Separator Owners':
                        socket=next(s for s in root.interface.items_tree
                            if s.item_type=='SOCKET' and s.in_out=='INPUT' and s.name==key)
                        root.interface.remove(socket)
                    else:root[key]=0
                    with self.assertRaises(ValueError):self.check_root(root)
                finally:bpy.data.node_groups.remove(root)

    def test_boundary_run_capability_has_reachable_native_construction(self):
        root = self.blocks[assets.ORIGINAL_SOLVER]
        if self.require_boundary_runs: self.assertEqual(root.get('blosm_boundary_runs'), 1)
        if not root.get('blosm_boundary_runs'): self.skipTest('No boundary-run capability')
        pending, names = [root], set()
        while pending:
            tree = pending.pop()
            if tree.name in names: continue
            names.add(tree.name)
            pending.extend(n.node_tree for n in tree.nodes if n.type == 'GROUP' and n.node_tree)
        for name in ('Category-bounded run assembly v1', 'Category-bounded profiled run v1',
                     'Category-bounded original span path v1'):
            if root.get('blosm_original_interpolation'):
                name+=(' / Exact or smooth v3 / Mixed seam frames' if name.endswith('span path v1')
                       else ' / Exact or smooth v1')
            if root.get('blosm_selection_boundaries') and name.startswith('Category-bounded profiled run v1'):
                name+=' / Declared selection trims v1'
            if root.get('blosm_interior_selection_cuts'):
                if name.startswith('Category-bounded profiled run v1'):name+=' / Paired selection cuts v1'
                elif name.startswith('Category-bounded original span path v1'):name+=' / Paired selection spans v1'
            self.assertIn('Prochitecture Streets / '+name, names)

    def test_missing_boundary_run_inputs_are_rejected(self):
        original = self.blocks[assets.ORIGINAL_SOLVER]
        if not original.get('blosm_boundary_runs'): self.skipTest('No boundary-run capability')
        for name in ('Boundary Runs', 'Boundary Run Count'):
            root = original.copy()
            try:
                socket = next(s for s in root.interface.items_tree
                              if s.item_type == 'SOCKET' and s.in_out == 'INPUT' and s.name == name)
                root.interface.remove(socket)
                with self.assertRaisesRegex(ValueError, name): self.check_root(root)
            finally: bpy.data.node_groups.remove(root)

    def test_boundary_run_capability_requires_its_prerequisites(self):
        original = self.blocks[assets.ORIGINAL_SOLVER]
        if not original.get('blosm_boundary_runs'): self.skipTest('No boundary-run capability')
        for name in ('blosm_section_scopes', 'blosm_standalone_streets', 'blosm_local_failures'):
            root = original.copy()
            try:
                root[name] = 0
                with self.assertRaises(ValueError): self.check_root(root)
            finally: bpy.data.node_groups.remove(root)

    def test_unknown_boundary_run_revision_is_rejected(self):
        root = self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            root['blosm_boundary_runs'] = 2
            with self.assertRaisesRegex(ValueError, 'category-bounded'): self.check_root(root)
        finally: bpy.data.node_groups.remove(root)

    def test_wrong_separator_socket_type_is_rejected(self):
        root = self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            socket = next(s for s in root.interface.items_tree
                          if s.item_type == 'SOCKET' and s.in_out == 'INPUT' and s.name == 'Separator Count')
            root.interface.remove(socket)
            root.interface.new_socket(name='Separator Count', in_out='INPUT', socket_type='NodeSocketFloat')
            with self.assertRaisesRegex(ValueError, 'Separator Count'):self.check_root(root)
        finally:bpy.data.node_groups.remove(root)

    def test_wrong_revision_is_rejected(self):
        root = self.blocks[assets.ORIGINAL_SOLVER].copy()
        try:
            root['blosm_street_asset_version'] = 'incompatible-test-revision'
            with self.assertRaisesRegex(ValueError, 'Incompatible'):self.check_root(root)
        finally:bpy.data.node_groups.remove(root)

    def local_failure_blocks(self):
        root=self.blocks[assets.ORIGINAL_SOLVER]
        if self.require_local_failures:self.assertEqual(root.get('blosm_local_failures'),1)
        if not root.get('blosm_local_failures'):self.skipTest('Bundled revision has no local failure capability')
        return root,self.blocks[assets.ORIGINAL_STRIP_PROVIDER]

    def test_local_failure_display_is_loaded_through_the_asset_service(self):
        root,strip=self.local_failure_blocks()
        assets.validate_original_failure_support(root,strip)
        loaded=self.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,assets.ORIGINAL_DIAGNOSTICS))
        assets.validate_original_failure_display(loaded[assets.ORIGINAL_FAILURE_DISPLAY])
        self.assertEqual(loaded[assets.ORIGINAL_DIAGNOSTICS].bl_rna.identifier,'Material')

    def test_incomplete_failure_transport_is_rejected(self):
        root,strip=self.local_failure_blocks();copy=strip.copy()
        try:
            del copy['blosm_retains_local_failures']
            with self.assertRaises(ValueError):assets.validate_original_failure_support(root,copy)
        finally:bpy.data.node_groups.remove(copy)

    def test_missing_partial_status_output_is_rejected(self):
        original,strip=self.local_failure_blocks();root=original.copy()
        try:
            socket=next(s for s in root.interface.items_tree if s.item_type=='SOCKET'
                        and s.in_out=='OUTPUT' and s.name=='Output Valid')
            root.interface.remove(socket)
            with self.assertRaises(ValueError):assets.validate_original_failure_support(root,strip)
        finally:bpy.data.node_groups.remove(root)

    def test_standalone_capability_requires_its_ownership_input(self):
        original=self.blocks[assets.ORIGINAL_SOLVER]
        if self.require_standalone:self.assertEqual(original.get('blosm_standalone_streets'),1)
        if not original.get('blosm_standalone_streets'):self.skipTest('No standalone capability declared')
        for name in ('Standalone Streets','Standalone Count'):
            with self.subTest(name=name):
                root=original.copy()
                try:
                    socket=next(s for s in root.interface.items_tree if s.item_type=='SOCKET'
                                and s.in_out=='INPUT' and s.name==name)
                    root.interface.remove(socket)
                    with self.assertRaisesRegex(ValueError,name):self.check_root(root)
                finally:bpy.data.node_groups.remove(root)

    def test_filtered_ports_use_the_junction_identity_registry(self):
        root=self.blocks[assets.ORIGINAL_SOLVER]
        if self.require_filtered_ports:self.assertEqual(root.get('blosm_filtered_physical_ports'),1)
        if not root.get('blosm_filtered_physical_ports'):self.skipTest('No filtered physical port capability')
        pending=[root];names=set()
        while pending:
            tree=pending.pop()
            if tree.name in names:continue
            names.add(tree.name)
            pending.extend(n.node_tree for n in tree.nodes if n.type=='GROUP' and n.node_tree)
        self.assertTrue(any(name.startswith('Prochitecture Streets / Compound original incidence registry v3') for name in names))

    def test_section_scopes_use_checked_adjacent_section_handoffs(self):
        root=self.blocks[assets.ORIGINAL_SOLVER]
        if self.require_section_scopes:self.assertEqual(root.get('blosm_section_scopes'),1)
        if not root.get('blosm_section_scopes'):self.skipTest('No Section scope capability')
        # The cached assembly consumes Strip Results. Adjacent-Section checks
        # execute in the shared strip provider; an unused copy in the assembly
        # is not evidence of this capability and may legitimately be removed.
        strip=self.blocks[assets.ORIGINAL_STRIP_PROVIDER]
        assets.validate_original_strip_provider(strip)
        pending=[strip];names=set()
        while pending:
            tree=pending.pop()
            if tree.name in names:continue
            names.add(tree.name)
            predecessors={node:[] for node in tree.nodes}
            for link in tree.links:
                if link.is_valid and not link.is_muted:
                    predecessors[link.to_node].append(link.from_node)
            for node in tree.nodes:
                paired=getattr(node,'paired_output',None)
                if paired is not None:
                    predecessors[node].append(paired);predecessors[paired].append(node)
            live=set();todo=[node for node in tree.nodes if node.type=='GROUP_OUTPUT']
            while todo:
                node=todo.pop()
                if node in live:continue
                live.add(node);todo.extend(predecessors[node])
            pending.extend(n.node_tree for n in live if n.type=='GROUP' and n.node_tree)
        self.assertTrue(any('Checked Section handoffs' in name or 'Checked feature handoffs' in name for name in names))
        if root.get('blosm_original_interpolation'):
            suffix=' / Authored through Connectors v1' if root.get('blosm_authored_through_connectors') else ''
            self.check_original_port_reader(names,' / Exact or smooth v3 / Mixed seam frames / Declared dimensions v1'+suffix)
        else:
            self.assertTrue(names.intersection({
                'Prochitecture Streets / Original Section-end port path v2 / Feature curves',
                'Prochitecture Streets / Original Section-end port path v3 / Feature curves and seams',
                'Prochitecture Streets / Original Section-end port path v4 / Profile admission',
            }))

    def test_internal_end_capabilities_must_agree_before_authoring(self):
        root=self.blocks[assets.ORIGINAL_SOLVER].copy()
        patch=self.blocks[assets.ORIGINAL_PATCH_PROVIDER].copy()
        try:
            for flags in ((0,0),(1,1),(1,0),(0,1),(2,2)):
                root['blosm_internal_logical_ends'],patch['blosm_internal_logical_ends']=flags
                with self.subTest(flags=flags):
                    if flags in ((0,0),(1,1)):
                        self.assertEqual(assets.validate_original_internal_ends(root,patch),flags==(1,1))
                    else:
                        with self.assertRaises(ValueError):assets.validate_original_internal_ends(root,patch)
        finally:
            bpy.data.node_groups.remove(root);bpy.data.node_groups.remove(patch)

    def test_internal_end_capability_contains_the_native_footprint_guard(self):
        root=self.blocks[assets.ORIGINAL_SOLVER];patch=self.blocks[assets.ORIGINAL_PATCH_PROVIDER]
        enabled=assets.validate_original_internal_ends(root,patch)
        if self.require_internal_ends:self.assertTrue(enabled)
        if not enabled:self.skipTest('No retained internal logical endpoint capability')
        suffix=' / Exact profile frames v1' if root.get('blosm_original_interpolation') else ''
        name='Prochitecture Streets / Original logical cluster boundary / Retained internal ends v2'+suffix
        # Pattern recipes copy boundary callers while retaining their native
        # guards. A Blender copy suffix is not a missing footprint check, and
        # an unused group with the expected name is not proof of an active one.
        copies=[tree for tree in self.reachable_groups(patch).values()
                if re.fullmatch(re.escape(name)+r'(?:\.\d{3,})?',tree.name)]
        self.assertTrue(copies,'Missing retained internal-footprint constructor')
        def connected_guard(tree):
            inputs=[node for node in tree.nodes if node.bl_idname=='GeometryNodeRepeatInput'
                    and node.label=='Preserve complete internal Street footprints']
            if len(inputs)!=1:return False
            guard=inputs[0].paired_output
            if guard is None:return False
            output=next(node for node in tree.nodes if node.type=='GROUP_OUTPUT' and node.is_active_output)
            for key in ('Valid','Boundary'):
                pending=[link.from_node for link in output.inputs[key].links if link.is_valid and not link.is_muted]
                seen=set()
                while pending:
                    node=pending.pop()
                    if node in seen:continue
                    seen.add(node)
                    pending.extend(link.from_node for socket in node.inputs for link in socket.links
                                   if link.is_valid and not link.is_muted)
                if guard not in seen:return False
            return True
        for tree in copies:
            with self.subTest(group=tree.name):
                self.assertTrue(connected_guard(tree),'Footprint guard must control validity and published boundary')
                self.assertTrue(any(node.bl_idname=='GeometryNodeInputNamedAttribute'
                                    and node.data_type=='INT' and node.inputs['Name'].default_value=='st_lc_footprint'
                                    for node in tree.nodes))
                broken=tree.copy()
                try:
                    guard=next(node.paired_output for node in broken.nodes
                               if node.bl_idname=='GeometryNodeRepeatInput'
                               and node.label=='Preserve complete internal Street footprints')
                    for link in list(guard.outputs['Valid'].links):broken.links.remove(link)
                    self.assertFalse(connected_guard(broken),'Disconnected guard must not qualify the library')
                finally:bpy.data.node_groups.remove(broken)

    def test_shared_owner_capabilities_and_program_inputs_must_agree(self):
        root=self.blocks[assets.ORIGINAL_SOLVER];patch=self.blocks[assets.ORIGINAL_PATCH_PROVIDER]
        enabled=assets.validate_original_shared_owners(root,patch)
        if self.require_shared_owners:self.assertTrue(enabled)
        if not enabled:self.skipTest('No shared logical owner capability')
        for target in (root,patch):
            copy=target.copy()
            try:
                for flag in (0,2):
                    copy['blosm_shared_logical_owners']=flag
                    with self.assertRaises(ValueError):assets.validate_original_shared_owners(
                        copy if target==root else root,copy if target==patch else patch)
                copy['blosm_shared_logical_owners']=1
                item=next(s for s in copy.interface.items_tree if s.item_type=='SOCKET'
                          and s.in_out=='INPUT' and s.name=='Chain Steps')
                copy.interface.remove(item)
                with self.assertRaisesRegex(ValueError,'boundary program'):
                    assets.validate_original_shared_owners(copy if target==root else root,copy if target==patch else patch)
            finally:bpy.data.node_groups.remove(copy)

    def test_original_seam_capabilities_must_agree_on_all_three_providers(self):
        from itertools import product
        original=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER))
        enabled=assets.validate_original_seams(*original)
        if self.require_seams:self.assertTrue(enabled)
        copies=tuple(g.copy() for g in original)
        try:
            for flags in product((0,1),repeat=3):
                for g,flag in zip(copies,flags):
                    g['blosm_original_seams']=flag;g['blosm_finite_seam_transitions']=flag
                if len(set(flags))==1:self.assertEqual(assets.validate_original_seams(*copies),flags[0]==1)
                else:
                    with self.assertRaisesRegex(ValueError,'seam'):assets.validate_original_seams(*copies)
            for g in copies:g['blosm_original_seams']=1;g['blosm_finite_seam_transitions']=1
            for g in copies:
                g['blosm_finite_seam_transitions']=0
                with self.assertRaisesRegex(ValueError,'finite-profile'):assets.validate_original_seams(*copies)
                g['blosm_finite_seam_transitions']=1
        finally:
            for g in copies:bpy.data.node_groups.remove(g)

    def test_original_seam_provider_contains_both_qualified_branches(self):
        blocks=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER))
        if not assets.validate_original_seams(*blocks):self.skipTest('No original seam capability')
        pending=[blocks[1]];names=set()
        while pending:
            tree=pending.pop()
            if tree.name in names:continue
            names.add(tree.name);pending.extend(n.node_tree for n in tree.nodes if n.type=='GROUP' and n.node_tree)
        for name in ('Matching original seam mouths v1','Finite original seam transition v1','Existing cubic interval envelope v1'):
            if blocks[0].get('blosm_original_interpolation'):
                name+=(' / Exact frames v1' if name=='Existing cubic interval envelope v1' else ' / Exact and mixed frames v2')
            self.assertIn('Prochitecture Streets / '+name,names)
        suffix=' / Exact or smooth v3 / Mixed seam frames / Declared dimensions v1' if blocks[0].get('blosm_original_interpolation') else ''
        if blocks[0].get('blosm_authored_through_connectors'):suffix+=' / Authored through Connectors v1'
        self.check_original_port_reader(names,suffix)

    def test_category_handoff_capabilities_and_native_branch(self):
        from itertools import product
        original=tuple(self.blocks[k] for k in (assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,assets.ORIGINAL_STRIP_PROVIDER))
        enabled=assets.validate_original_category_connections(*original)
        if self.require_category_connections:self.assertTrue(enabled)
        copies=tuple(g.copy() for g in original)
        try:
            for g in copies:
                g['blosm_section_scopes']=1;g['blosm_original_seams']=1;g['blosm_finite_seam_transitions']=1
            for flags in product((0,1),repeat=3):
                for g,flag in zip(copies,flags):g['blosm_category_connections']=flag
                if len(set(flags))==1:self.assertEqual(assets.validate_original_category_connections(*copies),flags[0]==1)
                else:
                    with self.assertRaisesRegex(ValueError,'category'):assets.validate_original_category_connections(*copies)
            copies[0]['blosm_section_scopes']=0
            with self.assertRaisesRegex(ValueError,'Section scopes'):assets.validate_original_category_connections(*copies)
        finally:
            for g in copies:bpy.data.node_groups.remove(g)
        if enabled:
            pending=[original[2]];names=set()
            while pending:
                tree=pending.pop()
                if tree.name in names:continue
                names.add(tree.name);pending.extend(n.node_tree for n in tree.nodes if n.type=='GROUP' and n.node_tree)
            self.assertTrue(any('Original Section spans v3 / Checked feature handoffs' in name for name in names))

    def test_stable_quad_energy_capability_matches_reachable_construction(self):
        root=self.blocks[assets.ORIGINAL_SOLVER];patch=self.blocks[assets.ORIGINAL_PATCH_PROVIDER]
        flags=tuple(g.get('blosm_stable_quad_energy',0) for g in (root,patch))
        self.assertIn(flags,((0,0),(2,2)))
        self.assertEqual(assets.validate_original_quad_energy(root,patch),flags==(2,2))
        copies=(root.copy(),patch.copy())
        try:
            for malformed in ((2,0),(0,2),(2,1),(1,1)):
                for tree,flag in zip(copies,malformed):tree['blosm_stable_quad_energy']=flag
                with self.assertRaisesRegex(ValueError,'stable quad energy'):
                    assets.validate_original_quad_energy(*copies)
        finally:
            for tree in copies:bpy.data.node_groups.remove(tree)
        if self.require_stable_quad_energy:self.assertEqual(flags,(2,2))
        if flags==(0,0):self.skipTest('No stable-energy capability in this revision')
        names=set();pending=[root,patch]
        while pending:
            tree=pending.pop()
            if tree.name in names:continue
            names.add(tree.name)
            pending.extend(n.node_tree for n in tree.nodes if n.type=='GROUP' and n.node_tree)
        relaxation=[name for name in names if name.startswith('Prochitecture Streets / Pinned curved quad descent v1')]
        self.assertTrue(relaxation)
        self.assertTrue(all('Stable energy change v2' in name for name in relaxation),relaxation)
        self.assertIn('Prochitecture Streets / Quad corner stable energy change v2',names)

    def test_shared_owner_asset_contains_final_mesh_footprint_guard(self):
        root=self.blocks[assets.ORIGINAL_SOLVER];patch=self.blocks[assets.ORIGINAL_PATCH_PROVIDER]
        if not assets.validate_original_shared_owners(root,patch):self.skipTest('No shared logical owner capability')
        pending=[patch];names=set()
        while pending:
            tree=pending.pop()
            if tree.name in names:continue
            names.add(tree.name)
            pending.extend(n.node_tree for n in tree.nodes if n.type=='GROUP' and n.node_tree)
        self.assertTrue(any('Original chain profiled patch' in name for name in names))
        suffix=' / Exact or smooth paths v1' if root.get('blosm_original_interpolation') else ''
        self.assertIn('Prochitecture Streets / Exact generated chain mesh boundary v1'+suffix,names)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--link-library', action='store_true',
                        help='Run contracts through the universal loader with linked asset datablocks')
    parser.add_argument('--library', type=Path, help='Optional library to qualify before publication')
    parser.add_argument('--revision', help='Expected revision of that candidate library')
    parser.add_argument('--require-local-failures',action='store_true')
    parser.add_argument('--require-standalone',action='store_true')
    parser.add_argument('--require-filtered-ports',action='store_true')
    parser.add_argument('--require-section-scopes',action='store_true')
    parser.add_argument('--require-internal-ends',action='store_true')
    parser.add_argument('--require-shared-owners',action='store_true')
    parser.add_argument('--require-seams',action='store_true')
    parser.add_argument('--require-category-connections',action='store_true')
    parser.add_argument('--require-stable-quad-energy',action='store_true')
    parser.add_argument('--require-free-ended',action='store_true')
    parser.add_argument('--require-boundary-runs',action='store_true')
    parser.add_argument('--require-chained',action='store_true')
    parser.add_argument('--require-free-chained',action='store_true')
    parser.add_argument('--require-crossing-regions',action='store_true')
    parser.add_argument('--require-continuation-regions',action='store_true')
    parser.add_argument('--require-internal-regions',action='store_true')
    parser.add_argument('--require-assembly-payload-cleanup',action='store_true')
    parser.add_argument('--require-interpolation',action='store_true')
    parser.add_argument('--require-cycle-cores',action='store_true')
    parser.add_argument('--require-shared-kernels',action='store_true')
    parser.add_argument('--require-strip-refinement',action='store_true')
    parser.add_argument('--require-guide-correction',action='store_true')
    parser.add_argument('--require-guarded-profile',action='store_true')
    parser.add_argument('--require-core-homotopy',action='store_true')
    parser.add_argument('--require-separator-parity',action='store_true')
    parser.add_argument('--require-scoped-separator-parity',action='store_true')
    parser.add_argument('--require-checked-chain-contacts',action='store_true')
    parser.add_argument('--require-selection-boundaries',action='store_true')
    parser.add_argument('--require-authored-through-connectors',action='store_true')
    parser.add_argument('--require-authored-dimension-runs',action='store_true')
    parser.add_argument('--require-authored-dimension-approaches',action='store_true')
    parser.add_argument('--require-interior-selection-cuts',action='store_true')
    parser.add_argument('--require-rejected-profiles',action='store_true')
    parser.add_argument('--require-longitudinal-separator-chains',action='store_true')
    parser.add_argument('--tests',nargs='+',choices=sorted(name for name in dir(OriginalSpanAssetTests)
                                                       if name.startswith('test_')),
                        help='Rerun named contracts; omission runs the complete suite')
    args = parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    OriginalSpanAssetTests.linked_library = args.link_library
    if bool(args.library) != bool(args.revision):
        parser.error('--library and --revision must be provided together')
    if args.library:
        assets.ORIGINAL_LIBRARY = replace(assets.ORIGINAL_LIBRARY, source='external',
            path=str(args.library.resolve()), revision=args.revision)
    # A shipped library must have all current capabilities. Optional external
    # revision checks can still exercise an older contract explicitly.
    OriginalSpanAssetTests.require_local_failures=args.require_local_failures or not args.library
    OriginalSpanAssetTests.require_standalone=args.require_standalone or not args.library
    OriginalSpanAssetTests.require_filtered_ports=args.require_filtered_ports or not args.library
    OriginalSpanAssetTests.require_section_scopes=args.require_section_scopes or not args.library
    OriginalSpanAssetTests.require_internal_ends=args.require_internal_ends or not args.library
    OriginalSpanAssetTests.require_shared_owners=args.require_shared_owners or not args.library
    OriginalSpanAssetTests.require_seams=args.require_seams or not args.library
    OriginalSpanAssetTests.require_category_connections=args.require_category_connections or not args.library
    OriginalSpanAssetTests.require_stable_quad_energy=args.require_stable_quad_energy or not args.library
    OriginalSpanAssetTests.require_free_ended=args.require_free_ended or not args.library
    OriginalSpanAssetTests.require_boundary_runs=args.require_boundary_runs or not args.library
    OriginalSpanAssetTests.require_chained=args.require_chained
    OriginalSpanAssetTests.require_free_chained=args.require_free_chained
    OriginalSpanAssetTests.require_crossing_regions=args.require_crossing_regions
    OriginalSpanAssetTests.require_continuation_regions=args.require_continuation_regions
    OriginalSpanAssetTests.require_internal_regions=args.require_internal_regions
    OriginalSpanAssetTests.require_assembly_payload_cleanup=args.require_assembly_payload_cleanup
    OriginalSpanAssetTests.require_interpolation=args.require_interpolation
    OriginalSpanAssetTests.require_cycle_cores=args.require_cycle_cores
    OriginalSpanAssetTests.require_shared_kernels=args.require_shared_kernels
    OriginalSpanAssetTests.require_strip_refinement=args.require_strip_refinement
    OriginalSpanAssetTests.require_guide_correction=args.require_guide_correction
    OriginalSpanAssetTests.require_guarded_profile=args.require_guarded_profile
    OriginalSpanAssetTests.require_core_homotopy=args.require_core_homotopy
    OriginalSpanAssetTests.require_separator_parity=args.require_separator_parity or args.require_scoped_separator_parity
    OriginalSpanAssetTests.require_scoped_separator_parity=args.require_scoped_separator_parity
    OriginalSpanAssetTests.require_checked_chain_contacts=args.require_checked_chain_contacts
    OriginalSpanAssetTests.require_selection_boundaries=args.require_selection_boundaries
    OriginalSpanAssetTests.require_authored_through_connectors=args.require_authored_through_connectors
    OriginalSpanAssetTests.require_authored_dimension_runs=args.require_authored_dimension_runs
    OriginalSpanAssetTests.require_authored_dimension_approaches=args.require_authored_dimension_approaches
    OriginalSpanAssetTests.require_interior_selection_cuts=args.require_interior_selection_cuts
    OriginalSpanAssetTests.require_rejected_profiles=args.require_rejected_profiles
    OriginalSpanAssetTests.require_longitudinal_separator_chains=args.require_longitudinal_separator_chains
    suite = (unittest.TestSuite(OriginalSpanAssetTests(name) for name in args.tests) if args.tests else
             unittest.defaultTestLoader.loadTestsFromTestCase(OriginalSpanAssetTests))
    if not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful():
        raise RuntimeError('Original-span asset contract regression failed')
