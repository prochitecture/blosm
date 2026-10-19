"""Small native asset-preflight tests; surface preservation has separate fixtures.

Run in factory-startup Blender. No complete street library is loaded here.
"""
from itertools import product
from pathlib import Path
import sys
import unittest

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from blosm.way.street_network import assets
from blosm.way.gn.editable_input.curved_port_nodes import protect_saved_profiles
sys.path.insert(0,str(Path(__file__).resolve().parent))
from test_original_span_assets import OriginalSpanAssetTests

PREFIX = 'Prochitecture Streets / Original Section-end port path '


class ProfileAssetTests(unittest.TestCase):
    def setUp(self):
        self.existing = set(bpy.data.node_groups)

    def tearDown(self):
        for tree in set(bpy.data.node_groups)-self.existing:
            bpy.data.node_groups.remove(tree)

    def tree(self, name):
        return bpy.data.node_groups.new(name, 'GeometryNodeTree')

    def roots(self):
        roots = [self.tree(name) for name in ('Network', 'Patch', 'Strip')]
        roots[0]['blosm_local_failures'] = 1
        roots[2]['blosm_retains_local_failures'] = 1
        for name, kind in (('Valid', 'NodeSocketBool'), ('Output Valid', 'NodeSocketBool'),
                           ('Failure Records', 'NodeSocketGeometry')):
            roots[0].interface.new_socket(name=name, in_out='OUTPUT', socket_type=kind)
        return roots

    def reader(self, name='v3 / Feature curves and seams'):
        tree = self.tree(PREFIX+name)
        for direction, key, kind in (('INPUT','Ports','NodeSocketGeometry'),
                                    ('INPUT','Index','NodeSocketInt'),
                                    ('OUTPUT','Geometry','NodeSocketGeometry'),
                                    ('OUTPUT','Valid','NodeSocketBool')):
            tree.interface.new_socket(name=key, in_out=direction, socket_type=kind)
        source = tree.nodes.new('NodeGroupInput')
        output = tree.nodes.new('NodeGroupOutput')
        valid = tree.nodes.new('FunctionNodeInputBool')
        valid.boolean = True
        tree.links.new(source.outputs['Ports'], output.inputs['Geometry'])
        tree.links.new(valid.outputs[0], output.inputs['Valid'])
        return tree

    @staticmethod
    def snapshot(tree):
        return (tuple((n.as_pointer(), n.bl_idname) for n in tree.nodes),
                tuple((l.from_socket.as_pointer(), l.to_socket.as_pointer()) for l in tree.links),
                tuple((s.identifier, s.name, s.in_out, s.socket_type)
                      for s in tree.interface.items_tree if s.item_type=='SOCKET'),
                dict(tree.items()))

    def test_capability_requires_all_three_consumers(self):
        roots = self.roots()
        self.assertFalse(assets.validate_original_rejected_profiles(*roots))
        self.assertFalse(assets.validate_original_rejected_profiles(roots[0],None,None))
        for flags in product((0,1), repeat=3):
            with self.subTest(flags=flags):
                for tree, flag in zip(roots, flags):tree['blosm_rejected_profiles']=flag
                if len(set(flags))==1:
                    self.assertEqual(assets.validate_original_rejected_profiles(*roots), bool(flags[0]))
                else:
                    with self.assertRaisesRegex(ValueError, 'rejected-profile'):
                        assets.validate_original_rejected_profiles(*roots)
        for at in range(3):
            for tree in roots:tree['blosm_rejected_profiles']=1
            roots[at]['blosm_rejected_profiles']=2
            with self.assertRaisesRegex(ValueError, 'rejected-profile'):
                assets.validate_original_rejected_profiles(*roots)

    def test_missing_partial_contract_is_rejected(self):
        roots = self.roots()
        for tree in roots:tree['blosm_rejected_profiles']=1
        for at, key in ((0,'blosm_local_failures'),(2,'blosm_retains_local_failures')):
            roots[at][key]=0
            with self.assertRaisesRegex(ValueError, 'local failure'):
                assets.validate_original_rejected_profiles(*roots)
            roots[at][key]=1
        root=roots[0]
        root.interface.remove(next(s for s in root.interface.items_tree if s.name=='Failure Records'))
        with self.assertRaisesRegex(ValueError, 'partial output'):
            assets.validate_original_rejected_profiles(*roots)

    def test_saved_reader_upgrade_preserves_schema_and_original_nodes(self):
        reader=self.reader();before=self.snapshot(reader)
        rows=protect_saved_profiles([reader])
        self.assertEqual(rows,[dict(group=reader.name,changed=True)])
        after=self.snapshot(reader)
        self.assertEqual(before[2],after[2])
        self.assertLessEqual(set(before[0]),set(after[0]))
        self.assertEqual(reader.get('st_saved_profile_admission'),1)
        self.assertEqual(protect_saved_profiles([reader]),[dict(group=reader.name,changed=False)])
        self.assertEqual(self.snapshot(reader),after)

    def test_reachable_reader_is_upgraded_once(self):
        root=self.tree('Shared assembly');reader=self.reader()
        for unused in range(2):
            root.nodes.new('GeometryNodeGroup').node_tree=reader
        self.assertEqual(len(protect_saved_profiles([root,reader])),1)

    def test_all_readers_are_preflighted_before_mutation(self):
        good=self.reader()
        unknown=self.reader('v99 / Unknown')
        before=self.snapshot(good)
        with self.assertRaisesRegex(ValueError, 'Unknown saved port-reader'):
            protect_saved_profiles([good,unknown])
        self.assertEqual(self.snapshot(good),before)
        unknown.name=PREFIX+'v3 / Z broken output'
        output=next(n for n in unknown.nodes if n.type=='GROUP_OUTPUT')
        unknown.links.remove(output.inputs['Valid'].links[0])
        with self.assertRaisesRegex(ValueError, 'linked original'):
            protect_saved_profiles([good,unknown])
        self.assertEqual(self.snapshot(good),before)

    def test_unknown_admission_revision_is_not_overwritten(self):
        reader=self.reader();reader['st_port_profile_admission']=2
        before=self.snapshot(reader)
        with self.assertRaisesRegex(ValueError, 'profile-admission revision'):
            protect_saved_profiles([reader])
        self.assertEqual(self.snapshot(reader),before)

    def test_feature_producer_is_not_a_port_consumer(self):
        reader=self.reader()
        producer=self.reader('v4 / Profile admission / Shared feature producer v1')
        before=self.snapshot(producer)
        self.assertEqual(len(protect_saved_profiles([producer,reader])),1)
        self.assertEqual(self.snapshot(producer),before)
        with self.assertRaisesRegex(ValueError, 'No original Section-end reader'):
            protect_saved_profiles([producer])

    def test_release_check_rejects_a_bypassed_admission_guard(self):
        reader=self.reader()
        checks=OriginalSpanAssetTests('test_rejected_profiles_require_consistent_native_consumers')
        protect_saved_profiles([reader])
        checks.check_profile_reader_admission(reader)
        output=next(n for n in reader.nodes if n.type=='GROUP_OUTPUT')
        source=next(n for n in reader.nodes if n.type=='GROUP_INPUT')
        valid=next(n for n in reader.nodes if n.bl_idname=='FunctionNodeInputBool')
        for name,old in (('Geometry',source.outputs['Ports']),('Valid',valid.outputs[0])):
            admitted=output.inputs[name].links[0].from_socket
            reader.links.new(old,output.inputs[name])
            with self.assertRaises(AssertionError):checks.check_profile_reader_admission(reader)
            reader.links.new(admitted,output.inputs[name])
        checks.check_profile_reader_admission(reader)


if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(ProfileAssetTests)
    if not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful():
        raise RuntimeError('Profile asset preflight regression')
