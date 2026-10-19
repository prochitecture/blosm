"""Bind shared optional complete/partial separator fills to native dependencies.

Both objects use the same provider program. Their inputs reference the original
owner outputs and metadata, and numeric controls follow native Blender drivers.
No callback computes geometry or responds to source vertex edits.
"""
import bpy
from ..gn.kernel.pipeline import set_input


class SeparatorFillProvider:
    result_inputs=('Complete Separator Fill','Partial Separator Fill')
    key_input='Separator Fill Key'
    policy_input='Separator Fill Policy'
    marker='blosm_separator_fill_provider'

    def __init__(self,solver):
        from .assets import validate_original_separator_fill_provider
        validate_original_separator_fill_provider(solver)
        self.solver=solver

    def install(self,controller,collection):
        from ..gn.editable_input.curved_patch_provider_bindings import input_values,set_inputs
        from .free_bundle_provider import FreeBundleProvider
        network=controller.modifiers[0].node_group
        if network.get('blosm_separator_fill_provider_input')!=1:
            raise ValueError('The shared solver has no separator fill provider inputs')
        if self.solver.get('blosm_street_asset_version')!=network.get('blosm_street_asset_version'):
            raise ValueError('Separator fill provider and shared network revisions must match')
        inputs=input_values(controller)
        for name in (*self.result_inputs,self.key_input):
            if name not in inputs:raise ValueError('Missing shared separator fill input: '+name)
        rows=inputs['Separators']
        if (inputs['Separator Count']<=0 or rows is None or rows.type!='MESH' or
                len(rows.data.vertices)!=inputs['Separator Count']):
            raise ValueError('Only consistent nonempty separator declarations need fill providers')
        if (any(inputs[name] is not None for name in self.result_inputs) or
                any(o.get(self.marker) and o.parent==controller for o in bpy.data.objects)):
            raise ValueError('Rebind or remove the existing separator fill providers first')
        required={s.name for s in self.solver.interface.items_tree if s.item_type=='SOCKET'
                  and s.in_out=='INPUT' and s.name!=self.policy_input}
        if not required<=inputs.keys():
            raise ValueError('Missing original fill controls: '+', '.join(sorted(required-inputs.keys())))
        used={o.get('blosm_separator_fill_key') for o in bpy.data.objects}
        key=1
        while key in used:key+=1
        if key>=2**24:raise ValueError('Separator fill provider identity capacity exceeded')
        children=[];meshes=[];old_key=controller.get('blosm_separator_fill_key')
        try:
            for policy,name in enumerate(self.result_inputs):
                data=bpy.data.meshes.new(controller.name+' / '+name+' host');meshes.append(data)
                child=bpy.data.objects.new(controller.name+' / '+name,data);children.append(child)
                collection.objects.link(child);child.parent=controller
                child.hide_render=True;child.hide_select=True;child.hide_set(True)
                child[self.marker]=1;child['blosm_separator_fill_policy']=policy
                modifier=child.modifiers.new('Shared native separator fill','NODES')
                modifier.node_group=self.solver;modifier.show_viewport=False
                set_inputs(child,dict(inputs,**{self.key_input:key,self.policy_input:policy}))
                FreeBundleProvider._drivers(controller,child,exclude=(self.policy_input,))
            controller['blosm_separator_fill_key']=key
            set_input(controller,self.key_input,key)
            for name,child in zip(self.result_inputs,children):set_input(controller,name,child)
            for child in children:child.modifiers[0].show_viewport=True
            return tuple(children)
        except BaseException:
            for name in (*self.result_inputs,self.key_input):set_input(controller,name,inputs[name])
            if old_key is None:
                if 'blosm_separator_fill_key' in controller:del controller['blosm_separator_fill_key']
            else:controller['blosm_separator_fill_key']=old_key
            for child in children:bpy.data.objects.remove(child,do_unlink=True)
            for data in meshes:
                if data.users==0:bpy.data.meshes.remove(data)
            raise
