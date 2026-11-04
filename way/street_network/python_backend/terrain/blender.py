"""Original-terrain sampling and explicit Blender publication adapter."""
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import time
import bpy
import numpy as np
from mathutils import Vector
from mathutils.bvhtree import BVHTree
from .alignment import PlacementSettings, fit, quad_samples
from .corridor import prepare as prepare_corridor
from ..model import Key, fingerprint
from ..mesh import Mesh, publication_coordinates
from ..providers import Part, AnchorTransform

RAW_Z = 'st_terrain_prepared_z'
MODE_ITEMS = (('PROFILE','Preserve profiles (Python)','Fit complete cross-sections and shared junction minima when Update is pressed'),
              ('DIRECT','Direct projection (GN)','Project each vertex and restore its relative height; terrain cross-slopes affect the profile'))


def require_original_terrain(terrain):
    """Reject explicit modifier/GN dependencies on street output or sources.

    Independent procedural terrain is allowed. Drivers and arbitrary extension
    code still need their own dependency policy; this inspects data references.
    """
    pending, seen = [terrain], set()
    while pending:
        value = pending.pop()
        if value is None or value.as_pointer() in seen:
            continue
        seen.add(value.as_pointer())
        if isinstance(value, bpy.types.Collection):
            pending.extend(value.all_objects)
        elif isinstance(value, bpy.types.Object):
            if value.get('st_py_controller') or value.get('st_py_source_data'):
                raise ValueError('Original terrain must not depend on generated streets or their source objects')
            for modifier in value.modifiers:
                for prop in modifier.bl_rna.properties:
                    if prop.type == 'POINTER':
                        item = getattr(modifier, prop.identifier, None)
                        if isinstance(item, (bpy.types.Object, bpy.types.Collection, bpy.types.GeometryNodeTree)):
                            pending.append(item)
                if modifier.type == 'NODES' and modifier.node_group:
                    for socket in modifier.node_group.interface.items_tree:
                        if (socket.item_type == 'SOCKET' and socket.in_out == 'INPUT'
                                and socket.socket_type in ('NodeSocketObject', 'NodeSocketCollection')):
                            pending.append(getattr(modifier.properties.inputs, socket.identifier).value)
        elif isinstance(value, bpy.types.GeometryNodeTree):
            for node in value.nodes:
                if getattr(node, 'node_tree', None):
                    pending.append(node.node_tree)
                for socket in node.inputs:
                    item = getattr(socket, 'default_value', None)
                    if isinstance(item, (bpy.types.Object, bpy.types.Collection)):
                        pending.append(item)


class TerrainSampler:
    def __init__(self, terrain, controller, scale):
        if terrain.type != 'MESH' or terrain == controller or terrain.get('st_py_controller'):
            raise ValueError('Choose an original terrain mesh, not a generated street object')
        require_original_terrain(terrain)
        matrix = controller.matrix_world
        rotation = np.array(matrix.to_3x3())
        if not np.allclose(rotation.T@rotation,np.eye(3),atol=1e-6) or not np.allclose(rotation[:,2],(0,0,1),atol=1e-6):
            raise ValueError('Terrain placement allows controller translation and rotation about Z; apply no scale or tilt')
        evaluated = terrain.evaluated_get(bpy.context.evaluated_depsgraph_get())
        mesh = evaluated.to_mesh()
        try:
            xyz=np.empty((len(mesh.vertices),3),np.float32); mesh.vertices.foreach_get('co',xyz.ravel())
            transform=np.array(matrix.inverted()@evaluated.matrix_world)
            xyz=(xyz@transform[:3,:3].T+transform[:3,3])*scale
            mesh.calc_loop_triangles()
            triangles=np.empty((len(mesh.loop_triangles),3),np.int32)
            mesh.loop_triangles.foreach_get('vertices',triangles.ravel())
            if not len(triangles) or not np.isfinite(xyz).all():
                raise ValueError('Original terrain needs a finite mesh surface')
            self.signature=hashlib.sha256(xyz.tobytes()+triangles.tobytes()).hexdigest()
            self.top=float(xyz[:,2].max())+10.
            self.reach=float(np.ptp(xyz[:,2]))+20.
            self.bvh=BVHTree.FromPolygons(xyz.tolist(),triangles.tolist(),all_triangles=True)
        finally:
            evaluated.to_mesh_clear()
        self.count=0

    def __call__(self, xy):
        result=np.empty(len(xy))
        for i,(x,y) in enumerate(xy):
            hit=self.bvh.ray_cast(Vector((float(x),float(y),self.top)),Vector((0,0,-1)),self.reach)[0]
            if hit is None:
                raise ValueError(f'Original terrain has no surface below corridor point ({x:.3f}, {y:.3f}); extend the terrain or revise the street footprint')
            result[i]=hit.z
        self.count+=len(xy)
        return result


def configure(controller, app):
    props=controller.blosm_street_terrain
    selection=getattr(app,'terrainObject',None)
    if isinstance(selection,str):
        # Blosm may create an auxiliary flat terrain for other features even
        # when no original terrain was selected. Do not implicitly enable
        # placement against that import-rectangle-sized helper.
        target=bpy.context.scene.objects.get(selection.strip()) if selection.strip() else None
    else:
        # Programmatic callers can explicitly supply the standard wrapper.
        target=getattr(getattr(app,'terrain',None),'terrain',None)
    if not isinstance(target,bpy.types.Object):
        target=None
    props.target=target
    props.mode=getattr(app,'streetTerrainMode','PROFILE')
    for field,app_name in (('clearance','streetTerrainClearance'),('maximum_grade','streetTerrainGrade'),
                            ('sample_spacing','streetTerrainSampling'),('maximum_lift','streetTerrainLift')):
        setattr(props,field,getattr(app,app_name,getattr(props,field)))


def hydrate(update, owned, scale):
    products=dict(update.products)
    for key,product in products.items():
        if not isinstance(product,Part): continue
        if product.mesh is None:
            obj=owned.get(key)
            if obj is None: raise ValueError('A saved street mesh is missing; use a full regeneration')
            xyz=np.empty((len(obj.data.vertices),3),np.float32)
            obj.data.vertices.foreach_get('co',xyz.ravel())
            xyz=xyz.astype(float)*scale
            attr=obj.data.attributes.get(RAW_Z)
            if attr:
                z=np.empty(len(xyz),np.float32); attr.data.foreach_get('value',z); xyz[:,2]=z
            ids=np.empty(len(obj.data.loops),np.int32); obj.data.loops.foreach_get('vertex_index',ids)
            role=product.metadata.get('terrain_material')
            if role is None:
                grass=obj.get('st_py_controller').get('st_py_grass_material')
                role='grass' if grass is not None and grass in obj.data.materials[:] else 'pavement'
            mesh=Mesh(xyz,ids.reshape(-1,4),role)
        else:
            mesh=Mesh(publication_coordinates(product.mesh.xyz,scale,physical=True),product.mesh.quads,
                      product.mesh.material,product.mesh.seams)
        products[key]=Part(product.owner,mesh,dict(product.metadata,terrain_material=mesh.material))
    return replace(update,products=products)


class ScenePlacement:
    def __init__(self,controller,network,update,owned,scale):
        self.controller,self.scale,self.raw=controller,scale,update
        self.datums={}; self.group=None; self.report=None
        props=getattr(controller,'blosm_street_terrain',None)
        self.target=props.target if props else None
        self.mode=props.mode if self.target else 'NONE'
        old=controller.get('st_py_terrain_signature','')
        self.settings=PlacementSettings(props.clearance,props.maximum_grade,props.control_spacing,
            props.sample_spacing,props.maximum_lift) if props else PlacementSettings()
        start=time.perf_counter()
        sampler=TerrainSampler(self.target,controller,scale) if self.target else None
        self.signature=fingerprint((sampler.signature if sampler else None,self.mode,asdict(self.settings),scale,'1'))
        self.update=update
        changed=any(isinstance(update.products[k],(Part,AnchorTransform)) for k in update.rebuilt)
        if self.mode=='NONE' and controller.get('st_py_terrain_mode','NONE')=='NONE': return
        if old==self.signature and not changed: return
        self.raw=hydrate(update,owned,scale)
        refs,components,report=(prepare_corridor(network,self.raw.products,self.settings)
                               if self.mode!='NONE' else ([],[],{}))
        self.datums={Key.parse(c.identity):c.datum for c in components}
        products=dict(self.raw.products)
        fitted=None
        if self.mode=='PROFILE':
            fitted=fit(refs,components,sampler,self.settings)
            report.update(fitted.report)
        elif self.mode=='DIRECT':
            report.update(protected_bundles=[],bundle_policy='Each vertex is projected independently')
            # A complete footprint preflight makes initial publication atomic.
            # GN remains live afterwards; Update rechecks changed terrain.
            for c in components:
                for _,points in quad_samples(c.xyz,c.quads,self.settings.sample_spacing):
                    sampler(points[:,:,:2].reshape(-1,2))
            self.group=load_projection_group()
        for c in components:
            key=Key.parse(c.identity); raw=products[key]
            xyz=fitted.coordinates[c.identity] if fitted else raw.mesh.xyz
            products[key]=Part(raw.owner,Mesh(xyz,raw.mesh.quads,raw.mesh.material,raw.mesh.seams),
                dict(raw.metadata,terrain_datum=c.datum,terrain_reference=c.reference))
        by_key={c.identity:c for c in components}
        placements={p.id:p for p in network.furniture}
        for key,p in list(products.items()):
            if not isinstance(p,AnchorTransform): continue
            source=placements[key.id]; datum=network.section_index[source.section_id].datum
            if self.mode=='DIRECT':
                delta=float(sampler(np.array([p.position[:2]]))[0])+self.settings.clearance-datum
            elif fitted:
                road=self.raw.products[Key('road',source.section_id)]
                component=by_key.get(Key('road',source.section_id).text())
                if component is None or not road.metadata['start']-1e-7<=source.station<=road.metadata['end']+1e-7:
                    owner=next((span['owner'] for name in ('compound_spans','junction_spans','transition_spans')
                        for span in road.metadata.get(name,()) if span['start']-1e-7<=source.station<=span['end']+1e-7),
                        road.metadata.get('composite_owner'))
                    owner_product=self.raw.products.get(Key.parse(owner)) if owner else None
                    if getattr(owner_product,'pavement',None) is not None:
                        owner=owner_product.pavement.owner.text()
                    component=next((c for c in components if self.raw.products[Key.parse(c.identity)].owner.text()==owner),None)
                if component is None:
                    raise ValueError('No terrain reference owns furniture '+key.text())
                curve=fitted.curves[component.reference]
                station=source.station if len(curve.stations)>1 else curve.stations[0]
                delta=float(curve(station))-datum
            else: delta=0.
            products[key]=replace(p,position=(*p.position[:2],p.position[2]+delta))
        rebuilt=set(update.rebuilt)
        for key,product in products.items():
            obj=owned.get(key)
            if isinstance(product,Part):
                equal=False
                if obj is not None and len(obj.data.vertices)==len(product.mesh.xyz):
                    existing=np.empty((len(obj.data.vertices),3),np.float32);obj.data.vertices.foreach_get('co',existing.ravel())
                    equal=np.array_equal(existing,publication_coordinates(product.mesh.xyz,scale))
                if not equal or old!=self.signature or key in update.rebuilt: rebuilt.add(key)
            elif isinstance(product,AnchorTransform):
                if obj is None or not np.allclose(np.array(obj.location)*scale,product.position,rtol=0,atol=1e-6): rebuilt.add(key)
        # On a geometric edit a propagated fit can affect adjacent, unchanged
        # outputs. Publish those too, while committing only raw solver products.
        self.report=dict(report,mode=self.mode,seconds=time.perf_counter()-start,
                         terrain_samples=sampler.count if sampler else 0)
        self.update=replace(update,products=products,rebuilt=tuple(sorted(rebuilt)))
        self.update.timings['terrain']=self.report

    def decorate(self,obj,key):
        raw=self.raw.products[key]
        if isinstance(raw,Part):
            obj.data.attributes.new(RAW_Z,'FLOAT','POINT').data.foreach_set('value',raw.mesh.xyz[:,2])
            if self.mode=='DIRECT':
                modifier=obj.modifiers.new('Street terrain placement','NODES'); modifier.node_group=self.group
                values={'Terrain':self.target,'Prepared Datum':self.datums[key]/self.scale,
                        'Clearance':self.settings.clearance/self.scale,'Ray Margin':10./self.scale}
                for socket in self.group.interface.items_tree:
                    if socket.item_type=='SOCKET' and socket.in_out=='INPUT' and socket.name in values:
                        item=getattr(modifier.properties.inputs,socket.identifier)
                        item.value=values[socket.name]
                        if hasattr(item,'type'):item.type='VALUE'
        elif isinstance(raw,AnchorTransform):
            obj['st_terrain_prepared_position']=list(raw.position)


def load_projection_group():
    from .....app.asset_library import BlenderAssetLoader
    from ...assets import TERRAIN_PROJECTION, terrain_catalog
    loader=BlenderAssetLoader(Path(__file__).resolve().parents[4]/'assets',terrain_catalog())
    group=loader.load(TERRAIN_PROJECTION)
    if group.get('st_terrain_asset_revision')!='1.0':
        raise ValueError('Incompatible direct terrain projection asset')
    return group
