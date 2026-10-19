"""Author ordinary editable Streets and private original-span GN inputs."""
import json
from dataclasses import asdict
import bpy
from mathutils import Vector
from .authoring import SourceMeshBuilder,metre_scale
from .port_source_layout import PortSourceLayout
from .physical_dimensions import PreparedPhysicalDimensions
from .original_surface_design import OriginalSurfaceDesign
from ..gn.editable_input import curved_port_bindings,curved_port_source_bindings
from ..gn.kernel.pipeline import set_input


class OriginalSpanNetworkBuilder:
    def __init__(self, design, solver, pavement, *, patch_solver=None,strip_solver=None,grass=None,
                 diagnostic_solver=None,diagnostic_material=None,surface_design=None,source_organizer=None,
                 free_bundle_provider=None, separator_ownership_provider=None, patch_programs=None,
                 separator_fill_provider=None, original_intake_provider=None):
        if patch_programs is not None:
            if patch_solver is not None:raise ValueError('Select a patch solver or a finite shared program set')
            patch_solver=patch_programs.representative
        self.patch_programs=patch_programs
        self.design=design;self.solver=solver;self.pavement=pavement;self.patch_solver=patch_solver
        from .assets import validate_original_rejected_profiles
        self.retain_unresolved_profiles=validate_original_rejected_profiles(solver,patch_solver,strip_solver)
        self.surface_design=surface_design if surface_design is not None else OriginalSurfaceDesign(
            design,retain_unresolved_profiles=self.retain_unresolved_profiles)
        self.source_organizer=source_organizer
        self.strip_solver=strip_solver
        self.grass=grass
        self.free_bundle_provider=free_bundle_provider
        # A callable resolves a retained optional library only after the final
        # separator table exists. Appending callers can provide a batched tree.
        self.separator_ownership_provider=separator_ownership_provider
        self.separator_fill_provider=separator_fill_provider
        self.original_intake_provider=original_intake_provider
        if solver.get('blosm_original_provider_input'):
            if original_intake_provider is None:
                raise ValueError('Original-span construction requires its declared shared input provider asset')
            if not callable(original_intake_provider) and not callable(getattr(original_intake_provider,'resolve',None)):
                from .assets import validate_original_intake_provider
                validate_original_intake_provider(original_intake_provider,solver)
        self.diagnostic_solver=diagnostic_solver;self.diagnostic_material=diagnostic_material
        if strip_solver is not None and patch_solver is None:raise ValueError('Strip providers require patch providers')
        from .assets import validate_original_shared_owners
        self.shared_owners=validate_original_shared_owners(solver,patch_solver) if patch_solver is not None else False
        from .assets import validate_original_internal_regions
        self.internal_regions=validate_original_internal_regions(solver,patch_solver)

    @staticmethod
    def coordinate_origin(inventory, plan, scale):
        """Return one representable world translation for sources and records."""
        location=(plan.sites[0].location if plan.sites else
                  next(i for i in inventory.streets[plan.source_streets[0]]['items'] if i['kind']=='Section')['points'][0])
        return Vector(tuple(v/scale for v in location[:2])+(0.,))

    def create(self, name, inventory, component, collection):
        if self.patch_programs is not None:self.patch_programs.validate_component(component)
        plan=component.ports;network=component.network
        free_bundles=component.free_bundles
        if free_bundles and self.solver.get('blosm_free_bundle_provider_input'):
            from .free_bundle_provider import FreeBundleProvider
            if self.free_bundle_provider is None:
                raise ValueError('The declared free Bundle requires its shared native provider asset')
            free_provider = FreeBundleProvider(self.free_bundle_provider)
        else:
            free_provider = None
        if free_bundles and (self.solver.get('blosm_free_bundle_surfaces') != 1
                            or self.solver.get('blosm_free_bundle_guards',0) < 3):
            raise ValueError('The loaded shared asset has no guarded free Bundle surface branch')
        roots=[self.solver,self.strip_solver]+(list(self.patch_programs.trees.values())
                if self.patch_programs is not None else [self.patch_solver])
        if any(port.source_kind=='THROUGH' and port.connector_id is not None
               for site in plan.sites for port in site.ports):
            if any(tree is not None and tree.get('blosm_authored_through_connectors') != 1 for tree in roots):
                raise ValueError('The loaded shared asset does not support authored through Connectors')
        if plan.selection_boundaries:
            if any(tree is not None and tree.get('blosm_selection_boundaries') != 1 for tree in roots):
                raise ValueError('The loaded shared asset does not support declared selection boundaries')
            if any(plan.selection_peer(row.endpoint) is not None for row in plan.selection_boundaries):
                if any(tree is not None and tree.get('blosm_interior_selection_cuts') != 1 for tree in roots):
                    raise ValueError('The loaded shared asset does not support paired interior selection cuts')
        declarations=PreparedPhysicalDimensions().read_spans(inventory,plan.spans,source_plan=plan)
        for row in declarations:
            routes=[r for r in network.routes if r.span.street_id==row.street_id
                    and (row.before_section,row.after_section) in tuple(zip(r.span.section_ids,r.span.section_ids[1:]))]
            if len(routes)!=1:raise ValueError('A dimensional transition requires one original surface owner')
            kind=routes[0].kind
            capability=('blosm_free_bundle_surfaces' if kind=='FREE_BUNDLE' else
                        'blosm_authored_dimension_runs' if kind in ('STANDALONE','BOUNDARY_ONLY') else
                        'blosm_authored_dimension_approaches' if kind in ('APPROACH','CATEGORY_APPROACH','SELECTION_APPROACH') else None)
            if capability:
                # A free Bundle is constructed wholly by the network root.
                # Patch/strip providers do not read its transition table and
                # correctly do not advertise this independent branch.
                capability_roots = [self.solver] if kind == 'FREE_BUNDLE' else roots
                if any(tree is not None and tree.get(capability)!=1 for tree in capability_roots):
                    raise ValueError('The loaded shared asset does not support dimensional '+kind.lower())
            elif row.source_kind=='AuthoredTransition':
                raise ValueError('The loaded shared asset does not support authored dimensional '+kind.lower())
        layout=PortSourceLayout.prepare(inventory,plan)
        scale=metre_scale()
        # Use the same representable controller translation for source objects
        # and private records; otherwise subtracting two different rounded
        # origins can introduce a false seam mismatch at large map coordinates.
        translation=self.coordinate_origin(inventory,plan,scale)
        origin=tuple(v*scale for v in translation)
        boundary_runs=bool(self.solver.get('blosm_boundary_runs'))
        if network.boundary_runs and not boundary_runs:
            raise ValueError('The loaded original-span asset has no category-bounded run branch')
        compiled=curved_port_bindings.compile(inventory,plan,name,origin=origin,endpoint_ownership=boundary_runs)
        port_designs=tuple(self.surface_design.ports(inventory,plan))
        approaches=curved_port_bindings.surface_approaches(plan,compiled,
            port_designs,name,retain_unresolved=self.retain_unresolved_profiles)
        dimensions=curved_port_bindings.dimension_transitions(plan,compiled,
            declarations,name)
        tables=curved_port_bindings.network_tables(plan,network,compiled,approaches,name,
            allow_unresolved=True,dimensions=dimensions,logical_clusters=True,shared_logical=self.shared_owners)
        if free_bundles:
            from ..gn.editable_input.free_bundle_bindings import path_rows
            from ..gn.editable_input.curved_junction_bindings import _mesh
            profiles={sid:self.surface_design.span(inventory,span) for span in plan.spans for sid in span.section_ids}
            clips={(row.endpoint.section_id,row.endpoint.end):row.trim_m for row in plan.selection_boundaries}
            rows=path_rows([inventory.streets[sid] for sid in plan.source_streets],free_bundles,profiles,
                           clips=clips,transitions=declarations)
            tables['Bundle Paths']=_mesh(name+' private free Bundle paths',[(0.,0.,0.)]*len(rows),[],[
                (key,'BOOLEAN' if isinstance(value,bool) else 'FLOAT' if isinstance(value,float) else 'INT',
                 'POINT',[row[key] for row in rows]) for key,value in rows[0].items()])
        free_routes=[route for route in network.routes if route.kind=='STANDALONE']
        accepts_free=any(s.item_type=='SOCKET' and s.in_out=='INPUT' and s.name=='Standalone Streets'
                         for s in self.solver.interface.items_tree)
        if free_routes and not accepts_free:
            raise ValueError('The loaded original-span asset has no standalone Street branch')
        if accepts_free:
            from ..gn.editable_input.curved_junction_bindings import _mesh
            rows=[]
            for route in free_routes:
                rise,datum=self.surface_design.span(inventory,route.span)
                rows.append(dict(span=plan.spans.index(route.span),street_id=route.span.street_id,
                                 rise=rise,datum=datum))
            tables['Standalone Streets']=_mesh(name+' private standalone ownership',[(0.,0.,0.)]*len(rows),[],[
                ('st_fs_'+key,kind,'POINT',[r[key] for r in rows]) for key,kind in
                (('span','INT'),('street_id','INT'),('rise','FLOAT'),('datum','FLOAT'))])
        if boundary_runs:
            from ..gn.editable_input.curved_junction_bindings import _mesh
            rows=[]
            indices={span.identity:i for i,span in enumerate(plan.spans)}
            for run in network.boundary_runs:
                span=plan.spans[indices[run.span_identity]]
                rise,datum=self.surface_design.span(inventory,span)
                row=dict(span=indices[run.span_identity],street_id=run.span_identity[0],
                         rise=rise,datum=datum)
                for end,side in enumerate(('start','end')):
                    row.update({side+'_role':run.roles[end],side+'_owner':run.owners[end],
                                side+'_next_section':run.next_sections[end]})
                rows.append(row)
            schema=[('span','INT'),('street_id','INT'),('rise','FLOAT'),('datum','FLOAT')]
            schema.extend((side+'_'+key,'INT') for side in ('start','end') for key in ('role','owner','next_section'))
            tables['Boundary Runs']=_mesh(name+' private category-bounded runs',[(0.,0.,0.)]*len(rows),[],[
                ('st_br_'+key,kind,'POINT',[r[key] for r in rows]) for key,kind in schema])
        separators=None;regions=None;continuations=None;continuation_metadata=[]
        cycles=None;cycle_metadata=[]
        if self.internal_regions:
            from .original_span_separator_cycles import InternalCycleRegionPlanner
            cycles=InternalCycleRegionPlanner(self.design.categories).prepare(inventory,component)
        if any(s.item_type=='SOCKET' and s.in_out=='INPUT' and s.name=='Separators' for s in self.solver.interface.items_tree):
            from .original_span_separators import OriginalSpanSeparatorPlanner
            from ..gn.editable_input.curved_separator_bindings import bind
            free_ends = bool(self.solver.get('blosm_free_ended_separators'))
            if self.solver.get('blosm_chained_separators'):
                from .original_span_separator_chains import ChainedSeparatorPlanner
                from ..gn.editable_input.curved_separator_chain_bindings import bind_mixed
                free_chains = bool(self.solver.get('blosm_free_chained_separators'))
                longitudinal = bool(self.solver.get('blosm_longitudinal_separator_chains'))
                separators=ChainedSeparatorPlanner(self.design.categories,longitudinal_chains=longitudinal).prepare(inventory,component,
                    free_ends=free_ends,free_chains=free_chains)
                if self.solver.get('blosm_crossing_separator_regions'):
                    from .original_span_separator_regions import CrossingSeparatorRegionPlanner
                    from ..gn.editable_input.curved_separator_region_bindings import bind_mixed as bind_regions
                    regions=CrossingSeparatorRegionPlanner(self.design.categories).prepare(inventory,component,
                        bundle_ids=(row['bundle_id'] for row in separators.pending))
                    tables['Separators'],tables['Separator Owners']=bind_regions(
                        separators,regions,tables,name,free_chains=free_chains,longitudinal_chains=longitudinal)
                    if self.solver.get('blosm_continuation_separator_regions'):
                        from .original_span_separator_continuations import CrossingFreeRegionPlanner
                        from ..gn.editable_input.curved_separator_continuation_bindings import bind as bind_continuations
                        continuations=CrossingFreeRegionPlanner(self.design.categories).prepare_continuations(
                            inventory,component,prefixes=regions.regions)
                        tables['Continuation Regions'],tables['Continuation Owners'],continuation_metadata=bind_continuations(
                            continuations,tables,name,original_bundle_ids=inventory.bundles.keys())
                else:
                    tables['Separators'],tables['Separator Owners']=bind_mixed(separators,tables,name,free_chains=free_chains,
                        longitudinal_chains=longitudinal)
            else:
                separators=OriginalSpanSeparatorPlanner(self.design.categories).prepare(inventory,component,free_ends=free_ends)
                tables['Separators']=(bind(separators,tables,name,free_ends=True) if free_ends
                                      else bind(separators,tables,name))
        if cycles is not None:
            from ..gn.editable_input.curved_internal_cycle_bindings import bind as bind_cycles
            from ..gn.editable_input.curved_separator_region_identity import RegionNamespaces
            from ..gn.editable_input.curved_patch_provider_bindings import values
            primary=[int(v) for v in values(tables['Separators'],'st_sb_bundle_id')]
            if set(primary).intersection(bid for row in cycles.regions for bid in row.bundle_ids):
                raise ValueError('A complete internal Bundle cannot also declare an external separator')
            namespaces=RegionNamespaces.allocate(inventory.bundles,primary,
                [asdict(row) for row in continuations.regions] if continuations is not None else (),
                internal_rows=[asdict(row) for row in cycles.regions])
            # The extended allocation preserves existing continuation IDs;
            # verify this before writing any additional region references.
            if any(row['region_id']!=namespaces.continuation(row) for row in continuation_metadata):
                raise ValueError('Shared namespace allocation changed an existing continuation identity')
            tables['Internal Cycles'],cycle_metadata=bind_cycles(cycles,tables,name,namespaces=namespaces)
        tables.pop('Sources')
        tables['Source Binding']=curved_port_source_bindings.bind(compiled,layout,name)
        private=bpy.data.collections.new(name+' / private topology');collection.children.link(private)
        sources=bpy.data.collections.new(name+' / editable Streets');collection.children.link(sources)
        obj=bpy.data.objects.new(name,bpy.data.meshes.new(name+' input'));collection.objects.link(obj)
        obj.location=translation
        bpy.context.view_layer.update()
        obj['blosm_original_span_network']=True
        obj['blosm_projection_coordinate_offset']=origin
        obj['blosm_source_schema']='original_spans.1'
        obj['blosm_logical_owner_ids']=json.dumps({str(c.logical_junction_id):c.owner_ids for c in network.logical_clusters})
        if free_bundles:
            obj['blosm_free_bundle_preparation']=json.dumps([asdict(row) for row in free_bundles])
        if separators is not None:
            declaration=asdict(separators)
            if not declaration['deferred_interiors']:declaration.pop('deferred_interiors')
            if regions is not None:
                declaration['crossing_regions']=[asdict(row) for row in regions.regions]
                declaration['crossing_region_pending']=list(regions.pending)
            if continuations is not None:
                declaration['continuation_regions']=continuation_metadata
                declaration['continuation_pending']=list(continuations.pending)
            if cycles is not None:
                declaration['internal_cycle_regions']=cycle_metadata
                declaration['internal_cycle_pending']=list(cycles.pending)
                if cycles.not_applicable:
                    declaration['internal_cycle_not_applicable']=list(cycles.not_applicable)
            obj['blosm_separator_preparation']=json.dumps(declaration)
        obj.modifiers.new('Shared original-span streets','NODES').node_group=self.solver
        obj.modifiers[0].show_viewport=False
        for key,mesh in tables.items():
            data=bpy.data.objects.new(name+' / '+key,mesh);private.objects.link(data)
            data.hide_set(True);data.hide_render=True;set_input(obj,key,data)
        vertex=edge=0
        for sid in plan.source_streets:
            sections=[dict(s,points=[tuple(p)+(0.,)*(3-len(p)) for p in s['points']],
                           source_way_ids=s.get('source_way_ids',()))
                      for s in inventory.streets[sid]['items'] if s['kind']=='Section']
            source=SourceMeshBuilder().create(name+' / Street '+str(sid),sections,sources,obj,self.design.interpolation)
            modes=[self.surface_design.interpolation(section) for section in sections]
            if any(type(mode) is not int or mode not in (0,3) for mode in modes):
                raise ValueError('Surface design must retain Exact or Smooth Section interpolation')
            source.data.attributes['st_interpolation'].data.foreach_set('value',
                [mode for section,mode in zip(sections,modes) for _ in range(len(section['points'])-1)])
            source['blosm_street_id']=sid
            for key,domain,start,count in (('st_live_vertex_index','POINT',vertex,len(source.data.vertices)),
                                          ('st_live_edge_index','EDGE',edge,len(source.data.edges))):
                source.data.attributes.new(key,'INT',domain).data.foreach_set('value',list(range(start,start+count)))
            vertex+=len(source.data.vertices);edge+=len(source.data.edges)
        assert (vertex,edge)==(layout.source_vertex_count,layout.source_edge_count)
        for key,value in (('Source Objects',sources),('Source Vertex Count',vertex),('Source Edge Count',edge),
            ('Metres to Blender Units',1./scale),('Cell Size',self.design.cell_size),
            ('Mesh Resolution',self.design.resolution),('Height Exaggeration',self.design.exaggeration),
            ('Pavement',self.pavement),('Samples Per Side',129)):
            set_input(obj,key,value)
        if separators is not None:
            set_input(obj,'Grass',self.grass)
            # The bound table includes every supported declaration kind. Using
            # only the closed subset would reject a mixed closed/free network.
            set_input(obj,'Separator Count',len(tables['Separators'].vertices))
        if continuations is not None:
            set_input(obj,'Continuation Count',len(tables['Continuation Regions'].vertices))
        if accepts_free:set_input(obj,'Standalone Count',len(free_routes))
        if free_bundles:
            set_input(obj,'Bundle Count',len(free_bundles))
            set_input(obj,'Grass',self.grass)
        if boundary_runs:set_input(obj,'Boundary Run Count',len(network.boundary_runs))
        topology=bpy.data.texts.new(name+' / input ownership')
        # This record is referenced by name, not by an ID pointer.
        topology.use_fake_user = True
        topology.write(json.dumps(dict(ports=asdict(plan),network=asdict(network),
            layout=asdict(layout),frontiers=[asdict(f) for f in component.frontiers],
            compound_decisions=component.compound_decisions,
            rejected_profiles=[asdict(d) for d in port_designs if getattr(d,'profile_error',None)]),indent=2))
        obj['blosm_original_span_ownership']=topology.name
        obj['blosm_rejected_profile_count']=sum(bool(getattr(d,'profile_error',None)) for d in port_designs)
        # Intermediate meshes are never authored inputs or retained solver state.
        for mesh in (compiled.sources,compiled.ports,approaches):
            if mesh.users==0:bpy.data.meshes.remove(mesh)
        # Provider dependency collections bind actual source objects. Organize
        # them before installing providers; regrouping afterward would leave
        # stale per-owner dependencies even if all live indices were retained.
        if self.source_organizer is not None:
            self.source_organizer(obj)
        if free_provider is not None:
            free_provider.install(obj, private)
        if self.patch_solver is not None:
            if self.patch_programs is not None:
                _,providers=self.patch_programs.bind(obj,network_tree=self.solver,collection=collection)
            else:
                from ..gn.editable_input.curved_patch_provider_bindings import install
                _,providers=install(obj,provider_tree=self.patch_solver,network_tree=self.solver,collection=collection)
            obj['blosm_cached_patch_owners']=len(providers)
        if self.strip_solver is not None:
            from ..gn.editable_input.curved_strip_provider_bindings import install
            _,providers=install(obj,provider_tree=self.strip_solver,network_tree=self.solver,collection=collection)
            obj['blosm_cached_strip_owners']=len(providers)
        if (self.solver.get('blosm_separator_ownership_provider_input') and separators is not None
                and len(tables['Separators'].vertices)):
            from .separator_ownership_provider import SeparatorOwnershipProvider
            provider = self.separator_ownership_provider
            if callable(provider): provider = provider()
            if provider is None:
                raise ValueError('Declared separators require the shared ownership provider asset')
            SeparatorOwnershipProvider(provider).install(obj, private)
        if (self.solver.get('blosm_separator_fill_provider_input') and separators is not None
                and len(tables['Separators'].vertices)):
            from .separator_fill_provider import SeparatorFillProvider
            provider=self.separator_fill_provider
            if callable(provider):provider=provider()
            if provider is None:
                raise ValueError('Declared separators require the shared fill provider asset')
            SeparatorFillProvider(provider).install(obj,private)
        if self.solver.get('blosm_original_provider_input'):
            from .original_intake_provider import OriginalIntakeProvider
            provider=self.original_intake_provider
            if callable(getattr(provider,'resolve',None)):provider=provider.resolve(obj)
            elif callable(provider):provider=provider()
            if provider is None:
                raise ValueError('Original-span construction requires its declared shared input provider asset')
            OriginalIntakeProvider(provider).install(obj,private)
        obj.modifiers[0].show_viewport=True
        if self.diagnostic_solver is not None:
            diagnostic=bpy.data.objects.new(name+' / Live failure labels',bpy.data.meshes.new(name+' diagnostic host'))
            collection.objects.link(diagnostic);diagnostic.parent=obj
            diagnostic.modifiers.new('Shared native failure labels','NODES').node_group=self.diagnostic_solver
            for key,value in (('Network',obj),('Material',self.diagnostic_material),('Metres to Blender Units',1./scale)):
                set_input(diagnostic,key,value)
            obj['blosm_failure_display']=diagnostic.name
        return obj
