"""Load and assemble the newer native solver through Blosm's feature lifecycle."""
from dataclasses import asdict
import json
import time
import bpy
import numpy as np
from . import assets
from .renderer import StreetFeature
from .physical_sources import PhysicalStreetInventory
from .original_span_plan import OriginalSpanPlanner
from .original_span_authoring import OriginalSpanNetworkBuilder
from .transaction import ImportTransaction
from .coverage import CoverageReport
from .evidence import evaluated_mesh
from .surface_validation import minimum_positive_quad_corner
from ..way_properties import wayCategoryProps


def surface_evidence(obj):
    """Check generated output before declaring any imported route generated."""
    bpy.context.view_layer.update()
    with evaluated_mesh(obj) as mesh:
        count=len(mesh.polygons)
        if not count:
            raise ValueError('Native original-span construction rejected this component; inspect the shared GN diagnostic outputs')
        sizes=np.empty(count,dtype=np.int32);mesh.polygons.foreach_get('loop_total',sizes)
        if np.any(sizes!=4):raise ValueError('Original-span output contains a non-quad face')
        xyz=np.empty(len(mesh.vertices)*3);mesh.vertices.foreach_get('co',xyz);xyz=xyz.reshape(-1,3)
        faces=np.empty(count*4,dtype=np.int32);mesh.loops.foreach_get('vertex_index',faces);faces=faces.reshape(-1,4)
        minimum_corner=minimum_positive_quad_corner(xyz,faces)
        evidence=dict(quads=count,other_faces=0,vertices=len(xyz),positive_quads=True,
                      minimum_corner_blender_units=minimum_corner)
        if obj.get('blosm_aggregate_network_version'):
            from .aggregate_coverage import AggregateSurfaceCoverage
            fields = []
            for name in ('_cn_face_kind', '_cn_face_owner'):
                attribute = mesh.attributes.get(name)
                if attribute is None or attribute.domain != 'FACE' or attribute.data_type != 'INT':
                    raise ValueError('Aggregate output lost its exact surface ownership')
                field = np.empty(count, dtype=np.int32); attribute.data.foreach_get('value', field); fields.append(field)
            evidence['surface_coverage'] = AggregateSurfaceCoverage.from_object(obj).evaluate(zip(*fields))
        elif obj.modifiers[0].node_group.get('blosm_local_failures'):
            from .original_span_coverage import OriginalSurfaceCoverage
            fields=[]
            for name in ('st_cn_kind','st_cn_owner'):
                attribute=mesh.attributes.get(name)
                if attribute is None or attribute.domain!='FACE' or attribute.data_type!='INT':
                    raise ValueError('Native partial output lost its exact surface ownership')
                field=np.empty(count,dtype=np.int32);attribute.data.foreach_get('value',field);fields.append(field)
            coverage=OriginalSurfaceCoverage.from_object(obj)
            region_evidence=()
            shared_regions = fields[0] == 7
            if coverage.regions or np.any(shared_regions):
                attribute=mesh.attributes.get('st_cn_region_id')
                if attribute is None or attribute.domain!='FACE' or attribute.data_type!='INT':
                    raise ValueError('Native separator output lost its region ownership')
                regions=np.empty(count,dtype=np.int32);attribute.data.foreach_get('value',regions)
                selected=fields[0]==4
                if coverage.regions: region_evidence=zip(fields[1][selected],regions[selected])
                if np.any(regions[shared_regions] != fields[1][shared_regions]):
                    raise ValueError('Shared separator face differs from its declared region owner')
            evidence['surface_coverage']=coverage.evaluate(zip(*fields),
                matching_connections=coverage.matching_connections(obj),accepted_regions=region_evidence)
        return evidence


class OriginalSpanCarriagewayFeature(StreetFeature):
    """An explicit newer solver option; no input file or fixture IDs in policy."""
    required_libraries=(assets.ORIGINAL_LIBRARY,)
    ground_passages=True

    def __init__(self, *, physical_ids=None, interior_ids=(), expand=True, compound_spans=()):
        self.selection=dict(physical_ids=physical_ids,interior_ids=tuple(interior_ids),
                            expand=expand,compound_spans=tuple(tuple(s) for s in compound_spans))

    def render(self, context):
        started=time.perf_counter()
        inventory=PhysicalStreetInventory.view(context.snapshot)
        from .patch_programs import OriginalPatchCatalog, OriginalPatchPrograms
        from .intake_programs import OriginalIntakeCatalog, OriginalIntakePrograms
        program_catalog=OriginalPatchCatalog.load_optional(context.loader)
        intake_catalog=OriginalIntakeCatalog.load_optional(context.loader)
        preparation=patch_programs=None;program_refs={};preparation_seconds=0.
        if program_catalog is not None:
            preparing=time.perf_counter()
            preparation=program_catalog.prepare(inventory,context.design.categories,
                category_properties=wayCategoryProps,**self.selection)
            kinds={kind for component in preparation.components
                   for kind in OriginalPatchPrograms.for_component(component).values()}
            # The compatibility representative is also required by the current
            # standalone-only authoring contract; it is never evaluated there.
            program_refs=program_catalog.references(kinds or {'single'})
            preparation_seconds=time.perf_counter()-preparing
        requested=(assets.ORIGINAL_SOLVER,assets.ORIGINAL_PATCH_PROVIDER,
                   assets.ORIGINAL_STRIP_PROVIDER,assets.ORIGINAL_PAVEMENT,assets.ORIGINAL_GRASS)
        if program_catalog is not None:
            requested=tuple(ref for ref in requested if ref!=assets.ORIGINAL_PATCH_PROVIDER)+tuple(program_refs.values())
        # Appending must batch optional programs to avoid duplicate child IDs.
        # Retained linking reuses children across loads, so it can defer this
        # program until preparation proves that a free Bundle actually needs it.
        retained = context.loader.uses_link(assets.ORIGINAL_LIBRARY.key)
        requested+=assets.optional_original_programs(context.loader,bundles=bool(inventory.bundles),
            original=intake_catalog is None,defer_linked=True)
        if intake_catalog is not None:requested+=intake_catalog.preload_references(context.loader)
        loaded=context.loader.load_many(requested)
        solver=loaded[assets.ORIGINAL_SOLVER]
        assets.validate_original_solver(solver,cached=True,cached_strips=True,separators=True)
        if intake_catalog is not None and solver.get('blosm_original_provider_input')!=1:
            raise ValueError('The intake catalog requires native original-span provider transport')
        if program_catalog is not None:
            patch_programs=OriginalPatchPrograms({kind:loaded[ref] for kind,ref in program_refs.items()})
            program_catalog.validate_loaded(solver,patch_programs)
            patch_solver=patch_programs.representative
        else:
            patch_solver=loaded[assets.ORIGINAL_PATCH_PROVIDER]
        assets.validate_original_patch_provider(patch_solver)
        assets.validate_original_internal_regions(solver,patch_solver)
        internal_ends=assets.validate_original_internal_ends(solver,patch_solver)
        shared_owners=assets.validate_original_shared_owners(solver,patch_solver)
        strip_solver=loaded[assets.ORIGINAL_STRIP_PROVIDER];assets.validate_original_strip_provider(strip_solver)
        assets.validate_original_interpolation(solver,patch_solver,strip_solver)
        original_seams=assets.validate_original_seams(solver,patch_solver,strip_solver)
        category_connections=assets.validate_original_category_connections(solver,patch_solver,strip_solver)
        assets.validate_original_quad_energy(solver,patch_solver)
        diagnostic_solver=diagnostic_material=None
        if solver.get('blosm_local_failures'):
            assets.validate_original_failure_support(solver,strip_solver)
            extra=context.loader.load_many((assets.ORIGINAL_FAILURE_DISPLAY,assets.ORIGINAL_DIAGNOSTICS))
            diagnostic_solver=extra[assets.ORIGINAL_FAILURE_DISPLAY]
            assets.validate_original_failure_display(diagnostic_solver)
            diagnostic_material=extra[assets.ORIGINAL_DIAGNOSTICS]
        context.timings['original_span_asset_loading_seconds']=time.perf_counter()-started-preparation_seconds
        if preparation is None:
            preparing=time.perf_counter()
            preparation=OriginalSpanPlanner(context.design.categories,category_properties=wayCategoryProps,
                defer_physical_categories=bool(solver.get('blosm_filtered_physical_ports')),
                section_scopes=bool(solver.get('blosm_section_scopes')),
                internal_logical_ends=internal_ends,shared_logical_owners=shared_owners,
                seam_connections=original_seams,category_connections=category_connections,
                boundary_runs=bool(solver.get('blosm_boundary_runs')),
                dimension_runs=bool(solver.get('blosm_authored_dimension_runs')),
                free_bundles=bool(solver.get('blosm_free_bundle_surfaces'))).prepare(
                inventory,standalone=bool(solver.get('blosm_standalone_streets')),**self.selection)
            preparation_seconds=time.perf_counter()-preparing
        context.timings['original_span_preparation_seconds']=preparation_seconds
        for row in preparation.issues:
            context.snapshot.issue(row['kind'],row.get('id',-1),'ORIGINAL_SPAN_SCOPE',row['reason'])
        selected={sid for c in preparation.components for sid in c.ports.source_streets}
        for sid,street in inventory.streets.items():
            excluded=self.selection['physical_ids'] is not None and sid not in selected
            deferred=(street.get('reason') or '').startswith('Reserved component/category:')
            street['surface_status']='DELIBERATELY_EXCLUDED' if excluded or deferred else 'UNSUPPORTED'
            street['surface_reason']=('Outside explicit original-site selection' if excluded else
                street.get('reason') or 'No accepted original-span surface owner')
        junction_selection = None
        if self.selection['physical_ids'] is not None:
            junction_selection = {site.junction_id for component in preparation.components
                                  for site in component.ports.sites}
            junction_selection.update(jid for component in preparation.components
                                      for cluster in component.network.logical_clusters for jid in cluster.owner_ids)
        CoverageReport.initialize_junctions(context.snapshot,junction_selection)
        free_bundle_provider = None
        if solver.get('blosm_free_bundle_provider_input') and any(c.free_bundles for c in preparation.components):
            free_bundle_provider = loaded.get(assets.ORIGINAL_FREE_BUNDLE_PROVIDER)
            if free_bundle_provider is None and retained:
                load_started = time.perf_counter()
                free_bundle_provider = context.loader.load(assets.ORIGINAL_FREE_BUNDLE_PROVIDER)
                context.timings['original_span_asset_loading_seconds'] += time.perf_counter()-load_started
            if free_bundle_provider is None:
                raise ValueError('Prepared free Bundles require the declared shared provider asset')
        def ownership_provider():
            reference = assets.ORIGINAL_SEPARATOR_OWNERSHIP_PROVIDER
            if reference not in loaded and retained:
                load_started = time.perf_counter()
                loaded[reference] = context.loader.load(reference)
                context.timings['original_span_asset_loading_seconds'] += time.perf_counter()-load_started
            if reference not in loaded:
                raise ValueError('Prepared separators require the declared shared ownership provider asset')
            return loaded[reference]
        def fill_provider():
            reference=assets.ORIGINAL_SEPARATOR_FILL_PROVIDER
            if reference not in loaded and retained:
                load_started=time.perf_counter()
                loaded[reference]=context.loader.load(reference)
                context.timings['original_span_asset_loading_seconds']+=time.perf_counter()-load_started
            if reference not in loaded:
                raise ValueError('Prepared separators require the declared shared fill provider asset')
            return loaded[reference]
        def intake_loaded(seconds):
            context.timings['original_span_asset_loading_seconds']+=seconds
        builder=OriginalSpanNetworkBuilder(context.design,solver,loaded[assets.ORIGINAL_PAVEMENT],
            patch_solver=patch_solver if patch_programs is None else None,patch_programs=patch_programs,
            strip_solver=strip_solver,grass=loaded[assets.ORIGINAL_GRASS],
            diagnostic_solver=diagnostic_solver,diagnostic_material=diagnostic_material,
            free_bundle_provider=free_bundle_provider,
            original_intake_provider=(OriginalIntakePrograms(intake_catalog,context.loader,loaded,on_load=intake_loaded)
                if intake_catalog is not None else loaded.get(assets.ORIGINAL_INTAKE_PROVIDER)),
            separator_ownership_provider=ownership_provider if solver.get('blosm_separator_ownership_provider_input') else None,
            separator_fill_provider=fill_provider if solver.get('blosm_separator_fill_provider_input') else None)
        # Keep later source allocations from repeatedly invalidating a large
        # already evaluated component during the same import.
        for component in sorted(preparation.components,key=lambda c:len(c.ports.source_streets)):
            ports=component.ports;network=component.network
            scope_issues=OriginalSpanPlanner.issues_for(inventory,component,preparation.issues)
            result=dict(provider='original_spans.1',streets=list(ports.source_streets),
                intersections=[s.junction_id for s in ports.sites],
                logical_intersections=[jid for c in network.logical_clusters for jid in c.owner_ids],
                routes=[dict(street_id=r.span.street_id,section_ids=r.span.section_ids,kind=r.kind) for r in network.routes],
                frontiers=[asdict(f) for f in component.frontiers],
                deferred_ports=[asdict(p) for p in ports.deferred],
                compound_decisions=component.compound_decisions,
                scope_issues=scope_issues,
                unresolved_routes=len(network.unresolved),status='FAILED')
            if ports.section_scopes:
                result['section_scopes']=[asdict(scope) for scope in ports.section_scopes]
            if component.free_bundles:
                result['free_bundle_surfaces']=[asdict(row) for row in component.free_bundles]
            if ports.category_connections:
                result['category_connections']=[asdict(boundary) for boundary in ports.category_connections]
            if network.boundary_runs:
                result['boundary_runs']=[asdict(run) for run in network.boundary_runs]
            before=time.perf_counter()
            try:
                # Keep authored inputs and diagnostics when GN rejects geometry;
                # roll back allocations only for an actual assembly exception.
                with ImportTransaction():
                    label=('Streets / original network '+str(min(result['intersections'])) if result['intersections'] else
                           ('Streets / category-bounded ' if network.boundary_runs else 'Streets / standalone ')
                           +str(min(result['streets'])))
                    obj=builder.create(label,
                                       inventory,component,context.collection)
                context.outputs.append(obj)
                if obj.get('blosm_rejected_profile_count',0):
                    record=bpy.data.texts[obj['blosm_original_span_ownership']]
                    result['rejected_profiles']=json.loads(record.as_string())['rejected_profiles']
                result['evidence']=surface_evidence(obj)
                geometry=result['evidence'].get('surface_coverage',{})
                affected=set(geometry.get('affected_streets',()))
                failed=set(geometry.get('failed_streets',()))
                missing_junctions=set(geometry.get('missing_junctions',()))
                separator_plan=json.loads(obj.get('blosm_separator_preparation','{"closed":[],"pending":[]}'))
                result['separator_preparation']=separator_plan
                missing = set(geometry.get('missing_bundles', ()))
                closed = {row['bundle_id'] for row in separator_plan['closed']}.difference(missing)
                free = {row['bundle_id'] for row in separator_plan.get('free_ended', ())}.difference(missing)
                chained = {row['bundle_id'] for row in separator_plan.get('chained', ())}.difference(missing)
                free_chained = {row['bundle_id'] for row in separator_plan.get('free_chained', ())}.difference(missing)
                filled = closed | free | chained | free_chained
                if component.free_bundles:
                    accepted={tuple(key) for key in geometry.get('accepted_owners',())}
                    whole={row.bundle_id for row in component.free_bundles if (8,row.bundle_id) in accepted}
                    filled.update(whole)
                    result['free_bundle_surface_ids']=sorted(whole)
                result['closed_separator_ids']=sorted(closed)
                result['free_ended_separator_ids']=sorted(free)
                if 'chained' in separator_plan:result['chained_separator_ids']=sorted(chained)
                if 'free_chained' in separator_plan:result['free_chained_separator_ids']=sorted(free_chained)
                if 'crossing_regions' in separator_plan:
                    # A rendered region is not proof of the whole Bundle. Its
                    # remaining original extents stay explicit in preparation,
                    # and this capability does not add its Bundle to `filled`.
                    result['crossing_separator_region_ids']=sorted({row['bundle_id'] for row in
                        separator_plan['crossing_regions']}.difference(missing))
                if 'continuation_regions' in separator_plan:
                    pending={row['bundle_id'] for row in separator_plan.get('continuation_pending',())}
                    completed={row['bundle_id'] for row in separator_plan['continuation_regions']}.difference(missing|pending)
                    # The continuation planner revalidates the complete original
                    # member inventories. Every required primary and extra region
                    # must have distinct accepted FACE evidence before completion.
                    if 'missing_regions' not in geometry:completed=set()
                    filled.update(completed)
                    result['continued_separator_ids']=sorted(completed)
                if 'internal_cycle_regions' in separator_plan:
                    accepted = {tuple(key) for key in geometry.get('accepted_owners', ())}
                    required = {}
                    completed_regions = []
                    for row in separator_plan['internal_cycle_regions']:
                        rid = row['region_id']
                        if (7, rid) in accepted: completed_regions.append(rid)
                        for bid in row['bundle_ids']: required.setdefault(bid, set()).add((7, rid))
                    completed = {bid for bid, keys in required.items() if keys <= accepted and bid not in missing}
                    filled.update(completed)
                    result['internal_separator_region_ids'] = sorted(completed_regions)
                    result['internal_separator_bundle_ids'] = sorted(completed)
                bundle_pending=any(inventory.streets[sid].get('bundle_id') is not None and
                    inventory.streets[sid].get('bundle_id') not in filled for sid in ports.source_streets)
                reasons=[]
                if affected:reasons.append('Some native surface owners failed; unaffected geometry was retained')
                if network.unresolved:reasons.append('Original connections remain unresolved')
                if scope_issues:reasons.append('Logical intersection ownership remains unresolved')
                if bundle_pending:reasons.append('Bundle separator coverage pending')
                result['status']='PARTIAL' if reasons else 'READY'
                result['reason']='; '.join(reasons) or None
                for sid in ports.source_streets:
                    source=inventory.streets[sid]
                    scoped=next((scope for scope in ports.section_scopes if scope.street_id==sid),None)
                    incomplete=any(r.span.street_id==sid for r in network.unresolved)
                    # Longitudinal pavement is not proof of a completed Bundle
                    # separator. Keep those sources partial until it is supplied.
                    bundle_pending=source.get('bundle_id') is not None and source.get('bundle_id') not in filled
                    source['surface_status']=('FAILED' if sid in failed else 'PARTIAL'
                        if sid in affected or incomplete or bundle_pending or scope_issues or scoped else 'GENERATED')
                    source['surface_reason']=('Native surface owner failed' if sid in affected else
                        'Original routes unresolved' if incomplete else
                        'Bundle separator coverage pending' if bundle_pending else
                        'Logical intersection ownership unresolved' if scope_issues else
                        'Excluded category Sections retained without carriageway surfaces' if scoped else None)
                    if scoped:
                        source['surface_sections']=[dict(section_id=section,category=item['category'],
                            status=('DELIBERATELY_EXCLUDED' if section not in scoped.included_sections else
                                    'FAILED' if sid in failed else 'PARTIAL' if sid in affected or incomplete else 'GENERATED'))
                            for item in source['items'] if item['kind']=='Section' for section in (item['id'],)]
                for jid in result['intersections']+result['logical_intersections']:
                    if jid in context.snapshot.intersections:
                        context.snapshot.intersections[jid]['surface_status']='FAILED' if jid in missing_junctions else 'GENERATED'
            except (ValueError,RuntimeError) as exc:
                result['reason']=str(exc)
                for sid in ports.source_streets:
                    inventory.streets[sid].update(surface_status='FAILED',surface_reason=str(exc))
            result['seconds']=time.perf_counter()-before
            context.results.append(result)
