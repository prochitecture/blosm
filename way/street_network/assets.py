"""Street asset declarations; library resolution/loading is feature-independent."""
from ...app.asset_library import AssetCatalog, AssetRef, LibrarySpec

LIBRARY = LibrarySpec('streets.native.v8', 'prochitecture_streets_v8.blend', revision='8.14')
SOLVER = AssetRef(LIBRARY.key, 'node_groups', 'Prochitecture Streets / Editable network v8')
PAVEMENT = AssetRef(LIBRARY.key, 'materials', 'Prochitecture Streets / Pavement')
GRASS = AssetRef(LIBRARY.key, 'materials', 'Prochitecture Streets / Grass')

ORIGINAL_LIBRARY = LibrarySpec('streets.original_spans.v1', 'prochitecture_streets_original_spans.blend', revision='1.14')
ORIGINAL_SOLVER = AssetRef(ORIGINAL_LIBRARY.key, 'node_groups', 'Prochitecture Streets / Original span network v1')
ORIGINAL_PATCH_PROVIDER = AssetRef(ORIGINAL_LIBRARY.key, 'node_groups', 'Prochitecture Streets / Original patch provider v1')
ORIGINAL_PATCH_CATALOG = AssetRef(ORIGINAL_LIBRARY.key, 'texts', 'Blosm / Street constructor catalog v1')
ORIGINAL_STRIP_PROVIDER = AssetRef(ORIGINAL_LIBRARY.key, 'node_groups', 'Prochitecture Streets / Original strip provider v1')
ORIGINAL_AGGREGATE_PROVIDER = AssetRef(ORIGINAL_LIBRARY.key, 'node_groups', 'Prochitecture Streets / Aggregate provider v1')
ORIGINAL_INTAKE_PROVIDER = AssetRef(ORIGINAL_LIBRARY.key, 'node_groups', 'Prochitecture Streets / Original intake provider v1')
ORIGINAL_INTAKE_CATALOG = AssetRef(ORIGINAL_LIBRARY.key, 'texts', 'Blosm / Street intake catalog v1')
ORIGINAL_FREE_BUNDLE_PROVIDER = AssetRef(ORIGINAL_LIBRARY.key, 'node_groups', 'Prochitecture Streets / Free Bundle provider v1')
ORIGINAL_SEPARATOR_OWNERSHIP_PROVIDER = AssetRef(ORIGINAL_LIBRARY.key, 'node_groups', 'Prochitecture Streets / Separator ownership provider v1')
ORIGINAL_SEPARATOR_FILL_PROVIDER = AssetRef(ORIGINAL_LIBRARY.key, 'node_groups', 'Prochitecture Streets / Separator fill provider v1')
ORIGINAL_PAVEMENT = AssetRef(ORIGINAL_LIBRARY.key, 'materials', PAVEMENT.name)
ORIGINAL_GRASS = AssetRef(ORIGINAL_LIBRARY.key, 'materials', GRASS.name)
ORIGINAL_FAILURE_DISPLAY = AssetRef(ORIGINAL_LIBRARY.key, 'node_groups', 'Prochitecture Streets / Original failure display v1')
ORIGINAL_DIAGNOSTICS = AssetRef(ORIGINAL_LIBRARY.key, 'materials', 'Prochitecture Streets / Diagnostics')


def catalog(external_path=None):
    from dataclasses import replace
    library = replace(LIBRARY, path=external_path, source='external') if external_path else LIBRARY
    original = replace(ORIGINAL_LIBRARY, path=external_path, source='external') if external_path else ORIGINAL_LIBRARY
    return AssetCatalog((library, original))


def optional_original_programs(loader, *, bundles=False, aggregate=False, original=True, defer_linked=False):
    """Batch optional roles so appended programs share their child datablocks.

    Retained-link consumers may resolve Bundle roles after final preparation.
    Appending consumers must include them in the common graph load; a second
    append can duplicate shared recipes. Availability/edition errors propagate.
    """
    requested=[]
    if original:requested.append(ORIGINAL_INTAKE_PROVIDER)
    if aggregate:requested.append(ORIGINAL_AGGREGATE_PROVIDER)
    if bundles and not (defer_linked and loader.uses_link(ORIGINAL_LIBRARY.key)):
        requested.extend((ORIGINAL_FREE_BUNDLE_PROVIDER,ORIGINAL_SEPARATOR_OWNERSHIP_PROVIDER,
                          ORIGINAL_SEPARATOR_FILL_PROVIDER))
    return tuple(ref for ref in requested if loader.contains(ref))


def validate_original_solver(group, *, cached=False,cached_strips=False,separators=False):
    _original_interpolation_revision(group)
    if (group.bl_idname != 'GeometryNodeTree' or
            group.get('blosm_street_asset_version') != ORIGINAL_LIBRARY.revision or
            group.get('blosm_street_input_schema') != 'original_spans.1'):
        raise ValueError('Incompatible original-span street GN library')
    required = {'Geometry':'NodeSocketGeometry','Source Objects':'NodeSocketCollection',
        'Source Binding':'NodeSocketObject','Source Vertex Count':'NodeSocketInt',
        'Source Edge Count':'NodeSocketInt','Cell Size':'NodeSocketFloat',
        'Mesh Resolution':'NodeSocketInt','Height Exaggeration':'NodeSocketFloat',
        'Metres to Blender Units':'NodeSocketFloat','Pavement':'NodeSocketMaterial'}
    required.update({k:'NodeSocketObject' for k in ('Mapping','Approaches','Junctions',
        'Corridors','External Approaches','Compounds','Transitions','Logical Clusters',
        'Cluster Sites','Boundary Links','Cluster Spans')})
    if cached or cached_strips:required['Patch Results']='NodeSocketCollection'
    if cached_strips:required['Strip Results']='NodeSocketCollection'
    if separators:
        required.update({'Separators':'NodeSocketObject','Separator Count':'NodeSocketInt',
                         'Grass':'NodeSocketMaterial'})
    free_ends = group.get('blosm_free_ended_separators', 0)
    chained = group.get('blosm_chained_separators', 0)
    if chained not in (0,1):raise ValueError('Unknown chained separator capability revision')
    if chained:
        if free_ends != 1 or group.get('blosm_local_failures') != 1:
            raise ValueError('Chained separators require free-ended dispatch and native local failure support')
        required['Separator Owners']='NodeSocketObject'
    free_chains = group.get('blosm_free_chained_separators', 0)
    if free_chains not in (0,1):raise ValueError('Unknown free-chain separator capability revision')
    if free_chains and chained != 1:
        raise ValueError('Free-chain separators require chained separator dispatch')
    longitudinal = group.get('blosm_longitudinal_separator_chains', 0)
    if longitudinal not in (0, 1):
        raise ValueError('Unknown longitudinal separator capability revision')
    if longitudinal and chained != 1:
        raise ValueError('Longitudinal separators require chained separator dispatch')
    regions = group.get('blosm_crossing_separator_regions', 0)
    if regions not in (0,1):raise ValueError('Unknown crossing-region capability revision')
    if regions and chained != 1:
        raise ValueError('Crossing regions require chained separator dispatch')
    continuations=group.get('blosm_continuation_separator_regions',0)
    if continuations not in (0,1):raise ValueError('Unknown continuation-region capability revision')
    if continuations:
        if regions!=1 or group.get('blosm_local_failures')!=1:
            raise ValueError('Continuation regions require crossing regions and local failure support')
        required.update({'Continuation Regions':'NodeSocketObject','Continuation Owners':'NodeSocketObject',
                         'Continuation Count':'NodeSocketInt'})
    internal_regions=group.get('blosm_internal_cycle_regions',0)
    if internal_regions not in (0,1):raise ValueError('Unknown internal-region capability revision')
    if internal_regions:
        if group.get('blosm_local_failures')!=1:
            raise ValueError('Internal regions require cached patches, separators and local failure support')
        required.update({'Internal Cycles':'NodeSocketObject','Patch Results':'NodeSocketCollection',
                         'Separators':'NodeSocketObject','Grass':'NodeSocketMaterial'})
    boundary_runs=group.get('blosm_boundary_runs',0)
    if boundary_runs not in (0,1):raise ValueError('Unknown category-bounded run capability revision')
    if boundary_runs:
        required.update({'Boundary Runs':'NodeSocketObject','Boundary Run Count':'NodeSocketInt'})
        if any(group.get(key)!=1 for key in ('blosm_section_scopes','blosm_standalone_streets','blosm_local_failures')):
            raise ValueError('Category-bounded runs require Section scopes, standalone and local failure support')
    if free_ends not in (0, 1):
        raise ValueError('Unknown free-ended separator capability revision')
    if free_ends:
        required.update({'Separators':'NodeSocketObject', 'Separator Count':'NodeSocketInt',
                         'Grass':'NodeSocketMaterial', 'Patch Results':'NodeSocketCollection',
                         'Strip Results':'NodeSocketCollection'})
        if group.get('blosm_local_failures') != 1:
            raise ValueError('Free-ended separators require native local failure support')
    if group.get('blosm_standalone_streets'):
        required['Standalone Streets']='NodeSocketObject'
        required['Standalone Count']='NodeSocketInt'
        if group.get('blosm_local_failures')!=1:
            raise ValueError('Standalone Streets require native local failure support')
    aggregate_provider = group.get('blosm_aggregate_provider_input', 0)
    original_provider = group.get('blosm_original_provider_input', 0)
    if original_provider not in (0, 1):
        raise ValueError('Unknown original-span provider transport revision')
    if original_provider:
        if group.get('blosm_prepared_network_contracts') != 'original_spans.1,aggregate_ports.1':
            raise ValueError('Original-span transport must retain both preparation contracts')
        required.update({'Original Span Result':'NodeSocketObject', 'Original Span Provider Key':'NodeSocketInt',
                         'Prepared Network Kind':'NodeSocketInt'})
    ownership_provider = group.get('blosm_separator_ownership_provider_input', 0)
    if ownership_provider not in (0, 1):
        raise ValueError('Unknown separator ownership provider transport revision')
    if ownership_provider:
        if chained != 1 or group.get('blosm_local_failures') != 1:
            raise ValueError('Ownership transport requires chained separators and local failure support')
        required.update({'Separator Ownership Result':'NodeSocketObject', 'Separator Ownership Key':'NodeSocketInt'})
    fill_provider=group.get('blosm_separator_fill_provider_input',0)
    if fill_provider not in (0,1):
        raise ValueError('Unknown separator fill provider transport revision')
    if fill_provider:
        if ownership_provider!=1 or group.get('blosm_local_failures')!=1:
            raise ValueError('Fill transport requires shared separator ownership and local failure support')
        required.update({'Complete Separator Fill':'NodeSocketObject','Partial Separator Fill':'NodeSocketObject',
                         'Separator Fill Key':'NodeSocketInt'})
    free_bundle_provider = group.get('blosm_free_bundle_provider_input', 0)
    if free_bundle_provider not in (0, 1):
        raise ValueError('Unknown free-Bundle provider transport revision')
    if free_bundle_provider:
        if group.get('blosm_free_bundle_surfaces') != 1 or group.get('blosm_free_bundle_guards', 0) < 3:
            raise ValueError('Free-Bundle provider transport requires the complete guarded branch')
        required.update({'Free Bundle Result':'NodeSocketObject', 'Free Bundle Provider Key':'NodeSocketInt'})
    if aggregate_provider not in (0, 1): raise ValueError('Unknown aggregate provider transport revision')
    if aggregate_provider:
        if not (group.get('blosm_prepared_network_dispatch') == 1 and group.get('st_aggregate_collection_inputs') == 1):
            raise ValueError('Aggregate provider input requires canonical prepared-network dispatch')
        required.update({'Aggregate Result':'NodeSocketObject', 'Aggregate Provider Key':'NodeSocketInt',
                         'Prepared Network Kind':'NodeSocketInt'})
    sockets = {s.name:s.socket_type for s in group.interface.items_tree
               if s.item_type=='SOCKET' and s.in_out=='INPUT'}
    missing = [name for name,kind in required.items() if sockets.get(name)!=kind]
    if missing:
        raise ValueError('Original-span GN interface mismatch: '+', '.join(missing))


def validate_original_intake_provider(group, network=None):
    """Validate the optional native packet producer against its common entry."""
    if (group is None or group.bl_idname != 'GeometryNodeTree'
            or group.get('blosm_street_asset_version') != ORIGINAL_LIBRARY.revision
            or group.get('blosm_street_input_schema') != 'original_intake_provider.1'
            or group.get('blosm_original_span_provider') != 1):
        raise ValueError('Incompatible native original-span provider asset')
    sockets = lambda tree, direction: {s.name:s.socket_type for s in tree.interface.items_tree
        if s.item_type == 'SOCKET' and s.in_out == direction}
    inputs = sockets(group, 'INPUT')
    if (inputs.get('Geometry') != 'NodeSocketGeometry'
            or inputs.get('Original Span Provider Key') != 'NodeSocketInt'
            or sockets(group, 'OUTPUT').get('Geometry') != 'NodeSocketGeometry'):
        raise ValueError('Original-span provider packet interface is incomplete')
    if network is not None:
        if (network.get('blosm_original_provider_input') != 1
                or network.get('blosm_street_asset_version') != group.get('blosm_street_asset_version')):
            raise ValueError('Original-span provider and common entry revisions must match')
        public = sockets(network, 'INPUT')
        if public.get('Original Span Result') != 'NodeSocketObject':
            raise ValueError('The common entry has no original-span provider input')
        if any(public.get(name) != kind for name, kind in inputs.items()):
            raise ValueError('Original-span provider controls differ from the common entry')


def validate_original_separator_fill_provider(group):
    if (group is None or group.bl_idname!='GeometryNodeTree' or
            group.get('blosm_street_asset_version')!=ORIGINAL_LIBRARY.revision or
            group.get('blosm_street_input_schema')!='separator_fill_provider.1' or
            group.get('blosm_separator_fill_provider')!=1):
        raise ValueError('Incompatible native separator fill provider asset')
    required={name:'NodeSocketObject' for name in ('Mapping','Approaches','Junctions','Compounds',
        'Logical Clusters','Cluster Sites','Boundary Links','Cluster Spans','Corridors',
        'External Approaches','Separators','Separator Owners','Separator Ownership Result')}
    required.update({'Patch Results':'NodeSocketCollection','Strip Results':'NodeSocketCollection',
        'Grass':'NodeSocketMaterial','Cell Size':'NodeSocketFloat','Metres to Blender Units':'NodeSocketFloat',
        'Mesh Resolution':'NodeSocketInt','Separator Count':'NodeSocketInt',
        'Separator Ownership Key':'NodeSocketInt','Separator Fill Key':'NodeSocketInt',
        'Separator Fill Policy':'NodeSocketInt'})
    inputs={s.name:s.socket_type for s in group.interface.items_tree
            if s.item_type=='SOCKET' and s.in_out=='INPUT'}
    outputs={s.name:s.socket_type for s in group.interface.items_tree
             if s.item_type=='SOCKET' and s.in_out=='OUTPUT'}
    missing=[name for name,kind in required.items() if inputs.get(name)!=kind]
    if missing or outputs!={'Geometry':'NodeSocketGeometry'}:
        raise ValueError('Separator fill provider contract mismatch: '+', '.join(missing))


def validate_original_separator_ownership_provider(group):
    if (group is None or group.bl_idname != 'GeometryNodeTree' or
            group.get('blosm_street_asset_version') != ORIGINAL_LIBRARY.revision or
            group.get('blosm_street_input_schema') != 'separator_ownership_provider.1' or
            group.get('blosm_separator_ownership_provider') != 1):
        raise ValueError('Incompatible native separator ownership provider asset')
    sockets = {s.name:s.socket_type for s in group.interface.items_tree
               if s.item_type == 'SOCKET' and s.in_out == 'INPUT'}
    required = {name:'NodeSocketObject' for name in ('Mapping','Approaches','Junctions','Compounds',
        'Logical Clusters','Cluster Sites','Boundary Links','Cluster Spans','Corridors',
        'External Approaches','Separators','Separator Owners')}
    required.update({'Patch Results':'NodeSocketCollection', 'Separator Count':'NodeSocketInt',
                     'Separator Ownership Key':'NodeSocketInt', 'Metres to Blender Units':'NodeSocketFloat'})
    missing = [name for name, kind in required.items() if sockets.get(name) != kind]
    if missing:
        raise ValueError('Separator ownership provider input contract mismatch: '+', '.join(missing))


def validate_original_free_bundle_provider(group):
    if (group is None or group.bl_idname != 'GeometryNodeTree' or
            group.get('blosm_street_asset_version') != ORIGINAL_LIBRARY.revision or
            group.get('blosm_street_input_schema') != 'free_bundle_provider.1' or
            group.get('blosm_free_bundle_provider') != 1 or group.get('blosm_free_bundle_guards', 0) < 3):
        raise ValueError('Incompatible native free-Bundle provider asset')
    sockets = {s.name:s.socket_type for s in group.interface.items_tree
               if s.item_type == 'SOCKET' and s.in_out == 'INPUT'}
    required = {'Source Objects':'NodeSocketCollection', 'Source Binding':'NodeSocketObject',
        'Source Vertex Count':'NodeSocketInt', 'Source Edge Count':'NodeSocketInt',
        'Bundle Paths':'NodeSocketObject', 'Bundle Count':'NodeSocketInt',
        'Cell Size':'NodeSocketFloat', 'Mesh Resolution':'NodeSocketInt',
        'Height Exaggeration':'NodeSocketFloat', 'Metres to Blender Units':'NodeSocketFloat',
        'Minimum Separator Width':'NodeSocketFloat', 'Free Bundle Provider Key':'NodeSocketInt'}
    missing = [name for name, kind in required.items() if sockets.get(name) != kind]
    if missing:
        raise ValueError('Free-Bundle provider input contract mismatch: '+', '.join(missing))


def validate_original_aggregate_provider(group):
    """Native provider protocol; formulas and geometry are separately qualified."""
    if (group is None or group.bl_idname != 'GeometryNodeTree' or
            group.get('blosm_street_asset_version') != ORIGINAL_LIBRARY.revision or
            group.get('blosm_street_input_schema') != 'aggregate_provider.1' or
            group.get('blosm_aggregate_provider') != 1 or group.get('st_aggregate_collection_inputs') != 1):
        raise ValueError('Incompatible native aggregate provider asset')
    inputs = {s.name:s.socket_type for s in group.interface.items_tree
              if s.item_type == 'SOCKET' and s.in_out == 'INPUT'}
    required = {'Geometry':'NodeSocketGeometry', 'Source Objects':'NodeSocketCollection',
        'Source Vertex Count':'NodeSocketInt', 'Source Edge Count':'NodeSocketInt',
        'Aggregate Provider Key':'NodeSocketInt', 'Mesh Resolution':'NodeSocketInt',
        'Separator Count':'NodeSocketInt', 'Cell Size':'NodeSocketFloat',
        'Metres to Blender Units':'NodeSocketFloat', 'Height Exaggeration':'NodeSocketFloat'}
    required.update({key:'NodeSocketObject' for key in
        ('Mapping','Approaches','Junctions','Corridors','External Approaches','Transitions','Separators')})
    missing = [key for key,kind in required.items() if inputs.get(key) != kind]
    if missing: raise ValueError('Aggregate provider input contract mismatch: '+', '.join(missing))


def validate_original_patch_provider(group):
    _original_interpolation_revision(group)
    if (group.bl_idname!='GeometryNodeTree' or
            group.get('blosm_street_asset_version')!=ORIGINAL_LIBRARY.revision or
            group.get('blosm_street_input_schema')!='original_patch_provider.1'):
        raise ValueError('Incompatible original patch provider library')
    required={'Geometry':'NodeSocketGeometry','Source Objects':'NodeSocketCollection',
        'Source Vertex Count':'NodeSocketInt','Source Edge Count':'NodeSocketInt',
        'Mesh Resolution':'NodeSocketInt','Owner Index':'NodeSocketInt',
        'Cell Size':'NodeSocketFloat','Height Exaggeration':'NodeSocketFloat',
        'Metres to Blender Units':'NodeSocketFloat','Samples Per Side':'NodeSocketInt',
        'Relax Iterations':'NodeSocketInt'}
    required.update({key:'NodeSocketObject' for key in ('Source Binding','Source Order',
        'Mapping','Approaches','Junctions','Compounds','Logical Clusters','Cluster Sites',
        'Boundary Links','Cluster Spans')})
    internal_regions=group.get('blosm_internal_cycle_regions',0)
    if internal_regions not in (0,1):raise ValueError('Unknown internal-region provider capability revision')
    if internal_regions:required['Internal Cycles']='NodeSocketObject'
    sockets={s.name:s.socket_type for s in group.interface.items_tree
        if s.item_type=='SOCKET' and s.in_out=='INPUT'}
    missing=[key for key,kind in required.items() if sockets.get(key)!=kind]
    if missing:raise ValueError('Original patch provider interface mismatch: '+', '.join(missing))


def _original_interpolation_revision(group):
    revision = group.get('blosm_original_interpolation', 0)
    if revision not in (0, 1):
        raise ValueError('Unknown original-span interpolation capability revision')
    return revision


def validate_original_interpolation(group, patch_provider, strip_provider):
    """Exact source admission and every cached producer must agree.

    Existing unmarked assets remain Smooth-only. A single upgraded wrapper
    cannot enable Exact input while another producer still smooths its paths.
    Geometric applicability, including mode changes inside one feature run,
    remains the responsibility of the native constructors and their diagnostics.
    """
    flags = tuple(_original_interpolation_revision(block)
                  for block in (group, patch_provider, strip_provider))
    if flags not in ((0, 0, 0), (1, 1, 1)):
        raise ValueError('Incomplete original-span interpolation support')
    return flags == (1, 1, 1)


def validate_original_internal_regions(group, patch_provider):
    """Surface promotion and patch construction must share one capability."""
    flags=(group.get('blosm_internal_cycle_regions',0),
           patch_provider.get('blosm_internal_cycle_regions',0) if patch_provider is not None else 0)
    if flags not in ((0,0),(1,1)):
        raise ValueError('Incomplete original-span internal-region support')
    if flags==(0,0):return False
    if group.get('blosm_local_failures')!=1:
        raise ValueError('Internal regions require native local failure support')
    for block in (group,patch_provider):
        sockets={s.name:s.socket_type for s in block.interface.items_tree
                 if s.item_type=='SOCKET' and s.in_out=='INPUT'}
        if sockets.get('Internal Cycles')!='NodeSocketObject':
            raise ValueError('Internal-region asset lacks original cycle references')
    return True


def validate_original_internal_ends(group, patch_provider):
    """Do not admit retained internal connectors with an older patch solver."""
    capabilities=(group.get('blosm_internal_logical_ends',0),
                  patch_provider.get('blosm_internal_logical_ends',0))
    if capabilities not in ((0,0),(1,1)):
        raise ValueError('Incomplete original-span internal logical endpoint support')
    return capabilities==(1,1)


def validate_original_shared_owners(group, patch_provider):
    """Admit shared physical ownership only with matching chain-capable assets."""
    capabilities=tuple(g.get('blosm_shared_logical_owners',0) for g in (group,patch_provider))
    if capabilities not in ((0,0),(1,1)):
        raise ValueError('Incomplete original-span shared logical ownership support')
    if capabilities==(0,0):return False
    if not validate_original_internal_ends(group,patch_provider):
        raise ValueError('Shared logical owners require internal endpoint support')
    for g in (group,patch_provider):
        sockets={s.name:s.socket_type for s in g.interface.items_tree
                 if s.item_type=='SOCKET' and s.in_out=='INPUT'}
        if any(sockets.get(key)!='NodeSocketObject' for key in ('Chain Links','Chain Steps')):
            raise ValueError('Shared logical owner asset lacks original boundary program inputs')
    return True


def validate_original_seams(group, patch_provider, strip_provider):
    """Enable original two-ended connections only with a complete asset set.

    An older library may omit this capability. A partially upgraded library
    must reject before authoring inputs instead of emitting two free caps.
    """
    blocks=(group,patch_provider,strip_provider)
    flags=tuple(g.get('blosm_original_seams',0) for g in blocks)
    transitions=tuple(g.get('blosm_finite_seam_transitions',0) for g in blocks)
    if flags not in ((0,0,0),(1,1,1)) or transitions!=flags:
        raise ValueError('Incomplete original-span seam and finite-profile support')
    return flags==(1,1,1)


def validate_original_category_connections(group, patch_provider, strip_provider):
    """Cross-Street feature boundaries require the checked role-five branch."""
    flags=tuple(g.get('blosm_category_connections',0) for g in (group,patch_provider,strip_provider))
    if flags not in ((0,0,0),(1,1,1)):
        raise ValueError('Incomplete original cross-Street category handoff support')
    if flags==(0,0,0):return False
    if not group.get('blosm_section_scopes') or not validate_original_seams(group,patch_provider,strip_provider):
        raise ValueError('Category handoffs require Section scopes and original seams')
    return True


def validate_original_quad_energy(group, patch_provider):
    """The assembly and cached patches must use one numerical revision."""
    flags=tuple(g.get('blosm_stable_quad_energy',0) for g in (group,patch_provider))
    if flags not in ((0,0),(2,2)):
        raise ValueError('Incomplete original-span stable quad energy support')
    return flags==(2,2)


def validate_original_failure_support(group, strip_provider):
    if group.get('blosm_local_failures')!=1 or strip_provider.get('blosm_retains_local_failures')!=1:
        raise ValueError('Incomplete original-span local failure support')
    outputs={s.name:s.socket_type for s in group.interface.items_tree
             if s.item_type=='SOCKET' and s.in_out=='OUTPUT'}
    if any(outputs.get(name)!=kind for name,kind in (
            ('Valid','NodeSocketBool'),('Output Valid','NodeSocketBool'),('Failure Records','NodeSocketGeometry'))):
        raise ValueError('Missing original-span partial output contract')


def validate_original_rejected_profiles(group, patch_provider, strip_provider):
    """Unresolved designs need the same explicit contract on every consumer."""
    roots=(group,patch_provider,strip_provider)
    flags=tuple(g.get('blosm_rejected_profiles',0) if g is not None else 0 for g in roots)
    if flags==(0,0,0):return False
    if flags!=(1,1,1):
        raise ValueError('Incomplete rejected-profile support in the original-span asset')
    validate_original_failure_support(group,strip_provider)
    return True


def validate_original_failure_display(group):
    if (group.bl_idname!='GeometryNodeTree' or
            group.get('blosm_street_asset_version')!=ORIGINAL_LIBRARY.revision or
            group.get('blosm_street_input_schema')!='original_failure_display.1'):
        raise ValueError('Incompatible original failure display library')
    inputs={s.name:s.socket_type for s in group.interface.items_tree
            if s.item_type=='SOCKET' and s.in_out=='INPUT'}
    if any(inputs.get(name)!=kind for name,kind in (
            ('Network','NodeSocketObject'),('Material','NodeSocketMaterial'),('Metres to Blender Units','NodeSocketFloat'))):
        raise ValueError('Original failure display interface mismatch')


def validate_original_strip_provider(group):
    _original_interpolation_revision(group)
    if (group.bl_idname!='GeometryNodeTree' or
            group.get('blosm_street_asset_version')!=ORIGINAL_LIBRARY.revision or
            group.get('blosm_street_input_schema')!='original_strip_provider.1'):
        raise ValueError('Incompatible original strip provider library')
    required={'Geometry':'NodeSocketGeometry','Source Objects':'NodeSocketCollection',
        'Patch Results':'NodeSocketCollection','Source Vertex Count':'NodeSocketInt',
        'Source Edge Count':'NodeSocketInt','Mesh Resolution':'NodeSocketInt',
        'Owner Index':'NodeSocketInt','Owner Kind':'NodeSocketInt','Cell Size':'NodeSocketFloat',
        'Height Exaggeration':'NodeSocketFloat','Metres to Blender Units':'NodeSocketFloat',
        'Samples Per Side':'NodeSocketInt','Relax Iterations':'NodeSocketInt'}
    required.update({key:'NodeSocketObject' for key in ('Source Binding','Source Order',
        'Mapping','Approaches','Corridors','External Approaches','Transitions')})
    sockets={s.name:s.socket_type for s in group.interface.items_tree
        if s.item_type=='SOCKET' and s.in_out=='INPUT'}
    missing=[key for key,kind in required.items() if sockets.get(key)!=kind]
    if missing:raise ValueError('Original strip provider interface mismatch: '+', '.join(missing))


def validate_solver(group):
    if group.bl_idname != 'GeometryNodeTree' or group.get('blosm_street_asset_version') != LIBRARY.revision:
        raise ValueError('Incompatible street GN library revision; rebuild or select a compatible library')
    required = {'Source Objects':'NodeSocketCollection', 'Junction Inputs':'NodeSocketCollection',
        'Path Inputs':'NodeSocketObject', 'Source Vertex Count':'NodeSocketInt', 'Source Edge Count':'NodeSocketInt',
        'Cell Size':'NodeSocketFloat', 'Mesh Resolution':'NodeSocketInt', 'Height Exaggeration':'NodeSocketFloat',
        'Enabled':'NodeSocketBool', 'Show Diagnostics':'NodeSocketBool',
        'Pavement':'NodeSocketMaterial', 'Grass':'NodeSocketMaterial', 'Metres to Blender Units':'NodeSocketFloat',
        'Minimum Separator Width':'NodeSocketFloat'}
    sockets = {s.name:s.socket_type for s in group.interface.items_tree
               if s.item_type=='SOCKET' and s.in_out=='INPUT'}
    mismatches = [name for name,kind in required.items() if sockets.get(name)!=kind]
    if mismatches:
        raise ValueError('Street GN interface mismatch: '+', '.join(mismatches))
