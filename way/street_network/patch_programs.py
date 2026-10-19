"""Shared native patch programs selected from prepared original topology.

Asset resolution remains in BlenderAssetLoader. This class owns only the street
program contract and routing; it never creates trees or evaluates geometry.
"""
from types import MappingProxyType
import json
import math
from ..gn.editable_input.patch_program_selection import owner_program_kinds, KINDS


def capabilities(tree):
    from ...app.asset_library import BlenderAssetLoader
    return {k:tree[k] for k in tree.keys() if k.startswith('blosm_')
            and k not in (BlenderAssetLoader.provenance_key,'blosm_builder_reference_manifest')}


class OriginalPatchCatalog:
    """Small library descriptor permits preparation before appending GN graphs.

    Descriptors contain data only. Resolution, editions, hashes and append/link
    reuse belong to the universal loader. Actual loaded roots must agree with
    the descriptor before their declared preparation can be used for authoring.
    """
    def __init__(self, data, revision):
        if (not isinstance(data,dict) or set(data)!={'version','revision','solver','patch','programs'}
                or type(data['version']) is not int or data['version']!=1 or data['revision']!=revision):
            raise ValueError('Incompatible street constructor catalog')
        for role in ('solver','patch'):
            block=data[role]
            if (not isinstance(block,dict) or any(not isinstance(k,str) or not k.startswith('blosm_')
                                                or type(v) not in (str,int,bool,float)
                                                or (type(v) is float and not math.isfinite(v))
                                                for k,v in block.items())
                    or block.get('blosm_street_asset_version')!=revision
                    or block.get('blosm_finite_owner_programs')!=1):
                raise ValueError('Invalid constructor catalog capabilities: '+role)
        programs=data['programs']
        if (not isinstance(programs,dict) or set(programs)!=set(KINDS)
                or any(not isinstance(name,str) or not name for name in programs.values())
                or len(set(programs.values()))!=len(KINDS)):
            raise ValueError('Require four distinct shared constructor entries')
        self.solver=MappingProxyType(dict(data['solver']))
        self.patch=MappingProxyType(dict(data['patch']))
        self.programs=MappingProxyType(dict(programs))

    @classmethod
    def load_optional(cls, loader):
        from . import assets
        ref=assets.ORIGINAL_PATCH_CATALOG
        if not loader.contains(ref):return None
        return cls(json.loads(loader.load(ref).as_string()),assets.ORIGINAL_LIBRARY.revision)

    def references(self, kinds):
        from ...app.asset_library import AssetRef
        from . import assets
        kinds=set(kinds)
        if not kinds<=set(KINDS):raise ValueError('Unknown required constructor kind')
        return {kind:AssetRef(assets.ORIGINAL_LIBRARY.key,'node_groups',self.programs[kind])
                for kind in KINDS if kind in kinds}

    def validate_solver(self, solver):
        """Check the common entry even when no original patch is requested."""
        if capabilities(solver)!=self.solver:
            raise ValueError('Loaded network capabilities differ from the constructor catalog')

    def validate_loaded(self, solver, programs):
        self.validate_solver(solver)
        if any(capabilities(tree)!=self.patch for tree in programs.trees.values()):
            raise ValueError('Loaded constructor capabilities differ from the catalog')

    def prepare(self, inventory, categories, *, category_properties=None, **selection):
        from .original_span_plan import OriginalSpanPlanner
        s,p=self.solver,self.patch
        return OriginalSpanPlanner(categories,category_properties=category_properties,
            defer_physical_categories=bool(s.get('blosm_filtered_physical_ports')),
            section_scopes=bool(s.get('blosm_section_scopes')),
            internal_logical_ends=bool(s.get('blosm_internal_logical_ends') and p.get('blosm_internal_logical_ends')),
            shared_logical_owners=bool(s.get('blosm_shared_logical_owners') and p.get('blosm_shared_logical_owners')),
            seam_connections=bool(s.get('blosm_original_seams')),
            category_connections=bool(s.get('blosm_category_connections')),
            boundary_runs=bool(s.get('blosm_boundary_runs')),
            dimension_runs=bool(s.get('blosm_authored_dimension_runs')),
            free_bundles=bool(s.get('blosm_free_bundle_surfaces'))).prepare(
                inventory,standalone=bool(s.get('blosm_standalone_streets')),**selection)


class OriginalPatchPrograms:
    @staticmethod
    def for_component(component):
        plan, network = component.ports, component.network
        index = {site.junction_id:i for i, site in enumerate(plan.sites)}
        if len(index) != len(plan.sites):
            raise ValueError('Repeated original site identity')
        return owner_program_kinds(len(plan.sites),
            (tuple(index[j] for j in pair.site_ids) for pair in network.compounds),
            (tuple(index[j] for j in cluster.physical_junction_ids) for cluster in network.logical_clusters),
            (i for i, site in enumerate(plan.sites)
             if len(site.ports) == 2 and site.ports[0].source_kind == 'SEAM'))

    def __init__(self, programs):
        from .assets import validate_original_patch_provider
        self.trees = MappingProxyType(dict(programs))
        if not self.trees or not set(self.trees) <= set(KINDS):
            raise ValueError('Require known finite shared patch programs')
        if len(set(self.trees.values())) != len(self.trees):
            raise ValueError('Each guarded kind requires its declared shared tree')
        self.representative = next(iter(self.trees.values()))
        def contract(tree):
            return [(s.identifier,s.name,s.in_out,s.socket_type) for s in tree.interface.items_tree
                    if s.item_type == 'SOCKET']
        for kind, tree in self.trees.items():
            validate_original_patch_provider(tree)
            family = 'logical' if kind == 'logical' else 'ordinary'
            if (tree.get('st_shared_constructor_family') != family or
                    tree.get('st_constructor_family_guarded') != 1 or
                    (kind != 'logical' and (tree.get('st_shared_owner_kind') != kind or
                                            tree.get('st_owner_kind_guarded') != 1))):
                raise ValueError('Missing native constructor applicability guard: '+kind)
            if (contract(tree) != contract(self.representative) or
                    capabilities(tree) != capabilities(self.representative)):
                raise ValueError('Shared constructor contracts or capabilities differ')

    def validate_component(self, component):
        missing = set(self.for_component(component).values()) - set(self.trees)
        if missing:
            raise ValueError('Missing required prepared constructor programs: '+', '.join(sorted(missing)))

    def bind(self, controller, *, network_tree, collection):
        from ..gn.editable_input.curved_patch_provider_bindings import install
        return install(controller, provider_programs=self.trees, network_tree=network_tree, collection=collection)
