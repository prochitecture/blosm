"""Host checks of finite program loading policy and original topology routing."""
from copy import deepcopy
import importlib
import json
from pathlib import Path
import sys
from types import ModuleType,SimpleNamespace
import unittest

ROOT=Path(__file__).resolve().parents[1]
PREFIX='_patch_program_tests'
for suffix,relative in (('', ''),('.app','app'),('.way','way'),('.way.street_network','way/street_network'),
                        ('.way.gn','way/gn'),('.way.gn.editable_input','way/gn/editable_input')):
    package=ModuleType(PREFIX+suffix);package.__path__=[str(ROOT/relative)];sys.modules[package.__name__]=package
m=importlib.import_module(PREFIX+'.way.street_network.patch_programs')
assets=importlib.import_module(PREFIX+'.way.street_network.assets')


def descriptor():
    common=dict(blosm_street_asset_version=assets.ORIGINAL_LIBRARY.revision,blosm_finite_owner_programs=1,
        blosm_filtered_physical_ports=1,blosm_section_scopes=1,blosm_internal_logical_ends=1,
        blosm_shared_logical_owners=1,blosm_original_seams=1,blosm_category_connections=1,
        blosm_boundary_runs=1,blosm_authored_dimension_runs=1,blosm_free_bundle_surfaces=1,
        blosm_standalone_streets=1)
    return dict(version=1,revision=assets.ORIGINAL_LIBRARY.revision,solver=common.copy(),patch=common.copy(),
                programs={kind:'Shared constructor / '+kind for kind in m.KINDS})


class Tree(dict):
    __hash__=object.__hash__
    bl_idname='GeometryNodeTree'


def tree(kind):
    block=Tree(blosm_street_asset_version=assets.ORIGINAL_LIBRARY.revision,
        blosm_street_input_schema='original_patch_provider.1',blosm_finite_owner_programs=1,
        st_shared_constructor_family='logical' if kind=='logical' else 'ordinary',st_constructor_family_guarded=1)
    if kind!='logical':block.update(st_shared_owner_kind=kind,st_owner_kind_guarded=1)
    inputs={'Geometry':'NodeSocketGeometry','Source Objects':'NodeSocketCollection'}
    inputs.update({name:'NodeSocketInt' for name in ('Source Vertex Count','Source Edge Count','Mesh Resolution',
        'Owner Index','Samples Per Side','Relax Iterations')})
    inputs.update({name:'NodeSocketFloat' for name in ('Cell Size','Height Exaggeration','Metres to Blender Units')})
    inputs.update({name:'NodeSocketObject' for name in ('Source Binding','Source Order','Mapping','Approaches',
        'Junctions','Compounds','Logical Clusters','Cluster Sites','Boundary Links','Cluster Spans')})
    block.interface=SimpleNamespace(items_tree=[SimpleNamespace(name=name,socket_type=kind,identifier=str(i),
        item_type='SOCKET',in_out='INPUT') for i,(name,kind) in enumerate(inputs.items())])
    return block


class PatchProgramCatalogTests(unittest.TestCase):
    def catalog(self):return m.OriginalPatchCatalog(descriptor(),assets.ORIGINAL_LIBRARY.revision)

    def test_batch_contains_only_required_kinds_in_same_library(self):
        refs=self.catalog().references(('paired','single','single'))
        self.assertEqual(tuple(refs),('single','paired'))
        self.assertTrue(all(ref.library==assets.ORIGINAL_LIBRARY.key for ref in refs.values()))
        self.assertEqual(self.catalog().references(()),{})

    def test_unknown_kind_and_incomplete_descriptor_reject(self):
        with self.assertRaises(ValueError):self.catalog().references(('unknown',))
        for change in (lambda d:d.update(version=True),lambda d:d.update(revision='wrong'),
                       lambda d:d['programs'].pop('logical'),lambda d:d['programs'].update(seam=d['programs']['single']),
                       lambda d:d['solver'].update(blosm_finite_owner_programs=0),
                       lambda d:d['patch'].update(blosm_external_code=['not data'])):
            data=descriptor();change(data)
            with self.assertRaises(ValueError):m.OriginalPatchCatalog(data,assets.ORIGINAL_LIBRARY.revision)

    def test_optional_descriptor_loads_only_text_before_graph_batch(self):
        calls=[]
        class Loader:
            def contains(self,ref):calls.append(('contains',ref));return True
            def load(self,ref):calls.append(('load',ref));return SimpleNamespace(as_string=lambda:json.dumps(descriptor()))
        self.assertIsNotNone(m.OriginalPatchCatalog.load_optional(Loader()))
        self.assertEqual([ref.kind for _,ref in calls],['texts','texts'])
        absent=SimpleNamespace(contains=lambda ref:False,load=lambda ref:self.fail('Absent descriptor loaded'))
        self.assertIsNone(m.OriginalPatchCatalog.load_optional(absent))

    def test_programs_retain_identity_and_ignore_loader_provenance(self):
        single,paired=tree('single'),tree('paired')
        single['blosm_asset_reference']='one';paired['blosm_asset_reference']='two'
        programs=m.OriginalPatchPrograms({'single':single,'paired':paired})
        self.assertIs(programs.trees['single'],single)
        self.assertIs(programs.trees['paired'],paired)
        with self.assertRaises(TypeError):programs.trees['seam']=tree('seam')

    def test_native_kind_guard_and_contract_cannot_be_substituted(self):
        single,paired=tree('single'),tree('paired')
        with self.assertRaises(ValueError):m.OriginalPatchPrograms({'single':single,'paired':single})
        paired['st_owner_kind_guarded']=0
        with self.assertRaises(ValueError):m.OriginalPatchPrograms({'paired':paired})
        paired=tree('paired');paired.interface.items_tree[0].identifier='changed'
        with self.assertRaises(ValueError):m.OriginalPatchPrograms({'single':single,'paired':paired})
        paired=tree('paired');paired['blosm_original_interpolation']=1
        with self.assertRaises(ValueError):m.OriginalPatchPrograms({'single':single,'paired':paired})

    def test_loaded_descriptor_is_checked_before_authoring(self):
        single=tree('single');programs=m.OriginalPatchPrograms({'single':single})
        data=descriptor();data['patch']=m.capabilities(single)
        catalog=m.OriginalPatchCatalog(data,assets.ORIGINAL_LIBRARY.revision)
        catalog.validate_loaded(data['solver'],programs)
        changed=dict(data['solver'],blosm_original_seams=0)
        with self.assertRaises(ValueError):catalog.validate_loaded(changed,programs)
        single['blosm_boundary_runs']=1
        with self.assertRaises(ValueError):catalog.validate_loaded(data['solver'],programs)

    def test_aggregate_only_entry_is_checked_without_loading_original_programs(self):
        catalog=self.catalog();data=descriptor()
        catalog.validate_solver(data['solver'])
        for change in ({'blosm_street_asset_version':'wrong'}, {'blosm_original_seams':0}):
            with self.subTest(change=change),self.assertRaises(ValueError):
                catalog.validate_solver(dict(data['solver'],**change))

    def test_finite_decimal_construction_settings_round_trip_and_match(self):
        single=tree('single')
        single.update(blosm_corner_radial_edge_budget=7.5,blosm_corner_refinement_limit=7.9,
                      blosm_corner_guide_blend=.25,blosm_corner_minimum_profile_depth=.75,
                      blosm_corner_profile_slope_reserve=.02)
        programs=m.OriginalPatchPrograms({'single':single})
        data=descriptor();data['patch']=m.capabilities(single)
        data['solver'].update({k:v for k,v in single.items() if type(v) is float})
        catalog=m.OriginalPatchCatalog(json.loads(json.dumps(data)),assets.ORIGINAL_LIBRARY.revision)
        catalog.validate_loaded(data['solver'],programs)
        changed=dict(data['solver'],blosm_corner_guide_blend=.3)
        with self.assertRaises(ValueError):catalog.validate_loaded(changed,programs)
        single['blosm_corner_minimum_profile_depth']=.5
        with self.assertRaises(ValueError):catalog.validate_loaded(data['solver'],programs)

    def test_nonfinite_decimal_settings_are_rejected_for_both_roles(self):
        for role in ('solver','patch'):
            for value in (float('nan'),float('inf'),float('-inf')):
                data=descriptor();data[role]['blosm_corner_refinement_limit']=value
                with self.subTest(role=role,value=value),self.assertRaises(ValueError):
                    m.OriginalPatchCatalog(data,assets.ORIGINAL_LIBRARY.revision)

    def test_original_seam_preparation_and_selection_survive(self):
        from test_physical_ports import seam_fixture
        source=seam_fixture();before=deepcopy(source)
        result=self.catalog().prepare(source,('residential',))
        self.assertFalse(result.issues);component,=result.components
        kinds=m.OriginalPatchPrograms.for_component(component)
        seam=next(i for i,s in enumerate(component.ports.sites) if s.junction_id==2000)
        self.assertEqual(kinds[seam],'seam');self.assertEqual(set(kinds.values()),{'single','seam'})
        self.assertEqual(source,before)
        with self.assertRaises(ValueError):m.OriginalPatchPrograms({'single':tree('single')}).validate_component(component)

    def test_shared_logical_cluster_has_one_program_and_original_members(self):
        from test_shared_port_clusters import fixture
        source=fixture();before=deepcopy(source)
        result=self.catalog().prepare(source,('residential',))
        self.assertFalse(result.issues);component,=result.components
        self.assertEqual(m.OriginalPatchPrograms.for_component(component),{0:'logical'})
        self.assertEqual(component.network.logical_clusters[0].owner_ids,(5000,5001))
        self.assertEqual(source,before)

    def test_category_handoff_does_not_become_a_road_seam(self):
        from test_physical_ports import seam_fixture
        source=seam_fixture();source.streets[8]['items'][0]['category']='footway'
        source.streets[8].update(status='UNSUPPORTED',reason='Reserved component/category: footway')
        result=self.catalog().prepare(source,('residential',))
        self.assertFalse(result.issues);component,=result.components
        self.assertNotIn(2000,{s.junction_id for s in component.ports.sites})
        self.assertEqual(next(r.kind for r in component.network.routes if r.span.street_id==4),'CATEGORY_APPROACH')

    def test_optional_bundle_batch_includes_separator_ownership(self):
        seen=[]
        loader=SimpleNamespace(uses_link=lambda key:False,contains=lambda ref:seen.append(ref) or True)
        expected=(assets.ORIGINAL_INTAKE_PROVIDER,assets.ORIGINAL_AGGREGATE_PROVIDER,assets.ORIGINAL_FREE_BUNDLE_PROVIDER,
                  assets.ORIGINAL_SEPARATOR_OWNERSHIP_PROVIDER,assets.ORIGINAL_SEPARATOR_FILL_PROVIDER)
        self.assertEqual(assets.optional_original_programs(loader,bundles=True,aggregate=True),expected)
        self.assertEqual(tuple(seen),expected)
        self.assertEqual(assets.optional_original_programs(loader),(assets.ORIGINAL_INTAKE_PROVIDER,))
        seen.clear()
        self.assertEqual(assets.optional_original_programs(loader,original=False,aggregate=True),
                         (assets.ORIGINAL_AGGREGATE_PROVIDER,))
        self.assertEqual(seen,[assets.ORIGINAL_AGGREGATE_PROVIDER])

    def test_retained_consumers_can_defer_bundles_but_appends_cannot(self):
        loader=SimpleNamespace(uses_link=lambda key:True,contains=lambda ref:True)
        self.assertEqual(assets.optional_original_programs(loader,bundles=True,aggregate=True,defer_linked=True),
                         (assets.ORIGINAL_INTAKE_PROVIDER,assets.ORIGINAL_AGGREGATE_PROVIDER))
        self.assertEqual(assets.optional_original_programs(loader,bundles=True),
                         (assets.ORIGINAL_INTAKE_PROVIDER,assets.ORIGINAL_FREE_BUNDLE_PROVIDER,assets.ORIGINAL_SEPARATOR_OWNERSHIP_PROVIDER,
                          assets.ORIGINAL_SEPARATOR_FILL_PROVIDER))
        loader.uses_link=lambda key:False
        self.assertEqual(assets.optional_original_programs(loader,bundles=True,defer_linked=True),
                         (assets.ORIGINAL_INTAKE_PROVIDER,assets.ORIGINAL_FREE_BUNDLE_PROVIDER,assets.ORIGINAL_SEPARATOR_OWNERSHIP_PROVIDER,
                          assets.ORIGINAL_SEPARATOR_FILL_PROVIDER))

    def test_legacy_library_without_optional_fill_keeps_its_original_roles(self):
        loader=SimpleNamespace(uses_link=lambda key:False,
            contains=lambda ref:ref not in (assets.ORIGINAL_SEPARATOR_FILL_PROVIDER,assets.ORIGINAL_INTAKE_PROVIDER))
        self.assertEqual(assets.optional_original_programs(loader,bundles=True),
                         (assets.ORIGINAL_FREE_BUNDLE_PROVIDER,assets.ORIGINAL_SEPARATOR_OWNERSHIP_PROVIDER))

    def test_optional_assets_do_not_suppress_edition_or_path_errors(self):
        loader=SimpleNamespace(uses_link=lambda key:False,
            contains=lambda ref:ref==assets.ORIGINAL_SEPARATOR_OWNERSHIP_PROVIDER)
        self.assertEqual(assets.optional_original_programs(loader,bundles=True),
                         (assets.ORIGINAL_SEPARATOR_OWNERSHIP_PROVIDER,))
        def unavailable(ref):raise RuntimeError('Edition or path unavailable')
        loader.contains=unavailable
        with self.assertRaisesRegex(RuntimeError,'Edition or path'):
            assets.optional_original_programs(loader,bundles=True)


if __name__=='__main__':unittest.main()
