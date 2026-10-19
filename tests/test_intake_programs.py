"""Final-declaration dispatch and shared loading without a Blender process."""
import importlib
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace as NS
import unittest


ROOT = Path(__file__).resolve().parents[1]
PREFIX = '_intake_program_tests'
for suffix, relative in (('', ''), ('.app', 'app'), ('.way', 'way'),
                         ('.way.street_network', 'way/street_network')):
    package = ModuleType(PREFIX+suffix)
    package.__path__ = [str(ROOT/relative)]
    sys.modules[package.__name__] = package
m = importlib.import_module(PREFIX+'.way.street_network.intake_programs')
assets = m.assets


def descriptor():
    return dict(version=1, revision=assets.ORIGINAL_LIBRARY.revision,
                programs=dict(complete='Shared complete intake', absent_optional_runs='Shared empty intake'))


def values():
    result = dict.fromkeys(m.TABLES, None)
    result.update(dict.fromkeys(m.COUNTS, 0))
    result['Prepared Network Kind'] = 0
    result['Original Span Provider Key'] = 1
    return result


def table(count=0, **kwargs):
    return NS(type='MESH', data=NS(vertices=[None]*count), modifiers=(), instance_type='NONE', **kwargs)


class Tree(dict):
    bl_idname = 'GeometryNodeTree'


def tree(empty=False):
    result = Tree(blosm_street_asset_version=assets.ORIGINAL_LIBRARY.revision,
        blosm_street_input_schema='original_intake_provider.1', blosm_original_span_provider=1)
    if empty: result['st_absent_optional_runs'] = 1
    kinds = {'Geometry': 'NodeSocketGeometry'}
    kinds.update(dict.fromkeys(m.TABLES, 'NodeSocketObject'))
    kinds.update(dict.fromkeys(m.COUNTS+('Prepared Network Kind', 'Original Span Provider Key'), 'NodeSocketInt'))
    result.interface = NS(items_tree=[NS(name=name, identifier='Socket_'+str(i),
        socket_type=kind, item_type='SOCKET', in_out='INPUT') for i, (name, kind) in enumerate(kinds.items())])
    result.interface.items_tree.append(NS(name='Geometry', identifier='Output_0',
        socket_type='NodeSocketGeometry', item_type='SOCKET', in_out='OUTPUT'))
    result.nodes = [NS(type='GROUP', label='Unchanged complete original-span construction',
                       node_tree={'st_absent_optional_runs_guard': 1} if empty else {})]
    return result


def controller():
    network = tree()
    network['blosm_original_provider_input'] = 1
    network.interface.items_tree.append(NS(name='Original Span Result', identifier='Result',
        socket_type='NodeSocketObject', item_type='SOCKET', in_out='INPUT'))
    props = NS()
    state = values()
    for socket in network.interface.items_tree:
        if socket.in_out == 'INPUT' and socket.socket_type != 'NodeSocketGeometry':
            setattr(props, socket.identifier, NS(value=state.get(socket.name)))
    return NS(modifiers=[NS(node_group=network, properties=NS(inputs=props))])


def change(obj, name, value):
    mod = obj.modifiers[0]
    socket, = [s for s in mod.node_group.interface.items_tree if s.in_out == 'INPUT' and s.name == name]
    getattr(mod.properties.inputs, socket.identifier).value = value


class Loader:
    def __init__(self, catalog, retained=True):
        self.retained = retained
        self.blocks = {catalog.reference('complete'): tree(), catalog.reference('absent_optional_runs'): tree(True)}
        self.calls = []

    def uses_link(self, key):
        return self.retained

    def load(self, ref):
        self.calls.append(ref)
        return self.blocks[ref]


class IntakeProgramsTests(unittest.TestCase):
    def catalog(self):
        return m.OriginalIntakeCatalog(descriptor(), assets.ORIGINAL_LIBRARY.revision)

    def test_empty_objects_and_unset_objects_use_same_empty_policy(self):
        self.assertEqual(m.declaration_kind(values()), 'absent_optional_runs')
        state = values()
        state.update({name: table() for name in m.TABLES})
        self.assertEqual(m.declaration_kind(state), 'absent_optional_runs')

    def test_every_actual_optional_table_can_require_complete_program(self):
        for name in m.TABLES:
            state = values(); state[name] = table(1)
            with self.subTest(name=name):
                self.assertEqual(m.declaration_kind(state), 'complete')

    def test_declared_nonzero_count_cannot_select_empty_table_variant(self):
        for name in m.COUNTS:
            state = values(); state[name] = 1
            self.assertEqual(m.declaration_kind(state), 'complete')

    def test_missing_or_invalid_counts_reject(self):
        for name in m.TABLES+m.COUNTS:
            state = values(); del state[name]
            with self.subTest(name=name), self.assertRaises(ValueError): m.declaration_kind(state)
        for value in (-1, .0, True, float('nan'), float('inf')):
            state = values(); state[m.COUNTS[0]] = value
            with self.subTest(value=value), self.assertRaises(ValueError): m.declaration_kind(state)

    def test_modified_instanced_and_nonmesh_tables_cannot_be_assumed_empty(self):
        for attribute, value in (('modifiers', [object()]), ('type', 'CURVE'), ('instance_type', 'COLLECTION')):
            obj = table(); setattr(obj, attribute, value)
            state = values(); state[m.TABLES[0]] = obj
            self.assertEqual(m.declaration_kind(state), 'complete')

    def test_catalog_is_data_only_and_revision_pinned(self):
        for edit in (lambda d: d.update(version=True), lambda d: d.update(revision='wrong'),
                     lambda d: d.update(code='run this'), lambda d: d['programs'].pop('complete'),
                     lambda d: d['programs'].update(complete=d['programs']['absent_optional_runs']),
                     lambda d: d['programs'].update(complete=' ')):
            data = descriptor(); edit(data)
            with self.assertRaises(ValueError): m.OriginalIntakeCatalog(data, assets.ORIGINAL_LIBRARY.revision)
        catalog = self.catalog()
        with self.assertRaises(TypeError): catalog.programs['complete'] = 'different'
        with self.assertRaises(ValueError): catalog.reference('unknown')

    def test_missing_catalog_keeps_legacy_path_without_loading_graphs(self):
        loader = NS(contains=lambda ref: False, load=lambda ref: self.fail('Absent catalog loaded'))
        self.assertIsNone(m.OriginalIntakeCatalog.load_optional(loader))

    def test_catalog_lookup_reads_only_text_and_preserves_access_errors(self):
        seen = []
        def load(ref):
            seen.append(ref)
            return NS(as_string=lambda: json.dumps(descriptor()))
        loader = NS(contains=lambda ref: True, load=load)
        self.assertIsNotNone(m.OriginalIntakeCatalog.load_optional(loader))
        self.assertEqual(seen, [assets.ORIGINAL_INTAKE_CATALOG])
        def unavailable(ref): raise RuntimeError('Edition unavailable')
        loader.contains = unavailable
        with self.assertRaisesRegex(RuntimeError, 'Edition unavailable'):
            m.OriginalIntakeCatalog.load_optional(loader)

    def test_linked_loading_is_lazy_shared_and_changes_with_final_tables(self):
        catalog = self.catalog(); loader = Loader(catalog)
        self.assertEqual(catalog.preload_references(loader), ())
        timings = []
        programs = m.OriginalIntakePrograms(catalog, loader, {}, on_load=timings.append)
        obj = controller()
        empty = programs.resolve(obj)
        self.assertEqual(loader.calls, [catalog.reference('absent_optional_runs')])
        self.assertIs(programs.resolve(controller()), empty)
        # Same controller now owns a nonempty optional table. Do not reuse the
        # old specialization merely because its original counts were zero.
        change(obj, 'Continuation Owners', table(1))
        self.assertIs(programs.resolve(obj), loader.blocks[catalog.reference('complete')])
        change(obj, 'Continuation Owners', None)
        self.assertIs(programs.resolve(obj), empty)
        self.assertEqual(len(loader.calls), 2)
        self.assertEqual(len(timings), 2)
        self.assertTrue(all(seconds >= 0 for seconds in timings))

    def test_append_loads_common_batch_and_forbids_late_dependency_duplication(self):
        catalog = self.catalog(); loader = Loader(catalog, retained=False)
        self.assertEqual(catalog.preload_references(loader), tuple(catalog.reference(k) for k in m.KINDS))
        programs = m.OriginalIntakePrograms(catalog, loader, loader.blocks)
        self.assertIs(programs.resolve(controller()), loader.blocks[catalog.reference('absent_optional_runs')])
        self.assertFalse(loader.calls)
        with self.assertRaisesRegex(ValueError, 'common asset batch'):
            m.OriginalIntakePrograms(catalog, loader, {}).resolve(controller())
        self.assertFalse(loader.calls)

    def test_actual_schema_and_native_guard_are_checked_before_binding(self):
        for mutate in (lambda t: t.pop('st_absent_optional_runs'),
                       lambda t: t.nodes.clear(),
                       lambda t: t.nodes[0].node_tree.clear(),
                       lambda t: t.update(blosm_street_asset_version='wrong'),
                       lambda t: t.interface.items_tree.pop(0)):
            catalog = self.catalog(); loader = Loader(catalog)
            mutate(loader.blocks[catalog.reference('absent_optional_runs')])
            with self.assertRaises(ValueError): m.OriginalIntakePrograms(catalog, loader, {}).resolve(controller())

    def test_complete_reference_cannot_silently_supply_small_program(self):
        catalog = self.catalog(); loader = Loader(catalog)
        loader.blocks[catalog.reference('complete')] = tree(True)
        obj = controller(); change(obj, 'Standalone Count', 1)
        with self.assertRaises(ValueError): m.OriginalIntakePrograms(catalog, loader, {}).resolve(obj)

    def test_aggregate_inputs_are_not_original_intake(self):
        for kind in (1, True, None):
            catalog = self.catalog(); loader = Loader(catalog)
            obj = controller(); change(obj, 'Prepared Network Kind', kind)
            with self.assertRaises(ValueError): m.OriginalIntakePrograms(catalog, loader, {}).resolve(obj)
            self.assertFalse(loader.calls)


if __name__ == '__main__': unittest.main()
