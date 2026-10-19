"""Select shared original-span intake programs from final input declarations.

This optional catalog is data, not executable builder code. The universal
asset loader retains responsibility for editions, paths, revisions and reuse.
No GN trees, source geometry or update handlers are created here. A library
without this catalog keeps the existing complete-program loading behavior.
"""
import json
import time
from types import MappingProxyType

from . import assets
from ...app.asset_library import AssetRef


KINDS = ('complete', 'absent_optional_runs')
TABLES = ('Standalone Streets', 'Boundary Runs', 'Continuation Regions', 'Continuation Owners')
COUNTS = ('Standalone Count', 'Boundary Run Count', 'Continuation Count')


def declaration_kind(values):
    """Use the smaller program only for demonstrably empty final declarations.

    Inspect both stored table sizes and declared counts. Tables with modifiers
    may evaluate to different geometry, so conservatively use the full program.
    The native empty-declaration guard is still mandatory after selection.
    """
    missing = set(TABLES+COUNTS)-values.keys()
    if missing:
        raise ValueError('Missing original intake declarations: '+', '.join(sorted(missing)))
    for name in COUNTS:
        if type(values[name]) is not int or values[name] < 0:
            raise ValueError('Invalid original intake count: '+name)
    if any(values[name] for name in COUNTS):
        return 'complete'
    for name in TABLES:
        obj = values[name]
        if obj is None:
            continue  # Native Object Info on an unset input has empty geometry.
        if (getattr(obj, 'type', None) != 'MESH' or getattr(obj, 'modifiers', ())
                or getattr(obj, 'instance_type', 'NONE') != 'NONE'):
            return 'complete'
        vertices = getattr(getattr(obj, 'data', None), 'vertices', None)
        if vertices is None or len(vertices):
            return 'complete'
    return 'absent_optional_runs'


class OriginalIntakeCatalog:
    """Finite complete/empty-declaration entries in one versioned library."""
    def __init__(self, data, revision):
        if (not isinstance(data, dict) or set(data) != {'version', 'revision', 'programs'}
                or type(data['version']) is not int or data['version'] != 1
                or data['revision'] != revision):
            raise ValueError('Incompatible original intake catalog')
        programs = data['programs']
        if (not isinstance(programs, dict) or set(programs) != set(KINDS)
                or any(not isinstance(name, str) or not name.strip() for name in programs.values())
                or len(set(programs.values())) != len(KINDS)):
            raise ValueError('Require distinct complete and guarded empty-declaration entries')
        self.revision = revision
        self.programs = MappingProxyType(dict(programs))

    @classmethod
    def load_optional(cls, loader):
        ref = assets.ORIGINAL_INTAKE_CATALOG
        if not loader.contains(ref):
            return None
        return cls(json.loads(loader.load(ref).as_string()), assets.ORIGINAL_LIBRARY.revision)

    def reference(self, kind):
        if kind not in KINDS:
            raise ValueError('Unknown original intake program: '+str(kind))
        return AssetRef(assets.ORIGINAL_LIBRARY.key, 'node_groups', self.programs[kind])

    def preload_references(self, loader):
        # Retained linking reuses dependencies when loading the selected entry
        # later. APPEND must batch both entries with other roles to avoid child
        # duplication. That compatibility path is not a memory optimization.
        if loader.uses_link(assets.ORIGINAL_LIBRARY.key):
            return ()
        return tuple(self.reference(kind) for kind in KINDS)

    def validate_loaded(self, kind, tree, network):
        if self.revision != assets.ORIGINAL_LIBRARY.revision:
            raise ValueError('Original intake catalog revision changed')
        assets.validate_original_intake_provider(tree, network)
        empty = tree.get('st_absent_optional_runs', 0)
        if kind not in KINDS or type(empty) is not int or empty != int(kind == 'absent_optional_runs'):
            raise ValueError('Loaded original intake has a different specialization')
        if empty:
            constructors = [n.node_tree for n in tree.nodes if n.type == 'GROUP'
                            and n.label == 'Unchanged complete original-span construction']
            if len(constructors) != 1 or constructors[0].get('st_absent_optional_runs_guard') != 1:
                raise ValueError('Empty-declaration intake lacks its native applicability guard')


class OriginalIntakePrograms:
    """Resolve a shared entry after the controller's final tables are bound.

    Resolution is an authoring/topology-operation step. Ordinary coordinates
    and widths keep their existing native dependencies. If topology introduces
    optional declarations, its managed command must resolve again and rebind;
    the native guard rejects stale selection in the meantime.
    """
    def __init__(self, catalog, loader, loaded, *, on_load=None):
        self.catalog, self.loader = catalog, loader
        self.on_load = on_load
        self.loaded = {catalog.reference(kind): loaded[catalog.reference(kind)]
                       for kind in KINDS if catalog.reference(kind) in loaded}
        self.retained = loader.uses_link(assets.ORIGINAL_LIBRARY.key)

    @staticmethod
    def controller_values(controller):
        mod = controller.modifiers[0]
        return {socket.name: getattr(mod.properties.inputs, socket.identifier).value
                for socket in mod.node_group.interface.items_tree
                if socket.item_type == 'SOCKET' and socket.in_out == 'INPUT'
                and socket.socket_type != 'NodeSocketGeometry'}

    def resolve(self, controller):
        values = self.controller_values(controller)
        if type(values.get('Prepared Network Kind')) is not int or values['Prepared Network Kind'] != 0:
            raise ValueError('Original intake selection requires an original-span controller')
        kind = declaration_kind(values)
        ref = self.catalog.reference(kind)
        if ref not in self.loaded:
            if not self.retained:
                raise ValueError('Append-mode intake programs must be loaded in the common asset batch')
            started = time.perf_counter()
            self.loaded[ref] = self.loader.load(ref)
            if self.on_load is not None:
                self.on_load(time.perf_counter()-started)
        tree = self.loaded[ref]
        self.catalog.validate_loaded(kind, tree, controller.modifiers[0].node_group)
        return tree
