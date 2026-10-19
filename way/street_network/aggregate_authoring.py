"""Author clean aggregate source polylines with one supplied shared GN asset.

Preparation supplies complete Street/Bundle ownership. This class only authors
source data and relation tables; it neither builds per-network solvers nor
computes surfaces, and requires no handler for coordinate/parameter edits.
"""
from dataclasses import asdict, dataclass
import json
import bpy
from ..gn.editable_input import aggregate_port_bindings, curved_port_bindings, curved_separator_bindings
from ..gn.kernel.pipeline import set_input
from .aggregate_coverage import AggregateSurfaceCoverage


@dataclass(frozen=True)
class AggregateAuthoredNetwork:
    controller: object
    tables: dict
    relation_collection: object
    source_objects: tuple = ()


class AggregateNetworkBuilder:
    def __init__(self, solver, *, pavement=None, grass=None, provider=None):
        if solver.get('st_aggregate_modifier_version') != 1:
            raise ValueError('A compatible shared aggregate modifier asset is required')
        self.solver = solver; self.pavement = pavement; self.grass = grass
        self.provider = provider
        if solver.get('blosm_aggregate_provider_input'):
            from .assets import validate_original_aggregate_provider
            validate_original_aggregate_provider(provider)
            if provider.get('blosm_street_asset_version') != solver.get('blosm_street_asset_version'):
                raise ValueError('Aggregate provider and entry revisions must match')

    def create(self, name, snapshot, plan, collection, *, section_modes,
               cell_size=.7, resolution=1, exaggeration=10., origin=(0., 0., 0.), source_groups=None):
        if plan.separators.pending:
            raise ValueError('Aggregate network has unresolved separator declarations')
        compiled, approaches = aggregate_port_bindings.compile(snapshot, plan.ports, name, origin=origin)
        modes = [section_modes[d.value] for d in compiled.sources.attributes['st_section_id'].data]
        compiled.sources.attributes.new('st_interpolation', 'INT', 'EDGE').data.foreach_set('value', modes)
        dimensions = curved_port_bindings.dimension_transitions(plan.ports, compiled, plan.dimensions, name)
        tables = curved_port_bindings.network_tables(plan.ports, plan.network, compiled, approaches, name, dimensions=dimensions)
        tables['Separators'] = curved_separator_bindings.bind(plan.separators, tables, name,
            free_ends=True, direct_sites=True, selection_ends=True)
        collection_inputs = self.solver.get('st_aggregate_collection_inputs') == 1
        if collection_inputs:
            from .authoring import metre_scale
            scale = metre_scale()
            for vertex in compiled.sources.vertices: vertex.co /= scale
            for key, domain, count in (('st_live_vertex_index', 'POINT', len(compiled.sources.vertices)),
                                       ('st_live_edge_index', 'EDGE', len(compiled.sources.edges))):
                compiled.sources.attributes.new(key, 'INT', domain).data.foreach_set('value', list(range(count)))
        else:
            scale = 1.
            if source_groups is not None: raise ValueError('Source groups need a collection-capable shared asset')
        private = bpy.data.collections.new(name+' / private relations'); collection.children.link(private)
        objects = {}
        for key, mesh in tables.items():
            if collection_inputs and key == 'Sources': mesh = bpy.data.meshes.new(name+' controller')
            obj = bpy.data.objects.new(name if key == 'Sources' else name+' / '+key, mesh)
            (collection if key == 'Sources' else private).objects.link(obj)
            obj.location = tuple(v/scale for v in origin); objects[key] = obj
            if key != 'Sources':
                obj.hide_render = True; obj.hide_set(True); obj.hide_select = True
        controller = objects['Sources']
        mod = controller.modifiers.new('Shared aggregate street network', 'NODES'); mod.node_group = self.solver
        mod.show_viewport = False
        mod.show_in_editmode = True
        for socket in self.solver.interface.items_tree:
            if socket.item_type == 'SOCKET' and socket.in_out == 'INPUT' and socket.socket_type == 'NodeSocketObject':
                if socket.name in objects: set_input(controller, socket.name, objects[socket.name])
        for key, value in (('Cell Size', cell_size), ('Mesh Resolution', resolution),
                           ('Height Exaggeration', exaggeration), ('Separator Count', len(tables['Separators'].vertices)),
                           ('Pavement', self.pavement), ('Grass', self.grass)):
            set_input(controller, key, value)
        sources = ()
        if collection_inputs:
            editable = bpy.data.collections.new(name+' / editable Streets'); collection.children.link(editable)
            groups = source_groups if source_groups is not None else tuple((name+' / Street '+str(s), (s,)) for s in plan.ports.source_streets)
            sources = self._sources(compiled.sources, groups, controller, editable)
            for key, value in (('Source Objects', editable), ('Source Vertex Count', len(compiled.sources.vertices)),
                               ('Source Edge Count', len(compiled.sources.edges)), ('Metres to Blender Units', 1./scale)):
                set_input(controller, key, value)
            bpy.data.meshes.remove(compiled.sources)
        if self.solver.get('blosm_prepared_network_dispatch'):
            set_input(controller, 'Prepared Network Kind', 1)
        for mesh in (compiled.ports, approaches):
            if mesh.users == 0: bpy.data.meshes.remove(mesh)
        controller['blosm_aggregate_network_version'] = 1
        controller['blosm_aggregate_streets'] = list(plan.ports.source_streets)
        controller['blosm_aggregate_surface_owners'] = AggregateSurfaceCoverage.prepare(plan, tables).encode()
        controller['blosm_projection_coordinate_offset'] = origin
        ownership = bpy.data.texts.new(name+' / aggregate input ownership')
        # Preserve the record across save/reopen despite its string reference.
        ownership.use_fake_user = True
        ownership.write(json.dumps(asdict(plan), indent=2))
        controller['blosm_aggregate_ownership'] = ownership.name
        if self.solver.get('blosm_aggregate_provider_input'):
            self._provider(controller, private)
        mod.show_viewport = True
        return AggregateAuthoredNetwork(controller, objects, private, tuple(sources))

    def _provider(self, controller, collection):
        """Native dependencies for geometry edits; topology rebinding stays explicit."""
        from ..gn.editable_input.curved_patch_provider_bindings import input_values, set_inputs
        # Share original input data, never the controller's evaluated surface.
        # This also keeps the normal Group Input Geometry path live without
        # introducing an Object Info cycle back to the output controller.
        child = bpy.data.objects.new(controller.name+' / aggregate provider', controller.data)
        collection.objects.link(child)
        child.parent = controller
        child.hide_render = True; child.hide_select = True; child.hide_set(True)
        child['blosm_aggregate_provider'] = 1
        mod = child.modifiers.new('Shared native aggregate construction', 'NODES')
        mod.node_group = self.provider; mod.show_viewport = False
        used = {o.get('blosm_aggregate_provider_key') for o in bpy.data.objects}
        key = 1
        while key in used: key += 1
        if key >= 2**24: raise ValueError('Aggregate provider identity capacity exceeded')
        # Native numeric drivers transport a float; use exactly representable
        # integers, and allocate uniquely within this saved Blender file.
        controller['blosm_aggregate_provider_key'] = key
        set_input(controller, 'Aggregate Provider Key', key)
        set_input(controller, 'Aggregate Result', child)
        set_inputs(child, input_values(controller))
        source = {s.name: getattr(controller.modifiers[0].properties.inputs, s.identifier)
            for s in self.solver.interface.items_tree if s.item_type == 'SOCKET' and s.in_out == 'INPUT'
            and s.socket_type in ('NodeSocketInt', 'NodeSocketFloat')}
        for socket in self.provider.interface.items_tree:
            if socket.item_type != 'SOCKET' or socket.in_out != 'INPUT' or socket.name not in source: continue
            prop = getattr(mod.properties.inputs, socket.identifier)
            curve = prop.driver_add('value'); curve.driver.type = 'AVERAGE'; curve.extrapolation = 'LINEAR'
            for point in curve.keyframe_points: point.interpolation = 'LINEAR'
            var = curve.driver.variables.new(); var.name = 'control'; var.type = 'SINGLE_PROP'
            var.targets[0].id = controller; var.targets[0].data_path = source[socket.name].path_from_id('value')
        mod.show_viewport = True

    @staticmethod
    def _sources(mesh, groups, controller, collection):
        """Preserve global source identities through grouped/per-Street objects."""
        streets = mesh.attributes['st_street_id'].data
        declared = {d.value for d in streets}; assigned = [sid for _, ids in groups for sid in ids]
        if len(assigned) != len(set(assigned)) or set(assigned) != declared:
            raise ValueError('Source object groups must cover each complete Street exactly once')
        result = []
        for name, ids in groups:
            edges = [e for e in mesh.edges if streets[e.index].value in ids]
            points = sorted({v for e in edges for v in e.vertices}); local = {v: i for i, v in enumerate(points)}
            data = bpy.data.meshes.new(name+' source')
            data.from_pydata([tuple(mesh.vertices[i].co) for i in points],
                             [tuple(local[v] for v in e.vertices) for e in edges], [])
            for attr in mesh.attributes:
                if attr.name == 'position' or attr.name.startswith('.') or attr.domain not in ('POINT', 'EDGE'): continue
                if attr.data_type not in ('INT', 'FLOAT', 'BOOLEAN'): continue
                indices = points if attr.domain == 'POINT' else [e.index for e in edges]
                data.attributes.new(attr.name, attr.data_type, attr.domain).data.foreach_set('value',
                    [attr.data[i].value for i in indices])
            obj = bpy.data.objects.new(name, data); collection.objects.link(obj)
            obj.location = controller.location; obj.display_type = 'WIRE'; obj.show_in_front = True
            if hasattr(obj, 'st_source_controller'): obj.st_source_controller = controller
            obj['blosm_aggregate_source'] = True
            result.append(obj)
        return result
