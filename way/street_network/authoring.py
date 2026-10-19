"""Author clean input geometry and persistent topology; GN owns final surfaces."""
from uuid import uuid4
import math
import bpy
from ..gn.editable_input import model, runtime


class SourceMeshBuilder:
    def create(self, name, sections, collection, controller=None, interpolation=3):
        vertices, edges = [], []
        scale = metre_scale()
        origin = tuple(controller.matrix_world.translation) if controller is not None else (0., 0., 0.)
        ids, widths, offsets, roles, modes = [], [], [], [], []
        for section in sections:
            points = section['points']
            if len(points) < 2:
                continue
            start = len(vertices)
            vertices.extend(tuple(p[k]/scale-origin[k] for k in range(3)) for p in points)
            for i in range(len(points)-1):
                edges.append((start+i, start+i+1))
                ids.append(section['id']); widths.append(section['width'] or 0.)
                offsets.append(section['offset'] or 0.); modes.append(interpolation)
            roles.extend([1]+[0]*(len(points)-2)+[2])
        mesh = bpy.data.meshes.new(name+' source')
        mesh.from_pydata(vertices, edges, [])
        for name_attr, kind, domain, values in (
            ('st_section_id','INT','EDGE',ids), ('st_width','FLOAT','EDGE',widths),
            ('st_offset','FLOAT','EDGE',offsets), ('st_interpolation','INT','EDGE',modes),
            ('st_endpoint_role','INT','POINT',roles)):
            mesh.attributes.new(name_attr, kind, domain).data.foreach_set('value', values)
        obj = bpy.data.objects.new(name, mesh)
        if controller is not None:
            obj.matrix_world = controller.matrix_world.copy()
            obj['blosm_projection_coordinate_offset'] = tuple(v*scale for v in origin)
        collection.objects.link(obj)
        obj.st_source_controller = controller
        obj.display_type = 'WIRE'; obj.show_in_front = True
        obj['blosm_osm_way_ids'] = ','.join(sorted({w for s in sections for w in s['source_way_ids']}))
        return obj


def metre_scale():
    """Prepared coordinates/dimensions are metres, Blender coordinates are not."""
    scale = float(bpy.context.scene.unit_settings.scale_length)
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError('Scene metre conversion must be finite and positive')
    return scale


class EditableNetworkBuilder:
    def __init__(self, design, assets):
        self.design = design
        self.assets = assets
        self.meshes = SourceMeshBuilder()

    def create(self, name, snapshot, street_ids, junction_ids, collection):
        model.register()
        mesh = bpy.data.meshes.new(name+' input')
        obj = bpy.data.objects.new(name, mesh)
        collection.objects.link(obj)
        points = [p for sid in street_ids for s in snapshot.streets[sid]['items']
                  if s['kind']=='Section' for p in s['points']]
        # Sources and controller share a local frame. Equal transforms cancel
        # exactly in Object Info, including after reversible gallery placement.
        # World coordinates and the saved offset retain the import projection.
        if points:
            scale = metre_scale()
            obj.matrix_world.translation = tuple((min(p[k] for p in points)+max(p[k] for p in points))*.5/scale
                                                   if k<2 else 0. for k in range(3))
        obj['st_network_v8'] = True
        props = obj.st_network
        props.enabled = True; props.network_uuid = str(uuid4())
        props.resolution = self.design.resolution
        props.exaggeration = self.design.exaggeration
        props.cell_size = self.design.cell_size
        props.solver_group = self.assets['solver']
        props.pavement_material = self.assets['pavement']
        props.grass_material = self.assets['grass']
        sources = bpy.data.collections.new(name+' / editable Streets')
        collection.children.link(sources)
        props.source_collection = sources
        with runtime.suspend():
            for street_id in sorted(street_ids):
                street = snapshot.streets[street_id]
                record = props.streets.add(); record.street_id = street_id
                sections = [i for i in street['items'] if i['kind'] == 'Section']
                for section in sections:
                    item = record.items.add(); item.section_id = section['id']; item.forward = True
                    s = props.sections.add(); s.section_id = section['id']; s.street_id = street_id
                    s.rise = self.design.profile_rise
                    s.preset = 'SMOOTH_CROWN' if s.rise else 'FLAT'
                    s.clip_start = section.get('clip_start', 0.)
                    s.clip_end = section.get('clip_end', 0.)
                self.meshes.create(f'{name} / Street {street_id}', sections, sources, obj, self.design.interpolation)
                for prepared in street['items']:
                    if not prepared.get('transition'):
                        continue
                    transition = props.transitions.add()
                    transition.transition_id = len(props.transitions)-1
                    transition.street_id = street_id
                    for key,value in prepared['transition'].items():
                        setattr(transition,key,value)
            bundles = {}
            for jid in sorted(junction_ids):
                junction = snapshot.junctions[jid]
                row = props.intersections.add(); row.intersection_id = jid
                row.position = tuple(junction.location[:2])+(0.,)
                row.datum = 0.; row.radius = 2.
                for index, c in enumerate(junction.connectors):
                    record = props.connectors.add(); record.connector_id = c.source_id
                    record.intersection_id = jid
                    record.radius_override = c.radius
                    record.radius_next_connector_id = junction.connectors[(index+1)%len(junction.connectors)].source_id
                    record.use_profile_override = True
                    profile = junction.profiles[c.profile_index]
                    record.rise = profile.rise; record.separator = profile.separator
                    record.preset = 'SMOOTH_CROWN' if profile.rise else 'FLAT'
                    if c.item_kind == 1:
                        record.bundle_id = c.bundle_source_id
                        bundles[c.bundle_source_id] = (junction, c)
                    for member in c.members:
                        street = junction.streets[member.street_index]
                        section = street.sections[0 if member.leaving else -1]
                        r = record.members.add(); r.section_id = section.source_id
                        r.section_end = 0 if member.leaving else 1
                        r.lateral_order = member.lateral_order
            for bid, (junction, connector) in bundles.items():
                record = props.bundles.add(); record.bundle_id = bid
                for m in connector.members:
                    r = record.members.add(); r.street_id = junction.streets[m.street_index].source_id
                    r.forward = m.leaving == connector.leaving
                    r.head_order = m.head_order; r.tail_order = m.tail_order
            for bid, bundle in snapshot.bundles.items():
                if bid in bundles or not set(bundle['head']) <= set(street_ids):
                    continue
                if bundle.get('status') != 'SUPPORTED':
                    raise ValueError(bundle.get('reason') or 'Unsupported Bundle')
                record = props.bundles.add(); record.bundle_id = bid
                for member in bundle['members']:
                    r = record.members.add()
                    for key in ('street_id', 'forward', 'head_order', 'tail_order'):
                        setattr(r, key, member[key])
            for entity in ('section','street','intersection','connector','bundle'):
                entries = getattr(props, entity+'s' if entity != 'vertex' else 'vertices')
                setattr(props, 'next_'+entity+'_id', max((getattr(v, entity+'_id') for v in entries),default=-1)+1)
        runtime.refresh(obj, force=True, raise_errors=True)
        return obj
