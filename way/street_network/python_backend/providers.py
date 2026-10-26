"""Composable numerical providers for paths, carriageways and source anchors."""
from dataclasses import dataclass
import math
import numpy as np
from .mesh import Mesh
from .model import Key
from .paths import Path, unit
from .updates import FeatureProvider, Task


@dataclass
class Part:
    owner: Key
    mesh: Mesh
    metadata: dict


@dataclass(frozen=True)
class RoadSpan:
    """Visible source span whose surface is owned by a composite feature."""
    metadata: dict


@dataclass(frozen=True)
class RoadSampling:
    """Retained row policy shared by pavement and downstream boundaries."""
    first: float
    last: float
    local_curvature: bool = False

    @classmethod
    def from_metadata(cls, metadata):
        return cls(metadata['start'], metadata['end'], metadata.get('local_curvature_sampling', False))

    def samples(self, path, settings):
        return path.samples(self.first, self.last, settings.spacing, local_curvature=self.local_curvature)


@dataclass(frozen=True)
class AnchorTransform:
    position: tuple
    heading: float
    asset: str


class PathProvider(FeatureProvider):
    source_kinds = frozenset(('section',))
    def tasks(self, network):
        for section in network.sections:
            neighbors = [j for j in (*network.intersections, *network.connections) if len(j.ends) == 2
                         and any(e.section_id == section.id for e in j.ends)]
            deps = tuple(sorted({section.key, *(Key('section', e.section_id) for j in neighbors for e in j.ends),
                                 *(j.key for j in neighbors)}))
            yield Task(Key('path', section.id), deps,
                lambda context, source=section, seams=neighbors: self.build(context, source, seams))

    @staticmethod
    def build(context, section, seams):
        endpoints = {}
        for seam in seams:
            own = next(e for e in seam.ends if e.section_id == section.id)
            partner = next(e for e in seam.ends if e != own)
            other = context.network.section_index[partner.section_id]
            def outward(source, end):
                points = np.asarray(source.points)[:, :2]
                return unit(points[1] - points[0]) if end == 0 else unit(points[-2] - points[-1])
            common = unit(outward(section, own.end) - outward(other, partner.end))
            endpoints[own.end] = common if own.end == 0 else -common
        return Path(section, context.network.settings.curve_tolerance, endpoints)


class CarriagewayProvider(FeatureProvider):
    source_kinds = frozenset(('section',))
    def tasks(self, network):
        from .profiles import memberships
        from .routes import RouteTopology
        from .ownership import PLAN, has_junctions
        membership = memberships(network)
        attachments = RouteTopology(network).attachments(network)
        for section in network.sections:
            junctions = attachments[section.id]
            deps = (section.key, Key('path', section.id), Key('profile', section.id),
                    *sorted({junction.layout_key for junction, end in junctions}),
                    *((Key('bundle_plan', membership[section.id][0].id),) if section.id in membership else ()),
                    *((PLAN,) if has_junctions(network) else ()))
            yield Task(Key('road', section.id), deps,
                       lambda context, source=section, joins=junctions: self.build(context, source, joins),
                       version='2',
                       signature=lambda context, dependencies=deps, joins=junctions, sid=section.id:
                           self.signature(context, dependencies, joins, sid))

    @staticmethod
    def signature(context, dependencies, attachments, section_id):
        from .ownership import PLAN
        plan = context.products.get(PLAN)
        layout_keys = {junction.layout_key for junction, end in attachments}
        consumed = []
        for key in dependencies:
            if key == PLAN or key in layout_keys:
                continue
            # A composite owns this road only according to its mode. Its full
            # geometry belongs to the separate bundle_pavement task.
            value = context[key].mode if key.kind == 'bundle_plan' else context.digest(key)
            consumed.append((key.text(), value))
        mouths = []
        for junction, end in attachments:
            layout = context[junction.layout_key]
            mouth = layout.mouths.get((section_id, end))
            if mouth is None:
                continue  # A potential route that does not reach this Section.
            owner = plan.owner(junction.key) if plan is not None else None
            mouths.append((junction.key.text(), end, mouth.station, mouth.transverse_intervals,
                           owner.text() if owner is not None else None,
                           layout.transition is not None, len(junction.ends) > 2))
        return consumed, mouths

    @staticmethod
    def span_metadata(context, section, attachments):
        path = context[Key('path', section.id)]
        first, last = 0., path.length
        from .ownership import PLAN
        plan = context.products.get(PLAN)
        covered, junction_spans, transition_spans = [], [], []
        for junction, end in attachments:
            mouth = context[junction.layout_key].mouths.get((section.id, end))
            if mouth is None:
                continue
            owner = plan.owner(junction.key) if plan is not None else None
            if owner is not None:
                a, b = (0., mouth.station) if end == 0 else (mouth.station, path.length)
                if b > a + 1e-8:
                    covered.append(dict(start=a, end=b, owner=owner.text()))
            elif context[junction.layout_key].transition is not None:
                from .transitions import product_key
                a, b = (0., mouth.station) if end == 0 else (mouth.station, path.length)
                transition_spans.append(dict(start=a, end=b, source_end=end, owner=product_key(junction.key).text()))
            elif len(junction.ends) > 2:
                a, b = (0., mouth.station) if end == 0 else (mouth.station, path.length)
                if b > a + 1e-8:
                    junction_spans.append(dict(start=a, end=b, owner=junction.key.text()))
            if end == 0:
                first = max(first, mouth.station)
            else:
                last = min(last, mouth.station)
        metadata = dict(start=first, end=last, compound_spans=covered,
                        junction_spans=junction_spans, transition_spans=transition_spans)
        refine = any(len(j.ends) == 2 and context[j.layout_key].transition is not None
            and 'curvature_refinement' in context[j.layout_key].transition.certification for j, _ in attachments)
        if refine:
            metadata['local_curvature_sampling'] = True
        return metadata

    @staticmethod
    def build(context, section, attachments):
        path = context[Key('path', section.id)]
        profile = context[Key('profile', section.id)]
        metadata = CarriagewayProvider.span_metadata(context, section, attachments)
        first, last = metadata['start'], metadata['end']
        covered, junction_spans, transition_spans = (metadata[name] for name in
            ('compound_spans', 'junction_spans', 'transition_spans'))
        settings = context.network.settings
        from .profiles import memberships
        member = memberships(context.network).get(section.id)
        if member is not None and context[Key('bundle_plan', member[0].id)].mode == 'union':
            return RoadSpan(dict(start=first, end=last, composite_owner=member[0].key.text(),
                compound_spans=covered, junction_spans=junction_spans, transition_spans=transition_spans))
        if last < first - 1e-7:
            owners = {span['owner'] for span in covered}
            if len(owners) == 1 and min(span['start'] for span in covered) <= 1e-8 and max(span['end'] for span in covered) >= path.length - 1e-8:
                return RoadSpan(dict(start=first, end=last, compound_spans=covered,
                                     consumed_by=sorted(owners), junction_spans=junction_spans, transition_spans=transition_spans))
            raise ValueError(f'Section {section.id} intersection trims overlap; compound junction ownership is required')
        if last <= first + 1e-7:
            owners = sorted({j.key.text() for j, end in attachments if len(j.ends) > 2
                             and (section.id, end) in context[j.layout_key].mouths})
            return RoadSpan(dict(start=first, end=last, consumed_by=owners, compound_spans=covered,
                                 junction_spans=junction_spans, transition_spans=transition_spans))
        # Tighten sampling only beside the newly certified Smooth-seam
        # fallback. Preserve every previously qualified road's row layout.
        refine_sampling = metadata.get('local_curvature_sampling', False)
        stations = RoadSampling(first, last, refine_sampling).samples(path, settings)
        from .intersections import half_intervals
        across = 2 * half_intervals(settings)
        # Path spacing is not assumed uniform; compulsory Exact control rows
        # and curvature samples remain where needed. The grid caches topology.
        grid = context.patterns.grid(len(stations) - 1, across)
        p, _, normal = path.frame(stations)
        q = np.linspace(-1., 1., across + 1)
        xy = p[:, None, :] + (section.offset + .5 * section.width * q)[None, :, None] * normal[:, None, :]
        xyz = np.empty((len(stations), across + 1, 3))
        xyz[..., :2] = xy
        xyz[..., 2] = profile.datum + profile.value(q) * settings.exaggeration
        # Increasing station followed by increasing rightward distance winds
        # clockwise in XY. Reverse every directly generated quad consistently.
        mesh = Mesh(xyz.reshape(-1, 3), grid.quads[:, ::-1])
        mesh.seams = {0: np.arange(across + 1),
                      1: np.arange((len(stations) - 1) * (across + 1), len(stations) * (across + 1))}
        quality = mesh.validate()
        metadata = dict(start=first, end=last, quality=quality,
                    transverse_intervals=across, profile_half=profile.half, compound_spans=covered,
                    junction_spans=junction_spans, transition_spans=transition_spans)
        if refine_sampling:
            metadata['local_curvature_sampling'] = True
        return Part(section.key, mesh, metadata)


class FurnitureProvider(FeatureProvider):
    source_kinds = frozenset(('furniture',))
    def tasks(self, network):
        from .ownership import SURFACES_READY, has_junctions
        from .transitions import required, product_key
        from .profiles import memberships
        membership = memberships(network)
        for item in network.furniture:
            transitions = tuple(product_key(j.key) for j in (*network.intersections, *network.connections)
                                if required(network, j) and any(e.section_id == item.section_id for e in j.ends))
            deps = (item.key, Key('path', item.section_id), Key('road', item.section_id), Key('profile', item.section_id),
                    *transitions,
                    *((Key('bundle_plan', membership[item.section_id][0].id),
                       Key('separator', membership[item.section_id][0].id)) if item.section_id in membership else ()),
                    *((SURFACES_READY,) if has_junctions(network) else ()))
            yield Task(Key('placement', item.id), deps,
                lambda context, source=item: self.build(context, source),
                version='2', signature=lambda context, source=item, dependencies=deps:
                    self.signature(context, source, dependencies))

    @staticmethod
    def compound_owner(context, placement):
        road = context[Key('road', placement.section_id)]
        return next((Key.parse(span['owner']) for span in (*road.metadata.get('compound_spans', ()), *road.metadata.get('junction_spans', ()))
                     if span['start'] - 1e-8 <= placement.station <= span['end'] + 1e-8), None)

    @classmethod
    def signature(cls, context, placement, dependencies):
        from .ownership import SURFACES_READY, surface_key
        owner = cls.compound_owner(context, placement)
        composite = context[Key('road', placement.section_id)].metadata.get('composite_owner')
        separator_height = False
        if composite:
            owner_key = Key.parse(composite)
            separator_height = getattr(context[Key('bundle_plan', owner_key.id)], 'requires_separator_query', False)
        return ([(key.text(), context.digest(key)) for key in dependencies if key != SURFACES_READY
                 and (key.kind != 'bundle_plan' or composite)
                 and (key.kind != 'separator' or separator_height)],
                (owner.text(), context.digest(surface_key(owner))) if owner is not None else None)

    @staticmethod
    def build(context, placement):
        section = context.network.section_index[placement.section_id]
        path = context[Key('path', placement.section_id)]
        road = context[Key('road', placement.section_id)]
        profile = context[Key('profile', placement.section_id)]
        transition = next((Key.parse(span['owner']) for span in road.metadata.get('transition_spans', ())
                           if span['start'] - 1e-8 <= placement.station <= span['end'] + 1e-8), None)
        if transition is not None:
            result = context[transition].anchor(placement, context.network.settings.exaggeration)
            composite = road.metadata.get('composite_owner')
            if composite:
                from dataclasses import replace
                bundle = Key.parse(composite)
                height = float(context[Key('bundle_plan', bundle.id)].height(np.asarray(result.position[:2])))
                result = replace(result, position=(*result.position[:2], height + placement.height))
            return result
        owner = FurnitureProvider.compound_owner(context, placement)
        if (placement.station > path.length + 1e-8 or
                (not road.metadata['start'] - 1e-8 <= placement.station <= road.metadata['end'] + 1e-8 and owner is None)):
            raise ValueError(f'Furniture {placement.id} station lies outside its generated road span')
        side = 1. if placement.side == 'RIGHT' else -1.
        p, tangent, normal = path.frame(placement.station)
        position = p + (section.offset + side * (.5 * section.width + placement.offset)) * normal
        q = np.clip(side * (1. + 2. * placement.offset / section.width), -1., 1.)
        height = profile.datum + float(profile.value(q)) * context.network.settings.exaggeration
        if road.metadata.get('composite_owner'):
            bundle_owner = Key.parse(road.metadata['composite_owner'])
            plan = context[Key('bundle_plan', bundle_owner.id)]
            height = float(plan.height(position))
        if owner is not None:
            from .ownership import surface_key
            junction = context[surface_key(owner)]
            # A merged radial mouth need not be perpendicular to either
            # original centreline. Its source station interval alone cannot
            # classify an offset furniture point's physical ownership.
            if (not road.metadata.get('composite_owner') or not hasattr(junction, 'contains')
                    or junction.contains(position[None])[0]):
                height = float(junction.height(position[None])[0])
        return AnchorTransform(tuple(map(float, position)) + (height + placement.height,),
            float(math.atan2(tangent[1], tangent[0])), placement.asset)
