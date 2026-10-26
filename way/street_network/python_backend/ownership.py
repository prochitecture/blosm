"""Shared physical junction ownership before per-feature surface dispatch."""
from dataclasses import dataclass, replace
import numpy as np
from .compounds import CompoundGroup, CompoundPlanner, rim_label, sample_footprint
from .compound_surfaces import CompoundSurfaceBuilder
from .conditions import BoundaryConditions
from .mesh import Assembler
from .model import Key
from .providers import Part
from .regions import RegionUnion, _classify
from .updates import Task

PLAN = Key('compound_plan', 0)
SURFACES_READY = Key('junction_surfaces_ready', 0)


class CurbChainUnavailable(ValueError):
    """A junction has no uninterrupted curb between two requested ports."""


def has_junctions(network):
    return any(len(j.ends) > 2 for j in network.intersections)


def footprint_key(owner):
    return Key('physical_' + owner.kind, owner.id)


def surface_key(owner):
    return Key('junction_surface', owner.id)


@dataclass
class OwnershipPlan:
    groups: tuple
    by_site: dict

    @classmethod
    def build(cls, context, junctions):
        layouts = {j.key: context[footprint_key(j.key)] for j in junctions}
        from .bundle_merged_entrances import preparation_key
        prepared = tuple(context.products.get(preparation_key(bundle)) for bundle in context.network.bundles)
        merged_groups = tuple(spec.group for spec in prepared if spec is not None)
        merged_sites = {site for group in merged_groups for site in group.members}
        ordinary = {key: layout for key, layout in layouts.items() if key not in merged_sites}
        groups = CompoundPlanner().build(ordinary)
        from .compound_ports import CompoundPortPlanner
        planner = CompoundPortPlanner(context, junctions, ordinary)
        groups = tuple(planner.build(group) for group in groups)+merged_groups
        planner.validate_deferred_profiles(groups)
        return cls(groups, {site: group for group in groups for site in group.members})

    def owner(self, site):
        group = self.by_site.get(site)
        return group.members[0] if group is not None else None


@dataclass
class CompoundProduct:
    owner: Key
    group: object
    field: object
    rims: dict
    kinds: dict
    peaks: dict
    port_points: dict
    pavement: Part
    grass: object
    exaggeration: float
    physical_rims: dict = None

    @classmethod
    def build(cls, context, owner, group):
        surface = CompoundSurfaceBuilder(context.patterns).build(group, context)
        group = surface.chart.partition.group
        layouts = {p.region.owner: p.layout for p in group.footprints}
        data = [BoundaryConditions().build(group, loop, layouts) for loop in group.loops]
        sources = {p.region.owner: p for p in group.footprints}
        rims, kinds, peaks, port_points, physical_rims = {}, {}, {}, {}, {}
        for li, ids in surface.mesh.seams.items():
            rims[li] = surface.mesh.xyz[ids].copy()
            physical_rims[li] = rims[li].copy()
            for row, vertex in enumerate(ids):
                ref = surface.boundary_refs[int(vertex)]
                physical_rims[li][row, 2] = data[ref.loop].rules[ref.edge].sample(ref.fraction)[0]
            labels = []
            for vertex in ids:
                ref = surface.boundary_refs[int(vertex)]
                edge = (ref.edge + int(ref.fraction >= 1. - 1e-9)) % len(group.loops[li].spans)
                labels.append(rim_label(group.loops[li].spans[edge], sources).kind)
            kinds[li] = tuple(labels)
            for port in data[li].ports:
                profile = layouts[port.owner].approaches[port.approach_index].profile
                u = .5 * (profile.peak + 1.)
                at = int(np.argmin(abs(np.asarray(port.outward_u) - u)))
                xy = data[li].xy[port.vertices[at]]
                index = int(np.argmin(np.linalg.norm(rims[li][:, :2] - xy, axis=1)))
                if np.linalg.norm(rims[li][index, :2] - xy) > 1e-7:
                    raise ValueError('Compound output lost an authoritative incoming profile peak')
                from .source_mouths import source_mouths
                aliases = source_mouths(layouts[port.owner].approaches[port.approach_index],
                                        layouts[port.owner].outward_stations[port.approach_index])
                for alias in aliases:
                    peaks[alias.endpoint.section_id, alias.endpoint.end] = (li, index, alias.station)
                for vertex, u in zip(port.vertices, port.outward_u):
                    xy = data[li].xy[vertex]
                    at = int(np.argmin(np.linalg.norm(rims[li][:, :2] - xy, axis=1)))
                    if np.linalg.norm(rims[li][at, :2] - xy) > 1e-7:
                        raise ValueError('Boundary output lost an incoming profile vertex')
                    for alias in aliases:
                        if alias.contains(u):
                            original_u = float(alias.outward_parameter(u))
                            port_points[alias.endpoint.section_id, alias.endpoint.end, round(original_u, 12)] = (li, at, alias.station)
        identity = Key('compound', owner.id) if len(group.members) > 1 else owner
        metadata = dict(surface.metadata, member_junctions=[key.text() for key in group.members],
                        consumed_sections=sorted({c.section_id for c in group.claims}),
                        physical_exterior_ports=group.physical_layout is not None,
                        default_provider_integration=True)
        from .corner_radii import fitted_radius_report
        metadata['fitted_corner_radii'] = [fit for key, layout in layouts.items()
                                           for fit in fitted_radius_report(key, layout)]
        pavement = Part(identity, surface.mesh, metadata)
        grass = None
        if surface.holes:
            assembly = Assembler()
            seams = {}
            for li, mesh in surface.holes.items():
                offset = len(assembly.xyz)
                assembly.add(mesh.xyz, mesh.quads)
                seams[li] = mesh.seams[0] + offset
            mesh = assembly.mesh('grass')
            mesh.seams = seams
            grass = Part(identity, mesh, dict(quality=mesh.validate(), holes=sorted(surface.holes),
                construction='compound separator holes', member_junctions=metadata['member_junctions']))
        return cls(owner, group, surface.field, rims, kinds, peaks, port_points, pavement, grass,
                   context.network.settings.exaggeration, physical_rims)

    def contains(self, xy):
        """Actual physical ownership, including its separator holes."""
        xy = np.asarray(xy)
        inside = np.zeros(len(xy), bool)
        for loop in self.group.loops:
            if not loop.hole:
                hit, boundary = _classify(xy, loop.xy, 1e-7)
                inside |= hit | (boundary >= 0)
        return inside

    def height(self, xy):
        xy = np.asarray(xy)
        inside = np.zeros(len(xy), bool)
        for loop in self.group.loops:
            if not loop.hole:
                hit, boundary = _classify(xy, loop.xy, 1e-7)
                inside |= hit | (boundary >= 0)
        holes = np.zeros(len(xy), bool)
        for loop in self.group.loops:
            if loop.hole:
                hit, boundary = _classify(xy, loop.xy, 1e-7)
                holes |= hit | (boundary >= 0)
        result = np.full(len(xy), self.field.datum)
        selected = inside & ~holes
        result[selected] = self.field.datum + self.exaggeration * (self.field(xy[selected]) - self.field.datum)
        result[holes] = self.field.datum + self.exaggeration * (self.field.ceiling - self.field.datum)
        return result

    def cap(self, first, last):
        return self._cap_points(self.peaks[first], self.peaks[last])

    def edge_cap(self, first, last, first_u, last_u):
        return self._cap_points(self.port_points[(*first, round(first_u, 12))],
                                self.port_points[(*last, round(last_u, 12))])

    def _cap_points(self, a, b):
        if a[0] != b[0]:
            raise CurbChainUnavailable('Bundle terminals are on different compound boundary loops')
        li, start, _ = a
        finish = b[1]
        size = len(self.rims[li])
        candidates = []
        for step in (1, -1):
            indices, at, valid = [start], start, True
            while at != finish:
                edge = at if step == 1 else (at - 1) % size
                if self.kinds[li][edge] != 'curb':
                    valid = False
                    break
                at = (at + step) % size
                indices.append(at)
            if valid:
                candidates.append(self.rims[li][indices])
        if len(candidates) != 1:
            raise CurbChainUnavailable('Bundle compound cap must follow one curb chain without crossing another mouth')
        return candidates[0]

    def release_geometry(self):
        # Reopening can reconstruct this query product while hydrated output
        # Parts retain their identities. Release both aliased and unused meshes.
        for part in (self.pavement, self.grass):
            if part is not None and part.mesh is not None:
                part.metadata['bounds'] = part.mesh.bounds
                part.mesh = None


@dataclass
class StandaloneProduct:
    pavement: Part
    grass: object = None
    field: object = None
    layout: object = None
    exaggeration: float = 1.

    def height(self, xy):
        if self.field is None:
            raise ValueError('Standalone pavement has no attachment height query')
        return self.field.evaluate(xy, default=self.field.datum)

    def cap(self, first, last):
        """Actual band boundary joining two inner peaks, with no intervening port."""
        if self.layout is None:
            raise CurbChainUnavailable('Standalone pavement has no curb layout')
        ports = {list(a.path.trims(s))[-1]: a for a, s in
                 zip(self.layout.approaches, self.layout.outward_stations)}
        return self.edge_cap(first, last, .5 * (ports[first].profile.peak + 1.),
                             .5 * (ports[last].profile.peak + 1.))

    def edge_cap(self, first, last, first_u, last_u):
        """One actual outer curb; full crowns keep their zero-height edges."""
        from .intersections import IntersectionProvider
        if self.layout is None:
            raise CurbChainUnavailable('Standalone pavement has no curb layout')
        ports = {list(a.path.trims(s))[-1]: a for a, s in
                 zip(self.layout.approaches, self.layout.outward_stations)}
        a, b = ports[first], ports[last]
        for side, corner in enumerate(self.layout.corners):
            forward = corner.first is a and corner.second is b and first_u == 0. and last_u == 1.
            reverse = corner.first is b and corner.second is a and last_u == 0. and first_u == 1.
            if forward or reverse:
                t = np.linspace(0., 1., self.layout.rows + 1)
                xy, z = IntersectionProvider.band_surface(self.layout, side, t, np.zeros_like(t))
                xyz = np.column_stack((xy, a.profile.datum + self.exaggeration * z))
                return xyz if forward else xyz[::-1]
        raise CurbChainUnavailable('No adjacent inner-profile curb joins these standalone ports')

    def release_geometry(self):
        if self.pavement is not None and self.pavement.mesh is not None:
            self.pavement.metadata['bounds'] = self.pavement.mesh.bounds
            self.pavement.mesh = None


class JunctionOwnership:
    @staticmethod
    def tasks(network):
        junctions = tuple(j for j in network.intersections if len(j.ends) > 2)
        if not junctions:
            return
        yield Task(PLAN, tuple(footprint_key(j.key) for j in junctions),
                   lambda context: OwnershipPlan.build(context, junctions))
        for junction in junctions:
            yield Task(surface_key(junction.key), (PLAN, junction.layout_key),
                lambda context, site=junction.key: JunctionOwnership.surface(context, site),
                signature=lambda context, site=junction.key: JunctionOwnership.surface_signature(context, site))
            yield Task(Key('junction_grass', junction.id), (surface_key(junction.key),),
                lambda context, site=junction.key: context[surface_key(site)].grass
                    if context[surface_key(site)] is not None else None)
        yield Task(SURFACES_READY, tuple(surface_key(j.key) for j in junctions), lambda context: None)

    @staticmethod
    def surface_signature(context, site):
        plan = context[PLAN]
        owner = plan.owner(site)
        if owner is None or owner != site:
            return owner.text() if owner is not None else context.digest(Key('mouths', site.id))
        return [(key.text(), context.digest(footprint_key(key))) for key in plan.by_site[site].members]

    @staticmethod
    def surface(context, site):
        plan = context[PLAN]
        owner = plan.owner(site)
        if owner is not None:
            return CompoundProduct.build(context, site, plan.by_site[site]) if owner == site else None
        from .intersections import IntersectionProvider
        junction = next(j for j in context.network.intersections if j.key == site)
        layout = context[junction.layout_key]
        if layout.mesh_strategy == 'bands':
            try:
                part = IntersectionProvider.build(context, junction)
                from .conditions import _profile_gradient
                incoming = max(float(np.linalg.norm(_profile_gradient(a, s, np.array((-1., 0., 1.))), axis=1).max())
                               for a, s in zip(layout.approaches, layout.outward_stations))
                if part.metadata['physical_slopes']['maximum'] <= max(.10, 1.1 * incoming):
                    from .surface_query import QuadHeightField
                    datum = layout.approaches[0].profile.datum
                    ceiling = datum + context.network.settings.exaggeration * max(a.profile.rise for a in layout.approaches)
                    return StandaloneProduct(part, field=QuadHeightField(part.mesh, datum=datum, ceiling=ceiling),
                                             layout=layout, exaggeration=context.network.settings.exaggeration)
            except ValueError:
                pass  # A valid physical rim can need another internal chart.
        footprint = sample_footprint(site, context[footprint_key(site)])
        loops = RegionUnion().build((footprint.region,))
        group = CompoundGroup((site,), (), (footprint,), loops)
        return CompoundProduct.build(context, site, group)

    @staticmethod
    def layout_signature(context, junction):
        plan = context[PLAN]
        owner = plan.owner(junction.key)
        group = plan.by_site.get(junction.key)
        if group is not None and group.physical_layout is not None:
            return ('exterior_ports', [(key.text(), context.digest(footprint_key(key))) for key in group.members])
        return context.digest(footprint_key(junction.key)), owner.text() if owner is not None else None

    @staticmethod
    def layout(context, junction):
        from .intersections import IntersectionProvider
        footprint = context[footprint_key(junction.key)]
        group = context[PLAN].by_site.get(junction.key)
        if group is not None and group.physical_layout is not None:
            return group.physical_layout
        if context[PLAN].owner(junction.key) is not None:
            return footprint
        try:
            return IntersectionProvider.layout(context, junction, footprint)
        except ValueError:
            return replace(footprint, mesh_strategy='boundary_chart')

    @staticmethod
    def pavement(context, junction):
        from .intersections import IntersectionProvider
        if len(junction.ends) > 2:
            product = context[surface_key(junction.key)]
            return product.pavement if product is not None else None
        return IntersectionProvider.build(context, junction)
