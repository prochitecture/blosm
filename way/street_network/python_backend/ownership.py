"""Shared physical junction ownership before per-feature surface dispatch."""
from dataclasses import dataclass, field, replace
import numpy as np
from .compounds import CompoundGroup, CompoundPlanner, rim_label, sample_footprint
from .compound_surfaces import CompoundSurfaceBuilder
from .conditions import BoundaryConditions
from .mesh import Assembler
from .model import Key
from .providers import Part
from .regions import RegionUnion, _classify
from .updates import CooperativeBuild, Task, drain_steps
from .coupled_source_owners import COUPLED_OWNER_POLICY
from .ordinary_mouth_clearance import ORDINARY_MOUTH_POLICY
from .source_strip_containment import SOURCE_RANGE_POLICY
from .source_curb_coordinates import SOURCE_CURB_COORDINATE_POLICY
from .source_bank_reconciliation import SOURCE_BANK_POLICY
from .retained_arc_fillets import SOURCE_ARC_POLICY
from .source_owned_connector_fillets import SOURCE_OWNED_CONNECTOR_POLICY
from .source_parent_chain_fillets import SOURCE_PARENT_CHAIN_POLICY
from .residual_road_gaps import SOURCE_RESIDUAL_GAP_POLICY
from .owned_bundle_normal_atlas import SOURCE_OWNED_BUNDLE_NORMAL_POLICY
from .owned_prefix_corner_authority import SOURCE_OWNED_PREFIX_CORNER_POLICY
from .source_root_templates import SOURCE_ROOT_POLICY

PLAN = Key('compound_plan', 0)
SURFACES_READY = Key('junction_surfaces_ready', 0)
# Projected ownership signatures intentionally omit the global plan digest.
# Keep its source-sampling epoch when retaining otherwise unchanged inputs.
SOURCE_SAMPLING_POLICY = '4-complete-future-grids-with-local-curvature-bounds'
SOURCE_HEIGHT_POLICY = '8-retained-crossfall-and-disjoint-paired-cores'


class CurbChainUnavailable(ValueError):
    """A junction has no uninterrupted curb between two requested ports."""


def has_junctions(network):
    return any(len(j.ends) > 2 for j in network.intersections)


def footprint_key(owner):
    return Key('physical_' + owner.kind, owner.id)


def surface_key(owner):
    return Key('junction_surface', owner.id)


def bundle_signature(context,group):
    """Local Bundle semantics used by a physical ownership construction."""
    sections={claim.section_id for claim in group.claims}
    for footprint in group.footprints:
        for approach in footprint.layout.approaches:
            sections.add(approach.source.id)
            sections.update(piece.section.id for piece in getattr(approach.path,'pieces',()))
    layout=group.physical_layout
    ids={claim.bundle for claim in getattr(layout,'partial_bundle_claims',())}
    ids.update(getattr(layout,'source_cycle_proofs',{}))
    from .bundle_merged_entrances import preparation_key
    result=[]
    if getattr(layout,'shared_bank_signature',()):
        result.append(('shared native bank sampling',layout.shared_bank_signature))
    if getattr(layout, 'shallow_bank_fillet', None):
        result.append(('shallow native bank fillet', layout.shallow_bank_fillet))
    for bundle in sorted(context.network.bundles,key=lambda b:b.id):
        if bundle.id not in ids and not sections.intersection(bundle.sections):
            continue
        signature=(bundle.key.text(),context.digest(bundle.key))
        prepared=context.products.get(preparation_key(bundle))
        if getattr(prepared,'source_phase_proof',None) is not None:
            # Joint sampling can change a local physical mouth when the
            # opposite member changes, although its raw local footprint and
            # Bundle membership are unchanged. Track that prepared resource.
            signature=(*signature,(preparation_key(bundle).text(),context.digest(preparation_key(bundle))))
        result.append(signature)
    derived = getattr(layout, 'paired_source_dependencies', ())
    if derived:
        result.append(('paired original sources', tuple((Key(kind, sid).text(), context.digest(Key(kind, sid)))
            for sid in derived for kind in ('path', 'profile'))))
    coupled = getattr(layout, 'coupled_source_dependencies', ())
    if coupled:
        result.append(('coupled original sources', tuple((Key(kind, sid).text(), context.digest(Key(kind, sid)))
            for sid in coupled for kind in ('path', 'profile'))))
    free_sources = getattr(layout, 'wholly_owned_free_source_dependencies', ())
    if free_sources:
        result.append(('whole owned free counterparts', tuple((Key(kind, sid).text(), context.digest(Key(kind, sid)))
            for sid in free_sources for kind in ('path', 'profile'))))
    from .ordinary_mouth_clearance import signature as ordinary_mouth_signature
    ordinary = ordinary_mouth_signature(group)
    if ordinary is not None:
        result.append(('complete ordinary source mouths', ordinary))
    from .paired_interval_atlas import signature as paired_prefix_signature
    prefix=paired_prefix_signature(context,group)
    if prefix is not None:
        result.append(('finite original paired-prefix charts',prefix))
    from .ordinary_source_prefix_atlas import signature as ordinary_prefix_signature
    ordinary_prefix=ordinary_prefix_signature(context,group)
    if ordinary_prefix is not None:
        result.append(('finite original ordinary-prefix charts',ordinary_prefix))
    from .source_bank_reconciliation import signature as bank_signature
    banks=bank_signature(group)
    if banks is not None:
        result.append(('original source-bank geometry',banks))
    arcs=getattr(layout,'source_retained_arc_policy',None)
    if arcs is not None:
        result.append(('retained original source-arc geometry',(SOURCE_ARC_POLICY,arcs)))
    joins=getattr(layout,'original_lead_joins',())
    if joins:
        result.append(('original auxiliary lead joins',tuple(record.signature() for record in joins)))
    bank_joins=getattr(layout,'original_bank_joins',())
    if bank_joins:
        result.append(('rounded exposed Section banks',tuple(record.signature() for record in bank_joins)))
    from .source_owned_connector_fillets import signature as connector_signature
    connectors=connector_signature(group)
    if connectors is not None:
        result.append(('finite owned original-source connector geometry',connectors))
    from .source_parent_chain_fillets import signature as parent_chain_signature
    parent_chain=parent_chain_signature(group)
    if parent_chain is not None:
        result.append(('finite certified original parent-chain geometry',parent_chain))
    from .residual_road_gaps import signature as residual_gap_signature
    residual_gaps=residual_gap_signature(group)
    if residual_gaps is not None:
        result.append(('certified residual road-gap fills',residual_gaps))
    from .source_root_templates import construction_signature as root_signature
    roots=root_signature(context,group)
    if roots is not None:
        result.append(('typed native-curvature original root templates',roots))
    if group.claims:
        from .owned_bundle_normal_atlas import signature as owned_normal_signature
        owned_normals=owned_normal_signature(context,group)
        if owned_normals is not None:
            result.append(('certified owned Bundle normal charts',owned_normals))
    # The original-body frame may use another visit on either declared Bundle
    # role, even when only a finite local prefix is owned by this junction.
    # Track those actual source inputs without acquiring their physical ranges.
    if getattr(layout,'source_cycle_full_width_proof',None) is not None:
        body_ids=set(sections)
        body_ids.update(row['section'] for row in layout.source_cycle_full_width_proof['sections'])
        local_body_ids=frozenset(body_ids)
        for bundle in context.network.bundles:
            if local_body_ids.intersection(bundle.sections):body_ids.update(bundle.sections)
        result.append(('original source-height inputs',tuple((Key(kind,sid).text(),context.digest(Key(kind,sid)))
            for sid in sorted(body_ids) for kind in ('path','profile'))))
    from .native_owned_cap_chains import construction_signature as cap_chain_signature
    from .ordinary_half_charts import applicable as half_chart_applies, NativeHalfJunction
    if half_chart_applies(group):
        result.append(('ordinary native half-profile chart', NativeHalfJunction.read(context, group).signature(context)))
    from .native_zero_bank_root import signature as endpoint_root_signature
    from .bundle_crest_roots import signature as crest_signature
    for label, construction in (('owned native cap chain', cap_chain_signature),
                                ('native Bundle crest', crest_signature),
                                ('native endpoint root', endpoint_root_signature)):
        signature = construction(context, group)
        if signature is not None:
            result.append((label, signature))
    return tuple(result)


def resolved_footprint(context, site):
    """Physical ownership replaces provisional inputs before geometry consumers."""
    plan = context.products.get(PLAN)
    group = plan.by_site.get(site) if plan is not None else None
    if group is not None and group.physical_layout is not None:
        return group.physical_layout
    footprint = context[footprint_key(site)]
    from .compound_preparation import DeferredCompoundFootprint
    if isinstance(footprint, DeferredCompoundFootprint):
        raise ValueError('Deferred junction footprint was consumed before physical ownership resolved it')
    return footprint


@dataclass
class OwnershipPlan:
    groups: tuple
    by_site: dict
    bundle_sources: dict = field(default_factory=dict)

    @classmethod
    def build(cls, context, junctions):
        return drain_steps(cls.build_steps(context, junctions))

    @classmethod
    def build_steps(cls, context, junctions):
        layouts = {j.key: context[footprint_key(j.key)] for j in junctions}
        from .bundle_merged_entrances import preparation_key
        prepared = tuple(context.products.get(preparation_key(bundle)) for bundle in context.network.bundles)
        merged_groups = tuple(group for spec in prepared if spec is not None
            for group in (spec.group,*(side.group for side in getattr(spec,'paired_sources',()))))
        from .finite_scope_ports import groups as finite_scope_groups
        merged_groups += finite_scope_groups(context, layouts)
        from .finite_free_caps import groups as finite_free_cap_groups
        merged_groups += finite_free_cap_groups(context, layouts)
        from .finite_scope_compounds import compose as compose_finite_scopes
        merged_groups = compose_finite_scopes(context,junctions,layouts,merged_groups)
        merged_sites = {site for group in merged_groups for site in group.members}
        ordinary = {key: layout for key, layout in layouts.items() if key not in merged_sites}
        from .compound_preparation import DeferredCompoundFootprint, DeferredCompoundPlanner
        qualified = {key:layout for key,layout in ordinary.items()
                     if not isinstance(layout,DeferredCompoundFootprint)}
        groups = CompoundPlanner().build(qualified)
        yield
        from .compound_port_fragments import normalize as normalize_fragments
        groups = tuple(normalize_fragments(context,group) for group in groups)
        groups = yield from DeferredCompoundPlanner(context,junctions,ordinary).build_steps(groups)
        yield
        from .compound_ports import CompoundPortPlanner
        planner = CompoundPortPlanner(context, junctions, ordinary)
        groups = tuple(planner.build(group) if group.physical_layout is None else group for group in groups)+merged_groups
        from .source_companions import prepare as prepare_source_companions
        groups, bundle_sources = prepare_source_companions(context, groups,
            prepared_bundle_ids=(spec.bundle_id for spec in prepared if spec is not None))
        from .coupled_source_owners import prepare_steps as prepare_coupled_source_owners
        groups = yield from prepare_coupled_source_owners(context, junctions, layouts, groups,
            excluded_bundle_ids=(*bundle_sources,
                *(spec.bundle_id for spec in prepared if spec is not None)))
        from .ordinary_mouth_clearance import finalize_steps
        yield
        groups = yield from finalize_steps(context, junctions, layouts, groups)
        from .source_bank_reconciliation import validate_prepared_neighbors
        groups = validate_prepared_neighbors(context,groups,layouts)
        planner.validate_deferred_profiles(groups)
        return cls(groups, {site: group for group in groups for site in group.members}, bundle_sources)

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
    derived_peaks: dict = None
    derived_port_points: dict = None
    physical_pavement: object = None
    physical_grass: object = None
    hole_fields: object = None

    @classmethod
    def build(cls, context, owner, group, *, source_heights=False):
        return drain_steps(cls.build_steps(context, owner, group, source_heights=source_heights))

    @classmethod
    def build_steps(cls, context, owner, group, *, source_heights=False):
        if type(source_heights) is not bool:
            raise ValueError('Original source heights require an explicit Boolean policy')
        builder = CompoundSurfaceBuilder(context.patterns, source_heights=source_heights)
        surface = yield from builder.build_steps(group, context)
        yield
        group = surface.chart.partition.group
        layouts = {p.region.owner: p.layout for p in group.footprints}
        data = surface.conditions
        if data is None:
            raise ValueError('Compound surface did not retain its actual boundary conditions')
        sources = {p.region.owner: p for p in group.footprints}
        rims, kinds, peaks, port_points, physical_rims = {}, {}, {}, {}, {}
        derived_peaks, derived_port_points = {}, {}
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
                    if alias.original_incidence:
                        peaks[alias.endpoint.section_id, alias.endpoint.end] = (li, index, alias.station)
                    else:
                        derived_peaks[owner, alias.endpoint.section_id, alias.endpoint.end] = (li, index, alias.station)
                for vertex, u in zip(port.vertices, port.outward_u):
                    xy = data[li].xy[vertex]
                    at = int(np.argmin(np.linalg.norm(rims[li][:, :2] - xy, axis=1)))
                    if np.linalg.norm(rims[li][at, :2] - xy) > 1e-7:
                        raise ValueError('Boundary output lost an incoming profile vertex')
                    for alias in aliases:
                        if alias.contains(u):
                            original_u = float(alias.outward_parameter(u))
                            if alias.original_incidence:
                                port_points[alias.endpoint.section_id, alias.endpoint.end, round(original_u, 12)] = (li, at, alias.station)
                            else:
                                derived_port_points[owner, alias.endpoint.section_id, alias.endpoint.end,
                                                    round(original_u, 12)] = (li, at, alias.station)
        identity = Key('compound', owner.id) if len(group.members) > 1 else owner
        metadata = dict(surface.metadata, member_junctions=[key.text() for key in group.members],
                        consumed_sections=sorted({c.section_id for c in group.claims}),
                        physical_exterior_ports=group.physical_layout is not None,
                        default_provider_integration=True)
        if source_heights:metadata['source_height_policy']=SOURCE_HEIGHT_POLICY
        metadata['finite_scope_ports'] = [port.report() for layout in layouts.values()
                                         for port in getattr(layout, 'finite_scope_ports', ())]
        metadata['finite_scope_source_proof'] = getattr(group.physical_layout,'finite_scope_source_proof',None)
        metadata['finite_scope_curb_sampling'] = getattr(group.physical_layout,'finite_scope_curb_sampling',())
        metadata['finite_scope_rounded_contacts'] = getattr(group.physical_layout,'finite_scope_rounded_contacts',())
        metadata['finite_free_pads'] = [claim.report() for layout in layouts.values()
                                      for claim in getattr(layout, 'finite_free_pads', ())]
        metadata['source_cycle_bundle_proofs'] = (getattr(group.physical_layout,'source_cycle_proofs',{})
                                                 if group.physical_layout is not None else {})
        metadata['source_cycle_covered_corners'] = (getattr(group.physical_layout,'source_cycle_covered_corners',[])
            if group.physical_layout is not None else [])
        metadata['source_cycle_inner_bundle_width_support'] = getattr(
            group.physical_layout,'source_cycle_inner_bundle_width_support',())
        metadata['source_cycle_covered_partial_corners'] = getattr(
            group.physical_layout,'source_cycle_covered_partial_corners',())
        metadata['source_cycle_separated_terminal_cuts'] = getattr(
            group.physical_layout,'source_cycle_separated_terminal_cuts',())
        metadata['owned_prefix_width_proof'] = getattr(group.physical_layout,'owned_prefix_width_proof',None)
        metadata['source_cycle_overlap_expansions'] = getattr(group.physical_layout,'source_cycle_overlap_expansions',[])
        metadata['source_cycle_displaced_contacts'] = getattr(group.physical_layout,'source_cycle_displaced_contacts',[])
        metadata['source_cycle_dimensional_contacts'] = getattr(group.physical_layout,'source_cycle_dimensional_contacts',[])
        metadata['source_cycle_root_envelopes'] = getattr(group.physical_layout,'source_cycle_root_envelopes',[])
        metadata['source_cycle_root_support_envelopes'] = getattr(group.physical_layout,'source_cycle_root_support_envelopes',[])
        metadata['source_cycle_source_contact_envelopes'] = getattr(group.physical_layout,'source_cycle_source_contact_envelopes',[])
        metadata['source_cycle_thin_median_gores'] = getattr(group.physical_layout,'source_cycle_thin_median_gores',[])
        metadata['extended_source_pool'] = getattr(group.physical_layout,'extended_source_pool',None)
        metadata['source_cycle_terminal_ports'] = getattr(group.physical_layout,'source_cycle_terminal_ports',[])
        metadata['source_cycle_finite_profile_search'] = getattr(group.physical_layout,'source_cycle_finite_profile_search',None)
        metadata['source_cycle_mouth_clearance'] = getattr(group.physical_layout,'source_cycle_mouth_clearance',())
        metadata['source_cycle_exposed_root_contact_retries'] = getattr(group.physical_layout,'source_cycle_exposed_root_contact_retries',())
        metadata['wholly_owned_free_counterparts'] = getattr(group.physical_layout,'wholly_owned_free_counterparts',None)
        metadata['source_cycle_mouth_expansions'] = getattr(group.physical_layout,'source_cycle_mouth_expansions',())
        metadata['source_rim_normalization'] = getattr(group.physical_layout,'source_rim_normalization',())
        metadata['source_cycle_full_width_proof'] = getattr(group.physical_layout,'source_cycle_full_width_proof',None)
        metadata['coupled_source_owner_proofs'] = getattr(group.physical_layout,'coupled_source_owner_proofs',())
        metadata['coupled_source_dependencies'] = getattr(group.physical_layout,'coupled_source_dependencies',())
        metadata['ordinary_complete_mouth_clearance'] = getattr(group.physical_layout,'ordinary_complete_mouth_clearance',None)
        metadata['partial_bundle_claims'] = [claim.report() for claim in
            getattr(group.physical_layout,'partial_bundle_claims',())]
        metadata['partial_bundle_width_proof'] = getattr(group.physical_layout,'partial_bundle_width_proof',None)
        metadata['internal_profile_fragments'] = getattr(group.physical_layout,'internal_profile_fragments',None)
        metadata['source_curb_owner_canonicalization'] = getattr(group.physical_layout,'source_curb_owner_canonicalization',())
        metadata['source_bank_reconciliation'] = getattr(group.physical_layout,'source_bank_reconciliation',None)
        bank_circles = getattr(group.physical_layout, 'native_original_circle_fillets', ())
        if bank_circles:
            metadata['bank_circle_fillets'] = [fillet.report() for fillet in bank_circles]
        metadata['source_retained_arc_policy'] = getattr(group.physical_layout,'source_retained_arc_policy',None)
        metadata['source_owned_connector_reconciliation'] = getattr(group.physical_layout,'source_owned_connector_reconciliation',None)
        metadata['source_parent_chain_reconciliation'] = getattr(group.physical_layout,'source_parent_chain_reconciliation',None)
        metadata['residual_road_gap_reconciliation'] = getattr(group.physical_layout,'residual_road_gap_reconciliation',None)
        metadata['source_root_prepared_neighbor_guard'] = getattr(group.physical_layout,'source_root_prepared_neighbor_guard',None)
        metadata['generated_curb_prepared_neighbor_guard'] = getattr(group.physical_layout,'generated_curb_prepared_neighbor_guard',None)
        from .corner_radii import fitted_radius_report
        metadata['fitted_corner_radii'] = [fit for key, layout in layouts.items()
            for fit in fitted_radius_report(getattr(layout,'original_junction',key), layout)]
        pavement = Part(identity, surface.mesh, metadata)
        grass = None
        physical_grass = None
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
            if context.network.settings.exaggeration == 1.:
                physical_grass = mesh
            else:
                physical_assembly = Assembler()
                for li, physical in surface.physical_holes.items():
                    physical_assembly.add(physical.xyz,physical.quads)
                physical_grass = physical_assembly.mesh('grass')
                physical_grass.seams = {li:ids.copy() for li,ids in seams.items()}
        return cls(owner, group, surface.field, rims, kinds, peaks, port_points, pavement, grass,
                   context.network.settings.exaggeration, physical_rims, derived_peaks, derived_port_points,
                   surface.physical_mesh,physical_grass,surface.hole_fields)

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
        hole_points = {}
        for li,loop in enumerate(self.group.loops):
            if loop.hole:
                rim = self.physical_rims[li][:,:2] if self.physical_rims is not None else loop.xy
                hit, boundary = _classify(xy, rim, 1e-7)
                selected_hole = hit | (boundary >= 0)
                hole_points[li] = selected_hole
                holes |= selected_hole
        result = np.full(len(xy), self.field.datum)
        selected = inside & ~holes
        result[selected] = self.field.datum + self.exaggeration * (self.field(xy[selected]) - self.field.datum)
        for li,selected_hole in hole_points.items():
            if not selected_hole.any():
                continue
            if self.hole_fields is None or li not in self.hole_fields:
                raise ValueError('Compound hole query is missing its actual source-derived height field')
            physical = self.hole_fields[li](xy[selected_hole])
            result[selected_hole] = self.field.datum + self.exaggeration*(physical-self.field.datum)
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
        self.physical_pavement = self.physical_grass = None


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
        source_ids = {visit.section_id for bundle in network.bundles if len(bundle.members) == 2
                      for member in bundle.members for visit in member.path}
        yield Task(PLAN, (*tuple(footprint_key(j.key) for j in junctions),
                         *tuple(bundle.key for bundle in network.bundles),
                         *tuple(Key(kind, sid) for sid in sorted(source_ids) for kind in ('path', 'profile'))),
                   CooperativeBuild(lambda context: OwnershipPlan.build(context, junctions),
                                    lambda context: OwnershipPlan.build_steps(context, junctions)),
                   version=SOURCE_SAMPLING_POLICY+':'+COUPLED_OWNER_POLICY+':'+ORDINARY_MOUTH_POLICY+':'+SOURCE_RANGE_POLICY+':'+SOURCE_CURB_COORDINATE_POLICY+':'+SOURCE_BANK_POLICY+':'+SOURCE_ARC_POLICY+':'+SOURCE_OWNED_CONNECTOR_POLICY+':'+SOURCE_PARENT_CHAIN_POLICY+':'+SOURCE_RESIDUAL_GAP_POLICY+':'+SOURCE_OWNED_BUNDLE_NORMAL_POLICY+':'+SOURCE_OWNED_PREFIX_CORNER_POLICY+':'+SOURCE_ROOT_POLICY)
        for junction in junctions:
            yield Task(surface_key(junction.key), (PLAN, junction.layout_key),
                CooperativeBuild(lambda context, site=junction.key: JunctionOwnership.surface(context, site),
                                 lambda context, site=junction.key: JunctionOwnership.surface_steps(context, site)),
                signature=lambda context, site=junction.key: JunctionOwnership.surface_signature(context, site),
                version='2-independent-vertex-quad-repair')
            yield Task(Key('junction_grass', junction.id), (surface_key(junction.key),),
                lambda context, site=junction.key: context[surface_key(site)].grass
                    if context[surface_key(site)] is not None else None)
        yield Task(SURFACES_READY, tuple(surface_key(j.key) for j in junctions), lambda context: None)

    @staticmethod
    def surface_signature(context, site):
        from .hole_surfaces import SOURCE_HOLE_HEIGHT_POLICY
        plan = context[PLAN]
        owner = plan.owner(site)
        if owner is None or owner != site:
            return SOURCE_SAMPLING_POLICY, SOURCE_RANGE_POLICY, SOURCE_CURB_COORDINATE_POLICY, SOURCE_BANK_POLICY, SOURCE_ARC_POLICY, SOURCE_OWNED_CONNECTOR_POLICY, SOURCE_PARENT_CHAIN_POLICY, SOURCE_RESIDUAL_GAP_POLICY, SOURCE_ROOT_POLICY, SOURCE_OWNED_BUNDLE_NORMAL_POLICY, SOURCE_OWNED_PREFIX_CORNER_POLICY, SOURCE_HEIGHT_POLICY, SOURCE_HOLE_HEIGHT_POLICY, (owner.text() if owner is not None else context.digest(Key('mouths', site.id)))
        group=plan.by_site[site]
        return (SOURCE_SAMPLING_POLICY, SOURCE_RANGE_POLICY, SOURCE_CURB_COORDINATE_POLICY, SOURCE_BANK_POLICY, SOURCE_ARC_POLICY, SOURCE_OWNED_CONNECTOR_POLICY, SOURCE_PARENT_CHAIN_POLICY, SOURCE_RESIDUAL_GAP_POLICY, SOURCE_ROOT_POLICY, SOURCE_OWNED_BUNDLE_NORMAL_POLICY, SOURCE_OWNED_PREFIX_CORNER_POLICY, SOURCE_HEIGHT_POLICY, SOURCE_HOLE_HEIGHT_POLICY, [(key.text(), context.digest(footprint_key(key))) for key in group.members],
                bundle_signature(context,group))

    @staticmethod
    def surface(context, site):
        return drain_steps(JunctionOwnership.surface_steps(context, site))

    @staticmethod
    def surface_steps(context, site):
        plan = context[PLAN]
        owner = plan.owner(site)
        if owner is not None:
            if owner != site:
                return None
            return (yield from CompoundProduct.build_steps(context, site, plan.by_site[site], source_heights=True))
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
        return (yield from CompoundProduct.build_steps(context, site, group, source_heights=True))

    @staticmethod
    def layout_signature(context, junction):
        plan = context[PLAN]
        owner = plan.owner(junction.key)
        group = plan.by_site.get(junction.key)
        if group is not None and group.physical_layout is not None:
            return (SOURCE_SAMPLING_POLICY, SOURCE_RANGE_POLICY, SOURCE_CURB_COORDINATE_POLICY, SOURCE_BANK_POLICY, SOURCE_ARC_POLICY, SOURCE_OWNED_CONNECTOR_POLICY, SOURCE_PARENT_CHAIN_POLICY, SOURCE_RESIDUAL_GAP_POLICY, SOURCE_ROOT_POLICY, SOURCE_OWNED_BUNDLE_NORMAL_POLICY, SOURCE_OWNED_PREFIX_CORNER_POLICY, 'exterior_ports', [(key.text(), context.digest(footprint_key(key))) for key in group.members],
                    bundle_signature(context,group))
        return SOURCE_SAMPLING_POLICY, SOURCE_RANGE_POLICY, SOURCE_CURB_COORDINATE_POLICY, SOURCE_BANK_POLICY, SOURCE_ARC_POLICY, SOURCE_OWNED_CONNECTOR_POLICY, SOURCE_PARENT_CHAIN_POLICY, SOURCE_RESIDUAL_GAP_POLICY, SOURCE_ROOT_POLICY, SOURCE_OWNED_BUNDLE_NORMAL_POLICY, SOURCE_OWNED_PREFIX_CORNER_POLICY, context.digest(footprint_key(junction.key)), owner.text() if owner is not None else None

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
