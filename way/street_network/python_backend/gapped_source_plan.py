"""Qualified union plan for actual separated side-junction source spans.

This module exposes an explicit last-family API. The provider chooses when to
call it; importing it neither discovers nor installs a new source preparation.
"""
import numpy as np
from .bundle_merged_plan import physical_context
from .gapped_source_inventory import GappedSourceInventory, GappedSourceUnavailable,source_ledger
from .gapped_source_body import GappedSourceBodyPiece
from .mesh import Assembler, Mesh
from .model import Key
from .ownership import surface_key
from .providers import Part, CarriagewayProvider
from .routes import RouteTopology
from .surface_query import QuadHeightField
from .gapped_runtime_admission import visible_overlap
from .bundle_axis import AxialChartUnavailable
from .ownership import CurbChainUnavailable


def _coverage(context, inventory):
    prefixes,_=source_ledger(context,inventory.run,inventory.spans_by_role,inventory.owned_spans)
    attachments = RouteTopology(context.network).attachments(context.network)
    active = {row['section']: [dict(row, owner='new_source_body')]
              for own in inventory.spans_by_role for span in own for row in span.source_intervals}
    for span in inventory.owned_spans:
        for row in span.intervals:
            active[row['section']] = [dict(row, owner=span.owner.text(), complete_owned_source=True)]
    for row in prefixes:
        active.setdefault(row['section'],[]).append(dict(section=row['section'],forward=row['forward'],
            first=row['first'],last=row['last'],owner=row['owner'].text(),
            complete_original_route_prefix=True))
    reports = []
    for member in inventory.run.members:
        for visit in member.visits:
            sid = visit.section_id; path = context[Key('path', sid)]
            metadata = CarriagewayProvider.span_metadata(context, context.network.section_index[sid], attachments[sid])
            owned = [dict(row, owner=kind) for kind in ('compound_spans', 'junction_spans', 'transition_spans')
                     for row in metadata.get(kind, ())]
            rows = active[sid]
            pieces = [(row['first'], row['last']) for row in rows]
            pieces.extend((row['start'], row['end']) for row in owned)
            last = 0.
            for first, tail in sorted(pieces):
                if tail <= first+1e-9: continue
                if not 0. <= first < tail <= path.length+1e-7 or first > last+1e-7:
                    raise ValueError('A gapped source plan drops an original finite source interval')
                last = max(last, tail)
            if abs(last-path.length) > 1e-7:
                raise ValueError('A gapped source plan does not cover its complete original source extent')
            reports.append(dict(section=sid, forward=visit.forward, body_or_whole_owner_intervals=rows,
                actual_original_owner_prefixes=owned, original_length=path.length,
                complete_original_extent=True,
                width_scope='visible original intervals, wholly owned Sections and certified original junction route prefixes'))
    if len(reports) != len({row['section'] for row in reports}):
        raise ValueError('A gapped source plan repeats an original source Section')
    return reports


class GappedSourcePlan:
    mode = 'union'
    requires_separator_query = True
    transition_owners = ()

    @classmethod
    def build(cls, context, bundle, *, require_visible_overlap=False):
        if type(require_visible_overlap) is not bool:
            raise ValueError('Gapped source admission flag must be Boolean')
        if len(bundle.members) != 2 or not all(member.complete for member in bundle.members):
            raise GappedSourceUnavailable('A gapped source plan needs two complete original members')
        sections = set(bundle.sections)
        required = {site.key for site in context.network.intersections
                    if any(end.section_id in sections for end in site.ends)}
        if not any(len(site.ends) >= 3 and sum(end.section_id in sections for end in site.ends) == 2
                   for site in context.network.intersections):
            raise GappedSourceUnavailable('A gapped source plan needs an actual original one-sided junction')
        try:
            physical = physical_context(context, required_owners=required)
            inventory = GappedSourceInventory.build(physical, bundle)
            overlap = visible_overlap(physical, inventory) if require_visible_overlap else None
        except (AxialChartUnavailable, CurbChainUnavailable) as unavailable:
            if not require_visible_overlap or isinstance(unavailable, GappedSourceUnavailable):
                raise
            # Discovery may decline a typed unsupported inventory/chart. Once
            # admitted, coverage, source certificates and actual meshing errors
            # remain fatal; their checks are deliberately outside this block.
            raise GappedSourceUnavailable('Gapped runtime source chart is unavailable: '+str(unavailable)) from unavailable
        if require_visible_overlap and overlap is None:
            raise GappedSourceUnavailable('Gapped runtime replacement needs certified visible member overlap')
        source_coverage = _coverage(physical, inventory)
        piece = GappedSourceBodyPiece()
        pavement, grass, report = piece.build(physical, inventory, namespace=(bundle.key, 'gapped source body'))
        return cls(context, physical, bundle, inventory, piece, pavement, grass,
                   dict(report, original_source_coverage=source_coverage,
                        runtime_discovery_enabled=require_visible_overlap,
                        runtime_overlap_admission=overlap))

    def __init__(self, context, physical, bundle, inventory, piece, pavement, grass, report):
        self.settings, self.bundle = context.network.settings, bundle
        self.field, self.frame = inventory.field, inventory.frame
        self.datum = self.field.datum
        self.source_ids = tuple(v.section_id for member in bundle.members for v in member.path)
        covered = [v.section_id for member in inventory.run.members for v in member.visits]
        if len(covered) != len(set(covered)) or set(covered) != set(self.source_ids):
            raise ValueError('A gapped union plan must retain every original member Section exactly once')
        self.transition_owners = tuple(sorted({owner for own in inventory.spans_by_role for span in own
                                               for owner in span.transition_owners}))
        owners = {node.owner for own in inventory.spans_by_role for span in own for node in (span.head, span.tail)}
        owners.update(span.owner for span in inventory.owned_spans)
        query = Assembler()
        query.add(pavement.xyz, pavement.quads); query.add(grass.xyz, grass.quads)
        ceiling = self.field.ceiling
        for owner in sorted(owners):
            product = physical[surface_key(owner)]
            if abs(product.field.datum-self.datum) > 1e-9:
                raise ValueError('A gapped source query needs its original common minimum datum')
            ceiling = max(ceiling, product.field.ceiling)
            for part in (product.pavement, product.grass):
                if part is None: continue
                if part.mesh is None: raise ValueError('A gapped source query lost an actual physical owner buffer')
                query.add(part.mesh.xyz, part.mesh.quads)
        self.query = QuadHeightField(query.mesh(), datum=self.datum, ceiling=ceiling, max_cells=200000)
        self.profile_ids = piece.profile_ids
        self.inventory = inventory
        metadata = dict(report, source_sections=list(self.source_ids),
            complete_selected_source_coverage=True,
            query_only_junction_owners=[owner.text() for owner in sorted(owners)],
            junction_geometry_republished=False, all_direct_quads=True)
        self._pavement = Part(bundle.key, self._display(pavement), dict(metadata, quality=pavement.validate()))
        self._separator = Part(bundle.key, self._display(grass), dict(metadata, quality=grass.validate()))

    def _display(self, mesh):
        xyz = mesh.xyz.copy()
        xyz[:, 2] = self.datum+self.settings.exaggeration*(xyz[:, 2]-self.datum)
        return Mesh(xyz, mesh.quads.copy(), mesh.material)

    def pavement(self, context): return self._pavement

    def separator(self, context): return self._separator

    def height(self, xy):
        xy = np.asarray(xy)
        heights = self.query.evaluate(xy.reshape(-1, 2), default=self.datum).reshape(xy.shape[:-1])
        return self.datum+self.settings.exaggeration*(heights-self.datum)

    def release_geometry(self):
        for part in (self._pavement, self._separator):
            if part.mesh is not None:
                part.metadata['bounds'] = part.mesh.bounds
                part.mesh = None
