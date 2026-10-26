"""Complete merged-entrance Bundle output in the normal preparation graph."""
from dataclasses import replace
import numpy as np
from .model import Key
from .updates import BuildContext
from .ownership import PLAN, CompoundProduct, surface_key
from .bundle_merged_entrances import PairedAngularEntrance
from .bundle_merged_bridge import MergedProfileBridge
from .bundle_merged_runs import MergedRunWithExternalCap
from .mesh import Assembler
from .providers import Part
from .surface_query import QuadHeightField


def physical_context(context):
    """Reuse prepared junctions and their exact physical rims, even at E0.

    No junction, path or source profile is recomputed. Only the small product
    views needed by cap sampling differ from the published inspection views.
    """
    products = dict(context.products)
    for key, product in products.items():
        if isinstance(product, CompoundProduct):
            if product.physical_rims is None:
                raise ValueError('Merged Bundle cap requires its prepared physical rim')
            products[key] = replace(product, rims=product.physical_rims, exaggeration=1.)
    return BuildContext(replace(context.network, settings=replace(context.network.settings, exaggeration=1.)),
                        products, context.patterns, context.task_digests, context.source_digests)


class MergedBundlePlan:
    mode = 'union'
    transition_owners = ()
    requires_separator_query = True

    def __init__(self, context, bundle, specification):
        physical = physical_context(context)
        owner = min(specification.crossing.sites)
        product = physical[surface_key(owner)]
        if not isinstance(product, CompoundProduct):
            raise ValueError('Merged Bundle needs its completed shared physical junction')
        self.bundle, self.settings = bundle, context.network.settings
        self.source_ids = tuple(visit.section_id for member in bundle.members for visit in member.path)
        self.datum = product.field.datum
        paving, grass, reports = Assembler(), Assembler(), []
        for index, approach in enumerate(specification.layout.approaches):
            if not isinstance(approach, PairedAngularEntrance):
                continue
            station = specification.layout.outward_stations[index]
            identities = {(end.section_id, end.end) for end in approach.original_ends}
            refs = sorted({(loop, node) for (sid, end, u), (loop, node, at) in product.port_points.items()
                           if (sid, end) in identities})
            mouth = np.asarray([product.rims[loop][node] for loop, node in refs])
            if len(mouth) < 3:
                raise ValueError('Merged Bundle lost its complete physical mouth')
            first = approach.point(station, -.5*approach.width)
            last = approach.point(station, .5*approach.width)
            mouth = mouth[np.argsort((mouth[:, :2]-first)@(last-first))]
            bridge = MergedProfileBridge()
            connector, connector_report = bridge.build(physical, approach, station, mouth)
            run = MergedRunWithExternalCap()
            pavement, separator, run_report = run.build(physical, bundle, approach, bridge.body_cut)
            joined = Assembler()
            keys = [('connector', i) for i in range(len(connector.xyz))]
            for row, node in enumerate(bridge.back_ids):
                keys[node] = ('merged cut', row)
            joined.add(connector.xyz, connector.quads, keys)
            joined.add(pavement.xyz, pavement.quads, run.pavement_keys)
            mesh = joined.mesh()
            paving.add(mesh.xyz, mesh.quads)
            grass.add(separator.xyz, separator.quads)
            reports.append(dict(source_sections=list(run.domain.run.source_sections),
                                connector=connector_report, body=run_report))
        if len(reports) != 2:
            raise ValueError('Merged Bundle needs its two complete source-owned runs')
        physical_paving, physical_grass = paving.mesh(), grass.mesh('grass')
        total = len(physical_paving.quads)+len(physical_grass.quads)
        if total > 16000*self.settings.resolution:
            raise ValueError('Merged Bundle exceeds its complete finite cell budget')
        query = Assembler()
        query.add(physical_paving.xyz, physical_paving.quads)
        query.add(physical_grass.xyz, physical_grass.quads)
        self.query = QuadHeightField(query.mesh(), datum=self.datum, ceiling=product.field.ceiling)
        for mesh in (physical_paving, physical_grass):
            mesh.xyz[:, 2] = self.datum+self.settings.exaggeration*(mesh.xyz[:, 2]-self.datum)
        metadata = dict(construction='source-owned merged Bundle entrances and complete rounded runs',
            source_sections=list(self.source_ids), member_junctions=[key.text() for key in specification.crossing.sites],
            consumed_sections=list(specification.crossing.sections), runs=reports,
            centerline_proof_checks=specification.centerline_proof_checks,
            complete_quads=total, all_direct_quads=True, default_provider_integration=True)
        self._pavement = Part(bundle.key, physical_paving, dict(metadata, quality=physical_paving.validate()))
        self._separator = Part(bundle.key, physical_grass, dict(metadata, quality=physical_grass.validate()))

    def pavement(self, context):
        return self._pavement

    def separator(self, context):
        return self._separator

    def height(self, xy):
        xy = np.asarray(xy)
        shape = xy.shape[:-1]
        value = self.query.evaluate(xy.reshape(-1, 2), default=self.datum)
        return (self.datum+self.settings.exaggeration*(value-self.datum)).reshape(shape)

    def release_geometry(self):
        for part in (self._pavement, self._separator):
            if part.mesh is not None:
                part.metadata['bounds'] = part.mesh.bounds
                part.mesh = None
