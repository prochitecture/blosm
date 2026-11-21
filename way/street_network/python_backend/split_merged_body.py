"""Two bounded direct-quad bodies sharing one computational merged fibre.

Only the original body's explicit combined-cell failure admits this family.
Both child bodies retain their existing eight-thousand-cell construction cap;
the complete source leaf and its connectors retain their sixteen-thousand cap.
"""
from types import SimpleNamespace
import numpy as np
from .bundle_common_fibres import CommonFibreCuts
from .bundle_merged_body import MergedRunBody, MergedBodyCellBudgetError
from .bundle_angular_open_islands import gap_intervals
from .drainage import BoundaryDrainage
from .four_side_charts import quality
from .mesh import Assembler, Mesh
from .model import Key
from .validation import validate_overlaps


class BoundedMergedRunBody(MergedRunBody):
    allow_separated_profiles = False
    def build_between(self, context, bundle, approach, cuts, *, namespace):
        try:
            return super().build_between(context, bundle, approach, cuts, namespace=namespace)
        except MergedBodyCellBudgetError as original:
            return self._split(context, bundle, approach, cuts, namespace, original)

    def _split(self, context, bundle, approach, cuts, namespace, original):
        frame, field = approach.frame_chart, approach.field
        prepared = CommonFibreCuts(frame, field, cuts, namespace)
        if not self.allow_separated_profiles and tuple(prepared.merged_ends) != (0, 1):
            raise ValueError('Bounded two-body recovery needs two actual merged outer cuts') from original
        first, last = prepared.interval
        domain = SimpleNamespace(frame=frame, field=field, rays=field.curves)
        regions, roots = gap_intervals(domain, first, last)
        reserve = 2.*context.network.settings.spacing
        allowed = [(index, a, b) for index, (a, b, positive) in enumerate(regions)
                   if not positive and b-a > 2.*reserve]
        if not allowed and self.allow_separated_profiles:
            allowed = [(index, a, b) for index, (a, b, positive) in enumerate(regions)
                       if b-a > 2.*reserve]
        if not allowed:
            raise ValueError('Two bounded source bodies need a finite merged computational fibre') from original
        fractions = (.5, .45, .55, .4, .6) if self.allow_separated_profiles else (.5,)
        tried, failures = set(), []
        for fraction in fractions:
            middle = first+fraction*(last-first) if fraction != .5 else .5*(first+last)
            candidates = [(abs(float(np.clip(middle, a+reserve, b-reserve))-middle), index,
                           float(np.clip(middle, a+reserve, b-reserve))) for index, a, b in allowed]
            station = min(candidates)[2]
            if station in tried:continue
            tried.add(station)
            try:
                result = self._at_station(context, bundle, approach, cuts, namespace,
                    original, first, last, reserve, station)
                if self.allow_separated_profiles:
                    result[2]['separated_source_computation'] = dict(
                        original_outer_merged_ends=list(prepared.merged_ends),
                        prior_budget_failures=failures, station_fraction=fraction)
                return result
            except MergedBodyCellBudgetError as failure:
                if not self.allow_separated_profiles:raise
                failures.append(dict(station=station, quads=failure.quads, limit=failure.limit))
        raise ValueError('No bounded two-body station meets the unchanged child cell budgets: '+str(failures)) from original

    def _at_station(self, context, bundle, approach, cuts, namespace, original,
                    first, last, reserve, station):
        frame, field = approach.frame_chart, approach.field
        if min(station-first, last-station) <= reserve:
            raise ValueError('Two bounded source bodies leave no positive computational subinterval') from original
        ao, ai, bi, bo = map(float, field.edges(station))
        high, low = map(float, field.peaks(station))
        merged = ai < bi and abs(high-low) <= 1e-9
        if not merged and not self.allow_separated_profiles:
            raise ValueError('The new computational cut would change a genuine separator') from original
        crests = [np.r_[frame.world(np.array((station, peak))), field.ceiling] for peak in (high, low)]
        if merged:crests[1] = crests[0]
        midpoint = {}
        for role, outer in enumerate((ao, bo)):
            count = max(len(cuts[role, end])-1 for end in (0, 1))
            values = np.linspace(outer, (high, low)[role], count+1)
            xy = frame.world(np.c_[np.full(len(values), station), values])
            midpoint[role] = np.c_[xy, field.values(xy)]
            midpoint[role][-1] = crests[role]
            midpoint[role][0, 2] = field.datum
        paving, grassy = Assembler(), Assembler()
        children, reports = [], []
        middle_keys = {}
        for role in (0, 1):
            middle_keys[role] = [(namespace, 'internal fibre', role, i) for i in range(len(midpoint[role]))]
            middle_keys[role][-1] = ((namespace, 'internal fibre crest') if merged
                                    else (namespace, 'internal fibre crest', role))
        grass_cut = None
        grass_keys = None
        for part in (0, 1):
            local = {(role, end): (cuts[role, 0] if part == end == 0 else
                     cuts[role, 1] if part == end == 1 else midpoint[role])
                     for role in (0, 1) for end in (0, 1)}
            body = MergedRunBody()
            pavement, grass, report = body.build_between(context, bundle, approach, local,
                                                       namespace=(namespace, 'bounded body', part))
            keys = list(body.pavement_keys)
            at = 1-part
            aliases = {}
            for role in (0, 1):
                for key, target in zip(body.cut_keys[role, at], middle_keys[role]): aliases[key] = target
            keys = [aliases.get(key, key) for key in keys]
            paving.add(pavement.xyz, pavement.quads, keys)
            if grass is not None:
                names = list(body.grass_keys)
                if not merged:
                    if at not in body.open_grass_cuts:
                        raise ValueError('A separated computational fibre lost its whole grass cross-section')
                    points, cut_names = body.open_grass_cuts[at]
                    if grass_cut is None:
                        grass_cut = points.copy()
                        grass_keys = [(namespace, 'internal grass fibre', i) for i in range(len(points))]
                        for i, point in enumerate(points):
                            for role in (0, 1):
                                if np.array_equal(point, midpoint[role][-1, :2]):
                                    grass_keys[i] = middle_keys[role][-1]
                        if len({middle_keys[role][-1] for role in (0, 1)} & set(grass_keys)) != 2:
                            raise ValueError('The computational grass fibre lost its actual two profile peaks')
                    if np.array_equal(points, grass_cut):
                        targets = grass_keys
                    elif np.array_equal(points, grass_cut[::-1]):
                        targets = grass_keys[::-1]
                    else:
                        raise ValueError('The two bounded bodies need one identical complete grass fibre')
                    aliases.update(zip(cut_names, targets))
                names = [aliases.get(key, key) for key in names]
                grassy.add(grass.xyz, grass.quads, names)
            children.append(body); reports.append(report)
        pavement = paving.mesh()
        grass = grassy.mesh('grass') if grassy.quads else None
        limit = 16000*context.network.settings.resolution
        total = len(pavement.quads)+(len(grass.quads) if grass is not None else 0)
        if total >= limit:
            raise ValueError('The two source bodies exceed the unchanged complete source-leaf budget')
        edge_ids = np.sort(np.stack((pavement.quads, np.roll(pavement.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edge_ids, axis=0, return_counts=True)
        incidence = {tuple(edge): int(count) for edge, count in zip(unique, counts)}
        shared = 0
        for names in middle_keys.values():
            nodes = [paving.lookup[key] for key in names]
            for a, b in zip(nodes, nodes[1:]):
                if incidence.get(tuple(sorted((a, b)))) != 2:
                    raise ValueError('The two bounded bodies lack one exact whole common-crown edge')
                shared += 1
        grass_shared = 0
        if grass_keys is not None:
            edge_ids = np.sort(np.stack((grass.quads, np.roll(grass.quads, -1, axis=1)),
                                        axis=-1).reshape(-1, 2), axis=1)
            unique, counts = np.unique(edge_ids, axis=0, return_counts=True)
            incidence = {tuple(edge): int(count) for edge, count in zip(unique, counts)}
            nodes = [grassy.lookup[key] for key in grass_keys]
            for a, b in zip(nodes, nodes[1:]):
                if incidence.get(tuple(sorted((a, b)))) != 2:
                    raise ValueError('The two bounded bodies lack one exact whole computational grass edge')
                grass_shared += 1
        checks = []
        payloads = {Key('bounded_source_body', 0): pavement}
        if grass is not None: payloads[Key('bounded_source_grass', 0)] = grass
        for mesh in payloads.values():
            stored = Mesh(mesh.xyz.astype(np.float32).astype(float), mesh.quads, mesh.material)
            q, published = quality(mesh), quality(stored)
            if (max(q['maximum_edge_ratio'], published['maximum_edge_ratio']) > 8.
                    or mesh.material != 'grass' and min(q['minimum_corner_sine'], published['minimum_corner_sine']) < .03):
                raise ValueError('The assembled bounded source bodies fail physical or stored grading')
            slope = max(mesh.slope_statistics()['maximum'], stored.slope_statistics()['maximum'])
            along = max(frame.longitudinal_bound(mesh), frame.longitudinal_bound(stored))
            if slope > .1+1e-6 or mesh.material != 'grass' and along > .03+1e-8:
                raise ValueError('The assembled bounded source bodies exceed physical slope limits')
            fixed = np.zeros(len(mesh.xyz), bool); fixed[np.unique(self._boundary_vertices(mesh))] = True
            drainage = BoundaryDrainage(max_cells=limit).fill(stored, fixed)[1]
            if drainage['maximum_fill'] > 1e-6: raise ValueError('The assembled bounded source bodies create a drainage depression')
            checks.append(dict(material=mesh.material, physical_quality=q, stored_quality=published,
                               maximum_slope=slope, longitudinal_bound=along, drainage=drainage))
        self.max_cells = limit
        self.pavement_keys = [None]*len(pavement.xyz)
        for key, at in paving.lookup.items(): self.pavement_keys[at] = key
        self.grass_keys = [None]*len(grass.xyz) if grass is not None else []
        for key, at in grassy.lookup.items(): self.grass_keys[at] = key
        self.cuts = {(role, end): children[end].cuts[role, end] for role in (0, 1) for end in (0, 1)}
        self.cut_keys = {(role, end): children[end].cut_keys[role, end] for role in (0, 1) for end in (0, 1)}
        self.open_grass_cuts = {end: children[end].open_grass_cuts[end] for end in (0, 1)
                               if end in children[end].open_grass_cuts}
        overlap = validate_overlaps(payloads, cell_size=2.)
        report = dict(construction='two bounded source bodies after one-body cell-budget failure',
            original_cell_budget=original.limit, original_candidate_quads=original.quads,
            child_cell_budgets=[child.max_cells for child in children], complete_leaf_limit=limit,
            complete_quads=total, computational_fibre=station, exact_shared_crown_edges=shared,
            physical_source_changed=False, physical_profiles_changed=False, source_intervals_changed=False,
            children=reports, stored_checks=checks, overlap=overlap)
        if self.allow_separated_profiles:
            report['computational_fibre_merged'] = merged
            report['exact_shared_grass_edges'] = grass_shared
        return pavement, grass, report


class BoundedGappedRunBody(BoundedMergedRunBody):
    """Only an explicit gapped inventory admits separated computational cuts."""
    allow_separated_profiles = True
