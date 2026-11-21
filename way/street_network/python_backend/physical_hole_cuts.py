"""Bounded computational cuts for a certified physical source rim.

Parallel pavement intervals open each hole on both sides. Visibility cuts then
divide the resulting simple blocks. Every endpoint addresses an unchanged
physical boundary chord; no source connection or replacement rim is inferred.
"""
import math
from numbers import Real
import numpy as np
from .compounds import rim_label
from .decomposition import VisibilitySplitter
from .mesh import cross
from .model import Key
from .partition import BoundaryRef, CompoundPartitioner, CutPolicy, JunctionCut


class _FixedCuts(CutPolicy):
    def __init__(self, values):
        self.values = tuple(values)

    def cuts(self, group, context):
        return self.values


class PhysicalHoleCutError(ValueError):
    """The finite computational cut families cannot partition this fixed rim."""


class PhysicalHoleCuts(CutPolicy):
    def __init__(self, *, spacing=.7, maximum_holes=8, maximum_cuts=128, source_corridors=False,
                 hole_bridges=False, existing_vertices=False, allow_filled=False, computational_extensions=False,
                 near_endpoint_fraction=0., source_ranges=None, cadence_multiplier=1,
                 maximum_unpaired_span_widths=None):
        if (type(maximum_holes) is not int or maximum_holes < 1
                or type(maximum_cuts) is not int or maximum_cuts < 2
                or not math.isfinite(spacing) or spacing <= 0.):
            raise ValueError('Invalid physical-hole cut budgets')
        self.maximum_holes, self.maximum_cuts = maximum_holes, maximum_cuts
        self.spacing = float(spacing)
        self.source_corridors = bool(source_corridors)
        self.hole_bridges = bool(hole_bridges)
        self.existing_vertices = bool(existing_vertices)
        if (isinstance(near_endpoint_fraction, bool) or not isinstance(near_endpoint_fraction, Real)
                or not math.isfinite(near_endpoint_fraction) or not 0. <= near_endpoint_fraction <= .125
                or self.existing_vertices and near_endpoint_fraction != 0.):
            raise ValueError('Invalid bounded near-end computational cut fraction')
        self.near_endpoint_fraction = float(near_endpoint_fraction)
        self.endpoint_reuse_reports = []
        self.allow_filled = bool(allow_filled)
        if type(computational_extensions) is not bool:
            raise ValueError('Computational cut extension needs an explicit boolean policy')
        self.computational_extensions=computational_extensions
        self.extension_reports=[]
        if type(cadence_multiplier) is not int or cadence_multiplier not in (1,2):
            raise ValueError('Source-normal cadence must use one of the two bounded policies')
        if (maximum_unpaired_span_widths is not None and
                (isinstance(maximum_unpaired_span_widths,bool)
                 or not isinstance(maximum_unpaired_span_widths,Real)
                 or not math.isfinite(maximum_unpaired_span_widths)
                 or not 0.<maximum_unpaired_span_widths<=8.)):
            raise ValueError('Unpaired computational spans need a finite at-most-eight-width bound')
        self.cadence_multiplier=cadence_multiplier
        self.maximum_unpaired_span_widths=maximum_unpaired_span_widths
        self.source_ranges=None
        if source_ranges is not None:
            if not isinstance(source_ranges,dict) or not 1<=len(source_ranges)<=20000:
                raise ValueError('Source-normal intervals need a bounded original Section mapping')
            owned={};count=0
            for sid,rows in source_ranges.items():
                if type(sid) is not int or sid<0:
                    raise ValueError('Source-normal intervals need original Section IDs')
                rows=tuple(tuple(row) for row in rows);count+=len(rows)
                if (not rows or count>20000 or any(len(row)!=2 for row in rows)
                        or any(not all(math.isfinite(v) for v in row) or not 0.<=row[0]<row[1] for row in rows)
                        or any(a[1]>b[0] for a,b in zip(rows,rows[1:]))):
                    raise ValueError('Source-normal intervals must be finite ordered disjoint original ranges')
                owned[sid]=rows
            self.source_ranges=owned
        self.source_cut_reports=[]
        if self.hole_bridges and not self.source_corridors:
            raise ValueError('Physical-hole source bridges require original source corridors')

    @staticmethod
    def eligible(group):
        return (group.physical_layout is not None and len(group.loops) > 1
                and sum(not loop.hole for loop in group.loops) == 1)

    @staticmethod
    def filled_eligible(group):
        layout = group.physical_layout
        if layout is None or len(group.loops) != 1 or group.loops[0].hole:
            return False
        fragments = getattr(layout, 'internal_profile_fragments', {})
        return (hasattr(layout, 'source_cycle_proofs')
                or (getattr(layout, 'source_region_union', False) and isinstance(fragments, dict)
                    and bool(fragments.get('complete_source_width'))))

    @staticmethod
    def reference(group, hit):
        _, loop, edge, fraction = hit
        return BoundaryRef(loop, edge, fraction)

    def endpoint(self, group, hit):
        if not self.existing_vertices:
            if self.near_endpoint_fraction == 0.:
                return self.reference(group, hit)
            _, loop, edge, fraction = hit
            # One eighth is the reciprocal of the unchanged 8:1 edge guard.
            # Only a generated cut contact moves along this same curb chord;
            # exact original vertices and all physical/profile samples stay.
            if not 0. < min(fraction, 1.-fraction) < self.near_endpoint_fraction:
                return self.reference(group, hit)
            polygon = group.loops[loop].xy
            target = (edge+int(fraction > .5)) % len(polygon)
            point = polygon[edge]+fraction*(polygon[(edge+1) % len(polygon)]-polygon[edge])
            self.endpoint_reuse_reports.append(dict(loop=loop, edge=edge, fraction=fraction,
                new_edge=target, displacement=float(np.linalg.norm(polygon[target]-point))))
        _, loop, edge, fraction = hit
        # Only the computational cut moves. The original physical polygon
        # and every prescribed profile vertex remain exactly unchanged.
        return BoundaryRef(loop, (edge+int(fraction > .5)) % len(group.loops[loop].xy), 0.)

    def cuts(self, group, context):
        if not (self.eligible(group) or (self.allow_filled and self.filled_eligible(group))):
            raise ValueError('Physical-hole cuts need one certified exterior and physical holes')
        holes = [loop for loop in group.loops if loop.hole]
        if len(holes) > self.maximum_holes:
            raise ValueError('Physical-hole partition exceeds its hole budget')
        exterior = next(loop.xy for loop in group.loops if not loop.hole)
        center = exterior.mean(axis=0)
        _, vectors = np.linalg.eigh((exterior-center).T @ (exterior-center))
        principal = vectors[:, -1]
        sources = {p.region.owner: p for p in group.footprints}
        labels = [[rim_label(span, sources) for span in loop.spans] for loop in group.loops]
        failures = []
        # A mouth row is prescribed. Alternative scan orientations can reach
        # owned curbs without inserting a new vertex into an incoming profile.
        for angle in (0., math.pi*.5, math.pi*.25, -math.pi*.25):
            c, s = math.cos(angle), math.sin(angle)
            axis = np.array(((c, -s), (s, c))) @ principal
            if axis[int(np.argmax(abs(axis)))] < 0.:
                axis = -axis
            try:
                self.extension_reports=[]
                self.endpoint_reuse_reports=[]
                return self._cuts(group, context, holes, center, axis, labels)
            except ValueError as exc:
                failures.append(str(exc))
        raise PhysicalHoleCutError('Physical-hole partition found no bounded compatible cuts: ' + failures[-1])

    def _cuts(self, group, context, holes, center, axis, labels):
        normal = np.array((-axis[1], axis[0]))
        projected = [(loop.xy-center) @ axis for loop in group.loops]
        stations = []
        exterior = next(along for loop, along in zip(group.loops, projected) if not loop.hole)
        low, high = float(exterior.min()), float(exterior.max())
        length = high-low
        width = sum(loop.area for loop in group.loops)/length
        bands = int(math.ceil(length/max(width, 4.*self.spacing)))
        if bands > 64:
            raise ValueError('Physical-hole partition exceeds its longitudinal band budget')
        if bands > 2 and not self.hole_bridges:
            stations.extend(np.linspace(low, high, bands+1)[1:-1])
        for loop in holes:
            along = (loop.xy-center) @ axis
            low, high = float(along.min()), float(along.max())
            if high-low <= 1e-7:
                raise ValueError('A physical hole has no resolved scan extent')
            stations.extend((low+.25*(high-low), low+.75*(high-low)))
        owner = group.members[0]
        cuts = []
        for station in sorted(set(stations)):
            hits = []
            for li, (loop, along) in enumerate(zip(group.loops, projected)):
                following = np.roll(along, -1)
                selected = ((along <= station) & (station < following)) | (
                            (following <= station) & (station < along))
                for edge in np.flatnonzero(selected):
                    fraction = float((station-along[edge])/(following[edge]-along[edge]))
                    point = loop.xy[edge] + fraction*(loop.xy[(edge+1) % len(loop.xy)]-loop.xy[edge])
                    hits.append((float((point-center) @ normal), li, int(edge), fraction))
            hits.sort()
            if len(hits) % 2:
                raise ValueError('Physical-hole scan has an unresolved crossing')
            for first, last in zip(hits[::2], hits[1::2]):
                if last[0]-first[0] <= 1e-7:
                    raise ValueError('Physical-hole scan has a collapsed pavement interval')
                if any(labels[hit[1]][hit[2]].kind != 'curb' for hit in (first, last)):
                    continue
                cuts.append(JunctionCut(None, owner, owner,
                    self.endpoint(group, first), self.endpoint(group, last)))
        if self.source_corridors:
            source = self._source_cuts(group, context, labels)
            cuts = cuts+source if self.hole_bridges else source
        # A physical hole can taper to a sharp tip. Bisect its pavement-side
        # reflex angle before charting, rather than make a radial ring collapse
        # against that tip. Intersections with other computational cuts remain
        # ordinary conforming partition nodes.
        for li, loop in enumerate(group.loops):
            if not loop.hole:
                continue
            edge = np.roll(loop.xy, -1, axis=0)-loop.xy
            unit = edge/np.linalg.norm(edge, axis=1)[:, None]
            turns = np.arctan2(cross(np.roll(unit, 1, axis=0), unit),
                              np.sum(np.roll(unit, 1, axis=0)*unit, axis=1))
            for vertex in np.flatnonzero(turns < -.7):
                tangent = unit[vertex-1]+unit[vertex]
                direction = np.array((-tangent[1], tangent[0]))
                size = np.linalg.norm(direction)
                if size <= 1e-8:
                    raise ValueError('Physical-hole tip has no resolved pavement bisector')
                direction /= size
                origin = loop.xy[vertex]
                hits = []
                for target, boundary in enumerate(group.loops):
                    vector = np.roll(boundary.xy, -1, axis=0)-boundary.xy
                    determinant = cross(direction, vector)
                    valid = abs(determinant) > 1e-12*np.linalg.norm(vector, axis=1)
                    delta = boundary.xy-origin
                    distance = np.divide(cross(delta, vector), determinant,
                        out=np.full(len(vector), np.inf), where=valid)
                    fraction = np.divide(cross(delta, direction), determinant,
                        out=np.full(len(vector), np.inf), where=valid)
                    for at in np.flatnonzero((distance > 1e-7) & (fraction >= 0.) & (fraction <= 1.)):
                        hits.append((float(distance[at]), target, int(at), float(fraction[at])))
                if not hits:
                    raise ValueError('Physical-hole tip has no visible pavement-side boundary')
                _, target, at, fraction = min(hits)
                if labels[target][at].kind == 'mouth':
                    from .reflex_cuts import ReflexCurbCuts
                    at = ReflexCurbCuts.profile_endpoint(group.loops[target].xy, labels[target], at, fraction)
                    fraction = 0.
                else:
                    reference = self.endpoint(group, (0., target, at, fraction))
                    at, fraction = reference.edge, reference.fraction
                cuts.append(JunctionCut(None, owner, owner, BoundaryRef(li, int(vertex), 0.),
                    BoundaryRef(target, at, fraction)))
        if len(cuts) > self.maximum_cuts:
            raise ValueError('Physical-hole partition exceeds its cut budget')
        # Two separate connections remove the repeated bridge vertex that a
        # single hole-to-exterior cut would leave in the half-edge face walk.
        provisional = CompoundPartitioner(policy=_FixedCuts(cuts),
            max_cuts=self.maximum_cuts, require_visibility=False).build(group, context)
        seen = set()
        for cut in cuts:
            seen.add(tuple(sorted(((cut.first.loop, cut.first.edge, cut.first.fraction),
                                   (cut.last.loop, cut.last.edge, cut.last.fraction)))))
        for block in provisional.blocks:
            indices = list(block.vertices)
            if len(indices) < 4:
                continue
            for refs in self._visible_diagonals(provisional,indices):
                token = tuple(sorted((ref.loop, ref.edge, ref.fraction) for ref in refs))
                if token not in seen:
                    seen.add(token)
                    cuts.append(JunctionCut(None, owner, owner, *refs))
                    if len(cuts) > self.maximum_cuts:
                        raise ValueError('Physical-hole partition exceeds its cut budget')
        return tuple(cuts)

    def _visible_diagonals(self,provisional,indices):
        policy=VisibilitySplitter(max_candidates=3072,max_cuts=32,minimum_clearance_ratio=.01)
        physical=[i for i,node in enumerate(indices) if node in provisional.boundary_refs]
        try:
            diagonals=policy.split(provisional.xy[indices],endpoints=physical)
        except ValueError as failure:
            if (not self.computational_extensions
                    or str(failure)!='Visibility decomposition needs another compatible cut family'):
                raise
            semantic=set(range(len(indices)))-set(physical)
            if not semantic:raise
            diagonals=policy.split(provisional.xy[indices],endpoints=None,semantic_endpoints=semantic)
        if all(indices[a] in provisional.boundary_refs and indices[b] in provisional.boundary_refs
               for a,b in diagonals):
            return tuple((provisional.boundary_refs[indices[a]],provisional.boundary_refs[indices[b]])
                         for a,b in diagonals)
        result=[];group=provisional.group
        sources={source.region.owner:source for source in group.footprints}
        curb_edges={(li,ei) for li,loop in enumerate(group.loops) for ei,span in enumerate(loop.spans)
                    if rim_label(span,sources).kind=='curb'}
        for a,b in diagonals:
            first,last=(indices[i] for i in (a,b))
            ar,br=(provisional.boundary_refs.get(i) for i in (first,last))
            if ar is not None and br is not None:
                result.append((ar,br));continue
            from .computational_cut_extension import extend
            refs,proof=extend(group,provisional.xy[first],provisional.xy[last],curb_edges=curb_edges,
                first_ref=ar,last_ref=br)
            self.extension_reports.append(proof);result.append(refs)
        return tuple(result)

    def _source_inventory(self, group, context):
        """Choose computational axes without changing physical source claims."""
        claimed = {claim.section_id for claim in group.claims}
        selected, paired = set(), set()
        proofs = getattr(group.physical_layout, 'source_cycle_proofs', {})
        for bundle in context.network.bundles:
            if bundle.id in proofs:
                paired.update(bundle.sections)
                selected.update(visit.section_id for visit in bundle.members[0].path)
        selected.update(claimed-paired)
        return selected, paired

    def _source_cuts(self, group, context, labels):
        """Keep cross-sections aligned with the original consumed corridors."""
        owner = group.members[0]
        cuts = []
        selected, paired = self._source_inventory(group, context)
        self.source_cut_reports=[]
        for sid in sorted(selected):
            path = context[Key('path', sid)]
            ranges=((0.,path.length),) if self.source_ranges is None else self.source_ranges.get(sid,())
            if not ranges or any(last>path.length for first,last in ranges):
                raise ValueError('Computational source normals escaped their finite original owned ranges')
            def intervals(station):
                point, tangent, normal = path.frame(station)
                hits = []
                for li, loop in enumerate(group.loops):
                    along = (loop.xy-point) @ tangent
                    following = np.roll(along, -1)
                    selected = ((along <= 0.) & (0. < following)) | ((following <= 0.) & (0. < along))
                    for edge in np.flatnonzero(selected):
                        fraction = float(-along[edge]/(following[edge]-along[edge]))
                        at = loop.xy[edge]+fraction*(loop.xy[(edge+1) % len(loop.xy)]-loop.xy[edge])
                        hits.append((float((at-point) @ normal), li, int(edge), fraction))
                hits.sort()
                if len(hits) % 2:
                    raise ValueError('Original corridor cut has an unresolved physical crossing')
                result = []
                for first, last in zip(hits[::2], hits[1::2]):
                    if sid not in paired and not first[0] < 0. < last[0]:
                        continue
                    if any(labels[hit[1]][hit[2]].kind != 'curb' for hit in (first, last)):
                        continue
                    result.append((first, last))
                return result
            bands=0
            for low,high in ranges:
                count = max(2, int(math.ceil((high-low)/max(
                    3.*path.section.width/self.cadence_multiplier,4.*self.spacing))))
                bands+=count
                if bands>64:
                    raise ValueError('Physical-hole source corridor exceeds its shared finite range band budget')
                for nominal in np.linspace(low,high,count+1)[1:-1]:
                    for first, last in intervals(float(nominal)):
                        if self.maximum_unpaired_span_widths is not None and sid not in paired:
                            def point(hit):
                                _,li,ei,fraction=hit;xy=group.loops[li].xy
                                return xy[ei]+fraction*(xy[(ei+1)%len(xy)]-xy[ei])
                            length=float(np.linalg.norm(point(last)-point(first)))
                            retain=length<=self.maximum_unpaired_span_widths*path.section.width
                            self.source_cut_reports.append(dict(original_section=sid,
                                owned_interval=(low,high),station=float(nominal),physical_span=length,
                                maximum_unpaired_span_widths=self.maximum_unpaired_span_widths,
                                retained=retain,computational_only=True))
                            if not retain:continue
                        cuts.append(JunctionCut(None, owner, owner,
                            self.endpoint(group, first), self.endpoint(group, last)))
                        if len(cuts) > self.maximum_cuts:
                            raise ValueError('Physical-hole partition exceeds its cut budget')
        return cuts
