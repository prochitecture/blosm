"""Complete physical entrances for a compound with consumed profile seams.

Logical junctions, Connections and Sections remain original source records.
Only an invalid compound exterior is replaced; valid footprint unions keep
their existing boundaries and charts. Corner radius ownership stays attached
to the original outgoing incidence, including strict/style override semantics.
"""
from dataclasses import dataclass, replace
import math
import numpy as np
from .compounds import rim_label, sample_footprint
from .compound_coverage import CompoundSourceCoverage
from .mesh import cross
from .model import End, Key
from .regions import RegionUnion, _classify


@dataclass(frozen=True)
class ExteriorSpecification:
    """Derived contact request, never a synthetic source Intersection."""
    id: int
    radii: dict

    def corner_radius(self, end):
        return self.radii[end]


class CompoundPortPlanner:
    maximum_support_samples = 20000
    support_chunk = 32

    def __init__(self, context, junctions, layouts):
        self.context = context
        self.sites = {site.key: site for site in junctions}
        self.layouts = layouts
        self._coverage = {}

    def coverage(self, group):
        if group.members not in self._coverage:
            self._coverage[group.members] = CompoundSourceCoverage(self.context, group, self.layouts)
        return self._coverage[group.members]

    def validate_deferred_profiles(self, groups):
        """A provisional height lead may disappear only inside its compound."""
        by_site = {key: group for group in groups for key in group.members}
        checked = set()
        for key, layout in self.layouts.items():
            for approach in layout.approaches:
                pending = getattr(approach.path, 'deferred_profile_reserve', None)
                if pending is None:
                    continue
                group = by_site.get(key)
                if group is None or pending.terminal_site not in group.members:
                    raise ValueError('Deferred profile reserve requires actual opposite-junction compound ownership')
                if not set(pending.section_ids) <= self.coverage(group).fully_owned:
                    raise ValueError('Deferred profile reserve has unowned original source intervals')
                if len(group.loops) != 1 or group.loops[0].hole:
                    raise ValueError('Deferred profile reserve needs another compound boundary policy')
                token = group.members, tuple(sorted(pending.section_ids))
                if token not in checked:
                    self.check_support(pending.section_ids, group.loops[0])
                    checked.add(token)

    @staticmethod
    def incomplete_ports(group):
        """Inspect the same physical boundary that sampling will consume."""
        sources = {p.region.owner: p for p in group.footprints}
        for loop in group.loops:
            labels = [rim_label(span, sources) for span in loop.spans]
            for i, span in enumerate(loop.spans):
                previous = loop.spans[i - 1]
                if previous.owner != span.owner and (
                        labels[i - 1].kind == 'mouth' or labels[i].kind == 'mouth'):
                    return True
            tokens = [(span.owner, label.index) if label.kind == 'mouth' else None
                      for span, label in zip(loop.spans, labels)]
            for i, token in enumerate(tokens):
                if token is None or token == tokens[i - 1]:
                    continue
                last = i
                for _ in range(len(tokens)):
                    at = (last + 1) % len(tokens)
                    if tokens[at] != token:
                        break
                    last = at
                else:
                    return True
                if abs(labels[i].first - 1.) > 1e-8 or abs(labels[last].last) > 1e-8:
                    return True
        return False

    def build(self, group):
        if not self.incomplete_ports(group):
            return group
        layouts = {key: self.layouts[key] for key in group.members}
        claimed = {}
        for key, layout in layouts.items():
            for site in layout.absorbed_transitions:
                claimed.setdefault(site, set()).add(key)
        consumed = self.coverage(group).internal_sections()
        seams = {site.key: site for site in self.context.network.connections}
        shared = [seams[key] for key, owners in claimed.items()
                  if key in seams and owners
                  and all(e.section_id in consumed for e in seams[key].ends)]
        if not shared:
            return group
        radii = {end: self.sites[key].corner_radius(end) for key in group.members
                 for end in self.sites[key].ends if end.section_id not in consumed}
        approaches = [a for layout in layouts.values() for a in layout.approaches
                      if End(a.source.id, a.end) in radii]
        identities = [End(a.source.id, a.end) for a in approaches]
        if len(identities) < 3 or len(set(identities)) != len(radii) or len(identities) != len(radii):
            raise ValueError('Compound exterior needs distinct complete source incidences')
        # Bundle half-profiles and separator holes need a compound boundary
        # policy retaining their inner plateau boundaries, not a full-crown
        # exterior substitution. Leave those established branches in place.
        if any(a.profile.half != 0 for a in approaches) or any(loop.hole for loop in group.loops):
            return group
        if max(a.profile.datum for a in approaches) - min(a.profile.datum for a in approaches) > 1e-8:
            raise ValueError('Compound exterior requires the original common minimum height')
        if any(abs(a.profile.datum - approaches[0].profile.datum) > 1e-7
               for layout in layouts.values() for a in layout.approaches):
            raise ValueError('Consumed compound approaches require the original common minimum height')
        from .intersections import IntersectionProvider
        specification = ExteriorSpecification(group.members[0].id, radii)
        physical = IntersectionProvider.curved_boundary(self.context, specification, list(approaches))
        mouths = {}
        for layout in layouts.values():
            for identity, mouth in layout.mouths.items():
                previous = mouths.get(identity)
                if previous is None or (mouth.station > previous.station if identity[1] == 0
                                        else mouth.station < previous.station):
                    mouths[identity] = mouth
        mouths.update(physical.mouths)
        physical = replace(physical, mouths=mouths,
            absorbed_transitions=tuple(sorted(claimed)))
        footprint = sample_footprint(group.members[0], physical)
        loops = RegionUnion().build((footprint.region,))
        if len(loops) != 1 or loops[0].hole:
            raise ValueError('Complete compound entrances need one simple physical exterior')
        self.check_support(consumed, loops[0])
        return replace(group, footprints=(footprint,), loops=loops, physical_layout=physical)

    def check_support(self, consumed, loop):
        """Keep complete consumed centerlines strictly inside the new exterior.

        Bound the chord approximation of Smooth paths, then test every segment
        against every boundary edge in small blocks. Endpoints alone would
        miss a segment leaving/re-entering a concave compound. The clearance
        margin includes both source and boundary approximation tolerances.
        """
        settings = self.context.network.settings
        tolerance = min(.001, settings.curve_tolerance * .1)
        clearance = tolerance + settings.curve_tolerance + 1e-7
        x, y = loop.xy, np.roll(loop.xy, -1, axis=0)
        edge = y - x
        edge2 = np.sum(edge * edge, axis=1)
        total = 0
        for sid in sorted(consumed):
            path = self.context[Key('path', sid)]
            spacing = min(.5, settings.spacing)
            if not path.exact and path.curvature_bound > 1e-12:
                reach = abs(path.section.offset) + .5 * path.section.width
                curvature = path.curvature_bound / (1. - reach * path.curvature_bound)
                speed = 1. + reach * path.curvature_bound
                spacing = min(spacing, math.sqrt(8. * tolerance / curvature) / speed,
                              math.pi / (2. * curvature * speed))
            if not math.isfinite(spacing) or spacing <= 0. or path.length / spacing + total > self.maximum_support_samples - 4:
                raise ValueError('Consumed compound centerline support exceeds its sample budget')
            stations = path.samples(0., path.length, spacing, tolerance=tolerance)
            total += len(stations)
            if total > self.maximum_support_samples:
                raise ValueError('Consumed compound centerline support exceeds its sample budget')
            points = path.frame(stations)[0]
            inside, boundary = _classify(points, x, clearance)
            if not np.all(inside & (boundary < 0)):
                raise ValueError('Complete compound exterior does not retain its consumed centerline')
            for start in range(0, len(points) - 1, self.support_chunk):
                a = points[start:start + self.support_chunk]
                b = points[start + 1:start + self.support_chunk + 1]
                a = a[:len(b)]
                d = b - a
                offset = x[None] - a[:, None]
                determinant = cross(d[:, None], edge[None])
                valid = abs(determinant) > 1e-12 * np.linalg.norm(d, axis=1)[:, None] * np.sqrt(edge2)[None]
                t = np.divide(cross(offset, edge[None]), determinant, out=np.zeros_like(determinant), where=valid)
                u = np.divide(cross(offset, d[:, None]), determinant, out=np.zeros_like(determinant), where=valid)
                if np.any(valid & (t >= 0.) & (t <= 1.) & (u >= 0.) & (u <= 1.)):
                    raise ValueError('A consumed centerline exits its compound exterior')
                # Without crossings, the minimum segment distance occurs at
                # an endpoint projected onto the other segment (also parallel).
                distances = []
                for p in (a, b):
                    v = p[:, None] - x[None]
                    f = np.clip(np.sum(v * edge[None], axis=2) / edge2[None], 0., 1.)
                    distances.append(np.linalg.norm(v - f[..., None] * edge[None], axis=2))
                d2 = np.sum(d * d, axis=1)
                for p in (x, y):
                    v = p[None] - a[:, None]
                    f = np.clip(np.sum(v * d[:, None], axis=2) / d2[:, None], 0., 1.)
                    distances.append(np.linalg.norm(v - f[..., None] * d[:, None], axis=2))
                if min(float(distance.min()) for distance in distances) <= clearance:
                    raise ValueError('Consumed centerline lacks certified clearance from its compound exterior')
