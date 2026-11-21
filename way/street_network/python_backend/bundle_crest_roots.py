"""Close an exposed outer root using two consecutive median-side Bundle banks."""
from dataclasses import dataclass
from hashlib import sha256
import math
from types import SimpleNamespace

import numpy as np

from .compound_root_envelopes import RootEnvelope
from .compounds import rim_label
from .generated_junction_curbs import GeneratedJunctionArea
from .model import End, Key
from .profiles import memberships
from .regions import RegionUnion
from .routes import route_key
from .source_root_ports import NativeRolePort, refine_native_contact
from .source_root_witnesses import _group_token


def _end(arm):
    return End(arm.source.id, arm.end)


def crest_bundle(network, corner):
    """Return the actual common Bundle, member and flat crest elevation."""
    first, second = corner.first, corner.second
    if (first.profile.half, second.profile.half) != (-1, 1):
        return None
    z0 = first.profile.datum + float(first.profile.value(-1.))
    z1 = second.profile.datum + float(second.profile.value(1.))
    if (z0 != z1 or not math.isfinite(z0) or z0 <= max(first.profile.datum, second.profile.datum)
            or first.profile.derivative(-1.) != 0. or second.profile.derivative(1.) != 0.):
        return None
    roles = memberships(network)
    a, b = roles.get(first.source.id), roles.get(second.source.id)
    if a is None or b is None or a[0] is not b[0] or a[1] != b[1]:
        return None
    bundle, role, _ = a
    walk = bundle.members[role].path
    ends = {_end(first), _end(second)}
    if any({End(left.section_id, int(left.forward)), End(right.section_id, int(not right.forward))} == ends
           for left, right in zip(walk, walk[1:])):
        return bundle, role, z0
    return None


@dataclass(frozen=True, eq=False)
class CrestSource:
    site: object
    group: object
    donor: object
    owner: object
    bundle: object
    role: int
    height: float
    bounds: object
    geometry: tuple
    objects: tuple

    @classmethod
    def build(cls, context, group, site, donor, owner, approaches, limits, bounds):
        from .cycle_dimensional_retry import certify_layout
        from .junction_bank_transitions import original_corner
        from .ordinary_singleton_sources import OriginalSingletonRootBounds
        from .ownership import footprint_key
        from .source_strip_containment import SourceStripContainmentError

        if (group is None or group.members != (site.key,) or group.claims or len(group.loops) != 1
                or group.loops[0].hole or type(bounds) is not OriginalSingletonRootBounds
                or not isinstance(bounds.failure, SourceStripContainmentError) or bounds.site != site.key
                or bounds.layout is not context.products.get(footprint_key(site.key))
                or limits != tuple(bounds.limits()[_end(a)] for a in (donor.first, donor.second))
                or donor.radius_mode != 'STYLE' or context.network.records().get(site.key) != site):
            raise ValueError('Bundle crest needs its original bounded singleton and style corner')
        layout = group.physical_layout
        arms = {_end(a): a for a in layout.approaches}
        inventory = {s.region.owner: s for s in group.area_footprints}
        resolver = SimpleNamespace(context=context, group=group, sources=inventory)
        if (set(arms) != set(site.ends) or len(arms) != len(approaches)
                or any(arms.get(_end(a)) is not a for a in approaches)
                or not original_corner(resolver, site, donor)):
            raise ValueError('Bundle crest needs the original cyclic incidence and corner')
        donors = []
        ports = []
        keys = set()
        for source in inventory.values():
            if len(source.layout.approaches) == 1:
                arm, = source.layout.approaches
                if arms.get(_end(arm)) is not arm or context[route_key(_end(arm))] is not arm.path:
                    raise ValueError('Bundle crest port changed its original approach')
                NativeRolePort.build(context, source, arm)
                ports.append(_end(arm))
                keys.add(route_key(_end(arm)))
                for piece in arm.path.pieces:
                    if context[Key('path', piece.section.id)] is not piece.path:
                        raise ValueError('Bundle crest lost an original Section path')
                    keys.update(Key(kind, piece.section.id) for kind in ('path', 'profile'))
            else:
                if (getattr(source.layout, 'original_junction', None) != site.key or not source.layout.corners
                        or any(not any(c is old for old in layout.corners) for c in source.layout.corners)):
                    raise ValueError('Bundle crest cannot borrow another corner')
                if any(c is donor for c in source.layout.corners):
                    donors.append(source.region.owner)
        if sorted(ports) != sorted(site.ends) or donors != [owner]:
            raise ValueError('Bundle crest needs the complete source inventory and exactly one donor')
        membership = crest_bundle(context.network, donor)
        if membership is None:
            raise ValueError('Bundle crest donor needs consecutive banks at equal crest height')
        bundle, role, height = membership
        certify_layout(context, group, ())
        for areas in (inventory.values(), [s for key, s in inventory.items() if key != owner]):
            RegionUnion(max_segments=20000).build(s.region for s in areas)
        return cls(site, group, donor, owner, bundle, role, height, bounds, _group_token(group),
                   tuple((k, context[k]) for k in sorted(keys)))

    def validate(self, context):
        from .ownership import footprint_key
        role = crest_bundle(context.network, self.donor)
        if (role is None or role[0] is not self.bundle or role[1:] != (self.role, self.height)
                or self.bounds.layout is not context.products.get(footprint_key(self.site.key))
                or self.site != context.network.records().get(self.site.key)
                or _group_token(self.group) != self.geometry or any(context[k] is not v for k, v in self.objects)):
            raise ValueError('Bundle crest lost its original profile, source geometry or mouth bounds')

    def refine(self, context, source, contact, maximum_samples=4096):
        arm, = source.layout.approaches
        if not any(s.region.owner == source.region.owner and len(s.layout.approaches) == 1
                   and s.layout.approaches[0] is arm for s in self.group.area_footprints):
            raise ValueError('Bundle crest contact belongs to a different source port')
        port = NativeRolePort.build(context, source, arm)
        refined, _ = refine_native_contact(context, port, contact, maximum_samples=maximum_samples)
        NativeRolePort.build(context, refined, arm)
        if not set(map(tuple, source.region.xy)).issubset(map(tuple, refined.region.xy)):
            raise ValueError('Bundle crest refinement moved an existing source point')
        return refined


class CrestSearch:
    """Cache candidate contact samples only for the duration of one search."""
    def __init__(self, source, context):
        self.source = source
        self.context = context
        self.donor = source.donor
        self.samples = {}

    def regions(self, candidate, requested):
        key = tuple(candidate.stations)
        if key not in self.samples:
            contacts = dict(zip((_end(candidate.first), _end(candidate.second)), key))
            fixed = [s for s in self.source.group.area_footprints if s.region.owner != self.source.owner]
            total = sum(len(s.region.xy) for s in fixed)
            regions = {}
            for source in fixed:
                arms = source.layout.approaches
                if len(arms) == 1 and _end(arms[0]) in contacts:
                    size = len(source.region.xy)
                    source = self.source.refine(self.context, source, contacts[_end(arms[0])],
                                                min(4096, (20000-total)//2))
                    total += len(source.region.xy) - size
                regions[source.region.owner] = source.region
            RegionUnion(max_segments=20000).build(regions.values())
            self.samples[key] = regions
        return tuple(self.samples[key][region.owner] for region in requested)


def _controls(curve):
    digest = sha256()
    for piece in curve.pieces:
        digest.update(piece.curve.continuation.controls.tobytes())
    return digest.hexdigest()


@dataclass(frozen=True)
class CrestEnvelope(RootEnvelope, GeneratedJunctionArea):
    reference: CrestSource
    controls: str
    refined: bool

    def validate(self, context, group):
        self.reference.validate(context)
        if (self.source.layout.generated_junction_area is not self
                or sum(b is self.curve for b in self.source.layout.boundaries) != 1
                or self.corner.stations != self.stations or _controls(self.curve) != self.controls
                or not any(c is self.corner for c in group.physical_layout.corners)
                or getattr(self.corner, 'root_support_envelope', None) is not self):
            raise ValueError('Bundle crest lost its original curve or source ownership')
        for a, s, q, t in zip((self.corner.first, self.corner.second), self.stations, (-1, 1), (0., 1.)):
            if (a.profile.datum + float(a.profile.value(q)) != self.reference.height
                    or not np.array_equal(a.physical_gradient(s, q), np.zeros(2))
                    or np.linalg.norm(self.curve(t)-a.point(s, a.offset+.5*a.width*q)) > 1e-8):
                raise ValueError('Bundle crest changed an original contact position or slope')


def prepare(owner, donor, settings, *, context, source_group, original_bounds, **geometry):
    from .updates import drain_steps
    return drain_steps(prepare_steps(owner,donor,settings,context=context,
        source_group=source_group,original_bounds=original_bounds,**geometry))


def prepare_steps(owner, donor, settings, *, context, source_group, original_bounds, **geometry):
    from .cyclic_root_support import _prepare_steps
    ref = CrestSource.build(context, source_group, geometry['site'], donor, owner, geometry['original_approaches'],
                            geometry['maximum_stations'], original_bounds)
    search = CrestSearch(ref, context)
    record = yield from _prepare_steps(owner, donor, settings, crest=search, **geometry)
    result = CrestEnvelope(record.source, record.curve, record.stations,
                           dict(record.report, bundle=ref.bundle.id, crest_height=ref.height),
                           record.corner, ref, _controls(record.curve), record.stations in search.samples)
    result.source.layout.generated_junction_area = result
    return result


def signature(context, group):
    entries = []
    for source in group.area_footprints:
        record = getattr(source.layout, 'generated_junction_area', None)
        if not isinstance(record, CrestEnvelope):
            continue
        record.validate(context, group)
        inventory = {s.region.owner:s for s in group.footprints}
        if record.refined:
            for arm, station in zip((record.corner.first, record.corner.second), record.stations):
                port, = (s for s in record.reference.group.area_footprints
                         if len(s.layout.approaches)==1 and s.layout.approaches[0] is arm)
                actual = inventory.get(port.region.owner)
                if actual is None:
                    raise ValueError('Final union lost a Bundle crest port')
                # Other awarded contacts can change the port's background grid.
                # Replaying this particular refinement must add zero samples.
                checked = record.reference.refine(context, actual, station)
                if (not np.array_equal(checked.region.xy,actual.region.xy)
                        or not np.array_equal(checked.layout.boundaries[0].source_coordinates.stations,
                                              actual.layout.boundaries[0].source_coordinates.stations)):
                    raise ValueError('Final union is missing accepted Bundle crest contact samples')
        labels = (rim_label(span, inventory) for loop in group.loops if not loop.hole
                  for span in loop.spans if span.owner == source.region.owner)
        if abs(sum(abs(label.last-label.first) for label in labels if label.kind=='curb')-1.) > 1e-8:
            raise ValueError('Final union clips the Bundle crest before its tangent endpoints')
        ref = record.reference
        entries.append((ref.site.key.text(), ref.bundle.id, ref.role, ref.height, record.stations,
                        record.controls, tuple((k.text(), context.digest(k)) for k, _ in ref.objects)))
    return tuple(entries) if entries else None
