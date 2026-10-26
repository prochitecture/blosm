"""Resolved full/half crowns over the shared authored cross-section contract."""
from dataclasses import dataclass
import numpy as np
from .model import Key
from .updates import FeatureProvider, Task


@dataclass(frozen=True)
class Crown:
    rise: float
    datum: float
    half: int = 0  # 0 full crown, +1 rising right, -1 rising left.

    @property
    def peak(self):
        return float(self.half)

    def value(self, q):
        distance = (np.asarray(q) - self.peak) / (1. + abs(self.half))
        return self.rise * (1. - distance * distance)

    def derivative(self, q):
        return -2. * self.rise * (np.asarray(q) - self.peak) / (1. + abs(self.half)) ** 2

    def oriented(self, end):
        return self if end == 0 else Crown(self.rise, self.datum, -self.half)

    def coordinate(self, side, v):
        return side + (self.peak - side) * np.asarray(v)

    def intervals(self, side, base):
        return int(abs(self.peak - side)) * base


def plateau_extension(xy, boundary, datum, height):
    """Bounded physical-distance extension of a prescribed separator rim.

    Invert the quintic smoothstep on the rim, then extend its values with the
    minimal sufficient Euclidean Lipschitz rate. The upper distance envelope
    reaches the exact plateau away from lowered end caps. Applying smoothstep
    again removes the derivative jump at the plateau/zero cutoff. The field
    can have distance-envelope seams; finite-quad slopes are measured separately.

    Only boundary points below the plateau can influence this upper extension.
    Scratch arrays are bounded by 128 rows, independent of surface resolution.
    """
    if height <= 1e-15:
        return np.full(len(xy), datum), dict(rate=0., continuous_slope_bound=0., active_samples=0)
    target = np.clip((boundary[:, 2] - datum) / height, 0., 1.)
    lo, hi = np.zeros(len(boundary)), np.ones(len(boundary))
    smoothstep = lambda t: t ** 3 * (10. + t * (-15. + 6. * t))
    for _ in range(48):
        middle = .5 * (lo + hi)
        lower = smoothstep(middle) < target
        lo, hi = np.where(lower, middle, lo), np.where(lower, hi, middle)
    values = .5 * (lo + hi)
    # Subtracting a translated datum loses a few ulps of the absolute height.
    # Inverse smoothstep is ill-conditioned at its endpoints: treating that
    # rounding as a real height difference can falsely reject coincident caps.
    scale = max(1., abs(datum), float(np.max(abs(boundary[:, 2]))))
    endpoint_tolerance = min(1e-10, max(1e-14, 8. * np.spacing(scale) / height))
    values[target <= endpoint_tolerance], values[target >= 1. - endpoint_tolerance] = 0., 1.
    active = np.flatnonzero(values < 1.)
    rate = 0.
    for start in range(0, len(active), 128):
        ids = active[start:start + 128]
        distance = np.linalg.norm(boundary[ids, None, :2] - boundary[None, :, :2], axis=-1)
        difference = abs(values[ids, None] - values[None, :])
        if np.any((distance <= 1e-10) & (difference > 1e-9)):
            raise ValueError('Separator boundary has conflicting heights at one physical point')
        ratios = np.divide(difference, distance, out=np.zeros_like(distance), where=distance > 1e-10)
        rate = max(rate, float(ratios.max()))
    rate *= 1.000001  # Numerical margin retains every prescribed rim value.
    result = np.full(len(xy), datum + height)
    if len(active) and rate > 0.:
        low, high = boundary[active, :2].min(axis=0) - 1. / rate, boundary[active, :2].max(axis=0) + 1. / rate
        affected = np.flatnonzero(np.all((xy >= low) & (xy <= high), axis=1))
        for start in range(0, len(affected), 128):
            ids = affected[start:start + 128]
            distance = np.linalg.norm(xy[ids, None, :] - boundary[None, active, :2], axis=-1)
            field = np.minimum(1., np.min(values[active] + rate * distance, axis=1))
            result[ids] = datum + height * smoothstep(field)
    return result, dict(rate=rate, continuous_slope_bound=height * 1.875 * rate, active_samples=len(active))


def memberships(network):
    result = {}
    for bundle in network.bundles:
        for role, member in enumerate(bundle.members):
            for traversal in member.path:
                if traversal.section_id in result:
                    raise ValueError('A complete Section has overlapping Bundle memberships')
                result[traversal.section_id] = (bundle, role, traversal.forward)
    return result


class ProfileProvider(FeatureProvider):
    source_kinds = frozenset(('section',))

    def tasks(self, network):
        ownership = memberships(network)
        for section in network.sections:
            member = ownership.get(section.id)
            deps = (section.key,) if member is None else tuple(sorted({section.key, member[0].key,
                *(Key('section', t.section_id) for m in member[0].members for t in m.path)}))
            yield Task(Key('profile', section.id), deps,
                lambda context, source=section, role=member: self.build(source, role, context.network.section_index),
                signature=lambda context, source=section, role=member: self.build(source, role, context.network.section_index))

    @staticmethod
    def build(section, membership, sections):
        half = 0
        rise = section.rise
        if membership is not None:
            bundle, role, forward = membership
            if len(bundle.members) != 2 or not all(m.complete for m in bundle.members):
                raise ValueError('Shared profile requires two complete lateral Bundle members')
            half = (1 if role == 0 else -1) * (1 if forward else -1)
            rise = bundle.profile_rise if bundle.profile_rise is not None else max(
                sections[t.section_id].rise for m in bundle.members for t in m.path)
        return Crown(rise, section.datum, half)
