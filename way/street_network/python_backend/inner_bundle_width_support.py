"""Finite plateau curb support at an exposed original Bundle width change.

Only a real hole's exposed inner root can request this last-family support.
The two approaches are consecutive original member traversals at the same
owned junction. Their profiles and the existing outer circular corners stay
unchanged; the enclosing source union must qualify the resulting boundary.
"""
from dataclasses import dataclass, replace
import math
import numpy as np

from .compounds import rim_label
from .intersections import Corner
from .model import End, Key


@dataclass(frozen=True)
class InnerBundleWidthSupport:
    source: object
    curve: object
    stations: tuple
    report: dict
    corner: object


def required(planner, loops, sources, retried, consumed, members, settings):
    inventory = {source.region.owner: source for source in sources}
    incidents = {end: site for site in planner.context.network.intersections for end in site.ends}
    exposed = {}
    for loop in loops:
        if not loop.hole:
            continue
        for span in loop.spans:
            if rim_label(span, inventory).kind != 'internal':
                continue
            source = inventory[span.owner]
            end = getattr(source.layout, 'original_root_caps', {}).get(span.segment)
            site = incidents.get(end)
            if site is None or site.key not in members or end.section_id in consumed:
                continue
            approach = next((a for a in planner.layouts[site.key].approaches
                if End(a.source.id, a.end) == end), None)
            if approach is None or approach.profile.half not in (-1, 1):
                continue
            xy = source.region.xy
            first, delta = xy[span.segment], xy[(span.segment+1) % len(xy)]-xy[span.segment]
            points = first+np.asarray((span.first, span.last))[:, None]*delta
            length = float(np.linalg.norm(points[1]-points[0]))
            if length == 0.:
                continue
            low = approach.point(0., approach.offset-.5*approach.width)
            high = approach.point(0., approach.offset+.5*approach.width)
            axis = high-low
            u = (points-low)@axis/(axis@axis)
            if (np.max(np.linalg.norm(points-(low+u[:, None]*axis), axis=1)) > 1e-7
                    or u.min() < -1e-8 or u.max() > 1.+1e-8):
                raise ValueError('Inner Bundle support lost its actual original root cross-section')
            if (approach.profile.half == -1 and u.max() <= .5+1e-8
                    or approach.profile.half == 1 and u.min() >= .5-1e-8):
                exposed[end] = dict(original_root_parameters=u.tolist(),
                    exposed_length=length)
    if not exposed:
        return None
    for bundle in planner.context.network.bundles:
        if len(bundle.members) != 2 or not all(member.complete for member in bundle.members):
            continue
        for role, member in enumerate(bundle.members):
            for before, after in zip(member.path, member.path[1:]):
                ends = (End(before.section_id, int(before.forward)),
                        End(after.section_id, int(not after.forward)))
                site = incidents.get(ends[0])
                token = bundle.id, role, site.key if site else None
                if (token in retried or not set(ends).intersection(exposed) or site is None
                        or not 3 <= len(site.ends) <= 8 or incidents.get(ends[1]) is not site
                        or site.key not in members or any(end.section_id in consumed for end in ends)):
                    continue
                by_end = {End(a.source.id, a.end): a for a in planner.layouts[site.key].approaches}
                if any(end not in by_end for end in ends):
                    continue
                pair = tuple(sorted((by_end[end] for end in ends), key=lambda a: a.profile.half))
                a, b = pair
                if (a.profile.half != -1 or b.profile.half != 1
                        or any(piece.section.id in consumed for approach in pair for piece in approach.path.pieces)
                        or abs(a.width-b.width) <= 1e-6):
                    continue
                requests = [site.corner_radius(End(approach.source.id, approach.end)) for approach in pair]
                if any(mode != 'STYLE' for radius, mode in requests):
                    continue
                gap = math.atan2(float(np.linalg.det(np.stack((a.direction, b.direction)))),
                                 float(a.direction@b.direction)) % math.tau
                if abs(gap-math.pi) > math.radians(15.):
                    continue
                radius = max(radius for radius, mode in requests)
                corner = Corner(a, b, (0., 0.), None, radius, 0., 0., gap,
                    continuation=True, requested_radius=radius, radius_mode='STYLE',
                    fit_policy='EXPOSED_ORIGINAL_INNER_BUNDLE_WIDTH_ROOT')
                from .compound_dimensional_corners import prepare
                record = prepare(Key('inner_bundle_support', 0), corner, settings,
                                 original_inner_half_profiles=True)
                record.report.update(policy='EXPOSED_ORIGINAL_INNER_BUNDLE_WIDTH_ROOT',
                    bundle=bundle.id, member=role, original_junction=site.key.text(),
                    original_member_pair=[dict(section=e.section_id, end=e.end) for e in ends],
                    actual_exposed_roots=[dict(section=e.section_id, end=e.end, **exposed[e])
                                          for e in ends if e in exposed],
                    physical_hole_only=True, outer_circular_corners_unchanged=True)
                corner=replace(corner,stations=record.stations)
                record.curve.corner=corner
                return token, InnerBundleWidthSupport(record.source, record.curve,
                    record.stations, record.report, corner)
    return None
