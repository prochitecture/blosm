"""Failed-only STYLE noses of genuine, narrow original-source median holes.

The ordinary radius floor is a style heuristic. Two complete original roads
can enclose a real separator narrower than that floor. This bounded fallback
uses their original finite banks, chooses a radius before any surface exists,
and still requires the complete source union to expose the whole new circle.
"""
from dataclasses import dataclass, replace
import math
import numpy as np
from .model import End, Key
from .regions import Region, RegionUnion
from .mesh import cross


@dataclass(frozen=True)
class ThinMedianGore:
    source_hole: np.ndarray
    report: dict


def radius_candidates(area, perimeter, old_minimum):
    """Keep the largest bounded nose below the characteristic mean width.

    2A/P is a characteristic width of a long separator. Requiring a nose's
    diameter to be no larger than that width avoids consuming most of a
    tapered lens merely because its widest cross-section admits a circle.
    This is an explicit style choice, not a curvature or containment proof.
    """
    if not all(np.isfinite(x) and x > 0. for x in (area, perimeter, old_minimum)):
        raise ValueError('A thin median needs positive finite area, perimeter and old radius')
    cap = area / perimeter
    return tuple(r for r in (.25, .125) if r <= cap and r < old_minimum)


def _strip_frame(points, direction):
    values = sorted(set(map(tuple, points)))
    def half(rows):
        result = []
        for q in rows:
            while len(result) >= 2 and cross(np.subtract(result[-1], result[-2]),
                                              np.subtract(q, result[-1])) <= 0.:
                result.pop()
            result.append(q)
        return result
    hull = np.array(half(values)[:-1] + half(values[::-1])[:-1])
    edges = np.roll(hull, -1, axis=0) - hull
    normals = np.c_[-edges[:, 1], edges[:, 0]] / np.linalg.norm(edges, axis=1)[:, None]
    if len(points)*len(normals)>8_000_000:
        raise ValueError('Thin median band search exceeds its bounded projection budget')
    widths = np.concatenate([np.ptp(points @ normals[i:i+128].T, axis=0)
                             for i in range(0, len(normals), 128)])
    normal = normals[int(np.argmin(widths))]
    axis = np.array((-normal[1], normal[0]))
    if axis @ direction < 0.: normal, axis = -normal, -axis
    return normal, axis


def prepare(planner, pair, consumed, members, *, exhausted):
    """Prepare one typed, exposed pair after all old contact families failed.

    Admission deliberately needs two full original Sections joining the SAME
    two authored physical junctions. Scope cuts, free ends, derived approaches
    and a hole made only by unrelated junction operands cannot qualify.
    """
    from .corner_radii import CornerRadiusPolicy
    from .intersections import _corner
    from .compound_contact_routes import regular_contact
    from .corner_patch import corner_patch
    if pair not in exhausted:
        raise ValueError('Thin median contacts require exhausted ordinary source contacts')
    first, last = pair
    ids = {first.section_id, last.section_id}
    if len(ids) != 2 or not ids <= set(consumed):
        raise ValueError('A thin median needs two complete consumed original roads')
    sites = planner.sites
    by_end = {e: key for key, site in sites.items() for e in site.ends}
    site = by_end.get(first)
    far = by_end.get(End(first.section_id, 1-first.end))
    if (site is None or far is None or site == far or site not in members or far not in members
            or by_end.get(last) != site or by_end.get(End(last.section_id, 1-last.end)) != far
            or len(sites[site].ends) < 3 or len(sites[far].ends) < 3):
        raise ValueError('A thin median needs both authored physical source junctions')
    radius, mode = sites[site].corner_radius(first)
    if mode != 'STYLE':
        raise ValueError('Only STYLE source radii admit a thin median nose')
    ordered = sorted(planner.layouts[site].approaches,
        key=lambda a: math.atan2(a.direction[1], a.direction[0]) % math.tau)
    endpoints = [End(a.source.id, a.end) for a in ordered]
    if first not in endpoints or endpoints[(endpoints.index(first)+1) % len(endpoints)] != last:
        raise ValueError('Thin median banks must be actual consecutive original incidences')
    a, b = (ordered[endpoints.index(e)] for e in pair)
    settings = planner.context.network.settings
    tolerance = min(.001, .05 * settings.curve_tolerance)
    paths = []; samples = []; regions = []; counts = {}; inner_sides = {}
    for approach, side in ((a, -1), (b, 1)):
        source = approach.source
        path = planner.context[Key('path', source.id)]
        pieces = approach.path.pieces
        if (source.selection_ends or source.scope_ends or approach.profile.half
                or source is not path.section or len(pieces) != 1 or pieces[0].path is not path
                or pieces[0].forward != (approach.end == 0) or path.exact
                or approach.path.terminal_free or getattr(approach.path, 'absorbed_sites', ())
                or hasattr(approach.path, 'mouth_piece') or hasattr(approach, 'original_approaches')
                or getattr(approach.path, 'deferred_profile_reserve', 0.)
                or abs(approach.offset) > 1e-8 or abs(source.offset) > 1e-8
                or abs(approach.width-source.width) > 1e-8
                or abs(approach.profile.rise-source.rise) > 1e-10
                or abs(approach.profile.datum-source.datum) > 1e-10
                or abs(float(approach.profile.value(side))) > 1e-10):
            raise ValueError('Thin median banks require unchanged full original smooth profiles')
        at = path.samples(0., path.length, min(.5, settings.spacing), tolerance=tolerance)
        counts[source.id] = len(at)
        if sum(counts.values()) * 2 > 20000:
            raise ValueError('Thin median original strips exceed their finite sampling budget')
        xy = np.vstack((path.point(at, .5*source.width), path.point(at, -.5*source.width)[::-1]))
        regions.append(Region(source.key, xy))
        paths.append(path); samples.append(at)
        inner_sides[source.key] = side * (1 if approach.end == 0 else -1)
    if abs(a.width-b.width) > 1e-8:
        raise ValueError('Thin median admission currently needs equal original carriageway widths')
    loops = RegionUnion(max_segments=20000).build(regions)
    holes = [loop for loop in loops if loop.hole]
    if len(holes) != 1 or sum(not loop.hole for loop in loops) != 1:
        raise ValueError('Two complete original strips must enclose one genuine median hole')
    hole = holes[0]
    for span in hole.spans:
        count = counts[span.owner.id]
        valid = (0 <= span.segment < count-1 if inner_sides[span.owner] == 1
                 else count <= span.segment < 2*count-1)
        if not valid:
            raise ValueError('A genuine median must be bounded only by the two original facing banks')
    if {span.owner.id for span in hole.spans} != ids:
        raise ValueError('A genuine median needs both original facing banks')
    normal, axis = _strip_frame(hole.xy, a.direction)
    proofs = []; bounds = []
    for approach, path, at in zip((a, b), paths, samples):
        c1, c2, c3 = (path.coefficients[:, i] for i in (1, 2, 3))
        derivative = np.stack((c1, c1+c2, c1+2*c2+3*c3), axis=1).reshape(-1, 2)
        derivative *= 1 if approach.end == 0 else -1
        axial = derivative @ axis
        scale = 1. - .5*path.section.width*path.curvature_bound
        if axial.min() <= 0. or scale <= 0.:
            raise ValueError('Complete original median banks need positive continuous axial tangents')
        cosine = 1. / np.sqrt(1. + float(np.max(abs(derivative @ normal) / axial))**2)
        side = inner_sides[path.section.key]
        projection = path.point(at, side*.5*path.section.width) @ normal
        bounds.append((float(projection.min()-tolerance), float(projection.max()+tolerance)))
        proofs.append(dict(section=path.section.id, original_side=side,
            minimum_axial_derivative=float(axial.min()), normal_projection_lower=cosine,
            offset_tangent_scale_lower=scale, source_chord_error=tolerance, intervals=len(at)-1))
    width = bounds[0][1] - bounds[1][0]
    old = tuple(CornerRadiusPolicy.contact_candidates(a, b, radius, mode))
    required = old[-1] * sum(row['normal_projection_lower'] for row in proofs)
    if width <= 0. or required <= width:
        raise ValueError('The certified original bank band must exclude the old minimum circle')
    area = -hole.area
    perimeter = float(np.linalg.norm(np.roll(hole.xy, -1, axis=0)-hole.xy, axis=1).sum())
    candidates = radius_candidates(area, perimeter, old[-1])
    for index, value in enumerate(candidates):
        try:
            corner = _corner(a, b, value, math.radians(settings.continuation_angle_degrees))
            if (corner.center is None or not all(regular_contact(p, s) for p, s in zip((a, b), corner.stations))):
                continue
            corner_patch(Key('thin_median_proof', 0), corner, settings.spacing, tolerance)
        except ValueError:
            continue
        corner = replace(corner, requested_radius=radius, radius_mode='STYLE',
            radius_attempts=(*old, *candidates[:index+1]), fit_policy='STYLE_ORIGINAL_THIN_MEDIAN')
        report = dict(policy=corner.fit_policy, source_sections=sorted(ids),
            junction=site.text(), opposite_junction=far.text(), original_hole_area=area,
            original_hole_perimeter=perimeter, characteristic_radius_cap=area/perimeter,
            radius_candidates=list(candidates), candidate_attempts=index+1, radius=value,
            requested_radius=radius, old_minimum_radius=old[-1],
            finite_source_stations=list(corner.stations), curvature_bound=1./value,
            available_bank_width=width, old_circle_required_projected_width=required,
            normal=normal.tolist(), axis=axis.tolist(), original_bank_proofs=proofs,
            original_sources_unchanged=True, source_curve_tolerance=tolerance)
        corner.thin_median_gore = ThinMedianGore(hole.xy.copy(), report)
        return corner
    raise ValueError('No bounded STYLE thin-median radius fits its characteristic width and finite banks')


def certify(corners, sources, loops):
    """Require each whole nose on the retained true hole after normal union."""
    from .compounds import rim_label
    from .finite_containment import contains_segment
    inventory = {source.region.owner: source for source in sources}
    for corner in corners:
        record = getattr(corner, 'thin_median_gore', None)
        if record is None: continue
        owners = {s.region.owner for s in sources if any(c is corner for c in s.layout.corners)}
        exposed = []
        for loop in loops:
            if not loop.hole: continue
            amount = sum(abs(rim_label(span, inventory).last-rim_label(span, inventory).first)
                for span in loop.spans if span.owner in owners and rim_label(span, inventory).kind == 'curb')
            if amount > 1e-8: exposed.append((loop, amount))
        if len(exposed) != 1 or abs(exposed[0][1]-1.) > 1e-8:
            raise ValueError('A thin median nose must expose its whole original tangent-to-tangent circle on one true hole')
        loop = exposed[0][0]
        tolerance = 2*record.report['source_curve_tolerance'] + 1e-7
        if not all(contains_segment(a, b, record.source_hole, tolerance=tolerance)
                   for a, b in zip(loop.xy, np.roll(loop.xy, -1, axis=0))):
            raise ValueError('A rounded thin median must retain the genuine original hole')
        record.report.update(retained_hole_area=-loop.area,
            retained_hole_fraction=-loop.area/record.report['original_hole_area'],
            whole_inner_circle_exposed=True, hole_containment_allowance=tolerance)
