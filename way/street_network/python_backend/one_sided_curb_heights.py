"""Exterior and separator curb conditions for a side street on a divided road."""
import numpy as np

from .compounds import RimSegment
from .intersections import BoundaryCurve
from .model import Key
from .terminal_curb_ancestry import OriginalTerminalBank


def separator_level(context, layout, pair):
    """Prove both halves enclose the same plateau on the two original members."""
    from .regions import signed_area
    from .source_curb_conditions import SourceCurbConditionUnavailable as Unsupported

    first, second = (layout.boundaries[i] for i in pair)
    if (not np.array_equal(first.points[-1], second.points[0])
            or not np.array_equal(second.points[-1], first.points[0])
            or signed_area(np.vstack((first.points[:-1], second.points[:-1]))) >= 0.):
        raise Unsupported('One-sided separator needs a closed clockwise rim')
    paths = {(role, path.section.id): path for approach in layout.approaches
             if hasattr(approach, 'near_role') for role, path in enumerate(approach.paths)}
    near_members = {approach.source.id: approach.near_role for approach in layout.approaches
                    if hasattr(approach, 'near_role')}
    used = set()
    levels = []
    checked_corners = set()
    for curve in (first, second):
        if len(curve.one_sided_edges) != len(curve.points)-1:
            raise Unsupported('One-sided separator lost an original edge')
        for edge in curve.one_sided_edges:
            if isinstance(edge, OriginalTerminalBank):
                key = edge.role, edge.section
                profile = context[Key('profile', edge.section)]
                if (paths.get(key) is not context[Key('path', edge.section)]
                        or profile.half != edge.side):
                    raise Unsupported('One-sided separator needs an original inner bank')
                if key not in used:
                    levels.append(float(profile.datum + profile.value(edge.side)))
                    used.add(key)
            elif isinstance(edge, RimSegment) and edge.kind == 'curb':
                if edge.index in checked_corners:
                    continue
                parent = layout.original_junction_layout.boundaries[edge.index]
                if not isinstance(parent, BoundaryCurve):
                    raise Unsupported('One-sided separator needs an original junction curb')
                for arm, side in ((parent.corner.first, -1), (parent.corner.second, 1)):
                    sid = arm.source.id
                    role = near_members.get(sid)
                    if (role is None or paths[role, sid] is not context[Key('path', sid)]
                            or arm.source is not paths[role, sid].section
                            or arm.profile != context[Key('profile', sid)].oriented(arm.end)
                            or arm.profile.half != side):
                        raise Unsupported('One-sided separator curb leaves the original plateau')
                    levels.append(float(arm.profile.datum + arm.profile.value(side)))
                    used.add((role, sid))
                checked_corners.add(edge.index)
            else:
                raise Unsupported('One-sided separator contains an unclassified edge')
    if {role for role, _ in used} != {0, 1} or max(levels)-min(levels) > 1e-9:
        raise Unsupported('One-sided separator needs two facing banks with a common plateau')
    return levels[0]


def exterior_profile(context, layout, index):
    """Return the common bank level, endpoint crossfalls and transition interval."""
    from .source_curb_conditions import SourceCurbConditionUnavailable as Unsupported

    approaches = layout.approaches
    if index >= len(approaches):
        raise Unsupported('An interior separator needs its own closed boundary profile')
    curve = layout.boundaries[index]
    ends = (index, (index+1) % len(approaches))
    heights, slopes = [], []
    for i, side, point in zip(ends, (-1., 1.), (curve.points[0], curve.points[-1])):
        arm, station = approaches[i], layout.outward_stations[i]
        if np.linalg.norm(point - arm.point(station, arm.offset + .5*side*arm.width)) > 1e-7:
            raise Unsupported('One-sided exterior curb does not connect its two actual mouths')
        heights.append(float(arm.profile.datum + arm.profile.value(side)))
        slopes.append(float(np.linalg.norm(arm.physical_gradient(station, side))))

    # Only a level outer edge is supported by this construction. A separator
    # plateau or a sloping curb requires its own height function.
    arc = []
    sections, corners = set(), set()
    if len(curve.one_sided_edges) != len(curve.points)-1:
        raise Unsupported('One-sided exterior needs one source for every retained edge')
    for i, edge in enumerate(curve.one_sided_edges):
        if isinstance(edge, OriginalTerminalBank):
            profile = context[Key('profile', edge.section)]
            if profile.half != -edge.side:
                raise Unsupported('One-sided exterior contains an inner carriageway bank')
            sections.add(edge.section)
            heights.append(float(profile.datum + profile.value(edge.side)))
            continue
        if not isinstance(edge, RimSegment) or edge.kind != 'curb' or edge.first >= edge.last:
            raise Unsupported('One-sided exterior lost its directed original corner interval')
        parent = layout.original_junction_layout.boundaries[edge.index]
        if not isinstance(parent, BoundaryCurve) or parent.corner.center is None:
            raise Unsupported('One-sided exterior needs a finite original circular corner')
        corners.add(edge.index)
        for arm, side in ((parent.corner.first, -1.), (parent.corner.second, 1.)):
            heights.append(float(arm.profile.datum + arm.profile.value(side)))
        lo, hi = parent.lengths[1:3]/parent.length
        a, b = max(edge.first, lo), min(edge.last, hi)
        if a < b:
            scale = (curve.stations[i+1]-curve.stations[i])/(edge.last-edge.first)
            arc.append(tuple((curve.stations[i]+(p-edge.first)*scale)/curve.length for p in (a,b)))
    if max(heights)-min(heights)>1e-9:
        raise Unsupported('One-sided exterior banks do not share the same level')
    if len(corners) == 1 and arc and len(sections) == 1:
        # Retain the incoming crossfalls along the leads. Reconcile them over
        # just the exposed part of the original circular bend.
        if any(abs(a[1]-b[0]) > 1e-10 for a, b in zip(arc, arc[1:])):
            raise Unsupported('One-sided corner has multiple disjoint exposed arcs')
        phase = arc[0][0], arc[-1][1]
    elif not corners and len(sections) == 1:
        # The uninterrupted opposite bank bridges two independently evaluated
        # paired profiles. Both mouths must identify this very same Section.
        sid, = sections
        if any(not hasattr(approaches[i], 'near_role') or
               approaches[i].paths[1-approaches[i].near_role].section.id != sid for i in ends):
            raise Unsupported('One-sided through bank has no matching paired mouths')
        phase = (0., 1.)
    else:
        raise Unsupported('One-sided exterior needs one continuous bank or one retained bend')
    return heights[0], slopes[0], slopes[1], phase
