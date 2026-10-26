"""Small input-driven regression networks, shared by Python/native checks."""
import math
from dataclasses import replace


def ownership_fixture_networks(model):
    pair = short_junction_pair(model)
    pair = replace(pair, uuid='compound_short_pair',
        furniture=(model.Placement(81, 1, 3., offset=-2.),))
    cross = compound_cross_network(model)
    section = next(s for s in cross.sections if s.id == 8)
    a, b = section.points
    edited = replace(section, points=(a, tuple(.5 * (x + y) for x, y in zip(a, b)), b))
    cross = replace(cross, uuid='compound_divided_cross',
        sections=tuple(edited if s.id == 8 else s for s in cross.sections),
        furniture=(model.Placement(81, 8, 5.),))
    return pair, cross, merging_bundle_network(model)


def height_fixture_networks(model):
    """A normal input network exercising projected junctions and attachments."""
    import numpy as np
    origin = np.array((464.9602355957031, 3.907151222229004, 0.))
    south = np.array(((455.3344421386719, -67.87256622314453, 0.),
                      (456.4639587402344, -49.12625503540039, 0.), origin)) - origin
    north = np.array((origin, (463.63800048828125, 17.97780990600586, 0.),
                      (458.25830078125, 31.246585845947266, 0.))) - origin
    east = np.array((origin, (472.9197082519531, 2.6054763793945312, 0.),
                     (486.00201416015625, .20298385620117188, 0.))) - origin
    delta = np.array((453.810302734375, 5.186278343200684, 0.)) - origin
    west = 20. * delta / np.linalg.norm(delta)
    down = 25. * (south[0] - south[1]) / np.linalg.norm(south[0] - south[1])
    up = 25. * (north[-1] - north[-2]) / np.linalg.norm(north[-1] - north[-2])
    # Two divided corridors meet a narrow crossroad. Their profiles retain the
    # half/full-crown combination of a real junction, with complete free ends.
    data = ((south, 5.4), (south + delta, 5.4), (north, 5.4), (north + delta, 5.4),
        ((delta, np.zeros(3)), 2.5), (east, 5.), ((delta, delta + west), 2.5),
        ((south[0] + delta, south[0]), 5.), ((south[0], south[0] + down), 5.4),
        ((south[0] + delta, south[0] + delta + down), 5.4),
        ((north[-1] + delta, north[-1]), 5.), ((north[-1], north[-1] + up), 5.4),
        ((north[-1] + delta, north[-1] + delta + up), 5.4))
    sections = tuple(model.Section(i + 1, i + 1, tuple(tuple(map(float, p)) for p in points), width=width)
                     for i, (points, width) in enumerate(data))
    incidence = (((1, 1), (3, 0), (5, 1), (6, 0)), ((2, 1), (4, 0), (5, 0), (7, 0)),
        ((1, 0), (8, 1), (9, 0)), ((2, 0), (8, 0), (10, 0)),
        ((3, 1), (11, 1), (12, 0)), ((4, 1), (11, 0), (13, 0)))
    junctions = tuple(model.Intersection(i + 1, tuple(model.End(*end) for end in ends))
                      for i, ends in enumerate(incidence))
    bundles = tuple(model.Bundle(i + 1, tuple(model.BundleMember((model.Traversal(sid),)) for sid in members))
                    for i, members in enumerate(((2, 1), (4, 3))))
    return (model.Network('profiled_corridor_crossing', sections, junctions, bundles,
                         furniture=(model.Placement(81, 3, 1., offset=-2.7),)),)


def transition_fixture_networks(model):
    """Ordinary and divided Streets with authored dimensional Section seams."""
    sections = tuple(model.Section(i + 1, i // 2 + 1,
        ((x0, y, 0.), (.5 * (x0 + x1), y, 0.), (x1, y, 0.)), width=width)
        for i, (x0, x1, y, width) in enumerate(((-40., 0., 8., 5.4), (0., 40., 8., 8.1),
                                                (-40., 0., -8., 5.4), (0., 40., -8., 8.1))))
    connections = tuple(model.Connection(6 + i, (model.End(2 * i + 1, 1), model.End(2 * i + 2, 0))) for i in range(2))
    single = model.Network('dimensional_seam', sections[:2], connections=connections[:1],
        furniture=(model.Placement(81, 1, 39., offset=-1.),))
    members = tuple(model.BundleMember(tuple(model.Traversal(sid) for sid in ids)) for ids in ((1, 2), (3, 4)))
    divided = model.Network('divided_dimensional_seams', sections, bundles=(model.Bundle(7, members),),
        connections=connections, furniture=(model.Placement(81, 1, 39.),))
    curved = model.Network('curved_dimensional_seam', (
        model.Section(1, 7, ((10.79864501953125, 63.49762940406799, 0.),
                            (1.279632568359375, 7.347206115722656, 0.), (0., 0., 0.)), width=5.4),
        model.Section(2, 7, ((0., 0., 0.), (-4.32623291015625, -26.082542419433594, 0.),
                            (-5.131927490234375, -30.947280883789062, 0.)), width=8.1)),
        connections=connections[:1])
    return single, divided, curved


def clearance_fixture_networks(model):
    """A junction close to a width transition, plus an unrelated control road."""
    sections = (model.Section(1, 1, ((0., 0., 0.), (3.75, 0., 0.), (7.5, 0., 0.)), width=5.4),
        model.Section(2, 1, ((7.5, 0., 0.), (20., 0., 0.), (40., 0., 0.)), width=8.1, offset=1.35),
        model.Section(3, 3, ((0., 0., 0.), (0., 12.5, 0.), (0., 25., 0.)), width=5.4),
        model.Section(4, 4, ((0., 0., 0.), (-12.5, 0., 0.), (-25., 0., 0.)), width=5.4),
        model.Section(5, 5, ((-25., -25., 0.), (0., -25., 0.), (25., -25., 0.))))
    return (model.Network('junction_dimensional_clearance', sections,
        (model.Intersection(1, (model.End(1, 0), model.End(3, 0), model.End(4, 0))),),
        connections=(model.Connection(6, (model.End(1, 1), model.End(2, 0))),),
        furniture=(model.Placement(81, 1, 6.8, offset=-1.),)),)


def absorbed_transition_fixture_networks(model):
    """A short inner Section changes width inside a rounded T intersection."""
    return (model.Network('junction_owned_width_seam', (
        model.Section(1, 1, ((0., 0., 0.), (1., 0., 0.), (2., 0., 0.))),
        model.Section(2, 1, ((2., 0., 0.), (20., 0., 0.), (40., 0., 0.)), width=6.05),
        model.Section(3, 3, ((0., 0., 0.), (0., 20., 0.), (0., 40., 0.))),
        model.Section(4, 4, ((0., 0., 0.), (-20., 0., 0.), (-40., 0., 0.))),
        model.Section(99, 99, ((-40., -40., 0.), (0., -40., 0.), (40., -40., 0.)))),
        (model.Intersection(1, (model.End(1, 0), model.End(3, 0), model.End(4, 0))),),
        connections=(model.Connection(6, (model.End(1, 1), model.End(2, 0))),),
        furniture=(model.Placement(81, 1, 1., offset=-1.),
                   model.Placement(82, 2, .5, offset=-1.), model.Placement(83, 2, 15., offset=-1.))),)


def absorbed_profile_fixture_networks(model):
    """The rounded junction owns a rise-only seam into a lower/higher crown."""
    base = absorbed_transition_fixture_networks(model)[0]
    return tuple(replace(base, uuid='junction_owned_' + name + '_crown',
        sections=tuple(replace(s, width=6., rise=rise) if s.id == 2 else s for s in base.sections))
        for name, rise in (('lower', .04), ('higher', .12)))


def side_junction_fixture_networks(model):
    """An outside turnoff interrupts one member but leaves the median whole."""
    sections = (model.Section(1, 1, ((-40., 8., 0.), (-20., 8., 0.), (0., 8., 0.)), width=5.4),
        model.Section(2, 1, ((0., 8., 0.), (20., 8., 0.), (40., 8., 0.)), width=5.4),
        model.Section(3, 2, ((-40., -8., 0.), (0., -8., 0.), (40., -8., 0.)), width=5.4),
        model.Section(4, 4, ((0., 8., 0.), (0., 22., 0.), (0., 36., 0.)), width=5.4))
    bundle = model.Bundle(7, (model.BundleMember((model.Traversal(1), model.Traversal(2))),
                             model.BundleMember((model.Traversal(3),))))
    outward = model.Network('bundle_outward_side_junction', sections,
        (model.Intersection(1, (model.End(1, 1), model.End(2, 0), model.End(4, 0))),),
        bundles=(bundle,))
    cross_sections = (*sections[:2],
        replace(sections[2], points=((-40., -8., 0.), (-20., -8., 0.), (0., -8., 0.))),
        replace(sections[3], points=((0., 8., 0.), (0., 0., 0.), (0., -8., 0.))),
        model.Section(5, 2, ((0., -8., 0.), (20., -8., 0.), (40., -8., 0.)), width=5.4))
    cross_bundle = model.Bundle(7, (bundle.members[0],
        model.BundleMember((model.Traversal(3), model.Traversal(5)))), interior_section_ids=(4,))
    crossing = model.Network('bundle_crossing_side_junctions', cross_sections,
        (outward.intersections[0], model.Intersection(2, (model.End(3, 1), model.End(5, 0), model.End(4, 1)))),
        bundles=(cross_bundle,))
    return outward, crossing


def interior_fixture_networks(model):
    """One member ends at an interior crossing; a separate Street continues."""
    network = side_junction_fixture_networks(model)[1]
    bundle = replace(network.bundles[0], members=(network.bundles[0].members[0],
        model.BundleMember((model.Traversal(3),))))
    sections = tuple(replace(s, street_id=5, rise=.12) if s.id == 5 else s for s in network.sections)
    remote = model.Section(99, 99, ((-40., 40., 0.), (0., 40., 0.), (40., 40., 0.)))
    return (replace(network, uuid='bundle_crossing_free_continuation', bundles=(bundle,),
        sections=(*sections, remote), furniture=(model.Placement(82, 5, 20.),)),)


def consumed_cap_fixture_networks(model):
    """A rounded cap uses an actual port beyond a fully consumed end Section."""
    sections = (
        model.Section(1, 1, ((-40., 8., 0.), (-20., 8., 0.), (0., 8., 0.)), width=5.4),
        model.Section(2, 2, ((-40., -8., 0.), (-20., -8., 0.), (-1., -8., 0.)), width=5.4),
        model.Section(3, 2, ((-1., -8., 0.), (0., -8., 0.)), width=5.4),
        model.Section(4, 4, ((0., 8., 0.), (0., 0., 0.), (0., -8., 0.)), width=5.4),
        model.Section(5, 5, ((0., -8., 0.), (20., -8., 0.)), width=5.4),
        model.Section(6, 6, ((0., 8., 0.), (20., 8., 0.)), width=5.4),
        model.Section(99, 99, ((-40., 40., 0.), (0., 40., 0.), (40., 40., 0.))))
    bundle = model.Bundle(7, (model.BundleMember((model.Traversal(1),)),
        model.BundleMember((model.Traversal(2), model.Traversal(3)))))
    junctions = (model.Intersection(1, (model.End(1, 1), model.End(4, 0), model.End(6, 0))),
                 model.Intersection(2, (model.End(3, 1), model.End(4, 1), model.End(5, 0))))
    return (model.Network('bundle_consumed_terminal_cap', sections, junctions, (bundle,),
        (model.Placement(83, 2, 20.),), connections=(model.Connection(3, (model.End(2, 1), model.End(3, 0))),)),)


def locality_fixture_networks(model):
    """A potential far-junction route does not trim an unchanged Bundle rim."""
    network = side_junction_fixture_networks(model)[0]
    sections = tuple(replace(s, points=(s.points[0], s.points[1])) if s.id == 4 else s
                     for s in network.sections)
    sections += (model.Section(6, 4, ((0., 22., 0.), (0., 36., 0.)), width=5.4),
                 model.Section(7, 4, ((0., 36., 0.), (0., 50., 0.)), width=5.4),
                 model.Section(8, 8, ((0., 36., 0.), (20., 36., 0.)), width=5.4))
    return (replace(network, uuid='bundle_remote_mouth_locality', sections=sections,
        intersections=(*network.intersections, model.Intersection(2,
            (model.End(6, 1), model.End(7, 0), model.End(8, 0)))),
        connections=(model.Connection(3, (model.End(4, 1), model.End(6, 0))),),
        furniture=(model.Placement(84, 3, 20., side='LEFT'),)),)


def merging_bundle_network(model):
    sections = (model.Section(1, 1, ((0., 0., 0.), (-5., 20., 0.), (-10., 40., 0.))),
                model.Section(2, 2, ((0., 0., 0.), (5., 20., 0.), (10., 40., 0.))),
                model.Section(3, 3, ((0., 0., 0.), (0., -35., 0.))))
    return model.Network('merging_bundle', sections,
        (model.Intersection(1, tuple(model.End(s.id, 0) for s in sections)),),
        bundles=(model.Bundle(1, (model.BundleMember((model.Traversal(1),)),
                                  model.BundleMember((model.Traversal(2),)))),))


def compound_cross_network(model):
    """Two divided roads crossing at four sites; preparation fixture only."""
    sections, vertices, street_sections = [], {}, {}
    for axis in (0, 1):
        for side in (-5., 5.):
            ids = []
            for a, b in ((-40., -5.), (-5., 5.), (5., 40.)):
                def point(value):
                    return (value, side, 0.) if axis == 0 else (side, value, 0.)
                sid = len(sections) + 1
                street = (axis, side)
                sections.append(model.Section(sid, len(street_sections) + 1, (point(a), point(b))))
                ids.append(sid)
                for end, p in enumerate((point(a), point(b))):
                    vertices.setdefault(p, []).append(model.End(sid, end))
            street_sections[(axis, side)] = tuple(model.Traversal(sid) for sid in ids)
    junctions = tuple(model.Intersection(i + 1, tuple(ends)) for i, (point, ends) in enumerate(sorted(vertices.items())) if len(ends) > 2)
    bundles = (model.Bundle(1, (model.BundleMember(street_sections[(1, -5.)]), model.BundleMember(street_sections[(1, 5.)]))),
               model.Bundle(2, (model.BundleMember(street_sections[(0, 5.)]), model.BundleMember(street_sections[(0, -5.)]))))
    return model.Network('compound_divided_cross_preparation', tuple(sections), junctions, bundles)


def short_junction_pair(model, length=6.):
    sections = (model.Section(1, 1, ((0., 0., 0.), (length, 0., 0.))),
        model.Section(2, 2, ((0., 0., 0.), (-30., 0., 0.))),
        model.Section(3, 3, ((0., 0., 0.), (0., 30., 0.))),
        model.Section(4, 4, ((length, 0., 0.), (length + 30., 0., 0.))),
        model.Section(5, 5, ((length, 0., 0.), (length, -30., 0.))))
    return model.Network('short_pair_preparation', sections,
        (model.Intersection(1, (model.End(1, 0), model.End(2, 0), model.End(3, 0))),
         model.Intersection(2, (model.End(1, 1), model.End(4, 0), model.End(5, 0)))))


def compound_seam_fixture(model, *, width=6.4, rise=.12, mode=3, length=6., split=3.):
    """Both close junctions consume the connecting Street's dimensional seam."""
    n = short_junction_pair(model, length=length)
    original = n.section_index[1]
    first = replace(original, points=((0., 0., 0.), (split / 2., 0., 0.), (split, 0., 0.)))
    second = replace(original, id=6, points=((split, 0., 0.), ((split + length) / 2., 0., 0.),
                                          (length, 0., 0.)), width=width, rise=rise)
    j2 = replace(n.intersections[1], ends=(model.End(6, 1), *n.intersections[1].ends[1:]))
    return replace(n, uuid='compound_owned_profile_seam',
        sections=tuple(replace(s, interpolation=mode) for s in
            (first, *n.sections[1:], second, model.Section(99, 99, ((-30., -40., 0.), (30., -40., 0.))))),
        intersections=(n.intersections[0], j2),
        connections=(model.Connection(6, (model.End(1, 1), model.End(6, 0))),),
        furniture=(model.Placement(81, 1, 1.5, offset=-1.), model.Placement(82, 6, 1.5, offset=1.)))


def compound_seam_fixture_networks(model):
    base = compound_seam_fixture(model)
    return (replace(base, uuid='compound_owned_width_seam', sections=tuple(
                replace(s, rise=.08) if s.id == 6 else s for s in base.sections)), base)


def compound_route_fixture_networks(model):
    return (replace(compound_seam_fixture(model, width=6., rise=.08), uuid='compound_split_same_profile'),
            replace(compound_seam_fixture(model, rise=.5), uuid='compound_split_internal_crown'),
            replace(compound_seam_fixture(model, rise=.5, length=10., split=2.), uuid='compound_short_prefix'))


def compound_internal_crown_fixture(model):
    """A high crown wholly inside a compound must not raise its exterior ceiling."""
    base = short_junction_pair(model)
    return replace(base, uuid='compound_internal_crown',
        sections=tuple(replace(s, rise=.18) if s.id == 1 else s for s in
            (*base.sections, model.Section(99, 99, ((-30., -40., 0.), (30., -40., 0.))))),
        furniture=(model.Placement(81, 1, 1.5, offset=-1.), model.Placement(82, 1, 4.5, offset=1.)))


def seamed_approach_network(model, *, reverse_first=False, reverse_next=False, curved=False):
    first = ((0., 0., 0.), (0., 2., 0.))
    next_points = ((0., 2., 0.), (1., 20., 0.), (0., 40., 0.)) if curved else ((0., 2., 0.), (0., 40., 0.))
    sections = (model.Section(1, 10, first[::-1] if reverse_first else first),
        model.Section(2, 20, ((0., 0., 0.), (-40., 0., 0.))),
        model.Section(3, 30, ((0., 0., 0.), (40., 0., 0.))),
        model.Section(4, 10, next_points[::-1] if reverse_next else next_points),
        model.Section(5, 50, ((60., 0., 0.), (80., 0., 0.))))
    return model.Network('seamed_approach', sections,
        (model.Intersection(10, (model.End(1, int(reverse_first)), model.End(2, 0), model.End(3, 0))),),
        connections=(model.Connection(20, (model.End(1, 0 if reverse_first else 1),
                                           model.End(4, 1 if reverse_next else 0))),),
        furniture=(model.Placement(81, 4, 15.),))


def divided_network(model, separation=20.):
    left, right = -.5 * separation, .5 * separation
    sections = (model.Section(1, 1, ((left, 0., 0.), (left, 22.5, 0.), (left, 45., 0.))),
        model.Section(2, 2, ((right, 0., 0.), (right, 22.5, 0.), (right, 45., 0.))),
        model.Section(3, 3, ((left, 0., 0.), (right, 0., 0.))),
        model.Section(4, 4, ((left, 0., 0.), (left - 35., 0., 0.))),
        model.Section(5, 5, ((right, 0., 0.), (right + 35., 0., 0.))))
    return model.Network('divided_T', sections,
        (model.Intersection(100, (model.End(1, 0), model.End(3, 0), model.End(4, 0))),
         model.Intersection(200, (model.End(2, 0), model.End(3, 1), model.End(5, 0)))),
        (model.Bundle(7, (model.BundleMember((model.Traversal(1),)),
                          model.BundleMember((model.Traversal(2),)))),))


def fixture_networks(model):
    cases = (('T', (0, 90, 180), (6, 6, 6)),
             ('three_way', (0, 120, 240), (6, 6, 6)),
             ('oblique', (0, 65, 180, 255), (6, 6, 6, 6)),
             ('five_way', (0, 72, 144, 216, 288), (6, 6, 6, 6, 6)),
             ('unequal_widths', (0, 90, 180, 270), (6, 12, 6, 12)))
    result = []
    for name, angles, widths in cases:
        sections = tuple(model.Section(i, i, ((0., 0., 0.),
            (45. * math.cos(math.radians(a)), 45. * math.sin(math.radians(a)), 0.)),
            width=width, rise=.08 if i % 2 == 0 else .06)
            for i, (a, width) in enumerate(zip(angles, widths)))
        result.append(model.Network(name, sections,
            (model.Intersection(1, tuple(model.End(s.id, 0) for s in sections)),)))
    result.append(divided_network(model))
    from _python_bundle_fixtures import bundle_networks
    # Keep Exact-mode qualification in the numerical suite. The combined
    # inspection scene uses the requested Smooth default on its authored inputs.
    result.extend(bundle_networks(model, interpolation=3))
    result.append(seamed_approach_network(model, curved=True))
    sections = (model.Section(1, 1, ((0., 0., 0.), (20., 2., 0.), (40., 0., 0.))),
        model.Section(2, 2, ((0., 0., 0.), (-30., 0., 0.))),
        model.Section(3, 3, ((0., 0., 0.), (0., 30., 0.))),
        model.Section(4, 4, ((40., 0., 0.), (70., 0., 0.))),
        model.Section(5, 5, ((40., 0., 0.), (40., -30., 0.))))
    result.append(model.Network('bent_two_intersections', sections,
        (model.Intersection(1, (model.End(1, 0), model.End(2, 0), model.End(3, 0))),
         model.Intersection(2, (model.End(1, 1), model.End(4, 0), model.End(5, 0)))),
        furniture=(model.Placement(81, 1, 18.), model.Placement(82, 1, 25., side='LEFT'))))
    return tuple(result)


def terminal_fixture_networks(model):
    sections = (model.Section(1, 1, ((0., 0., 0.), (4.968, 0., 0.)), width=5.4),
                model.Section(2, 2, ((0., 0., 0.), (0., 30., 0.)), width=5.4),
                model.Section(3, 3, ((0., 0., 0.), (-30., 0., 0.)), width=5.4))
    junction = model.Intersection(10, tuple(model.End(s.id, 0) for s in sections), radius=2.)
    placement = model.Placement(85, 1, 4.96, side='RIGHT', offset=-.5)
    return (model.Network('terminal_profiled_stub', sections, intersections=(junction,), furniture=(placement,)),)


def radius_fixture_networks(model):
    sections = tuple(model.Section(i, i, ((0., 0., 0.), end), width=6.) for i, end in
                     enumerate(((30., 0., 0.), (0., 8.5, 0.), (-30., 0., 0.)), 1))
    remote = model.Section(4, 4, ((60., 0., 0.), (80., 0., 0.)))
    junction = model.Intersection(1, tuple(model.End(s.id, 0) for s in sections), radius=10.,
                                  radius_modes=('STYLE',) * 3)
    return (model.Network('preferred_radius_short_T', (*sections, remote), intersections=(junction,)),)


def clipped_fixture_source(model, *, reverse=False, interpolation=3, curved=False):
    """A junction arm is cut at 2.5 m; authentic same-way support bends outside."""
    from types import SimpleNamespace
    tags = {'highway': 'residential', 'name': 'Original source road'}
    north = ((0., 0., 0.), (0., 2.5, 0.))
    sections = (model.Section(1, 1, ((0., 0., 0.), (30., 0., 0.)), interpolation=interpolation),
                model.Section(2, 2, north[::-1] if reverse else north, interpolation=interpolation),
                model.Section(3, 3, ((0., 0., 0.), (-30., 0., 0.)), interpolation=interpolation),
                model.Section(4, 4, ((60., 0., 0.), (80., 0., 0.)), interpolation=interpolation))
    ends = (model.End(1, 0), model.End(2, int(reverse)), model.End(3, 0))
    junction = model.Intersection(1, ends, radius=2., radius_modes=('STYLE',) * 3)
    network = model.Network('clipped_original_support_T', sections, intersections=(junction,),
                            furniture=(model.Placement(85, 2, 1., offset=-.5),))
    cut_end = 0 if reverse else 1
    interval = (.625, 0.) if reverse else (0., .625)
    streets = {s.street_id: dict(items=[dict(kind='Section', id=s.id, tags=tags,
        source_segments=[('way', 'root', 'outside', *interval)] if s.id == 2 else [])]) for s in sections}
    nodes = [dict(osm_node_id='outside', position=(0., 4.), topology_kind='CONTINUATION',
        incidences=[dict(source_way_id='way', source_nodes=('root', 'outside'), parameter=.625,
                        source_end=1, clipped_point=(0., 2.5))],
        open_ends=[dict(section_id=2, end=cut_end, street_id=2, owner_scope='ACTIVE',
                       source_way_id='way', source_nodes=('root', 'outside'), parameter=.625)],
        legs=[dict(source_way_id='way', source_nodes=('root', 'outside'), source_end=1,
                   other_node_id='root', other_position=(0., 0.), category='residential', tags=tags),
              dict(source_way_id='way', source_nodes=('outside', 'far'), source_end=0,
                   other_node_id='far', other_position=(2. if curved else 0., 30.), category='residential', tags=tags)])]
    return network, SimpleNamespace(streets=streets), dict(version=1, nodes=nodes, issues=[])


def attached_overlap_fixture_networks(model):
    upper = ((-40., 8., 0.), (-28., 8., 0.), (-16., 8., 0.), (0., 2., 0.),
             (16., 8., 0.), (28., 8., 0.), (40., 8., 0.))
    lower = tuple((x, -y, z) for x, y, z in upper)
    sections = [model.Section(1, 1, upper, width=6., rise=.12),
                model.Section(2, 2, lower, width=6., rise=.12)]
    for first, x, street in ((3, -40., 3), (6, 40., 4)):
        for i, ys in enumerate(((8., -8.), (8., 28.), (-8., -28.))):
            sections.append(model.Section(first + i, street, tuple((x, y, 0.) for y in ys), width=5.))
    sections.append(model.Section(99, 99, ((-45., 45., 0.), (45., 45., 0.)), width=5.))
    incidence = (((1, 0), (3, 0), (4, 0)), ((2, 0), (3, 1), (5, 0)),
                 ((1, 1), (6, 0), (7, 0)), ((2, 1), (6, 1), (8, 0)))
    sites = tuple(model.Intersection(i + 1, tuple(model.End(*v) for v in ends)) for i, ends in enumerate(incidence))
    bundle = model.Bundle(7, tuple(model.BundleMember((model.Traversal(i),)) for i in (1, 2)))
    return (model.Network('attached_curved_overlap', tuple(sections), sites, (bundle,),
        furniture=(model.Placement(81, 1, 40., offset=-2.5), model.Placement(82, 1, 2., 'LEFT', -3.),
                   model.Placement(83, 2, 40., 'LEFT', -2.5), model.Placement(99, 99, 30., offset=-1.))),)


def source_support_fixture_networks(model):
    import importlib
    cls = importlib.import_module(model.__package__ + '.source_support').ClippedApproachSupport
    network, source, context = clipped_fixture_source(model, reverse=True, curved=True)
    prepared = cls(source, context).prepare(network)
    if prepared.error:
        raise ValueError(prepared.error)
    return (prepared.network,)


def split_attached_fixture(model, *, short_ends=False):
    """The same Bundle with explicit Section seams, including consumed ends."""
    from dataclasses import replace
    network = attached_overlap_fixture_networks(model)[0]
    cut = 38. if short_ends else 28.
    added, heads, connections = [], [], []
    for role, original in enumerate(network.sections[:2]):
        first, middle, last = (1, 11, 12) if role == 0 else (2, 21, 22)
        y = 8. if role == 0 else -8.
        start, finish = (-cut, y, 0.), (cut, y, 0.)
        interior = tuple(p for p in original.points if -cut < p[0] < cut)
        heads.append(replace(original, points=(original.points[0], start)))
        added.extend((replace(original, id=middle, points=(start, *interior, finish)),
                      replace(original, id=last, points=(finish, original.points[-1]))))
        connections.extend((model.Connection(40 + 2 * role, (model.End(first, 1), model.End(middle, 0))),
                            model.Connection(41 + 2 * role, (model.End(middle, 1), model.End(last, 0)))))
    members = tuple(replace(member, path=tuple(model.Traversal(sid) for sid in ids))
        for member, ids in zip(network.bundles[0].members, ((1, 11, 12), (2, 21, 22))))
    sites = tuple(replace(j, ends=tuple(model.End(12 if e.section_id == 1 else 22, 1)
        if e.section_id in (1, 2) and e.end == 1 else e for e in j.ends)) for j in network.intersections)
    furniture = tuple(replace(p, section_id=11 if p.id == 81 else 21, station=p.station - (40. - cut))
        if p.id in (81, 83) else p for p in network.furniture)
    return replace(network, uuid='attached_consumed_section_chain' if short_ends else 'attached_section_chain',
        sections=(*heads, *network.sections[2:], *added), intersections=sites,
        connections=tuple(connections), bundles=(replace(network.bundles[0], members=members),), furniture=furniture)


def multisection_attached_fixture_networks(model):
    return tuple(split_attached_fixture(model, short_ends=value) for value in (False, True))


def skew_bundle_fixture(model, *, move=.2, split=False, free=False):
    from dataclasses import replace
    network = split_attached_fixture(model,short_ends=True) if split else attached_overlap_fixture_networks(model)[0]
    sid = 11 if split else 1
    sections = tuple(replace(s,points=tuple((x,y+move if i==1 else y,z)
        for i,(x,y,z) in enumerate(s.points))) if s.id==sid else s for s in network.sections)
    if free:
        assert not split
        sections = tuple(s for s in sections if s.id in (1,2,99))
    return replace(network,uuid='skew_free_overlap' if free else 'skew_attached_overlap',
        sections=sections,intersections=() if free else network.intersections,
        furniture=(*network.furniture,model.Placement(84,sid,.2 if free else 5.1,offset=-3.)))


def skew_bundle_fixture_networks(model):
    return tuple(skew_bundle_fixture(model,free=free) for free in (False,True))


def dimensional_bundle_fixture(model, *, bent=False, attached=False, skew=False):
    """Width transitions close a median; original Section IDs remain editable."""
    if attached:
        network = split_attached_fixture(model)
        return replace(network, uuid='attached_dimensional_bundle',
            sections=tuple(replace(s, width=8.) if s.id in (11, 21) else s for s in network.sections),
            furniture=(*network.furniture, model.Placement(85, 1, 11., 'RIGHT', -1.),
                       model.Placement(86, 11, 1., 'LEFT', -1.)))
    sections = []
    for member, y in enumerate((4., -4.)):
        for part, (start, end, width) in enumerate(((-40., 0., 6.), (0., 40., 10.))):
            sid = 1 + 2 * member + part
            xs = tuple(start + (end-start)*i/4. for i in range(5))
            points = tuple((float(x), y + ((.3 if member == 0 else .5) * (1. if part == 0 else -1.)
                if bent and abs(x) == 20. else 0.), 0.) for x in xs)
            sections.append(model.Section(sid, member+1, points, width=width))
    remote = model.Section(99, 99, ((-40., 35., 0.), (40., 35., 0.)))
    members = tuple(model.BundleMember(tuple(model.Traversal(sid) for sid in ids)) for ids in ((1, 2), (3, 4)))
    network = model.Network('bent_dimensional_bundle' if bent else 'dimensional_bundle', tuple(sections)+(remote,),
        bundles=(model.Bundle(7, members),),
        connections=(model.Connection(40, (model.End(1, 1), model.End(2, 0))),
                     model.Connection(41, (model.End(3, 1), model.End(4, 0)))),
        furniture=(model.Placement(81, 1, 39., 'RIGHT', -1.), model.Placement(82, 2, 1., 'LEFT', -1.),
                   model.Placement(99, 99, 20., 'RIGHT', -1.)))
    if skew:
        upper = (((-40.,4.,0.),(-28.,4.2,0.),(-16.,1.7,0.),(0.,1.7,0.)),
                 ((0.,1.7,0.),(16.,1.7,0.),(28.,4.,0.),(40.,4.,0.)))
        network = replace(network, uuid='skew_dimensional_bundle',
            sections=tuple(replace(s, points=tuple((x,y if s.id<3 else -y,z) for x,y,z in upper[(s.id-1)%2]),
                width=3. if s.id%2 else 5.) if s.id!=99 else s for s in network.sections),
            furniture=(model.Placement(81,1,39.,'RIGHT',-1.),model.Placement(82,2,1.,'LEFT',-1.),
                       model.Placement(99,99,20.,'RIGHT',-1.)))
    return network


def dimensional_bundle_fixture_networks(model):
    return (dimensional_bundle_fixture(model), dimensional_bundle_fixture(model, bent=True),
            dimensional_bundle_fixture(model, attached=True))


def one_sided_free_cap_fixture(model, *, move=.05, reversed_direction=False):
    """A captured valid edit: skew separated profiles opposite a merged free cap."""
    network=dimensional_bundle_fixture(model)
    network=replace(network,uuid='one_sided_free_cap',
        sections=tuple(replace(s,points=tuple((x,y+move if i==1 else y,z)
            for i,(x,y,z) in enumerate(s.points))) if s.id==1 else s for s in network.sections),
        furniture=(*network.furniture,model.Placement(84,1,.2,'RIGHT',-1.)))
    if reversed_direction:
        bundle=network.bundles[0]
        members=tuple(replace(m,path=tuple(replace(t,forward=not t.forward) for t in m.path[::-1]))
                      for m in bundle.members[::-1])
        network=replace(network,uuid='one_sided_free_cap_reversed',bundles=(replace(bundle,members=members),))
    return network


def one_sided_free_cap_fixture_networks(model):
    return (one_sided_free_cap_fixture(model),)


def nonparallel_terminal_fixture(model, *, move=.05):
    """Move a control at the merged free end without changing its endpoint."""
    network = one_sided_free_cap_fixture(model)
    return replace(network, uuid='nonparallel_merged_terminal',
        sections=tuple(replace(s, points=tuple((x, y+move if i==len(s.points)-2 else y, z)
            for i,(x,y,z) in enumerate(s.points))) if s.id==2 else s for s in network.sections),
        furniture=(*network.furniture, model.Placement(85, 2, 39.8, 'RIGHT', -2.)))


def higher_valence_fixture_networks(model):
    """Even/odd higher valence with independent widths and crown rises."""
    result = []
    for name, angles in (('six_way', tuple(range(0, 360, 60))),
                         ('seven_way', (0, 50, 105, 150, 204, 257, 307))):
        sections = tuple(model.Section(i+1, i+1,
            tuple((r*math.cos(math.radians(angle)), r*math.sin(math.radians(angle)), 0.)
                  for r in (0., 20., 40.)),
            width=(4., 6., 7.5, 5., 6., 4.5, 7.)[i], rise=.05+.005*(i % 4))
            for i, angle in enumerate(angles))
        result.append(model.Network(name, sections,
            intersections=(model.Intersection(1, tuple(model.End(s.id, 0) for s in sections)),)))
    return tuple(result)


def two_merged_terminal_fixture(model, *, bent_seams=False):
    """Two free union caps surround one island; keep the bent-seam capture."""
    sections, connections = [], []
    for member, y in enumerate((4., -4.)):
        for part, (a, b, width) in enumerate(((-40., -15., 10.), (-15., 15., 6.), (15., 40., 10.))):
            sid = member*3+part+1
            points = ((a, y, 0.), (.5*(a+b), y, 0.), (b, y, 0.))
            if member == 0 and part == 0:
                points = (((a, y, 0.), (.5*(a+b), y+.05, 0.), (b, y, 0.)) if bent_seams else
                          ((a, y, 0.), (a+5., y+.05, 0.), (a+10., y, 0.), (b, y, 0.)))
            if member == 1 and part == 2:
                points = (((a, y, 0.), (.5*(a+b), y+.05, 0.), (b, y, 0.)) if bent_seams else
                          ((a, y, 0.), (b-10., y, 0.), (b-5., y+.05, 0.), (b, y, 0.)))
            sections.append(model.Section(sid, member+1, points, width=width))
            if part:
                connections.append(model.Connection(sid, (model.End(sid-1, 1), model.End(sid, 0))))
    sections.append(model.Section(99, 99, ((60., -25., 0.), (80., -25., 0.), (100., -25., 0.))))
    bundle = model.Bundle(7, tuple(model.BundleMember(tuple(model.Traversal(i) for i in ids))
                                   for ids in ((1, 2, 3), (4, 5, 6))))
    return model.Network('two_merged_terminals_bent_seams' if bent_seams else 'two_merged_terminals',
        tuple(sections), bundles=(bundle,), connections=tuple(connections),
        furniture=(model.Placement(85, 1, .2, 'RIGHT', -2.),
                   model.Placement(86, 6, 24.8, 'LEFT', -2.), model.Placement(99, 99, 20., 'RIGHT', -2.),
                   *((model.Placement(87, 1, 24.5, 'RIGHT', -2.),) if bent_seams else ())))


def exact_bent_seam_fixture(model, *, offsets=(0.,0.), rises=(.08,.08), datum=0.):
    """Ordinary width seam with source bends inside both transition trims."""
    sources = (model.Section(1,10,((-30.,0.,0.),(-12.,-.2,0.),(-6.,.25,0.),(0.,0.,0.)),
                             width=6.,offset=offsets[0],rise=rises[0],datum=datum),
               model.Section(2,10,((0.,0.,0.),(5.,.4,0.),(12.,1.,0.),(30.,0.,0.)),
                             width=10.,offset=offsets[1],rise=rises[1],datum=datum),
               model.Section(99,99,((60.,-25.,0.),(100.,-25.,0.)),datum=datum))
    return model.Network('exact_bent_dimensional_seam',sources,
        connections=(model.Connection(8,(model.End(1,1),model.End(2,0))),),
        furniture=(model.Placement(85,1,29.5,'RIGHT',-1.),model.Placement(99,99,20.,'RIGHT',-1.)))


def sharp_exact_seam_fixture(model, *, length=30.):
    """Explicit Exact-mode regression; the normal Section default stays Smooth."""
    import math
    angle=math.radians(150.)
    sources=(model.Section(1,10,((-length,0.,0.),(0.,0.,0.)),width=6.,interpolation=0),
        model.Section(2,10,((0.,0.,0.),(length*math.cos(angle),length*math.sin(angle),0.)),width=10.,interpolation=0),
        model.Section(99,99,((70.,-25.,0.),(110.,-25.,0.))))
    return model.Network('sharp_exact_seam_'+format(length,'g'),sources,
        connections=(model.Connection(8,(model.End(1,1),model.End(2,0))),),
        furniture=(model.Placement(85,1,length-.5,'RIGHT',-1.),model.Placement(99,99,20.,'RIGHT',-1.)))


def sharp_smooth_seam_fixture(model):
    """The 60 m, 150-degree Smooth seam needs a local curvature proof."""
    original = sharp_exact_seam_fixture(model, length=60.)
    return replace(original, uuid='sharp_smooth_seam_60',
        sections=tuple(replace(s, interpolation=3) for s in original.sections))

