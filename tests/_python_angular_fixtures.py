"""Authored curved Bundle inputs shared by numerical and Blender checks.

The caller supplies the model module so fixtures use the same class identities
as the backend being exercised. No addon runtime imports this test helper.
"""
import numpy as np


def angular_network(model, sweep=180., variable=False, mode=3, *, widths=(6., 6.),
                    offsets=(0., 0.), rises=(.08, .08), datum=0., reversed_sources=(),
                    settings=None, displacement=0., rotation=0., translation=(0., 0.)):
    """Two curved members and an unrelated road, with source-relative anchors."""
    t = np.linspace(0., 1., 13)
    angles = np.deg2rad(sweep * (t - .5))
    gap = 4. + (4. * np.sin(np.pi * t) if variable else np.zeros_like(t))
    matrix = np.array(((np.cos(rotation), -np.sin(rotation)),
                       (np.sin(rotation), np.cos(rotation))))
    sections, members = [], []
    for i, (sid, radius) in enumerate(((17, np.full(len(t), 26.)), (29, 26. + gap))):
        points = np.c_[radius * np.cos(angles), radius * np.sin(angles)]
        if i == 1:
            points[6, 0] += displacement
        points = points @ matrix.T + np.asarray(translation)
        forward = sid not in reversed_sources
        if not forward:
            points = points[::-1]
        sections.append(model.Section(
            sid, 41 + 30 * i, tuple(map(tuple, np.c_[points, np.zeros(len(t))])),
            width=widths[i], offset=offsets[i] * (1 if forward else -1), rise=rises[i],
            datum=datum, interpolation=mode))
        members.append(model.BundleMember((model.Traversal(sid, forward),)))
    remote = np.array(((-80., -100.), (80., -100.))) @ matrix.T + np.asarray(translation)
    sections.append(model.Section(99, 99, tuple(map(tuple, np.c_[remote, np.zeros(2)])), datum=datum))
    return model.Network(
        f'angular_{sweep:g}_{int(variable)}_{mode}', tuple(sections),
        bundles=(model.Bundle(13, tuple(members)),),
        furniture=(model.Placement(901, 17, 30., side='LEFT', offset=-1.),
                   model.Placement(999, 99, 40.)), settings=settings or model.Settings())


def angular_fixture_networks(model):
    """The twelve accepted Exact/Smooth constant-overlap and island examples."""
    return tuple(angular_network(model, sweep, variable, mode)
                 for sweep in (180., 225., 270.)
                 for variable in (False, True) for mode in (0, 3))


def angular_chain_network(model, sweep=225., variable=True, mode=3, *,
                          cuts=((4, 8), (3, 7, 10)), reversed_sections=(),
                          **kwargs):
    """Actual Section chains with unequal partitions and explicit incidences.

    Split only the authored inputs. The runtime must retain the separately
    prepared paths; reconstructing one cubic from all controls is not allowed.
    First-piece IDs stay17/29 and later pieces use117/217 and129/229/329.
    """
    from dataclasses import replace
    network = angular_network(model, sweep, variable, mode, **kwargs)
    sections, members, connections, furniture = [], [], [], []
    for role, original in enumerate(network.sections[:2]):
        if not network.bundles[0].members[role].path[0].forward:
            original = replace(original, points=original.points[::-1], offset=-original.offset)
        if mode == 0:
            # Area-free Exact seams are collinear. Put the Section breaks
            # inside source segments; retain every original bend as a control.
            points, joins = [], []
            for i, point in enumerate(original.points):
                points.append(point)
                if i in cuts[role]:
                    joins.append(len(points))
                    points.append(tuple(.5*(np.asarray(point)+original.points[i+1])))
            original = replace(original, points=tuple(points))
            indices = (0, *joins, len(points)-1)
        else:
            indices = (0, *cuts[role], len(original.points)-1)
        visits = []
        for piece, (a, b) in enumerate(zip(indices, indices[1:])):
            sid = original.id + 100*piece
            forward = sid not in reversed_sections
            points = original.points[a:b+1]
            sections.append(replace(original, id=sid,
                points=points if forward else points[::-1],
                offset=original.offset if forward else -original.offset))
            visits.append(model.Traversal(sid, forward))
            if piece:
                previous = visits[-2]
                connections.append(model.Connection(1000+100*role+piece, (
                    model.End(previous.section_id, int(previous.forward)),
                    model.End(sid, int(not forward)))))
            if piece == 1:
                # Interior anchors exercise per-Section dependency mapping.
                furniture.append(model.Placement(910+role, sid, 5., side='LEFT'))
        members.append(model.BundleMember(tuple(visits)))
    furniture.extend((replace(network.furniture[0], station=5.), network.furniture[1]))
    return replace(network, uuid=f'angular_chain_{sweep:g}_{int(variable)}_{mode}',
        sections=(*sections, network.sections[-1]),
        bundles=(replace(network.bundles[0], members=tuple(members)),),
        connections=tuple(connections), furniture=tuple(furniture))


def angular_chain_fixture_networks(model):
    """Both interpolation modes and island states, with reversed middle pieces."""
    return tuple(angular_chain_network(model, sweep, variable, mode,
                                       reversed_sections=(117,229))
                 for sweep in (180.,225.,270.) for variable in (False,True) for mode in (0,3))


def angular_dimensional_network(model, sweep=225., variable=False, mode=3, *, both=True,
                                width_delta=1., offset_delta=0., **kwargs):
    """Original curved Sections with finite width/offset changes at four seams."""
    from dataclasses import replace
    kwargs.setdefault('reversed_sections',(117,229))
    network = angular_chain_network(model,sweep,variable,mode,**kwargs)
    changed = {117,229} if both else {117}
    directions = {v.section_id:v.forward for member in network.bundles[0].members for v in member.path}
    return replace(network,uuid=f'angular_dimensional_{sweep:g}_{int(variable)}_{mode}',
        sections=tuple(replace(s,width=s.width+width_delta,
            offset=s.offset+offset_delta*(1 if directions.get(s.id,True) else -1))
            if s.id in changed else s for s in network.sections))


def angular_island_chain_network(model, mode=3, *, float32=False, spacing=1.,
                                  islands=3, reversed_sources=()):
    """Noncrossing members with repeated overlap and closed grass islands."""
    angles = np.linspace(-.75*np.pi, .75*np.pi, 49)
    phase = 4.*angles if islands == 3 else np.linspace(-np.pi, (2*islands-1)*np.pi, len(angles))
    radii = (np.full(len(angles), 80.), 87.+4.*np.cos(phase))
    sections, members = [], []
    for index, radius in enumerate(radii):
        sid = 17+12*index
        points = np.c_[radius*np.cos(angles), radius*np.sin(angles), np.zeros(len(angles))]
        if float32:
            points = points.astype(np.float32).astype(float)
        forward = sid not in reversed_sources
        if not forward:
            points = points[::-1]
        sections.append(model.Section(sid, 41+30*index, tuple(map(tuple, points)), width=6.,
            rise=float(np.float32(.08)) if float32 else .08, interpolation=mode))
        members.append(model.BundleMember((model.Traversal(sid, forward),)))
    sections.append(model.Section(99, 99, ((-100., -140., 0.), (100., -140., 0.))))
    return model.Network(f'angular_island_chain_{islands}_{mode}', tuple(sections),
        bundles=(model.Bundle(13, tuple(members)),),
        furniture=(model.Placement(901, 17, 40., side='LEFT', offset=-1.),
                   model.Placement(999, 99, 40.)), settings=model.Settings(cell_size=spacing))



def attached_angular_network(model,mode=3,*,sweep=270.,end_gap=16.,samples=17,radius=40.):
    t=np.linspace(0.,1.,samples)
    angle=np.deg2rad(sweep*(t-.5))
    gap=4.+(end_gap-4.)*(2.*t-1.)**2
    sections=[]
    for index,member_radius in enumerate((np.full(samples,radius),radius+gap)):
        points=np.c_[member_radius*np.cos(angle),member_radius*np.sin(angle),np.zeros(samples)]
        sections.append(model.Section(17+12*index,41+30*index,tuple(map(tuple,points)),
                                      width=6.,rise=.08,interpolation=mode))
    intersections=[]
    for end in (0,1):
        inner=np.asarray(sections[0].points[0 if end==0 else -1])
        outer=np.asarray(sections[1].points[0 if end==0 else -1])
        radial=(outer-inner)/np.linalg.norm(outer-inner)
        inner_stub=22. if radius>=40. else min(10.,.5*radius)
        cuts=(inner-inner_stub*radial,inner,outer,outer+22.*radial)
        base=100+10*end
        for index,(first,last) in enumerate(zip(cuts,cuts[1:])):
            sections.append(model.Section(base+index,500+end,(tuple(first),tuple(last)),
                                          width=5.,rise=.08,interpolation=mode))
        intersections.extend((
            model.Intersection(200+10*end,(model.End(17,end),model.End(base,1),model.End(base+1,0)),radius=2.),
            model.Intersection(201+10*end,(model.End(29,end),model.End(base+1,1),model.End(base+2,0)),radius=2.)))
    sections.append(model.Section(99,99,((-100.,-150.,0.),(100.,-150.,0.)),width=6.,rise=.08))
    identity=f'attached_angular_{sweep:g}_{mode}'+(f'_radius_{radius:g}' if radius!=40. else '')
    return model.Network(identity,tuple(sections),tuple(intersections),
        bundles=(model.Bundle(13,(model.BundleMember((model.Traversal(17),)),
                                  model.BundleMember((model.Traversal(29),)))),),
        furniture=(model.Placement(901,17,30.,side='LEFT',offset=-1.),
                   model.Placement(902,29,30.,side='RIGHT',offset=-1.),
                   model.Placement(903,101,.5*end_gap),model.Placement(999,99,40.)),
        settings=model.Settings(cell_size=1.))



def attached_angular_state_network(model, mode=3, *, state='separated'):
    from dataclasses import replace
    source=attached_angular_network(model,mode)
    outer=source.section_index[29];t=np.linspace(0.,1.,len(outer.points))
    angle=np.deg2rad(270.*(t-.5))
    if state=='separated':gap=10.+6.*(2.*t-1.)**2
    elif state=='multiple':gap=4.+12.*(2.*t-1.)**2+6.*np.sin(3.*np.pi*t)**2
    elif state=='overlap':gap=4.+12.*(2.*t-1.)**2
    else:raise ValueError('Unknown authored Bundle state')
    points=np.c_[(40.+gap)*np.cos(angle),(40.+gap)*np.sin(angle),np.zeros(len(t))]
    # The existing external junction copies stay authoritative and unchanged.
    points[[0,-1]]=np.asarray(outer.points)[[0,-1]]
    outer=replace(outer,points=tuple(map(tuple,points)))
    source=replace(source,uuid=f'attached_angular_{state}_{mode}',
        sections=tuple(outer if section.id==29 else section for section in source.sections))
    return source


def internal_angular_state_network(model, mode=3, *, state='separated'):
    """The original curved Bundle, interrupted by one full transverse Street."""
    from dataclasses import replace
    source = attached_angular_state_network(model, mode, state=state)
    middle = len(source.section_index[17].points)//2
    sections = []
    for section in source.sections:
        if section.id not in (17, 29):
            sections.append(section)
            continue
        sections.extend((replace(section, points=section.points[:middle+1]),
                         replace(section, id=section.id+100, points=section.points[middle:])))
    inner = np.asarray(source.section_index[17].points[middle])
    outer = np.asarray(source.section_index[29].points[middle])
    normal = (outer-inner)/np.linalg.norm(outer-inner)
    cuts = (inner-22.*normal, inner, outer, outer+22.*normal)
    sections.extend(model.Section(300+i, 800, (tuple(a), tuple(b)), width=5., rise=.08,
        interpolation=mode) for i, (a, b) in enumerate(zip(cuts, cuts[1:])))
    junctions = [replace(site, ends=tuple(model.End(end.section_id+100, 1)
        if end.end == 1 and end.section_id in (17, 29) else end for end in site.ends))
        for site in source.intersections]
    junctions.extend((
        model.Intersection(500, (model.End(17, 1), model.End(117, 0), model.End(300, 1), model.End(301, 0)), radius=2.),
        model.Intersection(501, (model.End(29, 1), model.End(129, 0), model.End(301, 1), model.End(302, 0)), radius=2.)))
    members = tuple(model.BundleMember((model.Traversal(sid), model.Traversal(sid+100))) for sid in (17, 29))
    return replace(source, uuid=source.uuid+'_internal', sections=tuple(sections), intersections=tuple(junctions),
        bundles=(replace(source.bundles[0], members=members, interior_section_ids=(301,)),))
