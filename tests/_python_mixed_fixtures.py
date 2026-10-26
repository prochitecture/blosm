"""Model-independent authored mixed seams and the translated corner regression."""
from dataclasses import replace
import math
import numpy as np
from _python_angular_fixtures import internal_angular_state_network


def mixed_transition_networks(model):
    result=[]
    for modes in ((0,3),(3,0)):
        angle=math.radians(65.)
        sections=(model.Section(1,10,((-60.,0.,0.),(0.,0.,0.)),width=6.,interpolation=modes[0]),
            model.Section(2,10,((0.,0.,0.),(60.*math.cos(angle),60.*math.sin(angle),0.)),
                          width=10.,interpolation=modes[1]),
            model.Section(99,99,((1000.,-100.,0.),(1040.,-100.,0.))))
        result.append(model.Network(f'mixed_dimensional_{modes[0]}_{modes[1]}',sections,
            connections=(model.Connection(8,(model.End(1,1),model.End(2,0))),),
            furniture=(model.Placement(85,1,59.5,'RIGHT',-1.),model.Placement(99,99,20.,'RIGHT',-1.))))
    angle=math.radians(35.)
    base=result[0]
    result.append(replace(base,uuid='mixed_dimensional_unequal_rises',sections=(
        replace(base.sections[0],rise=.06,offset=.2),
        replace(base.sections[1],rise=.1,offset=-.1,
                points=((0.,0.,0.),(60.*math.cos(angle),60.*math.sin(angle),0.))),base.sections[2])))
    return tuple(result)


def uniform_mixed_networks(model):
    """Bent Exact/Smooth seams with unchanged dimensions on both sides."""
    return tuple(replace(network,uuid='uniform_'+network.uuid,
        sections=(network.sections[0],replace(network.sections[1],width=6.),network.sections[2]))
        for network in mixed_transition_networks(model)[:2])


def sharp_uniform_mixed_networks(model):
    """Short 110-degree mixed joins requiring bounded interior height changes."""
    angle=math.radians(110.)
    return tuple(replace(source,uuid=f'sharp_uniform_mixed_{modes[0]}_{modes[1]}',
        sections=(replace(source.sections[0],points=((-30.,0.,0.),(0.,0.,0.))),
            replace(source.sections[1],points=((0.,0.,0.),(30.*math.cos(angle),30.*math.sin(angle),0.))),
            source.sections[2]),
        furniture=(model.Placement(85,1,24.,'RIGHT',-2.),source.furniture[1]))
        for modes,source in zip(((0,3),(3,0)),uniform_mixed_networks(model)))


def uniform_mixed_bundle_network(model,modes=(0,3),*,overlap=False):
    from _python_street_fixtures import dimensional_bundle_fixture
    source=dimensional_bundle_fixture(model)
    sections=[]
    for section in source.sections:
        if section.id==99:
            sections.append(section)
            continue
        points=list(section.points)
        if section.id%2:
            x,y,z=points[-2];points[-2]=(x,y+.25,z)
        if overlap:
            shift=-1.5 if section.id<3 else 1.5
            points=[(x,y+shift,z) for x,y,z in points]
        sections.append(replace(section,points=tuple(points),width=6.,
                                interpolation=modes[0 if section.id%2 else 1]))
    return replace(source,uuid='uniform_mixed_bundle_overlap' if overlap else 'uniform_mixed_bundle',
                   sections=tuple(sections))


def translated_corner_seed_network(model):
    source=internal_angular_state_network(model,0,state='overlap')
    sections=[]
    for section in source.sections:
        # Match the accepted native fixture: translate the authored doubles,
        # then round once for Blender storage. Rounding before translation is
        # a different source and can select another valid quad partition.
        native=np.asarray(section.points,dtype=float)
        points=(native+(1000.,-500.,0.)).astype(np.float32).astype(float)
        sections.append(replace(section,points=tuple(map(tuple,points)),
            width=float(np.float32(section.width)),rise=float(np.float32(section.rise)),datum=1.7))
    return replace(source,uuid='translated_exact_corner_seed',sections=tuple(sections),
        settings=replace(source.settings,exaggeration=10.),furniture=(*source.furniture,
            model.Placement(904,301,2.,offset=-2.5),model.Placement(905,117,10.,offset=-1.)))
