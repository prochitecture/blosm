"""Three real source sites; a narrow consumed connector meets a wide bank.

Only the local exposed corner operands are needed for these discovery tests.
The source graph and its other incident arms remain explicit original inputs.
"""
from types import SimpleNamespace
import math
import numpy as np
from tested_python_streets.model import Section,End,Intersection,Network,Key,Settings
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route,Piece,RouteTopology
from tested_python_streets.profiles import Crown
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.intersections import Approach,JunctionLayout,Mouth
from tested_python_streets.corner_radii import CornerRadiusPolicy
from tested_python_streets.corner_patch import corner_patch
from tested_python_streets.compounds import CompoundGroup,SharedClaim
from tested_python_streets.compound_cycles import _port_footprint,_footprint
from tested_python_streets.source_curb_coordinates import CurbCoordinates
from tested_python_streets.regions import RegionUnion


def fixture(*,reverse=False,angle=0.,scope=False,through_radius=4.,through_mode='STYLE',consume=True):
    rotate=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    a=np.array((0.,0.));b=np.array((7.66,-6.42));c=np.array((2.13,8.83))
    controls={11:(a,a+(-7.,-43.),a+(-6.,-90.)),12:(a,b),13:(c,b),
        14:(b,b+(3.42,-22.02),b+(6.59,-172.89)),15:(a,a+(-20.,0.)),
        16:(c,c+(20.,0.)),17:(c,c+(0.,20.))}
    widths={11:5.4,12:2.7,13:10.8,14:10.8,15:5.4,16:10.8,17:10.8}
    sections=[];paths={};profiles={};arms={};grids={}
    for sid,points in controls.items():
        points=np.array(points)@rotate.T+(12.,-8.)
        if reverse:points=points[::-1]
        section=Section(sid,sid,tuple((float(x),float(y),0.) for x,y in points),width=widths[sid],
            interpolation=3,scope_ends=((0 if reverse else 1),) if scope and sid==12 else ())
        sections.append(section);path=Path(section,.005);paths[sid]=path
        profiles[sid]=Crown(.08,0.,(1 if reverse else -1) if sid in (13,14,16,17) else 0)
        grids[sid]=path.samples(0.,path.length,.25,tolerance=.00025,max_samples=20000)
        for end in (0,1):
            route=Route((Piece(path,end==0),));arms[sid,end]=Approach(route,end,profiles[sid].oriented(end))
    def end(sid,original):return End(sid,1-original if reverse else original)
    sites=(Intersection(1,(end(11,0),end(12,0),end(15,0)),2.,(2.,)*3,('STYLE',)*3),
        Intersection(2,(end(12,1),end(13,1),end(14,0)),through_radius,(through_radius,)*3,(through_mode,)*3),
        Intersection(3,(end(13,0),end(16,0),end(17,0)),4.,(4.,)*3,('STYLE',)*3))
    by_site={site.key:tuple(arms[e.section_id,e.end] for e in site.ends) for site in sites}
    first=arms[11,int(reverse)];connector=arms[12,int(reverse)]
    corner=CornerRadiusPolicy.solve(first,connector,2.,'STYLE')
    for arm,station in zip((first,connector),corner.stations):
        native=station if arm.end==0 else arm.length-station
        grids[arm.source.id]=np.unique(np.r_[grids[arm.source.id],native])
    source=corner_patch(Key('cycle_corner',0),corner,.7,.00025);source.layout.original_junction=sites[0].key
    sources=[source];external=(first,arms[14,int(reverse)])
    outward=(7.2,1.65);mouths={}
    for arm,at in zip(external,outward):
        sources.append(_port_footprint(Key('cycle_port',2*arm.source.id+arm.end),arm,at,8,.7,.00025,grids))
        mouths[arm.source.id,arm.end]=Mouth(arm.length-at if arm.end else at,16)
    for sid in (12,13):
        path=paths[sid];at=grids[sid];side=.5*path.section.width
        xy=np.vstack((path.point(at,side),path.point(at,-side)[::-1]))
        coordinates=CurbCoordinates(path,np.r_[at,at[::-1]],np.r_[np.full(len(at),side),np.full(len(at),-side)],
            np.r_[np.ones(len(at)),-np.ones(len(at))],maximum_chord_error=.00025)
        sources.append(_footprint(Key('cycle_strip',sid),xy,internal_edges=(len(at)-1,2*len(at)-1),source_coordinates=coordinates))
    layout=JunctionLayout(external,(corner,),(),outward,mouths,4,8,np.zeros(2),[0.]*2,mesh_strategy='boundary_chart')
    claims=tuple(SharedClaim(sid,sites[0].key if sid==12 else sites[2].key,sites[1].key,paths[sid].length,0.)
        for sid in ((12,13) if consume else (13,)))
    group=CompoundGroup(tuple(site.key for site in sites),claims,tuple(sources),
        RegionUnion(max_segments=20000).build(source.region for source in sources),layout)
    network=Network('owned-source-connector',tuple(sections),sites,settings=Settings(exaggeration=1.))
    products={**{Key('path',sid):path for sid,path in paths.items()},
              **{Key('profile',sid):profile for sid,profile in profiles.items()}}
    context=BuildContext(network,products,PatternLibrary())
    layouts={key:SimpleNamespace(approaches=value) for key,value in by_site.items()}
    planner=SimpleNamespace(context=context,sites={site.key:site for site in sites},layouts=layouts,topology=RouteTopology(network))
    return planner,group,grids
