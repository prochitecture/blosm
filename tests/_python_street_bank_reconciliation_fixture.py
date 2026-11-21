"""Input-driven sharp source union with a real three-arm junction."""
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
from tested_python_streets.compounds import CompoundGroup
from tested_python_streets.compound_cycles import _port_footprint
from tested_python_streets.regions import RegionUnion


def fixture(*,angle=0.,shift=(0.,0.),reverse=False,neighbor=None):
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    groups=[];sections=[];sites=[];paths={};profiles={};routes={};data=[]
    for ordinal,centre in enumerate(((0.,0.),) if neighbor is None else ((0.,0.),neighbor)):
        arms=[];grids={};outward=[];mouths={};ends=[];sources=[]
        for role,(direction,width) in enumerate((((1.,0.),5.),((-1.,0.),5.),((.7,.7),2.7))):
            sid=3*ordinal+role+1
            points=(np.array(((0.,0.),20.*np.array(direction)/np.linalg.norm(direction)))+centre)@rotation.T+shift
            if reverse:points=points[::-1]
            section=Section(sid,sid,tuple((float(x),float(y),0.) for x,y in points),width=width,interpolation=0)
            sections.append(section);path=Path(section,.005);paths[sid]=path;profiles[sid]=Crown(section.rise,section.datum)
            route=Route((Piece(path,not reverse),));routes[sid]=route
            arm=Approach(route,int(reverse),profiles[sid].oriented(int(reverse)));arms.append(arm)
            station=15.;outward.append(station);end=End(sid,int(reverse));ends.append(end)
            mouths[sid,int(reverse)]=Mouth(path.length-station if reverse else station,16)
            grids[sid]=np.linspace(0.,path.length,81)
        site=Intersection(ordinal+1,tuple(ends),2.,(2.,)*3,('STYLE',)*3);sites.append(site)
        layout=JunctionLayout(tuple(arms),(),(),tuple(outward),mouths,4,8,np.zeros(2),[0.]*3,
                              mesh_strategy='boundary_chart')
        for arm,station in zip(arms,outward):
            sources.append(_port_footprint(Key('cycle_port',2*arm.source.id+arm.end),arm,station,8,.7,.00025,grids))
        loops=RegionUnion(max_segments=20000).build(source.region for source in sources)
        group=CompoundGroup((site.key,),(),tuple(sources),loops,layout)
        groups.append(group);data.append(grids)
    network=Network('original-bank-reconciliation',tuple(sections),tuple(sites),
                    settings=Settings(exaggeration=1.))
    products={**{Key('path',sid):path for sid,path in paths.items()},
              **{Key('profile',sid):profile for sid,profile in profiles.items()}}
    context=BuildContext(network,products,PatternLibrary())
    planner=SimpleNamespace(context=context,sites={s.key:s for s in sites},
        layouts={s.key:g.physical_layout for s,g in zip(sites,groups)},topology=RouteTopology(network))
    return planner,groups[0],data[0],tuple(groups[1:])
