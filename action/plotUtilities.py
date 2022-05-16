import matplotlib.pyplot as plt
from itertools import cycle

from numpy import False_
from lib.pygeos.geom import GeometryFactory, CAP_STYLE
from mathutils import Vector

# --------------------------------------------------------------------------------
def plotPoly(polygon,vertsOrder,color='k',width=1.,order=100):
    count = 0
    for v1,v2 in zip(polygon[:-1],polygon[1:]):
        plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
        if vertsOrder:
            plt.text(v1[0],v1[1],str(count))
        count += 1
        # plt.plot(v1[0],v1[1],'kx')
    v1, v2 = polygon[-1], polygon[0]
    plt.plot([v1[0],v2[0]],[v1[1],v2[1]],color,linewidth=width,zorder=order)
    if vertsOrder:
        plt.text(v1[0],v1[1],str(count))

def plotGeom(geom,vertsOrder,color='k',width=1.,order=100):
    coords = [(v.x,v.y) for v in geom.coords]
    plotPoly(coords,vertsOrder,color,width,order)

def plotGeomList(geom,vertsOrder,color='k',width=1.,order=100):
    coords = [(v.x,v.y) for v in geom]
    plotPoly(coords,vertsOrder,color,width,order)

def plotGeomPoly(poly,vertsOrder,color='k',width=1.,order=100):
    plotGeom(poly.exterior,vertsOrder,color,width,order)
    for interior in poly.interiors:
        plotGeom(interior,vertsOrder,'r',width,order)

def plotFillMultiPolyList(polyList):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    iterColor = cycle(colors)
    for poly in polyList:
        plotGeomList(poly,False,'k')
        x = [v.x for v in poly]
        y = [v.y for v in poly]
        plt.fill(x,y,next(iterColor))

def plotGeosWithHoles(geosPoly,vertsOrder,color='k',width=1.,order=100):
    poly = [(c.x,c.y) for c in geosPoly.exterior.coords]
    plotPoly(poly,vertsOrder,color,width,order)
    for ring in geosPoly.interiors:
        p = [(c.x,c.y) for c in ring.coords]
        plotPoly(p,vertsOrder,'r',width,order)

def plotGeosAll(geosPoly,vertsOrder,color='k',width=1.,order=100):
    if geosPoly.geom_type == 'Polygon':
        plotGeosWithHoles(geosPoly,vertsOrder,color,width,order)
    else:
        for geom in geosPoly.geoms:
            plotGeosWithHoles(geom,vertsOrder,color,width,order)

def iterColorCycle():
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    return cycle(colors)
    # import matplotlib.cm as cm
    # import numpy as np
    # colors = cm.prism(np.linspace(0, 1, 50))
    # return cycle(colors)

def plotCycles(cycleSegs):
    colorCycler = iterColorCycle()
    for cNr,segList in enumerate(cycleSegs):
        nodes = [v for s in segList for v in s.path[:-1] ] + [segList[0].s]
        area = sum( (v2[0]-v1[0])*(v2[1]+v1[1]) for v1,v2 in zip(nodes[:-1],nodes[1:]+[nodes[0]]) )
        if area>0.:
            plotPoly(nodes,False,'r',5,50)
        # else:
        #     plotPoly(nodes,False,'k:',1,100)
        x = [v[0] for v in nodes]
        y = [v[1] for v in nodes]
        plt.fill(x,y,next(colorCycler),alpha=1.0,zorder=10)
        x0 = sum(x)/len(x)
        y0 = sum(y)/len(y)
        plt.text(x0,y0,str(cNr),zorder=12)

def plotGeosString(string,vertsOrder,color='k',width=1.,order=100):
    points = [(v.x,v.y) for v in string.coords]
    x = [n[0] for n in points]
    y = [n[1] for n in points]
    plt.plot(x,y,color=color,zorder=order)

def plotGeosPolyFill(geosPoly,color='#ff0000',width=1.,alpha = 0.2, order=100):
    poly = [(v.x,v.y) for v in geosPoly.coords[:-1]]
    x = [n[0] for n in poly]
    y = [n[1] for n in poly]
    plt.fill(x,y,color,alpha = alpha,zorder = 500)

globalColorCycler = iterColorCycle()
def plotGeosAllFill(geosPoly,vertsOrder,color='k',width=1.,order=100):
    global globalColorCycler
    color = next(globalColorCycler)
    if geosPoly.geom_type == 'Polygon':
        plotGeosPolyFill(geosPoly,vertsOrder,color,width,order)
    else:
        for geom in geosPoly.geoms:
            plotGeosPolyFill(geom,vertsOrder,color,width,order)

def iniPlot():
    fig = plt.figure()
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)


def plotEnd():
    plt.gca().axis('equal')
    plt.show()

def plotEndRange(xlim,ylim,setLim=False):
    if setLim:
        plt.xlim(xlim)
        plt.ylim(ylim)
    plt.gca().set_aspect(1)
    plt.show()

def plotCCWFlow(poly,dist,inner=True,color='r',whichInner = -1):
    gf = GeometryFactory()
    gSeg = gf.createLinearRing([gf.createCoordinate(v) for v in poly] + [gf.createCoordinate(poly[0])])
    buffer = gSeg.buffer(dist,resolution=1,cap_style=CAP_STYLE.square)
    if inner:
        if hasattr(buffer, 'interiors'):
            for nr,inters in enumerate(buffer.interiors):
                if whichInner < 0 or whichInner == nr:
                    buf = inters.coords
                    for v1c,v2c in zip(buf,buf[1:]):
                        v1 = Vector((v1c.x,v1c.y))
                        v2 = Vector((v2c.x,v2c.y))
                        plt.gca().annotate("", xy=(v2.x,v2.y), xytext=(v1.x,v1.y), arrowprops=dict(arrowstyle="->",color=color))
    else:
        if hasattr(buffer, 'exterior'):
            buf = buffer.exterior.coords
            for v1c,v2c in zip(buf,buf[1:]):
                v1 = Vector((v1c.x,v1c.y))
                v2 = Vector((v2c.x,v2c.y))
                plt.gca().annotate("", xy=(v2.x,v2.y), xytext=(v1.x,v1.y), arrowprops=dict(arrowstyle="->",color=color))
    
def plotPolyLine(line,halfWidth=1.5,color='#ff0000',alpha=0.3):
    from lib.pygeos.shared import Coordinate
    if halfWidth > 0.:
        geosF = GeometryFactory()
        coords = [Coordinate(v[0],v[1]) for v in line.verts]
        geosString = geosF.createLineString(coords)
        buffer = geosString.buffer(halfWidth,resolution=3,cap_style=CAP_STYLE.flat)
        if buffer.geom_type == 'Polygon':
            plotGeosPolyFill(buffer,color,1,alpha,800)
            plotGeom(buffer,False,'w',1.,800)
        else:
            for geom in buffer.geoms:
                plotGeosPolyFill(geom,color,1,alpha,800)
                plotGeom(geom,False,color,1.,100)

    else:
        for v1,v2 in line.segments():
            plt.plot([v1.x,v2.x],[v1.y,v2.y],color)

# code fragments ----------------------------------------------------------------------------

# at end of sliceLargeCycles()
        # iniPlot()
        # plotGeosWithHoles(cycle,False)
        # plotEndRange([],[],False)

        # iniPlot()
        # r = [10, 13, 14, 8, 9, 11, 12, 2, 3, 4, 7, 0, 1, 5, 6]
        # for i, slic in enumerate(finalSlices):
        #     n = r.index(i)+2#int((1+1)/4)*4 + (i+1)%4 + 1
        #     plt.subplot(4,4,n)
        #     plt.gca().set_axis_off()
        #     # plt.title(str(i))
        #     plotGeosWithHoles(slic,False)
        #     plt.gca().axis('equal')
        #     plt.tight_layout()
        # plotEndRange([],[],False)

        # iniPlot()
        # for i, slic in enumerate(finalSlices):
        #     plotGeosWithHoles(slic,False)
        #     plt.gca().axis('equal')
        # plotEndRange([],[],False)

# __init__() of GraphCycle, before cleanCycle()
        # nodes = [n for s in segList for n in s.path[:-1]]
        # if len(nodes) > 2:
        #     iniPlot()
        #     plotCCWFlow(nodes,3,True,'r',0)
        #     plotPoly(nodes,False,'k')
        #     plotEndRange([-560,-200],[-230,150],False)

# __init__() of GraphCycle, after cleanCycle()
        # iniPlot()
        # plotCCWFlow(boundary,6,True,'r',0)
        # plotPoly(boundary,False,'k')

        # if holes:
        #     for hole in holes:
        #         plotCCWFlow(hole,6,False,'r',-1)
        #         plotPoly(hole,False,'b')
    
        # plotEndRange([-560,-200],[-230,150],False)

        # iniPlot()
        # plotPoly(boundary,False,'k:')
        # if holes:
        #     for hole in holes:
        #         plotPoly(hole,False,'b:')

        # if spurs:
        #     for spur in spurs:
        #         v1,v2 = spur
        #         plt.plot([v1[0],v2[0]],[v1[1],v2[1]],'r')

        # plotEndRange([-560,-200],[-230,150],False)

# at begin of createCyclePolygons() in road_polygons.py
        # from action.plotUtilities import plotCCWFlow, iniPlot, plotEnd
        # for segList in cycleSegs:
        #     nodes = [n for s in segList for n in s.path[:-1]]
            # if len(nodes) > 2:
            #     iniPlot()
            #     plotCCWFlow(nodes,3,True,'r',0)
            #     plotPoly(nodes,False,'k')
            #     plt.xlim([-560,-200])
            #     plt.ylim([-230,150])
            #     plt.gca().set_aspect(1)
            #     plt.show()




    
