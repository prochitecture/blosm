from mathutils import Vector
from itertools import tee, islice, cycle
import matplotlib.pyplot as plt
from . import Renderer
from lib.CompGeom.BoolPolyOps import _isectSegSeg
from lib.CompGeom.PolyLine import PolyLine
from way.item import IntConnector, Intersection, Section, SideLane, SymLane

def cyclePair(iterable):
    # iterable -> (p0,p1), (p1,p2), (p2, p3), ..., (pn, p0)
    prevs, nexts = tee(iterable)
    prevs = islice(cycle(prevs), len(iterable) - 1, None)
    return zip(prevs,nexts)

class StreetRenderer(Renderer):

    def __init__(self, debug):
        super().__init__()
        self.debug = debug
    
    def prepare(self):
        return

    def render(self, manager, data):
        def isMinorCategory(section):
            return  section.category in  ['footway', 'cycleway','service'] or \
                    ('service' in section.tags and \
                    section.tags['service']=='driveway')

        for Id,isect in enumerate(manager.majorIntersections):
            p = isect.location
            plt.plot(p[0],p[1],'ro',markersize=6,zorder=999,markeredgecolor='red', markerfacecolor='orange')
            # if self.debug:
            #     plt.text(p[0]+2,p[1]-2,str(isect.id),color='r',fontsize=10,zorder=130,ha='left', va='top', clip_on=True)

        # for Id,isect in enumerate(manager.minorIntersections):
        #     p = isect.location
        #     plt.plot(p[0],p[1],'rv',markersize=6,zorder=999,markeredgecolor='green', markerfacecolor='none')
        #     isect.leaving.plot('c',1,True)
        #     isect.arriving.plot('m',1,False)
        #     plt.text(p[0],p[1],'  '+str(isect.id),color='r',fontsize=6,zorder=130,ha='left', va='top', clip_on=True)

        # for src, dst, multKey, street in manager.waymap.edges(data='object',keys=True):
        #     allVertices = []
        #     for item in street.iterItems():
        #         if isinstance(item,Section):
        #             allVertices.extend(item.centerline)
        #             item.polyline.plot('b',1,'dotted')
        #     if len(allVertices):
        #         c = sum(allVertices,Vector((0,0))) / len(allVertices)
        #         plt.text(c[0],c[1],'S '+str(street.id),color='blue',fontsize=6,zorder=130,ha='left', va='top', clip_on=True)
        # return
        for street in manager.iterStreets():
        # for src, dst, multKey, street in manager.waymap.edges(data='object',keys=True):
            allVertices = []
            streetIsMinor = False
            for item in street.iterItems():
                # if street.id == 525:
                #     continue
                if isinstance(item,Section):
                    section = item
                    allVertices.extend(section.centerline)
                    if section.valid:
                        color = 'b' if isMinorCategory(section) else 'r'
                        width = 1 if isMinorCategory(section) else 2
                        style = 'dotted' if isMinorCategory(section) else 'solid'
                        section.polyline.plot(color,width,style)
                        allColor = color
                        if isMinorCategory(section):
                            streetIsMinor = True
                if isinstance(item,Intersection):
                    if not item.isMinor:
                        print('False non-minor')
                        continue
                    p = item.location
                    plt.plot(p[0],p[1],'rv',markersize=6,zorder=999,markeredgecolor='green', markerfacecolor='none')
                    if self.debug:
                        plt.text(p[0]+2,p[1]-2,'M '+str(item.id),color='g',fontsize=10,zorder=130,ha='left', va='top', clip_on=True)

                    if self.debug:
                        for conn in Intersection.iterate_from(item.leftHead):
                            line = conn.item.head.polyline if conn.leaving else conn.item.tail.polyline
                            vec = line[1]-line[0] if conn.leaving else line[-2]-line[-1]
                            vec = vec/vec.length
                            p0 = line[0] if conn.leaving else line[-1]
                            p1 = p0 +3*vec
                            plt.plot([p0[0],p1[0]], [p0[1],p1[1]], 'g')
                            plt.text(p1[0]+2,p1[1]-2,'C '+str(conn.id),color='k',fontsize=8,zorder=130,ha='left', va='top', clip_on=True)
                        for conn in Intersection.iterate_from(item.rightHead):
                            line = conn.item.head.polyline if conn.leaving else conn.item.tail.polyline
                            vec = line[1]-line[0] if conn.leaving else line[-2]-line[-1]
                            vec = vec/vec.length
                            p0 = line[0] if conn.leaving else line[-1]
                            p1 = p0 +3*vec
                            plt.plot([p0[0],p1[0]], [p0[1],p1[1]], 'r')
                            plt.text(p1[0]+2,p1[1]-2,'C '+str(conn.id),color='k',fontsize=8,zorder=130,ha='left', va='top', clip_on=True)

                if isinstance(item,SideLane):
                    p = item.location
                    plt.plot(p[0],p[1],'rs',markersize=6,zorder=999,markeredgecolor='green', markerfacecolor='cyan')

                if isinstance(item,SymLane):
                    p = item.location
                    plt.plot(p[0],p[1],'rs',markersize=6,zorder=999,markeredgecolor='green', markerfacecolor='cyan')

            # if street.id == 525:
            #     allVertices = list(dict.fromkeys(allVertices))
            #     longPoly = PolyLine(allVertices)
            #     longPoly.plot('g',4)

            if self.debug:
                color = 'cornflowerblue' if streetIsMinor else 'crimson'
                width = 8 if streetIsMinor else 10
                if len(allVertices):
                    c = sum(allVertices,Vector((0,0))) / len(allVertices)
                    plt.text(c[0]+2,c[1]-2,'S '+str(street.id),color=color,fontsize=width,zorder=130,ha='left', va='top', clip_on=True)

def plotPolygon(poly,vertsOrder,lineColor='k',fillColor='k',width=1.,fill=False,alpha = 0.2,order=100):
    x = [n[0] for n in poly] + [poly[0][0]]
    y = [n[1] for n in poly] + [poly[0][1]]
    if fill:
        plt.fill(x[:-1],y[:-1],color=fillColor,alpha=alpha,zorder = order)
    plt.plot(x,y,lineColor,linewidth=width,zorder=order)
    if vertsOrder:
        for i,(xx,yy) in enumerate(zip(x[:-1],y[:-1])):
            plt.text(xx,yy,str(i),fontsize=12, clip_on=True)

def plotWay(way,vertsOrder,lineColor='k',width=1.,order=100):
    x = [n[0] for n in way]
    y = [n[1] for n in way]
    plt.plot(x,y,lineColor,linewidth=width,zorder=order)
    if vertsOrder:
        for i,(xx,yy) in enumerate(zip(x,y)):
            plt.text(xx,yy,str(i),fontsize=12, clip_on=True)

def plotSideLaneWay(smallWay,wideWay,color):
    from lib.CompGeom.PolyLine import PolyLine
    smallLine = PolyLine(smallWay.centerline)
    wideLine = PolyLine(wideWay.centerline)
    smallLine.plotWithArrows('g',2)
    wideLine.plotWithArrows('r',2)

    wideLine = wideLine.parallelOffset(-wideWay.offset)
    poly = wideLine.buffer(wideWay.width/2.,wideWay.width/2.)
    plotPolygon(poly,False,color,color,1.,True,0.3,120)
    poly = smallLine.buffer(smallWay.width/2.,smallWay.width/2.)
    plotPolygon(poly,False,color,color,1.,True,0.3,120)

def plotEnd():
    plt.gca().axis('equal')
    plt.show()