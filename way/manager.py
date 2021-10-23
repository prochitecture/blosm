from . import Way
from defs.base.polyline import Polyline
from defs.way import allWayCategories, facadeVisibilityWayCategories


class WayManager:
    
    def __init__(self, data, app):
        self.id = "ways"
        self.data = data
        self.app = app
        
        # use the default layer class in the <app>
        self.layerClass = None
        
        # don't accept broken multipolygons
        self.acceptBroken = False
        
        self.layers = dict((category, []) for category in allWayCategories)
        
        self.renderers = []
        
        self.actions = []
        
        # <self.networkGraph>, <self.waySectionGraph> and <self.junctions> are set in an action,
        # for example <action.way_clustering.Way>
        self.networkGraph = self.waySectionGraph = self.junctions = None
        
        app.addManager(self)

    def parseWay(self, element, elementId):
        self.createWay(element)
    
    def parseRelation(self, element, elementId):
        return
    
    def createWay(self, element):
        # create a wrapper for the OSM way <element>
        way = Way(element, self)
        self.layers[way.category].append(way)
    
    def getAllWays(self):
        return (way for category in allWayCategories for way in self.layers[category])
    
    def getFacadeVisibilityWays(self):
        return (
            way for category in facadeVisibilityWayCategories for way in self.layers[category] \
            if not way.bridge and not way.tunnel
        )

    def getSelectedWays(self, categories):
        return ( way for category in categories for way in self.layers[category] )
    
    def process(self):
        for way in self.getAllWays():
            # <osm.projection> may not be availabe in the constructor of a <Way>, but it's needed to
            # create way-segments. That's why an additional <init(..)> method is needed.
            way.init(self)
        for action in self.actions:
            action.do(self)
    
    def addRenderer(self, renderer):
        self.renderers.append(renderer)
        self.app.addRenderer(renderer)
    
    def render(self):
        for renderer in self.renderers:
            renderer.render(self, self.data)
    
    def addAction(self, action):
        action.app = self.app
        self.actions.append(action)


class RoadPolygonsManager:
    
    def __init__(self, data, app):
        self.id = "road_polygons"
        self.data = data
        self.polylines = []
        self.connectedManagers = []
        self.actions = []
        
        # don't accept broken multipolygons
        self.acceptBroken = False
        
        app.addManager(self)
    
    def process(self):
        for polyline in self.polylines:
            polyline.init(self)
            
        for connectedManager in self.connectedManagers:
            self.polylines.extend(connectedManager.getPolylines())
        
        for action in self.actions:
            action.do(self)

    def parseWay(self, element, elementId):
        self.polylines.append(Polyline(element))
    
    def parseRelation(self, element, elementId):
        return
    
    def render(self):
        RoadPolygonsRenderer().render(self)


from mpl.renderer import Renderer
class RoadPolygonsRenderer(Renderer):
    """
    A temporary class
    """
    def render(self, manager):
        for polyline in manager.polylines:
            for edge in polyline.edges:
                self.mpl.ax.plot(
                    (edge.v1[0], edge.v2[0]),
                    (edge.v1[1], edge.v2[1]),
                    color="brown", linewidth=1.,
                )