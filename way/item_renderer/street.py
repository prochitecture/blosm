from . import ItemRenderer


class Street(ItemRenderer):
    
    def init(self, globalRenderer):
        super().init(globalRenderer)
        
        self.intersectionRenderer = globalRenderer.itemRenderers["Intersection"]
    
    def requestNodeGroups(self, nodeGroupNames):
        return

    def setNodeGroups(self, nodeGroups):
        return
    
    def setPolyline1ParamsForCorner(self, modifier, connector):
        modifier["Socket_3"] = connector.item.obj
        modifier["Socket_4"] = connector.leaving

    def setPolyline2ParamsForCorner(self, modifier, connector):
        modifier["Socket_6"] = connector.item.obj
        modifier["Socket_7"] = connector.leaving
    
    def renderNeighborIntersection(self, intersection, connector, index, modifier):
        street = connector.item
        order = intersection.order
            
        modifier[ self.intersectionRenderer.inputCenterlines[order][index][0] ] = street.obj
        modifier[ self.intersectionRenderer.inputWidths[order][index][0] ] = street.head.width if connector.leaving else street.tail.width
        modifier[ self.intersectionRenderer.inputLocations[order][index][0] ] = connector.leaving