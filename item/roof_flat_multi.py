from .roof_flat import RoofFlat


class RoofFlatMulti(RoofFlat):
    
    def __init__(self, footprint, firstVertIndex):
        super().__init__(footprint, firstVertIndex)
        self.innerPolygons = []
    
    def init(self):
        super().init()
        self.innerPolygons.clear()