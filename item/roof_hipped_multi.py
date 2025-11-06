from .roof_hipped import RoofHipped


class RoofHippedMulti(RoofHipped):
    
    def __init__(self, footprint, firstVertIndex):
        super().__init__(footprint)